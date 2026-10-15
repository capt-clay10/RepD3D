# -*- coding: utf-8 -*-
"""
interpolation.py - Spatial and temporal interpolation routines for RepD3D v2.

Contains three categories of interpolation:

1. TEMPORAL (1D):
   - fill_nans_with_interpolation: Fill NaN gaps in time series via linear interp
   - fill_zeros_with_interpolation: Fill zero gaps (NaN proxy) via linear interp

2. SPATIAL - STRUCTURED GRIDS:
   - bilinear_interpolation: 4-point bilinear on regular lat/lon grids
   - inspect_bilinear_stencil: expose the raw bilinear stencil for QC
   - griddata_interpolation: scipy.griddata wrapper for 2D structured grids

3. SPATIAL - UNSTRUCTURED GRIDS:
   - build_kdtree: Build KD-tree from cell centres
   - kdtree_idw_interpolation: Inverse-distance weighted interpolation
"""

from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata as _scipy_griddata
from scipy.spatial import cKDTree


def fill_zeros_with_interpolation(arr: np.ndarray) -> np.ndarray:
    """
    Replace zero values with linearly interpolated values from non-zero neighbours.

    WARNING: This function treats all zeros as missing data. The v1 behaviour is
    preserved for backward compatibility.
    """
    non_zero_idx = np.where(arr != 0)[0]
    zero_idx = np.where(arr == 0)[0]

    if non_zero_idx.size == 0:
        return arr

    non_zero_vals = arr[non_zero_idx]
    interp_func = interp1d(
        non_zero_idx,
        non_zero_vals,
        kind='linear',
        fill_value='extrapolate',
    )

    arr[zero_idx] = interp_func(zero_idx)
    return arr


def fill_nans_with_interpolation(data_array: np.ndarray) -> np.ndarray:
    """
    Fill NaN values in a 1D array via linear interpolation (v1 legacy).

    Two-step approach preserved from v1:
    1. Replace NaN with 0
    2. Interpolate zeros with fill_zeros_with_interpolation

    WARNING: this EXTRAPOLATES linearly beyond the first/last valid sample.
    A stencil corner that is wet for a few hours of a storm and dry otherwise
    comes back as a ramp of tens of metres. bilinear_interpolation() no longer
    uses it (2026-10-13); it is kept only for backward compatibility of the
    public API.
    """
    filled = data_array.copy()
    filled = np.nan_to_num(filled)
    return fill_zeros_with_interpolation(filled)


def _bridge_gaps_without_extrapolation(
    series: np.ndarray,
    stats: Optional[Dict[str, int]] = None,
) -> np.ndarray:
    """
    Fill NaN gaps of a 1-D time series without ever extrapolating.

    Internal gaps are bridged linearly between the neighbouring valid samples;
    leading/trailing gaps are held at the nearest valid sample. An all-NaN or
    all-valid series is returned unchanged. ``stats`` (optional dict) receives
    the counts ``bridged_internal``, ``held_leading`` and ``held_trailing``.
    """
    values = np.asarray(series, dtype=float).copy()
    valid = ~np.isnan(values)
    n_valid = int(valid.sum())
    if n_valid == 0 or n_valid == values.size:
        return values
    first = int(np.argmax(valid))
    last = int(values.size - 1 - np.argmax(valid[::-1]))
    leading = first
    trailing = values.size - 1 - last
    internal = int((~valid[first:last + 1]).sum())
    idx = np.arange(values.size, dtype=float)
    # np.interp clamps outside [xp[0], xp[-1]] -> leading/trailing gaps are held.
    values[~valid] = np.interp(idx[~valid], idx[valid], values[valid])
    if stats is not None:
        stats['bridged_internal'] = stats.get('bridged_internal', 0) + internal
        stats['held_leading'] = stats.get('held_leading', 0) + leading
        stats['held_trailing'] = stats.get('held_trailing', 0) + trailing
    return values


def _masked_bilinear(
    corners: Dict[str, np.ndarray],
    normalized_weights: Dict[str, float],
    stats: Optional[Dict[str, int]] = None,
) -> np.ndarray:
    """
    Bilinear average that, per timestep, uses only the corners that hold a
    value (mask-normalised weights). Timesteps where every contributing corner
    is NaN are bridged in time by _bridge_gaps_without_extrapolation. The
    result can therefore never leave the range of the source values.
    """
    keys = ('Q11', 'Q12', 'Q21', 'Q22')
    weights = np.array([float(normalized_weights[k]) for k in keys])
    stack = np.vstack([np.asarray(corners[k], dtype=float) for k in keys])
    wet = ~np.isnan(stack)
    numerator = (np.where(wet, stack, 0.0) * weights[:, None]).sum(axis=0)
    denominator = (wet * weights[:, None]).sum(axis=0)
    covered = denominator > 0.0
    out = np.full(stack.shape[1], np.nan)
    out[covered] = numerator[covered] / denominator[covered]
    if stats is not None:
        stats['renormalised_steps'] = int((covered & (denominator < 1.0 - 1e-9)).sum())
        stats['all_dry_steps'] = int((~covered).sum())
    if not covered.any():
        raise ValueError(
            "No source value at any timestep: every contributing stencil corner "
            "is NaN (land or permanently dry) for the whole period."
        )
    out = _bridge_gaps_without_extrapolation(out, stats)
    lo = float(np.nanmin(stack))
    hi = float(np.nanmax(stack))
    tol = 1e-9 * max(1.0, abs(lo), abs(hi))
    if float(out.min()) < lo - tol or float(out.max()) > hi + tol:
        raise RuntimeError(
            "Internal error: interpolated series left the range of its source "
            f"corners ({lo:g}..{hi:g}); interpolation must never extrapolate."
        )
    return out


# ── Bounded-memory access to gridded time series ────────────────────────────
#
# TrilaWatt and EasyGSH files store one complete 2-D field per on-disk chunk
# (chunk shape (1, nlat, nlon), deflate-compressed), so reading ONE cell's
# time series costs one chunk decompression per timestep, whatever the size
# of the request. Reading the four corners of a stencil one at a time thus
# decompresses every chunk four times, and reading N boundary points one at
# a time decompresses the whole window N times. Measured 2026-10-15 on the
# 2017+2018 TrilaWatt wave files (27649 steps): ~50 s per corner series
# through dask, 17 min per strict-QC point, ~35 h for a 90-section wave run.
#
# The helpers below keep the arithmetic untouched (the arrays handed to the
# interpolation are the same cells, same dtype, same order) and only change
# how the cells are read: a stencil is read as one 2x2 block, per-cell
# validity statistics for QC come from one pass over the time axis, and the
# corner series of many points are prefetched in one pass per variable. All
# passes walk the time axis a block at a time, so memory never depends on
# the file size. Blocks are always sliced from the array in its native
# dimension order and transposed in memory afterwards: transposing a lazily
# indexed backend array first routes the read through xarray's vectorized
# indexing path, which reads and copies ~10x the block (EasyGSH files are
# stored (lon, lat, time)).

DEFAULT_BLOCK_BYTES = 64 * 1024 * 1024     # one in-memory time block
DEFAULT_CACHE_BYTES = 512 * 1024 * 1024    # resident prefetched corner series
MAX_CHUNKS_PER_BLOCK = 4                   # dask chunks materialised per block


def _valid_mask(values: np.ndarray) -> np.ndarray:
    """'Valid' exactly as inspect_bilinear_stencil sees it: not NaN (floats), not NaT (times); other dtypes are always valid."""
    if np.issubdtype(values.dtype, np.floating) or np.issubdtype(values.dtype, np.complexfloating):
        return ~np.isnan(values)
    if np.issubdtype(values.dtype, np.datetime64) or np.issubdtype(values.dtype, np.timedelta64):
        return ~np.isnat(values)
    return np.ones(values.shape, dtype=bool)


def _axis_order(dims, leading) -> list:
    """Axis permutation that puts ``leading`` first (in that order), remaining dims in native order."""
    dims = tuple(dims)
    rest = [d for d in dims if d not in leading]
    return [dims.index(d) for d in tuple(leading) + tuple(rest)]


def iter_time_blocks(da, time_coord: str, target_bytes: int = DEFAULT_BLOCK_BYTES, cancel_check=None, dims_out=None):
    """
    Yield ``(start, stop, values)`` blocks of ``da`` along ``time_coord``.

    Every block is a numpy array with time as its leading axis, followed by
    the dims named in ``dims_out`` (if given) and then any remaining dims in
    native order. Dask-backed arrays are read along their own chunk
    boundaries: chunks are merged while the block stays within
    ``target_bytes`` (at most MAX_CHUNKS_PER_BLOCK chunks per block) and
    chunks larger than the target are split. Plain lazily-indexed arrays are
    read in slices of the target size. ``cancel_check`` (zero-argument
    callable) runs before every block so a GUI cancel is honoured mid-pass.
    """
    dims = tuple(da.dims)
    time_axis = dims.index(time_coord)
    n_steps = int(da.sizes[time_coord])
    if n_steps == 0:
        return
    leading = (time_coord,) + tuple(d for d in (dims_out or ()) if d != time_coord)
    order = _axis_order(dims, leading)
    slice_bytes = max(1, int(da.dtype.itemsize) * int(da.size // n_steps))
    target_steps = max(1, min(n_steps, int(target_bytes // slice_bytes)))
    bounds = []
    chunks = getattr(da, 'chunks', None)
    if chunks:
        start = None
        stop = 0
        count = 0
        edge = 0
        for size in chunks[time_axis]:
            c0, c1 = edge, edge + int(size)
            edge = c1
            if start is not None and (c1 - start > target_steps or count >= MAX_CHUNKS_PER_BLOCK):
                bounds.append((start, stop))
                start = None
                count = 0
            if start is None:
                start = c0
            stop = c1
            count += 1
        if start is not None:
            bounds.append((start, stop))
    else:
        bounds = [(s, min(n_steps, s + target_steps)) for s in range(0, n_steps, target_steps)]
    for b0, b1 in bounds:
        for s in range(b0, b1, target_steps):       # split anything larger than the target
            e = min(b1, s + target_steps)
            if cancel_check is not None:
                cancel_check()
            values = np.asarray(da.isel({time_coord: slice(s, e)}).values)
            yield s, e, np.transpose(values, order)


def cell_validity_stats(
    da,
    time_coord: str,
    lat_coord: str,
    lon_coord: str,
    target_bytes: int = DEFAULT_BLOCK_BYTES,
    cancel_check=None,
) -> Optional[Tuple[int, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Per-cell validity statistics of a (time, lat, lon) field from ONE pass.

    Returns ``(n_steps, valid_count, first_valid, last_valid)`` with the three
    arrays shaped (nlat, nlon); ``first_valid`` is n_steps and ``last_valid``
    is -1 for cells that are never valid. These four numbers are all that
    inspect_bilinear_stencil() derives per corner (valid fraction, leading,
    trailing and internal NaN counts), so QC verdicts built from them are
    identical to verdicts built from the corner series themselves.

    Returns None when the field is not indexed by exactly {time, lat, lon}
    (callers then fall back to reading stencils directly). A field without
    a time dimension counts as a single timestep.
    """
    dims = set(da.dims)
    if lat_coord not in dims or lon_coord not in dims:
        return None
    if time_coord not in dims:
        if dims != {lat_coord, lon_coord}:
            return None
        values = np.transpose(np.asarray(da.values), _axis_order(da.dims, (lat_coord, lon_coord)))
        valid = _valid_mask(values)
        return 1, valid.astype(np.int64), np.where(valid, 0, 1), np.where(valid, 0, -1)
    if dims != {time_coord, lat_coord, lon_coord}:
        return None

    n_steps = int(da.sizes[time_coord])
    shape = (int(da.sizes[lat_coord]), int(da.sizes[lon_coord]))
    count = np.zeros(shape, dtype=np.int64)
    first = np.full(shape, n_steps, dtype=np.int64)
    last = np.full(shape, -1, dtype=np.int64)
    blocks = iter_time_blocks(da, time_coord, target_bytes, cancel_check, dims_out=(time_coord, lat_coord, lon_coord))
    for start, stop, block in blocks:
        valid = _valid_mask(block)
        count += valid.sum(axis=0)
        seen = valid.any(axis=0)
        block_first = np.argmax(valid, axis=0) + start
        block_last = (valid.shape[0] - 1 - np.argmax(valid[::-1], axis=0)) + start
        first = np.where(seen & (first == n_steps), block_first, first)
        last = np.where(seen, block_last, last)
    return n_steps, count, first, last


def _read_stencil_block(ds, variable, lat_coord, lon_coord, lat_indices, lon_indices) -> Dict[str, np.ndarray]:
    """Read the four stencil corners as ONE (…, 2, 2) block instead of four reads."""
    lat_idx1, lat_idx2 = lat_indices
    lon_idx1, lon_idx2 = lon_indices
    block = ds[variable].isel({
        lat_coord: slice(lat_idx1, lat_idx2 + 1),
        lon_coord: slice(lon_idx1, lon_idx2 + 1),
    })
    dims = tuple(block.dims)
    rest = [d for d in dims if d not in (lat_coord, lon_coord)]
    values = np.transpose(np.asarray(block.values), [dims.index(d) for d in rest + [lat_coord, lon_coord]])
    # At the grid edge inspect_bilinear_geometry clips idx2 onto idx1; the
    # slice then has length 1 and both corners are the same cell, exactly as
    # the former per-corner .sel() calls returned the same series twice.
    row2 = values.shape[-2] - 1
    col2 = values.shape[-1] - 1
    return {
        'Q11': values[..., 0, 0],
        'Q12': values[..., 0, col2],
        'Q21': values[..., row2, 0],
        'Q22': values[..., row2, col2],
    }


class CornerSeriesCache:
    """
    Prefetched stencil-corner time series for a known set of sampling points.

    Build it once per extraction from the (post-QC) sampling points and hand
    it to bilinear_interpolation(cache=...). The planned cells are loaded in
    batches sized so that ALL variables of a batch fit in ``max_bytes``
    together (the extractors loop point-major over several variables, so a
    per-variable sizing would thrash). Each (variable, batch) is one pass
    over the time axis: the bounding box of the batch is read block by block
    and the wanted cells gathered from every block, so N points cost a few
    passes over the window instead of N x 4 reads. Corners of points that
    were not planned are read together as one small ad-hoc unit (one pass
    for the four corners, the same cost as a direct block read). Resident
    data is capped at ``max_bytes``; the oldest units are dropped first.
    """

    def __init__(self, ds, points, variables, lat_coord: str, lon_coord: str, time_coord: str,
                 max_bytes: int = DEFAULT_CACHE_BYTES, block_bytes: int = DEFAULT_BLOCK_BYTES,
                 cancel_check=None, logger=None):
        self._ds = ds
        self._lat_coord = lat_coord
        self._lon_coord = lon_coord
        self._time_coord = time_coord
        self._max_bytes = int(max_bytes)
        self._block_bytes = int(block_bytes)
        self._cancel_check = cancel_check
        self._logger = logger
        self._variables = [v for v in dict.fromkeys(variables) if v and self._readable(v)]
        self._cells: list = []
        self._position: Dict[Tuple[int, int], int] = {}
        for lat, lon in points:
            geometry = inspect_bilinear_geometry(lat, lon, ds, lat_coord=lat_coord, lon_coord=lon_coord)
            lat_idx1, lat_idx2 = geometry['lat_indices']
            lon_idx1, lon_idx2 = geometry['lon_indices']
            for cell in ((lat_idx1, lon_idx1), (lat_idx1, lon_idx2), (lat_idx2, lon_idx1), (lat_idx2, lon_idx2)):
                self._plan(cell)
        self._n_steps = int(ds.sizes.get(time_coord, 1))
        bytes_per_cell_step = sum(int(ds[v].dtype.itemsize) for v in self._variables) or 8
        self._batch_size = max(1, self._max_bytes // (max(1, self._n_steps) * bytes_per_cell_step))
        self._data: Dict[Tuple[str, int, int], np.ndarray] = {}
        self._resident: list = []        # [(unit_key, variable, cells, nbytes)] oldest first
        self._resident_bytes = 0
        self.passes = 0

    # -- planning ---------------------------------------------------------
    def _readable(self, variable: str) -> bool:
        if variable not in self._ds:
            return False
        return set(self._ds[variable].dims) == {self._time_coord, self._lat_coord, self._lon_coord}

    def _plan(self, cell: Tuple[int, int]) -> int:
        cell = (int(cell[0]), int(cell[1]))
        if cell not in self._position:
            self._position[cell] = len(self._cells)
            self._cells.append(cell)
        return self._position[cell]

    @property
    def n_cells(self) -> int:
        return len(self._cells)

    @property
    def n_batches(self) -> int:
        return max(1, -(-len(self._cells) // self._batch_size))

    def has(self, variable: str) -> bool:
        return variable in self._variables

    # -- access -----------------------------------------------------------
    def get(self, variable: str, lat_index: int, lon_index: int) -> np.ndarray:
        cell = (int(lat_index), int(lon_index))
        key = (variable,) + cell
        if key not in self._data:
            if cell in self._position:
                batch_id = self._position[cell] // self._batch_size
                self._load(variable, ('batch', variable, batch_id),
                           self._cells[batch_id * self._batch_size:(batch_id + 1) * self._batch_size])
            else:
                self._load(variable, ('adhoc', variable, cell), [cell])
        return self._data[key]

    def get_corners(self, variable: str, lat_indices, lon_indices) -> Dict[str, np.ndarray]:
        """The four stencil corners; corners outside the plan are read together in one pass."""
        lat_idx1, lat_idx2 = (int(i) for i in lat_indices)
        lon_idx1, lon_idx2 = (int(j) for j in lon_indices)
        corners = {
            'Q11': (lat_idx1, lon_idx1),
            'Q12': (lat_idx1, lon_idx2),
            'Q21': (lat_idx2, lon_idx1),
            'Q22': (lat_idx2, lon_idx2),
        }
        missing = [c for c in dict.fromkeys(corners.values())
                   if (variable,) + c not in self._data and c not in self._position]
        if missing:
            self._load(variable, ('adhoc', variable, tuple(missing)), missing)
        return {key: self.get(variable, *cell) for key, cell in corners.items()}

    def _drop(self, entry) -> None:
        _unit_key, variable, cells, nbytes = entry
        for cell in cells:
            self._data.pop((variable,) + cell, None)
        self._resident.remove(entry)
        self._resident_bytes -= nbytes

    def _load(self, variable: str, unit_key, cells) -> None:
        for entry in [e for e in self._resident if e[0] == unit_key]:   # never keep two copies of a unit
            self._drop(entry)
        ii = np.array([c[0] for c in cells], dtype=int)
        jj = np.array([c[1] for c in cells], dtype=int)
        i0, j0 = int(ii.min()), int(jj.min())
        da = self._ds[variable].isel({
            self._lat_coord: slice(i0, int(ii.max()) + 1),
            self._lon_coord: slice(j0, int(jj.max()) + 1),
        })
        out = np.empty((self._n_steps, len(cells)), dtype=da.dtype)
        blocks = iter_time_blocks(da, self._time_coord, self._block_bytes, self._cancel_check,
                                  dims_out=(self._time_coord, self._lat_coord, self._lon_coord))
        for start, stop, block in blocks:
            out[start:stop] = block[:, ii - i0, jj - j0]
        self.passes += 1
        nbytes = int(out.nbytes)
        while self._resident and self._resident_bytes + nbytes > self._max_bytes:
            self._drop(self._resident[0])
        for k, cell in enumerate(cells):
            self._data[(variable,) + cell] = out[:, k]
        self._resident.append((unit_key, variable, list(cells), nbytes))
        self._resident_bytes += nbytes
        if self._logger is not None:
            if unit_key[0] == 'batch':
                where = f"batch {unit_key[2] + 1}/{self.n_batches}"
            else:
                where = "unplanned cell(s), read directly"
            self._logger.info(
                f"[prefetch]   {variable}: {len(cells)} cell(s) x {self._n_steps} steps read in one pass ({where})"
            )


def _build_bilinear_stencil(
    lat: float,
    lon: float,
    ds,
    variable: str,
    lat_coord: str = 'lat',
    lon_coord: str = 'lon',
    cache: Optional[CornerSeriesCache] = None,
) -> Dict[str, Any]:
    """Build the 4-corner stencil used by bilinear_interpolation()."""
    geometry = inspect_bilinear_geometry(
        lat, lon, ds, lat_coord=lat_coord, lon_coord=lon_coord
    )

    if cache is not None and cache.has(variable):
        corners = cache.get_corners(variable, geometry['lat_indices'], geometry['lon_indices'])
    elif lat_coord in ds[variable].dims and lon_coord in ds[variable].dims:
        corners = _read_stencil_block(
            ds, variable, lat_coord, lon_coord,
            geometry['lat_indices'], geometry['lon_indices'],
        )
    else:
        # Field not indexed by the lat/lon dimensions: the original per-corner
        # label selection (no shipped preset takes this path).
        lat1, lat2 = geometry['lat_coords']
        lon1, lon2 = geometry['lon_coords']
        corners = {
            'Q11': ds.sel({lat_coord: lat1, lon_coord: lon1})[variable].values,
            'Q12': ds.sel({lat_coord: lat1, lon_coord: lon2})[variable].values,
            'Q21': ds.sel({lat_coord: lat2, lon_coord: lon1})[variable].values,
            'Q22': ds.sel({lat_coord: lat2, lon_coord: lon2})[variable].values,
        }

    geometry = dict(geometry)
    geometry['corners'] = corners
    return geometry


def _normalized_bilinear_weights(
    weights: Dict[str, float],
    total_weight: float,
) -> Tuple[bool, Dict[str, float]]:
    """Normalize bilinear weights and special-case exact-grid sampling."""
    point_on_grid = bool(np.isclose(total_weight, 0.0))
    if point_on_grid:
        return point_on_grid, {'Q11': 1.0, 'Q12': 0.0, 'Q21': 0.0, 'Q22': 0.0}

    return point_on_grid, {
        key: abs(val / total_weight)
        for key, val in weights.items()
    }


def inspect_bilinear_geometry(
    lat: float,
    lon: float,
    ds,
    lat_coord: str = 'lat',
    lon_coord: str = 'lon',
) -> Dict[str, Any]:
    """Inspect the bilinear geometry without reading a variable time series."""
    lat_vals = ds.coords[lat_coord].values
    lon_vals = ds.coords[lon_coord].values

    lat_outside = bool(lat < lat_vals[0] or lat > lat_vals[-1])
    lon_outside = bool(lon < lon_vals[0] or lon > lon_vals[-1])

    lat_idx1 = int(np.clip(np.searchsorted(lat_vals, lat) - 1, 0, len(lat_vals) - 1))
    lat_idx2 = int(np.clip(lat_idx1 + 1, 0, len(lat_vals) - 1))
    lon_idx1 = int(np.clip(np.searchsorted(lon_vals, lon) - 1, 0, len(lon_vals) - 1))
    lon_idx2 = int(np.clip(lon_idx1 + 1, 0, len(lon_vals) - 1))

    lat1 = float(lat_vals[lat_idx1])
    lat2 = float(lat_vals[lat_idx2])
    lon1 = float(lon_vals[lon_idx1])
    lon2 = float(lon_vals[lon_idx2])

    weights = {
        'Q11': float((lat2 - lat) * (lon2 - lon)),
        'Q21': float((lat - lat1) * (lon2 - lon)),
        'Q12': float((lat2 - lat) * (lon - lon1)),
        'Q22': float((lat - lat1) * (lon - lon1)),
    }
    total_weight = float(sum(weights.values()))
    point_on_grid, normalized_weights = _normalized_bilinear_weights(
        weights, total_weight
    )

    return {
        'lat': float(lat),
        'lon': float(lon),
        'outside': lat_outside or lon_outside,
        'lat_outside': lat_outside,
        'lon_outside': lon_outside,
        'point_on_grid': point_on_grid,
        'lat_indices': (lat_idx1, lat_idx2),
        'lon_indices': (lon_idx1, lon_idx2),
        'lat_coords': (lat1, lat2),
        'lon_coords': (lon1, lon2),
        'weights': weights,
        'normalized_weights': normalized_weights,
        'total_weight': total_weight,
    }


def inspect_bilinear_stencil(
    lat: float,
    lon: float,
    ds,
    variable: str,
    lat_coord: str = 'lat',
    lon_coord: str = 'lon',
    cache: Optional[CornerSeriesCache] = None,
) -> Dict[str, Any]:
    """
    Inspect the raw bilinear interpolation stencil for a single point.

    Returns metadata about each source corner so callers can decide whether
    interpolation would rely on all-NaN cells or long edge extrapolation.
    """
    stencil = _build_bilinear_stencil(
        lat, lon, ds, variable, lat_coord=lat_coord, lon_coord=lon_coord, cache=cache
    )

    total_weight = stencil['total_weight']
    point_on_grid = stencil['point_on_grid']
    normalized = stencil['normalized_weights']

    corner_summary: Dict[str, Dict[str, Any]] = {}
    for key, arr in stencil['corners'].items():
        values = np.asarray(arr)
        n_steps = int(len(values))
        valid_mask = ~np.isnan(values)
        valid_count = int(valid_mask.sum())

        if valid_count == 0:
            first_valid = None
            last_valid = None
            leading = n_steps
            trailing = n_steps
            internal = 0
        else:
            valid_idx = np.where(valid_mask)[0]
            first_valid = int(valid_idx[0])
            last_valid = int(valid_idx[-1])
            leading = first_valid
            trailing = n_steps - last_valid - 1
            internal = int(np.isnan(values[first_valid:last_valid + 1]).sum())

        corner_summary[key] = {
            'n_steps': n_steps,
            'valid_count': valid_count,
            'valid_fraction': (valid_count / n_steps) if n_steps else 0.0,
            'all_nan': (valid_count == 0),
            'first_valid_index': first_valid,
            'last_valid_index': last_valid,
            'leading_nan_count': leading,
            'trailing_nan_count': trailing,
            'internal_nan_count': internal,
            'normalized_weight': float(normalized[key]),
        }

    stencil['point_on_grid'] = point_on_grid
    return {
        'lat': stencil['lat'],
        'lon': stencil['lon'],
        'outside': stencil['outside'],
        'lat_outside': stencil['lat_outside'],
        'lon_outside': stencil['lon_outside'],
        'point_on_grid': point_on_grid,
        'lat_indices': stencil['lat_indices'],
        'lon_indices': stencil['lon_indices'],
        'lat_coords': stencil['lat_coords'],
        'lon_coords': stencil['lon_coords'],
        'weights': dict(stencil['weights']),
        'total_weight': float(total_weight),
        'corners': corner_summary,
    }


def bilinear_interpolation(
    lat: float,
    lon: float,
    ds,
    variable: str,
    lat_coord: str = 'lat',
    lon_coord: str = 'lon',
    stats: Optional[Dict[str, int]] = None,
    cache: Optional[CornerSeriesCache] = None,
) -> np.ndarray:
    """
    Bilinear interpolation on a regular lat/lon grid.

    Finds the 4 surrounding grid cells and computes an area-weighted average.
    ``cache`` (optional CornerSeriesCache) serves prefetched corner series;
    without it the stencil is read as one 2x2 block from ``ds``.

    NaN handling (2026-10-13): a corner that is NaN at some timesteps (dry
    intertidal cell, storm-only wet cell, missing record) is simply left out
    of the weighting at those timesteps and the remaining corners are
    renormalised. Only when every contributing corner is NaN is the value
    bridged in time, and then without extrapolation (internal gaps linear,
    edges held). The v1 behaviour (fill_nans_with_interpolation, which
    extrapolated a storm-only corner into a ramp of tens of metres) is gone.
    Stencils without any NaN take exactly the v1 code path, so their output
    is byte-identical to before.

    ``stats`` (optional dict) receives renormalised_steps / all_dry_steps /
    bridged_internal / held_leading / held_trailing for the caller's log.
    """
    stencil = _build_bilinear_stencil(
        lat, lon, ds, variable, lat_coord=lat_coord, lon_coord=lon_coord, cache=cache
    )
    # Native dtype on purpose (TrilaWatt/EasyGSH are float32): the fully wet
    # path below must stay bit-identical to v1, which never upcast.
    Q11 = np.array(stencil['corners']['Q11'], copy=True)
    Q12 = np.array(stencil['corners']['Q12'], copy=True)
    Q21 = np.array(stencil['corners']['Q21'], copy=True)
    Q22 = np.array(stencil['corners']['Q22'], copy=True)

    if any(np.isnan(q).any() for q in (Q11, Q12, Q21, Q22)):
        return _masked_bilinear(
            {'Q11': Q11, 'Q12': Q12, 'Q21': Q21, 'Q22': Q22},
            stencil['normalized_weights'],
            stats,
        )

    w11 = stencil['weights']['Q11']
    w21 = stencil['weights']['Q21']
    w12 = stencil['weights']['Q12']
    w22 = stencil['weights']['Q22']
    total_weight = stencil['total_weight']

    if np.isclose(total_weight, 0.0):
        return Q11

    return (Q11 * w11 + Q21 * w21 + Q12 * w12 + Q22 * w22) / total_weight


def build_kdtree(
    source_lats: np.ndarray,
    source_lons: np.ndarray,
) -> cKDTree:
    """
    Build a KD-tree from unstructured grid cell centres.
    """
    coords = np.column_stack([source_lats, source_lons])
    return cKDTree(coords)


def kdtree_idw_interpolation(
    tree: cKDTree,
    source_values: np.ndarray,
    target_lats: np.ndarray,
    target_lons: np.ndarray,
    k: int = 6,
    power: float = 2.0,
) -> np.ndarray:
    """
    Inverse-distance weighted interpolation using KD-tree.
    """
    target_shape = target_lats.shape
    target_coords = np.column_stack([target_lats.ravel(), target_lons.ravel()])

    distances, indices = tree.query(target_coords, k=k)
    distances = np.maximum(distances, 1e-10)

    weights = 1.0 / distances ** power
    weights /= weights.sum(axis=1, keepdims=True)

    result = np.sum(weights * source_values[indices], axis=1)
    return result.reshape(target_shape)


def griddata_interpolation(
    source_x: np.ndarray,
    source_y: np.ndarray,
    source_values: np.ndarray,
    target_x: np.ndarray,
    target_y: np.ndarray,
    method: str = 'linear',
) -> np.ndarray:
    """
    Interpolate scattered 2D data onto a target grid using scipy.griddata.
    """
    return _scipy_griddata(
        (source_x.flatten(), source_y.flatten()),
        source_values.flatten(),
        (target_x, target_y),
        method=method,
    )


if __name__ == '__main__':
    print("=== interpolation.py self-test ===\n")
    all_pass = True

    def check(name, result, expected, tolerance=0.01):
        global all_pass
        if isinstance(result, (float, np.floating)):
            ok = abs(float(result) - float(expected)) < tolerance
        else:
            ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}] got {result!r}")
        else:
            print(f"  {name}: [{status}]")

    arr = np.array([1.0, 0.0, 3.0])
    filled = fill_zeros_with_interpolation(arr.copy())
    check("fill_zeros midpoint", filled[1], 2.0)

    nan_arr = np.array([1.0, np.nan, 3.0])
    nan_filled = fill_nans_with_interpolation(nan_arr)
    check("fill_nans midpoint", nan_filled[1], 2.0)

    try:
        import xarray as xr

        lats = np.array([53.0, 54.0, 55.0])
        lons = np.array([8.0, 9.0, 10.0])
        times = np.arange(5)

        wl_data = np.zeros((3, 3, 5))
        for i, lat_val in enumerate(lats):
            wl_data[i, :, :] = lat_val - 53.0

        ds = xr.Dataset(
            {'water_level': (['lat', 'lon', 'time'], wl_data)},
            coords={'lat': lats, 'lon': lons, 'time': times},
        )

        result_bilinear = bilinear_interpolation(53.5, 9.0, ds, 'water_level')
        check("bilinear at (53.5, 9.0)", result_bilinear[0], 0.5)
        check("bilinear length", len(result_bilinear), 5, tolerance=0)

        geometry = inspect_bilinear_geometry(53.5, 9.0, ds)
        check("geometry outside", geometry['outside'], False, tolerance=0)
        check("geometry Q12 weight", geometry['normalized_weights']['Q12'], 0.5)
        check("geometry Q22 weight", geometry['normalized_weights']['Q22'], 0.5)

        stencil = inspect_bilinear_stencil(53.5, 9.0, ds, 'water_level')
        check("stencil outside", stencil['outside'], False, tolerance=0)
        check("stencil Q11 valid", stencil['corners']['Q11']['all_nan'], False, tolerance=0)

        # ── NaN-aware bilinear: no extrapolation from partially dry corners ──
        # Point (53.5, 8.5): all four corners carry 25 % weight.
        tide = np.array([1.0, 0.5, -0.5, -1.0, 0.0, 1.0])
        field = np.zeros((3, 3, 6))
        field[:, :, :] = tide
        field[1, 0, :] = tide + 0.4                       # Q21 (lat 54, lon 8) slightly higher
        storm = np.full(6, np.nan)
        storm[3] = 2.0                                    # Q22 wet at ONE step only (storm)
        field[1, 1, :] = storm
        ds_nan = xr.Dataset(
            {'wl': (['lat', 'lon', 'time'], field)},
            coords={'lat': lats, 'lon': lons, 'time': np.arange(6)},
        )
        stats = {}
        got = bilinear_interpolation(53.5, 8.5, ds_nan, 'wl', stats=stats)
        expect_dry = (tide + tide + (tide + 0.4)) / 3.0   # three wet corners, renormalised
        check("storm-only corner ignored while dry (step 0)", got[0], expect_dry[0], tolerance=1e-12)
        check("storm-only corner ignored while dry (step 5)", got[5], expect_dry[5], tolerance=1e-12)
        check("storm-only corner used when wet (step 3)",
              got[3], (tide[3] + tide[3] + (tide[3] + 0.4) + 2.0) / 4.0, tolerance=1e-12)
        check("no value outside the source range", float(np.nanmax(np.abs(got))) <= 2.0, True, tolerance=0)
        check("stats: renormalised steps", stats['renormalised_steps'], 5, tolerance=0)
        check("stats: no all-dry steps", stats['all_dry_steps'], 0, tolerance=0)

        # Every corner dry at some steps -> bridged in time, edges held, never extrapolated.
        field2 = np.full((3, 3, 6), np.nan)
        field2[0, 0, :] = [np.nan, np.nan, 1.0, np.nan, 3.0, np.nan]
        field2[1, 0, :] = [np.nan, np.nan, 1.0, np.nan, 3.0, np.nan]
        ds_dry = xr.Dataset(
            {'wl': (['lat', 'lon', 'time'], field2)},
            coords={'lat': lats, 'lon': lons, 'time': np.arange(6)},
        )
        stats2 = {}
        got2 = bilinear_interpolation(53.5, 8.0, ds_dry, 'wl', stats=stats2)   # on the lon=8 grid line
        check("all-dry leading gap held (no backward ramp)", got2[0], 1.0, tolerance=1e-12)
        check("all-dry internal gap bridged linearly", got2[3], 2.0, tolerance=1e-12)
        check("all-dry trailing gap held (no forward ramp)", got2[5], 3.0, tolerance=1e-12)
        check("stats: all-dry steps counted", stats2['all_dry_steps'], 4, tolerance=0)
        check("stats: held edges counted", stats2['held_leading'] + stats2['held_trailing'], 3, tolerance=0)

        # Fully wet stencil: unchanged v1 arithmetic (byte-identical path).
        got3 = bilinear_interpolation(53.5, 9.0, ds, 'water_level')
        check("fully wet stencil unchanged", bool(np.array_equal(got3, result_bilinear)), True, tolerance=0)

        # All corners NaN for the whole period must fail loudly, not return NaN.
        try:
            bilinear_interpolation(54.5, 9.5, ds_dry, 'wl')
            check("all-NaN stencil raises", "no error", "ValueError", tolerance=0)
        except ValueError:
            check("all-NaN stencil raises", "ValueError", "ValueError", tolerance=0)

        # ── Bounded-memory access: one-pass statistics, block read, prefetch cache ──
        bm_rng = np.random.default_rng(7)
        bm_nt, bm_nlat, bm_nlon = 37, 5, 6
        bm_data = bm_rng.normal(size=(bm_nt, bm_nlat, bm_nlon)).astype('float32')
        bm_data[bm_rng.random(bm_data.shape) < 0.3] = np.nan
        bm_data[:, 0, 0] = np.nan                               # never valid
        bm_data[:5, 1, 1] = np.nan
        bm_data[-7:, 1, 1] = np.nan                             # leading + trailing gaps
        bm_data[:, 2, 2] = 1.0                                  # always valid
        bm_data[:-1, 3, 3] = np.nan
        bm_data[-1, 3, 3] = 1.0                                 # valid only at the last step
        bm_lat = np.linspace(0, 1, bm_nlat)
        bm_lon = np.linspace(0, 1, bm_nlon)
        bm_sets = {
            'numpy (time, lat, lon)': xr.Dataset(
                {'v': (['time', 'lat', 'lon'], bm_data)},
                coords={'time': np.arange(bm_nt), 'lat': bm_lat, 'lon': bm_lon}),
            'numpy (lon, lat, time)': xr.Dataset(                 # EasyGSH layout
                {'v': (['lon', 'lat', 'time'], bm_data.transpose(2, 1, 0))},
                coords={'time': np.arange(bm_nt), 'lat': bm_lat, 'lon': bm_lon}),
        }
        try:
            bm_sets['dask (time, lat, lon)'] = bm_sets['numpy (time, lat, lon)'].chunk({'time': 7})
        except Exception:   # dask not installed: the numpy-backed checks still run
            pass
        bm_points = [(0.12, 0.33), (0.0, 0.0), (1.0, 1.0), (0.5, 0.99), (0.37, 0.61), (1.3, 0.5), (0.25, 0.25)]

        def bm_brute(ds_bm):
            arr = ds_bm['v'].transpose('time', 'lat', 'lon').values
            valid = ~np.isnan(arr)
            idx = np.arange(arr.shape[0])[:, None, None]
            return (arr.shape[0], valid.sum(axis=0),
                    np.where(valid, idx, arr.shape[0]).min(axis=0), np.where(valid, idx, -1).max(axis=0))

        def bm_old_corners(ds_bm, la, lo):
            g = inspect_bilinear_geometry(la, lo, ds_bm)
            cells = [(a, b) for a in g['lat_coords'] for b in g['lon_coords']]
            return {k: ds_bm.sel(lat=a, lon=b)['v'].values for k, (a, b) in zip(('Q11', 'Q12', 'Q21', 'Q22'), cells)}

        for bm_label, ds_bm in bm_sets.items():
            brute = bm_brute(ds_bm)
            ok = True
            for block_bytes in (1, 4 * bm_nlat * bm_nlon * 10, 10 ** 9):
                got = cell_validity_stats(ds_bm['v'], 'time', 'lat', 'lon', target_bytes=block_bytes)
                ok = ok and got[0] == brute[0] and all(np.array_equal(a, b) for a, b in zip(got[1:], brute[1:]))
            check(f"one-pass statistics == brute force [{bm_label}]", ok, True, tolerance=0)
            cache = CornerSeriesCache(ds_bm, bm_points[:4], ['v'], 'lat', 'lon', 'time')
            tiny = CornerSeriesCache(ds_bm, bm_points[:4], ['v'], 'lat', 'lon', 'time', max_bytes=bm_nt * 4)
            ok = True
            for la, lo in bm_points:
                ref = bm_old_corners(ds_bm, la, lo)
                for corners in (_build_bilinear_stencil(la, lo, ds_bm, 'v')['corners'],
                                _build_bilinear_stencil(la, lo, ds_bm, 'v', cache=cache)['corners'],
                                _build_bilinear_stencil(la, lo, ds_bm, 'v', cache=tiny)['corners']):
                    for k in ref:
                        ok = ok and corners[k].dtype == ref[k].dtype and np.array_equal(corners[k], ref[k], equal_nan=True)
            check(f"block read / cache / evicting cache == per-corner .sel reads [{bm_label}]",
                  ok and tiny.passes > cache.passes, True, tolerance=0)

        bm_ds = bm_sets['numpy (time, lat, lon)']
        check("no-time field counts as one step", cell_validity_stats(bm_ds['v'].isel(time=0), 'time', 'lat', 'lon')[0], 1, tolerance=0)
        check("extra dimension -> None (caller falls back)", cell_validity_stats(bm_ds['v'].expand_dims(layer=2), 'time', 'lat', 'lon') is None, True, tolerance=0)
        bm_empty = bm_ds.isel(time=slice(0, 0))
        check("zero-timestep window: cache series empty like the direct read",
              CornerSeriesCache(bm_empty, bm_points[:2], ['v'], 'lat', 'lon', 'time').get('v', 0, 0).shape == (0,)
              and _build_bilinear_stencil(0.12, 0.33, bm_empty, 'v')['corners']['Q11'].shape == (0,), True, tolerance=0)
        check("unknown variable is not served by the cache",
              CornerSeriesCache(bm_ds, bm_points[:2], ['missing'], 'lat', 'lon', 'time').has('missing'), False, tolerance=0)
        bm_blocks = list(iter_time_blocks(bm_ds['v'], 'time', target_bytes=4 * bm_nlat * bm_nlon * 10))
        check("numpy time blocks of the target size cover the axis once",
              [b[0] for b in bm_blocks] == [0, 10, 20, 30] and bm_blocks[-1][1] == bm_nt, True, tolerance=0)
        if 'dask (time, lat, lon)' in bm_sets:
            bm_dask = bm_sets['dask (time, lat, lon)']['v']
            one = list(iter_time_blocks(bm_dask, 'time', target_bytes=4 * bm_nlat * bm_nlon * 10))     # target 10 < 2 chunks of 7
            two = list(iter_time_blocks(bm_dask, 'time', target_bytes=4 * bm_nlat * bm_nlon * 14))     # target 14 = 2 chunks
            # target 10: one 7-step chunk per block, the trailing 2-step chunk joins the last (9 <= 10);
            # target 14: two chunks per block.
            check("dask blocks follow chunk edges and never exceed the target",
                  [b[0] for b in one] == [0, 7, 14, 21, 28] and [b[0] for b in two] == [0, 14, 28]
                  and one[-1][1] == bm_nt and two[-1][1] == bm_nt
                  and max(b[1] - b[0] for b in one) <= 10 and max(b[1] - b[0] for b in two) <= 14, True, tolerance=0)

    except ImportError:
        print("  xarray not available - skipping bilinear tests")

    print(f"\n{'=' * 50}")
    print("ALL TESTS PASSED" if all_pass else "SOME TESTS FAILED")
    print(f"{'=' * 50}")
