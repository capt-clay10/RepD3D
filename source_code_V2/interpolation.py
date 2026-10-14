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


def _build_bilinear_stencil(
    lat: float,
    lon: float,
    ds,
    variable: str,
    lat_coord: str = 'lat',
    lon_coord: str = 'lon',
) -> Dict[str, Any]:
    """Build the 4-corner stencil used by bilinear_interpolation()."""
    geometry = inspect_bilinear_geometry(
        lat, lon, ds, lat_coord=lat_coord, lon_coord=lon_coord
    )
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
) -> Dict[str, Any]:
    """
    Inspect the raw bilinear interpolation stencil for a single point.

    Returns metadata about each source corner so callers can decide whether
    interpolation would rely on all-NaN cells or long edge extrapolation.
    """
    stencil = _build_bilinear_stencil(
        lat, lon, ds, variable, lat_coord=lat_coord, lon_coord=lon_coord
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
) -> np.ndarray:
    """
    Bilinear interpolation on a regular lat/lon grid.

    Finds the 4 surrounding grid cells and computes an area-weighted average.

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
        lat, lon, ds, variable, lat_coord=lat_coord, lon_coord=lon_coord
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

    except ImportError:
        print("  xarray not available - skipping bilinear tests")

    print(f"\n{'=' * 50}")
    print("ALL TESTS PASSED" if all_pass else "SOME TESTS FAILED")
    print(f"{'=' * 50}")
