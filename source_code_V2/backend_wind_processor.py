# -*- coding: utf-8 -*-
"""
backend_wind_processor.py — Wind/pressure field generation for Delft3D.

Supports three atmospheric reanalysis sources:
    1. COSMO REA6  — rotated grid, vector rotation required
    2. ICON-DREAM-EU — unstructured triangular, KDTree+IDW interpolation
    3. ERA5 — regular grid, longitude wrapping required

All three share the same output: .amu (x_wind), .amv (y_wind), .amp (air_pressure).

v1 source: cosmo_wind_file_generator.py (417 lines) — COSMO only.
v2 additions: ICON and ERA5 pipelines.

Usage:
    from backend_wind_processor import generate_wind_fields
"""

import hashlib
import logging
import os
import re
import shutil
import threading
import numpy as np
import scipy.io
from scipy.interpolate import griddata
from datetime import datetime
from pathlib import Path
from typing import Callable, List, Optional, Dict, Any, Tuple

import interpolation as interp
import coordinate_utils
import file_parsers
import presets as preset_module
import writer_meteo


class TaskCancelledError(Exception):
    """Raised when the backend observes a cancel_event set by the GUI."""
    pass


def _check_cancel(cancel_event: Optional[threading.Event]) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise TaskCancelledError("Cancelled by user")


# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

KELVIN_TO_CELSIUS = -273.15   # Add to K values to get °C
TEMPERATURE_NODATA = -999.000  # Matches example_temp_file.amt sentinel


def _cleanup_idx_files(data_path: str, logger: logging.Logger) -> None:
    """Remove .idx sidecar files created by cfgrib when reading GRIB data."""
    data_dir = Path(data_path)
    if data_dir.is_file():
        data_dir = data_dir.parent
    if not data_dir.is_dir():
        return
    removed = 0
    # Check the directory itself and immediate subdirectories (UV/, PS/ for COSMO)
    for idx_file in data_dir.glob('*.idx'):
        try:
            idx_file.unlink()
            removed += 1
        except OSError:
            pass
    for idx_file in data_dir.glob('*/*.idx'):
        try:
            idx_file.unlink()
            removed += 1
        except OSError:
            pass
    if removed:
        logger.info(f"Cleaned up {removed} .idx sidecar file(s)")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def generate_wind_fields(
    preset_name: str,
    wind_data_path: str,
    grid_mat_path: str,
    reference_time: str,
    output_dir: str,
    output_name: str,
    aux_file_path: Optional[str] = None,
    crs_type: str = 'Cartesian',
    utm_zone: int = 32,
    hemisphere: str = 'N',
    grid_dx: float = 6000.0,
    grid_dy: float = 6000.0,
    write_wnd: bool = True,
    validation_points_path: Optional[str] = None,
    validation_crs_type: str = 'Spherical',
    validation_points: Optional[List[Dict[str, Any]]] = None,
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Generate Delft3D wind/pressure field files from reanalysis data.

    Args:
        preset_name:     'COSMO REA6', 'ICON-DREAM-EU', 'ICON-DREAM-Global', or 'ERA5'.
        wind_data_path:  Path to wind data directory (COSMO/ICON) or file (ERA5).
        grid_mat_path:   Path to Delft3D equidistant grid .mat file.
        reference_time:  Reference datetime string 'YYYY-MM-DD HH:MM:SS'.
        output_dir:      Output directory for .amu/.amv/.amp files.
        output_name:     Base filename for output files.
        aux_file_path:   Auxiliary file (COSMO: latlon_rot.grb, ICON: grid.nc).
        crs_type:        'Cartesian' or 'Spherical'.
        utm_zone:        UTM zone (for Cartesian grids).
        grid_dx:         Grid spacing X (for file header).
        grid_dy:         Grid spacing Y (for file header).
        write_wnd:       Also write .xwind/.ywind files.
        validation_points_path: Optional path to validation points file
            (3-column: name, lat/northing, lon/easting).
        validation_crs_type: CRS of validation points ('Spherical' or 'Cartesian').
        validation_points: Optional inline validation point definitions
            [{'name': str, 'lat': float, 'lon': float}, ...].
        logger:          Logger instance.
        progress_callback: Progress callback (0.0–1.0).

    Returns:
        {'output_files': {'amu': path, 'amv': path, 'amp': path, ...},
         'n_timesteps': int,
         'validation_files': [path, ...] (if validation output requested)}
    """
    if logger is None:
        logger = logging.getLogger('repd3d.wind')

    preset = preset_module.get_wind_preset(preset_name)
    if preset.sampling_crs != 'geographic':
        raise NotImplementedError(
            f"Wind preset '{preset.name}' declares unsupported sampling CRS '{preset.sampling_crs}'."
        )
    logger.info(f"Wind field generation — source: {preset.name}, CRS: {crs_type}")

    REF = datetime.strptime(reference_time, '%Y-%m-%d %H:%M:%S')

    # Load the Delft3D equidistant grid
    logger.info(f"Loading target grid: {Path(grid_mat_path).name}")
    grid_ed = scipy.io.loadmat(grid_mat_path)
    X = grid_ed['data']['X'][0][0][:-1, :-1]
    Y = grid_ed['data']['Y'][0][0][:-1, :-1]

    # Replace zeros with NaN (v1 behaviour)
    X[X == 0] = np.nan
    Y[Y == 0] = np.nan

    n_rows, n_cols = X.shape
    # Delft3D meteo grid_unit accepts only 'm' or 'degree' (Delft3D-FLOW manual
    # §B.7.1); the abbreviation 'deg' is not a valid value.
    grid_unit = 'degree' if crs_type == 'Spherical' else 'm'

    grid_info = {
        'n_cols': n_cols,
        'n_rows': n_rows,
        # Values are sampled at the grid nodes, so the lower-left node is the
        # centre of the lower-left meteo cell -> written as x_llcenter/y_llcenter.
        # Writing it as x_llcorner shifted every field by dx/2, dy/2 (manual B.7.1).
        'x_llcenter': float(np.nanmin(X)),
        'y_llcenter': float(np.nanmin(Y)),
        'dx': grid_dx,
        'dy': grid_dy,
        'grid_unit': grid_unit,
    }

    logger.info(f"Target grid: {n_rows}×{n_cols}, unit={grid_unit}")
    progress_callback(0.05)

    query_grid = _prepare_geographic_query_grid(
        X, Y, crs_type, utm_zone, hemisphere, preset.name, logger,
    )

    # Dispatch to source-specific processor
    if preset.name == 'ERA5':
        # ERA5: streaming pipeline (chunked read → bilinear remap → write)
        n_ts, result, validation_files = _era5_streaming_pipeline_v2(
            wind_data_path, query_grid, preset, REF,
            grid_info, output_dir, output_name, write_wnd,
            X, Y, validation_points_path, validation_crs_type, validation_points,
            crs_type, utm_zone, hemisphere,
            logger, progress_callback, cancel_event,
        )
        progress_callback(1.0)
        _cleanup_idx_files(wind_data_path, logger)
        out = {'output_files': result, 'n_timesteps': n_ts}
        if validation_files:
            out['validation_files'] = validation_files
        return out

    # COSMO / ICON: existing batch approach
    if preset.name == 'COSMO REA6':
        time_series = _process_cosmo(
            wind_data_path, aux_file_path, query_grid,
            preset, REF, logger, progress_callback, cancel_event,
        )
    elif preset.name in ('ICON-DREAM-EU', 'ICON-DREAM-Global'):
        time_series = _process_icon(
            wind_data_path, aux_file_path, query_grid,
            preset, REF, logger, progress_callback, cancel_event,
        )
    else:
        raise ValueError(f"Unsupported wind preset: {preset.name}")

    _check_cancel(cancel_event)

    # Write output files
    logger.info(f"Writing {len(time_series)} timesteps to {output_dir}...")
    result = writer_meteo.write_meteo_files(
        output_dir, output_name, grid_info, time_series, REF,
        write_wnd=write_wnd, logger_ext=logger,
    )

    # Extract validation time series if requested
    validation_files = []
    if validation_points or validation_points_path:
        validation_files = _extract_validation_timeseries(
            time_series, X, Y, validation_points_path,
            validation_crs_type, crs_type, utm_zone, hemisphere,
            output_dir, logger, validation_points=validation_points,
        )

    progress_callback(1.0)
    _cleanup_idx_files(wind_data_path, logger)
    out = {
        'output_files': result,
        'n_timesteps': len(time_series),
    }
    if validation_files:
        out['validation_files'] = validation_files
    return out


# ─────────────────────────────────────────────────────────────────────────────
# TEMPERATURE (.amt) ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def generate_temperature_field(
    preset_name: str,
    temperature_data_path: str,
    grid_mat_path: str,
    reference_time: str,
    output_dir: str,
    output_name: str,
    aux_file_path: Optional[str] = None,
    crs_type: str = 'Cartesian',
    utm_zone: int = 32,
    hemisphere: str = 'N',
    grid_dx: float = 6000.0,
    grid_dy: float = 6000.0,
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
    **_unused,   # Swallow any extra kwargs the GUI sends (e.g. write_wnd)
) -> dict:
    """
    Generate a Delft3D .amt air-temperature field from reanalysis data.

    Supports: COSMO REA6, ICON-DREAM-EU, ICON-DREAM-Global, ERA5.

    Args mirror `generate_wind_fields` but only temperature is extracted.
    The temperature is converted from Kelvin to Celsius before writing.

    For ICON-DREAM-Global the temperature variable is multi-level
    (generalVerticalLayer); the lowest (near-surface) level is selected.
    """
    if logger is None:
        logger = logging.getLogger('repd3d.meteo')

    preset = preset_module.get_wind_preset(preset_name)
    if not preset.t_var:
        raise ValueError(
            f"Preset '{preset.name}' has no temperature variable (t_var) "
            f"configured. Temperature extraction is not supported for this source."
        )
    if preset.sampling_crs != 'geographic':
        raise NotImplementedError(
            f"Temperature preset '{preset.name}' declares unsupported sampling CRS '{preset.sampling_crs}'."
        )

    logger.info(
        f"Temperature field generation — source: {preset.name}, "
        f"variable: {preset.t_var}, CRS: {crs_type}"
    )

    REF = datetime.strptime(reference_time, '%Y-%m-%d %H:%M:%S')

    # Load the Delft3D equidistant grid
    logger.info(f"Loading target grid: {Path(grid_mat_path).name}")
    grid_ed = scipy.io.loadmat(grid_mat_path)
    X = grid_ed['data']['X'][0][0][:-1, :-1]
    Y = grid_ed['data']['Y'][0][0][:-1, :-1]
    X[X == 0] = np.nan
    Y[Y == 0] = np.nan

    n_rows, n_cols = X.shape
    grid_unit = 'degree' if crs_type == 'Spherical' else 'm'

    grid_info = {
        'n_cols': n_cols,
        'n_rows': n_rows,
        'x_llcenter': float(np.nanmin(X)),
        'y_llcenter': float(np.nanmin(Y)),
        'dx': grid_dx,
        'dy': grid_dy,
        'grid_unit': grid_unit,
    }
    logger.info(f"Target grid: {n_rows}×{n_cols}, unit={grid_unit}")
    progress_callback(0.05)

    query_grid = _prepare_geographic_query_grid(
        X, Y, crs_type, utm_zone, hemisphere, preset.name, logger,
    )

    output_path = str(Path(output_dir) / f'{output_name}.amt')
    if preset.name == 'ERA5':
        n_ts, final_path = _era5_temperature_streaming_pipeline(
            temperature_data_path, query_grid, preset, REF,
            grid_info, output_path, logger, progress_callback, cancel_event,
        )
        progress_callback(1.0)
        _cleanup_idx_files(temperature_data_path, logger)
        return {
            'output_files': {'amt': final_path},
            'n_timesteps': n_ts,
        }

    # Dispatch to source-specific temperature processor
    if preset.name == 'COSMO REA6':
        time_series = _process_cosmo_temperature(
            temperature_data_path, aux_file_path, query_grid,
            preset, REF, logger, progress_callback, cancel_event,
        )
    elif preset.name in ('ICON-DREAM-EU', 'ICON-DREAM-Global'):
        time_series = _process_icon_temperature(
            temperature_data_path, aux_file_path, query_grid,
            preset, REF, logger, progress_callback, cancel_event,
        )
    else:
        raise ValueError(f"Unsupported preset for temperature: {preset.name}")

    _check_cancel(cancel_event)

    # Write the .amt file
    logger.info(f"Writing {len(time_series)} timesteps to {output_path}...")
    writer_meteo.write_amt_file(
        output_path, grid_info, time_series, REF,
        nodata_value=TEMPERATURE_NODATA, logger_ext=logger,
    )

    progress_callback(1.0)
    _cleanup_idx_files(temperature_data_path, logger)
    return {
        'output_files': {'amt': output_path},
        'n_timesteps': len(time_series),
    }


def _celsius_and_fill(T: np.ndarray) -> np.ndarray:
    """Convert Kelvin → Celsius and replace NaN with TEMPERATURE_NODATA."""
    T_c = T + KELVIN_TO_CELSIUS
    T_c[np.isnan(T_c)] = TEMPERATURE_NODATA
    return T_c


# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION POINT EXTRACTION
# ─────────────────────────────────────────────────────────────────────────────

def _finite_bounds(values: np.ndarray, label: str) -> tuple:
    """Return (min, max) for finite values or raise if none exist."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError(f"{label} has no finite values.")
    return float(np.min(finite)), float(np.max(finite))


def _normalize_longitudes(lons: np.ndarray) -> np.ndarray:
    """Wrap longitudes to the [-180, 180) range."""
    out = np.asarray(lons, dtype=float).copy()
    finite = np.isfinite(out)
    out[finite] = ((out[finite] + 180.0) % 360.0) - 180.0
    return out


def _prepare_geographic_query_grid(
    X: np.ndarray,
    Y: np.ndarray,
    crs_type: str,
    utm_zone: int,
    hemisphere: str,
    source_name: str,
    logger: logging.Logger,
) -> Dict[str, np.ndarray]:
    """
    Build the geographic sampling grid used for all source lookups.

    The target Delft3D grid remains in its native CRS for output. Source
    processors should sample exclusively on the returned geographic query grid.
    """
    finite_mask = np.isfinite(X) & np.isfinite(Y)
    if not np.any(finite_mask):
        raise ValueError("Target grid has no finite cells to sample.")

    query_lats = np.full(X.shape, np.nan, dtype=float)
    query_lons = np.full(X.shape, np.nan, dtype=float)

    if str(crs_type).lower().startswith('cart'):
        lat, lon = coordinate_utils.utm_to_latlon(
            X[finite_mask], Y[finite_mask], utm_zone, hemisphere,
        )
        query_lats[finite_mask] = np.asarray(lat, dtype=float)
        query_lons[finite_mask] = _normalize_longitudes(np.asarray(lon, dtype=float))
    else:
        query_lats[finite_mask] = np.asarray(Y[finite_mask], dtype=float)
        query_lons[finite_mask] = _normalize_longitudes(
            np.asarray(X[finite_mask], dtype=float)
        )

    lat_min, lat_max = _finite_bounds(query_lats, 'Query latitude grid')
    lon_min, lon_max = _finite_bounds(query_lons, 'Query longitude grid')
    if lat_min < -90.0 or lat_max > 90.0:
        raise ValueError(
            f"{source_name}: geographic query latitude bounds [{lat_min:.4f}, {lat_max:.4f}] "
            f"are invalid. This usually means the target CRS '{crs_type}' or UTM settings "
            f"do not match the grid coordinates."
        )
    if lon_min < -360.0 or lon_max > 360.0:
        raise ValueError(
            f"{source_name}: geographic query longitude bounds [{lon_min:.4f}, {lon_max:.4f}] "
            f"are implausible. This usually means Cartesian metre coordinates were passed as "
            f"lat/lon or the target CRS settings are incorrect."
        )

    x_min, x_max = _finite_bounds(X, 'Target X grid')
    y_min, y_max = _finite_bounds(Y, 'Target Y grid')
    logger.info(
        f"Sampling setup for {source_name}: target CRS={crs_type}, sampling CRS=geographic, "
        f"target bounds X=[{x_min:.3f}, {x_max:.3f}] Y=[{y_min:.3f}, {y_max:.3f}], "
        f"query lat=[{lat_min:.4f}, {lat_max:.4f}] lon=[{lon_min:.4f}, {lon_max:.4f}]"
    )

    return {
        'target_x': X,
        'target_y': Y,
        'query_lats': query_lats,
        'query_lons': query_lons,
        'finite_mask': finite_mask,
    }


def _query_bounds(query_grid: Dict[str, np.ndarray]) -> tuple:
    """Return geographic query bounds (lat_min, lat_max, lon_min, lon_max)."""
    lat_min, lat_max = _finite_bounds(query_grid['query_lats'], 'Query latitude grid')
    lon_min, lon_max = _finite_bounds(query_grid['query_lons'], 'Query longitude grid')
    return lat_min, lat_max, lon_min, lon_max


def _find_structured_overlap_slices(
    source_lats: np.ndarray,
    source_lons: np.ndarray,
    query_grid: Dict[str, np.ndarray],
    source_name: str,
    buffer_cells: int = 2,
    buffer_degrees: float = 0.5,
) -> tuple:
    """Find 2D source slices covering the geographic query grid."""
    lat_min, lat_max, lon_min, lon_max = _query_bounds(query_grid)
    mask = (
        (source_lats >= lat_min - buffer_degrees)
        & (source_lats <= lat_max + buffer_degrees)
        & (source_lons >= lon_min - buffer_degrees)
        & (source_lons <= lon_max + buffer_degrees)
    )
    valid = np.where(mask)
    if valid[0].size == 0:
        raise ValueError(
            f"{source_name}: target grid query bounds lat=[{lat_min:.4f}, {lat_max:.4f}] "
            f"lon=[{lon_min:.4f}, {lon_max:.4f}] do not overlap the source domain."
        )

    mm1 = max(0, int(np.min(valid[0])) - buffer_cells)
    mm2 = min(source_lats.shape[0], int(np.max(valid[0])) + buffer_cells + 1)
    nn1 = max(0, int(np.min(valid[1])) - buffer_cells)
    nn2 = min(source_lats.shape[1], int(np.max(valid[1])) + buffer_cells + 1)
    return mm1, mm2, nn1, nn2


def _subset_regular_grid(
    source_lats: np.ndarray,
    source_lons: np.ndarray,
    query_grid: Dict[str, np.ndarray],
    source_name: str,
    buffer_degrees: float = 1.0,
) -> tuple:
    """Return boolean lat/lon masks covering the geographic query grid."""
    lat_min, lat_max, lon_min, lon_max = _query_bounds(query_grid)
    lon_mask = (
        (source_lons >= lon_min - buffer_degrees)
        & (source_lons <= lon_max + buffer_degrees)
    )
    lat_mask = (
        (source_lats >= lat_min - buffer_degrees)
        & (source_lats <= lat_max + buffer_degrees)
    )

    if not np.any(lon_mask) or not np.any(lat_mask):
        raise ValueError(
            f"{source_name}: target grid query bounds lat=[{lat_min:.4f}, {lat_max:.4f}] "
            f"lon=[{lon_min:.4f}, {lon_max:.4f}] do not overlap the source domain."
        )
    return lat_mask, lon_mask


def _griddata_geographic(
    source_lats: np.ndarray,
    source_lons: np.ndarray,
    source_values: np.ndarray,
    query_grid: Dict[str, np.ndarray],
    method: str = 'linear',
) -> np.ndarray:
    """Interpolate a geographic field onto the geographic query grid."""
    src_lats = np.asarray(source_lats, dtype=float).ravel()
    src_lons = np.asarray(source_lons, dtype=float).ravel()
    src_vals = np.asarray(source_values, dtype=float).ravel()
    valid = np.isfinite(src_lats) & np.isfinite(src_lons) & np.isfinite(src_vals)
    if not np.any(valid):
        return np.full(query_grid['query_lats'].shape, np.nan, dtype=float)
    return griddata(
        (src_lons[valid], src_lats[valid]),
        src_vals[valid],
        (query_grid['query_lons'], query_grid['query_lats']),
        method=method,
    )


# ─────────────────────────────────────────────────────────────────────────────
# ERA5 STRUCTURED BILINEAR REMAP
# ─────────────────────────────────────────────────────────────────────────────

def _build_era5_bilinear_weights(
    src_lats_1d: np.ndarray,
    src_lons_1d: np.ndarray,
    query_grid: Dict[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    """
    Precompute bilinear interpolation weights from an ERA5 regular lat/lon
    grid onto the geographic query grid.

    Exploits the regular spacing to compute bracket indices analytically
    rather than searching.  Called once per run; the returned dict is
    applied to every timestep via ``_apply_bilinear_remap``.

    Args:
        src_lats_1d: 1-D latitude array of the spatial subset (may be
                     ascending **or** descending — ERA5 is descending).
        src_lons_1d: 1-D longitude array of the spatial subset (ascending).
        query_grid:  Standard query-grid dict with ``query_lats``,
                     ``query_lons``, ``finite_mask``.

    Returns:
        dict with index arrays (j0, j1, i0, i1), weight arrays
        (w00, w01, w10, w11), and a boolean ``valid`` mask, all matching
        the shape of the query grid.
    """
    qlats = query_grid['query_lats']
    qlons = query_grid['query_lons']
    finite = query_grid['finite_mask']

    n_lat = len(src_lats_1d)
    n_lon = len(src_lons_1d)

    # Regular grid spacing (dlat is negative when lats descend)
    dlat = float(src_lats_1d[1] - src_lats_1d[0])
    dlon = float(src_lons_1d[1] - src_lons_1d[0])

    # Continuous fractional indices into the subset arrays
    j_float = (qlats - src_lats_1d[0]) / dlat
    i_float = (qlons - src_lons_1d[0]) / dlon

    j0 = np.floor(j_float).astype(int)
    i0 = np.floor(i_float).astype(int)
    j1 = j0 + 1
    i1 = i0 + 1

    fj = j_float - j0
    fi = i_float - i0

    valid = (
        finite
        & np.isfinite(qlats)
        & np.isfinite(qlons)
        & (j0 >= 0) & (j1 < n_lat)
        & (i0 >= 0) & (i1 < n_lon)
    )

    # Clamp so fancy-indexing never raises; invalid cells get NODATA anyway
    j0c = np.clip(j0, 0, n_lat - 1)
    j1c = np.clip(j1, 0, n_lat - 1)
    i0c = np.clip(i0, 0, n_lon - 1)
    i1c = np.clip(i1, 0, n_lon - 1)

    return {
        'j0': j0c, 'j1': j1c, 'i0': i0c, 'i1': i1c,
        'w00': (1.0 - fj) * (1.0 - fi),
        'w01': (1.0 - fj) * fi,
        'w10': fj * (1.0 - fi),
        'w11': fj * fi,
        'valid': valid,
        'shape': qlats.shape,
    }


def _apply_bilinear_remap(
    field_2d: np.ndarray,
    weights: Dict[str, np.ndarray],
    nodata: float = 9999.00,
) -> np.ndarray:
    """Apply precomputed bilinear weights to a 2-D source field."""
    w = weights
    result = (
        w['w00'] * field_2d[w['j0'], w['i0']]
        + w['w01'] * field_2d[w['j0'], w['i1']]
        + w['w10'] * field_2d[w['j1'], w['i0']]
        + w['w11'] * field_2d[w['j1'], w['i1']]
    )
    result[~w['valid']] = nodata
    result[~np.isfinite(result)] = nodata
    return result


def _kdtree_interpolate_query_grid(
    tree,
    source_values: np.ndarray,
    query_grid: Dict[str, np.ndarray],
    k: int = 6,
) -> np.ndarray:
    """Interpolate onto the geographic query grid while skipping invalid target cells."""
    result = np.full(query_grid['query_lats'].shape, np.nan, dtype=float)
    target_mask = (
        query_grid['finite_mask']
        & np.isfinite(query_grid['query_lats'])
        & np.isfinite(query_grid['query_lons'])
    )
    if not np.any(target_mask):
        return result

    source_values = np.asarray(source_values, dtype=float)
    if source_values.size == 0 or np.all(np.isnan(source_values)):
        return result

    sampled = interp.kdtree_idw_interpolation(
        tree,
        source_values,
        query_grid['query_lats'][target_mask],
        query_grid['query_lons'][target_mask],
        k=k,
    )
    result[target_mask] = sampled
    return result


def _parse_icon_nominal_year_month(path: str) -> Optional[tuple]:
    """Extract the nominal (year, month) from an ICON monthly filename."""
    match = re.search(r'_(\d{6})(?:_|\\.)', Path(path).name)
    if not match:
        return None
    token = match.group(1)
    return int(token[:4]), int(token[4:6])


def _collect_icon_time_records(ds, file_path: str, logger: logging.Logger) -> List[tuple]:
    """Flatten ICON time x step coordinates and trim records outside the nominal file month."""
    nominal = _parse_icon_nominal_year_month(file_path)
    records: List[tuple] = []
    trimmed = 0

    for t_idx in range(len(ds.time)):
        base_time = ds.time.values[t_idx].astype('datetime64[s]')
        for s_idx in range(len(ds.step)):
            valid_time = (base_time + ds.step.values[s_idx]).astype('datetime64[s]').tolist()
            if nominal is not None and (
                valid_time.year != nominal[0] or valid_time.month != nominal[1]
            ):
                trimmed += 1
                continue
            records.append((t_idx, s_idx, valid_time))

    if nominal is not None and trimmed:
        logger.info(
            f"{Path(file_path).name}: trimmed {trimmed} out-of-month ICON timestep(s) "
            f"outside {nominal[0]}-{nominal[1]:02d}"
        )
    return records


def _sort_and_deduplicate_time_series(
    time_series: List[tuple],
    logger: logging.Logger,
    label: str,
) -> List[tuple]:
    """Ensure timestamps are strictly increasing by dropping exact duplicates."""
    sorted_series = sorted(time_series, key=lambda item: item[0])
    deduped: List[tuple] = []
    duplicate_count = 0

    for item in sorted_series:
        if deduped and item[0] == deduped[-1][0]:
            duplicate_count += 1
            continue
        deduped.append(item)

    if duplicate_count:
        logger.warning(f"{label}: removed {duplicate_count} duplicate timestep(s) after concatenation.")
    return deduped


def _warn_if_globally_empty(
    logger: logging.Logger,
    label: str,
    timestamp: datetime,
    arrays: List[np.ndarray],
) -> None:
    """Log when a kept source timestep is entirely NaN."""
    if arrays and all(np.all(np.isnan(np.asarray(arr, dtype=float))) for arr in arrays):
        logger.warning(
            f"{label}: source timestep {timestamp:%Y-%m-%d %H:%M:%S} is entirely NaN; "
            f"keeping it as NODATA in the generated fields."
        )


def _find_nearest_grid_cell(
    target_x: float, target_y: float,
    X: np.ndarray, Y: np.ndarray,
) -> tuple:
    """
    Find the (row, col) of the nearest non-NaN grid cell to (target_x, target_y).

    Returns (row, col) index or None if no valid cell within reasonable distance.
    """
    valid = ~np.isnan(X) & ~np.isnan(Y)
    if not np.any(valid):
        return None

    rows, cols = np.where(valid)
    dist = np.sqrt(
        (X[valid] - target_x) ** 2 + (Y[valid] - target_y) ** 2
    )
    idx = np.argmin(dist)
    return int(rows[idx]), int(cols[idx])


def _resolve_validation_points(
    validation_points_path: Optional[str],
    validation_crs_type: str,
    utm_zone: int,
    hemisphere: str,
    validation_points: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Resolve inline or file-based validation point definitions."""
    if validation_points is not None:
        points: List[Dict[str, Any]] = []
        for idx, point in enumerate(validation_points, start=1):
            name = str(point.get('name', '')).strip()
            if not name:
                raise ValueError(f"Validation point {idx} is missing a name.")
            try:
                lat = float(point['lat'])
                lon = float(point['lon'])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    f"Validation point '{name}' must provide numeric lat/lon values."
                ) from exc
            if not -90.0 <= lat <= 90.0:
                raise ValueError(
                    f"Validation point '{name}' has latitude {lat}, but latitude must be between -90 and 90."
                )
            if not -180.0 <= lon <= 180.0:
                raise ValueError(
                    f"Validation point '{name}' has longitude {lon}, but longitude must be between -180 and 180."
                )
            points.append({
                'name': name,
                'lat': lat,
                'lon': lon,
                'easting': None,
                'northing': None,
            })
        return points

    if validation_points_path:
        return file_parsers.parse_validation_points(
            validation_points_path, validation_crs_type, utm_zone, hemisphere,
        )

    return []


def _extract_validation_timeseries(
    time_series: List,
    X: np.ndarray,
    Y: np.ndarray,
    validation_points_path: Optional[str],
    validation_crs_type: str,
    grid_crs_type: str,
    utm_zone: int,
    hemisphere: str,
    output_dir: str,
    logger: logging.Logger,
    validation_points: Optional[List[Dict[str, Any]]] = None,
) -> List[str]:
    """
    Extract single-point wind/pressure time series at validation points.

    For each point: finds the nearest grid cell, extracts u/v/p over all
    timesteps, computes wind speed and meteorological "from" direction,
    and writes one CSV to output_dir/validation/<name>_validation.csv.

    Returns list of written CSV paths.
    """
    NODATA = 9999.00

    points = _resolve_validation_points(
        validation_points_path, validation_crs_type, utm_zone, hemisphere,
        validation_points=validation_points,
    )
    if not points:
        return []
    logger.info(f"Extracting validation time series for {len(points)} point(s)")

    val_dir = Path(output_dir) / 'validation'
    val_dir.mkdir(parents=True, exist_ok=True)

    written = []

    for pt in points:
        # Convert validation lat/lon to grid coordinates
        if grid_crs_type.lower().startswith('cart'):
            pt_x, pt_y = coordinate_utils.latlon_to_utm(
                pt['lat'], pt['lon'], utm_zone, hemisphere,
            )
        else:
            pt_x = pt['lon']
            pt_y = pt['lat']

        cell = _find_nearest_grid_cell(float(pt_x), float(pt_y), X, Y)
        if cell is None:
            logger.warning(f"  {pt['name']}: no valid grid cell found — skipped")
            continue

        r, c = cell
        grid_dist = np.sqrt(
            (X[r, c] - float(pt_x)) ** 2 + (Y[r, c] - float(pt_y)) ** 2
        )
        logger.info(f"  {pt['name']}: nearest cell ({r},{c}), distance={grid_dist:.1f}")

        datetimes = []
        speeds = []
        dirs_from = []
        pressures = []

        for T, U, V, P in time_series:
            u_val = U[r, c]
            v_val = V[r, c]
            p_val = P[r, c]

            if u_val == NODATA or v_val == NODATA:
                speed = NODATA
                d_from = NODATA
            else:
                speed = np.sqrt(u_val ** 2 + v_val ** 2)
                d_from = (270.0 - np.degrees(np.arctan2(v_val, u_val))) % 360.0

            if p_val == NODATA:
                p_val = NODATA

            datetimes.append(T)
            speeds.append(float(speed))
            dirs_from.append(float(d_from))
            pressures.append(float(p_val))

        csv_path = val_dir / f"{pt['name']}_validation.csv"
        writer_meteo.write_validation_csv(
            str(csv_path), datetimes, speeds, dirs_from, pressures,
        )
        written.append(str(csv_path))
        logger.info(f"  Written: {csv_path.name} ({len(datetimes)} timesteps)")

    return written


# ─────────────────────────────────────────────────────────────────────────────
# COSMO REA6 PROCESSOR
# ─────────────────────────────────────────────────────────────────────────────

def calculate_theta(lat, lon, pollat=39.25, pollon=-162.0):
    """
    Calculate the angle that turns COSMO grid-relative winds into geographic winds.

    Rotated meridians are great circles through the rotated north pole, so grid
    north at a point is the bearing from that point to the pole. This is
    ``delta`` of DWD's ``uvrot2uv`` (COSMO utilities):

        zarg1 = cos(pollat) * sin(pollon - lon)
        zarg2 = sin(pollat) * cos(lat) - cos(pollat) * sin(lat) * cos(pollon - lon)
        delta = atan2(zarg1, zarg2)       (clockwise from true north)

    ``rotate_wind_vectors`` turns a vector counter-clockwise by theta, and
    grid-relative -> geographic is a clockwise turn by delta, so theta = -delta.

    Verified 2026-10-14 against the bearing of the grid axes measured from the
    exact grid coordinates in cosmo_latlon_rot.grb: max difference 0.0025 deg
    over the whole REA6 domain; theta = 7.84 deg in the German Bight.

    v1 (cosmo_wind_file_generator.py lines 27-53) used
    atan2(sin(lon - pollon), cos(pollat)*tan(lat) - sin(pollat)*cos(lon - pollon)),
    which is the bearing at the pole toward the point instead of at the point
    toward the pole: 1.9 deg short in the German Bight, 3.3 deg in the central
    North Sea, 13.6 deg near Iceland, 20 deg at the NW edge. Do not reintroduce it.

    IMPORTANT: The COSMO GRIB flag (GRIB_uvRelativeToGrid=0) says
    winds are geographic — BUT THIS IS WRONG. DWD documentation says
    winds are on the rotated grid. Rotation is ALWAYS applied.

    Args:
        lat: Geographic latitude (array, degrees). The cfgrib 'latitude'
             coordinate of the aux file is geographic, not rotated.
        lon: Geographic longitude (array, degrees).
        pollat: Rotated north pole latitude (COSMO-REA6: 39.25°).
        pollon: Rotated north pole longitude (COSMO-REA6: -162.0°).

    Returns:
        theta: Rotation angle (array, radians), same shape as lat/lon.
    """
    lat_rad = np.radians(lat)
    dlon_rad = np.radians(pollon - lon)
    pollat_rad = np.radians(pollat)

    zarg1 = np.cos(pollat_rad) * np.sin(dlon_rad)
    zarg2 = (
        np.sin(pollat_rad) * np.cos(lat_rad)
        - np.cos(pollat_rad) * np.sin(lat_rad) * np.cos(dlon_rad)
    )
    return -np.arctan2(zarg1, zarg2)


def rotate_wind_vectors(u, v, theta):
    """
    Rotate U/V wind components from COSMO rotated grid to true geographic.

    v1 source: cosmo_wind_file_generator.py lines 113-114.

    U_corrected = U * cos(θ) - V * sin(θ)
    V_corrected = U * sin(θ) + V * cos(θ)
    """
    return (
        u * np.cos(theta) - v * np.sin(theta),
        u * np.sin(theta) + v * np.cos(theta),
    )


def _process_cosmo(
    wind_dir, aux_path, query_grid, preset, REF, logger, progress_callback, cancel_event=None
):
    """
    Process COSMO REA6 data files.
    v1 source: cosmo_wind_file_generator.py lines 169-416.
    """
    import cfgrib
    import xarray as xr

    logger.info("Loading COSMO auxiliary file for coordinate rotation...")

    # Load rotated and true coordinates from aux file
    ds = xr.open_dataset(
        aux_path, engine="cfgrib",
        filter_by_keys={'typeOfLevel': 'surface'},
    )
    lat_true = ds['RLAT'].values
    lon_true = ds['RLON'].values
    # cfgrib's latitude/longitude are geographic (the same positions as
    # RLAT/RLON, computed exactly from the grid definition), not rotated.
    grid_lat = ds['latitude'].values
    grid_lon = ds['longitude'].values

    theta = calculate_theta(grid_lat, grid_lon)

    mm1, mm2, nn1, nn2 = _find_structured_overlap_slices(
        lat_true, _normalize_longitudes(lon_true), query_grid, 'COSMO REA6',
    )

    theta = theta[mm1:mm2, nn1:nn2]
    lat_sub = lat_true[mm1:mm2, nn1:nn2]
    lon_sub = _normalize_longitudes(lon_true[mm1:mm2, nn1:nn2])

    logger.info(f"COSMO overlap region: [{mm1}:{mm2}, {nn1}:{nn2}]")

    # Find and sort COSMO data files — support both UV/PS subdirectory
    # layout (DWD distribution) and a flat directory with all files together
    def extract_date(filename):
        return int(filename.split('.')[-2])

    uv_dir = os.path.join(wind_dir, 'UV')
    ps_dir = os.path.join(wind_dir, 'PS')

    if os.path.isdir(uv_dir) and os.path.isdir(ps_dir):
        u_search_dir = uv_dir
        v_search_dir = uv_dir
        p_search_dir = ps_dir
        logger.info("Using UV/ and PS/ subdirectory layout")
    else:
        u_search_dir = wind_dir
        v_search_dir = wind_dir
        p_search_dir = wind_dir
        logger.info("Using flat directory layout (all files in one folder)")

    u_files = sorted(
        [f for f in os.listdir(u_search_dir) if f.startswith('U_') and f.endswith('.grb')],
        key=extract_date,
    )
    v_files = sorted(
        [f for f in os.listdir(v_search_dir) if f.startswith('V_') and f.endswith('.grb')],
        key=extract_date,
    )
    p_files = sorted(
        [f for f in os.listdir(p_search_dir) if f.startswith('PS') and f.endswith('.grb')],
        key=extract_date,
    )

    if not u_files or not v_files or not p_files:
        raise FileNotFoundError(
            f"COSMO REA6: expected U_*.grb, V_*.grb, and PS*.grb files in "
            f"{wind_dir} (or UV/ and PS/ subdirectories), but found "
            f"{len(u_files)} U, {len(v_files)} V, and {len(p_files)} P files."
        )

    logger.info(f"Found {len(u_files)} U files, {len(v_files)} V files, {len(p_files)} P files")

    # Process each monthly file
    time_series = []
    total_files = len(u_files)

    for fi in range(total_files):
        _check_cancel(cancel_event)
        WU = cfgrib.open_dataset(os.path.join(u_search_dir, u_files[fi]))
        WV = cfgrib.open_dataset(os.path.join(v_search_dir, v_files[fi]))
        PS = cfgrib.open_dataset(os.path.join(p_search_dir, p_files[fi]))

        for kk in range(len(WU.time)):
            _check_cancel(cancel_event)
            DATA_U = WU[preset.u_var].values[kk]
            DATA_V = WV[preset.v_var].values[kk]
            DATA_P = PS[preset.p_var].values[kk]

            # Slice to overlap region
            wnd_u = DATA_U[mm1:mm2, nn1:nn2]
            wnd_v = DATA_V[mm1:mm2, nn1:nn2]

            # Rotate wind vectors
            U_rot, V_rot = rotate_wind_vectors(wnd_u, wnd_v, theta)

            # Interpolate onto the geographic sampling grid
            U_interp = _griddata_geographic(lat_sub, lon_sub, U_rot, query_grid)
            V_interp = _griddata_geographic(lat_sub, lon_sub, V_rot, query_grid)
            P_interp = _griddata_geographic(
                lat_sub, lon_sub, DATA_P[mm1:mm2, nn1:nn2], query_grid,
            ) / 100.0

            # Replace NaN with NODATA
            U_interp[np.isnan(U_interp)] = 9999.00
            V_interp[np.isnan(V_interp)] = 9999.00
            P_interp[np.isnan(P_interp)] = 9999.00

            T = WU.time.values[kk].astype('datetime64[s]').tolist()
            time_series.append((T, U_interp, V_interp, P_interp))

        progress_callback(0.05 + 0.9 * (fi + 1) / total_files)
        logger.info(f"  Processed file {fi + 1}/{total_files}: {u_files[fi]}, {v_files[fi]}, {p_files[fi]}")

    return _sort_and_deduplicate_time_series(time_series, logger, f"{preset.name} wind")


_ERA5_CHUNK_SIZE = 168   # timesteps per GRIB read-chunk (1 week hourly)
# Cache lives in a user-writable location so it works inside a frozen .exe
# (where __file__ resolves into the read-only PyInstaller _MEIPASS temp dir).
_ERA5_NETCDF_CACHE_DIR = (
    Path(os.getenv('LOCALAPPDATA') or Path.home() / '.cache')
    / 'RepD3D' / 'era5_netcdf'
)
_ERA5_GRIB_SHORTNAMES = {
    'u10': '10u',
    'v10': '10v',
    'sp': 'sp',
    't2m': '2t',
}


def _collect_era5_input_paths(data_path: str, label: str) -> List[str]:
    """Resolve ERA5 NetCDF input as either a single file or a directory of files."""
    if os.path.isdir(data_path):
        candidates = [
            name for name in os.listdir(data_path)
            if name.lower().endswith('.nc')
        ]
        preferred = [name for name in candidates if 'era5' in name.lower()]
        selected = preferred or candidates
        paths = sorted(str(Path(data_path) / name) for name in selected)
        if not paths:
            raise FileNotFoundError(
                f"No ERA5 NetCDF .nc files found in {data_path}. "
                "Use pre-converted/downloaded ERA5 NetCDF files; GRIB input requires ecCodes."
            )
        return paths
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"{label}: file not found: {data_path}")
    if not str(data_path).lower().endswith('.nc'):
        raise ValueError(
            f"{label}: ERA5 input must be a NetCDF .nc file. "
            "Use pre-converted/downloaded ERA5 NetCDF files; GRIB input requires ecCodes."
        )
    return [data_path]


def _build_era5_cache_key(source_path: str) -> str:
    """Build a stable cache key from path + size + mtime."""
    resolved = Path(source_path).resolve()
    stat = resolved.stat()
    payload = f"{resolved}|{stat.st_size}|{stat.st_mtime_ns}"
    return hashlib.sha1(payload.encode('utf-8')).hexdigest()[:16]


def _ensure_era5_cache_dir() -> Path:
    """Return the ERA5 NetCDF cache directory, creating it on demand."""
    _ERA5_NETCDF_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return _ERA5_NETCDF_CACHE_DIR


def _resolve_grib_to_netcdf_binary() -> str:
    """Locate the ecCodes GRIB->NetCDF converter."""
    exe = shutil.which('grib_to_netcdf.exe') or shutil.which('grib_to_netcdf')
    if exe:
        return exe
    raise RuntimeError(
        "ERA5 NetCDF conversion requires 'grib_to_netcdf.exe' (ecCodes) on PATH."
    )


def _resolve_grib_copy_binary() -> str:
    """Locate the ecCodes GRIB split/copy tool."""
    exe = shutil.which('grib_copy.exe') or shutil.which('grib_copy')
    if exe:
        return exe
    raise RuntimeError(
        "ERA5 NetCDF conversion requires 'grib_copy.exe' (ecCodes) on PATH."
    )


def _ensure_era5_netcdf_cache(source_path: str, logger: logging.Logger) -> str:
    """
    Convert an ERA5 monthly GRIB to a cached NetCDF, or reuse the existing cache.

    If the input is already a NetCDF file, it is returned unchanged.
    """
    source = Path(source_path)
    if source.suffix.lower() == '.nc':
        logger.info(f"  ERA5 NetCDF input: {source.name}")
        return str(source)

    cache_dir = _ensure_era5_cache_dir()
    cache_key = _build_era5_cache_key(str(source))
    cached_path = cache_dir / f"{source.stem}__{cache_key}.nc"
    tmp_path = cache_dir / f"{source.stem}__{cache_key}.tmp.nc"

    if cached_path.exists() and cached_path.stat().st_size > 0:
        logger.info(f"  ERA5 cache reuse: {source.name} -> {cached_path.name}")
        return str(cached_path)

    if tmp_path.exists():
        tmp_path.unlink()

    cmd = [
        _resolve_grib_to_netcdf_binary(),
        '-k', '4',
        '-d', '1',
        '-D', 'NC_FLOAT',
        '-o', str(tmp_path),
        str(source),
    ]

    logger.info(f"  ERA5 cache convert: {source.name} -> {cached_path.name}")
    import subprocess
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0 or not tmp_path.exists():
        if tmp_path.exists():
            tmp_path.unlink()
        detail = (proc.stderr or proc.stdout or '').strip()
        raise RuntimeError(
            f"ERA5 NetCDF conversion failed for {source.name}: {detail or 'unknown error'}"
        )

    os.replace(tmp_path, cached_path)
    return str(cached_path)


def _ensure_era5_variable_netcdf_cache(
    source_path: str,
    requested_var: str,
    logger: logging.Logger,
) -> str:
    """
    Ensure there is a NetCDF cache for one ERA5 variable from the source file.

    GRIB inputs are split to a single shortName first so `grib_to_netcdf`
    never has to process the >4 GB combined monthly bundle in one shot.
    """
    source = Path(source_path)
    if source.suffix.lower() == '.nc':
        return str(source)

    cache_dir = _ensure_era5_cache_dir()
    cache_key = _build_era5_cache_key(str(source))
    short_name = _ERA5_GRIB_SHORTNAMES.get(requested_var, requested_var)
    split_grib_path = cache_dir / f"{source.stem}__{cache_key}__{short_name}.grib"
    cached_netcdf_path = cache_dir / f"{source.stem}__{cache_key}__{requested_var}.nc"
    tmp_nc_path = cache_dir / f"{source.stem}__{cache_key}__{requested_var}.tmp.nc"

    if cached_netcdf_path.exists() and cached_netcdf_path.stat().st_size > 0:
        logger.info(
            f"  ERA5 cache reuse: {source.name} [{requested_var}] -> {cached_netcdf_path.name}"
        )
        return str(cached_netcdf_path)

    import subprocess

    if not split_grib_path.exists() or split_grib_path.stat().st_size == 0:
        split_cmd = [
            _resolve_grib_copy_binary(),
            '-w', f'shortName={short_name}',
            str(source),
            str(split_grib_path),
        ]
        logger.info(
            f"  ERA5 cache split: {source.name} [{short_name}] -> {split_grib_path.name}"
        )
        split_proc = subprocess.run(split_cmd, capture_output=True, text=True, check=False)
        if split_proc.returncode != 0 or not split_grib_path.exists():
            detail = (split_proc.stderr or split_proc.stdout or '').strip()
            raise RuntimeError(
                f"ERA5 variable split failed for {source.name} [{short_name}]: "
                f"{detail or 'unknown error'}"
            )

    if tmp_nc_path.exists():
        tmp_nc_path.unlink()

    convert_cmd = [
        _resolve_grib_to_netcdf_binary(),
        '-k', '4',
        '-d', '1',
        '-D', 'NC_FLOAT',
        '-o', str(tmp_nc_path),
        str(split_grib_path),
    ]
    logger.info(
        f"  ERA5 cache convert: {source.name} [{requested_var}] -> {cached_netcdf_path.name}"
    )
    convert_proc = subprocess.run(convert_cmd, capture_output=True, text=True, check=False)
    if convert_proc.returncode != 0 or not tmp_nc_path.exists():
        if tmp_nc_path.exists():
            tmp_nc_path.unlink()
        detail = (convert_proc.stderr or convert_proc.stdout or '').strip()
        raise RuntimeError(
            f"ERA5 NetCDF conversion failed for {source.name} [{requested_var}]: "
            f"{detail or 'unknown error'}"
        )

    os.replace(tmp_nc_path, cached_netcdf_path)

    if split_grib_path.exists():
        split_grib_path.unlink()

    return str(cached_netcdf_path)


def _coerce_datetime64_seconds(value: Any) -> np.datetime64:
    """Convert datetime-like values from netCDF readers to datetime64[s]."""
    if isinstance(value, np.datetime64):
        return value.astype('datetime64[s]')
    return np.datetime64(
        datetime(
            value.year,
            value.month,
            value.day,
            value.hour,
            value.minute,
            value.second,
            getattr(value, 'microsecond', 0),
        ),
        's',
    )


def _resolve_era5_variable_name(variable_names: List[str], requested: str) -> str:
    """Map logical ERA variable names to the actual NetCDF variable names."""
    aliases = {
        'u10': ('u10', '10u'),
        'v10': ('v10', '10v'),
        't2m': ('t2m', '2t'),
        'sp': ('sp',),
    }
    for candidate in aliases.get(requested, (requested,)):
        if candidate in variable_names:
            return candidate
    raise KeyError(
        f"ERA5 variable '{requested}' was not found. Available variables: {', '.join(sorted(variable_names))}"
    )


def _resolve_era5_coord_name(
    variable_names: List[str],
    candidates: Tuple[str, ...],
    label: str,
    netcdf_path: str,
) -> str:
    """Resolve common ERA5 coordinate aliases used by cfgrib/CDS NetCDF files."""
    for candidate in candidates:
        if candidate in variable_names:
            return candidate
    raise KeyError(
        f"Missing ERA5 {label} coordinate in {Path(netcdf_path).name}. "
        f"Tried: {', '.join(candidates)}"
    )


def _read_era5_netcdf_metadata(
    netcdf_path: str,
    required_vars: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Read coordinates, times, and resolved variable names from a NetCDF-backed ERA file."""
    import netCDF4

    with netCDF4.Dataset(netcdf_path, 'r') as ds:
        variable_names = list(ds.variables.keys())
        lat_coord = _resolve_era5_coord_name(
            variable_names, ('latitude', 'lat'), 'latitude', netcdf_path,
        )
        lon_coord = _resolve_era5_coord_name(
            variable_names, ('longitude', 'lon'), 'longitude', netcdf_path,
        )
        time_coord = _resolve_era5_coord_name(
            variable_names, ('time', 'valid_time'), 'time', netcdf_path,
        )

        lats = np.asarray(ds.variables[lat_coord][:], dtype=float)
        lons = np.asarray(ds.variables[lon_coord][:], dtype=float)

        time_var = ds.variables[time_coord]
        time_values = time_var[:]
        time_array = np.asarray(time_values)
        if np.issubdtype(time_array.dtype, np.datetime64):
            times = time_array.astype('datetime64[s]')
        else:
            time_units = getattr(time_var, 'units', None)
            if not time_units:
                raise ValueError(
                    f"ERA5 NetCDF file {Path(netcdf_path).name} has no time units "
                    f"on coordinate '{time_coord}'."
                )
            time_calendar = getattr(time_var, 'calendar', 'standard')
            datetimes = netCDF4.num2date(
                time_values,
                units=time_units,
                calendar=time_calendar,
                only_use_cftime_datetimes=False,
            )
            times = np.array(
                [_coerce_datetime64_seconds(dt) for dt in datetimes],
                dtype='datetime64[s]',
            )

        resolved_vars: Dict[str, str] = {}
        for required in required_vars or []:
            if required:
                resolved_vars[required] = _resolve_era5_variable_name(variable_names, required)

    return {
        'lats': lats,
        'lons': lons,
        'times': times,
        'var_names': resolved_vars,
    }


def _prepare_era5_source_entries(
    data_path: str,
    query_grid: Dict[str, np.ndarray],
    label: str,
    logger: logging.Logger,
    required_vars: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """
    Gather metadata for one or more ERA5 files, sorted chronologically.

    Each returned entry contains source coordinates, timestamps, subset
    indices, and precomputed bilinear weights for the geographic query grid.
    """
    entries: List[Dict[str, Any]] = []
    for fpath in _collect_era5_input_paths(data_path, label):
        try:
            if required_vars:
                cache_paths: Dict[str, str] = {}
                resolved_vars: Dict[str, str] = {}
                meta = None
                for required in required_vars:
                    if not required:
                        continue
                    variable_nc = _ensure_era5_variable_netcdf_cache(fpath, required, logger)
                    variable_meta = _read_era5_netcdf_metadata(
                        variable_nc, required_vars=[required],
                    )
                    cache_paths[required] = variable_nc
                    resolved_vars[required] = variable_meta['var_names'][required]
                    if meta is None:
                        meta = variable_meta
                    else:
                        if not np.array_equal(meta['times'], variable_meta['times']):
                            raise ValueError(
                                f"{Path(fpath).name}: ERA5 variable time axes do not match."
                            )
                        if not np.array_equal(meta['lats'], variable_meta['lats']):
                            raise ValueError(
                                f"{Path(fpath).name}: ERA5 variable latitude axes do not match."
                            )
                        if not np.array_equal(meta['lons'], variable_meta['lons']):
                            raise ValueError(
                                f"{Path(fpath).name}: ERA5 variable longitude axes do not match."
                            )
                if meta is None:
                    continue
                netcdf_path = next(iter(cache_paths.values()))
                meta['var_names'] = resolved_vars
            else:
                netcdf_path = _ensure_era5_netcdf_cache(fpath, logger)
                meta = _read_era5_netcdf_metadata(netcdf_path, required_vars=required_vars)
                cache_paths = {'_all': netcdf_path}
        except Exception as exc:
            logger.warning(
                f"{label}: skipping {Path(fpath).name} because the cached NetCDF is not usable: {exc}"
            )
            continue

        lons = meta['lons']
        lats = meta['lats']
        times = meta['times']

        if len(times) == 0:
            logger.warning(f"{label}: skipping empty ERA5 file {Path(fpath).name}")
            continue

        if lons.max() > 180:
            lons = np.where(lons > 180, lons - 360, lons)
        lons = _normalize_longitudes(lons)

        lat_mask, lon_mask = _subset_regular_grid(
            lats, lons, query_grid, f"{label} {Path(fpath).name}",
        )
        lat_idx = np.where(lat_mask)[0]
        lon_idx = np.where(lon_mask)[0]
        lat_slice = slice(int(lat_idx[0]), int(lat_idx[-1]) + 1)
        lon_slice = slice(int(lon_idx[0]), int(lon_idx[-1]) + 1)
        src_lats_sub = lats[lat_mask]
        src_lons_sub = lons[lon_mask]
        weights = _build_era5_bilinear_weights(src_lats_sub, src_lons_sub, query_grid)

        start_time = times[0].astype('datetime64[s]').tolist()
        end_time = times[-1].astype('datetime64[s]').tolist()
        n_valid = int(np.sum(weights['valid']))
        logger.info(
            f"  {Path(fpath).name}: cached as {Path(netcdf_path).name}, "
            f"{len(times)} timestep(s) from {start_time} to {end_time}, "
            f"subset {len(lat_idx)}x{len(lon_idx)}, "
            f"mapped {n_valid}/{int(np.prod(weights['shape']))} cells"
        )

        entries.append({
            'path': fpath,
            'cache_path': netcdf_path,
            'cache_paths': cache_paths,
            'times': times,
            'n_times': len(times),
            'lat_idx': lat_idx,
            'lon_idx': lon_idx,
            'lat_slice': lat_slice,
            'lon_slice': lon_slice,
            'weights': weights,
            'start_time': start_time,
            'end_time': end_time,
            'var_names': meta['var_names'],
        })

    if not entries:
        raise FileNotFoundError(f"{label}: no readable ERA5 files found.")

    entries.sort(key=lambda entry: entry['start_time'])
    logger.info(
        f"{label}: prepared {len(entries)} file(s), "
        f"{sum(entry['n_times'] for entry in entries)} timestep(s) total"
    )
    return entries


def _load_era5_variable_chunk(
    source: Any,
    var_name: str,
    time_slice: slice,
    lat_slice: slice,
    lon_slice: slice,
    logger: Optional[logging.Logger] = None,
) -> np.ndarray:
    """
    Load one ERA5 variable chunk from a NetCDF-backed ERA dataset.

    ``source`` may be an open ``netCDF4.Dataset`` or a path to the cached
    NetCDF file. The returned chunk is a plain float ndarray with masked
    values converted to NaN.
    """
    import netCDF4

    close_after = False
    if isinstance(source, (str, Path, os.PathLike)):
        ds = netCDF4.Dataset(str(source), 'r')
        close_after = True
        source_label = Path(source).name
    else:
        ds = source
        source_label = getattr(ds, 'filepath', lambda: '<memory>')()

    try:
        if var_name not in ds.variables:
            available = ', '.join(sorted(ds.variables.keys()))
            raise KeyError(
                f"{Path(str(source_label)).name}: variable '{var_name}' not found "
                f"in ERA5 NetCDF. Available variables: {available}"
            )
        raw = ds.variables[var_name][time_slice, lat_slice, lon_slice]
        return np.asarray(np.ma.filled(raw, np.nan), dtype=float)
    except Exception:
        if logger is not None:
            logger.exception(
                f"Failed to read ERA5 NetCDF chunk {var_name} "
                f"{time_slice.start}:{time_slice.stop} from {source_label}"
            )
        raise
    finally:
        if close_after:
            ds.close()


def _era5_streaming_pipeline(
    wind_path, query_grid, preset, REF,
    grid_info, output_dir, output_name, write_wnd,
    X, Y, val_path, val_crs, val_points,
    grid_crs, utm_zone, hemisphere,
    logger, progress_callback, cancel_event,
):
    """
    Legacy cfgrib-based ERA5 wind pipeline kept for reference only.

    The active production ERA path is `_era5_streaming_pipeline_v2()`,
    which uses cached NetCDF files instead of cfgrib chunk reads.

    ERA5 streaming pipeline: chunked GRIB read → bilinear remap →
    streaming file write → streaming validation extraction.

    Returns ``(n_timesteps, output_files_dict, validation_file_list)``.
    """
    import cfgrib
    import gc

    logger.info("Processing ERA5 data (streaming pipeline)...")

    # Phase 1 — open dataset to extract metadata, then close immediately
    # to avoid cfgrib/eccodes memory accumulation.
    ds = cfgrib.open_dataset(wind_path)

    lons = ds.longitude.values.copy()
    lats = ds.latitude.values.copy()
    times_all = ds.time.values.copy()
    n_times = len(times_all)

    ds.close()
    del ds
    gc.collect()

    # ERA5 uses 0-360° longitude — wrap for European domains
    if lons.max() > 180:
        lons = np.where(lons > 180, lons - 360, lons)
    lons = _normalize_longitudes(lons)

    lat_mask, lon_mask = _subset_regular_grid(lats, lons, query_grid, 'ERA5')
    lat_idx = np.where(lat_mask)[0]
    lon_idx = np.where(lon_mask)[0]
    src_lats_sub = lats[lat_mask]
    src_lons_sub = lons[lon_mask]

    logger.info(
        f"  Spatial subset: {len(lat_idx)} lats x {len(lon_idx)} lons"
    )

    # Precompute bilinear remap weights — applied to every timestep
    weights = _build_era5_bilinear_weights(src_lats_sub, src_lons_sub, query_grid)
    n_valid = int(np.sum(weights['valid']))
    logger.info(
        f"  Bilinear remap: {n_valid}/{int(np.prod(weights['shape']))} target cells mapped"
    )

    # Resolve validation points once
    val_cells = None
    val_accum = None
    if val_points or val_path:
        points = _resolve_validation_points(
            val_path, val_crs, utm_zone, hemisphere,
            validation_points=val_points,
        )
        if points:
            val_cells = []
            val_accum = []
            for pt in points:
                if grid_crs.lower().startswith('cart'):
                    pt_x, pt_y = coordinate_utils.latlon_to_utm(
                        pt['lat'], pt['lon'], utm_zone, hemisphere,
                    )
                else:
                    pt_x, pt_y = pt['lon'], pt['lat']
                cell = _find_nearest_grid_cell(float(pt_x), float(pt_y), X, Y)
                if cell is None:
                    logger.warning(f"  {pt['name']}: no valid grid cell — skipped")
                    continue
                r, c = cell
                val_cells.append((pt['name'], r, c))
                val_accum.append({
                    'name': pt['name'],
                    'datetimes': [],
                    'speeds': [],
                    'dirs': [],
                    'pressures': [],
                })
            logger.info(
                f"  Validation: {len(val_cells)} point(s) resolved"
            )

    # Open streaming writers
    writer_ctx = writer_meteo.open_meteo_writers(
        output_dir, output_name, grid_info,
        write_wnd=write_wnd, logger_ext=logger,
    )

    # Phase 2 — Chunked read → remap → write.
    # Reopen the GRIB file for each chunk and close immediately after
    # reading.  cfgrib / eccodes accumulate internal message handles
    # that are not freed until the Dataset is destroyed; reopening per
    # chunk keeps peak memory bounded regardless of file length.
    n_written = 0
    last_ts = None
    NODATA = 9999.00

    for chunk_start in range(0, n_times, _ERA5_CHUNK_SIZE):
        _check_cancel(cancel_event)
        chunk_end = min(chunk_start + _ERA5_CHUNK_SIZE, n_times)
        chunk_range = list(range(chunk_start, chunk_end))

        # Fresh dataset for this chunk
        ds = cfgrib.open_dataset(wind_path)
        u_chunk = (
            ds[preset.u_var]
            .isel(time=chunk_range, latitude=lat_idx, longitude=lon_idx)
            .values.copy()
        )
        v_chunk = (
            ds[preset.v_var]
            .isel(time=chunk_range, latitude=lat_idx, longitude=lon_idx)
            .values.copy()
        )
        p_chunk = (
            ds[preset.p_var]
            .isel(time=chunk_range, latitude=lat_idx, longitude=lon_idx)
            .values.copy()
        )
        ds.close()
        del ds

        p_chunk = p_chunk / 100.0   # Pa → mbar
        times_chunk = times_all[chunk_start:chunk_end]
        gc.collect()

        for k in range(len(chunk_range)):
            T = times_chunk[k].astype('datetime64[s]').tolist()

            # Skip exact duplicates, fail on decreasing timestamps
            if last_ts is not None:
                if T == last_ts:
                    continue
                if T < last_ts:
                    raise ValueError(
                        f"ERA5 timestamps are not monotonic: {last_ts} → {T}"
                    )
            last_ts = T

            U_remap = _apply_bilinear_remap(u_chunk[k], weights, NODATA)
            V_remap = _apply_bilinear_remap(v_chunk[k], weights, NODATA)
            P_remap = _apply_bilinear_remap(p_chunk[k], weights, NODATA)

            writer_meteo.append_meteo_timestep(
                writer_ctx, T, U_remap, V_remap, P_remap, REF,
            )

            # Accumulate validation samples
            if val_cells:
                for vi, (_, r, c) in enumerate(val_cells):
                    u_val = float(U_remap[r, c])
                    v_val = float(V_remap[r, c])
                    p_val = float(P_remap[r, c])
                    if u_val == NODATA or v_val == NODATA:
                        spd = NODATA
                        d_from = NODATA
                    else:
                        spd = np.sqrt(u_val ** 2 + v_val ** 2)
                        d_from = (270.0 - np.degrees(np.arctan2(v_val, u_val))) % 360.0
                    val_accum[vi]['datetimes'].append(T)
                    val_accum[vi]['speeds'].append(spd)
                    val_accum[vi]['dirs'].append(d_from)
                    val_accum[vi]['pressures'].append(p_val if p_val != NODATA else NODATA)

            n_written += 1

        del u_chunk, v_chunk, p_chunk
        gc.collect()

        pct = 0.05 + 0.90 * chunk_end / n_times
        progress_callback(pct)
        logger.info(
            f"  ERA5 chunk {chunk_start}-{chunk_end - 1} of {n_times} done "
            f"({n_written} timesteps written)"
        )

    # Close meteo files
    output_files = writer_meteo.close_meteo_writers(writer_ctx, logger_ext=logger)

    # Write validation CSVs
    validation_files: List[str] = []
    if val_accum:
        val_dir = Path(output_dir) / 'validation'
        val_dir.mkdir(parents=True, exist_ok=True)
        for acc in val_accum:
            csv_path = val_dir / f"{acc['name']}_validation.csv"
            writer_meteo.write_validation_csv(
                str(csv_path),
                acc['datetimes'], acc['speeds'], acc['dirs'], acc['pressures'],
            )
            validation_files.append(str(csv_path))
            logger.info(f"  Written: {csv_path.name} ({len(acc['datetimes'])} timesteps)")

    return n_written, output_files, validation_files


def _era5_streaming_pipeline_v2(
    wind_path, query_grid, preset, REF,
    grid_info, output_dir, output_name, write_wnd,
    X, Y, val_path, val_crs, val_points,
    grid_crs, utm_zone, hemisphere,
    logger, progress_callback, cancel_event,
):
    """
    ERA5 streaming pipeline: chunked GRIB read -> bilinear remap ->
    streaming file write -> streaming validation extraction.

    Supports a single GRIB file or a directory of monthly GRIB files.
    """
    import gc
    import netCDF4

    logger.info("Processing ERA5 data (streaming pipeline)...")
    entries = _prepare_era5_source_entries(
        wind_path,
        query_grid,
        f"{preset.name} wind",
        logger,
        required_vars=[preset.u_var, preset.v_var, preset.p_var],
    )
    total_source_steps = sum(entry['n_times'] for entry in entries)

    val_cells = None
    val_accum = None
    if val_points or val_path:
        points = _resolve_validation_points(
            val_path, val_crs, utm_zone, hemisphere,
            validation_points=val_points,
        )
        if points:
            val_cells = []
            val_accum = []
            for pt in points:
                if grid_crs.lower().startswith('cart'):
                    pt_x, pt_y = coordinate_utils.latlon_to_utm(
                        pt['lat'], pt['lon'], utm_zone, hemisphere,
                    )
                else:
                    pt_x, pt_y = pt['lon'], pt['lat']
                cell = _find_nearest_grid_cell(float(pt_x), float(pt_y), X, Y)
                if cell is None:
                    logger.warning(f"  {pt['name']}: no valid grid cell -> skipped")
                    continue
                r, c = cell
                val_cells.append((pt['name'], r, c))
                val_accum.append({
                    'name': pt['name'],
                    'datetimes': [],
                    'speeds': [],
                    'dirs': [],
                    'pressures': [],
                })
            logger.info(f"  Validation: {len(val_cells)} point(s) resolved")

    writer_ctx = writer_meteo.open_meteo_writers(
        output_dir, output_name, grid_info,
        write_wnd=write_wnd, logger_ext=logger,
    )

    n_written = 0
    n_processed = 0
    last_ts = None
    nodata = 9999.00

    for file_idx, entry in enumerate(entries, start=1):
        logger.info(
            f"  ERA5 file {file_idx}/{len(entries)}: {Path(entry['path']).name} "
            f"({entry['start_time']} -> {entry['end_time']})"
        )
        n_times = entry['n_times']
        u_var_name = entry['var_names'][preset.u_var]
        v_var_name = entry['var_names'][preset.v_var]
        p_var_name = entry['var_names'][preset.p_var]
        u_cache = entry['cache_paths'][preset.u_var]
        v_cache = entry['cache_paths'][preset.v_var]
        p_cache = entry['cache_paths'][preset.p_var]

        with (
            netCDF4.Dataset(u_cache, 'r') as u_ds,
            netCDF4.Dataset(v_cache, 'r') as v_ds,
            netCDF4.Dataset(p_cache, 'r') as p_ds,
        ):
            for chunk_start in range(0, n_times, _ERA5_CHUNK_SIZE):
                _check_cancel(cancel_event)
                chunk_end = min(chunk_start + _ERA5_CHUNK_SIZE, n_times)
                time_slice = slice(chunk_start, chunk_end)
                u_chunk = _load_era5_variable_chunk(
                    u_ds, u_var_name, time_slice, entry['lat_slice'], entry['lon_slice'],
                    logger=logger,
                )
                v_chunk = _load_era5_variable_chunk(
                    v_ds, v_var_name, time_slice, entry['lat_slice'], entry['lon_slice'],
                    logger=logger,
                )
                p_chunk = _load_era5_variable_chunk(
                    p_ds, p_var_name, time_slice, entry['lat_slice'], entry['lon_slice'],
                    logger=logger,
                ) / 100.0
                times_chunk = entry['times'][chunk_start:chunk_end]

                for k in range(len(times_chunk)):
                    timestamp = times_chunk[k].astype('datetime64[s]').tolist()

                    if last_ts is not None:
                        if timestamp == last_ts:
                            continue
                        if timestamp < last_ts:
                            raise ValueError(
                                f"ERA5 timestamps are not monotonic: {last_ts} -> {timestamp}"
                            )
                    last_ts = timestamp

                    u_remap = _apply_bilinear_remap(u_chunk[k], entry['weights'], nodata)
                    v_remap = _apply_bilinear_remap(v_chunk[k], entry['weights'], nodata)
                    p_remap = _apply_bilinear_remap(p_chunk[k], entry['weights'], nodata)

                    writer_meteo.append_meteo_timestep(
                        writer_ctx, timestamp, u_remap, v_remap, p_remap, REF,
                    )

                    if val_cells:
                        for vi, (_, r, c) in enumerate(val_cells):
                            u_val = float(u_remap[r, c])
                            v_val = float(v_remap[r, c])
                            p_val = float(p_remap[r, c])
                            if u_val == nodata or v_val == nodata:
                                speed = nodata
                                direction_from = nodata
                            else:
                                speed = np.sqrt(u_val ** 2 + v_val ** 2)
                                direction_from = (
                                    270.0 - np.degrees(np.arctan2(v_val, u_val))
                                ) % 360.0
                            val_accum[vi]['datetimes'].append(timestamp)
                            val_accum[vi]['speeds'].append(speed)
                            val_accum[vi]['dirs'].append(direction_from)
                            val_accum[vi]['pressures'].append(
                                p_val if p_val != nodata else nodata
                            )

                    n_written += 1

                del u_chunk, v_chunk, p_chunk
                gc.collect()

                n_processed += len(times_chunk)
                pct = 0.05 + 0.90 * n_processed / max(total_source_steps, 1)
                progress_callback(pct)
                logger.info(
                    f"  ERA5 chunk {chunk_start}-{chunk_end - 1} of {n_times} "
                    f"from {Path(entry['path']).name} done ({n_written} timesteps written)"
                )
        gc.collect()

    output_files = writer_meteo.close_meteo_writers(writer_ctx, logger_ext=logger)

    validation_files: List[str] = []
    if val_accum:
        val_dir = Path(output_dir) / 'validation'
        val_dir.mkdir(parents=True, exist_ok=True)
        for acc in val_accum:
            csv_path = val_dir / f"{acc['name']}_validation.csv"
            writer_meteo.write_validation_csv(
                str(csv_path),
                acc['datetimes'], acc['speeds'], acc['dirs'], acc['pressures'],
            )
            validation_files.append(str(csv_path))
            logger.info(f"  Written: {csv_path.name} ({len(acc['datetimes'])} timesteps)")

    return n_written, output_files, validation_files


def _era5_temperature_streaming_pipeline(
    data_path, query_grid, preset, REF,
    grid_info, output_path,
    logger, progress_callback, cancel_event=None,
):
    """
    ERA5 temperature streaming pipeline: chunked GRIB read -> bilinear remap ->
    streaming .amt write.

    Supports a single GRIB file or a directory of monthly GRIB files.
    """
    import gc
    import netCDF4

    logger.info("Processing ERA5 temperature data (streaming pipeline)...")
    entries = _prepare_era5_source_entries(
        data_path,
        query_grid,
        f"{preset.name} temperature",
        logger,
        required_vars=[preset.t_var],
    )
    total_source_steps = sum(entry['n_times'] for entry in entries)

    writer_ctx = writer_meteo.open_amt_writer(
        output_path,
        grid_info,
        nodata_value=TEMPERATURE_NODATA,
        logger_ext=logger,
    )

    n_written = 0
    n_processed = 0
    last_ts = None

    for file_idx, entry in enumerate(entries, start=1):
        logger.info(
            f"  ERA5 temperature file {file_idx}/{len(entries)}: {Path(entry['path']).name} "
            f"({entry['start_time']} -> {entry['end_time']})"
        )
        n_times = entry['n_times']
        t_var_name = entry['var_names'][preset.t_var]
        t_cache = entry['cache_paths'][preset.t_var]

        with netCDF4.Dataset(t_cache, 'r') as ds:
            for chunk_start in range(0, n_times, _ERA5_CHUNK_SIZE):
                _check_cancel(cancel_event)
                chunk_end = min(chunk_start + _ERA5_CHUNK_SIZE, n_times)
                time_slice = slice(chunk_start, chunk_end)
                t_chunk = _load_era5_variable_chunk(
                    ds, t_var_name, time_slice, entry['lat_slice'], entry['lon_slice'],
                    logger=logger,
                )
                times_chunk = entry['times'][chunk_start:chunk_end]

                for k in range(len(times_chunk)):
                    timestamp = times_chunk[k].astype('datetime64[s]').tolist()

                    if last_ts is not None:
                        if timestamp == last_ts:
                            continue
                        if timestamp < last_ts:
                            raise ValueError(
                                f"ERA5 temperature timestamps are not monotonic: {last_ts} -> {timestamp}"
                            )
                    last_ts = timestamp

                    t_remap = _apply_bilinear_remap(t_chunk[k], entry['weights'], nodata=np.nan)
                    t_c = _celsius_and_fill(t_remap)
                    writer_meteo.append_amt_timestep(writer_ctx, timestamp, t_c, REF)
                    n_written += 1

                del t_chunk
                gc.collect()

                n_processed += len(times_chunk)
                pct = 0.05 + 0.90 * n_processed / max(total_source_steps, 1)
                progress_callback(pct)
                logger.info(
                    f"  ERA5 T chunk {chunk_start}-{chunk_end - 1} of {n_times} "
                    f"from {Path(entry['path']).name} done ({n_written} timesteps written)"
                )
        gc.collect()

    final_path = writer_meteo.close_amt_writer(writer_ctx, logger_ext=logger)
    return n_written, final_path


def _process_icon(
    wind_dir, grid_path, query_grid, preset, REF, logger, progress_callback, cancel_event=None
):
    """
    Process ICON-DREAM-EU data files.
    NEW in v2 — does not exist in v1.
    """
    import cfgrib
    import xarray as xr

    logger.info("Loading ICON grid file...")

    # Load grid cell centres (in radians!)
    grid_ds = xr.open_dataset(grid_path)
    clat = np.degrees(grid_ds['clat'].values)
    clon = _normalize_longitudes(np.degrees(grid_ds['clon'].values))

    # Build KDTree (one-time operation)
    logger.info(f"Building KDTree from {len(clat)} ICON cells...")

    tree = interp.build_kdtree(clat, clon)
    progress_callback(0.1)

    # Find data files
    p_token = 'PS' if preset.name in ('ICON-DREAM-EU', 'ICON-DREAM-Global') else '_P_'
    u_files = sorted([f for f in os.listdir(wind_dir) if 'U_10M' in f and f.endswith('.grb')])
    v_files = sorted([f for f in os.listdir(wind_dir) if 'V_10M' in f and f.endswith('.grb')])
    p_files = sorted([f for f in os.listdir(wind_dir) if p_token in f and f.endswith('.grb')])

    if not u_files or not v_files or not p_files:
        raise FileNotFoundError(
            f"{preset.name}: expected U/V/P monthly GRIB files in {wind_dir}, "
            f"but found {len(u_files)} U, {len(v_files)} V, and {len(p_files)} P files."
        )
    if not (len(u_files) == len(v_files) == len(p_files)):
        raise ValueError(
            f"{preset.name}: monthly U/V/P file counts do not match "
            f"({len(u_files)}, {len(v_files)}, {len(p_files)})."
        )

    logger.info(f"Found {len(u_files)} monthly file sets")

    time_series = []
    total = len(u_files)

    for fi in range(total):
        _check_cancel(cancel_event)
        WU = cfgrib.open_dataset(os.path.join(wind_dir, u_files[fi]))
        WV = cfgrib.open_dataset(os.path.join(wind_dir, v_files[fi]))
        PS = cfgrib.open_dataset(os.path.join(wind_dir, p_files[fi]))

        # ICON has time(249) × step(3) — flatten to continuous hourly
        for t_idx, s_idx, valid_time in _collect_icon_time_records(
            WU, os.path.join(wind_dir, u_files[fi]), logger,
        ):
            _check_cancel(cancel_event)
            u_flat = WU[preset.u_var].values[t_idx, s_idx, :]
            v_flat = WV[preset.v_var].values[t_idx, s_idx, :]
            p_flat = PS[preset.p_var].values[t_idx, s_idx, :] / 100.0

            _warn_if_globally_empty(
                logger, preset.name, valid_time, [u_flat, v_flat, p_flat],
            )

            u_interp = _kdtree_interpolate_query_grid(tree, u_flat, query_grid, k=6)
            v_interp = _kdtree_interpolate_query_grid(tree, v_flat, query_grid, k=6)
            p_interp = _kdtree_interpolate_query_grid(tree, p_flat, query_grid, k=6)

            u_interp[np.isnan(u_interp)] = 9999.00
            v_interp[np.isnan(v_interp)] = 9999.00
            p_interp[np.isnan(p_interp)] = 9999.00

            time_series.append((valid_time, u_interp, v_interp, p_interp))

        progress_callback(0.1 + 0.85 * (fi + 1) / total)
        logger.info(f"  Processed ICON file {fi + 1}/{total}")

    return _sort_and_deduplicate_time_series(
        time_series, logger, f"{preset.name} wind",
    )


# ─────────────────────────────────────────────────────────────────────────────
# TEMPERATURE PROCESSORS (source-specific)
# ─────────────────────────────────────────────────────────────────────────────

def _process_cosmo_temperature(
    data_dir, aux_path, query_grid,
    preset, REF, logger, progress_callback, cancel_event=None,
):
    """
    Process COSMO REA6 2m-temperature files onto the target grid.

    Expects `data_dir` to contain T_2M.2D.YYYYMM.grb files (single variable t2m).
    Reuses the rotation+overlap mapping logic from `_process_cosmo` but only
    extracts the scalar temperature field (no vector rotation needed).
    """
    import cfgrib
    import xarray as xr

    logger.info("Loading COSMO auxiliary file for coordinate alignment...")
    ds = xr.open_dataset(
        aux_path, engine="cfgrib",
        filter_by_keys={'typeOfLevel': 'surface'},
    )
    lat_true = ds['RLAT'].values
    lon_true = _normalize_longitudes(ds['RLON'].values)

    mm1, mm2, nn1, nn2 = _find_structured_overlap_slices(
        lat_true, lon_true, query_grid, 'COSMO REA6 temperature',
    )
    lat_sub = lat_true[mm1:mm2, nn1:nn2]
    lon_sub = lon_true[mm1:mm2, nn1:nn2]

    def extract_date(filename):
        return int(filename.split('.')[-2])

    t_files = sorted(
        [f for f in os.listdir(data_dir)
         if f.startswith('T_2M') and f.endswith('.grb')],
        key=extract_date,
    )
    if not t_files:
        raise FileNotFoundError(
            f"No T_2M*.grb files found in {data_dir}. "
            f"Place COSMO T_2M.2D.YYYYMM.grb files in that directory."
        )
    logger.info(f"Found {len(t_files)} COSMO T_2M files")

    time_series = []
    for fi, fname in enumerate(t_files):
        _check_cancel(cancel_event)
        T_ds = cfgrib.open_dataset(os.path.join(data_dir, fname))
        t_var = T_ds[preset.t_var]

        for kk in range(len(T_ds.time)):
            _check_cancel(cancel_event)
            data_T = t_var.values[kk][mm1:mm2, nn1:nn2]
            T_interp = _griddata_geographic(lat_sub, lon_sub, data_T, query_grid)
            T_c = _celsius_and_fill(T_interp)
            timestamp = T_ds.time.values[kk].astype('datetime64[s]').tolist()
            time_series.append((timestamp, T_c))

        progress_callback(0.05 + 0.9 * (fi + 1) / len(t_files))
        logger.info(f"  Processed {fi + 1}/{len(t_files)}: {fname}")

    return time_series


def _process_era5_temperature(
    data_path, query_grid, preset, REF, logger, progress_callback, cancel_event=None,
):
    """
    Legacy cfgrib-based ERA5 temperature extractor kept for reference only.

    The active production ERA path is `_era5_temperature_streaming_pipeline()`,
    which uses cached NetCDF files instead of cfgrib reads.

    Process ERA5 2m-temperature. `data_path` may be a file or directory.

    Uses chunked GRIB reading + bilinear remap (same strategy as the
    wind streaming pipeline).
    """
    import cfgrib
    import gc

    logger.info("Processing ERA5 temperature data...")
    if os.path.isdir(data_path):
        files = sorted(
            f for f in os.listdir(data_path) if f.endswith(('.grib', '.grb', '.nc'))
        )
        if not files:
            raise FileNotFoundError(f"No .grib/.grb/.nc files found in {data_path}")
        paths = [os.path.join(data_path, f) for f in files]
    else:
        paths = [data_path]

    time_series = []

    for pi, fpath in enumerate(paths):
        ds = cfgrib.open_dataset(fpath)
        lons = ds.longitude.values.copy()
        lats = ds.latitude.values.copy()
        if lons.max() > 180:
            lons = np.where(lons > 180, lons - 360, lons)
        lons = _normalize_longitudes(lons)

        lat_mask, lon_mask = _subset_regular_grid(
            lats, lons, query_grid, 'ERA5 temperature',
        )
        lat_idx = np.where(lat_mask)[0]
        lon_idx = np.where(lon_mask)[0]

        weights = _build_era5_bilinear_weights(
            lats[lat_mask], lons[lon_mask], query_grid,
        )

        n_times = len(ds.time)
        for chunk_start in range(0, n_times, _ERA5_CHUNK_SIZE):
            _check_cancel(cancel_event)
            chunk_end = min(chunk_start + _ERA5_CHUNK_SIZE, n_times)
            chunk_range = list(range(chunk_start, chunk_end))

            t_chunk = (
                ds[preset.t_var]
                .isel(time=chunk_range, latitude=lat_idx, longitude=lon_idx)
                .values.copy()
            )
            times_chunk = ds.time.values[chunk_start:chunk_end].copy()

            for k in range(len(chunk_range)):
                T_remap = _apply_bilinear_remap(t_chunk[k], weights, nodata=np.nan)
                T_c = _celsius_and_fill(T_remap)
                timestamp = times_chunk[k].astype('datetime64[s]').tolist()
                time_series.append((timestamp, T_c))

            del t_chunk, times_chunk
            gc.collect()

            pct = 0.05 + 0.9 * chunk_end / n_times
            progress_callback(pct)
            logger.info(f"  ERA5 T chunk {chunk_start}-{chunk_end - 1}/{n_times}")

        ds.close()
        del ds
        gc.collect()

    return time_series


def _process_icon_temperature(
    data_dir, grid_path, query_grid,
    preset, REF, logger, progress_callback, cancel_event=None,
):
    """
    Process ICON-DREAM-EU or ICON-DREAM-Global temperature.

    EU:     single-level 't2m' (shape time × step × values).
    Global: multi-level 't' (shape time × step × generalVerticalLayer × values);
            the lowest vertical layer (highest layer index = near-surface) is used.
    """
    import cfgrib
    import xarray as xr

    logger.info("Loading ICON grid file...")
    grid_ds = xr.open_dataset(grid_path)
    clat = np.degrees(grid_ds['clat'].values)
    clon = _normalize_longitudes(np.degrees(grid_ds['clon'].values))

    logger.info(f"Building KDTree from {len(clat)} ICON cells...")
    tree = interp.build_kdtree(clat, clon)
    progress_callback(0.1)

    # ICON-DREAM-EU: T_2M files. Global: T_hourly files.
    token = 'T_2M' if preset.name == 'ICON-DREAM-EU' else '_T_'
    t_files = sorted(
        f for f in os.listdir(data_dir)
        if token in f and f.endswith(('.grb', '.grib'))
    )
    # Fallback search if naming differs slightly
    if not t_files:
        t_files = sorted(
            f for f in os.listdir(data_dir)
            if 'T' in f and f.endswith(('.grb', '.grib'))
            and 'U_10M' not in f and 'V_10M' not in f
            and 'PS' not in f and '_P_' not in f
        )
    if not t_files:
        raise FileNotFoundError(
            f"No temperature GRIB files found in {data_dir} for {preset.name}"
        )
    logger.info(f"Found {len(t_files)} ICON temperature files")

    time_series = []
    total = len(t_files)
    for fi, fname in enumerate(t_files):
        _check_cancel(cancel_event)
        T_ds = cfgrib.open_dataset(os.path.join(data_dir, fname))
        t_var = T_ds[preset.t_var]

        for t_idx, s_idx, valid_time in _collect_icon_time_records(
            T_ds, os.path.join(data_dir, fname), logger,
        ):
            _check_cancel(cancel_event)
            if preset.t_multilevel:
                # generalVerticalLayer dim present; pick the lowest near-surface level.
                n_levels = t_var.shape[2]
                t_flat = t_var.values[t_idx, s_idx, n_levels - 1, :]
            else:
                t_flat = t_var.values[t_idx, s_idx, :]

            _warn_if_globally_empty(logger, f"{preset.name} temperature", valid_time, [t_flat])

            T_interp = _kdtree_interpolate_query_grid(tree, t_flat, query_grid, k=6)
            T_c = _celsius_and_fill(T_interp)
            time_series.append((valid_time, T_c))

        progress_callback(0.1 + 0.85 * (fi + 1) / total)
        logger.info(f"  Processed ICON T file {fi + 1}/{total}: {fname}")

    return _sort_and_deduplicate_time_series(
        time_series, logger, f"{preset.name} temperature",
    )


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== backend_wind_processor.py self-test ===\n")
    all_pass = True

    def check(name, result, expected, tol=0.001):
        global all_pass
        if isinstance(result, (float, np.floating)):
            ok = abs(result - expected) < tol
        elif isinstance(result, np.ndarray):
            ok = np.allclose(result, expected, atol=tol)
        else:
            ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]  got {result!r}")
        else:
            print(f"  {name}: [{status}]")

    # Test calculate_theta
    print("calculate_theta:")
    theta = calculate_theta(np.array([[50.0]]), np.array([[10.0]]))
    check("theta is 2D", theta.shape, (1, 1))
    check("theta is finite", np.isfinite(theta[0, 0]), True)

    def _bearing_to_pole_3d(lat, lon, pollat=39.25, pollon=-162.0):
        """Grid-north bearing from 3-D vectors (independent of the closed form)."""
        def unit(la, lo):
            la, lo = np.radians(la), np.radians(lo)
            return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
        p, q = unit(lat, lon), unit(pollat, pollon)
        tangent = q - np.dot(p, q) * p
        east = np.array([-np.sin(np.radians(lon)), np.cos(np.radians(lon)), 0.0])
        north = np.cross(p, east)
        return np.degrees(np.arctan2(np.dot(tangent, east), np.dot(tangent, north)))

    pts = [(54.2, 7.9), (56.0, 3.0), (50.0, -2.0), (64.0, -20.0),
           (63.7, -37.3), (45.0, 30.0), (35.0, 15.0), (70.0, 55.0)]
    worst = max(
        abs(-np.degrees(calculate_theta(la, lo)) - _bearing_to_pole_3d(la, lo))
        for la, lo in pts
    )
    check("theta = -bearing to rotated pole (3-D check)", worst, 0.0, tol=1e-9)
    check("theta = 0 on the central meridian (18E)",
          float(np.degrees(calculate_theta(np.array([30.0, 55.0, 70.0]), 18.0)).max()),
          0.0, tol=1e-9)
    # Regression value (v1 gave ~5.9 deg here). Matches the 7.836 deg bearing
    # measured from cosmo_latlon_rot.grb at the nearest grid cell (54.22N, 7.88E).
    check("theta German Bight (54.2N, 7.9E) = 7.823 deg",
          float(np.degrees(calculate_theta(54.2, 7.9))), 7.823, tol=0.001)

    # End-to-end: a true 20 m/s westerly made grid-relative with the 3-D bearing
    # must come back out of rotate_wind_vectors(calculate_theta) unchanged.
    e2e = 0.0
    for la, lo in pts:
        d = np.radians(_bearing_to_pole_3d(la, lo))
        ur, vr = 20.0 * np.cos(d), 20.0 * np.sin(d)
        ug, vg = rotate_wind_vectors(ur, vr, calculate_theta(la, lo))
        e2e = max(e2e, abs(ug - 20.0), abs(vg))
    check("grid-relative westerly recovered (u=20, v=0)", float(e2e), 0.0, tol=1e-9)

    # Test rotation with known values
    print("\nrotate_wind_vectors:")
    u = np.array([[5.0]])
    v = np.array([[0.0]])
    theta_90 = np.array([[np.pi / 2]])  # 90° rotation
    u_rot, v_rot = rotate_wind_vectors(u, v, theta_90)
    check("90° rotation U→0", u_rot[0, 0], 0.0, tol=0.01)
    check("90° rotation V→5", v_rot[0, 0], 5.0, tol=0.01)

    # Identity rotation (theta=0)
    theta_0 = np.array([[0.0]])
    u_id, v_id = rotate_wind_vectors(u, v, theta_0)
    check("0° rotation U unchanged", u_id[0, 0], 5.0, tol=0.01)
    check("0° rotation V unchanged", v_id[0, 0], 0.0, tol=0.01)

    # Test _find_nearest_grid_cell
    print("\n_find_nearest_grid_cell:")
    X_test = np.array([[1.0, 2.0, np.nan], [1.0, 2.0, 3.0]])
    Y_test = np.array([[10.0, 10.0, np.nan], [20.0, 20.0, 20.0]])
    cell = _find_nearest_grid_cell(2.1, 19.9, X_test, Y_test)
    check("nearest cell to (2.1,19.9)", cell, (1, 1))
    cell2 = _find_nearest_grid_cell(3.0, 20.0, X_test, Y_test)
    check("nearest cell to (3,20)", cell2, (1, 2))

    # Test _extract_validation_timeseries with synthetic data
    print("\n_extract_validation_timeseries:")
    import tempfile
    import shutil

    # Create a 3x3 grid (lon, lat)
    lons = np.array([[7.0, 8.0, 9.0]] * 3)
    lats = np.array([[53.0] * 3, [54.0] * 3, [55.0] * 3])

    # 2 timesteps: constant u=3, v=4 -> speed=5, dir_from = (270 - atan2(4,3)*180/pi) % 360
    ts_test = [
        (datetime(2015, 1, 1, 0), np.full((3, 3), 3.0), np.full((3, 3), 4.0), np.full((3, 3), 1013.0)),
        (datetime(2015, 1, 1, 1), np.full((3, 3), 0.0), np.full((3, 3), 5.0), np.full((3, 3), 1010.0)),
    ]

    # Write a temp validation points file
    tmp_dir = tempfile.mkdtemp(dir=os.getcwd())
    vp_file = os.path.join(tmp_dir, 'val_pts.txt')
    with open(vp_file, 'w') as f:
        f.write("TestPt 54.0 8.0\n")

    written = _extract_validation_timeseries(
        ts_test, lons, lats, vp_file,
        'Spherical', 'Spherical', 32, 'N', tmp_dir,
        logging.getLogger('test'),
    )
    check("1 csv written", len(written), 1)
    check("csv filename", os.path.basename(written[0]), 'TestPt_validation.csv')

    with open(written[0], 'r') as f:
        csv_lines = f.readlines()
    check("csv has header + 2 rows", len(csv_lines), 3)

    # u=3, v=4 -> speed=5.0
    row1_parts = csv_lines[1].strip().split(',')
    check("speed t0 = 5.0", abs(float(row1_parts[1]) - 5.0) < 0.01, True)
    # direction_from for u=3,v=4: (270 - atan2(4,3)*180/pi) % 360 = (270 - 53.13) % 360 = 216.87
    check("dir_from t0 ~ 216.87", abs(float(row1_parts[2]) - 216.87) < 0.1, True)

    # u=0, v=5 -> speed=5.0, dir_from = (270-90) % 360 = 180
    row2_parts = csv_lines[2].strip().split(',')
    check("speed t1 = 5.0", abs(float(row2_parts[1]) - 5.0) < 0.01, True)
    check("dir_from t1 = 180", abs(float(row2_parts[2]) - 180.0) < 0.01, True)

    written_inline = _extract_validation_timeseries(
        ts_test, lons, lats, None,
        'Spherical', 'Spherical', 32, 'N', tmp_dir,
        logging.getLogger('test'),
        validation_points=[{'name': 'InlinePt', 'lat': 54.0, 'lon': 8.0}],
    )
    check("1 inline csv written", len(written_inline), 1)
    check("inline csv filename", os.path.basename(written_inline[0]), 'InlinePt_validation.csv')

    print("\nCRS-safe query grid helpers:")
    helgo_e, helgo_n = coordinate_utils.latlon_to_utm(54.1832, 7.8856, 32, 'N')
    query_grid_test = _prepare_geographic_query_grid(
        np.array([[helgo_e, np.nan]]),
        np.array([[helgo_n, np.nan]]),
        'Cartesian', 32, 'N', 'self-test', logging.getLogger('test'),
    )
    check("query grid lat", round(float(query_grid_test['query_lats'][0, 0]), 4), 54.1832, tol=0.01)
    check("query grid lon", round(float(query_grid_test['query_lons'][0, 0]), 4), 7.8856, tol=0.01)

    invalid_crs_raised = False
    try:
        _prepare_geographic_query_grid(
            np.array([[5000000.0]]), np.array([[6000000.0]]),
            'Spherical', 32, 'N', 'self-test', logging.getLogger('test'),
        )
    except ValueError:
        invalid_crs_raised = True
    check("invalid geographic bounds fail fast", invalid_crs_raised, True)

    print("\nICON month filtering:")
    class _Coord:
        def __init__(self, values):
            self.values = np.array(values)

        def __len__(self):
            return len(self.values)

    class _DummyIconDataset:
        def __init__(self, times, steps):
            self.time = _Coord(np.array(times, dtype='datetime64[s]'))
            self.step = _Coord(np.array(steps, dtype='timedelta64[h]'))

    icon_ds = _DummyIconDataset(
        ['2018-03-31T21:00:00'],
        [1, 2, 3],
    )
    records = _collect_icon_time_records(
        icon_ds, 'ICON-DREAM-EU_201804_U_10M_hourly.grb', logging.getLogger('test'),
    )
    check("trimmed out-of-month ICON steps", len(records), 1)
    check("kept first in-month ICON time", records[0][2], datetime(2018, 4, 1, 0, 0, 0))

    print("\nKDTree target masking:")
    tree = interp.build_kdtree(np.array([54.0, 55.0]), np.array([8.0, 9.0]))
    masked_query = {
        'query_lats': np.array([[54.1, np.nan]]),
        'query_lons': np.array([[8.1, np.nan]]),
        'finite_mask': np.array([[True, False]]),
    }
    kd_result = _kdtree_interpolate_query_grid(
        tree, np.array([1.0, 3.0]), masked_query, k=2,
    )
    check("kdtree valid target finite", np.isfinite(kd_result[0, 0]), True)
    check("kdtree invalid target stays NaN", np.isnan(kd_result[0, 1]), True)

    print("\nBilinear remap (ERA5 structured grid):")
    # Source: 5×5 regular grid covering lat [50,54], lon [6,10]
    src_lats = np.array([54.0, 53.0, 52.0, 51.0, 50.0])  # descending like ERA5
    src_lons = np.array([6.0, 7.0, 8.0, 9.0, 10.0])

    # Target: single point at (52.5, 8.5) — exactly between four source cells
    bl_qgrid = {
        'query_lats': np.array([[52.5]]),
        'query_lons': np.array([[8.5]]),
        'finite_mask': np.array([[True]]),
    }
    bl_w = _build_era5_bilinear_weights(src_lats, src_lons, bl_qgrid)
    check("bilinear weight sum = 1", float(
        bl_w['w00'][0, 0] + bl_w['w01'][0, 0] + bl_w['w10'][0, 0] + bl_w['w11'][0, 0]
    ), 1.0, tol=1e-12)

    # Constant field → remap preserves constant
    const_field = np.full((5, 5), 7.0)
    bl_out = _apply_bilinear_remap(const_field, bl_w)
    check("bilinear constant field", float(bl_out[0, 0]), 7.0, tol=1e-12)

    # Linear field: f(lat,lon) = lat + lon → interpolated value at (52.5,8.5) = 61.0
    lin_field = np.zeros((5, 5))
    for r in range(5):
        for c in range(5):
            lin_field[r, c] = src_lats[r] + src_lons[c]
    bl_lin = _apply_bilinear_remap(lin_field, bl_w)
    check("bilinear linear field", float(bl_lin[0, 0]), 52.5 + 8.5, tol=0.01)

    # Out-of-bounds target → NODATA
    oob_qgrid = {
        'query_lats': np.array([[49.0]]),
        'query_lons': np.array([[8.0]]),
        'finite_mask': np.array([[True]]),
    }
    oob_w = _build_era5_bilinear_weights(src_lats, src_lons, oob_qgrid)
    oob_out = _apply_bilinear_remap(const_field, oob_w, nodata=9999.0)
    check("bilinear out-of-bounds → NODATA", float(oob_out[0, 0]), 9999.0)

    print("\nDuplicate time cleanup:")
    deduped = _sort_and_deduplicate_time_series(
        [
            (datetime(2018, 4, 1, 1), 2),
            (datetime(2018, 4, 1, 0), 1),
            (datetime(2018, 4, 1, 1), 3),
        ],
        logging.getLogger('test'),
        'self-test',
    )
    check("dedupe keeps 2 timesteps", len(deduped), 2)
    check("dedupe sorts ascending", deduped[0][0], datetime(2018, 4, 1, 0))

    shutil.rmtree(tmp_dir)

    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
