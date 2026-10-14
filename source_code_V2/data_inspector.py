# -*- coding: utf-8 -*-
"""
data_inspector.py — Inspect NetCDF, GRIB, and MAT files for RepD3D v2.

Used by the GUI for:
    1. Listing available variables when a user loads a data file
    2. Auto-detecting which preset matches based on variable names
    3. Showing dimension shapes, coordinate ranges, and metadata
    4. Validating that required variables exist before extraction

This module opens files read-only and closes them immediately after
inspection — it does NOT load the full data into memory.

Usage:
    from data_inspector import (
        inspect_file, auto_detect_bc_preset, validate_variables,
    )
"""

import numpy as np
from pathlib import Path
from typing import Dict, Optional, Any, List
import logging

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# UNIFIED FILE INSPECTION
# ─────────────────────────────────────────────────────────────────────────────

def inspect_file(path: str) -> Dict[str, Any]:
    """
    Inspect a data file and return its metadata.

    Auto-detects file type from extension and delegates to the appropriate
    reader. Supports .nc (NetCDF), .grb/.grib (GRIB), and .mat (MATLAB).

    Args:
        path: Path to the data file.

    Returns:
        {
            'file_type':    'netcdf' | 'grib' | 'matlab',
            'file_size_mb': float,
            'variables':    {name: {'shape': tuple, 'units': str, 'long_name': str, 'standard_name': str}, ...},
            'coords':       {name: {'shape': tuple, 'min': ..., 'max': ..., 'units': str}, ...},
            'dimensions':   {name: int, ...},
            'global_attrs': dict,  (NetCDF only)
        }

    Raises:
        ValueError: If file type is not supported.
        FileNotFoundError: If file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    file_size_mb = path.stat().st_size / (1024 * 1024)
    suffix = path.suffix.lower()

    if suffix in ('.nc', '.nc4', '.netcdf'):
        result = _inspect_netcdf(str(path))
    elif suffix in ('.grb', '.grib', '.grb2', '.grib2'):
        result = _inspect_grib(str(path))
    elif suffix == '.mat':
        result = _inspect_mat(str(path))
    else:
        raise ValueError(
            f"Unsupported file type: {suffix}. "
            f"Supported: .nc, .grb, .grib, .mat"
        )

    result['file_size_mb'] = round(file_size_mb, 1)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# NETCDF INSPECTION
# ─────────────────────────────────────────────────────────────────────────────

def _inspect_netcdf(path: str) -> Dict[str, Any]:
    """Inspect a NetCDF file."""
    import xarray as xr

    ds = xr.open_dataset(path)

    result = {
        'file_type': 'netcdf',
        'variables': {},
        'coords': {},
        'dimensions': dict(ds.sizes),
        'global_attrs': dict(ds.attrs),
    }

    for name, var in ds.data_vars.items():
        result['variables'][name] = {
            'shape': var.shape,
            'units': var.attrs.get('units', ''),
            'long_name': var.attrs.get('long_name', ''),
            'standard_name': var.attrs.get('standard_name', ''),
        }

    for name, coord in ds.coords.items():
        vals = coord.values
        coord_info = {
            'shape': vals.shape,
            'units': coord.attrs.get('units', ''),
            'long_name': coord.attrs.get('long_name', ''),
            'standard_name': coord.attrs.get('standard_name', ''),
        }
        if np.issubdtype(vals.dtype, np.number):
            coord_info['min'] = float(np.nanmin(vals))
            coord_info['max'] = float(np.nanmax(vals))
            coord_info['dtype'] = str(vals.dtype)
        elif np.issubdtype(vals.dtype, np.datetime64):
            coord_info['min'] = str(vals[0])
            coord_info['max'] = str(vals[-1])
            coord_info['dtype'] = 'datetime64'
        else:
            coord_info['dtype'] = str(vals.dtype)

        result['coords'][name] = coord_info

    ds.close()
    return result


# ─────────────────────────────────────────────────────────────────────────────
# GRIB INSPECTION
# ─────────────────────────────────────────────────────────────────────────────

def _inspect_grib(path: str) -> Dict[str, Any]:
    """Inspect a GRIB file using cfgrib."""
    import cfgrib

    result = {
        'file_type': 'grib',
        'variables': {},
        'coords': {},
        'dimensions': {},
        'global_attrs': {},
    }

    try:
        ds = cfgrib.open_dataset(path)
    except Exception:
        # Some GRIB files need specific filter_by_keys
        ds = cfgrib.open_dataset(
            path, backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}}
        )

    result['dimensions'] = dict(ds.sizes)

    for name, var in ds.data_vars.items():
        result['variables'][name] = {
            'shape': var.shape,
            'units': var.attrs.get('units', var.attrs.get('GRIB_units', '')),
            'long_name': var.attrs.get(
                'long_name', var.attrs.get('GRIB_name', '')
            ),
            'standard_name': var.attrs.get(
                'standard_name', var.attrs.get('GRIB_cfName', '')
            ),
        }

    for name, coord in ds.coords.items():
        vals = coord.values
        coord_info = {
            'shape': vals.shape,
            'units': coord.attrs.get('units', coord.attrs.get('GRIB_units', '')),
            'long_name': coord.attrs.get(
                'long_name', coord.attrs.get('GRIB_name', '')
            ),
            'standard_name': coord.attrs.get(
                'standard_name', coord.attrs.get('GRIB_cfName', '')
            ),
        }
        if np.issubdtype(vals.dtype, np.number) and vals.size > 0:
            coord_info['min'] = float(np.nanmin(vals))
            coord_info['max'] = float(np.nanmax(vals))
        result['coords'][name] = coord_info

    ds.close()
    return result


# ─────────────────────────────────────────────────────────────────────────────
# MATLAB .mat INSPECTION
# ─────────────────────────────────────────────────────────────────────────────

def _inspect_mat(path: str) -> Dict[str, Any]:
    """
    Inspect a MATLAB .mat file (for Delft3D equidistant grid).

    Typical structure: mat['data']['X'][0][0] and mat['data']['Y'][0][0]
    containing 2D grid coordinate arrays.
    """
    import scipy.io

    mat = scipy.io.loadmat(path)

    result = {
        'file_type': 'matlab',
        'variables': {},
        'coords': {},
        'dimensions': {},
        'global_attrs': {},
    }

    for key in mat:
        if key.startswith('_'):
            continue
        val = mat[key]
        if hasattr(val, 'shape'):
            var_info = {'shape': val.shape, 'dtype': str(val.dtype)}
            # Check for struct arrays (common in Delft3D grid .mat files)
            if val.dtype.names:
                var_info['fields'] = list(val.dtype.names)
                # Try to extract sub-field shapes
                for field_name in val.dtype.names:
                    try:
                        sub = val[field_name][0][0]
                        if hasattr(sub, 'shape'):
                            var_info[f'{field_name}_shape'] = sub.shape
                            if np.issubdtype(sub.dtype, np.number):
                                var_info[f'{field_name}_min'] = float(np.nanmin(sub))
                                var_info[f'{field_name}_max'] = float(np.nanmax(sub))
                    except (IndexError, TypeError):
                        pass
            result['variables'][key] = var_info

    return result


# ─────────────────────────────────────────────────────────────────────────────
# PRESET AUTO-DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def auto_detect_bc_preset(file_info: Dict[str, Any]) -> Optional[str]:
    """
    Try to match a file's variables to a known BC preset.

    Checks for signature variables unique to each data source.

    Args:
        file_info: Output from inspect_file().

    Returns:
        Preset name string (e.g., 'EasyGSH') or None if no match.
    """
    var_names = set(file_info.get('variables', {}).keys())

    # EasyGSH: German variable names with Mesh2_face_ prefix
    if 'Mesh2_face_Wasserstand_2d' in var_names:
        return 'EasyGSH'
    if 'Mesh2_face_signifikante_Wellenhoehe_2d' in var_names:
        return 'EasyGSH'

    # TrilaWatt: English variable names
    if 'sea_surface_height_2d' in var_names:
        return 'TrilaWatt'
    if 'significant_wave_height' in var_names:
        return 'TrilaWatt'

    # Copernicus SSH/Currents
    if 'zos' in var_names:
        return 'Copernicus SSH/Currents'
    if 'uo' in var_names and 'vo' in var_names:
        return 'Copernicus SSH/Currents'

    # Copernicus Waves
    if 'VHM0' in var_names:
        return 'Copernicus Waves'

    return None


def auto_detect_wind_preset(file_info: Dict[str, Any]) -> Optional[str]:
    """
    Try to match a GRIB file's variables to a known Wind preset.

    This is less reliable than BC detection since variable names (u10, v10, sp)
    are shared across sources. Detection relies on grid structure clues.

    Args:
        file_info: Output from inspect_file().

    Returns:
        Preset name string or None.
    """
    dims = file_info.get('dimensions', {})
    coords = file_info.get('coords', {})

    # ICON: has 'values' dimension with ~659k points
    if 'values' in dims and dims.get('values', 0) > 500000:
        return 'ICON-DREAM-EU'

    # COSMO: has 2D latitude/longitude arrays and rotated pole
    var_names = set(file_info.get('variables', {}).keys())
    if 'RLAT' in var_names or 'RLON' in var_names:
        return 'COSMO REA6'

    # ERA5: regular lat/lon with 721 × 1440 grid
    if dims.get('latitude', 0) == 721 and dims.get('longitude', 0) == 1440:
        return 'ERA5 (.nc)'

    return None


# ─────────────────────────────────────────────────────────────────────────────
# VARIABLE VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

def validate_variables(
    file_info: Dict[str, Any],
    required_vars: List[str],
) -> Dict[str, bool]:
    """
    Check which required variables exist in the file.

    Args:
        file_info:     Output from inspect_file().
        required_vars: List of variable names to check.

    Returns:
        Dict mapping each variable name to True/False.

    Example:
        >>> validate_variables(info, ['zos', 'uo', 'vo'])
        {'zos': True, 'uo': True, 'vo': False}
    """
    var_names = set(file_info.get('variables', {}).keys())
    return {var: var in var_names for var in required_vars}


def get_required_vars_for_extraction(preset_name: str, extract_type: str) -> List[str]:
    """
    Get the list of variable names needed for a specific extraction type.

    Args:
        preset_name:  Name of the BC preset.
        extract_type: 'water_level', 'waves', 'currents', 'sediments', or 'salinity'.

    Returns:
        List of variable name strings.
    """
    from presets import get_bc_preset

    preset = get_bc_preset(preset_name)
    required = []

    if extract_type == 'water_level':
        if preset.water_level_var:
            required.append(preset.water_level_var)

    elif extract_type == 'waves':
        if preset.wave_hs_var:
            required.append(preset.wave_hs_var)
        if preset.wave_tp_var:
            required.append(preset.wave_tp_var)
        if preset.wave_dir_format == 'vector':
            if preset.wave_dir_x_var:
                required.append(preset.wave_dir_x_var)
            if preset.wave_dir_y_var:
                required.append(preset.wave_dir_y_var)
        elif preset.wave_dir_format == 'degrees':
            if preset.wave_dir_var:
                required.append(preset.wave_dir_var)
        if preset.wave_dir_spread_var:
            required.append(preset.wave_dir_spread_var)

    elif extract_type == 'currents':
        if preset.current_u_var:
            required.append(preset.current_u_var)
        if preset.current_v_var:
            required.append(preset.current_v_var)

    elif extract_type == 'sediments':
        required.extend(preset.sediment_conc_vars)

    elif extract_type == 'salinity':
        if preset.salinity_var:
            required.append(preset.salinity_var)

    return required


# ─────────────────────────────────────────────────────────────────────────────
# HUMAN-READABLE SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

def format_file_summary(file_info: Dict[str, Any]) -> str:
    """
    Format file inspection results as a human-readable string for the GUI log.

    Args:
        file_info: Output from inspect_file().

    Returns:
        Multi-line string summary.
    """
    lines = []
    lines.append(f"File type: {file_info.get('file_type', 'unknown')}")
    lines.append(f"Size: {file_info.get('file_size_mb', '?')} MB")

    dims = file_info.get('dimensions', {})
    if dims:
        dim_str = ', '.join(f"{k}={v}" for k, v in dims.items())
        lines.append(f"Dimensions: {dim_str}")

    coords = file_info.get('coords', {})
    if coords:
        lines.append("Coordinates:")
        for name, info in coords.items():
            parts = [f"  {name}: shape={info.get('shape', '?')}"]
            if 'min' in info:
                parts.append(f"range=[{info['min']}, {info['max']}]")
            lines.append(' '.join(parts))

    variables = file_info.get('variables', {})
    if variables:
        lines.append(f"Variables ({len(variables)}):")
        for name, info in variables.items():
            parts = [f"  {name}: shape={info.get('shape', '?')}"]
            if info.get('units'):
                parts.append(f"[{info['units']}]")
            if info.get('long_name'):
                parts.append(f"— {info['long_name']}")
            lines.append(' '.join(parts))

    return '\n'.join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== data_inspector.py self-test ===\n")
    all_pass = True

    def check(name, result, expected):
        global all_pass
        ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]  got {result!r}, expected {expected!r}")
        else:
            print(f"  {name}: [{status}]")

    # ── auto_detect_bc_preset ────────────────────────────────────────────
    print("auto_detect_bc_preset:")

    check("EasyGSH WL", auto_detect_bc_preset(
        {'variables': {'Mesh2_face_Wasserstand_2d': {}}}
    ), 'EasyGSH')

    check("EasyGSH waves", auto_detect_bc_preset(
        {'variables': {'Mesh2_face_signifikante_Wellenhoehe_2d': {}}}
    ), 'EasyGSH')

    check("TrilaWatt WL", auto_detect_bc_preset(
        {'variables': {'sea_surface_height_2d': {}}}
    ), 'TrilaWatt')

    check("TrilaWatt waves", auto_detect_bc_preset(
        {'variables': {'significant_wave_height': {}}}
    ), 'TrilaWatt')

    check("Copernicus SSH", auto_detect_bc_preset(
        {'variables': {'zos': {}, 'uo': {}, 'vo': {}}}
    ), 'Copernicus SSH/Currents')

    check("Copernicus Waves", auto_detect_bc_preset(
        {'variables': {'VHM0': {}, 'VTPK': {}, 'VMDR': {}}}
    ), 'Copernicus Waves')

    check("Unknown", auto_detect_bc_preset(
        {'variables': {'temperature': {}}}
    ), None)

    # ── auto_detect_wind_preset ──────────────────────────────────────────
    print("\nauto_detect_wind_preset:")

    check("ICON (large values dim)", auto_detect_wind_preset(
        {'dimensions': {'values': 659156}, 'variables': {}, 'coords': {}}
    ), 'ICON-DREAM-EU')

    check("COSMO (RLAT/RLON)", auto_detect_wind_preset(
        {'dimensions': {}, 'variables': {'RLAT': {}, 'RLON': {}}, 'coords': {}}
    ), 'COSMO REA6')

    check("ERA5 (721x1440)", auto_detect_wind_preset(
        {'dimensions': {'latitude': 721, 'longitude': 1440}, 'variables': {}, 'coords': {}}
    ), 'ERA5 (.nc)')

    check("Unknown wind", auto_detect_wind_preset(
        {'dimensions': {}, 'variables': {'u10': {}}, 'coords': {}}
    ), None)

    # ── validate_variables ───────────────────────────────────────────────
    print("\nvalidate_variables:")
    info = {'variables': {'zos': {}, 'uo': {}}}
    val = validate_variables(info, ['zos', 'uo', 'vo'])
    check("zos present", val['zos'], True)
    check("uo present", val['uo'], True)
    check("vo missing", val['vo'], False)

    # ── get_required_vars_for_extraction ──────────────────────────────────
    print("\nget_required_vars_for_extraction:")
    req_wl = get_required_vars_for_extraction('EasyGSH', 'water_level')
    check("EasyGSH WL vars", req_wl, ['Mesh2_face_Wasserstand_2d'])

    req_waves = get_required_vars_for_extraction('EasyGSH', 'waves')
    check("EasyGSH wave vars count", len(req_waves), 5)
    check("EasyGSH wave Hs in list", 'Mesh2_face_signifikante_Wellenhoehe_2d' in req_waves, True)

    req_cop_waves = get_required_vars_for_extraction('Copernicus Waves', 'waves')
    check("Copernicus wave dir", 'VMDR' in req_cop_waves, True)
    # Copernicus has no dir spread, so it shouldn't be in the list
    check("Copernicus no spread var", len([v for v in req_cop_waves if 'spread' in v.lower()]), 0)

    req_curr = get_required_vars_for_extraction('TrilaWatt', 'currents')
    check("TrilaWatt current vars", req_curr, ['current_velocity_2d_x', 'current_velocity_2d_y'])

    # ── format_file_summary ──────────────────────────────────────────────
    print("\nformat_file_summary:")
    summary = format_file_summary({
        'file_type': 'netcdf',
        'file_size_mb': 3400.5,
        'dimensions': {'lon': 372, 'lat': 218, 'time': 26280},
        'coords': {'lat': {'shape': (218,), 'min': 53.21, 'max': 55.14}},
        'variables': {'zos': {'shape': (26280, 218, 372), 'units': 'm', 'long_name': 'sea surface height'}},
    })
    check("summary has file type", 'netcdf' in summary, True)
    check("summary has size", '3400.5 MB' in summary, True)
    check("summary has variable", 'zos' in summary, True)

    # ── Unsupported file type ────────────────────────────────────────────
    print("\nError handling:")
    try:
        inspect_file('test.xyz')
        check("bad extension raises", False, True)
    except (ValueError, FileNotFoundError):
        check("bad extension raises", True, True)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
