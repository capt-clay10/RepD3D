# -*- coding: utf-8 -*-
"""
writer_meteo.py — Write Delft3D meteo-on-equidistant-grid files.

Writes:
    .amu — x-component wind field
    .amv — y-component wind field
    .amp — air pressure field
    .amt — air temperature field (quantity = air_temperature, unit = Celsius)

FileVersion 1.03 format, shared by COSMO, ICON, and ERA5 outputs.
The format is identical for both Delft3D-4 and Delft3D-FM.

v1 source: cosmo_wind_file_generator.py lines 311-416 (headers + write_data).

Usage:
    from writer_meteo import write_meteo_files, write_grid_data, build_meteo_header
"""

import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Callable, Optional
import os
import logging

logger = logging.getLogger(__name__)


def _time_offset_hours(timestamp: datetime, reference_time: datetime):
    """
    Hours between ``timestamp`` and ``reference_time`` for a meteo TIME header.

    Returns an ``int`` when the offset lands on a whole hour (the common case —
    keeps the header identical to hourly reanalysis output), and a ``float``
    otherwise so sub-hourly timesteps are not collapsed onto the same integer
    hour. Consecutive TIME values must stay strictly increasing (Delft3D-FLOW
    manual §B.7.1); rounding a sub-hourly cadence to whole hours would violate
    that and misread the series.
    """
    diff = (timestamp - reference_time).total_seconds() / 3600.0
    rounded = round(diff)
    if abs(diff - rounded) <= 1e-9:
        return int(rounded)
    return diff


# ─────────────────────────────────────────────────────────────────────────────
# HEADER BUILDER
# ─────────────────────────────────────────────────────────────────────────────

def _grid_origin(grid_info: Dict[str, any]) -> Tuple[float, float]:
    """
    Return the (x, y) coordinate of the lower-left data point of the meteo grid.

    The field values are sampled at the target grid nodes, so this coordinate is
    the *centre* of the lower-left meteo cell and must be written with the
    ``x_llcenter`` / ``y_llcenter`` keywords (Delft3D-FLOW manual B.7.1,
    D-Flow FM manual C.13.1: with ``x_llcorner`` the first value is placed at
    x_llcorner + dx/2, y_llcorner + dy/2, i.e. half a cell north-east of
    where it was sampled).

    The legacy ``x_llcorner`` / ``y_llcorner`` keys are still accepted from
    older callers; they always held this same sample coordinate.
    """
    x = grid_info.get('x_llcenter', grid_info.get('x_llcorner'))
    y = grid_info.get('y_llcenter', grid_info.get('y_llcorner'))
    if x is None or y is None:
        raise ValueError("grid_info must contain x_llcenter/y_llcenter")
    return float(x), float(y)


def build_meteo_header(
    n_cols: int,
    n_rows: int,
    x_llcenter: float,
    y_llcenter: float,
    dx: float,
    dy: float,
    grid_unit: str,
    quantity: str,
    unit: str,
) -> List[str]:
    """
    Build header lines for a meteo-on-equidistant-grid file.

    v1 source: cosmo_wind_file_generator.py lines 328-403.
    v1 had 5 nearly-identical copies of this header (amu, amv, amp, xwind, ywind).
    v2 generates it once with parameters.

    The origin is written as ``x_llcenter`` / ``y_llcenter``: the data values
    are sampled at the grid nodes, so the lower-left node is the centre of the
    lower-left meteo cell.  v1 (and v2 until 2026-10-14) wrote the same number
    as ``x_llcorner``, which made Delft3D place every field dx/2, dy/2 to the
    north-east of the sample points (Delft3D-FLOW manual B.7.1).

    Args:
        n_cols:     Number of columns in the grid.
        n_rows:     Number of rows in the grid.
        x_llcenter: X coordinate of the lower-left data point (cell centre).
        y_llcenter: Y coordinate of the lower-left data point (cell centre).
        dx:         Grid spacing in X direction.
        dy:         Grid spacing in Y direction.
        grid_unit:  'degree' for spherical or 'm' for Cartesian. The manual
                    (§B.7.1) accepts only these two values; the legacy
                    abbreviation 'deg' is normalized to 'degree' here.
        quantity:   'x_wind', 'y_wind', or 'air_pressure'.
        unit:       'm s-1' for wind, 'mbar' for pressure.

    Returns:
        List of header line strings.
    """
    if grid_unit == 'deg':
        grid_unit = 'degree'
    return [
        'FileVersion = 1.03',
        'Filetype = meteo_on_equidistant_grid',
        'NODATA_value = 9999.00',
        f'n_cols = {n_cols}',
        f'n_rows = {n_rows}',
        f'grid_unit = {grid_unit}',
        f'x_llcenter = {x_llcenter}',
        f'y_llcenter = {y_llcenter}',
        f'dx = {dx}',
        f'dy = {dy}',
        'n_quantity = 1',
        f'quantity1 = {quantity}',
        f'unit1 = {unit}',
    ]


# ─────────────────────────────────────────────────────────────────────────────
# GRID DATA WRITER
# ─────────────────────────────────────────────────────────────────────────────

def write_grid_data(fid, data_2d: np.ndarray, n_rows: int, n_cols: int):
    """
    Write a single timestep of 2D grid data to an open file handle.

    Format: fixed-width columns.
        Column 1: 6 chars wide (%6.2f)
        Columns 2..N-1: 8 chars wide ( %7.2f)
        Last column: 8 chars wide + newline ( %7.2f\\n)

    v1 source: cosmo_wind_file_generator.py lines 148-167.
    v1 had identical format strings for pressure and wind (despite having
    separate branches). v2 uses a single code path.

    Args:
        fid:     Open file handle.
        data_2d: 2D numpy array (n_rows × n_cols).
        n_rows:  Number of rows.
        n_cols:  Number of columns.
    """
    fmt_first_col = '%6.2f'
    fmt_mid_col = ' %7.2f'
    fmt_last_col = ' %7.2f\n'

    fmt_line = fmt_first_col + fmt_mid_col * (n_cols - 2) + fmt_last_col

    for jj in range(n_rows):
        fid.write(fmt_line % tuple(data_2d[jj, :]))


# ─────────────────────────────────────────────────────────────────────────────
# FULL METEO FILE WRITER
# ─────────────────────────────────────────────────────────────────────────────

def write_meteo_files(
    output_dir: str,
    filename_base: str,
    grid_info: Dict[str, any],
    time_series: List[Tuple[datetime, np.ndarray, np.ndarray, np.ndarray]],
    reference_time: datetime,
    write_wnd: bool = True,
    logger_ext: Optional[logging.Logger] = None,
    progress_callback: Optional[Callable[[float], None]] = None,
) -> Dict[str, str]:
    """
    Write complete .amu, .amv, .amp file set.

    Args:
        output_dir:     Output directory.
        filename_base:  Base filename (e.g., 'wind_2015').
        grid_info:      {
                            'n_cols': int,
                            'n_rows': int,
                            'x_llcenter': float,   # lower-left data point
                            'y_llcenter': float,
                            'dx': float,
                            'dy': float,
                            'grid_unit': str,  # 'm' or 'degree'
                        }
        time_series:    List of (timestamp, U_2d, V_2d, P_2d) tuples.
                        P_2d should already be in mbar.
                        NaN values should already be replaced with 9999.00.
        reference_time: Reference datetime for time headers.
        write_wnd:      Also write .xwind and .ywind files (default True).
        logger_ext:     External logger (optional, uses module logger if None).
        progress_callback: Progress callback (optional).

    Returns:
        {'amu': path, 'amv': path, 'amp': path}
        (plus 'xwind' and 'ywind' if write_wnd=True)
    """
    log = logger_ext or logger
    output_dir = Path(output_dir)
    n_cols = grid_info['n_cols']
    n_rows = grid_info['n_rows']
    x_llcenter, y_llcenter = _grid_origin(grid_info)

    # Build file paths
    file_defs = {
        'amu': (f'{filename_base}.amu', 'x_wind', 'm s-1'),
        'amv': (f'{filename_base}.amv', 'y_wind', 'm s-1'),
        'amp': (f'{filename_base}.amp', 'air_pressure', 'mbar'),
    }
    if write_wnd:
        file_defs['xwind'] = (f'xwind_{filename_base}.wnd', 'x_wind', 'm s-1')
        file_defs['ywind'] = (f'ywind_{filename_base}.wnd', 'y_wind', 'm s-1')

    # Open all files
    files = {}
    for key, (fname, quantity, unit) in file_defs.items():
        fpath = output_dir / fname
        fid = open(fpath, 'w')
        files[key] = {'fid': fid, 'path': str(fpath)}

        # Write header
        header = build_meteo_header(
            n_cols, n_rows,
            x_llcenter, y_llcenter,
            grid_info['dx'], grid_info['dy'],
            grid_info['grid_unit'], quantity, unit,
        )
        for line in header:
            fid.write(line + '\n')

    # Write timesteps
    total = len(time_series)
    for idx, (timestamp, U, V, P) in enumerate(time_series):
        # Time header
        diff_hours = _time_offset_hours(timestamp, reference_time)
        time_header = (
            f'TIME = {diff_hours} hours since '
            f'{reference_time.strftime("%Y-%m-%d %H:%M:%S")} +00:00 '
        )

        # Write time header to all files
        for key in files:
            files[key]['fid'].write(f'{time_header}\n')

        # Write data
        write_grid_data(files['amu']['fid'], U, n_rows, n_cols)
        write_grid_data(files['amv']['fid'], V, n_rows, n_cols)
        write_grid_data(files['amp']['fid'], P, n_rows, n_cols)
        if write_wnd:
            write_grid_data(files['xwind']['fid'], U, n_rows, n_cols)
            write_grid_data(files['ywind']['fid'], V, n_rows, n_cols)

        if progress_callback:
            progress_callback((idx + 1) / total)

    # Close all files
    result = {}
    for key in files:
        files[key]['fid'].close()
        result[key] = files[key]['path']

    log.info(f"Written {len(file_defs)} meteo files: {total} timesteps each")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# STREAMING METEO WRITER — open once, append per timestep, close at end
# ─────────────────────────────────────────────────────────────────────────────

def open_meteo_writers(
    output_dir: str,
    filename_base: str,
    grid_info: Dict[str, any],
    write_wnd: bool = True,
    logger_ext: Optional[logging.Logger] = None,
) -> dict:
    """
    Open .amu/.amv/.amp files, write headers, and return a context dict.

    The caller writes timesteps via ``append_meteo_timestep()`` and
    finishes with ``close_meteo_writers()``.
    """
    log = logger_ext or logger
    output_dir = Path(output_dir)
    n_cols = grid_info['n_cols']
    n_rows = grid_info['n_rows']
    x_llcenter, y_llcenter = _grid_origin(grid_info)

    file_defs = {
        'amu': (f'{filename_base}.amu', 'x_wind', 'm s-1'),
        'amv': (f'{filename_base}.amv', 'y_wind', 'm s-1'),
        'amp': (f'{filename_base}.amp', 'air_pressure', 'mbar'),
    }
    if write_wnd:
        file_defs['xwind'] = (f'xwind_{filename_base}.wnd', 'x_wind', 'm s-1')
        file_defs['ywind'] = (f'ywind_{filename_base}.wnd', 'y_wind', 'm s-1')

    files = {}
    for key, (fname, quantity, unit) in file_defs.items():
        fpath = output_dir / fname
        fid = open(fpath, 'w')
        files[key] = {'fid': fid, 'path': str(fpath)}
        header = build_meteo_header(
            n_cols, n_rows,
            x_llcenter, y_llcenter,
            grid_info['dx'], grid_info['dy'],
            grid_info['grid_unit'], quantity, unit,
        )
        for line in header:
            fid.write(line + '\n')

    log.info(f"Opened {len(file_defs)} meteo files for streaming write")
    return {
        'files': files,
        'n_cols': n_cols,
        'n_rows': n_rows,
        'write_wnd': write_wnd,
        'n_written': 0,
    }


def append_meteo_timestep(
    ctx: dict,
    timestamp: datetime,
    U: np.ndarray,
    V: np.ndarray,
    P: np.ndarray,
    reference_time: datetime,
) -> None:
    """Write one timestep to all open meteo file handles."""
    files = ctx['files']
    n_rows = ctx['n_rows']
    n_cols = ctx['n_cols']

    diff_hours = _time_offset_hours(timestamp, reference_time)
    time_header = (
        f'TIME = {diff_hours} hours since '
        f'{reference_time.strftime("%Y-%m-%d %H:%M:%S")} +00:00 '
    )

    for key in files:
        files[key]['fid'].write(f'{time_header}\n')

    write_grid_data(files['amu']['fid'], U, n_rows, n_cols)
    write_grid_data(files['amv']['fid'], V, n_rows, n_cols)
    write_grid_data(files['amp']['fid'], P, n_rows, n_cols)
    if ctx['write_wnd']:
        write_grid_data(files['xwind']['fid'], U, n_rows, n_cols)
        write_grid_data(files['ywind']['fid'], V, n_rows, n_cols)

    ctx['n_written'] += 1


def close_meteo_writers(
    ctx: dict,
    logger_ext: Optional[logging.Logger] = None,
) -> Dict[str, str]:
    """Close all open meteo file handles. Returns ``{key: path}`` dict."""
    log = logger_ext or logger
    result = {}
    for key in ctx['files']:
        ctx['files'][key]['fid'].close()
        result[key] = ctx['files'][key]['path']
    log.info(f"Closed meteo files: {ctx['n_written']} timesteps written")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# .amt FILE WRITER — Air Temperature on Equidistant Grid
# ─────────────────────────────────────────────────────────────────────────────
#
# The .amt file shares the meteo_on_equidistant_grid filetype with .amu/.amv/
# .amp, but differs from our internal .amu/.amv/.amp convention in several
# small ways that match the reference example_temp_file.amt:
#
#   - Header is delimited by `### START OF HEADER` / `### END OF HEADER`.
#   - Uses `x_llcenter` / `y_llcenter` (same as the .amu/.amv/.amp writers).
#   - NODATA_value = -999.000 (a common Delft3D temperature sentinel).
#   - grid_unit = 'degree' (long form, not 'deg').
#   - quantity1 = air_temperature, unit1 = Celsius.
#   - Data: free-form space-separated with 6 decimal places (no fixed width).

def build_amt_header(
    n_cols: int,
    n_rows: int,
    x_llcenter: float,
    y_llcenter: float,
    dx: float,
    dy: float,
    grid_unit: str = 'degree',
    nodata_value: float = -999.000,
) -> List[str]:
    """
    Build header lines for a .amt (air temperature) file.

    Matches example_temp_file.amt. Header is written between
    `### START OF HEADER` and `### END OF HEADER` lines by the caller.

    Args:
        n_cols, n_rows:     Grid dimensions.
        x_llcenter, y_llcenter: Pixel-centre of the lower-left cell.
        dx, dy:             Grid spacing.
        grid_unit:          'degree' (spherical) or 'm' (Cartesian).
        nodata_value:       Missing-data sentinel.

    Returns:
        List of header-body line strings (between the START/END markers).
    """
    return [
        'FileVersion = 1.03',
        'filetype = meteo_on_equidistant_grid',
        f'NODATA_value = {nodata_value:.3f}',
        f'n_cols = {n_cols}',
        f'n_rows = {n_rows}',
        f'grid_unit = {grid_unit}',
        f'x_llcenter = {x_llcenter}',
        f'y_llcenter = {y_llcenter}',
        f'dx = {dx}',
        f'dy = {dy}',
        'n_quantity = 1',
        'quantity1 = air_temperature',
        'unit1 = Celsius',
    ]


def write_amt_file(
    output_path: str,
    grid_info: Dict[str, any],
    time_series: List,
    reference_time: datetime,
    nodata_value: float = -999.000,
    logger_ext: Optional[logging.Logger] = None,
    progress_callback: Optional[Callable[[float], None]] = None,
) -> str:
    """
    Write a Delft3D .amt air-temperature file.

    The format is the same meteo_on_equidistant_grid envelope used by
    .amu/.amv/.amp but with the header/data styling of the reference
    example_temp_file.amt (see module docstring above).

    Args:
        output_path:    Full path to the output .amt file.
        grid_info:      {
                            'n_cols', 'n_rows',
                            'x_llcenter', 'y_llcenter',  (or 'x_llcorner'/'y_llcorner')
                            'dx', 'dy',
                            'grid_unit' ('degree' or 'm'),
                        }
        time_series:    List of (timestamp, T_2d) tuples. T_2d is the 2D
                        temperature field in **Celsius** with missing cells
                        already replaced by nodata_value.
        reference_time: Reference datetime for the TIME headers.
        nodata_value:   Missing-data sentinel (default -999.000, matching
                        the reference file).
        logger_ext:     External logger (optional).
        progress_callback: Progress callback (optional).

    Returns:
        The output_path written.
    """
    log = logger_ext or logger
    output_path = str(output_path)

    x_llcenter, y_llcenter = _grid_origin(grid_info)

    # The .amt example uses 'degree' for spherical grids rather than 'deg'.
    raw_unit = grid_info.get('grid_unit', 'degree')
    grid_unit = 'degree' if raw_unit in ('deg', 'degree') else raw_unit

    header_body = build_amt_header(
        grid_info['n_cols'], grid_info['n_rows'],
        x_llcenter, y_llcenter,
        grid_info['dx'], grid_info['dy'],
        grid_unit=grid_unit,
        nodata_value=nodata_value,
    )

    total = len(time_series)
    with open(output_path, 'w') as fid:
        fid.write('### START OF HEADER\n')
        for line in header_body:
            fid.write(line + '\n')
        fid.write('### END OF HEADER\n')

        for idx, (timestamp, T) in enumerate(time_series):
            diff_hours = (timestamp - reference_time).total_seconds() / 3600.0
            # Match the reference file's "TIME =   H.H hours since …" style.
            fid.write(
                f'TIME = {diff_hours:5.1f} hours since '
                f'{reference_time.strftime("%Y-%m-%d %H:%M:%S")} +00:00\n'
            )

            # Free-form space-separated values, 6 decimal places, one row per line.
            for jj in range(grid_info['n_rows']):
                row = T[jj, :]
                fid.write(' '.join(f'{v:.6f}' for v in row) + '\n')

            if progress_callback:
                progress_callback((idx + 1) / total)

    log.info(f"Written .amt file: {output_path} ({total} timesteps)")
    return output_path


def open_amt_writer(
    output_path: str,
    grid_info: Dict[str, any],
    nodata_value: float = -999.000,
    logger_ext: Optional[logging.Logger] = None,
) -> dict:
    """
    Open a .amt writer, emit the header once, and return a context dict.

    Timesteps are written with ``append_amt_timestep()`` and the file is
    finalized with ``close_amt_writer()``.
    """
    log = logger_ext or logger
    output_path = str(output_path)

    x_llcenter, y_llcenter = _grid_origin(grid_info)

    raw_unit = grid_info.get('grid_unit', 'degree')
    grid_unit = 'degree' if raw_unit in ('deg', 'degree') else raw_unit

    header_body = build_amt_header(
        grid_info['n_cols'], grid_info['n_rows'],
        x_llcenter, y_llcenter,
        grid_info['dx'], grid_info['dy'],
        grid_unit=grid_unit,
        nodata_value=nodata_value,
    )

    fid = open(output_path, 'w')
    fid.write('### START OF HEADER\n')
    for line in header_body:
        fid.write(line + '\n')
    fid.write('### END OF HEADER\n')

    log.info("Opened .amt file for streaming write")
    return {
        'fid': fid,
        'path': output_path,
        'n_rows': grid_info['n_rows'],
        'n_written': 0,
    }


def append_amt_timestep(
    ctx: dict,
    timestamp: datetime,
    T: np.ndarray,
    reference_time: datetime,
) -> None:
    """Append one temperature timestep to an open .amt stream."""
    fid = ctx['fid']
    diff_hours = (timestamp - reference_time).total_seconds() / 3600.0
    fid.write(
        f'TIME = {diff_hours:5.1f} hours since '
        f'{reference_time.strftime("%Y-%m-%d %H:%M:%S")} +00:00\n'
    )
    for jj in range(ctx['n_rows']):
        row = T[jj, :]
        fid.write(' '.join(f'{v:.6f}' for v in row) + '\n')
    ctx['n_written'] += 1


def close_amt_writer(
    ctx: dict,
    logger_ext: Optional[logging.Logger] = None,
) -> str:
    """Close an open .amt stream and return the written path."""
    log = logger_ext or logger
    ctx['fid'].close()
    log.info(f"Closed .amt file: {ctx['n_written']} timesteps written")
    return ctx['path']


# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION POINT CSV WRITER
# ─────────────────────────────────────────────────────────────────────────────

def write_validation_csv(
    output_path: str,
    datetimes: List[datetime],
    wind_speed: List[float],
    direction_from: List[float],
    pressure: List[float],
) -> str:
    """
    Write a single-point validation time series CSV.

    Columns: datetime, wind_speed_m_s, direction_from_deg_N, pressure_mbar

    Args:
        output_path:    Full path for the output CSV.
        datetimes:      List of datetime objects.
        wind_speed:     Wind speed in m/s.
        direction_from: Meteorological "wind from" direction in degrees (0=N).
        pressure:       Air pressure in mbar.

    Returns:
        The output_path written.
    """
    output_path = str(output_path)
    with open(output_path, 'w') as f:
        f.write('datetime,wind_speed_m_s,direction_from_deg_N,pressure_mbar\n')
        for dt, ws, wd, p in zip(datetimes, wind_speed, direction_from, pressure):
            f.write(f'{dt.strftime("%Y-%m-%d %H:%M:%S")},{ws:.4f},{wd:.2f},{p:.2f}\n')
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .wnd FILE WRITER — Single-point wind time series (speed + direction)
# ─────────────────────────────────────────────────────────────────────────────

def write_wnd_point(
    output_path: str,
    times_min: List[float],
    wind_speed: List[float],
    wind_direction: List[float],
) -> str:
    """
    Write a single-point Delft3D .wnd time series.

    Format (matches point_timeseries_wind.wnd):
        - 3 columns: time [min since reference], wind speed [m/s],
          wind direction [deg from North, nautical convention]
        - Scientific notation, precision 7, exponent 3 digits
        - Positive values padded with a leading space
        - Double-space column separator
        - CRLF line endings
    """
    from formatting_utils import list_to_sci, pad_positive

    if not (len(times_min) == len(wind_speed) == len(wind_direction)):
        raise ValueError("times, speed, direction must have equal length")

    t_str = list_to_sci(times_min, precision=7, exp_digits=3)
    s_str = list_to_sci(wind_speed, precision=7, exp_digits=3)
    d_str = list_to_sci(wind_direction, precision=7, exp_digits=3)

    with open(output_path, 'w', newline='') as f:
        for t, s, d in zip(t_str, s_str, d_str):
            line = f'{pad_positive(t)} {pad_positive(s)} {pad_positive(d)}\r\n'
            f.write(line)

    logger.info(f"Written .wnd file: {output_path} ({len(times_min)} rows)")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile

    print("=== writer_meteo.py self-test ===\n")
    all_pass = True

    def check(name, result, expected):
        global all_pass
        ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]  got {result!r}")
        else:
            print(f"  {name}: [{status}]")

    # ── build_meteo_header ───────────────────────────────────────────────
    print("build_meteo_header:")
    header = build_meteo_header(48, 37, 269913, 5899220, 6000, 6000, 'm', 'x_wind', 'm s-1')
    check("header line count", len(header), 13)
    check("has FileVersion", header[0], 'FileVersion = 1.03')
    check("has n_cols", 'n_cols = 48' in header[3], True)
    check("has quantity", 'quantity1 = x_wind' in header[11], True)
    check("origin written as llcenter", header[6], 'x_llcenter = 269913')
    check("origin y written as llcenter", header[7], 'y_llcenter = 5899220')
    check("legacy corner keys still map to the sample point",
          _grid_origin({'x_llcorner': 1, 'y_llcorner': 2}), (1.0, 2.0))

    # grid_unit 'deg' must be normalized to the manual-valid 'degree'
    header_deg = build_meteo_header(48, 37, 0, 0, 0.1, 0.1, 'deg', 'x_wind', 'm s-1')
    check("normalizes deg -> degree", 'grid_unit = degree' in header_deg[5], True)
    header_m = build_meteo_header(48, 37, 0, 0, 6000, 6000, 'm', 'x_wind', 'm s-1')
    check("keeps grid_unit m", 'grid_unit = m' in header_m[5], True)

    # sub-hourly TIME offsets must not collapse onto whole hours
    check("whole-hour offset is int", _time_offset_hours(datetime(2015, 1, 1, 6), datetime(2015, 1, 1)), 6)
    check("half-hour offset stays fractional",
          abs(_time_offset_hours(datetime(2015, 1, 1, 0, 30), datetime(2015, 1, 1)) - 0.5) < 1e-9, True)

    # ── write_grid_data ──────────────────────────────────────────────────
    print("\nwrite_grid_data:")
    import io
    buf = io.StringIO()
    test_data = np.array([[1.23, 4.56, 7.89], [10.0, 20.0, 30.0]])
    write_grid_data(buf, test_data, 2, 3)
    output = buf.getvalue()
    lines = output.split('\n')
    lines = [l for l in lines if l]  # Remove empty lines only
    check("grid data 2 rows", len(lines), 2)
    # First value should be 6 chars wide (format %6.2f → '  1.23')
    check("first val formatted", lines[0][:6], '  1.23')

    # ── write_meteo_files ────────────────────────────────────────────────
    print("\nwrite_meteo_files:")
    tmpdir = tempfile.mkdtemp(dir=os.getcwd())
    ref = datetime(2015, 1, 1)

    # Create 3 timesteps of 4x5 grid data
    ts_data = []
    for h in range(3):
        t = datetime(2015, 1, 1, h)
        U = np.ones((4, 5)) * (h + 1.0)
        V = np.ones((4, 5)) * (h + 0.5)
        P = np.ones((4, 5)) * (1013.0 + h)
        ts_data.append((t, U, V, P))

    grid = {
        'n_cols': 5, 'n_rows': 4,
        'x_llcenter': 269913, 'y_llcenter': 5899220,
        'dx': 6000, 'dy': 6000, 'grid_unit': 'm',
    }

    result = write_meteo_files(tmpdir, 'test', grid, ts_data, ref, write_wnd=False)
    check("returns amu path", 'amu' in result, True)
    check("returns amv path", 'amv' in result, True)
    check("returns amp path", 'amp' in result, True)

    # Verify .amu content
    with open(result['amu'], 'r') as f:
        amu_content = f.read()
    check("amu has FileVersion", 'FileVersion = 1.03' in amu_content, True)
    check("amu has x_wind", 'x_wind' in amu_content, True)
    check("amu has TIME headers", amu_content.count('TIME ='), 3)
    check("amu has reference", '2015-01-01 00:00:00' in amu_content, True)
    check("amu has x_llcenter", 'x_llcenter = 269913.0' in amu_content, True)
    check("amu has no llcorner keyword", 'llcorner' in amu_content, False)

    # Verify .amp content
    with open(result['amp'], 'r') as f:
        amp_content = f.read()
    check("amp has air_pressure", 'air_pressure' in amp_content, True)
    check("amp has mbar unit", 'mbar' in amp_content, True)

    # Cleanup
    import shutil
    shutil.rmtree(tmpdir)

    # ── write_amt_file ───────────────────────────────────────────────────
    print("\nwrite_amt_file:")
    tmp_amt = tempfile.mkdtemp(dir=os.getcwd())
    ref_amt = datetime(2023, 8, 1)
    t_ts = []
    for h in (0, 3, 6):
        t = datetime(2023, 8, 1, h)
        T = np.ones((3, 4)) * (15.0 + h * 0.5)
        t_ts.append((t, T))

    grid_amt = {
        'n_cols': 4, 'n_rows': 3,
        'x_llcenter': 8.875, 'y_llcenter': 53.375,
        'dx': 0.25, 'dy': 0.25, 'grid_unit': 'degree',
    }
    amt_path = os.path.join(tmp_amt, 'test.amt')
    write_amt_file(amt_path, grid_amt, t_ts, ref_amt)

    with open(amt_path, 'r') as f:
        amt_content = f.read()
    check("amt has START OF HEADER", '### START OF HEADER' in amt_content, True)
    check("amt has END OF HEADER", '### END OF HEADER' in amt_content, True)
    check("amt has air_temperature quantity", 'quantity1 = air_temperature' in amt_content, True)
    check("amt has Celsius unit", 'unit1 = Celsius' in amt_content, True)
    check("amt has x_llcenter", 'x_llcenter = 8.875' in amt_content, True)
    check("amt has -999.000 NODATA", 'NODATA_value = -999.000' in amt_content, True)
    check("amt has grid_unit = degree", 'grid_unit = degree' in amt_content, True)
    check("amt has 3 TIME headers", amt_content.count('TIME ='), 3)
    check("amt has correct reference", '2023-08-01 00:00:00' in amt_content, True)
    # Value format should have 6 decimals
    check("amt has 6-decimal values", '.000000' in amt_content, True)
    shutil.rmtree(tmp_amt)

    print("\nstreaming .amt writer:")
    tmp_amt_stream = tempfile.mkdtemp(dir=os.getcwd())
    stream_amt_path = os.path.join(tmp_amt_stream, 'stream_test.amt')
    ctx_amt = open_amt_writer(stream_amt_path, grid_amt)
    for t, T in t_ts:
        append_amt_timestep(ctx_amt, t, T, ref_amt)
    close_amt_writer(ctx_amt)
    with open(stream_amt_path, 'r') as f:
        stream_amt_content = f.read()
    check("streaming amt has 3 TIME headers", stream_amt_content.count('TIME ='), 3)
    check("streaming amt has END OF HEADER", '### END OF HEADER' in stream_amt_content, True)
    shutil.rmtree(tmp_amt_stream)

    # ── write_validation_csv ─────────────────────────────────────────────
    print("\nwrite_validation_csv:")
    tmpdir2 = tempfile.mkdtemp(dir=os.getcwd())
    csv_path = os.path.join(tmpdir2, 'P1_validation.csv')
    dts = [datetime(2015, 1, 1, h) for h in range(3)]
    ws = [5.0, 6.5, 3.2]
    wd = [180.0, 270.0, 45.0]
    pr = [1013.25, 1010.0, 1015.5]
    write_validation_csv(csv_path, dts, ws, wd, pr)
    with open(csv_path, 'r') as f:
        csv_lines = f.readlines()
    check("csv header", csv_lines[0].strip(),
          'datetime,wind_speed_m_s,direction_from_deg_N,pressure_mbar')
    check("csv row count", len(csv_lines), 4)  # header + 3 data rows
    check("csv row1 datetime", csv_lines[1].startswith('2015-01-01 00:00:00'), True)
    check("csv row2 has ws 6.5", '6.5000' in csv_lines[2], True)
    shutil.rmtree(tmpdir2)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
