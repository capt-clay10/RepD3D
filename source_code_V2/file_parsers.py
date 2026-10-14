# -*- coding: utf-8 -*-
"""
file_parsers.py — Read Delft3D-4 and Delft3D-FM model definition files.

Consolidates and fixes parsing code from:
    - extract_from_d3d_files.py (184 lines in v1)
    - Inline .mdf parsing in bct_generator.py (lines 201-236)
    - Inline .mdf parsing in main.py (lines 87-99, 152-176)
    - Inline .mdw parsing in mdw_writer.py

New in v2:
    - parse_mdu(): Delft3D-FM .mdu master definition file (INI format)
    - parse_pli(): Delft3D-FM .pli polyline boundary file
    - extract_mdf_timing(): Derives simulation timing from parsed .mdf
    - All parsers use `with open()` and pathlib (no leaked file handles)

Bugs fixed from v1:
    1. `break` before `file.close()` — dead code (bct_generator.py line 213)
    2. False error printed for every non-matching line in value_from_txt_file
    3. CRS detection only worked with \\r\\n line endings (extract_crs_type)
    4. .grd parser relied on \\r\\n — now handles both \\n and \\r\\n

Usage:
    from file_parsers import (
        parse_mdf, extract_mdf_timing,
        parse_mdu,
        parse_bnd, parse_grd, extract_boundary_coordinates,
        parse_pli, parse_mdw,
    )
"""

import math
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Dict, Optional, Tuple, Any
import logging

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D-4: .mdf (Master Definition Flow)
# ─────────────────────────────────────────────────────────────────────────────

def parse_mdf(path: str) -> Dict[str, str]:
    """
    Parse a Delft3D-4 .mdf master definition file.

    The .mdf format is key = value pairs, one per line. Values may be
    wrapped in quotes or hash marks. Comment lines start with '*'.

    Returns:
        dict mapping stripped keys to stripped values.

    Example .mdf content:
        Ident  = #Delft3D#
        Itdate = #2015-01-01#
        Tstart = 0.0000000e+000
        Tstop  = 4.3200000e+004

    Example return:
        {'Ident': '#Delft3D#', 'Itdate': '#2015-01-01#', 'Tstart': '0.0000000e+000', ...}

    v1 source: bct_generator.py lines 45-57 (value_from_txt_file) opened
    the file 3 separate times for Tstart, Tstop, Itdate. This version
    reads the file once and returns all key-value pairs.

    FIX: v1 had `break` before `file1.close()` (dead code) and printed
    a false error for every non-matching line.
    """
    result = {}
    path = Path(path)
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('*'):
                continue
            if '=' in line:
                key, _, val = line.partition('=')
                result[key.strip()] = val.strip()
    return result


def extract_mdf_timing(mdf_data: Dict[str, str]) -> Dict[str, Any]:
    """
    Extract simulation timing parameters from a parsed .mdf dictionary.

    Handles the Itdate format variants found in practice:
        - #2015-01-01#      (hash-wrapped with hyphens)
        - '2015-01-01'      (quote-wrapped)
        - 2015-01-01        (bare)

    v1 source: bct_generator.py lines 233-261, main.py lines 152-177.
    v1 used: ref_time_unedited[1:11] to strip the '#' and then
    ref_time.replace('-', '') to get the reference_time for bct headers.

    Args:
        mdf_data: Dict from parse_mdf().

    Returns:
        {
            'reference_date':     '20150101'              (for .bct header)
            'reference_date_iso': '2015-01-01'            (ISO format)
            'tstart':             0.0                     (minutes since ref)
            'tstop':              43200.0                 (minutes since ref)
            'start_datetime':     datetime(2015, 1, 1, 0) (absolute start)
            'end_datetime':       datetime(2015, 1, 31, 0)(absolute end)
        }
    """
    tstart = float(mdf_data['Tstart'])
    tstop = float(mdf_data['Tstop'])

    # Parse Itdate — strip quotes, hashes, whitespace
    itdate_raw = mdf_data['Itdate']
    itdate_clean = itdate_raw.strip().strip("'").strip('#').strip()
    # Extract just the date portion (first 10 chars: YYYY-MM-DD)
    ref_date_iso = itdate_clean[:10]
    reference_date = ref_date_iso.replace('-', '')

    # Compute absolute start/end datetimes
    ref_dt = datetime.strptime(ref_date_iso, '%Y-%m-%d')
    start_dt = ref_dt + timedelta(minutes=tstart)
    end_dt = ref_dt + timedelta(minutes=tstop)

    return {
        'reference_date': reference_date,
        'reference_date_iso': ref_date_iso,
        'tstart': tstart,
        'tstop': tstop,
        'start_datetime': start_dt,
        'end_datetime': end_dt,
    }


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D-4: .bnd (Boundary Definition)
# ─────────────────────────────────────────────────────────────────────────────

def parse_bnd(path: str) -> List[Dict[str, Any]]:
    """
    Parse a Delft3D-4 .bnd boundary definition file.

    The .bnd format uses fixed-width columns:
        Columns  1-20: Boundary name (padded with spaces)
        Column  21+:   Space-separated fields:
                       [type] [data_type] [m_a] [n_a] [m_b] [n_b] [refl_coeff] ...

    Returns:
        List of dicts, each with:
            'bnd_name': str  — boundary name
            'm_a', 'n_a':    int — start point grid indices
            'm_b', 'n_b':    int — end point grid indices

    v1 source: extract_from_d3d_files.py lines 1-32 (extract_bnd_grd_indices).
    """
    result = []
    path = Path(path)
    with open(path, 'r') as f:
        content = f.read()

    for record in content.split('\n'):
        record = record.rstrip('\r')  # Handle \r\n and \n
        if not record.strip():
            continue

        # First 20 chars = boundary name (fixed-width), rest = fields
        bnd_name = record[:20].strip()
        fields = record[20:].split()

        if len(fields) < 6:
            logger.warning(f"Skipping malformed .bnd line: {record!r}")
            continue

        # fields: [type, data_type, m_a, n_a, m_b, n_b, ...]
        result.append({
            'bnd_name': bnd_name,
            'bnd_type': fields[0] if len(fields) > 0 else '',
            'data_type': fields[1] if len(fields) > 1 else '',
            'm_a': int(fields[2]),
            'n_a': int(fields[3]),
            'm_b': int(fields[4]),
            'n_b': int(fields[5]),
        })

    return result


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D-4: .grd (Grid File)
# ─────────────────────────────────────────────────────────────────────────────

def parse_grd(path: str) -> Dict[str, Any]:
    """
    Parse a Delft3D-4 .grd structured grid file.

    The .grd format consists of:
        - Header lines (Coordinate System, Missing Value, grid dimensions)
        - X coordinates per ETA line (N blocks of M values)
        - Y coordinates per ETA line (N blocks of M values)

    Blocks are separated by 'ETA=' markers.

    Returns:
        {
            'x_values':      list of lists (one per ETA line, each containing M+1 strings)
            'y_values':      list of lists (same structure)
            'm':             int (grid dimension M)
            'n':             int (grid dimension N)
            'missing_value': str (e.g., '0.0000000e+000' or '0')
            'crs_type':      'Cartesian' or 'Spherical'
        }

    v1 source: extract_from_d3d_files.py lines 35-78 (read_grd) and
    lines 80-88 (extract_crs_type).

    FIX: v1 split header on '\\r\\n' only, failing on Unix files.
    This version handles both \\r\\n and \\n.
    """
    path = Path(path)
    with open(path, 'r') as f:
        content = f.read()

    # ── Detect CRS type from header ──────────────────────────────────────
    crs_type = 'Cartesian'
    # Normalise line endings for reliable header parsing
    header_lines = content.split('\n')
    for line in header_lines[:10]:
        line_clean = line.strip().rstrip('\r')
        if 'Spherical' in line_clean:
            crs_type = 'Spherical'
            break
        if 'Cartesian' in line_clean:
            crs_type = 'Cartesian'
            break

    # ── Split on ETA= markers ────────────────────────────────────────────
    records = content.split('ETA=')

    x_values = []
    y_values = []
    m = 0
    n = 0
    missing_value = '0'

    for i, record in enumerate(records):
        if i == 0:
            # Header block — extract missing value and grid dimensions
            # Normalise line endings
            header_split = record.replace('\r\n', '\n').split('\n')
            # Remove empty lines
            header_split = [line for line in header_split if line.strip()]

            # Look for 'Missing Value' line
            for line in header_split:
                parts = line.split()
                if parts and parts[0] == 'Missing':
                    if len(parts) >= 4:
                        missing_value = parts[3]
                    break

            # Grid dimensions line has exactly 2 positive integers: "  M  N"
            # Must distinguish from other numeric lines like "0 0 0"
            for line in reversed(header_split):
                parts = line.split()
                if len(parts) == 2:
                    try:
                        m_cand = int(parts[0])
                        n_cand = int(parts[1])
                        if m_cand > 0 and n_cand > 0:
                            m = m_cand
                            n = n_cand
                            break
                    except ValueError:
                        continue

        elif i <= n:
            # X values (ETA blocks 1 through N)
            x_values.append(record.split())
        elif i <= 2 * n:
            # Y values (ETA blocks N+1 through 2N)
            y_values.append(record.split())

    logger.debug(f"Grid dimensions: M = 1...{m}, N = 1...{n}, CRS = {crs_type}")

    return {
        'x_values': x_values,
        'y_values': y_values,
        'm': m,
        'n': n,
        'missing_value': missing_value,
        'crs_type': crs_type,
    }


def extract_boundary_coordinates(
    path_grd: str,
    bnd_list: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Extract x/y coordinates for boundary start/end points from a .grd file
    using the m/n indices from a .bnd file.

    For each boundary in bnd_list, this function:
    1. Looks up the grid coordinates at (m_a, n_a) and (m_b, n_b)
    2. Handles boundaries that extend 1 cell beyond the grid edge
    3. Falls back to adjacent cells if coordinates are 0 (cut grids)

    Args:
        path_grd: Path to .grd file.
        bnd_list: Output from parse_bnd() — list of dicts with m/n indices.

    Returns:
        Input list enriched with 'x_a', 'y_a', 'x_b', 'y_b' keys.

    Raises:
        ValueError: If any boundary index is more than 1 cell outside grid range.

    v1 source: extract_from_d3d_files.py lines 90-184.
    FIX: v1 silently set inputs_valid = True even for out-of-range indices
    (line 122). v2 raises a clear error.
    FIX: v1 only checked == 0.0 for missing coordinates. v2 reads the
    grid's declared missing value (e.g. -999.999) and uses axis-aware,
    pair-safe fallback that picks x and y from the same candidate cell.
    """
    grd = parse_grd(path_grd)
    x_values = grd['x_values']
    y_values = grd['y_values']
    m = grd['m']
    n = grd['n']

    try:
        miss_val = float(grd['missing_value'])
    except (ValueError, TypeError):
        miss_val = 0.0

    def _is_missing_coord(val: float) -> bool:
        if val == 0.0:
            return True
        if math.isclose(val, miss_val, abs_tol=1.0):
            return True
        return False

    def _clamp_index(idx: int, max_val: int, label: str) -> int:
        """Clamp boundary index that is 1 beyond grid to the grid edge."""
        if idx < 1 or idx > max_val + 1:
            raise ValueError(
                f"Boundary index {label}={idx} is outside valid range "
                f"(1 to {max_val + 1})"
            )
        if idx == max_val + 1:
            logger.debug(f"Clamping {label}={idx} to grid max {max_val}")
            return max_val
        return idx

    def _grid_xy(mi: int, ni: int) -> Tuple[float, float]:
        """Read (x, y) from the grid arrays at 1-based indices (m, n)."""
        return float(x_values[ni - 1][mi]), float(y_values[ni - 1][mi])

    def _find_valid_cell(
        mi: int, ni: int, m_a: int, m_b: int, n_a: int, n_b: int,
        name: str, point: str,
    ) -> Tuple[float, float]:
        """Search inward for the nearest cell where both x and y are valid."""
        max_search = 6
        constant_m = (m_a == m_b)
        constant_n = (n_a == n_b)

        if constant_m:
            for step in range(1, max_search):
                cm = mi - step
                if cm < 1:
                    break
                cx, cy = _grid_xy(cm, ni)
                if not _is_missing_coord(cx) and not _is_missing_coord(cy):
                    logger.debug(
                        f"Fallback for {name}.{point}: constant-m, "
                        f"used m={cm}, n={ni} → ({cx}, {cy})"
                    )
                    return cx, cy

        if constant_n:
            for step in range(1, max_search):
                cn = ni - step
                if cn < 1:
                    break
                cx, cy = _grid_xy(mi, cn)
                if not _is_missing_coord(cx) and not _is_missing_coord(cy):
                    logger.debug(
                        f"Fallback for {name}.{point}: constant-n, "
                        f"used m={mi}, n={cn} → ({cx}, {cy})"
                    )
                    return cx, cy

        for step in range(1, max_search):
            for dm, dn in [(-step, 0), (0, -step), (-step, -step),
                           (step, 0), (0, step)]:
                cm = mi + dm
                cn = ni + dn
                if cm < 1 or cn < 1 or cm > m or cn > n:
                    continue
                cx, cy = _grid_xy(cm, cn)
                if not _is_missing_coord(cx) and not _is_missing_coord(cy):
                    logger.debug(
                        f"Fallback for {name}.{point}: neighborhood, "
                        f"used m={cm}, n={cn} → ({cx}, {cy})"
                    )
                    return cx, cy

        raise ValueError(
            f"No valid grid coordinates found near {name}.{point} "
            f"(m={mi}, n={ni}). All nearby cells contain the missing "
            f"value ({miss_val})."
        )

    for bnd in bnd_list:
        name = bnd['bnd_name']
        m_a = _clamp_index(bnd['m_a'], m, f"{name}.m_a")
        m_b = _clamp_index(bnd['m_b'], m, f"{name}.m_b")
        n_a = _clamp_index(bnd['n_a'], n, f"{name}.n_a")
        n_b = _clamp_index(bnd['n_b'], n, f"{name}.n_b")

        for point, mi, ni in [('a', m_a, n_a), ('b', m_b, n_b)]:
            x, y = _grid_xy(mi, ni)
            if _is_missing_coord(x) or _is_missing_coord(y):
                x, y = _find_valid_cell(mi, ni, m_a, m_b, n_a, n_b,
                                        name, point)
            bnd[f'x_{point}'] = x
            bnd[f'y_{point}'] = y

    return bnd_list


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D-FM: .mdu (Master Definition Unstructured)
# ─────────────────────────────────────────────────────────────────────────────

def parse_mdu(path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a Delft3D-FM .mdu master definition file (INI format).

    The .mdu format uses [Section] headers with key = value pairs.

    Returns:
        dict of dicts: {section_name: {key: value, ...}, ...}

    Example .mdu content:
        [General]
        Program        = D-Flow FM
        Version        = 1.2.100

        [time]
        RefDate        = 20170101
        Tunit          = S
        DtUser         = 300.0
        TStart         = 0.0
        TStop          = 86400.0

    Example return:
        {
            'General': {'Program': 'D-Flow FM', 'Version': '1.2.100'},
            'time': {'RefDate': '20170101', 'Tunit': 'S', ...},
        }
    """
    result = {}
    current_section = None
    path = Path(path)

    with open(path, 'r') as f:
        for line in f:
            line = line.strip().rstrip('\r')

            # Skip empty lines and comments
            if not line or line.startswith('#') or line.startswith('*'):
                continue

            # Section header
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                result[current_section] = {}
            elif '=' in line and current_section is not None:
                key, _, val = line.partition('=')
                # Handle inline comments (# after value)
                if '#' in val:
                    val = val[:val.index('#')]
                result[current_section][key.strip()] = val.strip()

    return result


def extract_mdu_timing(mdu_data: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
    """
    Extract simulation timing from a parsed .mdu dictionary.

    .mdu stores time differently from .mdf:
        - RefDate: YYYYMMDD (no hyphens, no hashes)
        - Tunit: time unit ('S' for seconds, 'M' for minutes, 'H' for hours)
        - TStart / TStop: in whatever Tunit specifies

    Args:
        mdu_data: Dict-of-dicts from parse_mdu().

    Returns:
        Same structure as extract_mdf_timing():
        {
            'reference_date':     '20170101'
            'reference_date_iso': '2017-01-01'
            'tstart':             0.0              (always in minutes for consistency)
            'tstop':              1440.0            (always in minutes)
            'start_datetime':     datetime object
            'end_datetime':       datetime object
        }
    """
    # Find the time section (case-insensitive search)
    time_section = None
    for key in mdu_data:
        if key.lower() == 'time':
            time_section = mdu_data[key]
            break

    if time_section is None:
        raise ValueError("No [time] section found in .mdu file")

    ref_date_str = time_section.get('RefDate', '').strip()
    tunit = time_section.get('Tunit', 'S').strip().upper()
    tstart_raw = float(time_section.get('TStart', '0'))
    tstop_raw = float(time_section.get('TStop', '0'))

    # Convert to minutes (our internal standard)
    unit_to_minutes = {'S': 1 / 60.0, 'M': 1.0, 'H': 60.0}
    factor = unit_to_minutes.get(tunit, 1 / 60.0)
    tstart = tstart_raw * factor
    tstop = tstop_raw * factor

    # Parse reference date
    if len(ref_date_str) == 8:
        ref_date_iso = f"{ref_date_str[:4]}-{ref_date_str[4:6]}-{ref_date_str[6:8]}"
    elif '-' in ref_date_str:
        ref_date_iso = ref_date_str[:10]
    else:
        raise ValueError(f"Cannot parse RefDate: {ref_date_str!r}")

    reference_date = ref_date_iso.replace('-', '')
    ref_dt = datetime.strptime(ref_date_iso, '%Y-%m-%d')
    start_dt = ref_dt + timedelta(minutes=tstart)
    end_dt = ref_dt + timedelta(minutes=tstop)

    return {
        'reference_date': reference_date,
        'reference_date_iso': ref_date_iso,
        'tstart': tstart,
        'tstop': tstop,
        'start_datetime': start_dt,
        'end_datetime': end_dt,
    }


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D-FM: .pli (Polyline Boundary)
# ─────────────────────────────────────────────────────────────────────────────

def parse_pli(path: str) -> List[Dict[str, Any]]:
    """
    Parse a Delft3D-FM .pli polyline boundary file.

    The .pli format defines named polylines that specify boundary locations:

        BND_1
        2    2
        6.50000    54.00000
        6.50000    54.50000
        BND_2
        3    2
        7.00000    53.00000
        7.50000    53.50000
        8.00000    54.00000

    Format per boundary:
        Line 1: Name
        Line 2: N_points  N_dimensions
        Lines 3..N+2: coordinate values (x y [z])

    Returns:
        List of dicts:
        [
            {'name': 'BND_1', 'points': [(6.5, 54.0), (6.5, 54.5)]},
            {'name': 'BND_2', 'points': [(7.0, 53.0), (7.5, 53.5), (8.0, 54.0)]},
        ]

    Points are (x, y) = (longitude, latitude) in geographic coordinates.
    """
    result = []
    path = Path(path)

    with open(path, 'r') as f:
        lines = [line.strip().rstrip('\r') for line in f if line.strip()]

    i = 0
    while i < len(lines):
        # Line 1: Name
        name = lines[i]
        i += 1

        if i >= len(lines):
            logger.warning(f"Unexpected end of .pli file after name '{name}'")
            break

        # Line 2: N_points  N_dimensions
        dims = lines[i].split()
        n_points = int(dims[0])
        i += 1

        # Lines 3..N+2: coordinates
        points = []
        for _ in range(n_points):
            if i >= len(lines):
                logger.warning(f"Unexpected end of .pli file in boundary '{name}'")
                break
            coords = lines[i].split()
            points.append((float(coords[0]), float(coords[1])))
            i += 1

        result.append({'name': name, 'points': points})

    return result


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D WAVE: .mdw (Master Definition Wave)
# ─────────────────────────────────────────────────────────────────────────────

def parse_mdw(path: str) -> Dict[str, Any]:
    """
    Parse a Delft3D WAVE .mdw file (INI-like format with [Boundary] blocks).

    Returns:
        {
            'header_lines': list of raw lines before first [Boundary],
            'header_end_line': int — line number where [Boundary] blocks start,
            'boundaries': [
                {'Name': 'BND_1', 'Definition': 'xy-coordinates', ...},
                ...
            ],
        }

    v1 source: mdw_writer.py lines 16-74.
    """
    result = {
        'header_lines': [],
        'header_end_line': 0,
        'boundaries': [],
    }
    path = Path(path)

    with open(path, 'r') as f:
        lines = f.readlines()

    in_boundary = False
    current_bnd = {}

    for i, line in enumerate(lines):
        stripped = line.strip().rstrip('\r')

        if stripped == '[Boundary]':
            # Save previous boundary if exists
            if current_bnd:
                result['boundaries'].append(current_bnd)
            current_bnd = {}
            in_boundary = True

            # Record where the first boundary starts
            if result['header_end_line'] == 0:
                result['header_end_line'] = i

        elif in_boundary and '=' in stripped:
            key, _, val = stripped.partition('=')
            current_bnd[key.strip()] = val.strip()

        elif not in_boundary:
            result['header_lines'].append(line)

    # Don't forget the last boundary
    if current_bnd:
        result['boundaries'].append(current_bnd)

    return result


# ─────────────────────────────────────────────────────────────────────────────
# BOUNDARY COORDINATE CSV (intermediate format used in v1 workflow)
# ─────────────────────────────────────────────────────────────────────────────

def write_boundary_csv(
    bnd_data_list: List[Dict[str, Any]],
    output_path: str,
) -> str:
    """
    Write boundary coordinates to a CSV file.

    This creates the intermediate CSV that v1 used as input for BC extraction.
    In v2 this is optional (backends work directly with bnd_data_list),
    but is preserved for compatibility and debugging.

    Format (no header):
        BND_1_a,566463.0,5935270.0
        BND_1_b,566500.0,5935300.0

    v1 source: output_methods.py lines 5-17.
    """
    import csv
    output_path = str(output_path)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        for bnd in bnd_data_list:
            writer.writerow([f"{bnd['bnd_name']}_a", bnd['x_a'], bnd['y_a']])
            writer.writerow([f"{bnd['bnd_name']}_b", bnd['x_b'], bnd['y_b']])

    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# WIND VALIDATION POINTS
# ─────────────────────────────────────────────────────────────────────────────

def parse_validation_points(
    path: str,
    crs_type: str = 'Spherical',
    utm_zone: int = 32,
    hemisphere: str = 'N',
) -> List[Dict[str, Any]]:
    """
    Parse a 3-column validation-points file for wind time-series extraction.

    Accepts whitespace- or comma-separated files with the columns:
        name   y   x

    Where (y, x) is interpreted depending on `crs_type`:
        - 'Spherical' → (lat, lon) in degrees (WGS84)
        - 'Cartesian' → (northing, easting) in metres for the given UTM zone

    Header lines starting with '#' are ignored. The first row is also skipped
    if its second column is non-numeric (treated as a header row).

    Returns:
        [
            {'name': 'P1', 'lat': 54.0, 'lon': 8.0,
             'easting': <metres or None>, 'northing': <metres or None>},
            ...
        ]
        Geographic lat/lon are always populated; for Cartesian inputs they
        are computed via coordinate_utils.utm_to_latlon.

    Raises:
        FileNotFoundError, ValueError on parse failure.
    """
    from coordinate_utils import utm_to_latlon

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Validation points file not found: {p}")

    points: List[Dict[str, Any]] = []
    with open(p, 'r') as f:
        raw_lines = f.readlines()

    rows: List[List[str]] = []
    for line in raw_lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if ',' in line:
            parts = [c.strip() for c in line.split(',')]
        else:
            parts = line.split()
        if len(parts) < 3:
            continue
        rows.append(parts[:3])

    if not rows:
        raise ValueError(f"No data rows found in validation points file: {p}")

    # Skip header row if its second column is non-numeric
    try:
        float(rows[0][1])
    except ValueError:
        rows = rows[1:]

    if not rows:
        raise ValueError(f"Validation points file has only a header row: {p}")

    for parts in rows:
        name = parts[0]
        try:
            y_val = float(parts[1])
            x_val = float(parts[2])
        except ValueError as e:
            raise ValueError(
                f"Could not parse numeric coords for point '{name}' in {p}: {e}"
            )

        if crs_type.lower().startswith('cart'):
            northing = y_val
            easting = x_val
            lat, lon = utm_to_latlon(easting, northing, utm_zone, hemisphere)
            points.append({
                'name': name,
                'lat': float(lat),
                'lon': float(lon),
                'easting': easting,
                'northing': northing,
            })
        else:
            lat = y_val
            lon = x_val
            points.append({
                'name': name,
                'lat': lat,
                'lon': lon,
                'easting': None,
                'northing': None,
            })

    return points


# ─────────────────────────────────────────────────────────────────────────────
# DELFT3D-4: .bct / .bcc Section Headers (for synthetic BC header replay)
# ─────────────────────────────────────────────────────────────────────────────

def _parse_bc_section_table(path: str) -> Dict[str, Any]:
    """
    Shared parser for .bct / .bcc files that use the "table-name" header
    convention. Header length varies with parameter count (11 lines for
    water level / salinity / concentration, 13 for currents), so this
    scans each section until the `records-in-table` marker.

    Returns:
        {
          'reference_date': '20150101',
          'sections': [
            {
              'boundary_name': str,
              'parameters': [{'label': str, 'unit': str}, ...],  # excludes the 'time' row
              'records_in_table': int,
              'n_data_columns': int,   # number of value columns excluding time
              'header_lines': [...],   # raw header lines (verbatim for replay)
            }, ...
          ],
        }
    """
    lines = Path(path).read_text().splitlines()
    sections: List[Dict[str, Any]] = []
    top_ref_date: Optional[str] = None

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lstrip().startswith('table-name'):
            header_end = i
            while header_end < len(lines):
                if lines[header_end].lstrip().startswith('records-in-table'):
                    break
                header_end += 1
            if header_end >= len(lines):
                raise ValueError(
                    f"Unterminated header block at line {i + 1} in {path}"
                )

            header_block = lines[i:header_end + 1]

            boundary_name = ''
            parameters: List[Dict[str, str]] = []
            records_in_table = 0
            ref_date = ''

            for hl in header_block:
                hls = hl.strip()
                if hls.startswith('location'):
                    boundary_name = _extract_quoted(hls)
                elif hls.startswith('reference-time'):
                    ref_date = hls.split(None, 1)[1].strip()
                elif hls.startswith('records-in-table'):
                    try:
                        records_in_table = int(hls.split()[1])
                    except (IndexError, ValueError):
                        records_in_table = 0
                elif hls.startswith('parameter'):
                    label = _extract_quoted(hls)
                    unit = _extract_unit(hls)
                    parameters.append({'label': label, 'unit': unit})

            data_params = [p for p in parameters
                           if 'time' not in p['label'].lower()]

            if not top_ref_date and ref_date:
                top_ref_date = ref_date

            sections.append({
                'boundary_name': boundary_name,
                'parameters': data_params,
                'records_in_table': records_in_table,
                'n_data_columns': len(data_params),
                'header_lines': header_block,
            })

            i = header_end + 1 + max(records_in_table, 0)
            continue
        i += 1

    return {
        'reference_date': top_ref_date or '',
        'sections': sections,
    }


def _extract_quoted(line: str) -> str:
    start = line.find("'")
    if start < 0:
        return ''
    end = line.find("'", start + 1)
    if end < 0:
        return ''
    return line[start + 1:end].strip()


def _extract_unit(line: str) -> str:
    key = 'unit'
    idx = line.rfind(key)
    if idx < 0:
        return ''
    return _extract_quoted(line[idx:])


def parse_bct_sections(path: str) -> Dict[str, Any]:
    """
    Parse a Delft3D-4 .bct file's section headers (water level or currents).

    Each .bct section has an 11-line header block followed by `records-in-table`
    data rows. This parser skips the data rows and returns the header metadata
    plus the verbatim header lines so the synthetic BC module can replay them
    when the user supplies an existing file instead of .mdf+.bnd+.grd.

    Returns:
        See _parse_bc_section_table(). Each section has n_data_columns == 2
        for water level, 4 for currents.
    """
    return _parse_bc_section_table(path)


def parse_bcc_sections(path: str) -> Dict[str, Any]:
    """
    Parse a Delft3D-4 .bcc file's section headers (sediment or salinity).

    Same 11-line header envelope as .bct. For sediments the section list
    repeats per fraction — fraction identity is *positional*, not stored
    in the header. Length of sections = n_fractions * n_boundaries.
    """
    return _parse_bc_section_table(path)


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile
    import os

    print("=== file_parsers.py self-test ===\n")
    all_pass = True

    def check(name, result, expected):
        global all_pass
        ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]")
            print(f"    got:      {result!r}")
            print(f"    expected: {expected!r}")
        else:
            print(f"  {name}: [{status}]")

    # ── parse_mdf ────────────────────────────────────────────────────────
    print("parse_mdf:")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.mdf', delete=False) as f:
        f.write("* Comment line\n")
        f.write("Ident  = #Delft3D-FLOW#\n")
        f.write("Itdate = #2015-01-01#\n")
        f.write("Tstart = 0.0000000e+000\n")
        f.write("Tstop  = 4.3200000e+004\n")
        f.write("Dt     = 1.0000000e+000\n")
        mdf_path = f.name

    mdf = parse_mdf(mdf_path)
    check("Itdate", mdf['Itdate'], '#2015-01-01#')
    check("Tstart", mdf['Tstart'], '0.0000000e+000')
    check("Tstop", mdf['Tstop'], '4.3200000e+004')
    check("Comment skipped", 'Comment' not in mdf, True)
    os.unlink(mdf_path)

    # ── extract_mdf_timing ───────────────────────────────────────────────
    print("\nextract_mdf_timing:")
    timing = extract_mdf_timing(mdf)
    check("reference_date", timing['reference_date'], '20150101')
    check("reference_date_iso", timing['reference_date_iso'], '2015-01-01')
    check("tstart", timing['tstart'], 0.0)
    check("tstop", timing['tstop'], 43200.0)
    check("start_datetime", timing['start_datetime'], datetime(2015, 1, 1))
    check("end_datetime", timing['end_datetime'], datetime(2015, 1, 31))

    # Test with quoted Itdate
    mdf_quoted = {'Itdate': "'2017-06-15'", 'Tstart': '0.0', 'Tstop': '1440.0'}
    timing_q = extract_mdf_timing(mdf_quoted)
    check("quoted Itdate", timing_q['reference_date'], '20170615')

    # ── parse_bnd ────────────────────────────────────────────────────────
    print("\nparse_bnd:")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bnd', delete=False) as f:
        # 20-char fixed-width name field
        f.write("South               Z Z   1   1   1  50   0.0000000e+000\n")
        f.write("West                Z Z  50   1 100   1   0.0000000e+000\n")
        bnd_path = f.name

    bnd = parse_bnd(bnd_path)
    check("bnd count", len(bnd), 2)
    check("bnd[0] name", bnd[0]['bnd_name'], 'South')
    check("bnd[0] m_a", bnd[0]['m_a'], 1)
    check("bnd[0] n_a", bnd[0]['n_a'], 1)
    check("bnd[0] m_b", bnd[0]['m_b'], 1)
    check("bnd[0] n_b", bnd[0]['n_b'], 50)
    check("bnd[1] name", bnd[1]['bnd_name'], 'West')
    check("bnd[1] m_a", bnd[1]['m_a'], 50)
    os.unlink(bnd_path)

    # ── parse_grd ────────────────────────────────────────────────────────
    print("\nparse_grd:")
    # Create a minimal .grd file with 3x2 grid (M=3, N=2)
    grd_content = """Coordinate System = Cartesian
         0         0         0
 Missing Value =   0.0000000e+000
         3         2
 0 0 0
ETA=  1     100.0  200.0  300.0
ETA=  2     110.0  210.0  310.0
ETA=  1    1000.0 2000.0 3000.0
ETA=  2    1100.0 2100.0 3100.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_content)
        grd_path = f.name

    grd = parse_grd(grd_path)
    check("grd m", grd['m'], 3)
    check("grd n", grd['n'], 2)
    check("grd crs", grd['crs_type'], 'Cartesian')
    check("grd x count", len(grd['x_values']), 2)
    check("grd y count", len(grd['y_values']), 2)
    # First ETA x values: ['1', '100.0', '200.0', '300.0']
    # Element [0] is the ETA number, [1..M] are coordinates
    check("grd x[0][1]", grd['x_values'][0][1], '100.0')
    check("grd y[1][2]", grd['y_values'][1][2], '2100.0')
    os.unlink(grd_path)

    # Test Spherical
    grd_sph_content = """Coordinate System = Spherical
         0         0         0
 Missing Value =   0.0000000e+000
         2         2
 0 0 0
ETA=  1     6.500  7.500
ETA=  2     6.500  7.500
ETA=  1    53.000 53.000
ETA=  2    54.000 54.000
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_sph_content)
        grd_sph_path = f.name

    grd_sph = parse_grd(grd_sph_path)
    check("grd spherical crs", grd_sph['crs_type'], 'Spherical')
    os.unlink(grd_sph_path)

    # ── extract_boundary_coordinates ─────────────────────────────────────
    print("\nextract_boundary_coordinates:")

    # Test 1: legacy zero-missing fallback (missing value declared as 0.0)
    grd_zero_content = """Coordinate System = Cartesian
         0         0         0
 Missing Value =   0.0000000e+000
         4         3
 0 0 0
ETA=  1     100.0  200.0  300.0    0.0
ETA=  2     110.0  210.0  310.0    0.0
ETA=  3     120.0  220.0  320.0    0.0
ETA=  1    1000.0 2000.0 3000.0    0.0
ETA=  2    1100.0 2100.0 3100.0    0.0
ETA=  3    1200.0 2200.0 3200.0    0.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_zero_content)
        grd_zero_path = f.name

    bnd_zero = [
        {'bnd_name': 'North', 'm_a': 4, 'n_a': 1, 'm_b': 4, 'n_b': 3,
         'type1': 'Z', 'type2': 'Z', 'reflection': 0.0},
    ]
    result = extract_boundary_coordinates(grd_zero_path, bnd_zero)
    check("zero-miss: x_a != 0", result[0]['x_a'] != 0.0, True)
    check("zero-miss: y_a != 0", result[0]['y_a'] != 0.0, True)
    check("zero-miss: x_a from m=3", result[0]['x_a'], 300.0)
    check("zero-miss: y_a from m=3", result[0]['y_a'], 3000.0)
    check("zero-miss: x_b from m=3", result[0]['x_b'], 320.0)
    check("zero-miss: y_b from m=3", result[0]['y_b'], 3200.0)
    os.unlink(grd_zero_path)

    # Test 2: declared -999.999 missing value at edge
    grd_neg_content = """Coordinate System = Cartesian
         0         0         0
 Missing Value =  -9.99999000000000024E+002
         4         3
 0 0 0
ETA=  1     100.0  200.0  300.0 -9.99999000000000024E+002
ETA=  2     110.0  210.0  310.0 -9.99999000000000024E+002
ETA=  3     120.0  220.0  320.0 -9.99999000000000024E+002
ETA=  1    1000.0 2000.0 3000.0 -9.99999000000000024E+002
ETA=  2    1100.0 2100.0 3100.0 -9.99999000000000024E+002
ETA=  3    1200.0 2200.0 3200.0 -9.99999000000000024E+002
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_neg_content)
        grd_neg_path = f.name

    bnd_neg = [
        {'bnd_name': 'North', 'm_a': 4, 'n_a': 1, 'm_b': 4, 'n_b': 3,
         'type1': 'Z', 'type2': 'Z', 'reflection': 0.0},
    ]
    result = extract_boundary_coordinates(grd_neg_path, bnd_neg)
    check("neg-miss: x_a fallback", result[0]['x_a'], 300.0)
    check("neg-miss: y_a fallback", result[0]['y_a'], 3000.0)
    check("neg-miss: x_b fallback", result[0]['x_b'], 320.0)
    check("neg-miss: y_b fallback", result[0]['y_b'], 3200.0)
    os.unlink(grd_neg_path)

    # Test 3: pair-safe — x and y come from the same cell
    grd_pair_content = """Coordinate System = Cartesian
         0         0         0
 Missing Value =  -9.99900000000000000E+002
         3         2
 0 0 0
ETA=  1     100.0  200.0 -999.0
ETA=  2     110.0  210.0 -999.0
ETA=  1    1000.0 2000.0 -999.0
ETA=  2    1100.0 2100.0 -999.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_pair_content)
        grd_pair_path = f.name

    bnd_pair = [
        {'bnd_name': 'East', 'm_a': 3, 'n_a': 1, 'm_b': 3, 'n_b': 2,
         'type1': 'Z', 'type2': 'Z', 'reflection': 0.0},
    ]
    result = extract_boundary_coordinates(grd_pair_path, bnd_pair)
    check("pair: x_a=200, y_a=2000 (same cell m=2)", result[0]['x_a'], 200.0)
    check("pair: y_a=2000", result[0]['y_a'], 2000.0)
    check("pair: x_b=210, y_b=2100 (same cell m=2)", result[0]['x_b'], 210.0)
    check("pair: y_b=2100", result[0]['y_b'], 2100.0)
    os.unlink(grd_pair_path)

    # Test 4: max+1 clamp still works
    grd_clamp_content = """Coordinate System = Cartesian
         0         0         0
 Missing Value =   0.0000000e+000
         3         2
 0 0 0
ETA=  1     100.0  200.0  300.0
ETA=  2     110.0  210.0  310.0
ETA=  1    1000.0 2000.0 3000.0
ETA=  2    1100.0 2100.0 3100.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_clamp_content)
        grd_clamp_path = f.name

    bnd_clamp = [
        {'bnd_name': 'West', 'm_a': 1, 'n_a': 1, 'm_b': 1, 'n_b': 3,
         'type1': 'Z', 'type2': 'Z', 'reflection': 0.0},
    ]
    result = extract_boundary_coordinates(grd_clamp_path, bnd_clamp)
    check("clamp: n_b=3 clamped to 2, x_b", result[0]['x_b'], 110.0)
    check("clamp: n_b=3 clamped to 2, y_b", result[0]['y_b'], 1100.0)
    os.unlink(grd_clamp_path)

    # Test 5: no valid fallback raises ValueError
    grd_allbad_content = """Coordinate System = Cartesian
         0         0         0
 Missing Value =  -9.99900000000000000E+002
         2         2
 0 0 0
ETA=  1    -999.0 -999.0
ETA=  2    -999.0 -999.0
ETA=  1    -999.0 -999.0
ETA=  2    -999.0 -999.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.grd', delete=False) as f:
        f.write(grd_allbad_content)
        grd_allbad_path = f.name

    bnd_allbad = [
        {'bnd_name': 'Bad', 'm_a': 1, 'n_a': 1, 'm_b': 2, 'n_b': 2,
         'type1': 'Z', 'type2': 'Z', 'reflection': 0.0},
    ]
    try:
        extract_boundary_coordinates(grd_allbad_path, bnd_allbad)
        check("allbad: should raise ValueError", False, True)
    except ValueError as e:
        check("allbad: raises ValueError", 'No valid grid coordinates' in str(e), True)
    os.unlink(grd_allbad_path)

    # ── parse_mdu ────────────────────────────────────────────────────────
    print("\nparse_mdu:")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.mdu', delete=False) as f:
        f.write("[General]\n")
        f.write("Program   = D-Flow FM\n")
        f.write("Version   = 1.2.100\n")
        f.write("\n")
        f.write("[time]\n")
        f.write("RefDate   = 20170101\n")
        f.write("Tunit     = S\n")
        f.write("TStart    = 0.0\n")
        f.write("TStop     = 86400.0\n")
        f.write("DtUser    = 300.0  # timestep in seconds\n")
        mdu_path = f.name

    mdu = parse_mdu(mdu_path)
    check("mdu sections", sorted(mdu.keys()), ['General', 'time'])
    check("mdu Program", mdu['General']['Program'], 'D-Flow FM')
    check("mdu RefDate", mdu['time']['RefDate'], '20170101')
    check("mdu DtUser (inline comment stripped)", mdu['time']['DtUser'], '300.0')
    os.unlink(mdu_path)

    # ── extract_mdu_timing ───────────────────────────────────────────────
    print("\nextract_mdu_timing:")
    mdu_timing = extract_mdu_timing(mdu)
    check("mdu ref_date", mdu_timing['reference_date'], '20170101')
    check("mdu ref_date_iso", mdu_timing['reference_date_iso'], '2017-01-01')
    check("mdu tstart (minutes)", mdu_timing['tstart'], 0.0)
    check("mdu tstop (minutes)", mdu_timing['tstop'], 1440.0)  # 86400s = 1440min
    check("mdu start_dt", mdu_timing['start_datetime'], datetime(2017, 1, 1))
    check("mdu end_dt", mdu_timing['end_datetime'], datetime(2017, 1, 2))

    # ── parse_pli ────────────────────────────────────────────────────────
    print("\nparse_pli:")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pli', delete=False) as f:
        f.write("BND_1\n")
        f.write("2    2\n")
        f.write("6.50000    54.00000\n")
        f.write("6.50000    54.50000\n")
        f.write("BND_2\n")
        f.write("3    2\n")
        f.write("7.00000    53.00000\n")
        f.write("7.50000    53.50000\n")
        f.write("8.00000    54.00000\n")
        pli_path = f.name

    pli = parse_pli(pli_path)
    check("pli count", len(pli), 2)
    check("pli[0] name", pli[0]['name'], 'BND_1')
    check("pli[0] points", len(pli[0]['points']), 2)
    check("pli[0] point0", pli[0]['points'][0], (6.5, 54.0))
    check("pli[1] name", pli[1]['name'], 'BND_2')
    check("pli[1] points", len(pli[1]['points']), 3)
    check("pli[1] point2", pli[1]['points'][2], (8.0, 54.0))
    os.unlink(pli_path)

    # ── parse_mdw ────────────────────────────────────────────────────────
    print("\nparse_mdw:")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.mdw', delete=False) as f:
        f.write("[General]\n")
        f.write("   Program = Delft3D-WAVE\n")
        f.write("[Boundary]\n")
        f.write("   Name                 = BND_1\n")
        f.write("   Definition           = xy-coordinates\n")
        f.write("   StartCoordX          = 1.0000000e+005\n")
        f.write("[Boundary]\n")
        f.write("   Name                 = BND_2\n")
        f.write("   Definition           = xy-coordinates\n")
        mdw_path = f.name

    mdw = parse_mdw(mdw_path)
    check("mdw boundary count", len(mdw['boundaries']), 2)
    check("mdw bnd[0] Name", mdw['boundaries'][0]['Name'], 'BND_1')
    check("mdw bnd[1] Name", mdw['boundaries'][1]['Name'], 'BND_2')
    check("mdw header_end_line > 0", mdw['header_end_line'] > 0, True)
    os.unlink(mdw_path)

    # ── write_boundary_csv ───────────────────────────────────────────────
    print("\nwrite_boundary_csv:")
    bnd_data = [
        {'bnd_name': 'BND_1', 'x_a': 100.0, 'y_a': 200.0, 'x_b': 110.0, 'y_b': 210.0},
        {'bnd_name': 'BND_2', 'x_a': 300.0, 'y_a': 400.0, 'x_b': 310.0, 'y_b': 410.0},
    ]
    csv_path = tempfile.mktemp(suffix='.csv')
    write_boundary_csv(bnd_data, csv_path)
    with open(csv_path, 'r') as f:
        csv_lines = f.readlines()
    check("csv line count", len(csv_lines), 4)
    check("csv line 0", csv_lines[0].strip(), 'BND_1_a,100.0,200.0')
    check("csv line 1", csv_lines[1].strip(), 'BND_1_b,110.0,210.0')
    os.unlink(csv_path)

    # ── parse_validation_points ──────────────────────────────────────────
    print("\nparse_validation_points (Spherical):")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("# name lat lon\n")
        f.write("name lat lon\n")  # header row, should be skipped (non-numeric col 2)
        f.write("P1 54.0 8.0\n")
        f.write("P2, 53.5, 7.5\n")
        vp_path = f.name

    vp = parse_validation_points(vp_path, crs_type='Spherical')
    check("vp count", len(vp), 2)
    check("vp[0] name", vp[0]['name'], 'P1')
    check("vp[0] lat", vp[0]['lat'], 54.0)
    check("vp[0] lon", vp[0]['lon'], 8.0)
    check("vp[1] name (csv)", vp[1]['name'], 'P2')
    check("vp[1] lat (csv)", vp[1]['lat'], 53.5)
    os.unlink(vp_path)

    print("\nparse_validation_points (Cartesian):")
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write("P1 5935270.0 566463.0\n")  # northing, easting (Hamburg ish)
        vp_path2 = f.name

    vp_cart = parse_validation_points(
        vp_path2, crs_type='Cartesian', utm_zone=32, hemisphere='N'
    )
    check("vp_cart count", len(vp_cart), 1)
    check("vp_cart[0] name", vp_cart[0]['name'], 'P1')
    # lat ≈ 53.55, lon ≈ 10.0 — within tolerance
    lat_ok = abs(vp_cart[0]['lat'] - 53.55) < 0.1
    lon_ok = abs(vp_cart[0]['lon'] - 10.0) < 0.1
    check("vp_cart[0] lat ≈ 53.55", lat_ok, True)
    check("vp_cart[0] lon ≈ 10.0", lon_ok, True)
    os.unlink(vp_path2)

    # ── parse_bct_sections ───────────────────────────────────────────────
    print("\nparse_bct_sections:")
    from writer_delft3d4 import write_bct, write_bct_currents
    bct_round_path = tempfile.mktemp(suffix='.bct')
    bnd_rt = {
        'South': {'wl_a': [0.5, 0.6], 'wl_b': [0.4, 0.5]},
        'West':  {'wl_a': [1.0, 1.1], 'wl_b': [0.9, 1.0]},
    }
    write_bct(bct_round_path, bnd_rt, '20150101', 0.0, 20.0, 20.0)
    bct_parsed = parse_bct_sections(bct_round_path)
    check("bct ref_date", bct_parsed['reference_date'], '20150101')
    check("bct section count", len(bct_parsed['sections']), 2)
    check("bct sec0 name", bct_parsed['sections'][0]['boundary_name'], 'South')
    check("bct sec1 name", bct_parsed['sections'][1]['boundary_name'], 'West')
    check("bct sec0 n_data_cols", bct_parsed['sections'][0]['n_data_columns'], 2)
    check("bct sec0 unit[0]", bct_parsed['sections'][0]['parameters'][0]['unit'], '[m]')
    os.unlink(bct_round_path)

    # Currents .bct → 4 data columns
    bct_cur_path = tempfile.mktemp(suffix='.bct')
    cur_rt = {'South': {'u_a':[0.1,0.2], 'v_a':[-0.1,-0.2],
                         'u_b':[0.15,0.25], 'v_b':[-0.05,-0.15]}}
    write_bct_currents(bct_cur_path, cur_rt, '20150101', 0.0, 20.0, 20.0)
    cur_parsed = parse_bct_sections(bct_cur_path)
    check("currents n_data_cols", cur_parsed['sections'][0]['n_data_columns'], 4)
    os.unlink(bct_cur_path)

    # ── parse_bcc_sections ───────────────────────────────────────────────
    print("\nparse_bcc_sections:")
    from writer_delft3d4 import write_bcc_salinity
    bcc_path = tempfile.mktemp(suffix='.bcc')
    sal_rt = {
        'South': {'sal_a':[35.0,35.1], 'sal_b':[34.9,35.0]},
        'West':  {'sal_a':[34.5,34.6], 'sal_b':[34.4,34.5]},
    }
    write_bcc_salinity(bcc_path, sal_rt, '20150101', 0.0, 20.0, 20.0)
    bcc_parsed = parse_bcc_sections(bcc_path)
    check("bcc ref_date", bcc_parsed['reference_date'], '20150101')
    check("bcc section count", len(bcc_parsed['sections']), 2)
    check("bcc sec0 unit[0]", bcc_parsed['sections'][0]['parameters'][0]['unit'], '[ppt]')
    os.unlink(bcc_path)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
