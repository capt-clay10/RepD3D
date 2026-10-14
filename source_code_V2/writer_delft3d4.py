# -*- coding: utf-8 -*-
"""
writer_delft3d4.py — Write Delft3D-4 boundary condition files.

Writes:
    .bct — Water level / current time series (table-based ASCII)
    .bcc — Sediment concentration / salinity time series
    .bcw — Wave time series (SWAN boundary format)

The file formats are extremely sensitive to whitespace, column alignment,
and number formatting. Every format detail is preserved from v1.

v1 source:
    .bct: bct_generator.py lines 347-394
    .bcw: bcw_generator.py lines 414-467

Usage:
    from writer_delft3d4 import write_bct, write_bcw
"""

from pathlib import Path
from typing import Dict, List, Any
import csv
import os
import logging

from formatting_utils import (
    float_to_sci, list_to_sci, pad_positive,
    float_to_fixed, list_to_fixed, pad_to_integral_width,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# .bct parameter records (manual A.2.15)
#
# The quantity label and the 'end A'/'end B' suffix share one 26-character
# quoted field. Keeping the labels here means a mixed water-level/Neumann file
# cannot drift into two different column layouts.
# ─────────────────────────────────────────────────────────────────────────────

_BCT_PARAM_TIME = (
    "parameter            'time                '                     unit '[min]'"
)
_BCT_PARAMS = {
    'wl': (
        "parameter            'water elevation (z)  end A'               unit '[m]'",
        "parameter            'water elevation (z)  end B'               unit '[m]'",
    ),
    # EXACT Delft3D-FLOW GUI spelling: 'neumann' + 9 spaces + '(n)' + 2 spaces + 'end A'.
    # The parser matches this fixed-format string; a hand-edited variant with
    # 10 spaces and 1 space crashed FLOW at initialisation with
    # "Parameter name neumann (n) end a in (BCT, BCC or DIS) file not expected"
    # (NESH-HPC, 2026-10-13). Verified against a GUI-generated test file.
    'neumann': (
        "parameter            'neumann         (n)  end A'               unit '[-]'",
        "parameter            'neumann         (n)  end B'               unit '[-]'",
    ),
}


# ─────────────────────────────────────────────────────────────────────────────
# .bct FILE WRITER — Water Level / Current Time Series
# ─────────────────────────────────────────────────────────────────────────────

def write_bct(
    output_path: str,
    boundary_data: Dict[str, Dict[str, list]],
    reference_date: str,
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write a Delft3D-4 .bct (boundary condition time series) file.

    Args:
        output_path:    Path for the output .bct file.
        boundary_data:  Dict keyed by boundary name (without _a/_b suffix).
                        Each value is a dict:
                        {
                            'wl_a': [float, ...],   # water levels at end A
                            'wl_b': [float, ...],   # water levels at end B
                        }
        reference_date: Reference date string for header, e.g. '20150101'.
        tstart:         Simulation start time in minutes since reference.
        tstop:          Simulation stop time in minutes since reference.
        timestep:       Output timestep in minutes.

    Returns:
        The output file path.

    File format (preserved from v1 bct_generator.py lines 357-393):
        - Header block per boundary section (11 lines)
        - Data rows: time, wl_a, wl_b in scientific notation
        - Line ending: \\r\\n (CRLF — Windows convention for Delft3D)
        - Positive values have leading space for alignment
    """
    output_path = str(output_path)

    # Remove existing file (v1 behaviour)
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    # Generate time array
    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    section_number = 1
    for bnd_name, data in boundary_data.items():
        wl_a_list = list_to_sci(data['wl_a'], precision=7, exp_digits=3)
        wl_b_list = list_to_sci(data['wl_b'], precision=7, exp_digits=3)

        # Ensure all arrays have the same length
        n = min(len(time_list), len(wl_a_list), len(wl_b_list))
        # records-in-table must equal the number of data rows actually written
        # (manual A.2.15: it is the compulsory block terminator). Deriving it
        # from n rather than the requested time array keeps header and body in
        # agreement even if a series is short.
        records_in_table = n

        header_lines = [
            f"table-name           'Boundary Section : {section_number}'",
            "contents             'Uniform             '",
            f"location             '{bnd_name:<20s}'",
            "time-function        'non-equidistant'",
            f"reference-time       {reference_date}",
            "time-unit            'minutes'",
            "interpolation        'linear'",
            _BCT_PARAM_TIME,
            _BCT_PARAMS['wl'][0],
            _BCT_PARAMS['wl'][1],
            f"records-in-table     {records_in_table}",
        ]

        with open(output_path, 'a', newline='') as f:
            for line in header_lines:
                f.write(line + '\r\n')

            writer = csv.writer(f)
            for i in range(n):
                time_val = pad_positive(time_list[i])
                wl_a = pad_positive(wl_a_list[i])
                wl_b = pad_positive(wl_b_list[i])
                row_str = f'{time_val} {wl_a} {wl_b}'
                writer.writerow([row_str])

        section_number += 1

    logger.info(f"Written .bct file: {output_path} ({section_number - 1} sections)")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .bct FILE WRITER — Neumann (alongshore water-level gradient, type N)
# ─────────────────────────────────────────────────────────────────────────────

def write_bct_neumann(
    output_path: str,
    boundary_data: Dict[str, Dict[str, list]],
    reference_date: str,
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write a Delft3D-4 .bct file for Neumann (water-level gradient) sections.

    Same table structure as write_bct (Delft3D-FLOW manual A.2.15), but the
    parameter records carry the Neumann quantity with dimensionless unit
    (manual §9.4.1 Table 9.4: gradient amplitudes are [-]). The matching
    .bnd records must be typed 'N' and carry NO reflection coefficient
    (manual A.2.10).

    Args:
        output_path:    Path for the output .bct file.
        boundary_data:  Dict keyed by boundary name:
                        {'g_a': [...], 'g_b': [...]} gradients at ends A/B.
        reference_date: Reference date string, e.g. '20180101'.
        tstart:         Simulation start time in minutes since reference.
        tstop:          Simulation stop time in minutes since reference.
        timestep:       Output timestep in minutes.
    """
    output_path = str(output_path)
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    section_number = 1
    for bnd_name, data in boundary_data.items():
        g_a_list = list_to_sci(data['g_a'], precision=7, exp_digits=3)
        g_b_list = list_to_sci(data['g_b'], precision=7, exp_digits=3)
        n = min(len(time_list), len(g_a_list), len(g_b_list))
        records_in_table = n  # match actual rows (manual A.2.15 terminator)

        header_lines = [
            f"table-name           'Boundary Section : {section_number}'",
            "contents             'Uniform             '",
            f"location             '{bnd_name:<20s}'",
            "time-function        'non-equidistant'",
            f"reference-time       {reference_date}",
            "time-unit            'minutes'",
            "interpolation        'linear'",
            _BCT_PARAM_TIME,
            _BCT_PARAMS['neumann'][0],
            _BCT_PARAMS['neumann'][1],
            f"records-in-table     {records_in_table}",
        ]

        with open(output_path, 'a', newline='') as f:
            for line in header_lines:
                f.write(line + '\r\n')
            writer = csv.writer(f)
            for i in range(n):
                time_val = pad_positive(time_list[i])
                g_a = pad_positive(g_a_list[i])
                g_b = pad_positive(g_b_list[i])
                writer.writerow([f'{time_val} {g_a} {g_b}'])

        section_number += 1

    logger.info(
        f"Written Neumann .bct file: {output_path} ({section_number - 1} sections)"
    )
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .bct FILE WRITER — Mixed water-level (Z) + Neumann (N) sections
# ─────────────────────────────────────────────────────────────────────────────

def write_bct_mixed(
    output_path: str,
    sections: List[Dict[str, Any]],
    reference_date: str,
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write one Delft3D-4 .bct holding both water-level (.bnd type Z) and
    Neumann gradient (.bnd type N) sections, in the caller's order.

    Delft3D takes a table's *meaning* from the .bnd type letter, not from the
    .bct parameter label, and binds tables to sections by the 'location' name.
    A model whose .bnd mixes Z and N boundaries therefore needs a single .bct
    containing every section, each table carrying the quantity matching its
    .bnd type. Writing water-level values into an N section would have them
    read as gradients — metres where ~1e-6 is expected.

    Section order must follow the .bnd: 'table-name Boundary Section : k' is
    positional.

    Args:
        output_path:    Path for the output .bct file.
        sections:       Ordered list of dicts, one per .bnd section:
                        {
                            'name': str,             # .bnd section name
                            'kind': 'wl' | 'neumann',
                            'a':    [float, ...],    # series at end A
                            'b':    [float, ...],    # series at end B
                        }
        reference_date: Reference date string, e.g. '20180101'.
        tstart:         Simulation start time in minutes since reference.
        tstop:          Simulation stop time in minutes since reference.
        timestep:       Output timestep in minutes.

    Returns:
        The output file path.
    """
    output_path = str(output_path)

    for idx, section in enumerate(sections):
        missing = {'name', 'kind', 'a', 'b'} - set(section)
        if missing:
            raise ValueError(
                f"Section {idx} is missing key(s) {sorted(missing)} "
                "(need 'name', 'kind', 'a', 'b')."
            )
        if section['kind'] not in _BCT_PARAMS:
            raise ValueError(
                f"Section '{section['name']}' has unknown kind "
                f"'{section['kind']}' (use 'wl' or 'neumann')."
            )

    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    kind_counts: Dict[str, int] = {}
    section_number = 1
    for section in sections:
        a_list = list_to_sci(section['a'], precision=7, exp_digits=3)
        b_list = list_to_sci(section['b'], precision=7, exp_digits=3)
        n = min(len(time_list), len(a_list), len(b_list))
        records_in_table = n  # match actual rows (manual A.2.15 terminator)
        param_a, param_b = _BCT_PARAMS[section['kind']]

        header_lines = [
            f"table-name           'Boundary Section : {section_number}'",
            "contents             'Uniform             '",
            f"location             '{section['name']:<20s}'",
            "time-function        'non-equidistant'",
            f"reference-time       {reference_date}",
            "time-unit            'minutes'",
            "interpolation        'linear'",
            _BCT_PARAM_TIME,
            param_a,
            param_b,
            f"records-in-table     {records_in_table}",
        ]

        with open(output_path, 'a', newline='') as f:
            for line in header_lines:
                f.write(line + '\r\n')
            writer = csv.writer(f)
            for i in range(n):
                time_val = pad_positive(time_list[i])
                val_a = pad_positive(a_list[i])
                val_b = pad_positive(b_list[i])
                writer.writerow([f'{time_val} {val_a} {val_b}'])

        kind_counts[section['kind']] = kind_counts.get(section['kind'], 0) + 1
        section_number += 1

    breakdown = ", ".join(f"{count} {kind}" for kind, count in sorted(kind_counts.items()))
    logger.info(
        f"Written mixed .bct file: {output_path} "
        f"({section_number - 1} sections: {breakdown})"
    )
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .bct FILE WRITER — Currents (velocity u/v at end A and end B)
# ─────────────────────────────────────────────────────────────────────────────

def write_bct_currents(
    output_path: str,
    boundary_data: Dict[str, Dict[str, list]],
    reference_date: str,
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write a Delft3D-4 .bct file containing velocity boundary conditions.

    Spec: bc_file_format_spec.md §1.3 — same .bct envelope as water level,
    but four data columns (u@A, v@A, u@B, v@B) instead of two.

    Args:
        output_path:    Path for the output .bct file.
        boundary_data:  Dict keyed by boundary name (no _a/_b suffix).
                        Each value is a dict:
                        {
                            'u_a': [...], 'v_a': [...],
                            'u_b': [...], 'v_b': [...],
                        }
        reference_date: 'YYYYMMDD' string matching .mdf Itdate.
        tstart, tstop:  Simulation window in minutes since reference.
        timestep:       Output timestep in minutes.
    """
    output_path = str(output_path)

    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    section_number = 1
    for bnd_name, data in boundary_data.items():
        u_a = list_to_sci(data['u_a'], precision=7, exp_digits=3)
        v_a = list_to_sci(data['v_a'], precision=7, exp_digits=3)
        u_b = list_to_sci(data['u_b'], precision=7, exp_digits=3)
        v_b = list_to_sci(data['v_b'], precision=7, exp_digits=3)

        n = min(len(time_list), len(u_a), len(v_a), len(u_b), len(v_b))
        records_in_table = n  # match actual rows (manual A.2.15 terminator)

        header_lines = [
            f"table-name           'Boundary Section : {section_number}'",
            "contents             'Uniform             '",
            f"location             '{bnd_name:<20s}'",
            "time-function        'non-equidistant'",
            f"reference-time       {reference_date}",
            "time-unit            'minutes'",
            "interpolation        'linear'",
            "parameter            'time                '                     unit '[min]'",
            "parameter            'velocity (u)  end A'                      unit '[m/s]'",
            "parameter            'velocity (v)  end A'                      unit '[m/s]'",
            "parameter            'velocity (u)  end B'                      unit '[m/s]'",
            "parameter            'velocity (v)  end B'                      unit '[m/s]'",
            f"records-in-table     {records_in_table}",
        ]

        with open(output_path, 'a', newline='') as f:
            for line in header_lines:
                f.write(line + '\r\n')

            writer = csv.writer(f)
            for i in range(n):
                t = pad_positive(time_list[i])
                ua = pad_positive(u_a[i])
                va = pad_positive(v_a[i])
                ub = pad_positive(u_b[i])
                vb = pad_positive(v_b[i])
                writer.writerow([f'{t} {ua} {va} {ub} {vb}'])

        section_number += 1

    logger.info(
        f"Written .bct (currents) file: {output_path} "
        f"({section_number - 1} sections)"
    )
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .bcc FILE WRITER — Sediment Concentration Boundary Conditions
# ─────────────────────────────────────────────────────────────────────────────

def write_bcc(
    output_path: str,
    boundary_data: Dict[str, Dict[str, Dict[str, list]]],
    fraction_names: List[str],
    reference_date: str,
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write a Delft3D-4 .bcc (sediment concentration) file.

    Block ordering is **boundary-major, fraction-minor**, matching both the
    manual (A.2.17 example 1: each boundary section is followed by all of its
    constituents before the next section) and a Delft3D-NESTHD-generated
    reference .bcc: every sediment fraction of a boundary is emitted before
    moving to the next boundary, and `table-name 'Boundary Section : N'` is
    keyed to the boundary (N is constant across that boundary's fractions and
    increments per boundary). Each parameter is labelled
    '<fraction> end A uniform' / '<fraction> end B uniform' to pair with the
    'Uniform' contents.

    Args:
        output_path:   Output .bcc path.
        boundary_data: Nested dict
                       {
                           fraction_name: {
                               boundary_name: {
                                   'conc_a': [...], 'conc_b': [...],
                               },
                           },
                       }
                       Fraction key order MUST match the .sed fraction order;
                       boundary key order (within each fraction) MUST match the
                       .bnd open-boundary section order.
        fraction_names: Ordered list of fraction names (the .sed order); the
                        inner loop over each boundary's fractions.
        reference_date: 'YYYYMMDD'.
        tstart, tstop:  Simulation window in minutes since reference.
        timestep:       Output timestep in minutes.
    """
    output_path = str(output_path)

    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    for fraction_name in fraction_names:
        if fraction_name not in boundary_data:
            raise ValueError(
                f"Fraction '{fraction_name}' missing from boundary_data"
            )

    # Boundary section order is taken from the first fraction; every fraction
    # must cover the same boundary sections (one .sed run over one .bnd).
    boundary_order = list(boundary_data[fraction_names[0]].keys())
    for fraction_name in fraction_names[1:]:
        if list(boundary_data[fraction_name].keys()) != boundary_order:
            raise ValueError(
                f"Fraction '{fraction_name}' has a different set/order of "
                f"boundary sections than '{fraction_names[0]}'; all fractions "
                f"must share the same .bnd sections in the same order."
            )

    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    for section_number, bnd_name in enumerate(boundary_order, start=1):
        for fraction_name in fraction_names:
            data = boundary_data[fraction_name][bnd_name]
            conc_a = list_to_sci(data['conc_a'], precision=7, exp_digits=3)
            conc_b = list_to_sci(data['conc_b'], precision=7, exp_digits=3)
            n = min(len(time_list), len(conc_a), len(conc_b))
            records_in_table = n  # match actual rows (manual A.2.17 terminator)

            label_a = f"{fraction_name:<20s} end A uniform"
            label_b = f"{fraction_name:<20s} end B uniform"

            header_lines = [
                f"table-name           'Boundary Section : {section_number}'",
                "contents             'Uniform   '",
                f"location             '{bnd_name:<20s}'",
                "time-function        'non-equidistant'",
                f"reference-time       {reference_date}",
                "time-unit            'minutes'",
                "interpolation        'linear'",
                "parameter            'time                '  unit '[min]'",
                f"parameter            '{label_a}'       unit '[kg/m3]'",
                f"parameter            '{label_b}'       unit '[kg/m3]'",
                f"records-in-table     {records_in_table}",
            ]

            with open(output_path, 'a', newline='') as f:
                for line in header_lines:
                    f.write(line + '\r\n')
                writer = csv.writer(f)
                for i in range(n):
                    t = pad_positive(time_list[i])
                    ca = pad_positive(conc_a[i])
                    cb = pad_positive(conc_b[i])
                    writer.writerow([f'{t} {ca} {cb}'])

    logger.info(
        f"Written .bcc file: {output_path} "
        f"({len(fraction_names)} fractions)"
    )
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .bcc FILE WRITER — Salinity Boundary
# ─────────────────────────────────────────────────────────────────────────────

def write_bcc_salinity(
    output_path: str,
    boundary_data: Dict[str, Dict[str, list]],
    reference_date: str,
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write a Delft3D-4 .bcc file containing salinity boundary conditions.

    Spec: bc_file_format_spec.md §3.5. Same envelope as .bct, with
    'Salinity end A/B uniform' parameter labels. Single variable (no
    fractions), one table block per boundary section.

    Args:
        output_path:   Output .bcc path.
        boundary_data: {boundary_name: {'sal_a': [...], 'sal_b': [...]}}
        reference_date: 'YYYYMMDD'.
        tstart, tstop:  Simulation window in minutes since reference.
        timestep:       Output timestep in minutes.
    """
    output_path = str(output_path)

    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    section_number = 1
    for bnd_name, data in boundary_data.items():
        sal_a = list_to_sci(data['sal_a'], precision=7, exp_digits=3)
        sal_b = list_to_sci(data['sal_b'], precision=7, exp_digits=3)
        n = min(len(time_list), len(sal_a), len(sal_b))
        records_in_table = n  # match actual rows (manual A.2.17 terminator)

        header_lines = [
            f"table-name           'Boundary Section : {section_number}'",
            "contents             'Uniform             '",
            f"location             '{bnd_name:<20s}'",
            "time-function        'non-equidistant'",
            f"reference-time       {reference_date}",
            "time-unit            'minutes'",
            "interpolation        'linear'",
            "parameter            'time                '                     unit '[min]'",
            "parameter            'Salinity end A uniform'                   unit '[ppt]'",
            "parameter            'Salinity end B uniform'                   unit '[ppt]'",
            f"records-in-table     {records_in_table}",
        ]

        with open(output_path, 'a', newline='') as f:
            for line in header_lines:
                f.write(line + '\r\n')
            writer = csv.writer(f)
            for i in range(n):
                t = pad_positive(time_list[i])
                sa = pad_positive(sal_a[i])
                sb = pad_positive(sal_b[i])
                writer.writerow([f'{t} {sa} {sb}'])

        section_number += 1

    logger.info(f"Written .bcc (salinity) file: {output_path}")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# HEADER-REPLAY WRITER — Reuse an existing .bct/.bcc file's headers verbatim
# ─────────────────────────────────────────────────────────────────────────────

def write_bc_with_header_replay(
    output_path: str,
    parsed_sections: List[Dict[str, Any]],
    section_values: List[List[list]],
    tstart: float,
    tstop: float,
    timestep: float,
) -> str:
    """
    Write a .bct or .bcc file by replaying parsed section headers verbatim
    and substituting fresh time-series data.

    Used by the synthetic BC module when the user supplies an existing file
    of the same type (instead of .mdf+.bnd+.grd). The headers — parameter
    labels, units, location strings — are preserved exactly; only the data
    rows and `records-in-table` are refreshed for the new timeline.

    Args:
        output_path:     Target path.
        parsed_sections: List of section dicts from parse_bct_sections() /
                         parse_bcc_sections(). Each has 'header_lines' and
                         'n_data_columns'.
        section_values:  Parallel list to parsed_sections. Each element is a
                         list of value-column lists, length n_data_columns,
                         each inner list length == number of timesteps.
        tstart, tstop:   Time window in minutes since the reference date.
        timestep:        Output timestep in minutes.

    Returns:
        output_path.
    """
    if len(parsed_sections) != len(section_values):
        raise ValueError(
            f"parsed_sections ({len(parsed_sections)}) and section_values "
            f"({len(section_values)}) must have matching length"
        )

    output_path = str(output_path)
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    import numpy as np
    float_range = np.arange(tstart, tstop + timestep, timestep).tolist()
    time_list = list_to_sci(float_range, precision=7, exp_digits=3)

    with open(output_path, 'a', newline='') as f:
        for section, columns in zip(parsed_sections, section_values):
            expected = section['n_data_columns']
            if len(columns) != expected:
                raise ValueError(
                    f"Section '{section['boundary_name']}' expects "
                    f"{expected} data columns; got {len(columns)}"
                )
            sci_columns = [
                list_to_sci(col, precision=7, exp_digits=3) for col in columns
            ]

            # records-in-table must equal the rows actually written; clamp to
            # the shortest column so header and body agree and no short column
            # is indexed out of range (manual A.2.15 terminator).
            records_in_table = min(
                [len(time_list)] + [len(c) for c in sci_columns]
            )

            for line in section['header_lines']:
                if line.lstrip().startswith('records-in-table'):
                    f.write(f"records-in-table     {records_in_table}\r\n")
                else:
                    f.write(line + '\r\n')

            writer = csv.writer(f)
            for i in range(records_in_table):
                parts = [pad_positive(time_list[i])]
                for col in sci_columns:
                    parts.append(pad_positive(col[i]))
                writer.writerow([' '.join(parts)])

    logger.info(
        f"Written BC file with replayed header: {output_path} "
        f"({len(parsed_sections)} sections)"
    )
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .bcw FILE WRITER — SWAN Wave Boundary
# ─────────────────────────────────────────────────────────────────────────────

def write_bcw(
    output_path: str,
    boundary_data: Dict[str, Dict[str, list]],
    reference_date: str,
) -> str:
    """
    Write a Delft3D-4 .bcw (SWAN wave boundary) file.

    Args:
        output_path:    Path for the output .bcw file.
        boundary_data:  Dict keyed by boundary name (single point per boundary).
                        Each value is a dict:
                        {
                            'time':   [float, ...],  # time in minutes
                            'hs':     [float, ...],  # significant wave height (m)
                            'tp':     [float, ...],  # peak period (s)
                            'dir':    [float, ...],  # direction (degrees nautical)
                            'spread': [float, ...],  # directional spreading (degrees)
                        }
        reference_date: Reference date string, e.g. '20150101'.

    Returns:
        The output file path.

    File format (preserved from v1 bcw_generator.py lines 421-467):
        - Header per boundary (10 lines — NO table-name, NO records-in-table)
        - Data rows: time, Hs, Tp, Dir, DirSpreading in fixed-point format
        - Line ending: \\n (LF — NOT CRLF! Different from .bct!)
        - Columns padded to specific integral widths
    """
    output_path = str(output_path)

    # Remove existing file
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    for bnd_name, data in boundary_data.items():
        time_vals = list_to_fixed(data['time'], decimal_digits=2)
        hs_vals = list_to_fixed(data['hs'], decimal_digits=4)
        tp_vals = list_to_fixed(data['tp'], decimal_digits=4)
        dir_vals = list_to_fixed(data['dir'], decimal_digits=4)
        spread_vals = list_to_fixed(data['spread'], decimal_digits=4)

        n = min(len(time_vals), len(hs_vals), len(tp_vals),
                len(dir_vals), len(spread_vals))

        # Calculate integral width for time column
        max_time = data['time'][-1] if data['time'] else 0
        time_integral_width = len(str(int(max_time))) + 1  # +1 for padding

        header_lines = [
            f"location             '{bnd_name:<23s}'",
            "time-function        'non-equidistant'",
            f"reference-time       {reference_date}",
            "time-unit            'minutes'",
            "interpolation        'linear'",
            "parameter            'time                  '                   unit '[min]'",
            "parameter            'WaveHeight'                               unit '[m]'",
            "parameter            'Period'                                   unit '[s]'",
            "parameter            'Direction'                                unit '[N^o]'",
            "parameter            'DirSpreading'                             unit '[deg]'",
        ]

        with open(output_path, 'a', newline='') as f:
            for line in header_lines:
                f.write(line + '\n')

            writer = csv.writer(f, lineterminator='\n')
            for i in range(n):
                t = pad_to_integral_width(time_vals[i], time_integral_width)
                hs = pad_to_integral_width(hs_vals[i], 2)
                tp = pad_to_integral_width(tp_vals[i], 2)
                dr = pad_to_integral_width(dir_vals[i], 4)
                sp = pad_to_integral_width(spread_vals[i], 3)
                row_str = f'{t} {hs} {tp} {dr} {sp}'
                writer.writerow([row_str])

    logger.info(f"Written .bcw file: {output_path}")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile
    import numpy as np

    print("=== writer_delft3d4.py self-test ===\n")
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

    # ── write_bct ────────────────────────────────────────────────────────
    print("write_bct:")
    bct_path = tempfile.mktemp(suffix='.bct')
    bnd_data = {
        'South': {
            'wl_a': [0.5, 0.6, 0.7],
            'wl_b': [0.4, 0.5, 0.6],
        },
        'West': {
            'wl_a': [1.0, 1.1, 1.2],
            'wl_b': [0.9, 1.0, 1.1],
        },
    }
    result_path = write_bct(bct_path, bnd_data, '20150101', 0.0, 40.0, 20.0)
    check("returns path", result_path, bct_path)

    with open(bct_path, 'r') as f:
        content = f.read()
        lines = content.split('\n')

    check("has Section 1", "Boundary Section : 1" in content, True)
    check("has Section 2", "Boundary Section : 2" in content, True)
    check("has reference time", "20150101" in content, True)
    check("has records-in-table", "records-in-table" in content, True)
    check("has time parameter", "'time" in content, True)
    check("has water elevation", "water elevation" in content, True)

    # Check data format — first data line should be scientific notation
    # Find first line after header
    for line in lines:
        stripped = line.strip()
        if stripped and 'e+' in stripped.lower() and 'parameter' not in stripped.lower():
            # This is a data line
            parts = stripped.split()
            check("data has 3 columns", len(parts), 3)
            check("time starts with sci notation",
                  'e+' in parts[0].lower() or 'e-' in parts[0].lower(), True)
            break

    os.unlink(bct_path)

    # ── write_bcw ────────────────────────────────────────────────────────
    print("\nwrite_bcw:")
    bcw_path = tempfile.mktemp(suffix='.bcw')
    wave_data = {
        'South': {
            'time': [0.0, 20.0, 40.0],
            'hs': [1.5, 1.6, 1.4],
            'tp': [8.0, 8.5, 7.8],
            'dir': [270.0, 265.5, 280.0],
            'spread': [30.0, 28.5, 32.0],
        },
    }
    result_bcw = write_bcw(bcw_path, wave_data, '20150101')
    check("returns path", result_bcw, bcw_path)

    with open(bcw_path, 'r') as f:
        bcw_content = f.read()
        bcw_lines = bcw_content.split('\n')

    check("has location", "location" in bcw_content, True)
    check("has WaveHeight param", "WaveHeight" in bcw_content, True)
    check("has Direction param", "Direction" in bcw_content, True)
    check("has DirSpreading param", "DirSpreading" in bcw_content, True)
    check("NO table-name (bcw format)", "table-name" not in bcw_content, True)
    check("NO records-in-table (bcw format)", "records-in-table" not in bcw_content, True)

    # Check data format — fixed-point
    for line in bcw_lines:
        stripped = line.strip()
        if stripped and stripped[0].isdigit():
            parts = stripped.split()
            check("bcw data has 5 columns", len(parts), 5)
            check("bcw time is fixed-point", '.' in parts[0], True)
            check("bcw no sci notation", 'e' not in parts[1].lower(), True)
            break

    # Check line endings are LF not CRLF
    with open(bcw_path, 'rb') as f:
        raw = f.read()
    check("bcw uses LF (not CRLF)", b'\r\n' not in raw, True)

    os.unlink(bcw_path)

    # ── write_bct_mixed ──────────────────────────────────────────────────
    print("\nwrite_bct_mixed:")
    mixed_path = tempfile.mktemp(suffix='.bct')
    write_bct_mixed(
        mixed_path,
        [
            {'name': 'North01', 'kind': 'neumann',
             'a': [1.0e-5, 2.0e-5, 3.0e-5], 'b': [0.0, -1.0e-5, -2.0e-5]},
            {'name': 'West_1', 'kind': 'wl',
             'a': [0.5, 0.6, 0.7], 'b': [0.4, 0.5, 0.6]},
            {'name': 'South01', 'kind': 'neumann',
             'a': [3.0e-6, 4.0e-6, 5.0e-6], 'b': [3.0e-6, 4.0e-6, 5.0e-6]},
        ],
        '20180101', 0.0, 40.0, 20.0,
    )
    with open(mixed_path, 'rb') as f:
        mixed_raw = f.read()
    mixed_text = mixed_raw.decode('ascii')
    check("mixed: 3 sections", mixed_text.count("table-name"), 3)
    check("mixed: 2 neumann tables", mixed_text.count("'neumann         (n)"), 4)
    # Regression: the Delft3D-FLOW parser accepts exactly the GUI spelling
    # (9 spaces, 2 spaces); 10/1 crashed FLOW (see _BCT_PARAMS). Count, don't eyeball.
    import re as _re
    _label = _re.search(r"'neumann( +)\(n\)( +)end A'", mixed_text)
    check("neumann label: 9 spaces before (n)", len(_label.group(1)) if _label else None, 9)
    check("neumann label: 2 spaces before end A", len(_label.group(2)) if _label else None, 2)
    check("neumann label aligned with water-level label (same record length)",
          len(_BCT_PARAMS['neumann'][0]), len(_BCT_PARAMS['wl'][0]))
    check("mixed: 1 water-level table", mixed_text.count("'water elevation (z)"), 2)
    check("mixed: CRLF only",
          mixed_raw.count(b'\n') == mixed_raw.count(b'\r\n'), True)

    # Section numbering must follow the caller's (.bnd) order, not group by kind
    order = [
        mixed_text.index("'North01"),
        mixed_text.index("'West_1"),
        mixed_text.index("'South01"),
    ]
    check("mixed: preserves .bnd order", order == sorted(order), True)
    check("mixed: North01 is Section 1",
          mixed_text.index("Boundary Section : 1") < mixed_text.index("'North01"), True)
    check("mixed: South01 is Section 3",
          mixed_text.index("Boundary Section : 3") < mixed_text.index("'South01"), True)

    # Regression guard: mixed output must be byte-identical to the dedicated
    # writers when every section is of one kind.
    only_wl_mixed = tempfile.mktemp(suffix='.bct')
    only_wl_plain = tempfile.mktemp(suffix='.bct')
    write_bct_mixed(
        only_wl_mixed,
        [{'name': 'South', 'kind': 'wl', 'a': [0.5, 0.6, 0.7], 'b': [0.4, 0.5, 0.6]},
         {'name': 'West', 'kind': 'wl', 'a': [1.0, 1.1, 1.2], 'b': [0.9, 1.0, 1.1]}],
        '20150101', 0.0, 40.0, 20.0,
    )
    write_bct(only_wl_plain, bnd_data, '20150101', 0.0, 40.0, 20.0)
    check("mixed(all wl) == write_bct byte-for-byte",
          open(only_wl_mixed, 'rb').read() == open(only_wl_plain, 'rb').read(), True)

    only_n_mixed = tempfile.mktemp(suffix='.bct')
    only_n_plain = tempfile.mktemp(suffix='.bct')
    write_bct_mixed(
        only_n_mixed,
        [{'name': 'WestLat', 'kind': 'neumann',
          'a': [1.0e-5, 2.0e-5, 3.0e-5], 'b': [0.0, -1.0e-5, -2.0e-5]}],
        '20180101', 0.0, 40.0, 20.0,
    )
    write_bct_neumann(
        only_n_plain,
        {'WestLat': {'g_a': [1.0e-5, 2.0e-5, 3.0e-5], 'g_b': [0.0, -1.0e-5, -2.0e-5]}},
        '20180101', 0.0, 40.0, 20.0,
    )
    check("mixed(all neumann) == write_bct_neumann byte-for-byte",
          open(only_n_mixed, 'rb').read() == open(only_n_plain, 'rb').read(), True)

    for bad, label in (
        ({'name': 'X', 'kind': 'gradient', 'a': [1.0], 'b': [1.0]}, "unknown kind"),
        ({'name': 'X', 'kind': 'wl', 'a': [1.0]}, "missing key"),
    ):
        try:
            write_bct_mixed(tempfile.mktemp(suffix='.bct'), [bad],
                            '20180101', 0.0, 0.0, 20.0)
            check(f"mixed: rejects {label}", "no error", "ValueError")
        except ValueError:
            check(f"mixed: rejects {label}", "ValueError", "ValueError")

    for path in (mixed_path, only_wl_mixed, only_wl_plain, only_n_mixed, only_n_plain):
        os.remove(path)

    # ── write_bct_currents ───────────────────────────────────────────────
    print("\nwrite_bct_currents:")
    bct_cur_path = tempfile.mktemp(suffix='.bct')
    cur_data = {
        'South': {
            'u_a': [0.1, 0.2, 0.3], 'v_a': [-0.1, -0.2, -0.3],
            'u_b': [0.15, 0.25, 0.35], 'v_b': [-0.05, -0.15, -0.25],
        },
    }
    write_bct_currents(bct_cur_path, cur_data, '20150101', 0.0, 40.0, 20.0)
    with open(bct_cur_path) as f:
        cur_content = f.read()
    check("currents has u end A", "velocity (u)  end A" in cur_content, True)
    check("currents has v end B", "velocity (v)  end B" in cur_content, True)
    check("currents has 5 cols",
          all(len(l.split()) == 5
              for l in cur_content.split('\n')
              if l.strip() and 'e+' in l.lower() and 'parameter' not in l), True)
    os.unlink(bct_cur_path)

    # ── write_bcc ────────────────────────────────────────────────────────
    print("\nwrite_bcc:")
    bcc_path = tempfile.mktemp(suffix='.bcc')
    sed_data = {
        'sand': {
            'South': {'conc_a': [0.001, 0.002], 'conc_b': [0.0011, 0.0021]},
            'West':  {'conc_a': [0.003, 0.004], 'conc_b': [0.0031, 0.0041]},
        },
        'mud': {
            'South': {'conc_a': [0.005, 0.006], 'conc_b': [0.0051, 0.0061]},
            'West':  {'conc_a': [0.007, 0.008], 'conc_b': [0.0071, 0.0081]},
        },
    }
    write_bcc(bcc_path, sed_data, ['sand', 'mud'], '20150101', 0.0, 20.0, 20.0)
    with open(bcc_path) as f:
        bcc_content = f.read()
    check("bcc labels fraction 'sand' end A uniform",
          "'sand                 end A uniform'" in bcc_content, True)
    check("bcc labels fraction 'mud' end B uniform",
          "'mud                  end B uniform'" in bcc_content, True)
    check("bcc has no bare duplicated concentration label",
          "'concentration       '" in bcc_content, False)
    check("bcc has 4 sections (2 boundaries x 2 fractions)",
          bcc_content.count("table-name"), 4)
    check("bcc has kg/m3 unit", "[kg/m3]" in bcc_content, True)
    # Boundary-major ordering: section number is keyed to the boundary, so the
    # two blocks of the first boundary both read 'Boundary Section : 1'.
    check("bcc section number keyed to boundary",
          bcc_content.count("Boundary Section : 1"), 2)
    check("bcc has a second boundary section",
          bcc_content.count("Boundary Section : 2"), 2)
    # First two table blocks must be the SAME boundary (South) with the two
    # fractions in .sed order, proving boundary-major / fraction-minor.
    blocks = [b for b in bcc_content.split("table-name") if b.strip()]
    check("block 1 is South/sand",
          "'South" in blocks[0] and "'sand" in blocks[0], True)
    check("block 2 is South/mud (same boundary, next fraction)",
          "'South" in blocks[1] and "'mud" in blocks[1], True)
    check("block 3 moves to West/sand",
          "'West" in blocks[2] and "'sand" in blocks[2], True)
    with open(bcc_path, 'rb') as f:
        check("bcc uses CRLF", b'\r\n' in f.read(), True)
    os.unlink(bcc_path)

    # ── write_bcc_salinity ──────────────────────────────────────────────
    print("\nwrite_bcc_salinity:")
    bcc_sal_path = tempfile.mktemp(suffix='.bcc')
    sal_data = {
        'South': {'sal_a': [35.0, 35.1, 35.2], 'sal_b': [34.9, 35.0, 35.1]},
        'West':  {'sal_a': [34.5, 34.6, 34.7], 'sal_b': [34.4, 34.5, 34.6]},
    }
    write_bcc_salinity(bcc_sal_path, sal_data, '20150101', 0.0, 40.0, 20.0)
    with open(bcc_sal_path) as f:
        sal_content = f.read()
    check("salinity has Salinity param",
          "Salinity end A uniform" in sal_content, True)
    check("salinity has end B param",
          "Salinity end B uniform" in sal_content, True)
    check("salinity has ppt unit", "[ppt]" in sal_content, True)
    check("salinity has 2 sections (2 boundaries)",
          sal_content.count("table-name"), 2)
    with open(bcc_sal_path, 'rb') as f:
        check("salinity bcc uses CRLF", b'\r\n' in f.read(), True)
    os.unlink(bcc_sal_path)

    # ── write_bc_with_header_replay ─────────────────────────────────────
    print("\nwrite_bc_with_header_replay:")
    from file_parsers import parse_bct_sections
    src_path = tempfile.mktemp(suffix='.bct')
    write_bct(src_path, {
        'South': {'wl_a':[0.5,0.6], 'wl_b':[0.4,0.5]},
        'West':  {'wl_a':[1.0,1.1], 'wl_b':[0.9,1.0]},
    }, '20150101', 0.0, 20.0, 20.0)
    parsed = parse_bct_sections(src_path)
    replay_path = tempfile.mktemp(suffix='.bct')
    new_values = [
        [[9.0, 9.1, 9.2], [8.0, 8.1, 8.2]],   # South: wl_a, wl_b
        [[7.0, 7.1, 7.2], [6.0, 6.1, 6.2]],   # West:  wl_a, wl_b
    ]
    write_bc_with_header_replay(
        replay_path, parsed['sections'], new_values, 0.0, 40.0, 20.0,
    )
    with open(replay_path) as f:
        rp_content = f.read()
    check("replay has South location", "'South" in rp_content, True)
    check("replay has West location", "'West" in rp_content, True)
    check("replay has updated records-in-table",
          "records-in-table     3" in rp_content, True)
    with open(replay_path, 'rb') as f:
        check("replay uses CRLF", b'\r\n' in f.read(), True)
    os.unlink(src_path)
    os.unlink(replay_path)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
