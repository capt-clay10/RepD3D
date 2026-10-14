# -*- coding: utf-8 -*-
"""
backend_sea_level.py — Add sea level change to existing .bct files.

Reads an existing .bct file, modifies the water level columns (end A and
end B), and writes a new file with '_wl_changed.bct' suffix.

Two modes:
    1. Linear:   Linearly interpolate from 0 to target SLR across the simulation
    2. Constant: Add a fixed offset to all water levels

v1 source: sea_level_change.py (166 lines).
v2 improvements:
    - Uses formatting_utils instead of 6th inline copy of sci notation code
    - Uses logger instead of no logging
    - Returns output path
    - Accepts logger + progress_callback for GUI integration
    - Cleaner file structure parsing (no hardcoded line offsets)

Usage:
    from backend_sea_level import add_sea_level_change
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Callable, Optional
import os
import threading
import logging

from formatting_utils import list_to_sci


class TaskCancelledError(Exception):
    """Raised when the backend observes a cancel_event set by the GUI."""
    pass


def _check_cancel(cancel_event: Optional[threading.Event]) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise TaskCancelledError("Cancelled by user")


def add_sea_level_change(
    bct_file: str,
    sea_level_change_m: float,
    constant: bool = False,
    output_path: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Add sea level change to an existing .bct file.

    Args:
        bct_file:           Path to the input .bct file.
        sea_level_change_m: Sea level change in metres (positive = rise).
        constant:           If True, add constant offset. If False, linear 0→target.
        output_path:        Output file path. Defaults to '{name}_wl_changed.bct'.
        logger:             Logger for status messages.
        progress_callback:  Progress callback (0.0–1.0).

    Returns:
        {'output_file': str, 'n_sections': int, 'mode': str}
    """
    if logger is None:
        logger = logging.getLogger('repd3d.sea_level')

    bct_file = str(bct_file)

    # Default output path
    if output_path is None:
        stem = Path(bct_file).stem
        parent = Path(bct_file).parent
        output_path = str(parent / f'{stem}_wl_changed.bct')

    mode = 'constant' if constant else 'linear'
    logger.info(f"Adding {sea_level_change_m:.3f}m sea level change ({mode}) to {Path(bct_file).name}")

    # ── Parse the .bct file structure ────────────────────────────────────
    # .bct files have a repeating structure:
    #   11 header lines per section
    #   N data lines per section (from records-in-table)
    #   This repeats for each boundary section.

    # Read the entire file
    with open(bct_file, 'r') as f:
        all_lines = f.readlines()

    # Find records-in-table from the first section header (line 11, 0-indexed 10)
    data_lines_per_section = _extract_records_in_table(all_lines)
    header_lines_per_section = 11

    # Count total sections
    n_sections = sum(1 for line in all_lines if 'interpolation' in line.lower())

    logger.info(f"Found {n_sections} boundary sections, {data_lines_per_section} records each")

    # Generate the sea level change array
    if constant:
        wl_change = np.full(data_lines_per_section, sea_level_change_m)
    else:
        # Linear interpolation from 0 to target
        wl_change = np.linspace(0, sea_level_change_m, data_lines_per_section)

    # ── Write the new file ───────────────────────────────────────────────
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    block_size = header_lines_per_section + data_lines_per_section

    try:
        with open(output_path, 'w') as f_out:
            for section_idx in range(n_sections):
                _check_cancel(cancel_event)

                block_start = section_idx * block_size

                # Write header lines unchanged
                header_start = block_start
                header_end = block_start + header_lines_per_section
                for line in all_lines[header_start:header_end]:
                    f_out.write(line)

                # Read and modify data lines
                data_start = header_end
                data_end = data_start + data_lines_per_section

                for line_idx, line in enumerate(all_lines[data_start:data_end]):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        time_val = parts[0]
                        wl_a = float(parts[1]) + wl_change[line_idx]
                        wl_b = float(parts[2]) + wl_change[line_idx]

                        # Convert back to scientific notation
                        wl_a_str = list_to_sci([wl_a], precision=7, exp_digits=3)[0]
                        wl_b_str = list_to_sci([wl_b], precision=7, exp_digits=3)[0]

                        # Pad positive values
                        if not wl_a_str.startswith('-'):
                            wl_a_str = ' ' + wl_a_str
                        if not wl_b_str.startswith('-'):
                            wl_b_str = ' ' + wl_b_str
                        if not time_val.startswith('-'):
                            time_val = ' ' + time_val

                        f_out.write(f' {time_val} {wl_a_str} {wl_b_str}\n')
                    else:
                        f_out.write(line)

                progress_callback((section_idx + 1) / n_sections)
    except TaskCancelledError:
        logger.info("Sea level change cancelled — removing partial output")
        try:
            os.remove(output_path)
        except OSError:
            pass
        raise

    logger.info(f"Written: {output_path}")

    return {
        'output_file': output_path,
        'n_sections': n_sections,
        'mode': mode,
    }


def _extract_records_in_table(lines: list) -> int:
    """Extract records-in-table value from .bct header lines."""
    for line in lines[:15]:  # Check first 15 lines
        if 'records-in-table' in line.lower():
            parts = line.split()
            return int(parts[-1])
    raise ValueError("Could not find 'records-in-table' in .bct file header")


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile
    from writer_delft3d4 import write_bct

    print("=== backend_sea_level.py self-test ===\n")
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

    # Create a test .bct file first
    bct_path = tempfile.mktemp(suffix='.bct')
    bnd_data = {
        'South': {'wl_a': [0.5, 0.6, 0.7, 0.8, 0.9], 'wl_b': [0.4, 0.5, 0.6, 0.7, 0.8]},
        'West': {'wl_a': [1.0, 1.1, 1.2, 1.3, 1.4], 'wl_b': [0.9, 1.0, 1.1, 1.2, 1.3]},
    }
    write_bct(bct_path, bnd_data, '20150101', 0.0, 80.0, 20.0)

    # Test constant SLR
    print("Constant sea level change (+0.5m):")
    result = add_sea_level_change(bct_path, 0.5, constant=True)
    check("returns output_file", 'output_file' in result, True)
    check("mode is constant", result['mode'], 'constant')
    check("n_sections", result['n_sections'], 2)

    # Verify the output values increased
    with open(result['output_file'], 'r') as f:
        out_lines = f.readlines()
    # Find first data line (after 11 header lines)
    first_data = out_lines[11].strip().split()
    original_wl_a = 0.5
    modified_wl_a = float(first_data[1])
    check("WL increased", abs(modified_wl_a - (original_wl_a + 0.5)) < 0.001, True)

    os.unlink(result['output_file'])

    # Test linear SLR
    print("\nLinear sea level change (0→1.0m):")
    result_lin = add_sea_level_change(bct_path, 1.0, constant=False)
    check("mode is linear", result_lin['mode'], 'linear')

    with open(result_lin['output_file'], 'r') as f:
        lin_lines = f.readlines()
    first_data_lin = lin_lines[11].strip().split()
    last_data_lin = lin_lines[15].strip().split()
    first_change = float(first_data_lin[1]) - 0.5
    last_change = float(last_data_lin[1]) - 0.9
    check("first timestep ≈ no change", abs(first_change) < 0.01, True)
    check("last timestep ≈ +1.0m", abs(last_change - 1.0) < 0.01, True)

    os.unlink(result_lin['output_file'])
    os.unlink(bct_path)

    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
