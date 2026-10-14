# -*- coding: utf-8 -*-
"""
writer_mdw.py — Generate/update Delft3D WAVE .mdw definition files.

Reads an existing .mdw file (which must have at least one [Boundary]
placeholder), replaces the [Boundary] blocks with actual boundary
coordinates extracted from the .grd + .bnd files.

v1 source: mdw_writer.py (110 lines).
v2 improvements:
    - Uses formatting_utils instead of inline sci notation copy
    - Uses file_parsers.parse_mdw for reading
    - Proper file handle management

Usage:
    from writer_mdw import write_mdw
"""

from pathlib import Path
from typing import List, Dict, Any
import os
import logging

from formatting_utils import float_to_sci

logger = logging.getLogger(__name__)


def write_mdw(
    template_mdw_path: str,
    boundary_coordinates: List[Dict[str, Any]],
    output_path: str,
) -> str:
    """
    Write an updated .mdw file with boundary coordinates.

    Takes the header portion of an existing .mdw file (everything before
    the first [Boundary]) and appends new [Boundary] blocks with the
    correct coordinates from the grid/boundary file extraction.

    Args:
        template_mdw_path:    Path to existing .mdw file (used as template).
        boundary_coordinates: List of boundary dicts, each with:
                              {
                                  'name':  str,   # boundary name
                                  'x_a':   float, # start X coordinate
                                  'y_a':   float, # start Y coordinate
                                  'x_b':   float, # end X coordinate
                                  'y_b':   float, # end Y coordinate
                              }
        output_path:          Path for the output .mdw file.

    Returns:
        The output file path.

    v1 source: mdw_writer.py lines 60-109.
    v1 used inline convert_flt_to_sci_not (7th copy of that function).
    v2 uses formatting_utils.float_to_sci.
    """
    from file_parsers import parse_mdw

    output_path = str(output_path)

    # Remove existing output
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    # Parse template to get header lines
    mdw_data = parse_mdw(template_mdw_path)
    header_lines = mdw_data['header_lines']

    # Write header
    with open(output_path, 'w') as f:
        for line in header_lines:
            f.write(line)  # header_lines already include newlines

    # Write boundary blocks
    with open(output_path, 'a') as f:
        for bnd in boundary_coordinates:
            name = bnd['name']
            x_a = float_to_sci(bnd['x_a'], precision=7, exp_digits=3)
            y_a = float_to_sci(bnd['y_a'], precision=7, exp_digits=3)
            x_b = float_to_sci(bnd['x_b'], precision=7, exp_digits=3)
            y_b = float_to_sci(bnd['y_b'], precision=7, exp_digits=3)

            bnd_block = [
                '[Boundary]',
                f'   Name                 = {name:<30s}',
                '   Definition           = xy-coordinates               ',
                f'   StartCoordX          =  {x_a}              ',
                f'   EndCoordX            =  {x_b}              ',
                f'   StartCoordY          =  {y_a}              ',
                f'   EndCoordY            =  {y_b}              ',
                '   SpectrumSpec         = parametric                   ',
                '   SpShapeType          = jonswap                      ',
                '   PeriodType           = peak                         ',
                '   DirSpreadType        = degrees                        ',
                '   PeakEnhanceFac       =  3.3000000e+000              ',
                '   GaussSpread          =  9.9999998e-003              ',
            ]

            for line in bnd_block:
                f.write(line + '\n')

    logger.info(f"Written .mdw file: {output_path} ({len(boundary_coordinates)} boundaries)")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile

    print("=== writer_mdw.py self-test ===\n")
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

    # Create a template .mdw
    template_path = tempfile.mktemp(suffix='.mdw')
    with open(template_path, 'w') as f:
        f.write("[General]\n")
        f.write("   Program = Delft3D-WAVE\n")
        f.write("   Version = 3.04\n")
        f.write("[Boundary]\n")
        f.write("   Name = placeholder\n")

    # Write new .mdw
    output_path = tempfile.mktemp(suffix='_new.mdw')
    boundaries = [
        {'name': 'BND_1', 'x_a': 100000.0, 'y_a': 5900000.0, 'x_b': 110000.0, 'y_b': 5910000.0},
        {'name': 'BND_2', 'x_a': 200000.0, 'y_a': 5950000.0, 'x_b': 210000.0, 'y_b': 5960000.0},
    ]

    result = write_mdw(template_path, boundaries, output_path)
    check("returns path", result, output_path)

    with open(output_path, 'r') as f:
        content = f.read()

    # Verify header preserved
    check("has [General]", "[General]" in content, True)
    check("has Program", "Delft3D-WAVE" in content, True)

    # Verify old placeholder removed
    check("no placeholder", "placeholder" not in content, True)

    # Verify new boundaries
    check("has 2 [Boundary] blocks", content.count("[Boundary]"), 2)
    check("has BND_1", "BND_1" in content, True)
    check("has BND_2", "BND_2" in content, True)
    check("has sci notation coords", "e+" in content, True)
    check("has jonswap", "jonswap" in content, True)
    check("has PeakEnhanceFac", "3.3000000e+000" in content, True)

    os.unlink(template_path)
    os.unlink(output_path)

    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
