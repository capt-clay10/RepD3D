# -*- coding: utf-8 -*-
"""
writer_delft3dfm.py — Write Delft3D-FM boundary condition files.

Writes:
    .bc      — Universal INI-style forcing file (water level, velocity, etc.)
    .pli     — Polyline boundary location file
    _bnd.ext — External forcings specification linking .pli to .bc

NEW in v2 — these formats did not exist in v1.
Specifications from D-Flow FM User Manual and handoff document Section 4.

Usage:
    from writer_delft3dfm import write_bc, write_pli, write_ext
"""

from pathlib import Path
from typing import Dict, List, Any, Optional
import os
import logging

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# .bc FILE WRITER — Universal Forcing File
# ─────────────────────────────────────────────────────────────────────────────

def write_bc(
    output_path: str,
    boundary_data: Dict[str, Dict[str, list]],
    quantity: str,
    unit: str,
    reference_datetime: str,
    time_unit: str = "minutes",
    vector_name: Optional[str] = None,
) -> str:
    """
    Write a Delft3D-FM .bc forcing file.

    Args:
        output_path:        Path for the output .bc file.
        boundary_data:      Dict keyed by boundary point name.
                            Each value is a dict:
                            {
                                'time': [float, ...],     # time values
                                'values': [float, ...],   # data values
                            }
                            For vector quantities (velocity), use:
                            {
                                'time': [float, ...],
                                'values_x': [float, ...],
                                'values_y': [float, ...],
                            }
        quantity:            Forcing quantity name, e.g.:
                            'waterlevelbnd', 'velocitybnd',
                            'dischargebnd', 'salinitybnd'
        unit:               Unit string, e.g. 'm', 'm/s', 'm3/s'
        reference_datetime: Reference datetime string,
                            e.g. '2017-01-01 00:00:00'
        time_unit:          Time unit, e.g. 'minutes', 'seconds', 'hours'

    Per-block overrides:
        Each ``data`` entry may carry optional keys that override the
        file-wide arguments for that one forcing block — needed when a single
        .bc collects several quantities at the same support point (e.g. one
        sediment fraction per block, all named after the same polyline point):
            'name':     value for the INI ``Name`` key (default: the dict key)
            'quantity': scalar quantity name, OR a ``(qx, qy)`` pair for a
                        vector block (see below)
            'unit':     unit string, OR a ``(unit_x, unit_y)`` pair
        For a vector forcing, provide 'values_x'/'values_y'. The two data
        columns are labelled from the ``quantity`` tuple when one is given
        (e.g. ``('ux', 'uy')`` for ``uxuyadvectionvelocitybnd``); a bare
        string falls back to the legacy ``<quantity>_x`` / ``<quantity>_y``.

    Vector declaration:
        A vector boundary condition needs a ``Vector`` line binding the
        aggregate quantity to its component columns (manual C.5 example,
        p.487): ``Vector = uxuyadvectionvelocitybnd:ux,uy``. Pass the aggregate
        name via ``vector_name``; the component names come from the vector
        ``quantity`` tuple.

    Returns:
        The output file path.

    File format (D-Flow FM User Manual, Appendix C.5):
        [General]
        fileVersion                     = 1.01
        fileType                        = boundConds

        [forcing]
        Name                            = BND_1_0001
        Function                        = timeseries
        timeInterpolation               = linear
        Quantity                        = time
        Unit                            = minutes since 2017-01-01 00:00:00
        Quantity                        = waterlevelbnd
        Unit                            = m
        0.0    1.234
        20.0   1.245
        ...
    """
    output_path = str(output_path)

    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass

    with open(output_path, 'w', newline='\n') as f:
        # Mandatory [General] block (manual Table C.1) — without it D-Flow FM
        # does not recognise the file as a boundary-conditions forcing file.
        f.write("[General]\n")
        f.write("fileVersion                     = 1.01\n")
        f.write("fileType                        = boundConds\n")
        f.write("\n")

        for point_name, data in boundary_data.items():
            block_name = data.get('name', point_name)
            block_quantity = data.get('quantity', quantity)
            block_unit = data.get('unit', unit)

            f.write("[forcing]\n")
            f.write(f"Name                            = {block_name}\n")
            f.write("Function                        = timeseries\n")
            f.write("timeInterpolation               = linear\n")
            f.write("Quantity                        = time\n")
            f.write(f"Unit                            = {time_unit} since {reference_datetime}\n")

            is_vector = 'values_x' in data and 'values_y' in data

            if is_vector:
                if isinstance(block_quantity, (list, tuple)):
                    qx, qy = block_quantity
                else:
                    qx, qy = f"{block_quantity}_x", f"{block_quantity}_y"
                if isinstance(block_unit, (list, tuple)):
                    unit_x, unit_y = block_unit
                else:
                    unit_x = unit_y = block_unit
                block_vector_name = data.get('vector_name', vector_name)
                if block_vector_name:
                    f.write(
                        f"Vector                          = "
                        f"{block_vector_name}:{qx},{qy}\n"
                    )
                f.write(f"Quantity                        = {qx}\n")
                f.write(f"Unit                            = {unit_x}\n")
                f.write(f"Quantity                        = {qy}\n")
                f.write(f"Unit                            = {unit_y}\n")
            else:
                f.write(f"Quantity                        = {block_quantity}\n")
                f.write(f"Unit                            = {block_unit}\n")

            time_vals = data['time']

            if is_vector:
                vals_x = data['values_x']
                vals_y = data['values_y']
                n = min(len(time_vals), len(vals_x), len(vals_y))
                for i in range(n):
                    f.write(f"{time_vals[i]:.6f}    {vals_x[i]:.6f}    {vals_y[i]:.6f}\n")
            else:
                values = data['values']
                n = min(len(time_vals), len(values))
                for i in range(n):
                    f.write(f"{time_vals[i]:.6f}    {values[i]:.6f}\n")

            f.write("\n")  # Blank line between forcing blocks

    logger.info(f"Written .bc file: {output_path} ({len(boundary_data)} points)")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# .pli FILE WRITER — Polyline Boundary Locations
# ─────────────────────────────────────────────────────────────────────────────

def write_pli(
    output_path: str,
    boundaries: List[Dict[str, Any]],
) -> str:
    """
    Write a Delft3D-FM .pli polyline file.

    Args:
        output_path: Path for the output .pli file.
        boundaries:  List of boundary dicts:
                     [
                         {
                             'name': 'BND_1',
                             'points': [(lon1, lat1), (lon2, lat2), ...],
                         },
                         ...
                     ]

    Returns:
        The output file path.

    File format:
        BND_1
        2    2
        6.50000    54.00000
        6.50000    54.50000
    """
    output_path = str(output_path)

    with open(output_path, 'w', newline='\n') as f:
        for bnd in boundaries:
            name = bnd['name']
            points = bnd['points']
            n_points = len(points)
            n_dims = 2  # Always 2D (x, y)

            f.write(f"{name}\n")
            f.write(f"{n_points}    {n_dims}\n")
            for x, y in points:
                f.write(f"{x:.5f}    {y:.5f}\n")

    logger.info(f"Written .pli file: {output_path} ({len(boundaries)} boundaries)")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# _bnd.ext FILE WRITER — External Forcings Specification
# ─────────────────────────────────────────────────────────────────────────────

def write_ext(
    output_path: str,
    quantity: str,
    pli_file: str,
    bc_file: str,
) -> str:
    """
    Write a Delft3D-FM _bnd.ext external forcings file.

    Args:
        output_path: Path for the output .ext file.
        quantity:     Forcing quantity, e.g. 'waterlevelbnd'.
        pli_file:     Filename of the .pli file (relative, not absolute).
        bc_file:      Filename of the .bc file (relative, not absolute).

    Returns:
        The output file path.

    File format (D-Flow FM User Manual, Appendix C.6.2 — new-style external
    forcings; the [General] block is mandatory):
        [General]
        fileVersion       = 2.02
        fileType          = extForce

        [Boundary]
        quantity          = waterlevelbnd
        locationFile      = boundaries.pli
        forcingFile       = boundaries.bc
    """
    output_path = str(output_path)

    with open(output_path, 'w', newline='\n') as f:
        f.write("[General]\n")
        f.write("fileVersion       = 2.02\n")
        f.write("fileType          = extForce\n")
        f.write("\n")
        f.write("[Boundary]\n")
        f.write(f"quantity          = {quantity}\n")
        f.write(f"locationFile      = {pli_file}\n")
        f.write(f"forcingFile       = {bc_file}\n")

    logger.info(f"Written .ext file: {output_path}")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile

    print("=== writer_delft3dfm.py self-test ===\n")
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

    # ── write_bc (scalar) ────────────────────────────────────────────────
    print("write_bc (scalar — water level):")
    bc_path = tempfile.mktemp(suffix='.bc')
    bc_data = {
        'BND_1_0001': {'time': [0.0, 20.0, 40.0], 'values': [1.0, 1.1, 1.2]},
        'BND_1_0002': {'time': [0.0, 20.0, 40.0], 'values': [0.9, 1.0, 1.1]},
    }
    write_bc(bc_path, bc_data, 'waterlevelbnd', 'm', '2017-01-01 00:00:00')

    with open(bc_path, 'r') as f:
        content = f.read()
    check("has [General]", "[General]" in content, True)
    check("has fileVersion", "fileVersion                     = 1.01" in content, True)
    check("has fileType boundConds", "fileType                        = boundConds" in content, True)
    check("has [forcing]", "[forcing]" in content, True)
    check("has Name", "BND_1_0001" in content, True)
    check("has timeseries", "timeseries" in content, True)
    check("has timeInterpolation", "timeInterpolation" in content, True)
    check("no hyphenated Time-interpolation", "Time-interpolation" in content, False)
    check("has quantity", "waterlevelbnd" in content, True)
    check("has time unit", "minutes since 2017-01-01" in content, True)
    check("has 2 forcing blocks", content.count("[forcing]"), 2)
    os.unlink(bc_path)

    # ── write_bc (vector) ────────────────────────────────────────────────
    print("\nwrite_bc (vector — velocity):")
    bc_vec_path = tempfile.mktemp(suffix='.bc')
    bc_vec_data = {
        'BND_1_0001': {
            'time': [0.0, 20.0],
            'values_x': [0.5, 0.6],
            'values_y': [0.1, 0.2],
        },
    }
    write_bc(bc_vec_path, bc_vec_data, 'velocitybnd', 'm/s', '2017-01-01 00:00:00')

    with open(bc_vec_path, 'r') as f:
        vec_content = f.read()
    check("vector fallback x quantity", "velocitybnd_x" in vec_content, True)
    check("vector fallback y quantity", "velocitybnd_y" in vec_content, True)

    # Check data line has 3 columns
    for line in vec_content.split('\n'):
        stripped = line.strip()
        if stripped and stripped[0].isdigit():
            parts = stripped.split()
            check("vector data has 3 cols", len(parts), 3)
            break
    os.unlink(bc_vec_path)

    # ── write_bc (vector with explicit ux/uy column names + Vector line) ──
    print("\nwrite_bc (vector — uxuyadvectionvelocitybnd ux/uy):")
    bc_uxuy_path = tempfile.mktemp(suffix='.bc')
    write_bc(bc_uxuy_path, bc_vec_data, ('ux', 'uy'), 'm/s', '2017-01-01 00:00:00',
             vector_name='uxuyadvectionvelocitybnd')
    with open(bc_uxuy_path, 'r') as f:
        uxuy_content = f.read()
    check("uxuy has Vector declaration",
          "Vector                          = uxuyadvectionvelocitybnd:ux,uy" in uxuy_content, True)
    check("uxuy Vector precedes ux quantity",
          uxuy_content.index("Vector") < uxuy_content.index("Quantity                        = ux"), True)
    check("uxuy has ux column", "Quantity                        = ux" in uxuy_content, True)
    check("uxuy has uy column", "Quantity                        = uy" in uxuy_content, True)
    check("uxuy has no _x suffix", "velocitybnd_x" in uxuy_content, False)
    os.unlink(bc_uxuy_path)

    # Vector line only emitted when a vector_name is supplied
    bc_novec_path = tempfile.mktemp(suffix='.bc')
    write_bc(bc_novec_path, bc_vec_data, ('ux', 'uy'), 'm/s', '2017-01-01 00:00:00')
    with open(bc_novec_path, 'r') as f:
        novec_content = f.read()
    check("no Vector line without vector_name", "Vector" in novec_content, False)
    os.unlink(bc_novec_path)

    # ── write_bc (per-block quantity/name override — sediment fractions) ──
    print("\nwrite_bc (per-block quantity/name override):")
    bc_sed_path = tempfile.mktemp(suffix='.bc')
    sed_bc_data = {
        'BND_1_0001__sand': {
            'name': 'BND_1_0001', 'quantity': 'sedfracbndsand', 'unit': 'kg/m3',
            'time': [0.0, 20.0], 'values': [0.01, 0.02],
        },
        'BND_1_0001__mud': {
            'name': 'BND_1_0001', 'quantity': 'sedfracbndmud', 'unit': 'kg/m3',
            'time': [0.0, 20.0], 'values': [0.03, 0.04],
        },
    }
    write_bc(bc_sed_path, sed_bc_data, 'sedfracbnd', 'kg/m3', '2017-01-01 00:00:00')
    with open(bc_sed_path, 'r') as f:
        sed_content = f.read()
    check("both blocks share point Name", sed_content.count("Name                            = BND_1_0001"), 2)
    check("sand quantity present", "sedfracbndsand" in sed_content, True)
    check("mud quantity present", "sedfracbndmud" in sed_content, True)
    check("no fraction-prefixed Name leaked", "BND_1_0001__sand" in sed_content, False)
    os.unlink(bc_sed_path)

    # ── write_pli ────────────────────────────────────────────────────────
    print("\nwrite_pli:")
    pli_path = tempfile.mktemp(suffix='.pli')
    pli_data = [
        {'name': 'BND_1', 'points': [(6.5, 54.0), (6.5, 54.5)]},
        {'name': 'BND_2', 'points': [(7.0, 53.0), (7.5, 53.5), (8.0, 54.0)]},
    ]
    write_pli(pli_path, pli_data)

    with open(pli_path, 'r') as f:
        pli_content = f.read()
        pli_lines = [l for l in pli_content.split('\n') if l.strip()]
    check("pli has BND_1", "BND_1" in pli_content, True)
    check("pli has 2 2", "2    2" in pli_content, True)
    check("pli has 3 2", "3    2" in pli_content, True)
    check("pli has coordinate", "6.50000" in pli_content, True)
    os.unlink(pli_path)

    # ── write_ext ────────────────────────────────────────────────────────
    print("\nwrite_ext:")
    ext_path = tempfile.mktemp(suffix='.ext')
    write_ext(ext_path, 'waterlevelbnd', 'boundaries.pli', 'boundaries.bc')

    with open(ext_path, 'r') as f:
        ext_content = f.read()
    check("ext has [General]", "[General]" in ext_content, True)
    check("ext has fileType extForce", "fileType          = extForce" in ext_content, True)
    check("ext has [Boundary]", "[Boundary]" in ext_content, True)
    check("ext has quantity", "waterlevelbnd" in ext_content, True)
    check("ext has locationFile", "boundaries.pli" in ext_content, True)
    check("ext has forcingFile", "boundaries.bc" in ext_content, True)
    os.unlink(ext_path)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
