# -*- coding: utf-8 -*-
"""
backend_synthetic_bc.py — Generate hand-authored boundary condition files.

Given a small set of user-entered control dates and values, produce a full
equidistant-time Delft3D BC file by linear- or step-interpolating between
control points. Reuses the existing writer and parser primitives; no GUI
imports.

Entry point:
    generate_synthetic_bc(file_type, ...) -> dict

Supported file_type values:
    'waterlevel_bct'   -> writer_delft3d4.write_bct
    'currents_bct'     -> writer_delft3d4.write_bct_currents
    'salinity_bcc'     -> writer_delft3d4.write_bcc_salinity
    'sediments_bcc'    -> writer_delft3d4.write_bcc
    'waves_bcw'        -> writer_delft3d4.write_bcw (+ writer_mdw.write_mdw)
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Callable, Any, Tuple
import logging
import threading

import numpy as np

import file_parsers as fp
from writer_delft3d4 import (
    write_bct, write_bct_currents, write_bcc, write_bcc_salinity, write_bcw,
    write_bc_with_header_replay,
)
from writer_mdw import write_mdw
from writer_meteo import write_wnd_point


FILE_TYPES = [
    'waterlevel_bct',
    'currents_bct',
    'waves_bcw',
    'sediments_bcc',
    'salinity_bcc',
    'point_wnd',
]

# Column layout consumed by the GUI. Each tuple is (parameter_key, unit_label).
COLUMN_SPECS: Dict[str, List[Tuple[str, str]]] = {
    'waterlevel_bct': [('wl_a', '[m]'), ('wl_b', '[m]')],
    'currents_bct':   [('u_a',  '[m/s]'), ('v_a',  '[m/s]'),
                       ('u_b',  '[m/s]'), ('v_b',  '[m/s]')],
    'salinity_bcc':   [('sal_a', '[ppt]'),  ('sal_b', '[ppt]')],
    'sediments_bcc':  [('conc_a', '[kg/m3]'), ('conc_b', '[kg/m3]')],
    'waves_bcw':      [('hs', '[m]'), ('tp', '[s]'),
                       ('dir', '[deg]'), ('spread', '[deg]')],
    'point_wnd':      [('ws', '[m/s]'), ('wdir', '[deg]')],
}

WAVE_COL_ORDER = ('hs', 'tp', 'dir', 'spread')


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def _parse_dt(s: str) -> datetime:
    """Parse 'YYYY-MM-DD HH:MM:SS' or 'YYYY-MM-DD'."""
    s = s.strip()
    for fmt in ('%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H:%M', '%Y-%m-%d'):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"Unparseable datetime: {s!r}")


def _minutes_since(dt: datetime, ref: datetime) -> float:
    return (dt - ref).total_seconds() / 60.0


def _ref_from_yyyymmdd(ref_str: str) -> datetime:
    """Parse reference-date strings like '20150101' or '2015-01-01'."""
    s = ref_str.strip().strip("'").strip('#').replace('-', '')
    return datetime.strptime(s[:8], '%Y%m%d')


def _interp_series(
    timeline_min: np.ndarray,
    control_min: np.ndarray,
    control_vals: np.ndarray,
    mode: str,
) -> np.ndarray:
    """Interpolate a single column. mode ∈ {'linear', 'constant'}."""
    if mode == 'linear':
        return np.interp(timeline_min, control_min, control_vals)
    # 'constant' = step-hold: value at time t equals the most recent control
    # value at or before t.
    idx = np.searchsorted(control_min, timeline_min, side='right') - 1
    idx = np.clip(idx, 0, len(control_vals) - 1)
    return control_vals[idx]


def _validate_control_dates(
    dates: List[datetime], start: datetime, end: datetime
) -> None:
    if len(dates) < 2:
        raise ValueError("Need at least start and end control dates.")
    if dates[0] != start or dates[-1] != end:
        raise ValueError("First and last control dates must match start/end.")
    for earlier, later in zip(dates, dates[1:]):
        if later <= earlier:
            raise ValueError(
                f"Control dates must be strictly increasing; "
                f"{earlier} is not before {later}"
            )


def _preview_lines(path: str, first_n: int = 1, last_n: int = 5) -> List[str]:
    """Return first N + last M lines of a text file (no duplication)."""
    lines = Path(path).read_text().splitlines()
    if len(lines) <= first_n + last_n:
        return lines
    head = lines[:first_n]
    tail = lines[-last_n:]
    return head + ['...'] + tail


# ─────────────────────────────────────────────────────────────────────────────
# Context builders
# ─────────────────────────────────────────────────────────────────────────────

def _boundary_names_from_bnd(bnd_path: str) -> List[str]:
    return [b['bnd_name'] for b in fp.parse_bnd(bnd_path)]


def _bnd_pairs_from_bnd_grd(bnd_path: str, grd_path: str) -> List[Dict[str, Any]]:
    bnd_indices = fp.parse_bnd(bnd_path)
    return fp.extract_boundary_coordinates(grd_path, bnd_indices)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def generate_synthetic_bc(
    file_type: str,
    source_mode: str,
    start_date: str,
    end_date: str,
    timestep_min: float,
    control_dates: List[str],
    values: Dict[str, Any],
    interpolation: str,
    output_path: str,
    mdf_path: Optional[str] = None,
    mdw_path: Optional[str] = None,
    bnd_path: Optional[str] = None,
    grd_path: Optional[str] = None,
    existing_file_path: Optional[str] = None,
    fraction_names: Optional[List[str]] = None,
    write_companion_mdw: bool = True,
    mdw_output_path: Optional[str] = None,
    reference_datetime: Optional[str] = None,
    user_reference_date: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> Dict[str, Any]:
    """
    Generate a hand-authored Delft3D boundary condition file.

    Args:
        file_type: One of FILE_TYPES.
        source_mode: 'model_files' or 'existing_file'. Waves must be
            'model_files'.
        start_date, end_date, timestep_min: Output time window (datetime
            strings) and equidistant step in minutes.
        control_dates: Full list of control datetimes as strings, length
            2 + n_midpoints. First must equal start_date, last must equal
            end_date.
        values: Control-point values keyed in a shape that matches file_type:
            - waterlevel_bct:  {bnd_name: {'wl_a':[...], 'wl_b':[...]}}
            - currents_bct:    {bnd_name: {'u_a':[...], 'v_a':[...],
                                           'u_b':[...], 'v_b':[...]}}
            - salinity_bcc:    {bnd_name: {'sal_a':[...], 'sal_b':[...]}}
            - sediments_bcc:   {frac: {bnd_name: {'conc_a':[...], 'conc_b':[...]}}}
            - waves_bcw:       {bnd_name: {'hs':[...], 'tp':[...],
                                           'dir':[...], 'spread':[...]}}
            Each inner list has len(control_dates) entries.
        interpolation: 'linear' or 'constant' (step-hold).
        output_path: Target BC file path.
        mdf_path, bnd_path, grd_path: required for non-wave model_files mode.
        mdw_path, bnd_path, grd_path: required for waves.
        existing_file_path: required for non-wave existing_file mode.
        fraction_names: sediments_bcc only; fraction order for the writer.
        write_companion_mdw: waves only; also append boundary blocks to .mdw.
        mdw_output_path: waves only; defaults to <output_path>.mdw sibling.

    Returns:
        {'output_path': str, 'n_timesteps': int, 'preview_lines': [...],
         'mdw_path': Optional[str]}
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if file_type not in FILE_TYPES:
        raise ValueError(f"Unknown file_type '{file_type}'. Expected one of {FILE_TYPES}")
    if interpolation not in ('linear', 'constant'):
        raise ValueError(f"interpolation must be 'linear' or 'constant'; got {interpolation!r}")

    is_wave = file_type == 'waves_bcw'
    is_point_wnd = file_type == 'point_wnd'
    if is_wave and source_mode != 'model_files':
        raise ValueError("Waves (.bcw) only supports 'model_files' source mode.")
    if not is_point_wnd and source_mode not in ('model_files', 'existing_file'):
        raise ValueError(f"source_mode must be 'model_files' or 'existing_file'")

    # ── Parse dates ─────────────────────────────────────────────────────
    start_dt = _parse_dt(start_date)
    end_dt = _parse_dt(end_date)
    control_dt = [_parse_dt(s) for s in control_dates]
    _validate_control_dates(control_dt, start_dt, end_dt)
    if timestep_min <= 0:
        raise ValueError(f"timestep_min must be positive; got {timestep_min}")

    # ── Point .wnd: no boundaries, no model files — short-circuit ───────
    if is_point_wnd:
        if reference_datetime and reference_datetime.strip():
            ref_dt = _parse_dt(reference_datetime)
        else:
            ref_dt = start_dt.replace(hour=0, minute=0, second=0, microsecond=0)
        tstart = _minutes_since(start_dt, ref_dt)
        tstop = _minutes_since(end_dt, ref_dt)
        if tstop <= tstart:
            raise ValueError("end_date must be strictly after start_date")
        timeline = np.arange(tstart, tstop + timestep_min, timestep_min)
        control_min = np.array([_minutes_since(dt, ref_dt) for dt in control_dt])
        logger.info(
            f"Timeline: {len(timeline)} steps, "
            f"{tstart:.1f} → {tstop:.1f} min, dt={timestep_min}"
        )
        progress_callback(0.35)

        ws_ctrl = np.asarray(values['ws'], dtype=float)
        wd_ctrl = np.asarray(values['wdir'], dtype=float)
        if ws_ctrl.shape != control_min.shape or wd_ctrl.shape != control_min.shape:
            raise ValueError(
                f"Expected {len(control_min)} control values for ws/wdir"
            )
        ws = _interp_series(timeline, control_min, ws_ctrl, interpolation).tolist()
        wd = _interp_series(timeline, control_min, wd_ctrl, interpolation).tolist()

        progress_callback(0.7)
        write_wnd_point(output_path, timeline.tolist(), ws, wd)
        preview = _preview_lines(output_path)
        for line in preview:
            logger.info(line)
        progress_callback(1.0)
        logger.info(f"Done: {output_path}")
        return {
            'output_path': output_path,
            'n_timesteps': len(timeline),
            'preview_lines': preview,
            'mdw_path': None,
        }

    # ── Build context (reference_date, boundary list / section list) ────
    parsed_sections: Optional[List[Dict[str, Any]]] = None
    bnd_pairs: Optional[List[Dict[str, Any]]] = None

    if source_mode == 'existing_file':
        if not existing_file_path:
            raise ValueError("existing_file_path required for 'existing_file' mode")
        logger.info(f"Parsing existing file header: {Path(existing_file_path).name}")
        if file_type == 'waterlevel_bct' or file_type == 'currents_bct':
            parsed = fp.parse_bct_sections(existing_file_path)
        else:
            parsed = fp.parse_bcc_sections(existing_file_path)
        parsed_sections = parsed['sections']
        ref_str = parsed['reference_date']
        ref_dt = _ref_from_yyyymmdd(ref_str) if ref_str else start_dt.replace(
            hour=0, minute=0, second=0, microsecond=0
        )
        reference_date = ref_str or ref_dt.strftime('%Y%m%d')
        boundary_names = [s['boundary_name'] for s in parsed_sections]
        logger.info(
            f"Existing file: {len(parsed_sections)} sections, "
            f"reference-date {reference_date}"
        )
    else:  # model_files
        if is_wave:
            if not (mdw_path and bnd_path and grd_path):
                raise ValueError("Waves require mdw_path, bnd_path, grd_path")
            # .mdw has no reference time; use user's start date midnight.
            ref_dt = start_dt.replace(hour=0, minute=0, second=0, microsecond=0)
            reference_date = ref_dt.strftime('%Y%m%d')
            bnd_pairs = _bnd_pairs_from_bnd_grd(bnd_path, grd_path)
            boundary_names = [bp['bnd_name'] if 'bnd_name' in bp
                              else bp['name'] for bp in bnd_pairs]
        else:
            if not (mdf_path and bnd_path):
                raise ValueError("Non-wave model_files mode requires mdf_path and bnd_path")
            mdf = fp.parse_mdf(mdf_path)
            timing = fp.extract_mdf_timing(mdf)
            reference_date = timing['reference_date']
            ref_dt = timing['start_datetime'].replace(
                hour=0, minute=0, second=0, microsecond=0
            ) if timing.get('start_datetime') else _ref_from_yyyymmdd(reference_date)
            # Use the reference_date's midnight as the zero of the minutes axis
            ref_dt = _ref_from_yyyymmdd(reference_date)
            boundary_names = _boundary_names_from_bnd(bnd_path)

    # Override reference date with user-supplied value if provided
    if user_reference_date and user_reference_date.strip():
        reference_date = user_reference_date.strip()
        ref_dt = _ref_from_yyyymmdd(reference_date)
        logger.info(f"Using user-specified reference date: {reference_date}")

    progress_callback(0.15)
    if cancel_event is not None and cancel_event.is_set():
        raise RuntimeError("Cancelled")

    # ── Build timeline + tstart/tstop in minutes since reference ────────
    tstart = _minutes_since(start_dt, ref_dt)
    tstop = _minutes_since(end_dt, ref_dt)
    if tstop <= tstart:
        raise ValueError("end_date must be strictly after start_date")

    timeline = np.arange(tstart, tstop + timestep_min, timestep_min)
    control_min = np.array([_minutes_since(dt, ref_dt) for dt in control_dt])
    logger.info(
        f"Timeline: {len(timeline)} steps, "
        f"{tstart:.1f} → {tstop:.1f} min, dt={timestep_min}"
    )

    # ── Interpolate values per file type ────────────────────────────────
    progress_callback(0.35)

    def _interp_col(vals: List[float]) -> List[float]:
        arr = np.asarray(vals, dtype=float)
        if arr.shape != control_min.shape:
            raise ValueError(
                f"Expected {len(control_min)} control values, got {arr.shape[0]}"
            )
        return _interp_series(timeline, control_min, arr, interpolation).tolist()

    # ── Dispatch to writer ──────────────────────────────────────────────
    progress_callback(0.55)
    mdw_out: Optional[str] = None

    if file_type == 'waterlevel_bct':
        boundary_data = {
            bnd: {
                'wl_a': _interp_col(values[bnd]['wl_a']),
                'wl_b': _interp_col(values[bnd]['wl_b']),
            } for bnd in boundary_names
        }
        if parsed_sections:
            section_values = [
                [boundary_data[s['boundary_name']]['wl_a'],
                 boundary_data[s['boundary_name']]['wl_b']]
                for s in parsed_sections
            ]
            write_bc_with_header_replay(
                output_path, parsed_sections, section_values,
                tstart, tstop, timestep_min,
            )
        else:
            write_bct(output_path, boundary_data, reference_date,
                      tstart, tstop, timestep_min)

    elif file_type == 'currents_bct':
        boundary_data = {
            bnd: {
                'u_a': _interp_col(values[bnd]['u_a']),
                'v_a': _interp_col(values[bnd]['v_a']),
                'u_b': _interp_col(values[bnd]['u_b']),
                'v_b': _interp_col(values[bnd]['v_b']),
            } for bnd in boundary_names
        }
        if parsed_sections:
            section_values = [
                [boundary_data[s['boundary_name']]['u_a'],
                 boundary_data[s['boundary_name']]['v_a'],
                 boundary_data[s['boundary_name']]['u_b'],
                 boundary_data[s['boundary_name']]['v_b']]
                for s in parsed_sections
            ]
            write_bc_with_header_replay(
                output_path, parsed_sections, section_values,
                tstart, tstop, timestep_min,
            )
        else:
            write_bct_currents(output_path, boundary_data, reference_date,
                               tstart, tstop, timestep_min)

    elif file_type == 'salinity_bcc':
        boundary_data = {
            bnd: {
                'sal_a': _interp_col(values[bnd]['sal_a']),
                'sal_b': _interp_col(values[bnd]['sal_b']),
            } for bnd in boundary_names
        }
        if parsed_sections:
            section_values = [
                [boundary_data[s['boundary_name']]['sal_a'],
                 boundary_data[s['boundary_name']]['sal_b']]
                for s in parsed_sections
            ]
            write_bc_with_header_replay(
                output_path, parsed_sections, section_values,
                tstart, tstop, timestep_min,
            )
        else:
            write_bcc_salinity(output_path, boundary_data, reference_date,
                               tstart, tstop, timestep_min)

    elif file_type == 'sediments_bcc':
        if not fraction_names:
            raise ValueError("sediments_bcc requires fraction_names")
        boundary_data = {}
        for frac in fraction_names:
            if frac not in values:
                raise ValueError(f"values missing fraction '{frac}'")
            boundary_data[frac] = {
                bnd: {
                    'conc_a': _interp_col(values[frac][bnd]['conc_a']),
                    'conc_b': _interp_col(values[frac][bnd]['conc_b']),
                } for bnd in boundary_names
            }
        if parsed_sections:
            # Existing file may have n_fractions * n_boundaries sections in
            # fraction-major order.
            n_bnd = len(boundary_names)
            section_values = []
            for sec_idx, sec in enumerate(parsed_sections):
                frac_idx = sec_idx // n_bnd
                if frac_idx >= len(fraction_names):
                    raise ValueError(
                        f"Section count mismatch: existing file has "
                        f"{len(parsed_sections)} sections but only "
                        f"{len(fraction_names)} fractions x {n_bnd} boundaries "
                        f"supplied"
                    )
                frac = fraction_names[frac_idx]
                bnd = sec['boundary_name']
                section_values.append([
                    boundary_data[frac][bnd]['conc_a'],
                    boundary_data[frac][bnd]['conc_b'],
                ])
            write_bc_with_header_replay(
                output_path, parsed_sections, section_values,
                tstart, tstop, timestep_min,
            )
        else:
            write_bcc(output_path, boundary_data, fraction_names,
                      reference_date, tstart, tstop, timestep_min)

    elif file_type == 'waves_bcw':
        time_col = timeline.tolist()
        boundary_data = {}
        for bnd in boundary_names:
            boundary_data[bnd] = {
                'time': time_col,
                'hs':     _interp_col(values[bnd]['hs']),
                'tp':     _interp_col(values[bnd]['tp']),
                'dir':    _interp_col(values[bnd]['dir']),
                'spread': _interp_col(values[bnd]['spread']),
            }
        write_bcw(output_path, boundary_data, reference_date)

        if write_companion_mdw and mdw_path and bnd_pairs:
            mdw_out = mdw_output_path or str(
                Path(output_path).with_suffix('.mdw')
            )
            # write_mdw expects the bnd_pairs shape with keys name,x_a,y_a,x_b,y_b
            mdw_bnd_pairs = [
                {'name': bp.get('bnd_name', bp.get('name')),
                 'x_a': bp['x_a'], 'y_a': bp['y_a'],
                 'x_b': bp['x_b'], 'y_b': bp['y_b']}
                for bp in bnd_pairs
            ]
            write_mdw(mdw_path, mdw_bnd_pairs, mdw_out)
            logger.info(f"Companion .mdw written: {mdw_out}")

    else:
        raise ValueError(f"Unhandled file_type: {file_type}")

    progress_callback(0.9)

    preview = _preview_lines(output_path)
    for line in preview:
        logger.info(line)

    progress_callback(1.0)
    logger.info(f"Done: {output_path}")

    return {
        'output_path': output_path,
        'n_timesteps': len(timeline),
        'preview_lines': preview,
        'mdw_path': mdw_out,
    }


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import tempfile
    import os

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    print("=== backend_synthetic_bc.py self-test ===\n")
    all_pass = True

    def check(name, result, expected, tol=None):
        global all_pass
        if tol is not None:
            ok = abs(result - expected) < tol
        else:
            ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]  got {result!r} expected {expected!r}")
        else:
            print(f"  {name}: [{status}]")

    # Build a tiny .mdf + .bnd + .grd for model_files mode
    tmp = Path(tempfile.mkdtemp())
    mdf_p = tmp / 't.mdf'
    bnd_p = tmp / 't.bnd'
    grd_p = tmp / 't.grd'

    mdf_p.write_text(
        "Ident = #test#\n"
        "Itdate = #2015-01-01#\n"
        "Tstart = 0.0000000e+000\n"
        "Tstop  = 2.8800000e+003\n"
    )
    bnd_p.write_text(
        "South               Z Z   1   1   1  50   0.0\n"
        "West                Z Z  50   1 100   1   0.0\n"
    )
    # Minimal 3x2 grid
    grd_p.write_text(
        "Coordinate System = Cartesian\n"
        "         0         0         0\n"
        " Missing Value =   0.0000000e+000\n"
        "       100        50\n"
        " 0 0 0\n"
        + ''.join(
            f"ETA=  {r}  " + ' '.join(f"{100.0*c:.2f}" for c in range(1, 101)) + "\n"
            for r in range(1, 51)
        )
        + ''.join(
            f"ETA=  {r}  " + ' '.join(f"{1000.0*c:.2f}" for c in range(1, 101)) + "\n"
            for r in range(1, 51)
        )
    )

    # ── Test 1: linear interpolation waterlevel ──────────────────────────
    print("linear waterlevel bct:")
    out_p = tmp / 'out_linear.bct'
    res = generate_synthetic_bc(
        file_type='waterlevel_bct',
        source_mode='model_files',
        start_date='2015-01-01 00:00:00',
        end_date='2015-01-01 04:00:00',
        timestep_min=60.0,
        control_dates=['2015-01-01 00:00:00',
                       '2015-01-01 02:00:00',
                       '2015-01-01 04:00:00'],
        values={
            'South': {'wl_a': [0.0, 2.0, 0.0], 'wl_b': [0.0, 1.0, 0.0]},
            'West':  {'wl_a': [1.0, 1.0, 1.0], 'wl_b': [2.0, 2.0, 2.0]},
        },
        interpolation='linear',
        output_path=str(out_p),
        mdf_path=str(mdf_p), bnd_path=str(bnd_p), grd_path=str(grd_p),
    )
    check("n_timesteps", res['n_timesteps'], 5)
    check("preview is list", isinstance(res['preview_lines'], list), True)

    # Inspect South wl_a: should be [0, 1, 2, 1, 0] after linear interp
    content = out_p.read_text()
    lines = content.splitlines()
    south_data = [l.strip() for l in lines
                  if l.strip() and 'e+' in l.lower()
                  and 'parameter' not in l.lower()
                  and 'Boundary' not in l
                  and 'reference' not in l
                  and 'time-' not in l
                  and 'interpolation' not in l
                  and 'location' not in l
                  and 'records' not in l
                  and 'contents' not in l
                  and 'table-' not in l][:5]  # first 5 data rows = South
    wl_a_row2 = float(south_data[1].split()[1])
    wl_a_row3 = float(south_data[2].split()[1])
    check("linear mid wl_a row1 ≈ 1.0", wl_a_row2, 1.0, tol=1e-6)
    check("linear peak wl_a row2 ≈ 2.0", wl_a_row3, 2.0, tol=1e-6)

    # ── Test 2: constant (step-hold) interpolation ───────────────────────
    print("\nconstant waterlevel bct:")
    out_p2 = tmp / 'out_const.bct'
    generate_synthetic_bc(
        file_type='waterlevel_bct',
        source_mode='model_files',
        start_date='2015-01-01 00:00:00',
        end_date='2015-01-01 04:00:00',
        timestep_min=60.0,
        control_dates=['2015-01-01 00:00:00',
                       '2015-01-01 02:00:00',
                       '2015-01-01 04:00:00'],
        values={
            'South': {'wl_a': [0.0, 2.0, 5.0], 'wl_b': [0.0, 1.0, 0.0]},
            'West':  {'wl_a': [1.0, 1.0, 1.0], 'wl_b': [2.0, 2.0, 2.0]},
        },
        interpolation='constant',
        output_path=str(out_p2),
        mdf_path=str(mdf_p), bnd_path=str(bnd_p), grd_path=str(grd_p),
    )
    content2 = out_p2.read_text()
    lines2 = content2.splitlines()
    south2 = [l.strip() for l in lines2
              if l.strip() and 'e+' in l.lower()
              and not any(k in l.lower() for k in
                          ('parameter', 'boundary', 'reference', 'time-',
                           'interpolation', 'location', 'records', 'contents',
                           'table-'))][:5]
    # At t=60 (< 120), still control[0]=0.0
    # At t=120 (==control[1]), value = 2.0
    v_t60 = float(south2[1].split()[1])
    v_t120 = float(south2[2].split()[1])
    check("constant t=60 holds 0.0", v_t60, 0.0, tol=1e-6)
    check("constant t=120 jumps to 2.0", v_t120, 2.0, tol=1e-6)

    # ── Test 3: validation rejects bad control dates ─────────────────────
    print("\nvalidation:")
    try:
        generate_synthetic_bc(
            file_type='waterlevel_bct', source_mode='model_files',
            start_date='2015-01-01 00:00:00', end_date='2015-01-01 04:00:00',
            timestep_min=60.0,
            control_dates=['2015-01-01 00:00:00', '2015-01-01 04:00:00',
                           '2015-01-01 02:00:00'],
            values={'South': {'wl_a':[0,0,0], 'wl_b':[0,0,0]},
                    'West':  {'wl_a':[0,0,0], 'wl_b':[0,0,0]}},
            interpolation='linear',
            output_path=str(tmp / 'bad.bct'),
            mdf_path=str(mdf_p), bnd_path=str(bnd_p), grd_path=str(grd_p),
        )
        check("non-monotonic dates raise", False, True)
    except ValueError:
        check("non-monotonic dates raise", True, True)

    # ── Test 4: existing-file header replay ──────────────────────────────
    print("\nexisting-file replay:")
    seed_path = tmp / 'seed.bct'
    from writer_delft3d4 import write_bct as _wbct
    _wbct(str(seed_path),
          {'South': {'wl_a':[0.1, 0.2], 'wl_b':[0.3, 0.4]},
           'West':  {'wl_a':[0.5, 0.6], 'wl_b':[0.7, 0.8]}},
          '20200515', 0.0, 60.0, 60.0)
    out_replay = tmp / 'replay.bct'
    res_r = generate_synthetic_bc(
        file_type='waterlevel_bct',
        source_mode='existing_file',
        start_date='2020-05-15 00:00:00',
        end_date='2020-05-15 02:00:00',
        timestep_min=60.0,
        control_dates=['2020-05-15 00:00:00',
                       '2020-05-15 01:00:00',
                       '2020-05-15 02:00:00'],
        values={
            'South': {'wl_a':[9.0, 9.5, 10.0], 'wl_b':[1.0, 1.5, 2.0]},
            'West':  {'wl_a':[0.0, 0.5, 1.0], 'wl_b':[2.0, 2.5, 3.0]},
        },
        interpolation='linear',
        output_path=str(out_replay),
        existing_file_path=str(seed_path),
    )
    check("replay n_timesteps", res_r['n_timesteps'], 3)
    replayed = out_replay.read_text()
    check("replay reuses 20200515 ref", '20200515' in replayed, True)
    check("replay has South", "'South" in replayed, True)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
