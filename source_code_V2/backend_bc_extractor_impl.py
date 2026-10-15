# -*- coding: utf-8 -*-
"""
backend_bc_extractor.py — Boundary condition extraction engine.

This module merges 4 v1 files into one unified engine:
    - bct_generator.py              (single-year water level)
    - bct_year_overlap_file_generator.py (multi-year water level)
    - bcw_generator.py              (single-year waves)
    - bcw_year_overlap_file_generator.py (multi-year waves)

The multi-file/single-file distinction disappears: xr.open_mfdataset
handles both seamlessly (this needs dask — see _open_dataset).

The data-source distinction disappears: presets provide variable names,
coordinate names, time offsets, etc.

Flow:
    1. Load model definition (.mdf or .mdu) → extract timing
    2. Load boundary definitions (.bnd+.grd for D3D4, .pli for FM)
    3. Convert boundary coordinates to lat/lon
    4. Load data file(s) using preset configuration
    5. Apply time offset (EasyGSH: -12h)
    6. For each boundary point: bilinear interpolation across time
    7. For waves: convert direction vectors → nautical degrees
    8. Organize data into writer-ready format
    9. Write output via writer_delft3d4 or writer_delft3dfm

Usage:
    from backend_bc_extractor import extract_water_level, extract_waves
"""

import math
import logging
import threading
import time
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from datetime import timedelta
from typing import Callable, List, Optional, Dict, Any, Tuple

import file_parsers
import interpolation as interp
import coordinate_utils
import presets as preset_module
import writer_delft3d4
import writer_delft3dfm
import writer_mdw


class TaskCancelledError(Exception):
    """Raised when the backend observes a cancel_event set by the GUI."""


def _check_cancel(cancel_event: Optional[threading.Event]) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise TaskCancelledError("Cancelled by user")


QC_MIN_VALID_FRACTION = 0.75
QC_MAX_BAD_WEIGHT = 0.25
QC_MAX_SINGLE_BAD_CORNER_WEIGHT = 0.05
QC_EDGE_GAP_FRACTION = 0.02
QC_EDGE_GAP_CAP = 240
QC_INTERNAL_GAP_FRACTION = 0.01
QC_INTERNAL_GAP_CAP = 120
REPAIR_BINARY_STEPS = 6


def _gap_limit(n_steps: int, fraction: float, cap: int) -> int:
    if n_steps <= 0:
        return 0
    return max(12, min(cap, int(n_steps * fraction)))


def _assess_stencil_quality(stencil: Dict[str, Any]) -> Dict[str, Any]:
    """Decide whether a bilinear stencil is reliable enough for extraction."""

    corner_values = list(stencil['corners'].values())
    n_steps = corner_values[0]['n_steps'] if corner_values else 0
    edge_limit = _gap_limit(n_steps, QC_EDGE_GAP_FRACTION, QC_EDGE_GAP_CAP)
    internal_limit = _gap_limit(n_steps, QC_INTERNAL_GAP_FRACTION, QC_INTERNAL_GAP_CAP)

    bad_weight = 0.0
    max_bad_corner_weight = 0.0
    issues = []

    if stencil['outside']:
        issues.append('point lies outside the dataset extent')

    for corner_name, corner in stencil['corners'].items():
        weight = float(corner['normalized_weight'])
        reasons = []
        if corner['all_nan']:
            reasons.append('all-NaN')
        else:
            if corner['valid_fraction'] < QC_MIN_VALID_FRACTION:
                reasons.append(f"valid fraction {corner['valid_fraction']:.1%}")
            if corner['leading_nan_count'] > edge_limit:
                reasons.append(f"leading gap {corner['leading_nan_count']}")
            if corner['trailing_nan_count'] > edge_limit:
                reasons.append(f"trailing gap {corner['trailing_nan_count']}")
            if corner['internal_nan_count'] > internal_limit:
                reasons.append(f"internal gap {corner['internal_nan_count']}")
        if not reasons:
            continue
        bad_weight += weight
        max_bad_corner_weight = max(max_bad_corner_weight, weight)
        issues.append(
            f"{corner_name} ({weight:.1%} weight): " + ', '.join(reasons)
        )

    reliable = (
        not stencil['outside']
        and bad_weight <= QC_MAX_BAD_WEIGHT
        and max_bad_corner_weight <= QC_MAX_SINGLE_BAD_CORNER_WEIGHT
    )

    return {
        'reliable': reliable,
        'bad_weight': bad_weight,
        'max_bad_corner_weight': max_bad_corner_weight,
        'issues': issues,
        'edge_gap_limit': edge_limit,
        'internal_gap_limit': internal_limit,
    }


def _assess_derived_stencil(stencil: Dict[str, Any]) -> Dict[str, Any]:
    """
    Reliability of a derived quantity (wave peak period, directional
    spreading; see backend_bc_extractor._wave_derived_variables): such a
    field is undefined wherever there are no waves, so gaps are expected and
    the only disqualifying stencil is one without a value at any timestep in
    any contributing corner (the mask-normalised interpolation could then not
    produce a series).
    """
    has_value = any(
        float(corner['normalized_weight']) > 0.0 and not corner['all_nan']
        for corner in stencil['corners'].values()
    )
    reliable = has_value and not stencil['outside']
    issues = []
    if stencil['outside']:
        issues.append('point lies outside the dataset extent')
    if not has_value:
        issues.append('no source value at any timestep in any contributing corner')
    return {
        'reliable': reliable,
        'bad_weight': 0.0 if reliable else 1.0,
        'max_bad_corner_weight': 0.0 if reliable else 1.0,
        'issues': issues,
        'edge_gap_limit': 0,
        'internal_gap_limit': 0,
    }


def _assessment_summary(variable: str, assessment: Dict[str, Any]) -> str:
    """Format the most useful stencil issues into one short log line."""

    details = assessment['issues'][:3]
    msg = '; '.join(details) if details else 'insufficient valid source coverage'
    if len(assessment['issues']) > 3:
        msg += f"; +{len(assessment['issues']) - 3} more"
    return (
        f"{variable}: {msg} (bad weight "
        f"{assessment['bad_weight']:.1%})"
    )


def _point_cache_key(lat: float, lon: float) -> Tuple[float, float]:
    """Stable cache key for boundary point QC lookups."""

    return (round(float(lat), 10), round(float(lon), 10))


def _build_qc_context(ds, variables: List[str], preset, logger=None, cancel_event=None) -> Dict[str, Any]:
    """Precompute per-variable validity masks and per-cell statistics for QC.

    One pass over each required variable (block by block along the time axis,
    see interpolation.cell_validity_stats) yields, per grid cell, how many
    timesteps are valid and where the first and last valid samples sit. The
    fast-QC masks are derived from the counts exactly as before (valid at
    every step / at any step), and the strict and intertidal checks take
    their per-corner gap statistics from the same arrays instead of pulling
    four corner time series off the disk per variable and point. On
    one-field-per-chunk files (TrilaWatt, EasyGSH) such a read decompresses
    every chunk of the window — ~50 s per corner series, 17 min per
    strict-QC point for the 2017+2018 TrilaWatt waves (27649 steps,
    measured 2026-10-15) — so a 90-section wave run spent ~20 h in QC alone.
    """

    required_vars = list(dict.fromkeys(v for v in variables if v))
    metrics = {}
    cancel_check = (lambda: _check_cancel(cancel_event)) if cancel_event is not None else None
    if required_vars and logger is not None:
        logger.info(
            f"Scanning {len(required_vars)} variable(s) for valid-data coverage "
            f"(one pass over the source window)..."
        )

    for variable in required_vars:
        started = time.perf_counter()
        stats = interp.cell_validity_stats(
            ds[variable], preset.time_coord, preset.lat_coord, preset.lon_coord,
            cancel_check=cancel_check,
        )
        if stats is None:
            # Unexpected dimension layout: the original eager masks, and the
            # strict/intertidal checks read their stencils as before.
            if logger is not None:
                logger.warning(
                    f"  {variable}: dimensions {tuple(ds[variable].dims)} are not "
                    f"({preset.time_coord}, {preset.lat_coord}, {preset.lon_coord}); "
                    f"QC statistics unavailable, stencils will be read per point (slow)."
                )
            valid = ds[variable].notnull()
            if preset.time_coord in valid.dims:
                full_valid = valid.all(dim=preset.time_coord)
                any_valid = valid.any(dim=preset.time_coord)
            else:
                full_valid = valid
                any_valid = valid

            metrics[variable] = {
                'full_valid_mask': np.asarray(
                    full_valid.transpose(preset.lat_coord, preset.lon_coord).values,
                    dtype=bool,
                ),
                'any_valid_mask': np.asarray(
                    any_valid.transpose(preset.lat_coord, preset.lon_coord).values,
                    dtype=bool,
                ),
            }
            continue

        n_steps, count, first, last = stats
        metrics[variable] = {
            'full_valid_mask': count == n_steps,
            'any_valid_mask': count > 0,
            'n_steps': int(n_steps),
            'valid_count': count,
            'first_valid_index': first,
            'last_valid_index': last,
        }
        if logger is not None:
            n_full = int((count == n_steps).sum())
            n_none = int((count == 0).sum())
            logger.info(
                f"  {variable}: {n_full} of {count.size} cells valid at every step, "
                f"{count.size - n_full - n_none} with gaps, {n_none} never valid "
                f"({time.perf_counter() - started:.1f} s)"
            )

    return {
        'required_vars': tuple(required_vars),
        'metrics': metrics,
        'strict_assessment_cache': {},
        'point_result_cache': {},
        'stats': {'strict_point_calls': 0},
    }


def _fast_assess_variable_quality(
    variable: str,
    geometry: Dict[str, Any],
    qc_context: Dict[str, Any],
) -> Dict[str, Any]:
    """Classify one variable as fast-pass, hard-fail, or ambiguous."""

    masks = qc_context['metrics'][variable]
    full_valid_mask = masks['full_valid_mask']
    any_valid_mask = masks['any_valid_mask']
    lat_idx1, lat_idx2 = geometry['lat_indices']
    lon_idx1, lon_idx2 = geometry['lon_indices']
    corner_indices = {
        'Q11': (lat_idx1, lon_idx1),
        'Q12': (lat_idx1, lon_idx2),
        'Q21': (lat_idx2, lon_idx1),
        'Q22': (lat_idx2, lon_idx2),
    }

    all_contributing_full = True
    bad_weight = 0.0
    max_bad_corner_weight = 0.0
    issues = []

    for corner_name, (lat_idx, lon_idx) in corner_indices.items():
        weight = float(geometry['normalized_weights'][corner_name])
        if weight <= 0.0:
            continue

        full_valid = bool(full_valid_mask[lat_idx, lon_idx])
        any_valid = bool(any_valid_mask[lat_idx, lon_idx])
        if not full_valid:
            all_contributing_full = False
        if any_valid:
            continue
        bad_weight += weight
        max_bad_corner_weight = max(max_bad_corner_weight, weight)
        issues.append(f"{corner_name} ({weight:.1%} weight): all-NaN")

    if all_contributing_full:
        decision = 'accept'
    elif bad_weight > QC_MAX_BAD_WEIGHT or max_bad_corner_weight > QC_MAX_SINGLE_BAD_CORNER_WEIGHT:
        decision = 'hard_fail'
    else:
        decision = 'ambiguous'

    return {
        'decision': decision,
        'bad_weight': bad_weight,
        'max_bad_corner_weight': max_bad_corner_weight,
        'issues': issues,
    }


def _fast_evaluate_point_quality(
    lat: float,
    lon: float,
    ds,
    variables: List[str],
    preset,
    qc_context: Dict[str, Any],
) -> Dict[str, Any]:
    """Quickly screen obvious good and obvious bad points from mask coverage."""

    geometry = interp.inspect_bilinear_geometry(
        lat,
        lon,
        ds,
        lat_coord=preset.lat_coord,
        lon_coord=preset.lon_coord,
    )

    if geometry['outside']:
        return {
            'decision': 'hard_fail',
            'issues': ['point lies outside the dataset extent'],
        }

    all_accept = True
    hard_fail_issues = []
    for variable in variables:
        assessment = _fast_assess_variable_quality(variable, geometry, qc_context)
        if assessment['decision'] != 'accept':
            all_accept = False
        if assessment['decision'] == 'hard_fail':
            hard_fail_issues.append(_assessment_summary(variable, assessment))

    if all_accept:
        return {'decision': 'accept', 'issues': []}
    if hard_fail_issues:
        return {'decision': 'hard_fail', 'issues': hard_fail_issues}
    return {'decision': 'ambiguous', 'issues': []}


def _stencil_from_stats(lat, lon, ds, variable, preset, metrics) -> Dict[str, Any]:
    """
    interpolation.inspect_bilinear_stencil() rebuilt from the per-cell
    statistics of the QC context — same keys, same numbers, no file access.
    """
    geometry = interp.inspect_bilinear_geometry(
        lat, lon, ds, lat_coord=preset.lat_coord, lon_coord=preset.lon_coord
    )
    n_steps = int(metrics['n_steps'])
    count = metrics['valid_count']
    first = metrics['first_valid_index']
    last = metrics['last_valid_index']
    lat_idx1, lat_idx2 = geometry['lat_indices']
    lon_idx1, lon_idx2 = geometry['lon_indices']
    corner_cells = {
        'Q11': (lat_idx1, lon_idx1),
        'Q12': (lat_idx1, lon_idx2),
        'Q21': (lat_idx2, lon_idx1),
        'Q22': (lat_idx2, lon_idx2),
    }
    corners = {}
    for key, (i, j) in corner_cells.items():
        valid_count = int(count[i, j])
        if valid_count == 0:
            first_valid = None
            last_valid = None
            leading = n_steps
            trailing = n_steps
            internal = 0
        else:
            first_valid = int(first[i, j])
            last_valid = int(last[i, j])
            leading = first_valid
            trailing = n_steps - last_valid - 1
            internal = (last_valid - first_valid + 1) - valid_count
        corners[key] = {
            'n_steps': n_steps,
            'valid_count': valid_count,
            'valid_fraction': (valid_count / n_steps) if n_steps else 0.0,
            'all_nan': (valid_count == 0),
            'first_valid_index': first_valid,
            'last_valid_index': last_valid,
            'leading_nan_count': leading,
            'trailing_nan_count': trailing,
            'internal_nan_count': internal,
            'normalized_weight': float(geometry['normalized_weights'][key]),
        }
    return {
        'lat': geometry['lat'],
        'lon': geometry['lon'],
        'outside': geometry['outside'],
        'lat_outside': geometry['lat_outside'],
        'lon_outside': geometry['lon_outside'],
        'point_on_grid': geometry['point_on_grid'],
        'lat_indices': geometry['lat_indices'],
        'lon_indices': geometry['lon_indices'],
        'lat_coords': geometry['lat_coords'],
        'lon_coords': geometry['lon_coords'],
        'weights': dict(geometry['weights']),
        'total_weight': float(geometry['total_weight']),
        'corners': corners,
    }


def _inspect_stencil(lat, lon, ds, variable, preset, qc_context=None) -> Dict[str, Any]:
    """Stencil summary for QC: from the context's per-cell statistics when present, else read from disk."""
    metrics = ((qc_context or {}).get('metrics') or {}).get(variable)
    if metrics and 'valid_count' in metrics:
        return _stencil_from_stats(lat, lon, ds, variable, preset, metrics)
    return interp.inspect_bilinear_stencil(
        lat,
        lon,
        ds,
        variable=variable,
        lat_coord=preset.lat_coord,
        lon_coord=preset.lon_coord,
    )


def _evaluate_point_quality(
    lat: float,
    lon: float,
    ds,
    variables: List[str],
    preset,
    qc_context: Dict[str, Any],
) -> Dict[str, Any]:
    """Evaluate a boundary point against all variables needed for one extraction."""

    cache = qc_context['strict_assessment_cache']
    qc_context['stats']['strict_point_calls'] += 1
    derived = set(qc_context.get('derived_vars', ()))
    assessments = {}
    issues = []
    for variable in variables:
        key = (variable, *_point_cache_key(lat, lon))
        if key not in cache:
            stencil = _inspect_stencil(lat, lon, ds, variable, preset, qc_context)
            if variable in derived:
                cache[key] = _assess_derived_stencil(stencil)
            else:
                cache[key] = _assess_stencil_quality(stencil)
        assessment = cache[key]
        assessments[variable] = assessment
        if not assessment['reliable']:
            issues.append(_assessment_summary(variable, assessment))

    return {
        'reliable': not issues,
        'assessments': assessments,
        'issues': issues,
    }


def _evaluate_point_status(
    lat: float,
    lon: float,
    ds,
    variables: List[str],
    preset,
    qc_context: Dict[str, Any],
) -> Dict[str, Any]:
    """Run fast QC first, then strict QC only for ambiguous points."""

    point_cache = qc_context['point_result_cache']
    key = _point_cache_key(lat, lon)
    if key in point_cache:
        return point_cache[key]

    fast = _fast_evaluate_point_quality(
        lat, lon, ds, variables, preset, qc_context
    )

    if fast['decision'] == 'accept':
        result = {
            'reliable': True,
            'path': 'fast',
            'issues': [],
        }
    elif fast['decision'] == 'hard_fail':
        result = {
            'reliable': False,
            'path': 'fast',
            'issues': fast['issues'],
        }
    else:
        strict = _evaluate_point_quality(
            lat, lon, ds, variables, preset, qc_context
        )
        result = {
            'reliable': strict['reliable'],
            'path': 'strict',
            'issues': strict['issues'],
        }

    point_cache[key] = result
    return result


def _count_qc_pass(summary: Dict[str, Any], quality: Dict[str, Any]) -> None:
    """Track how each original point cleared QC."""

    if not quality['reliable']:
        return
    if quality['path'] == 'fast':
        summary['fast_pass_points'] += 1
    else:
        summary['strict_pass_points'] += 1


def _distance_meters(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Approximate distance between two geographic points in metres."""

    rad = math.pi / 180.0
    dlat = (lat2 - lat1) * rad
    dlon = (lon2 - lon1) * rad
    a = (
        math.sin(dlat / 2.0) ** 2
        + math.cos(lat1 * rad) * math.cos(lat2 * rad) * math.sin(dlon / 2.0) ** 2
    )
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(max(1.0 - a, 0.0)))
    return 6371000.0 * c


def _repair_point_toward_targets(
    start_lat: float,
    start_lon: float,
    targets: List[Tuple[str, float, float]],
    ds,
    variables: List[str],
    preset,
    qc_context: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    """
    Binary-search from a bad point toward already-good neighbours.
    """

    best = None
    for target_name, target_lat, target_lon in targets:
        target_quality = _evaluate_point_status(
            target_lat, target_lon, ds, variables, preset, qc_context
        )

        if not target_quality['reliable']:
            continue

        best_candidate = {
            'lat': target_lat,
            'lon': target_lon,
            'fraction': 1.0,
            'distance_m': _distance_meters(
                start_lat, start_lon, target_lat, target_lon
            ),
            'target_name': target_name,
        }

        low = 0.0
        high = 1.0

        for _ in range(REPAIR_BINARY_STEPS):
            fraction = (low + high) / 2.0
            cand_lat = start_lat + fraction * (target_lat - start_lat)
            cand_lon = start_lon + fraction * (target_lon - start_lon)
            quality = _evaluate_point_status(
                cand_lat, cand_lon, ds, variables, preset, qc_context
            )

            if quality['reliable']:
                best_candidate = {
                    'lat': cand_lat,
                    'lon': cand_lon,
                    'fraction': fraction,
                    'distance_m': _distance_meters(
                        start_lat, start_lon, cand_lat, cand_lon
                    ),
                    'target_name': target_name,
                }
                high = fraction
            else:
                low = fraction

        if best is None or best_candidate['distance_m'] < best['distance_m']:
            best = best_candidate

    return best


def _prepare_boundary_sampling(
    model_type: str,
    extraction_label: str,
    ds,
    bnd_loc_geo: pd.DataFrame,
    bnd_pairs: List[Dict[str, Any]],
    variables: List[str],
    preset,
    qc_context: Dict[str, Any],
    logger,
    waves_use_a_only: bool = False,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> Tuple[pd.DataFrame, List[Dict[str, Any]], Dict[str, Any]]:
    """
    Validate boundary sampling points against the source dataset and repair
    obvious edge cases before extraction begins.

    progress_callback receives values in [0.0, 1.0] for the QC phase only.
    """

    required_vars = list(qc_context.get('required_vars', ()))
    sample_geo = bnd_loc_geo.copy()
    mdw_pairs = list(bnd_pairs)
    summary = {
        'repaired_points': [],
        'rejected_boundaries': [],
        'fast_pass_points': 0,
        'strict_pass_points': 0,
    }

    if not required_vars or sample_geo.empty:
        return sample_geo, mdw_pairs, summary

    logger.info(
        f"QC boundary sampling points for {extraction_label} against "
        f"{len(required_vars)} variable(s)..."
    )

    if model_type == 'delft3d4':
        original = {
            row['boundaries']: (float(row['lat']), float(row['lon']))
            for _, row in bnd_loc_geo.iterrows()
        }

        rejected = set()
        boundary_names = sorted({name[:-2] for name in original})
        n_boundaries = len(boundary_names)
        active_suffixes = ('a',) if waves_use_a_only else ('a', 'b')

        for bnd_idx, boundary_name in enumerate(boundary_names):
            _check_cancel(cancel_event)
            point_results = {}
            pending_updates = []
            failing_points = []

            for suffix in active_suffixes:
                point_name = f"{boundary_name}_{suffix}"
                lat, lon = original[point_name]
                quality = _evaluate_point_status(
                    lat, lon, ds, required_vars, preset, qc_context
                )

                point_results[point_name] = quality
                _count_qc_pass(summary, quality)

            for suffix in active_suffixes:
                point_name = f"{boundary_name}_{suffix}"
                quality = point_results[point_name]
                if quality['reliable']:
                    continue
                other_suffix = 'b' if suffix == 'a' else 'a'
                target_name = f"{boundary_name}_{other_suffix}"
                target_quality = point_results.get(target_name)
                if target_quality is None:
                    target_lat, target_lon = original[target_name]
                    target_quality = _evaluate_point_status(
                        target_lat, target_lon, ds, required_vars, preset, qc_context
                    )

                    point_results[target_name] = target_quality

                lat, lon = original[point_name]
                repair = None
                if target_quality['reliable']:
                    target_lat, target_lon = original[target_name]
                    repair = _repair_point_toward_targets(
                        lat,
                        lon,
                        [(target_name, target_lat, target_lon)],
                        ds,
                        required_vars,
                        preset,
                        qc_context,
                    )
                if repair is None:
                    failing_points.append((point_name, quality['issues']))
                else:
                    pending_updates.append((point_name, repair, quality['issues']))

            if failing_points:
                rejected.add(boundary_name)
                summary['rejected_boundaries'].append(boundary_name)
                bad_points = ', '.join(name for name, _ in failing_points)
                logger.warning(
                    f"[qc] {extraction_label}: rejected boundary {boundary_name} because "
                    f"{bad_points} still have no reliable source data."
                )

                for point_name, issues in failing_points:
                    logger.warning(
                        f"[qc]   {point_name}: " + ' | '.join(issues[:3])
                    )

                progress_callback((bnd_idx + 1) / n_boundaries)
                continue

            for point_name, repair, issues in pending_updates:
                sample_geo.loc[
                    sample_geo['boundaries'] == point_name, ['lat', 'lon']
                ] = [repair['lat'], repair['lon']]
                msg = (
                    f"{point_name} moved {repair['distance_m']:.0f} m toward "
                    f"{repair['target_name']} for {extraction_label}"
                )
                summary['repaired_points'].append(msg)
                logger.warning(
                    f"[qc] {msg}. Reason: {' | '.join(issues[:2])}"
                )

            progress_callback((bnd_idx + 1) / n_boundaries)
            if (bnd_idx + 1) % max(1, n_boundaries // 5) == 0:
                logger.info(
                    f"  QC checked {bnd_idx + 1}/{n_boundaries} boundaries"
                )

        if rejected:
            sample_geo = sample_geo[
                ~sample_geo['boundaries'].str[:-2].isin(rejected)
            ].reset_index(drop=True)

            if waves_use_a_only:
                mdw_pairs = [pair for pair in mdw_pairs if pair['name'] not in rejected]

    else:
        original = {
            row['boundaries']: (float(row['lat']), float(row['lon']))
            for _, row in bnd_loc_geo.iterrows()
        }

        grouped = {}
        for point_name in original:
            base_name, idx_str = point_name.rsplit('_', 1)
            grouped.setdefault(base_name, []).append((int(idx_str), point_name))

        point_results = {}
        unresolved = []
        pending_updates = []
        n_fm_points = len(original)
        checked = 0
        for base_name, points in grouped.items():
            ordered = sorted(points)
            for pos, (_, point_name) in enumerate(ordered):
                _check_cancel(cancel_event)
                lat, lon = original[point_name]
                quality = _evaluate_point_status(
                    lat, lon, ds, required_vars, preset, qc_context
                )

                point_results[point_name] = quality
                _count_qc_pass(summary, quality)
                checked += 1
                if n_fm_points > 0:
                    progress_callback(checked / n_fm_points)
                if checked % max(1, n_fm_points // 5) == 0:
                    logger.info(f"  QC checked {checked}/{n_fm_points} FM points")

        for base_name, points in grouped.items():
            ordered = sorted(points)
            for pos, (_, point_name) in enumerate(ordered):
                quality = point_results[point_name]
                if quality['reliable']:
                    continue

                lat, lon = original[point_name]
                targets = []
                if pos > 0:
                    prev_name = ordered[pos - 1][1]
                    if point_results[prev_name]['reliable']:
                        prev_lat, prev_lon = original[prev_name]
                        targets.append((prev_name, prev_lat, prev_lon))
                if pos + 1 < len(ordered):
                    next_name = ordered[pos + 1][1]
                    if point_results[next_name]['reliable']:
                        next_lat, next_lon = original[next_name]
                        targets.append((next_name, next_lat, next_lon))

                repair = (
                    _repair_point_toward_targets(
                        lat, lon, targets, ds, required_vars, preset, qc_context
                    )
                    if targets else None
                )

                if repair is None:
                    unresolved.append((base_name, point_name, quality['issues']))
                else:
                    pending_updates.append((point_name, repair, quality['issues']))

        if unresolved:
            for base_name, point_name, issues in unresolved:
                logger.error(
                    f"[qc] {extraction_label}: FM point {point_name} on {base_name} has no reliable source data. "
                    f"{' | '.join(issues[:3])}"
                )

            bad_names = ', '.join(point_name for _, point_name, _ in unresolved[:6])
            if len(unresolved) > 6:
                bad_names += ', ...'
            raise ValueError(
                f"One or more FM boundary points have no reliable source data for "
                f"{extraction_label}: {bad_names}. Adjust the .pli geometry or use a source file with valid coverage."
            )

        for point_name, repair, issues in pending_updates:
            sample_geo.loc[
                sample_geo['boundaries'] == point_name, ['lat', 'lon']
            ] = [repair['lat'], repair['lon']]

            msg = (
                f"{point_name} moved {repair['distance_m']:.0f} m toward "
                f"{repair['target_name']} for {extraction_label}"
            )
            summary['repaired_points'].append(msg)
            logger.warning(
                f"[qc] {msg}. Reason: {' | '.join(issues[:2])}"
            )

    if sample_geo.empty:
        raise ValueError(
            f"No reliable boundary points remain for {extraction_label}. Extraction aborted."
        )

    if summary['repaired_points'] or summary['rejected_boundaries']:
        logger.warning(
            f"[qc] {extraction_label}: {len(summary['repaired_points'])} point(s) repaired, "
            f"{len(summary['rejected_boundaries'])} boundary section(s) rejected."
        )
    else:
        logger.info(f"QC passed for {extraction_label} - no repairs needed.")

    logger.info(
        f"[qc] {extraction_label}: fast-pass {summary['fast_pass_points']}, strict-pass "
        f"{summary['strict_pass_points']}, repaired {len(summary['repaired_points'])}, rejected "
        f"{len(summary['rejected_boundaries'])}"
    )

    return sample_geo, mdw_pairs, summary


def extract_water_level(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: Optional[str] = None,
    data_files: List[str] = None,
    preset_name: str = 'EasyGSH',
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = 'N',
    output_path: Optional[str] = None,
    output_dir: str = '.',
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Extract water level boundary conditions and write output files.

    Supports both Delft3D-4 (.bct) and Delft3D-FM (.bc) output.
    Supports single-file and multi-file (year-overlap) input seamlessly.

    Args:
        model_type:     'delft3d4' or 'delft3dfm'
        model_def_path: .mdf (D3D4) or .mdu (FM) file path
        boundary_path:  .bnd (D3D4) or .pli (FM) file path
        grid_path:      .grd file path (D3D4 only, None for FM)
        data_files:     List of NetCDF file paths
        preset_name:    Preset name from presets.py
        timestep:       Output timestep in minutes
        utm_zone:       UTM zone for Cartesian grids
        hemisphere:     'N' or 'S'
        output_path:    Output file path (auto-generated if None)
        output_dir:     Output directory
        logger:         Logger instance
        progress_callback: Progress callback (0.0–1.0)

    Returns:
        {'output_files': [...], 'n_boundaries': int, 'n_timesteps': int}
    """

    if logger is None:
        logger = logging.getLogger('repd3d.bc')

    preset = preset_module.get_bc_preset(preset_name)
    logger.info(f"Starting water level extraction — preset: {preset.name}")

    _check_cancel(cancel_event)
    timing, crs_type = _load_model_definition(model_type, model_def_path, grid_path, logger)

    _check_cancel(cancel_event)
    bnd_loc_geo, _bnd_native, _bnd_pairs = _load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )

    _check_cancel(cancel_event)
    ds = _open_dataset(data_files, preset, timing, logger)
    _check_cancel(cancel_event)
    qc_context = _build_qc_context(ds, [preset.water_level_var], preset)
    bnd_loc_geo, _unused_mdw_pairs, qc_summary = _prepare_boundary_sampling(
        model_type,
        'water level',
        ds,
        bnd_loc_geo,
        _bnd_pairs,
        [preset.water_level_var],
        preset,
        qc_context,
        logger,
        progress_callback=lambda p: progress_callback(p * 0.1),
        cancel_event=cancel_event,
    )
    progress_callback(0.1)

    logger.info(f"Extracting water level for {len(bnd_loc_geo)} boundary points...")
    output_dict = {}
    total_points = len(bnd_loc_geo)

    for idx, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        _check_cancel(cancel_event)
        wl = interp.bilinear_interpolation(
            row['lat'], row['lon'], ds,
            variable=preset.water_level_var,
            lat_coord=preset.lat_coord,
            lon_coord=preset.lon_coord,
        )
        output_dict[row['boundaries']] = wl
        progress_callback(0.1 + 0.7 * (idx + 1) / total_points)
        if (idx + 1) % max(1, total_points // 5) == 0:
            logger.info(f"  Processed {idx + 1}/{total_points} points")

    _check_cancel(cancel_event)
    logger.info("Water level extraction complete")

    if model_type == 'delft3d4':
        boundary_data = _assemble_bct_data(output_dict)
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '.bct'))
        writer_delft3d4.write_bct(
            output_path, boundary_data,
            timing['reference_date'],
            timing['tstart'], timing['tstop'], timestep,
        )
    else:
        boundary_data = _assemble_fm_data(output_dict)
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '.bc'))
        writer_delft3dfm.write_bc(
            output_path, boundary_data,
            'waterlevelbnd', 'm',
            str(timing['start_datetime']),
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")

    return {
        'output_files': [output_path],
        'n_boundaries': len(boundary_data),
        'n_timesteps': len(next(iter(output_dict.values()))),
        'repaired_points': qc_summary['repaired_points'],
        'rejected_boundaries': qc_summary['rejected_boundaries'],
    }


def extract_waves(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: Optional[str] = None,
    mdw_template_path: Optional[str] = None,
    data_files: List[str] = None,
    preset_name: str = 'EasyGSH',
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = 'N',
    output_path: Optional[str] = None,
    output_dir: str = '.',
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Extract wave boundary conditions and write output files.

    Extracts: Hs, Tp, direction, directional spreading.
    For vector-format direction (EasyGSH, TrilaWatt): converts x/y → degrees.
    For degree-format direction (Copernicus): uses directly.

    Output: .bcw (D3D4) or .bc (FM)

    Args: Same as extract_water_level.

    Returns:
        {'output_files': [...], 'n_boundaries': int, 'n_timesteps': int}
    """

    if logger is None:
        logger = logging.getLogger('repd3d.bc')

    preset = preset_module.get_bc_preset(preset_name)
    logger.info(f"Starting wave extraction — preset: {preset.name}")

    _check_cancel(cancel_event)
    timing, crs_type = _load_model_definition(model_type, model_def_path, grid_path, logger)
    _check_cancel(cancel_event)
    bnd_loc_geo, _bnd_native, bnd_pairs = _load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )

    _check_cancel(cancel_event)
    ds = _open_dataset(data_files, preset, timing, logger)
    _check_cancel(cancel_event)
    wave_qc_vars = [preset.wave_hs_var, preset.wave_tp_var]
    if preset.wave_dir_format == 'vector':
        wave_qc_vars.extend([preset.wave_dir_x_var, preset.wave_dir_y_var])
    else:
        wave_qc_vars.append(preset.wave_dir_var)
    if preset.wave_dir_spread_var and not preset.wave_dir_spread_use_default:
        wave_qc_vars.append(preset.wave_dir_spread_var)
    qc_context = _build_qc_context(ds, wave_qc_vars, preset)
    bnd_loc_geo, bnd_pairs, qc_summary = _prepare_boundary_sampling(
        model_type,
        'waves',
        ds,
        bnd_loc_geo,
        bnd_pairs,
        wave_qc_vars,
        preset,
        qc_context,
        logger,
        waves_use_a_only=True,
        progress_callback=lambda p: progress_callback(p * 0.1),
        cancel_event=cancel_event,
    )
    progress_callback(0.1)

    # Waves only use _a points (SWAN format requirement)
    bnd_loc_a = bnd_loc_geo[bnd_loc_geo['boundaries'].str.endswith('_a')].copy()
    if bnd_loc_a.empty:
        # FM boundaries don't have _a/_b suffix — use all points
        bnd_loc_a = bnd_loc_geo

    total_points = len(bnd_loc_a)
    logger.info(f"Extracting wave data for {total_points} boundary points...")

    wave_data = {}
    for idx, (_, row) in enumerate(bnd_loc_a.iterrows()):
        _check_cancel(cancel_event)
        bnd_name = row['boundaries']
        point_data = {}

        point_data['hs'] = interp.bilinear_interpolation(
            row['lat'], row['lon'], ds,
            variable=preset.wave_hs_var,
            lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
        ).tolist()

        point_data['tp'] = interp.bilinear_interpolation(
            row['lat'], row['lon'], ds,
            variable=preset.wave_tp_var,
            lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
        ).tolist()

        if preset.wave_dir_format == 'vector':
            dir_x = interp.bilinear_interpolation(
                row['lat'], row['lon'], ds,
                variable=preset.wave_dir_x_var,
                lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
            )
            dir_y = interp.bilinear_interpolation(
                row['lat'], row['lon'], ds,
                variable=preset.wave_dir_y_var,
                lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
            )
            point_data['dir'] = _wave_direction_from_vectors(dir_x, dir_y).tolist()
        else:
            point_data['dir'] = interp.bilinear_interpolation(
                row['lat'], row['lon'], ds,
                variable=preset.wave_dir_var,
                lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
            ).tolist()

        if preset.wave_dir_spread_var and not preset.wave_dir_spread_use_default:
            point_data['spread'] = interp.bilinear_interpolation(
                row['lat'], row['lon'], ds,
                variable=preset.wave_dir_spread_var,
                lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
            ).tolist()
        else:
            default_spread = preset.wave_dir_spread_default
            point_data['spread'] = [default_spread] * len(point_data['hs'])
            logger.warning(
                "No directional spreading variable — using default "
                f"{default_spread}°"
            )

        n_steps = len(point_data['hs'])
        point_data['time'] = np.arange(0, n_steps * timestep, timestep).tolist()

        clean_name = bnd_name[:-2] if bnd_name.endswith('_a') else bnd_name
        wave_data[clean_name] = point_data

        progress_callback(0.1 + 0.7 * (idx + 1) / total_points)
        if (idx + 1) % max(1, total_points // 5) == 0:
            logger.info(f"  Processed {idx + 1}/{total_points} points")

    _check_cancel(cancel_event)
    logger.info("Wave extraction complete")

    output_files = []
    if model_type == 'delft3d4':
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '.bcw'))
        writer_delft3d4.write_bcw(output_path, wave_data, timing['reference_date'])
        output_files.append(output_path)

        if mdw_template_path:
            mdw_out = str(
                Path(output_dir) / (Path(mdw_template_path).stem + '_updated.mdw')
            )
            writer_mdw.write_mdw(mdw_template_path, bnd_pairs, mdw_out)
            output_files.append(mdw_out)
            logger.info(f"Updated .mdw written: {mdw_out}")
        else:
            logger.warning(
                "No .mdw template provided — skipping .mdw update. The .bcw boundary names will not be linked to a wave model definition."
            )
    else:
        fm_data = {}
        for name, wd in wave_data.items():
            fm_data[name] = {'time': wd['time'], 'values': wd['hs']}
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '_wave.bc'))
        writer_delft3dfm.write_bc(
            output_path, fm_data, 'significantwaveheight', 'm',
            str(timing['start_datetime']),
        )
        output_files.append(output_path)

    progress_callback(1.0)
    logger.info(f"Output written: {', '.join(output_files)}")

    return {
        'output_files': output_files,
        'n_boundaries': len(wave_data),
        'n_timesteps': len(next(iter(wave_data.values()))['hs']),
        'repaired_points': qc_summary['repaired_points'],
        'rejected_boundaries': qc_summary['rejected_boundaries'],
    }


def extract_currents(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: Optional[str] = None,
    data_files: List[str] = None,
    preset_name: str = 'TrilaWatt',
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = 'N',
    output_path: Optional[str] = None,
    output_dir: str = '.',
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Extract depth-averaged currents (u, v) and write boundary file.

    D3D4 → .bct with velocity (u/v) at end A and end B (spec §1.3).
    FM   → .bc with vector quantity 'velocitybnd' at each .pli point.
    """

    if logger is None:
        logger = logging.getLogger('repd3d.bc')

    preset = preset_module.get_bc_preset(preset_name)
    if not preset.current_u_var or not preset.current_v_var:
        raise ValueError(
            f"Preset '{preset.name}' has no current_u_var/current_v_var. Use a preset that supports currents (e.g. TrilaWatt, Copernicus SSH/Currents) or fill in via Custom mode."
        )

    logger.info(f"Starting currents extraction — preset: {preset.name}")
    _check_cancel(cancel_event)
    timing, crs_type = _load_model_definition(model_type, model_def_path, grid_path, logger)
    _check_cancel(cancel_event)
    bnd_loc_geo, _native, _pairs = _load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )

    _check_cancel(cancel_event)
    ds = _open_dataset(data_files, preset, timing, logger)
    _check_cancel(cancel_event)
    qc_context = _build_qc_context(
        ds, [preset.current_u_var, preset.current_v_var], preset
    )
    bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
        model_type,
        'currents',
        ds,
        bnd_loc_geo,
        _pairs,
        [preset.current_u_var, preset.current_v_var],
        preset,
        qc_context,
        logger,
        progress_callback=lambda p: progress_callback(p * 0.1),
        cancel_event=cancel_event,
    )
    progress_callback(0.1)

    u_dict, v_dict = {}, {}
    total = len(bnd_loc_geo)
    for idx, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        _check_cancel(cancel_event)
        u = interp.bilinear_interpolation(
            row['lat'], row['lon'], ds,
            variable=preset.current_u_var,
            lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
        )
        v = interp.bilinear_interpolation(
            row['lat'], row['lon'], ds,
            variable=preset.current_v_var,
            lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
        )
        u_dict[row['boundaries']] = u
        v_dict[row['boundaries']] = v
        progress_callback(0.1 + 0.7 * (idx + 1) / total)
        if (idx + 1) % max(1, total // 5) == 0:
            logger.info(f"  Processed {idx + 1}/{total} points")

    if model_type == 'delft3d4':
        boundary_data = {}
        for key in u_dict:
            bnd_name = key[:-2]
            suffix = key[-1]
            boundary_data.setdefault(bnd_name, {})
            boundary_data[bnd_name][f"u_{suffix}"] = u_dict[key].tolist()
            boundary_data[bnd_name][f"v_{suffix}"] = v_dict[key].tolist()

        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '_currents.bct'))
        writer_delft3d4.write_bct_currents(
            output_path, boundary_data,
            timing['reference_date'],
            timing['tstart'], timing['tstop'], timestep,
        )
    else:
        fm_data = {}
        for key in u_dict:
            n = len(u_dict[key])
            fm_data[key] = {
                'time': list(range(n)),
                'values_x': u_dict[key].tolist(),
                'values_y': v_dict[key].tolist(),
            }
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '_currents.bc'))
        writer_delft3dfm.write_bc(
            output_path, fm_data, 'velocitybnd', 'm/s',
            str(timing['start_datetime']),
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")

    return {
        'output_files': [output_path],
        'n_boundaries': len(u_dict) // (2 if model_type == 'delft3d4' else 1),
        'n_timesteps': len(next(iter(u_dict.values()))),
        'repaired_points': qc_summary['repaired_points'],
        'rejected_boundaries': qc_summary['rejected_boundaries'],
    }


def extract_sediments(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: Optional[str] = None,
    data_files: List[str] = None,
    preset_name: str = 'Copernicus BGC',
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = 'N',
    output_path: Optional[str] = None,
    output_dir: str = '.',
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Extract sediment concentrations and write a .bcc (D3D4) or .bc (FM) file.

    The preset declares one variable per fraction (sediment_conc_vars) plus an
    optional split ratio (sediment_split_ratios). If two fractions reference
    the same NetCDF variable, the split ratios partition the extracted value.
    Spec §3.
    """

    if logger is None:
        logger = logging.getLogger('repd3d.bc')

    preset = preset_module.get_bc_preset(preset_name)
    if not preset.sediment_conc_vars:
        raise ValueError(
            f"Preset '{preset.name}' has no sediment_conc_vars. Use 'Copernicus BGC' or 'Custom' with sediment variables filled in."
        )

    n_frac = len(preset.sediment_conc_vars)
    fraction_names = (
        preset.sediment_fraction_names
        if len(preset.sediment_fraction_names) == n_frac
        else [f"fraction_{i + 1}" for i in range(n_frac)]
    )
    split_ratios = (
        preset.sediment_split_ratios
        if len(preset.sediment_split_ratios) == n_frac
        else [1.0] * n_frac
    )
    unique_vars = list(dict.fromkeys(preset.sediment_conc_vars))

    logger.info(
        f"Starting sediment extraction — preset: {preset.name}, "
        f"{n_frac} fraction(s)"
    )

    _check_cancel(cancel_event)
    timing, crs_type = _load_model_definition(model_type, model_def_path, grid_path, logger)
    _check_cancel(cancel_event)
    bnd_loc_geo, _native, _pairs = _load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )

    _check_cancel(cancel_event)
    ds = _open_dataset(data_files, preset, timing, logger)
    _check_cancel(cancel_event)
    qc_context = _build_qc_context(ds, unique_vars, preset)
    bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
        model_type,
        'sediments',
        ds,
        bnd_loc_geo,
        _pairs,
        unique_vars,
        preset,
        qc_context,
        logger,
        progress_callback=lambda p: progress_callback(p * 0.05),
        cancel_event=cancel_event,
    )
    progress_callback(0.05)

    raw = {}
    total_steps = len(unique_vars) * len(bnd_loc_geo)
    step = 0
    for var_name in unique_vars:
        raw[var_name] = {}
        for _, row in bnd_loc_geo.iterrows():
            _check_cancel(cancel_event)
            arr = interp.bilinear_interpolation(
                row['lat'], row['lon'], ds,
                variable=var_name,
                lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
            )
            raw[var_name][row['boundaries']] = arr
            step += 1
            progress_callback(0.05 + 0.75 * step / total_steps)

    if model_type == 'delft3d4':
        bcc_data = {}
        for f_idx, fname in enumerate(fraction_names):
            src_var = preset.sediment_conc_vars[f_idx]
            ratio = split_ratios[f_idx]
            bcc_data[fname] = {}
            for key in raw[src_var]:
                bnd_name = key[:-2]
                suffix = key[-1]
                bcc_data[fname].setdefault(bnd_name, {})
                bcc_data[fname][bnd_name][f"conc_{suffix}"] = (
                    raw[src_var][key] * ratio
                ).tolist()

        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '.bcc'))
        writer_delft3d4.write_bcc(
            output_path, bcc_data, fraction_names,
            timing['reference_date'],
            timing['tstart'], timing['tstop'], timestep,
        )
    else:
        fm_data = {}
        for f_idx, fname in enumerate(fraction_names):
            src_var = preset.sediment_conc_vars[f_idx]
            ratio = split_ratios[f_idx]
            for key in raw[src_var]:
                n = len(raw[src_var][key])
                fm_data[f"{fname}_{key}"] = {
                    'time': list(range(n)),
                    'values': (raw[src_var][key] * ratio).tolist(),
                }
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '_sediments.bc'))
        writer_delft3dfm.write_bc(
            output_path, fm_data, 'tracerbnd', 'kg/m3',
            str(timing['start_datetime']),
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")

    return {
        'output_files': [output_path],
        'n_boundaries': len(bnd_loc_geo) // (2 if model_type == 'delft3d4' else 1),
        'n_fractions': n_frac,
        'n_timesteps': len(next(iter(next(iter(raw.values())).values()))),
        'repaired_points': qc_summary['repaired_points'],
        'rejected_boundaries': qc_summary['rejected_boundaries'],
    }


def extract_salinity(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: Optional[str] = None,
    data_files: List[str] = None,
    preset_name: str = 'Custom',
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = 'N',
    output_path: Optional[str] = None,
    output_dir: str = '.',
    logger: Optional[logging.Logger] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event: Optional[threading.Event] = None,
) -> dict:
    """
    Extract salinity and write a .bcc (D3D4) or .bc (FM) file.

    Uses the preset's salinity_var field. Output format follows
    bc_file_format_spec.md §3.5 — 'Salinity end A/B uniform' parameters
    in [ppt]. Single variable, no fraction logic.
    """

    if logger is None:
        logger = logging.getLogger('repd3d.bc')

    preset = preset_module.get_bc_preset(preset_name)
    if not preset.salinity_var:
        raise ValueError(
            f"Preset '{preset.name}' has no salinity_var. Use 'Custom' mode and select a salinity variable."
        )

    logger.info(f"Starting salinity extraction — preset: {preset.name}")
    _check_cancel(cancel_event)
    timing, crs_type = _load_model_definition(model_type, model_def_path, grid_path, logger)
    _check_cancel(cancel_event)
    bnd_loc_geo, _native, _pairs = _load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )

    _check_cancel(cancel_event)
    ds = _open_dataset(data_files, preset, timing, logger)
    _check_cancel(cancel_event)
    qc_context = _build_qc_context(ds, [preset.salinity_var], preset)
    bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
        model_type,
        'salinity',
        ds,
        bnd_loc_geo,
        _pairs,
        [preset.salinity_var],
        preset,
        qc_context,
        logger,
        progress_callback=lambda p: progress_callback(p * 0.05),
        cancel_event=cancel_event,
    )
    progress_callback(0.05)

    raw = {}
    total_steps = len(bnd_loc_geo)
    for step_i, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        _check_cancel(cancel_event)
        arr = interp.bilinear_interpolation(
            row['lat'], row['lon'], ds,
            variable=preset.salinity_var,
            lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
        )
        raw[row['boundaries']] = arr
        progress_callback(0.05 + 0.75 * (step_i + 1) / total_steps)

    if model_type == 'delft3d4':
        sal_data = {}
        for key in raw:
            bnd_name = key[:-2]
            suffix = key[-1]
            sal_data.setdefault(bnd_name, {})
            sal_data[bnd_name][f"sal_{suffix}"] = raw[key].tolist()

        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '_salinity.bcc'))
        writer_delft3d4.write_bcc_salinity(
            output_path, sal_data,
            timing['reference_date'],
            timing['tstart'], timing['tstop'], timestep,
        )
    else:
        fm_data = {}
        for key in raw:
            n = len(raw[key])
            fm_data[key] = {
                'time': list(range(n)),
                'values': raw[key].tolist(),
            }
        if output_path is None:
            output_path = str(Path(output_dir) / (Path(model_def_path).stem + '_salinity.bc'))
        writer_delft3dfm.write_bc(
            output_path, fm_data, 'salinitybnd', 'ppt',
            str(timing['start_datetime']),
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")

    return {
        'output_files': [output_path],
        'n_boundaries': len(bnd_loc_geo) // (2 if model_type == 'delft3d4' else 1),
        'n_timesteps': len(next(iter(raw.values()))),
        'repaired_points': qc_summary['repaired_points'],
        'rejected_boundaries': qc_summary['rejected_boundaries'],
    }


def _load_model_definition(model_type, model_def_path, grid_path, logger):
    """Load model definition and extract timing + CRS."""

    logger.info(f"Loading model definition: {Path(model_def_path).name}")

    if model_type == 'delft3d4':
        mdf_data = file_parsers.parse_mdf(model_def_path)
        timing = file_parsers.extract_mdf_timing(mdf_data)
        if grid_path:
            grd = file_parsers.parse_grd(grid_path)
            crs_type = grd['crs_type']
        else:
            crs_type = 'Cartesian'
    else:
        mdu_data = file_parsers.parse_mdu(model_def_path)
        timing = file_parsers.extract_mdu_timing(mdu_data)
        crs_type = 'Spherical'

    logger.info(
        f"Simulation: {timing['start_datetime']} to {timing['end_datetime']} (ref: "
        f"{timing['reference_date_iso']}, CRS: {crs_type})"
    )

    return timing, crs_type


def _load_boundary_locations(model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger):
    """
    Load boundary coordinates and convert to lat/lon for interpolation.

    Returns:
        (bnd_loc_geo, bnd_loc_native, bnd_pairs)
            bnd_loc_geo:    DataFrame with lat/lon/boundaries (interpolation)
            bnd_loc_native: DataFrame with easting/northing/boundary
                            (projected source coords; needed by .mdw writer)
            bnd_pairs:      List of {'name','x_a','y_a','x_b','y_b'} dicts in
                            native coords. Empty for FM.
    """

    logger.info(f"Loading boundaries: {Path(boundary_path).name}")

    bnd_pairs = []

    if model_type == 'delft3d4':
        bnd_indices = file_parsers.parse_bnd(boundary_path)
        bnd_data = file_parsers.extract_boundary_coordinates(grid_path, bnd_indices)

        rows = []
        for bnd in bnd_data:
            rows.append({
                'boundary': f"{bnd['bnd_name']}_a",
                'easting': bnd['x_a'],
                'northing': bnd['y_a'],
            })
            rows.append({
                'boundary': f"{bnd['bnd_name']}_b",
                'easting': bnd['x_b'],
                'northing': bnd['y_b'],
            })
            bnd_pairs.append({
                'name': bnd['bnd_name'],
                'x_a': bnd['x_a'], 'y_a': bnd['y_a'],
                'x_b': bnd['x_b'], 'y_b': bnd['y_b'],
            })

        bnd_loc = pd.DataFrame(rows)
    else:
        pli_data = file_parsers.parse_pli(boundary_path)
        rows = []
        for pli in pli_data:
            # D-Flow FM numbers polyline support points 1-based in the .bc Name
            # (<polyline label>_<0001..>, manual C.5.1.2 example left01_0001);
            # a 0-based index would leave every block unbindable.
            for i, (lon, lat) in enumerate(pli['points']):
                rows.append({
                    'boundary': f"{pli['name']}_{i + 1:04d}",
                    'easting': lon,
                    'northing': lat,
                })

        bnd_loc = pd.DataFrame(rows)

    bnd_loc_geo = coordinate_utils.convert_boundary_coordinates(
        bnd_loc, crs_type, utm_zone, hemisphere
    )

    logger.info(f"Loaded {len(bnd_loc_geo)} boundary points (geographic coords)")
    return bnd_loc_geo, bnd_loc, bnd_pairs


GRID_MATCH_TOLERANCE_DEG = 1e-6   # coordinate agreement needed to call two grids identical
REGRID_MIN_WET_WEIGHT = 0.5       # regridded cell keeps a value only if >= this weight is wet
REGRID_TIME_CHUNK = 16            # timesteps per dask task when a file is regridded: interp holds several target-grid intermediates per task and the threaded scheduler runs many tasks at once (16 -> ~0.4-0.6 GB peak on a 200x200 target, 32 -> ~1 GB, 128 -> ~3 GB; measured 2026-10-15)
LAZY_CHUNK_TARGET_BYTES = 64 * 1024 * 1024   # size of one lazy (dask) chunk along time
MAX_DEDUP_SLICES = 64             # overlapping files: dedup by slicing up to this many runs, else by index


def _describe_grid(path, lat, lon):
    """One log line describing a file's lat/lon grid."""
    def spacing(values):
        return float(np.median(np.abs(np.diff(values)))) if values.size > 1 else float('nan')

    return (
        f"{Path(path).name}: {lon.size} x {lat.size} cells, spacing "
        f"{spacing(lon):.5f} x {spacing(lat):.5f} deg, lon {lon.min():.3f}..{lon.max():.3f}, "
        f"lat {lat.min():.3f}..{lat.max():.3f}"
    )


def _same_grid(lat_a, lon_a, lat_b, lon_b):
    return (
        lat_a.size == lat_b.size
        and lon_a.size == lon_b.size
        and np.allclose(lat_a, lat_b, rtol=0, atol=GRID_MATCH_TOLERANCE_DEG)
        and np.allclose(lon_a, lon_b, rtol=0, atol=GRID_MATCH_TOLERANCE_DEG)
    )


def _regrid_dataset(ds, lat_name, lon_name, target_lat, target_lon, logger):
    """Interpolate every lat/lon field of one source file onto the target axes.

    Everything stays lazy: each dask task interpolates one on-disk chunk (one
    time step for TrilaWatt/EasyGSH), so memory does not depend on file size.
    Floating fields use a mask-normalised bilinear interpolation — dry/missing
    cells are left out of the weights instead of turning their neighbours NaN,
    and a target cell is kept only when at least REGRID_MIN_WET_WEIGHT of its
    weight comes from wet cells. Non-float fields (flags, ids) use nearest.
    """
    targets = {lat_name: target_lat, lon_name: target_lon}
    kept = {}
    dropped = []
    for name, var in ds.data_vars.items():
        has_lat = lat_name in var.dims
        has_lon = lon_name in var.dims
        if not (has_lat or has_lon):
            kept[name] = var
            continue
        if not (has_lat and has_lon):
            dropped.append(name)   # would carry the old axis into the merge
            continue
        if var.chunks is not None:
            # Whole fields per task (interp needs the full axes), but no more
            # than REGRID_TIME_CHUNK steps per task: the interpolation holds
            # a few target-grid-sized intermediates per task, so the lazy
            # time chunking of _open_dataset would be too coarse here.
            var = var.chunk({
                lat_name: -1,
                lon_name: -1,
                **{d: REGRID_TIME_CHUNK for d in var.dims if d not in (lat_name, lon_name)},
            })
        if np.issubdtype(var.dtype, np.floating):
            wet = var.notnull().astype('float32')
            weighted = var.fillna(0.0).interp(targets, method='linear')
            weight = wet.interp(targets, method='linear')
            out = weighted / weight.where(weight > 0)
            out = out.where(weight >= REGRID_MIN_WET_WEIGHT).astype(var.dtype)
        else:
            out = var.interp(targets, method='nearest')
        out.attrs = dict(var.attrs)
        kept[name] = out

    # Build a fresh Dataset so the source file's own lat/lon axes cannot
    # survive into the merge (they would re-create the union grid).
    result = xr.Dataset(kept, attrs=dict(ds.attrs))
    source = Path(str(ds.encoding.get('source', 'dataset'))).name
    logger.warning(
        f"[grid]   regridded {source}: {ds.sizes[lon_name]} x {ds.sizes[lat_name]} "
        f"-> {target_lon.size} x {target_lat.size} cells"
        + (f" (dropped 1-D fields {', '.join(dropped)})" if dropped else "")
    )
    return result


def _grid_harmoniser(data_files, preset, logger):
    """Return an open_mfdataset preprocess that puts every file on one grid, or None.

    Fail-safe for annual files delivered at different resolutions (seen with
    TrilaWatt: tides_2017 at 0.0089 deg next to tides_2018 at 0.0177 deg).
    Merging those as-is is silently destructive: join='outer' builds the union
    of both coordinate sets, so the coarser year has data on one cell in four
    and NaN elsewhere, every bilinear stencil then carries a year-long NaN
    tail (the fill routine extrapolates that into hundreds of metres), and the
    reindex materialises the whole union grid — 46 GB resident for a 12 GB +
    3 GB pair.

    When all files share one grid this returns None and nothing changes.
    Otherwise it logs the grids loudly and returns a preprocess that regrids
    the coarser files onto the finest one (see _regrid_dataset). For nested
    grids (coarse nodes on every n-th fine node) the boundary values come out
    identical to a single-file run on the coarse file; for unrelated grids
    they are smoothed by one coarse cell.
    """
    lat_name, lon_name = preset.lat_coord, preset.lon_coord
    grids = []
    for path in data_files:
        try:
            with xr.open_dataset(path) as probe:
                if lat_name not in probe.coords or lon_name not in probe.coords:
                    return None
                lat = np.asarray(probe[lat_name].values, dtype=float)
                lon = np.asarray(probe[lon_name].values, dtype=float)
        except Exception as exc:   # never let the fail-safe itself break a run
            logger.info(
                f"[grid] Could not read the grid of {Path(path).name} ({exc}); "
                f"merging files as-is."
            )
            return None
        if lat.ndim != 1 or lon.ndim != 1 or lat.size < 2 or lon.size < 2:
            return None
        grids.append((path, lat, lon))

    first_lat, first_lon = grids[0][1], grids[0][2]
    if all(_same_grid(first_lat, first_lon, lat, lon) for _, lat, lon in grids[1:]):
        return None

    def fineness(entry):
        _, lat, lon = entry
        cell = float(np.median(np.abs(np.diff(lat)))) * float(np.median(np.abs(np.diff(lon))))
        return (cell, -(lat.size * lon.size))

    target_path, target_lat, target_lon = min(grids, key=fineness)

    logger.warning("[grid] Source files are on DIFFERENT grids:")
    for path, lat, lon in grids:
        marker = "  <- finest, used as target" if path == target_path else ""
        logger.warning(f"[grid]   {_describe_grid(path, lat, lon)}{marker}")
    logger.warning(
        "[grid] Merging them as-is would leave the coarser file(s) with data on only a "
        "fraction of the merged cells and give every boundary stencil year-long gaps. "
        f"Regridding the other file(s) onto the {Path(target_path).name} grid "
        "(bilinear, dry/missing cells excluded) before merging. Where the grids are "
        "nested this reproduces a single-file run of each file exactly; otherwise the "
        "regridded years are smoothed by one coarse cell. Check the seam at the year "
        "boundary in the output."
    )

    def preprocess(ds):
        lat = np.asarray(ds[lat_name].values, dtype=float)
        lon = np.asarray(ds[lon_name].values, dtype=float)
        if _same_grid(lat, lon, target_lat, target_lon):
            return ds
        return _regrid_dataset(ds, lat_name, lon_name, target_lat, target_lon, logger)

    return preprocess


def _lazy_time_chunks(paths, preset, logger):
    """Dask chunking for the multi-file open: whole 2-D fields, many steps per chunk.

    Without a chunks= argument open_mfdataset adopts the on-disk chunking, and
    TrilaWatt/EasyGSH files hold one 2-D field per chunk, so every lazy array
    carries one dask task per timestep: a 13-month window is ~27650 tasks and
    reading one cell's series costs ~1.8 ms per step of scheduler overhead
    (50 s per corner series, measured 2026-10-15) against ~0.2 ms per step for
    the decompression itself. Chunks of about LAZY_CHUNK_TARGET_BYTES along the
    time axis keep everything just as lazy and bounded, and remove that
    overhead. The chunk length is sized from the LARGEST per-step field over
    all files (mixed-resolution pairs may sort coarse-first). Dimensions other
    than time keep their on-disk chunking. Returns None (= on-disk chunking,
    the old behaviour) if any file cannot be sized.
    """
    time_name = preset.time_coord
    if isinstance(paths, (str, Path)):
        paths = [paths]
    slice_bytes = 0
    stored_steps = 1          # largest on-disk chunk length along time (1 for TrilaWatt/EasyGSH)
    for path in paths:
        try:
            with xr.open_dataset(path) as probe:
                if time_name not in probe.dims:
                    return None
                for var in probe.data_vars.values():
                    if time_name not in var.dims:
                        continue
                    per_step = int(var.size // max(1, int(var.sizes[time_name])))
                    slice_bytes = max(slice_bytes, int(var.dtype.itemsize) * per_step)
                    preferred = var.encoding.get('preferred_chunks')
                    if isinstance(preferred, dict) and time_name in preferred:
                        stored_steps = max(stored_steps, int(preferred[time_name]))
                    elif var.encoding.get('chunksizes'):
                        stored_steps = max(stored_steps, int(var.encoding['chunksizes'][var.dims.index(time_name)]))
        except Exception as exc:   # never let chunk sizing break a run
            logger.info(
                f"[open] Could not size lazy chunks for {Path(path).name} ({exc}); "
                f"using the on-disk chunking."
            )
            return None
    if slice_bytes <= 0:
        return None
    steps = int(max(1, min(4096, LAZY_CHUNK_TARGET_BYTES // slice_bytes)))
    if stored_steps > 1:      # whole stored chunks per dask chunk, never a chunk boundary inside one
        steps = max(stored_steps, (steps // stored_steps) * stored_steps)
    return {time_name: steps}


def _open_dataset(data_files, preset, timing, logger):
    """Open one or more NetCDF files and slice to simulation period.

    The multi-file branch REQUIRES dask, and that is deliberate — do not
    "simplify" it to open_dataset + combine_by_coords. open_mfdataset keeps the
    concatenation lazy, so the .sel() below reads only the simulation window,
    and the downstream per-point bilinear interpolation then touches just the
    handful of cells it needs. An eager combine materialises the entire merged
    grid instead: for a 370-day EasyGSH run that is 17.3 GB resident to serve
    ~73 MB of boundary data, a 236x overshoot that would exhaust a typical
    16 GB user machine.

    Nothing here imports dask statically, so PyInstaller cannot infer it from
    the source. Frozen builds must therefore collect it explicitly (see
    --collect-all dask in auto_py_to_exe_settings.json); without it the app
    builds and runs fine until a user selects a second file and then fails with
    "chunk manager 'dask' is not available".

    Fail-safe: files on different lat/lon grids are detected before the merge
    and the coarser ones are regridded onto the finest (see _grid_harmoniser).
    Files on one grid take exactly the path they always did.
    """

    logger.info(f"Opening {len(data_files)} data file(s)...")

    if len(data_files) == 1:
        ds = xr.open_dataset(data_files[0])
    else:
        files = sorted(data_files)
        ds = xr.open_mfdataset(
            files,
            combine='by_coords',
            join='outer',
            compat='no_conflicts',
            preprocess=_grid_harmoniser(files, preset, logger),
            chunks=_lazy_time_chunks(files, preset, logger),
        )

    start_dt = timing['start_datetime'] + timedelta(hours=preset.time_offset_hours)
    end_dt = timing['end_datetime'] + timedelta(hours=preset.time_offset_hours)

    ds = ds.sel({preset.time_coord: slice(str(start_dt), str(end_dt))})

    if len(data_files) > 1:
        ds = _drop_duplicate_times(ds, preset, logger)

    n_times = ds.sizes.get(preset.time_coord, ds.sizes.get('time', 0))
    logger.info(f"Data sliced: {start_dt} to {end_dt} ({n_times} timesteps)")

    return ds


def _drop_duplicate_times(ds, preset, logger):
    """Remove repeated timestamps produced by source files that overlap in time.

    Annual products are not always disjoint. The TrilaWatt CAU files repeat the
    year boundary — tides_2018 ends at 2019-01-01 00:00 and tides_2019 starts
    there — and the two copies of that instant disagree (measured on
    sea_surface_height_2d: 0.06 m mean, 1.21 m max, 94% of wet cells differing
    by >0.01 m, i.e. comparable to a full 20-minute tidal change). EasyGSH
    annual files are disjoint and unaffected.

    Left in place the extra sample is silently destructive rather than merely
    redundant: the writers pair the data column to a regular arange() time
    column positionally, not by timestamp (writer_delft3d4.write_bct), so one
    surplus sample shifts every later row one timestep late and the overflow is
    dropped by the trailing min(). The output stays well-formed and no error is
    raised, so it would only ever be caught as a phase error in the model.

    The first occurrence of each timestamp is kept. Which copy to keep is
    arbitrary — either choice leaves exactly one discontinuity at the seam — so
    the magnitude is logged for the user to judge.
    """
    time_coord = preset.time_coord
    if time_coord not in ds.dims:
        return ds

    times = ds[time_coord].values
    unique_times, first_index = np.unique(times, return_index=True)
    n_duplicates = times.size - unique_times.size
    if n_duplicates == 0:
        return ds

    repeated = unique_times[np.bincount(
        np.searchsorted(unique_times, times), minlength=unique_times.size) > 1]
    logger.warning(
        f"[overlap] Source files overlap in time: {n_duplicates} duplicate "
        f"timestamp(s) at {', '.join(str(pd.Timestamp(t)) for t in repeated[:5])}"
        f"{' ...' if repeated.size > 5 else ''}. Keeping the first occurrence "
        f"of each — without this the series would be one sample too long and "
        f"every value after the overlap would be written one timestep late."
    )

    # Report how much the duplicated copies actually disagree, so a benign
    # repeat is distinguishable from two genuinely different runs.
    for var_name, var in ds.data_vars.items():
        if time_coord not in var.dims:
            continue
        try:
            for stamp in repeated[:1]:
                positions = np.flatnonzero(times == stamp)
                if positions.size < 2:
                    continue
                first = np.asarray(
                    var.isel({time_coord: int(positions[0])}).values, dtype=float)
                second = np.asarray(
                    var.isel({time_coord: int(positions[1])}).values, dtype=float)
                both = np.isfinite(first) & np.isfinite(second)
                if not both.any():
                    continue
                diff = np.abs(first[both] - second[both])
                logger.warning(
                    f"[overlap]   {var_name} at {pd.Timestamp(stamp)}: copies differ by "
                    f"mean {diff.mean():.4g}, max {diff.max():.4g} (n={both.sum()} cells)"
                )
        except Exception as exc:  # never let a diagnostic break extraction
            logger.info(f"[overlap]   ({var_name}: could not compare copies — {exc})")

    keep = np.sort(first_index)
    # Prefer slices to an index array: the kept indices form a few contiguous
    # runs (one per seam), and slicing keeps the lazy chunk layout intact. An
    # index-array take re-chunks the window (dask's shuffle-take) so every
    # chunk straddles two source chunks — on the regrid path that doubles the
    # interpolation work and memory per read.
    breaks = np.flatnonzero(np.diff(keep) != 1) + 1
    starts = np.concatenate(([0], breaks))
    ends = np.concatenate((breaks, [keep.size]))
    runs = [(int(keep[s]), int(keep[e - 1]) + 1) for s, e in zip(starts, ends)]
    if len(runs) > MAX_DEDUP_SLICES:
        return ds.isel({time_coord: keep})
    parts = [ds.isel({time_coord: slice(a, b)}) for a, b in runs]
    if len(parts) == 1:
        return parts[0]
    return xr.concat(parts, dim=time_coord, data_vars='minimal', coords='minimal', compat='override')


def _assemble_bct_data(output_dict: dict) -> Dict[str, Dict[str, list]]:
    """
    Assemble extracted water levels into writer_delft3d4.write_bct format.

    Input: {'BND_1_a': array, 'BND_1_b': array, 'BND_2_a': array, ...}
    Output: {'BND_1': {'wl_a': list, 'wl_b': list}, 'BND_2': {...}}

    v1 source: bct_generator.py lines 324-337.
    """

    boundary_data = {}
    for key, values in output_dict.items():
        # Key format: 'BND_NAME_a' or 'BND_NAME_b'
        bnd_name = key[:-2]
        suffix = key[-1]

        if bnd_name not in boundary_data:
            boundary_data[bnd_name] = {}

        boundary_data[bnd_name][f"wl_{suffix}"] = values.tolist()

    return boundary_data


def _assemble_fm_data(output_dict: dict) -> Dict[str, Dict[str, list]]:
    """
    Assemble extracted data into writer_delft3dfm.write_bc format.

    Input: {'BND_1_0001': array, 'BND_1_0002': array, ...}
    Output: {'BND_1_0001': {'time': [...], 'values': [...]}, ...}
    """

    boundary_data = {}
    for key, values in output_dict.items():
        n = len(values)
        boundary_data[key] = {
            'time': list(range(n)),
            'values': values.tolist(),
        }

    return boundary_data


def _wave_direction_from_vectors(
    x_component: np.ndarray,
    y_component: np.ndarray,
) -> np.ndarray:
    """
    Convert wave direction from x/y vector components to nautical degrees.

    Uses the quadrant-based atan2 logic from v1 (bcw_generator.py lines 323-346).
    PRESERVED EXACTLY — this is scientifically verified.

    Convention: "wave FROM direction" in nautical degrees (0=N, 90=E, 180=S, 270=W).
    The +180° correction converts from "wave TO" to "wave FROM".

    Args:
        x_component: 1D array of x-direction values (-1 to 1).
        y_component: 1D array of y-direction values (-1 to 1).

    Returns:
        1D array of nautical direction in degrees (0–360).
    """

    directions = np.zeros(len(x_component))

    for i in range(len(x_component)):
        x = float(x_component[i])
        y = float(y_component[i])

        if x > 0 and y > 0:
            rad = math.atan2(x, y)
            deg = math.degrees(rad)
        elif x > 0 and y < 0:
            rad = math.atan2(-y, x)
            deg = math.degrees(rad) + 90
        elif x < 0 and y < 0:
            rad = math.atan2(-x, -y)
            deg = math.degrees(rad) + 180
        elif x < 0 and y > 0:
            rad = math.atan2(y, -x)
            deg = math.degrees(rad) + 270
        else:
            deg = 0.0

        directions[i] = deg

    directions = (directions + 180) % 360
    return directions


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== backend_bc_extractor.py self-test ===\n")

    global all_pass
    all_pass = True

    def check(name, result, expected, tol=0.01):
        global all_pass
        safe_name = str(name).encode('ascii', 'replace').decode('ascii')
        if isinstance(result, (float, np.floating)):
            ok = abs(result - expected) < tol
        else:
            ok = result == expected
        status = 'PASS' if ok else 'FAIL'
        if not ok:
            all_pass = False
            print(f"  {safe_name}: [{status}]  got {result!r}")
        else:
            print(f"  {safe_name}: [{status}]")

    # ── Wave direction from vectors ──────────────────────────────────────
    print("Wave direction from vectors:")

    d = _wave_direction_from_vectors(np.array([0.001, 1.0]), np.array([1.0, 0.001]))

    check('Near-north vector → ~180° (FROM south)', d[0], 180.0, tol=5.0)

    d2 = _wave_direction_from_vectors(np.array([1.0]), np.array([0.001]))
    check('East vector → ~270° (FROM west)', d2[0], 270.0, tol=5.0)

    d3 = _wave_direction_from_vectors(np.array([-1.0]), np.array([0.001]))
    check('West vector → ~90° (FROM east)', d3[0], 90.0, tol=5.0)

    print("\nAssemble BCT data:")
    test_dict = {
        'South_a': np.array([1.0, 2.0]),
        'South_b': np.array([1.1, 2.1]),
        'West_a': np.array([3.0, 4.0]),
        'West_b': np.array([3.1, 4.1]),
    }
    assembled = _assemble_bct_data(test_dict)
    check('2 boundaries', len(assembled), 2)
    check('South has wl_a', 'wl_a' in assembled['South'], True)
    check('South has wl_b', 'wl_b' in assembled['South'], True)
    check('South wl_a[0]', assembled['South']['wl_a'][0], 1.0)
    check('West wl_b[1]', assembled['West']['wl_b'][1], 4.1)

    print("\nAssemble FM data:")
    fm_dict = {
        'BND_0001': np.array([1.0, 2.0, 3.0]),
    }
    fm_assembled = _assemble_fm_data(fm_dict)
    check('FM has values', 'values' in fm_assembled['BND_0001'], True)
    check('FM values length', len(fm_assembled['BND_0001']['values']), 3)

    print("\nTwo-stage QC:")






    class _Preset:
        lat_coord = 'lat'
        lon_coord = 'lon'
        time_coord = 'time'

    preset = _Preset()
    test_logger = logging.getLogger('repd3d.bc.selftest')
    test_logger.handlers = [logging.NullHandler()]
    test_logger.setLevel(logging.CRITICAL)

    times = np.arange(40)
    lats = np.array([0.0, 1.0])
    lons = np.array([0.0, 1.0])

    ds_full = xr.Dataset(
        {'wl': (['time', 'lat', 'lon'], np.ones((40, 2, 2)))},
        coords={'time': times, 'lat': lats, 'lon': lons},
    )

    partial = np.ones((40, 2, 2))
    partial[:20, 0, 0] = np.nan
    ds_partial = xr.Dataset(
        {'wl': (['time', 'lat', 'lon'], partial)},
        coords={'time': times, 'lat': lats, 'lon': lons},
    )

    land = np.ones((40, 2, 2))
    land[:, :, 0] = np.nan
    ds_land = xr.Dataset(
        {'wl': (['time', 'lat', 'lon'], land)},
        coords={'time': times, 'lat': lats, 'lon': lons},
    )

    qc_full = _build_qc_context(ds_full, ['wl'], preset)
    fast_good = _evaluate_point_status(0.25, 0.75, ds_full, ['wl'], preset, qc_full)
    check('offshore point fast-pass', fast_good['path'], 'fast', tol=0)
    check('offshore point reliable', fast_good['reliable'], True, tol=0)
    check(
        'offshore point strict calls',
        qc_full['stats']['strict_point_calls'],
        0,
        tol=0,
    )

    qc_outside = _build_qc_context(ds_full, ['wl'], preset)
    outside = _evaluate_point_status(-0.1, 0.5, ds_full, ['wl'], preset, qc_outside)
    check('outside point reliable', outside['reliable'], False, tol=0)
    check('outside point uses fast path', outside['path'], 'fast', tol=0)

    qc_partial_fail = _build_qc_context(ds_partial, ['wl'], preset)
    fast_partial_fail = _fast_evaluate_point_quality(
        0.0, 0.5, ds_partial, ['wl'], preset, qc_partial_fail
    )
    partial_fail = _evaluate_point_status(
        0.0, 0.5, ds_partial, ['wl'], preset, qc_partial_fail
    )
    strict_partial_fail = _evaluate_point_quality(
        0.0, 0.5, ds_partial, ['wl'], preset, _build_qc_context(ds_partial, ['wl'], preset)
    )
    check('mixed coastal point flagged ambiguous', fast_partial_fail['decision'], 'ambiguous', tol=0)
    check('mixed coastal point strict-fail', partial_fail['reliable'], False, tol=0)
    check('mixed coastal parity fail', partial_fail['reliable'], strict_partial_fail['reliable'], tol=0)

    qc_partial_pass = _build_qc_context(ds_partial, ['wl'], preset)
    fast_partial_pass = _fast_evaluate_point_quality(
        0.0, 0.98, ds_partial, ['wl'], preset, qc_partial_pass
    )
    partial_pass = _evaluate_point_status(
        0.0, 0.98, ds_partial, ['wl'], preset, qc_partial_pass
    )
    strict_partial_pass = _evaluate_point_quality(
        0.0, 0.98, ds_partial, ['wl'], preset, _build_qc_context(ds_partial, ['wl'], preset)
    )
    check('mixed coastal strict-pass still ambiguous first', fast_partial_pass['decision'], 'ambiguous', tol=0)
    check('mixed coastal strict-pass reliable', partial_pass['reliable'], True, tol=0)
    check('mixed coastal strict-pass path', partial_pass['path'], 'strict', tol=0)
    check('mixed coastal parity pass', partial_pass['reliable'], strict_partial_pass['reliable'], tol=0)

    d3d4_repair = pd.DataFrame([
        {'boundaries': 'Edge_a', 'lat': 0.0, 'lon': 0.05},
        {'boundaries': 'Edge_b', 'lat': 0.0, 'lon': 1.0},
    ])

    d3d4_pairs = [{'name': 'Edge', 'x_a': 0.05, 'y_a': 0.0, 'x_b': 1.0, 'y_b': 0.0}]
    repaired_geo, repaired_pairs, repaired_summary = _prepare_boundary_sampling(
        'delft3d4',
        'water level',
        ds_land,
        d3d4_repair,
        d3d4_pairs,
        ['wl'],
        preset,
        _build_qc_context(ds_land, ['wl'], preset),
        test_logger,
    )
    edge_a_lon = float(
        repaired_geo.loc[repaired_geo['boundaries'] == 'Edge_a', 'lon'].iloc[0]
    )

    check('D3D4 repaired one point', len(repaired_summary['repaired_points']), 1, tol=0)
    check('D3D4 rejected none', len(repaired_summary['rejected_boundaries']), 0, tol=0)
    check('D3D4 repair keeps mdw pair', len(repaired_pairs), 1, tol=0)
    check('D3D4 repair uses binary search', edge_a_lon > 0.9525, True, tol=0)

    d3d4_reject = pd.DataFrame([
        {'boundaries': 'Bad_a', 'lat': 0.0, 'lon': 0.05},
        {'boundaries': 'Bad_b', 'lat': 0.0, 'lon': 0.1},
        {'boundaries': 'Good_a', 'lat': 0.0, 'lon': 1.0},
        {'boundaries': 'Good_b', 'lat': 0.0, 'lon': 1.0},
    ])

    reject_geo, _, reject_summary = _prepare_boundary_sampling(
        'delft3d4',
        'water level',
        ds_land,
        d3d4_reject,
        [
            {'name': 'Bad', 'x_a': 0.05, 'y_a': 0.0, 'x_b': 0.1, 'y_b': 0.0},
            {'name': 'Good', 'x_a': 1.0, 'y_a': 0.0, 'x_b': 1.0, 'y_b': 0.0},
        ],
        ['wl'],
        preset,
        _build_qc_context(ds_land, ['wl'], preset),
        test_logger,
    )
    check('D3D4 rejected bad boundary', reject_summary['rejected_boundaries'][0], 'Bad', tol=0)
    check('D3D4 keeps good boundary points', len(reject_geo), 2, tol=0)

    fm_repair = pd.DataFrame([
        {'boundaries': 'Line_0000', 'lat': 0.0, 'lon': 0.05},
        {'boundaries': 'Line_0001', 'lat': 0.0, 'lon': 1.0},
    ])

    fm_geo, _, fm_summary = _prepare_boundary_sampling(
        'delft3dfm',
        'water level',
        ds_land,
        fm_repair,
        [],
        ['wl'],
        preset,
        _build_qc_context(ds_land, ['wl'], preset),
        test_logger,
    )
    fm_lon = float(fm_geo.loc[fm_geo['boundaries'] == 'Line_0000', 'lon'].iloc[0])
    check('FM repaired one point', len(fm_summary['repaired_points']), 1, tol=0)
    check('FM repair moves inland point offshore', fm_lon > 0.9525, True, tol=0)

    fm_error = None
    try:
        _prepare_boundary_sampling(
            'delft3dfm',
            'water level',
            ds_land,
            pd.DataFrame([
                {'boundaries': 'Fail_0000', 'lat': 0.0, 'lon': 0.05},
                {'boundaries': 'Fail_0001', 'lat': 0.0, 'lon': 0.1},
            ]),
            [],
            ['wl'],
            preset,
            _build_qc_context(ds_land, ['wl'], preset),
            test_logger,
        )
    except ValueError as exc:
        fm_error = str(exc)
    check('FM unresolved raises', fm_error is not None, True, tol=0)
    check('FM unresolved mentions bad point', 'Fail_0000' in (fm_error or ''), True, tol=0)

    perf_points = pd.DataFrame([
        {'boundaries': 'Perf_0000', 'lat': 0.25, 'lon': 0.25},
        {'boundaries': 'Perf_0001', 'lat': 0.25, 'lon': 0.75},
        {'boundaries': 'Perf_0002', 'lat': 0.75, 'lon': 0.25},
        {'boundaries': 'Perf_0003', 'lat': 0.75, 'lon': 0.75},
    ])

    perf_context = _build_qc_context(ds_full, ['wl'], preset)
    _prepare_boundary_sampling(
        'delft3dfm',
        'water level',
        ds_full,
        perf_points,
        [],
        ['wl'],
        preset,
        perf_context,
        test_logger,
    )

    legacy_calls = len(perf_points)
    new_calls = perf_context['stats']['strict_point_calls']
    check('strict QC call count reduced', new_calls < legacy_calls, True, tol=0)

    print("\nStats-based stencil summary (strict/intertidal QC without file access):")
    for ds_name, ds_case in (('full', ds_full), ('partial', ds_partial), ('land', ds_land)):
        ctx = _build_qc_context(ds_case, ['wl'], preset)
        ok = 'valid_count' in ctx['metrics']['wl']
        for la, lo in ((0.25, 0.75), (0.0, 0.5), (0.0, 0.98), (0.5, 0.02), (1.0, 1.0), (0.3, 0.6), (1.2, 0.5)):
            ok = ok and _inspect_stencil(la, lo, ds_case, 'wl', preset, ctx) == interp.inspect_bilinear_stencil(
                la, lo, ds_case, 'wl', lat_coord='lat', lon_coord='lon')
        check(f'stats stencil == disk stencil ({ds_name})', ok, True, tol=0)
    check('stats absent -> disk stencil',
          _inspect_stencil(0.25, 0.75, ds_full, 'wl', preset, {'metrics': {'wl': {}}})
          == interp.inspect_bilinear_stencil(0.25, 0.75, ds_full, 'wl'), True, tol=0)

    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print(f"{'=' * 50}")
