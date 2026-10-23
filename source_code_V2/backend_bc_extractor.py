# -*- coding: utf-8 -*-
"""Wrapper around the BC extractor implementation.

Loads backend_bc_extractor_impl.py (reconstructed source, preferred) and
falls back to backend_bc_extractor_impl.pyc (legacy compiled module) if the
source file is absent. Both were verified byte-for-byte equivalent on real
EasyGSH extractions (water level + waves) in July 2026.
"""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime
import importlib.machinery
import importlib.util
import inspect
import io
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple


_THIS_DIR = Path(__file__).resolve().parent
_IMPL_SRC_PATH = _THIS_DIR / "backend_bc_extractor_impl.py"
_IMPL_PYC_PATH = _THIS_DIR / "backend_bc_extractor_impl.pyc"


def _load_impl():
    # Normal import first: this is what PyInstaller/auto-py-to-exe can follow,
    # so the engine gets bundled into frozen builds automatically.
    try:
        import backend_bc_extractor_impl as impl
        return impl
    except ImportError:
        pass

    # Path-based fallbacks for unusual sys.path setups.
    if _IMPL_SRC_PATH.exists():
        impl_path = _IMPL_SRC_PATH
        spec = importlib.util.spec_from_file_location(
            "_backend_bc_extractor_impl",
            str(_IMPL_SRC_PATH),
        )
    elif _IMPL_PYC_PATH.exists():
        impl_path = _IMPL_PYC_PATH
        loader = importlib.machinery.SourcelessFileLoader(
            "_backend_bc_extractor_impl",
            str(_IMPL_PYC_PATH),
        )
        spec = importlib.util.spec_from_loader(loader.name, loader)
    else:
        raise FileNotFoundError(
            "Backend implementation not found: "
            f"{_IMPL_SRC_PATH} or {_IMPL_PYC_PATH}"
        )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load backend implementation spec from {impl_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_IMPL = _load_impl()

_ORIG_PREPARE_BOUNDARY_SAMPLING = _IMPL._prepare_boundary_sampling

INTERTIDAL_MIN_VALID_FRACTION = 0.20
INTERTIDAL_MAX_EDGE_FRACTION = 0.10

# Physical plausibility bounds per preset variable field. Interpolation never
# extrapolates any more, so a sampled series outside these bounds can only come
# from the source data itself (wrong datum, units or fill value) — stop before
# it reaches a model. Variables not listed (wave direction components,
# spreading) are not range-checked.
PLAUSIBLE_RANGES: Dict[str, Tuple[float, float]] = {
    "water_level_var": (-15.0, 15.0),     # m, any datum
    "current_u_var": (-10.0, 10.0),       # m/s
    "current_v_var": (-10.0, 10.0),       # m/s
    "salinity_var": (0.0, 50.0),          # ppt / PSU
    "wave_hs_var": (0.0, 30.0),           # m
    "wave_tp_var": (0.0, 40.0),           # s
}
PLAUSIBLE_SEDIMENT_RANGE: Tuple[float, float] = (0.0, 1.0e4)   # any concentration unit; catches negatives/fill values
WAVE_TIMING_SOURCE_AUTO = "auto"
WAVE_TIMING_SOURCE_MANUAL = "manual"
_MANUAL_WAVE_DATE_FORMATS = ("%Y%m%d", "%Y-%m-%d")
_MANUAL_WAVE_DATETIME_FORMATS = (
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%dT%H:%M:%S",
    "%Y-%m-%dT%H:%M",
    "%Y%m%d %H:%M:%S",
    "%Y%m%d %H:%M",
)

_intertidal_tolerance_enabled = False


def _clone(value: Any) -> Any:
    return deepcopy(value)


def set_intertidal_tolerance(enabled: bool) -> None:
    global _intertidal_tolerance_enabled
    _intertidal_tolerance_enabled = enabled


def _point_result_cache_key(lat: float, lon: float, variables: List[str], allow_intertidal: bool = True) -> Tuple[Any, ...]:
    # The verdict depends on whether the intertidal path may be used, so the
    # flag is part of the key: a Neumann probe screened strictly must not reuse
    # a verdict cached for the same coordinate as an intertidal-tolerant
    # water-level endpoint (shared section corners make this a real case).
    return (_IMPL._point_cache_key(lat, lon), tuple(variables), bool(allow_intertidal))


class _SpacedQCLogger:
    def __init__(self, logger):
        self._logger = logger

    def info(self, msg, *args, **kwargs):
        if isinstance(msg, str) and msg.startswith("[qc] ") and not msg.startswith("[qc]   "):
            self._logger.info("", *args, **kwargs)
        self._logger.info(msg, *args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._logger, name)


def _format_intertidal_reasons(details: List[Dict[str, Any]]) -> str:
    if not details:
        return "intertidal tolerance"
    parts = []
    for detail in details:
        parts.append(
            f"{detail['variable']}/{detail['corner']} ({detail['weight']:.1%} weight): "
            f"valid {detail['valid_fraction']:.1%}, {detail['internal_nan_count']} internal gaps"
        )
    return "; ".join(parts)


def _format_fast_assessment_summary(variable: str, assessment: Dict[str, Any]) -> str:
    """Format fast-QC failures defensively without assuming every key exists."""
    issues = assessment.get("issues", [])
    if not issues:
        return f"{variable}: insufficient valid source coverage"
    if "bad_weight" in assessment:
        return _IMPL._assessment_summary(variable, assessment)
    return f"{variable}: {'; '.join(issues[:3])}"


def _register_intertidal_acceptance(qc_context: Dict[str, Any], lat: float, lon: float, details: List[Dict[str, Any]]) -> None:
    if details:
        qc_context.setdefault("_intertidal_accepts", {})[_IMPL._point_cache_key(lat, lon)] = _clone(details)


def _check_intertidal_tolerance(lat: float, lon: float, ds, variables: List[str], preset, qc_context=None) -> Dict[str, Any]:
    details: List[Dict[str, Any]] = []
    inspect_stencil = getattr(_IMPL, "_inspect_stencil", None)
    derived = set((qc_context or {}).get("derived_vars", ()))
    for variable in variables:
        # Per-corner gap statistics come from the QC context's one-pass
        # per-cell statistics (no file access); without them (or with a
        # legacy implementation module) the stencil is read from disk.
        if inspect_stencil is not None:
            stencil = inspect_stencil(lat, lon, ds, variable, preset, qc_context)
        else:
            stencil = _IMPL.interp.inspect_bilinear_stencil(
                lat, lon, ds, variable=variable, lat_coord=preset.lat_coord, lon_coord=preset.lon_coord,
            )
        if stencil["outside"]:
            return {"reliable": False, "issues": ["outside dataset extent"]}
        if variable in derived:
            if not _derived_stencil_has_value(stencil["corners"]):
                return {"reliable": False, "issues": [f"{variable}: {DERIVED_ISSUE}"]}
            continue
        for corner_name, corner in stencil["corners"].items():
            weight = corner["normalized_weight"]
            if weight <= 0.0:
                continue
            n_steps = corner["n_steps"]
            if n_steps == 0:
                return {"reliable": False, "issues": [f"{corner_name}: no timesteps"]}
            if not corner["all_nan"] and corner["valid_fraction"] >= _IMPL.QC_MIN_VALID_FRACTION:
                edge_frac = (corner["leading_nan_count"] + corner["trailing_nan_count"]) / n_steps
                if edge_frac <= INTERTIDAL_MAX_EDGE_FRACTION:
                    continue
            if corner["all_nan"]:
                return {"reliable": False, "issues": [f"{corner_name} ({weight:.1%} weight): all-NaN (land)"]}
            if corner["valid_fraction"] < INTERTIDAL_MIN_VALID_FRACTION:
                return {
                    "reliable": False,
                    "issues": [
                        f"{corner_name} ({weight:.1%} weight): valid fraction "
                        f"{corner['valid_fraction']:.1%} below intertidal minimum "
                        f"{INTERTIDAL_MIN_VALID_FRACTION:.0%}"
                    ],
                }
            edge_frac = (corner["leading_nan_count"] + corner["trailing_nan_count"]) / n_steps
            if edge_frac > INTERTIDAL_MAX_EDGE_FRACTION:
                return {
                    "reliable": False,
                    "issues": [
                        f"{corner_name} ({weight:.1%} weight): edge gap {edge_frac:.1%} - "
                        "data does not span simulation period"
                    ],
                }
            details.append(
                {
                    "variable": variable,
                    "corner": corner_name,
                    "weight": weight,
                    "valid_fraction": corner["valid_fraction"],
                    "internal_nan_count": corner["internal_nan_count"],
                }
            )
    return {"reliable": True, "path": "intertidal", "issues": [], "details": details}


DERIVED_ISSUE = "no source value at any timestep in any contributing corner"


def _wave_derived_variables(preset, variables) -> Tuple[str, ...]:
    """
    Wave quantities that the source leaves undefined wherever there are no
    waves: peak period and directional spreading. TrilaWatt sets them NaN
    exactly when Hs <= 0.05 m (checked 2026-10-15 on the 2017/2018 files:
    never while Hs > 0.1 m, and only 11 of 11639 wet cells are gap-free over
    13 months). Applying the strict gap rules to them rejected 20 of 90 wave
    sections for calm weather. They are QC'd for "has a value at some
    timestep in some contributing corner" only; the calm-period gaps are
    bridged by the mask-normalised interpolation at extraction and reported
    in the [fill] log lines.
    """
    return tuple(v for v in (preset.wave_tp_var, preset.wave_dir_spread_var) if v and v in variables)


def _log_derived_policy(logger, derived_vars) -> None:
    if derived_vars:
        logger.info(
            f"[qc] {', '.join(derived_vars)}: derived wave quantities (undefined while Hs ~ 0). QC only "
            f"requires a value at some timestep; calm-period gaps are bridged at extraction (see [fill] lines)."
        )


def _derived_stencil_has_value(corners: Dict[str, Dict[str, Any]]) -> bool:
    return any(float(c["normalized_weight"]) > 0.0 and not c["all_nan"] for c in corners.values())


def _fast_assess_derived_quality(variable: str, geometry: Dict[str, Any], masks: Dict[str, Any]) -> Dict[str, Any]:
    """Fast-path verdict for a derived quantity: a value somewhere in the stencil suffices."""
    any_valid_mask = masks["any_valid_mask"]
    lat_idx1, lat_idx2 = geometry["lat_indices"]
    lon_idx1, lon_idx2 = geometry["lon_indices"]
    corners = {"Q11": (lat_idx1, lon_idx1), "Q12": (lat_idx1, lon_idx2), "Q21": (lat_idx2, lon_idx1), "Q22": (lat_idx2, lon_idx2)}
    contributing = [name for name in corners if float(geometry["normalized_weights"][name]) > 0.0]
    has_value = any(bool(any_valid_mask[corners[name]]) for name in contributing)
    ok = has_value and not geometry["outside"]
    return {
        "decision": "accept" if ok else "hard_fail",
        "issues": [] if ok else [DERIVED_ISSUE],
        "bad_weight": 0.0 if ok else 1.0,
        "max_bad_corner_weight": 0.0 if ok else 1.0,
        "contributing_corner_count": len(contributing),
        "full_valid_contributing_count": 0,
        "all_contributing_full": False,
    }


def _fast_assess_variable_quality(variable: str, geometry: Dict[str, Any], qc_context: Dict[str, Any]) -> Dict[str, Any]:
    masks = qc_context["metrics"][variable]
    if variable in qc_context.get("derived_vars", ()):
        return _fast_assess_derived_quality(variable, geometry, masks)
    full_valid_mask = masks["full_valid_mask"]
    any_valid_mask = masks["any_valid_mask"]
    lat_idx1, lat_idx2 = geometry["lat_indices"]
    lon_idx1, lon_idx2 = geometry["lon_indices"]
    indices = {"Q11": (lat_idx1, lon_idx1), "Q12": (lat_idx1, lon_idx2), "Q21": (lat_idx2, lon_idx1), "Q22": (lat_idx2, lon_idx2)}

    contributing_corner_count = 0
    full_valid_contributing_count = 0
    all_contributing_full = True
    bad_weight = 0.0
    max_bad_corner_weight = 0.0
    issues: List[str] = []

    for corner_name, (lat_idx, lon_idx) in indices.items():
        weight = float(geometry["normalized_weights"][corner_name])
        if weight <= 0.0:
            continue
        contributing_corner_count += 1
        full_valid = bool(full_valid_mask[lat_idx, lon_idx])
        any_valid = bool(any_valid_mask[lat_idx, lon_idx])
        if full_valid:
            full_valid_contributing_count += 1
            continue
        all_contributing_full = False
        if any_valid:
            continue
        bad_weight += weight
        max_bad_corner_weight = max(max_bad_corner_weight, weight)
        issues.append(f"{corner_name} ({weight:.1%} weight): all-NaN")

    # Fast acceptance ONLY when every contributing corner is valid at every
    # timestep. The former "3-of-4 override" (accept when three corners are
    # fully valid, whatever the fourth looks like) let a 51 %-weight corner
    # that was wet for 11 storm timesteps through unchecked; the strict
    # assessment below applies the gap/weight limits to such corners instead.
    if all_contributing_full:
        return {
            "decision": "accept",
            "issues": [],
            "bad_weight": bad_weight,
            "max_bad_corner_weight": max_bad_corner_weight,
            "contributing_corner_count": contributing_corner_count,
            "full_valid_contributing_count": full_valid_contributing_count,
            "all_contributing_full": all_contributing_full,
        }
    if geometry["outside"] or bad_weight > _IMPL.QC_MAX_BAD_WEIGHT or max_bad_corner_weight > _IMPL.QC_MAX_SINGLE_BAD_CORNER_WEIGHT:
        return {
            "decision": "hard_fail",
            "issues": issues or ["outside dataset extent"],
            "bad_weight": bad_weight,
            "max_bad_corner_weight": max_bad_corner_weight,
            "contributing_corner_count": contributing_corner_count,
            "full_valid_contributing_count": full_valid_contributing_count,
            "all_contributing_full": all_contributing_full,
        }
    return {
        "decision": "ambiguous",
        "issues": issues,
        "bad_weight": bad_weight,
        "max_bad_corner_weight": max_bad_corner_weight,
        "contributing_corner_count": contributing_corner_count,
        "full_valid_contributing_count": full_valid_contributing_count,
        "all_contributing_full": all_contributing_full,
    }


def _fast_evaluate_point_quality(lat: float, lon: float, ds, variables: List[str], preset, qc_context: Dict[str, Any]) -> Dict[str, Any]:
    geometry = _IMPL.interp.inspect_bilinear_geometry(lat, lon, ds, lat_coord=preset.lat_coord, lon_coord=preset.lon_coord)
    if geometry["outside"]:
        return {"decision": "hard_fail", "assessments": {}, "issues": ["outside dataset extent"]}

    assessments: Dict[str, Dict[str, Any]] = {}
    issues: List[str] = []
    saw_hard_fail = False

    for variable in variables:
        assessment = _fast_assess_variable_quality(variable, geometry, qc_context)
        assessments[variable] = assessment
        if assessment["decision"] == "accept":
            continue
        issues.extend(
            _format_fast_assessment_summary(variable, assessment).split(" | ")
            if assessment["issues"] else []
        )
        saw_hard_fail = saw_hard_fail or assessment["decision"] == "hard_fail"

    if all(assessment["decision"] == "accept" for assessment in assessments.values()):
        return {"decision": "accept", "assessments": assessments, "issues": []}
    if saw_hard_fail:
        return {"decision": "hard_fail", "assessments": assessments, "issues": issues}
    return {"decision": "ambiguous", "assessments": assessments, "issues": issues}


def _evaluate_point_status(
    lat: float,
    lon: float,
    ds,
    variables: List[str],
    preset,
    qc_context: Dict[str, Any],
    allow_intertidal: Optional[bool] = None,
) -> Dict[str, Any]:
    """
    Fast QC, then strict QC for ambiguous stencils, then (water-level style
    boundaries only) the intertidal tolerance. ``allow_intertidal=None`` follows
    the GUI setting stored in the qc_context; ``False`` forces strict-only
    screening, which the Neumann probe loop uses because a probe on a drying
    flat cannot feed a gradient fit.
    """
    if allow_intertidal is None:
        allow_intertidal = bool(qc_context.get("_intertidal_tolerance", _intertidal_tolerance_enabled))
    cache = qc_context["point_result_cache"]
    key = _point_result_cache_key(lat, lon, variables, allow_intertidal)
    if key in cache:
        return _clone(cache[key])

    fast = _fast_evaluate_point_quality(lat, lon, ds, variables, preset, qc_context)
    if fast["decision"] == "accept":
        result = {"reliable": True, "path": "fast", "assessments": fast["assessments"], "issues": []}
    else:
        result = None
        strict = None
        if fast["decision"] == "ambiguous":
            strict = _IMPL._evaluate_point_quality(lat, lon, ds, variables, preset, qc_context)
            if strict["reliable"]:
                result = {"reliable": True, "path": "strict", "assessments": strict["assessments"], "issues": []}
        if result is None and allow_intertidal:
            intertidal = _check_intertidal_tolerance(lat, lon, ds, variables, preset, qc_context)
            if intertidal["reliable"]:
                _register_intertidal_acceptance(qc_context, lat, lon, intertidal["details"])
                result = {"reliable": True, "path": "intertidal", "assessments": {}, "issues": []}
        if result is None:
            source = strict if strict is not None else fast
            result = {
                "reliable": False,
                "path": "strict" if strict is not None else "fast",
                "assessments": source.get("assessments", {}),
                "issues": source.get("issues", []),
            }

    cache[key] = _clone(result)
    return _clone(result)


def _prepare_boundary_sampling(model_type, extraction_label, ds, bnd_loc_geo, bnd_pairs, variables, preset, qc_context, logger, waves_use_a_only=False, progress_callback=lambda x: None, cancel_event=None):
    spaced_logger = _SpacedQCLogger(logger)
    qc_context.pop("_intertidal_accepts", None)
    qc_context["_intertidal_tolerance"] = _intertidal_tolerance_enabled
    sample_geo, updated_pairs, summary = _ORIG_PREPARE_BOUNDARY_SAMPLING(
        model_type, extraction_label, ds, bnd_loc_geo, bnd_pairs, variables, preset, qc_context, spaced_logger,
        waves_use_a_only=waves_use_a_only, progress_callback=progress_callback, cancel_event=cancel_event,
    )
    for registry_name, prefix, formatter in (
        ("_intertidal_accepts", "accepted by intertidal tolerance", _format_intertidal_reasons),
    ):
        registry = qc_context.get(registry_name, {})
        names = []
        seen = set()
        for _, row in bnd_loc_geo.iterrows():
            point_name = row.get("boundaries", row.get("name"))
            if point_name in seen:
                continue
            seen.add(point_name)
            details = registry.get(_IMPL._point_cache_key(float(row["lat"]), float(row["lon"])))
            if not details:
                continue
            names.append(point_name)
            spaced_logger.info(f"[qc] {point_name} {prefix} for {extraction_label}. Reason: {formatter(details)}")
        if registry_name == "_intertidal_accepts":
            summary["intertidal_points"] = names
    summary.setdefault("intertidal_points", [])
    return sample_geo, updated_pairs, summary


def _normalize_wave_timing_source(value: Any) -> str:
    return WAVE_TIMING_SOURCE_MANUAL if str(value or "").strip().lower() == WAVE_TIMING_SOURCE_MANUAL else WAVE_TIMING_SOURCE_AUTO


def _parse_manual_wave_value(raw_value: Any, label: str, formats: Tuple[str, ...]) -> datetime:
    value = str(raw_value or "").strip()
    if not value:
        raise ValueError(f"{label} is required for manual wave timing.")
    for fmt in formats:
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            continue
    raise ValueError(f"{label} '{value}' has an invalid format. Expected one of: {', '.join(formats)}.")


def _build_manual_wave_timing(reference_date_value: Any, start_datetime_value: Any, end_datetime_value: Any) -> Dict[str, Any]:
    reference_midnight = _parse_manual_wave_value(reference_date_value, "Wave reference date", _MANUAL_WAVE_DATE_FORMATS)
    start_dt = _parse_manual_wave_value(start_datetime_value, "Wave simulation start", _MANUAL_WAVE_DATETIME_FORMATS)
    end_dt = _parse_manual_wave_value(end_datetime_value, "Wave simulation end", _MANUAL_WAVE_DATETIME_FORMATS)
    if start_dt < reference_midnight:
        raise ValueError("Wave simulation start must be on or after the wave reference date.")
    if end_dt < start_dt:
        raise ValueError("Wave simulation end must be on or after the wave simulation start.")
    return {
        "reference_date": reference_midnight.strftime("%Y%m%d"),
        "reference_date_iso": reference_midnight.strftime("%Y-%m-%d"),
        "tstart": (start_dt - reference_midnight).total_seconds() / 60.0,
        "tstop": (end_dt - reference_midnight).total_seconds() / 60.0,
        "start_datetime": start_dt,
        "end_datetime": end_dt,
    }


def _resolve_wave_timing_and_crs(model_type: str, model_def_path: str, grid_path: str | None, logger, wave_timing_source: Any = WAVE_TIMING_SOURCE_AUTO, wave_reference_date: Any = None, wave_start_datetime: Any = None, wave_end_datetime: Any = None) -> Tuple[Dict[str, Any], str]:
    if _normalize_wave_timing_source(wave_timing_source) == WAVE_TIMING_SOURCE_MANUAL:
        if model_type != "delft3d4":
            raise ValueError("Manual wave timing is only supported for Delft3D-4 wave extraction.")
        if not grid_path:
            raise ValueError("Delft3D-4 wave extraction requires a .grd file when manual timing is used.")
        logger.info(f"Loading grid definition: {Path(grid_path).name}")
        crs_type = _IMPL.file_parsers.parse_grd(grid_path)["crs_type"]
        timing = _build_manual_wave_timing(wave_reference_date, wave_start_datetime, wave_end_datetime)
        logger.info(f"Simulation: {timing['start_datetime']} to {timing['end_datetime']} (ref: {timing['reference_date_iso']}, CRS: {crs_type})")
        return timing, crs_type
    if not model_def_path:
        raise ValueError("Auto wave timing requires a model definition file (.mdf/.mdu).")
    return _IMPL._load_model_definition(model_type, model_def_path, grid_path, logger)


def _build_wave_time_values(timing: Dict[str, Any], n_steps: int, timestep: float) -> List[float]:
    if n_steps <= 0:
        return []
    return (float(timing["tstart"]) + (_IMPL.np.arange(n_steps, dtype=float) * float(timestep))).tolist()


def _build_wave_output_minutes(timing: Dict[str, Any], timestep: float) -> _IMPL.np.ndarray:
    """Build the requested output timeline in minutes since the reference date."""
    step = float(timestep)
    if step <= 0:
        raise ValueError("Wave timestep must be positive.")
    tstart = float(timing["tstart"])
    tstop = float(timing["tstop"])
    return _IMPL.np.arange(tstart, tstop + step, step, dtype=float)


def _source_times_on_model_clock(ds, preset) -> _IMPL.pd.DatetimeIndex:
    """Source timestamps moved onto the model clock.

    _open_dataset slices the source at model time + preset.time_offset_hours
    (EasyGSH: -12 h, as v1 did), so the sliced samples still carry source-clock
    stamps. The same offset must be undone here: the resamplers below align by
    timestamp, so without it they cancel the correction (output 12 h early, as
    every v2 build up to 20261015 wrote for EasyGSH) and np.interp clamps the
    last |offset| hours of the window to the final sample.
    """
    time_values = _IMPL.pd.to_datetime(ds[preset.time_coord].values)
    offset_hours = float(getattr(preset, "time_offset_hours", 0.0) or 0.0)
    return time_values - _IMPL.pd.Timedelta(hours=offset_hours)


def _native_time_minutes(ds, preset, timing: Dict[str, Any]) -> _IMPL.np.ndarray:
    """Convert the sliced dataset time coordinate to minutes since the reference date."""
    time_values = _source_times_on_model_clock(ds, preset)
    ref_dt = _IMPL.pd.Timestamp(timing["reference_date_iso"])
    return ((time_values - ref_dt) / _IMPL.pd.Timedelta(minutes=1)).to_numpy(dtype=float)


def _time_axis_reference_datetime(model_type: str, timing: Dict[str, Any]) -> _IMPL.pd.Timestamp:
    if model_type == "delft3d4":
        return _IMPL.pd.Timestamp(timing["reference_date_iso"])
    return _IMPL.pd.Timestamp(timing["start_datetime"])


def _build_output_minutes(model_type: str, timing: Dict[str, Any], timestep: float) -> _IMPL.np.ndarray:
    step = float(timestep)
    if step <= 0.0:
        raise ValueError("Output timestep must be positive.")
    if model_type == "delft3d4":
        start_minutes = float(timing["tstart"])
        stop_minutes = float(timing["tstop"])
    else:
        start_dt = _IMPL.pd.Timestamp(timing["start_datetime"])
        end_dt = _IMPL.pd.Timestamp(timing["end_datetime"])
        start_minutes = 0.0
        stop_minutes = float((end_dt - start_dt) / _IMPL.pd.Timedelta(minutes=1))
    return _IMPL.np.arange(start_minutes, stop_minutes + step, step, dtype=float)


def _native_time_minutes_for_model(model_type: str, ds, preset, timing: Dict[str, Any]) -> _IMPL.np.ndarray:
    time_values = _source_times_on_model_clock(ds, preset)
    ref_dt = _time_axis_reference_datetime(model_type, timing)
    return ((time_values - ref_dt) / _IMPL.pd.Timedelta(minutes=1)).to_numpy(dtype=float)


def _series_matches_output_time(native_minutes: _IMPL.np.ndarray, output_minutes: _IMPL.np.ndarray) -> bool:
    return len(native_minutes) == len(output_minutes) and _IMPL.np.allclose(native_minutes, output_minutes, atol=1e-6)


def _can_subsample_wave_series(native_minutes: _IMPL.np.ndarray, output_minutes: _IMPL.np.ndarray) -> bool:
    if len(native_minutes) < 2 or len(output_minutes) == 0:
        return False
    native_step = native_minutes[1] - native_minutes[0]
    if native_step <= 0:
        return False
    requested_step = output_minutes[1] - output_minutes[0] if len(output_minutes) > 1 else native_step
    ratio = requested_step / native_step
    if not _IMPL.np.isclose(ratio, round(ratio), atol=1e-8):
        return False
    step_multiple = int(round(ratio))
    if step_multiple < 1:
        return False
    offset_ratio = (output_minutes[0] - native_minutes[0]) / native_step
    if not _IMPL.np.isclose(offset_ratio, round(offset_ratio), atol=1e-8):
        return False
    nearest_indices = _IMPL.np.round((output_minutes - native_minutes[0]) / native_step).astype(int)
    if (nearest_indices < 0).any() or (nearest_indices >= len(native_minutes)).any():
        return False
    return _IMPL.np.allclose(native_minutes[nearest_indices], output_minutes, atol=1e-6)


def _resample_linear_series(native_minutes: _IMPL.np.ndarray, values: Any, output_minutes: _IMPL.np.ndarray) -> _IMPL.np.ndarray:
    series = _IMPL.np.asarray(values, dtype=float)
    if len(series) != len(native_minutes):
        raise ValueError("Wave series length does not match native time coordinate length.")
    if len(output_minutes) == 0:
        return _IMPL.np.asarray([], dtype=float)
    if _series_matches_output_time(native_minutes, output_minutes):
        return series.copy()
    if _can_subsample_wave_series(native_minutes, output_minutes):
        native_step = native_minutes[1] - native_minutes[0]
        indices = _IMPL.np.round((output_minutes - native_minutes[0]) / native_step).astype(int)
        return series[indices]
    return _IMPL.np.interp(output_minutes, native_minutes, series)


def _resample_circular_degrees(native_minutes: _IMPL.np.ndarray, values: Any, output_minutes: _IMPL.np.ndarray) -> _IMPL.np.ndarray:
    degrees = _IMPL.np.asarray(values, dtype=float)
    radians = _IMPL.np.deg2rad(degrees)
    sin_resampled = _resample_linear_series(native_minutes, _IMPL.np.sin(radians), output_minutes)
    cos_resampled = _resample_linear_series(native_minutes, _IMPL.np.cos(radians), output_minutes)
    return (_IMPL.np.degrees(_IMPL.np.arctan2(sin_resampled, cos_resampled)) + 360.0) % 360.0


def _default_output_path(output_path: str | None, output_dir: str, model_def_path: str, suffix: str) -> str:
    if output_path is not None:
        return output_path
    return str(Path(output_dir) / f"{Path(model_def_path).stem}{suffix}")


def _point_label(row) -> str:
    name = row.get("boundaries") if hasattr(row, "get") else None
    if name is None or (isinstance(name, float) and _IMPL.np.isnan(name)):
        try:
            return f"point ({float(row['lat']):.5f}, {float(row['lon']):.5f})"
        except Exception:
            return "point"
    return str(name)


def _plausible_range_for(preset, variable: str) -> Optional[Tuple[float, float]]:
    for field_name, bounds in PLAUSIBLE_RANGES.items():
        if getattr(preset, field_name, None) == variable:
            return bounds
    if variable in (getattr(preset, "sediment_conc_vars", None) or []):
        return PLAUSIBLE_SEDIMENT_RANGE
    return None


def _check_plausible_range(label: str, variable: str, values: Any, preset) -> None:
    """Abort on physically impossible sampled values (see PLAUSIBLE_RANGES)."""
    bounds = _plausible_range_for(preset, variable)
    if bounds is None:
        return
    arr = _IMPL.np.asarray(values, dtype=float)
    finite = arr[_IMPL.np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError(f"{label}: '{variable}' has no finite value after sampling.")
    lo, hi = float(finite.min()), float(finite.max())
    if lo < bounds[0] or hi > bounds[1]:
        raise ValueError(
            f"{label}: sampled '{variable}' spans {lo:.4g}..{hi:.4g}, outside the plausible "
            f"range {bounds[0]:g}..{bounds[1]:g}. Sampling never extrapolates, so this comes "
            "from the source data itself (datum, units or an unmasked fill value) — "
            "inspect the file before writing boundary conditions."
        )


def _log_fill_stats(logger, label: str, variable: str, stats: Dict[str, int]) -> None:
    """Report, per point, how a partially dry stencil was handled."""
    renormalised = int(stats.get("renormalised_steps", 0))
    all_dry = int(stats.get("all_dry_steps", 0))
    if not renormalised and not all_dry:
        return
    msg = (
        f"[fill] {label} '{variable}': {renormalised} timestep(s) had a dry/missing stencil "
        "corner (weights renormalised over the wet corners)"
    )
    if all_dry:
        held = int(stats.get("held_leading", 0)) + int(stats.get("held_trailing", 0))
        bridged = int(stats.get("bridged_internal", 0))
        msg += (
            f"; {all_dry} timestep(s) had every corner dry: {bridged} bridged linearly in time, "
            f"{held} held at the nearest valid value. Nothing is extrapolated."
        )
        logger.warning(msg)
    else:
        logger.info(msg)


def _resample_scalar_variable(
    row,
    ds,
    variable: str,
    preset,
    native_minutes: _IMPL.np.ndarray,
    output_minutes: _IMPL.np.ndarray,
    logger=None,
    cache=None,
) -> List[float]:
    stats: Dict[str, int] = {}
    native_series = _IMPL.interp.bilinear_interpolation(
        row["lat"],
        row["lon"],
        ds,
        variable=variable,
        lat_coord=preset.lat_coord,
        lon_coord=preset.lon_coord,
        stats=stats,
        cache=cache,
    )
    label = _point_label(row)
    _check_plausible_range(label, variable, native_series, preset)
    if logger is not None:
        _log_fill_stats(logger, label, variable, stats)
    return _resample_linear_series(native_minutes, native_series, output_minutes).tolist()


def _qc_context(ds, variables, preset, logger=None, cancel_event=None):
    """
    _IMPL._build_qc_context with the logger/cancel keywords when the
    implementation accepts them; a legacy compiled implementation module
    (backend_bc_extractor_impl.pyc fallback) takes the three positional
    arguments only.
    """
    builder = _IMPL._build_qc_context
    try:
        params = inspect.signature(builder).parameters
        accepts = "cancel_event" in params or any(
            p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()
        )
    except (TypeError, ValueError):
        accepts = False
    if accepts:
        return builder(ds, variables, preset, logger=logger, cancel_event=cancel_event)
    return builder(ds, variables, preset)


def _sampling_points(rows) -> List[Tuple[float, float]]:
    """(lat, lon) tuples of a boundary-location DataFrame, in extraction order."""
    return [(float(row["lat"]), float(row["lon"])) for _, row in rows.iterrows()]


def _prefetch_corner_series(ds, points, variables, preset, logger, cancel_event=None):
    """
    CornerSeriesCache for the extraction loop (see interpolation.CornerSeriesCache).

    Reading each point's stencil corners one by one decompresses the whole
    window once per corner on one-field-per-chunk files; the cache reads the
    corner cells of all points in one pass per variable instead. Returns None
    when there is nothing to prefetch; the loop then reads per point.
    """
    points = list(points)
    variables = [v for v in dict.fromkeys(variables) if v]
    if not points or not variables:
        return None
    try:
        cache = _IMPL.interp.CornerSeriesCache(
            ds, points, variables, preset.lat_coord, preset.lon_coord, preset.time_coord,
            cancel_check=lambda: _IMPL._check_cancel(cancel_event), logger=logger,
        )
    except Exception as exc:   # an optimisation must never break an extraction
        logger.info(f"[prefetch] Not prefetching stencil series ({exc}); reading per point instead.")
        return None
    if not any(cache.has(v) for v in variables):
        logger.warning(
            f"[prefetch] None of {', '.join(variables)} is laid out as "
            f"({preset.time_coord}, {preset.lat_coord}, {preset.lon_coord}); reading per point instead (slow)."
        )
        return None
    logger.info(
        f"[prefetch] {len(points)} point(s) x {len(variables)} variable(s) = {cache.n_cells} stencil cell(s); "
        f"read in {cache.n_batches} pass(es) per variable instead of {4 * len(points)} reads per variable."
    )
    return cache


def _prefetch_source_caches(sources, rows, point_sources, logger, cancel_event=None):
    """One CornerSeriesCache per source group, planned from the points that group serves."""
    caches = []
    for src_idx, src in enumerate(sources):
        served = rows[rows["boundaries"].map(lambda name: point_sources.get(name) == src_idx)]
        caches.append(
            _prefetch_corner_series(src["ds"], _sampling_points(served), src["variables"], src["preset"], logger, cancel_event)
            if len(served) else None
        )
    return caches


def _log_progress(logger, idx: int, total: int, noun: str = "points") -> None:
    if total <= 0:
        return
    if (idx + 1) % max(1, total // 5) == 0:
        logger.info(f"  Processed {idx + 1}/{total} {noun}")


# ─── Multi-source (cross-border) extraction ──────────────────────────────────
#
# A source group is one {"data_files": [...], "preset_name": "..."} entry.
# When extract_* receives several groups, every boundary point is QC-evaluated
# against each group on its own native grid — the files are never merged into
# one dataset, so disjoint or phase-offset regional grids (e.g. the TrilaWatt
# Germany and Denmark cutouts) cannot poison each other with union-grid NaNs.
# Each point is served by the first group (in list order) that passes QC, and
# the choice is logged per point in a coverage report.

JUNCTION_PAIR_MAX_M = 300.0        # endpoints of different sections this close form a junction
SEAM_WARN_THRESHOLD = 0.10         # junction series differing by more than this (variable units) get a warning
REPAIR_DISTANCE_ADVISORY_M = 1000.0


def _quantity_variables(preset, quantity: str) -> List[str]:
    """Variable names one source must provide for a given extraction quantity."""
    if quantity == "water_level":
        if not preset.water_level_var:
            raise ValueError(f"Preset '{preset.name}' has no water_level_var.")
        return [preset.water_level_var]
    if quantity == "currents":
        if not preset.current_u_var or not preset.current_v_var:
            raise ValueError(
                f"Preset '{preset.name}' has no current_u_var/current_v_var. "
                "Use a preset that supports currents (e.g. TrilaWatt, Copernicus SSH/Currents) or fill in via Custom mode."
            )
        return [preset.current_u_var, preset.current_v_var]
    if quantity == "salinity":
        if not preset.salinity_var:
            raise ValueError(
                f"Preset '{preset.name}' has no salinity_var. Use 'Custom' mode and select a salinity variable."
            )
        return [preset.salinity_var]
    if quantity == "sediments":
        if not preset.sediment_conc_vars:
            raise ValueError(
                f"Preset '{preset.name}' has no sediment_conc_vars. Use 'Copernicus BGC' or 'Custom' with sediment variables filled in."
            )
        return list(dict.fromkeys(preset.sediment_conc_vars))
    if quantity == "waves":
        variables = [preset.wave_hs_var, preset.wave_tp_var]
        if preset.wave_dir_format == "vector":
            variables.extend([preset.wave_dir_x_var, preset.wave_dir_y_var])
        else:
            variables.append(preset.wave_dir_var)
        if preset.wave_dir_spread_var and not preset.wave_dir_spread_use_default:
            variables.append(preset.wave_dir_spread_var)
        if not all(variables):
            raise ValueError(f"Preset '{preset.name}' is missing wave variables required for wave extraction.")
        return variables
    raise ValueError(f"Unknown extraction quantity: {quantity}")


def _build_boundary_sources(
    source_groups: List[Dict[str, Any]],
    quantity: str,
    model_type: str,
    timing: Dict[str, Any],
    logger,
    cancel_event=None,
) -> List[Dict[str, Any]]:
    """Open each source group and precompute its QC context and time axis."""
    if not source_groups:
        raise ValueError("source_groups must contain at least one {'data_files', 'preset_name'} group.")
    sources: List[Dict[str, Any]] = []
    for index, group in enumerate(source_groups):
        data_files = list(group.get("data_files") or [])
        preset_name = group.get("preset_name")
        if not data_files or not preset_name:
            raise ValueError(f"Source group {index + 1} needs both 'data_files' and 'preset_name'.")
        # Clone so per-group overrides (Custom mappings) never mutate the
        # shared preset singletons in presets.py.
        preset = _clone(_IMPL.preset_module.get_bc_preset(preset_name))
        for field_name, field_value in (group.get("preset_overrides") or {}).items():
            if not hasattr(preset, field_name):
                raise ValueError(
                    f"Source group {index + 1}: unknown preset field '{field_name}' in preset_overrides."
                )
            setattr(preset, field_name, field_value)
        label = group.get("label") or f"{index + 1}:{preset.name}/{Path(data_files[0]).stem}"
        variables = _quantity_variables(preset, quantity)
        logger.info(f"[src] Opening source {label}...")
        ds = _IMPL._open_dataset(data_files, preset, timing, logger)
        _IMPL._check_cancel(cancel_event)
        qc_context = _qc_context(ds, variables, preset, logger, cancel_event)
        qc_context["_intertidal_tolerance"] = _intertidal_tolerance_enabled
        if quantity == "waves":
            qc_context["derived_vars"] = _wave_derived_variables(preset, variables)
            _log_derived_policy(logger, qc_context["derived_vars"])
        if quantity == "waves":
            native_minutes = _native_time_minutes(ds, preset, timing)
        else:
            native_minutes = _native_time_minutes_for_model(model_type, ds, preset, timing)
        sources.append(
            {
                "label": label,
                "preset": preset,
                "variables": variables,
                "ds": ds,
                "qc_context": qc_context,
                "native_minutes": native_minutes,
            }
        )
    if quantity == "sediments":
        n_frac = len(sources[0]["preset"].sediment_conc_vars)
        for source in sources[1:]:
            if len(source["preset"].sediment_conc_vars) != n_frac:
                raise ValueError(
                    "All source groups must define the same number of sediment fractions "
                    f"({sources[0]['label']} has {n_frac}, {source['label']} has "
                    f"{len(source['preset'].sediment_conc_vars)})."
                )
    return sources


def _junction_pairs(
    original: Dict[str, Tuple[float, float]],
    waves_use_a_only: bool = False,
) -> List[Tuple[str, str, float]]:
    """Pairs of endpoints from different sections that are close enough to form a corner."""
    names = [n for n in original if not waves_use_a_only or n.endswith("_a")]
    pairs = []
    for i, name_1 in enumerate(names):
        lat_1, lon_1 = original[name_1]
        for name_2 in names[i + 1:]:
            if name_2[:-2] == name_1[:-2]:
                continue
            lat_2, lon_2 = original[name_2]
            distance_m = _IMPL._distance_meters(lat_1, lon_1, lat_2, lon_2)
            if distance_m <= JUNCTION_PAIR_MAX_M:
                pairs.append((name_1, name_2, distance_m))
    return pairs


def _status_verdict(status: Dict[str, Any]) -> str:
    if status["reliable"]:
        return {"fast": "ok", "strict": "ok*", "intertidal": "tidal"}.get(status.get("path"), "ok")
    if any("outside" in issue for issue in status.get("issues", [])):
        return "out"
    return "FAIL"


def _log_coverage_report(
    extraction_label: str,
    original: Dict[str, Tuple[float, float]],
    statuses: Dict[str, List[Dict[str, Any]]],
    assignments: Dict[str, int],
    labels: List[str],
    logger,
) -> None:
    logger.info(f"[src] Coverage report for {extraction_label} (ok=fast QC, ok*=strict, tidal=intertidal, out=outside grid):")
    col = max(10, max((len(lbl) for lbl in labels), default=10) + 2)
    header = "point".ljust(16) + "".join(lbl.ljust(col) for lbl in labels) + "-> serves"
    logger.info(f"[src]   {header}")
    for point_name in original:
        verdicts = "".join(_status_verdict(s).ljust(col) for s in statuses[point_name])
        if point_name in assignments:
            serves = labels[assignments[point_name]]
            logger.info(f"[src]   {point_name.ljust(16)}{verdicts}-> {serves}")
        else:
            logger.warning(f"[src]   {point_name.ljust(16)}{verdicts}-> NO SOURCE")


def _prepare_multi_source_sampling(
    model_type: str,
    extraction_label: str,
    sources: List[Dict[str, Any]],
    bnd_loc_geo,
    bnd_pairs: List[Dict[str, Any]],
    logger,
    waves_use_a_only: bool = False,
    max_repair_distance_m: Optional[float] = None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event=None,
) -> Tuple[Any, List[Dict[str, Any]], Dict[str, Any], Dict[str, int]]:
    """
    Multi-source counterpart of _prepare_boundary_sampling: assign every
    boundary point to the first source group that passes QC on its own grid,
    repair unassigned points within a source toward their reliable partner,
    and reject (D3D4) or abort (FM) with a coverage report when no source can
    serve a point. Returns the per-point source assignment alongside the usual
    (sample_geo, pairs, summary) triple.
    """
    spaced_logger = _SpacedQCLogger(logger)
    sample_geo = bnd_loc_geo.copy()
    mdw_pairs = list(bnd_pairs)
    labels = [source["label"] for source in sources]
    summary: Dict[str, Any] = {
        "repaired_points": [],
        "rejected_boundaries": [],
        "fast_pass_points": 0,
        "strict_pass_points": 0,
        "intertidal_points": [],
        "source_assignments": {},
        "seam_pairs": [],
    }
    original = {
        row["boundaries"]: (float(row["lat"]), float(row["lon"]))
        for _, row in bnd_loc_geo.iterrows()
    }
    if not original:
        return sample_geo, mdw_pairs, summary, {}

    spaced_logger.info(
        "[src] Source priority for "
        f"{extraction_label}: " + "  ".join(f"{i + 1}) {lbl}" for i, lbl in enumerate(labels))
    )

    statuses: Dict[str, List[Dict[str, Any]]] = {}
    assignments: Dict[str, int] = {}
    for point_idx, (point_name, (lat, lon)) in enumerate(original.items()):
        _IMPL._check_cancel(cancel_event)
        per_source = [
            _evaluate_point_status(lat, lon, src["ds"], src["variables"], src["preset"], src["qc_context"])
            for src in sources
        ]
        statuses[point_name] = per_source
        for src_idx, status in enumerate(per_source):
            if status["reliable"]:
                assignments[point_name] = src_idx
                break
        progress_callback(0.8 * (point_idx + 1) / len(original))

    def _try_repair(point_name: str, targets_by_source: List[List[Tuple[str, float, float]]]):
        """First source (in priority order) whose repair search succeeds within the cap."""
        lat, lon = original[point_name]
        for src_idx, src in enumerate(sources):
            targets = targets_by_source[src_idx]
            if not targets:
                continue
            candidate = _IMPL._repair_point_toward_targets(
                lat, lon, targets, src["ds"], src["variables"], src["preset"], src["qc_context"]
            )
            if candidate is None:
                continue
            if max_repair_distance_m is not None and candidate["distance_m"] > max_repair_distance_m:
                spaced_logger.warning(
                    f"[qc] {point_name}: repair in {labels[src_idx]} needs "
                    f"{candidate['distance_m']:.0f} m (> cap {max_repair_distance_m:.0f} m) - not applied."
                )
                continue
            return candidate, src_idx
        return None, None

    def _apply_repair(point_name: str, repair: Dict[str, Any], src_idx: int) -> None:
        assignments[point_name] = src_idx
        sample_geo.loc[sample_geo["boundaries"] == point_name, ["lat", "lon"]] = [repair["lat"], repair["lon"]]
        msg = (
            f"{point_name} moved {repair['distance_m']:.0f} m toward {repair['target_name']} "
            f"in {labels[src_idx]} for {extraction_label}"
        )
        summary["repaired_points"].append(msg)
        reasons = statuses[point_name][src_idx]["issues"][:2]
        spaced_logger.warning(f"[qc] {msg}. Reason: {' | '.join(reasons) or 'no source passed QC at the original location'}")
        if repair["distance_m"] > REPAIR_DISTANCE_ADVISORY_M:
            spaced_logger.warning(
                f"[qc] {point_name}: relocation exceeds the {REPAIR_DISTANCE_ADVISORY_M:.0f} m advisory distance - "
                "review the seam consistency report before trusting this boundary."
            )

    if model_type == "delft3d4":
        boundary_names = sorted({name[:-2] for name in original})
        active_suffixes = ("a",) if waves_use_a_only else ("a", "b")
        rejected = set()

        for boundary_name in boundary_names:
            _IMPL._check_cancel(cancel_event)
            failing = []
            for suffix in active_suffixes:
                point_name = f"{boundary_name}_{suffix}"
                if point_name not in original or point_name in assignments:
                    continue
                other_name = f"{boundary_name}_{'b' if suffix == 'a' else 'a'}"
                targets_by_source: List[List[Tuple[str, float, float]]] = []
                for src_idx in range(len(sources)):
                    if other_name in original and statuses[other_name][src_idx]["reliable"]:
                        targets_by_source.append([(other_name, *original[other_name])])
                    else:
                        targets_by_source.append([])
                repair, repair_src = _try_repair(point_name, targets_by_source)
                if repair is None:
                    failing.append(point_name)
                else:
                    _apply_repair(point_name, repair, repair_src)

            if failing:
                rejected.add(boundary_name)
                summary["rejected_boundaries"].append(boundary_name)
                spaced_logger.warning(
                    f"[qc] {extraction_label}: rejected boundary {boundary_name} because "
                    f"{', '.join(failing)} have no reliable data in any source."
                )
                for point_name in failing:
                    for src_idx, status in enumerate(statuses[point_name]):
                        issues = " | ".join(status["issues"][:2]) or "unreliable"
                        spaced_logger.warning(f"[qc]   {point_name} vs {labels[src_idx]}: {issues}")

        # Junction guard: near-coincident endpoints of different sections must not
        # jump between sources when a single source can serve both — a cross-source
        # step at a shared corner shows up in the model as an artificial
        # along-boundary gradient (spurious velocities near the corner).
        for name_1, name_2, distance_m in _junction_pairs(original, waves_use_a_only):
            if name_1 not in assignments or name_2 not in assignments:
                continue
            if name_1[:-2] in rejected or name_2[:-2] in rejected:
                continue
            idx_1, idx_2 = assignments[name_1], assignments[name_2]
            if idx_1 == idx_2:
                continue
            common = next(
                (
                    i
                    for i in range(len(sources))
                    if statuses[name_1][i]["reliable"] and statuses[name_2][i]["reliable"]
                ),
                None,
            )
            if common is not None:
                assignments[name_1] = assignments[name_2] = common
                spaced_logger.info(
                    f"[src] junction {name_1}/{name_2} ({distance_m:.0f} m apart): both endpoints pinned to "
                    f"{labels[common]} to avoid a cross-source jump at the corner."
                )
            else:
                summary["seam_pairs"].append((name_1, name_2))
                spaced_logger.warning(
                    f"[src] junction {name_1}/{name_2} ({distance_m:.0f} m apart) must be served by different "
                    f"sources ({labels[idx_1]} vs {labels[idx_2]}). Check the seam consistency report - a large "
                    "difference here creates an artificial along-boundary gradient in the model."
                )

        if not waves_use_a_only:
            for boundary_name in boundary_names:
                if boundary_name in rejected:
                    continue
                name_a, name_b = f"{boundary_name}_a", f"{boundary_name}_b"
                if (
                    name_a in assignments
                    and name_b in assignments
                    and assignments[name_a] != assignments[name_b]
                ):
                    spaced_logger.info(
                        f"[src] section {boundary_name} straddles sources: end a <- {labels[assignments[name_a]]}, "
                        f"end b <- {labels[assignments[name_b]]} (Delft3D interpolates linearly along the section)."
                    )

        _log_coverage_report(extraction_label, original, statuses, assignments, labels, spaced_logger)

        if rejected:
            sample_geo = sample_geo[~sample_geo["boundaries"].str[:-2].isin(rejected)].reset_index(drop=True)
            if waves_use_a_only:
                mdw_pairs = [pair for pair in mdw_pairs if pair["name"] not in rejected]

    else:
        grouped: Dict[str, List[Tuple[int, str]]] = {}
        for point_name in original:
            base_name, idx_str = point_name.rsplit("_", 1)
            grouped.setdefault(base_name, []).append((int(idx_str), point_name))

        unresolved = []
        for base_name, points in grouped.items():
            ordered = sorted(points)
            for pos, (_, point_name) in enumerate(ordered):
                _IMPL._check_cancel(cancel_event)
                if point_name in assignments:
                    continue
                neighbour_names = []
                if pos > 0:
                    neighbour_names.append(ordered[pos - 1][1])
                if pos + 1 < len(ordered):
                    neighbour_names.append(ordered[pos + 1][1])
                targets_by_source = [
                    [
                        (neighbour, *original[neighbour])
                        for neighbour in neighbour_names
                        if statuses[neighbour][src_idx]["reliable"]
                    ]
                    for src_idx in range(len(sources))
                ]
                repair, repair_src = _try_repair(point_name, targets_by_source)
                if repair is None:
                    unresolved.append(point_name)
                else:
                    _apply_repair(point_name, repair, repair_src)

        _log_coverage_report(extraction_label, original, statuses, assignments, labels, spaced_logger)

        if unresolved:
            bad_names = ", ".join(unresolved[:6]) + (", ..." if len(unresolved) > 6 else "")
            raise ValueError(
                f"One or more FM boundary points have no reliable data in any source group for "
                f"{extraction_label}: {bad_names}. Extend a source's coverage or adjust the .pli geometry."
            )

        previous_idx = None
        for base_name, points in grouped.items():
            for _, point_name in sorted(points):
                current_idx = assignments.get(point_name)
                if previous_idx is not None and current_idx != previous_idx:
                    spaced_logger.info(f"[src] source switches to {labels[current_idx]} at {point_name}.")
                previous_idx = current_idx

    for point_name, src_idx in assignments.items():
        status = statuses[point_name][src_idx]
        if not status["reliable"]:
            continue
        if status["path"] == "fast":
            summary["fast_pass_points"] += 1
        elif status["path"] == "strict":
            summary["strict_pass_points"] += 1
        elif status["path"] == "intertidal":
            summary["intertidal_points"].append(point_name)

    for src_idx, source in enumerate(sources):
        for registry_name, prefix, formatter in (
            ("_intertidal_accepts", "accepted by intertidal tolerance", _format_intertidal_reasons),
        ):
            registry = source["qc_context"].get(registry_name, {})
            for point_name, assigned_idx in assignments.items():
                if assigned_idx != src_idx or point_name not in original:
                    continue
                details = registry.get(_IMPL._point_cache_key(*original[point_name]))
                if details:
                    spaced_logger.info(
                        f"[qc] {point_name} {prefix} in {labels[src_idx]} for {extraction_label}. "
                        f"Reason: {formatter(details)}"
                    )

    summary["source_assignments"] = {name: labels[idx] for name, idx in assignments.items()}

    if sample_geo.empty:
        raise ValueError(f"No reliable boundary points remain for {extraction_label}. Extraction aborted.")

    if summary["repaired_points"] or summary["rejected_boundaries"]:
        spaced_logger.warning(
            f"[qc] {extraction_label}: {len(summary['repaired_points'])} point(s) repaired, "
            f"{len(summary['rejected_boundaries'])} boundary section(s) rejected."
        )
    else:
        spaced_logger.info(f"QC passed for {extraction_label} - no repairs needed.")

    spaced_logger.info(
        f"[qc] {extraction_label}: fast-pass {summary['fast_pass_points']}, strict-pass "
        f"{summary['strict_pass_points']}, repaired {len(summary['repaired_points'])}, rejected "
        f"{len(summary['rejected_boundaries'])}"
    )
    progress_callback(1.0)
    return sample_geo, mdw_pairs, summary, assignments


def _report_seam_consistency(
    extraction_label: str,
    series_by_point: Dict[str, List[float]],
    original: Dict[str, Tuple[float, float]],
    assignments: Dict[str, int],
    labels: List[str],
    logger,
    threshold: float = SEAM_WARN_THRESHOLD,
    waves_use_a_only: bool = False,
) -> List[Tuple[str, str, float]]:
    """
    Compare extracted series at junction corners (near-coincident endpoints of
    adjacent sections). A relocation or a cross-source assignment can leave the
    two series at one physical corner disagreeing — which the model would turn
    into an artificial along-boundary gradient — so quantify and report it.
    """
    flagged: List[Tuple[str, str, float]] = []
    checked = 0
    worst = 0.0
    for name_1, name_2, distance_m in _junction_pairs(original, waves_use_a_only):
        series_1 = series_by_point.get(name_1)
        series_2 = series_by_point.get(name_2)
        if series_1 is None or series_2 is None:
            continue
        arr_1 = _IMPL.np.asarray(series_1, dtype=float)
        arr_2 = _IMPL.np.asarray(series_2, dtype=float)
        n = min(len(arr_1), len(arr_2))
        if n == 0:
            continue
        diff = arr_1[:n] - arr_2[:n]
        valid = ~_IMPL.np.isnan(diff)
        if not valid.any():
            continue
        max_diff = float(_IMPL.np.abs(diff[valid]).max())
        rms = float(_IMPL.np.sqrt((diff[valid] ** 2).mean()))
        checked += 1
        worst = max(worst, max_diff)
        idx_1, idx_2 = assignments.get(name_1), assignments.get(name_2)
        source_note = (
            f" [{labels[idx_1]} vs {labels[idx_2]}]"
            if idx_1 is not None and idx_2 is not None and idx_1 != idx_2
            else ""
        )
        if max_diff > threshold:
            flagged.append((name_1, name_2, max_diff))
            logger.warning(
                f"[seam] {name_1} vs {name_2} ({distance_m:.0f} m apart{source_note}): "
                f"max diff {max_diff:.3f}, RMS {rms:.3f} - exceeds the {threshold:g} threshold. "
                "This corner will show an artificial along-boundary gradient; consider extending the "
                "higher-priority source or moving the section break."
            )
        else:
            logger.info(
                f"[seam] {name_1} vs {name_2} ({distance_m:.0f} m apart{source_note}): "
                f"max diff {max_diff:.3f}, RMS {rms:.3f} - OK"
            )
    if checked:
        logger.info(
            f"[seam] {extraction_label}: {checked} junction(s) checked, "
            f"{len(flagged)} above threshold, worst max diff {worst:.3f}."
        )
    return flagged


def neumann_section_names(model_type: str, boundary_path: str) -> set:
    """
    Names of .bnd sections typed 'N' (Neumann water-level gradient).

    Delft3D-4 only: the FM .pli carries no type letter, so FM always returns
    an empty set.
    """
    if model_type != "delft3d4":
        return set()
    return {
        b["bnd_name"]
        for b in _IMPL.file_parsers.parse_bnd(boundary_path)
        if str(b.get("bnd_type", "")).upper() == "N"
    }


def _drop_neumann_points(bnd_loc_geo, bnd_pairs, neumann_names: set, logger):
    """
    Remove N-typed sections from a boundary point set before water-level
    sampling.

    Delft3D reads a table's meaning from the .bnd type letter, not from the
    .bct parameter label, so a water-level table written for an N section is
    interpreted as a gradient: metres where ~1e-6 is expected. Sampling them
    here and letting them reach the writer is the dangerous outcome, so they
    are dropped loudly instead.
    """
    if not neumann_names:
        return bnd_loc_geo, bnd_pairs

    keep = bnd_loc_geo["boundaries"].map(
        lambda name: str(name)[:-2] not in neumann_names
    )
    n_dropped = int((~keep).sum())
    if not n_dropped:
        return bnd_loc_geo, bnd_pairs

    logger.warning(
        f"[water level] Skipping {len(neumann_names)} Neumann section(s) typed 'N' "
        f"in the .bnd ({n_dropped} boundary point(s)): "
        f"{', '.join(sorted(neumann_names))}."
    )
    logger.warning(
        "[water level] Delft3D reads a boundary's meaning from the .bnd type letter, "
        "so a water-level table on an N section would be applied as a gradient "
        "(metres instead of ~1e-6). Generate these sections with the "
        "'WL Gradient (Neumann)' extraction type, which can write the whole "
        "mixed .bct in one run."
    )
    filtered_pairs = [p for p in (bnd_pairs or []) if p["name"] not in neumann_names]
    return bnd_loc_geo[keep].reset_index(drop=True), filtered_pairs


def extract_water_level(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: str | None = None,
    data_files: List[str] = None,
    preset_name: str = "EasyGSH",
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = "N",
    output_path: str | None = None,
    output_dir: str = ".",
    source_groups: Optional[List[Dict[str, Any]]] = None,
    max_repair_distance_m: Optional[float] = None,
    logger=None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event=None,
) -> Dict[str, Any]:
    if logger is None:
        logger = _IMPL.logging.getLogger("repd3d.bc")
    multi_source = bool(source_groups)
    if multi_source:
        logger.info(f"Starting water level extraction — {len(source_groups)} source group(s)")
    else:
        preset = _IMPL.preset_module.get_bc_preset(preset_name)
        logger.info(f"Starting water level extraction — preset: {preset.name}")
    _IMPL._check_cancel(cancel_event)
    timing, crs_type = _IMPL._load_model_definition(model_type, model_def_path, grid_path, logger)
    _IMPL._check_cancel(cancel_event)
    bnd_loc_geo, _bnd_native, _bnd_pairs = _IMPL._load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )
    bnd_loc_geo, _bnd_pairs = _drop_neumann_points(
        bnd_loc_geo, _bnd_pairs,
        neumann_section_names(model_type, boundary_path), logger,
    )
    if len(bnd_loc_geo) == 0:
        raise ValueError(
            "Every section in the .bnd is typed 'N' (Neumann) — there are no "
            "water-level boundaries to extract. Use the 'WL Gradient (Neumann)' "
            "extraction type instead."
        )
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        sources = _build_boundary_sources(source_groups, "water_level", model_type, timing, logger, cancel_event=cancel_event)
        original_coords = {
            row["boundaries"]: (float(row["lat"]), float(row["lon"]))
            for _, row in bnd_loc_geo.iterrows()
        }
        bnd_loc_geo, _unused_pairs, qc_summary, point_sources = _prepare_multi_source_sampling(
            model_type,
            "water level",
            sources,
            bnd_loc_geo,
            _bnd_pairs,
            logger,
            max_repair_distance_m=max_repair_distance_m,
            progress_callback=lambda p: progress_callback(0.1 * p),
            cancel_event=cancel_event,
        )
    else:
        ds = _IMPL._open_dataset(data_files, preset, timing, logger)
        _IMPL._check_cancel(cancel_event)
        qc_context = _qc_context(ds, [preset.water_level_var], preset, logger, cancel_event)
        bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
            model_type,
            "water level",
            ds,
            bnd_loc_geo,
            _bnd_pairs,
            [preset.water_level_var],
            preset,
            qc_context,
            logger,
            progress_callback=lambda p: progress_callback(0.1 * p),
            cancel_event=cancel_event,
        )
    progress_callback(0.1)
    if not multi_source:
        native_minutes = _native_time_minutes_for_model(model_type, ds, preset, timing)
    output_minutes = _build_output_minutes(model_type, timing, timestep)
    output_time_list = output_minutes.tolist()
    logger.info(f"Extracting water level for {len(bnd_loc_geo)} boundary points...")
    if multi_source:
        caches = _prefetch_source_caches(sources, bnd_loc_geo, point_sources, logger, cancel_event)
    else:
        corner_cache = _prefetch_corner_series(
            ds, _sampling_points(bnd_loc_geo), [preset.water_level_var], preset, logger, cancel_event
        )
    output_dict: Dict[str, List[float]] = {}
    total_points = len(bnd_loc_geo)
    for idx, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        _IMPL._check_cancel(cancel_event)
        if multi_source:
            src_idx = point_sources[row["boundaries"]]
            src = sources[src_idx]
            output_dict[row["boundaries"]] = _resample_scalar_variable(
                row, src["ds"], src["preset"].water_level_var, src["preset"],
                src["native_minutes"], output_minutes, logger=logger, cache=caches[src_idx],
            )
        else:
            output_dict[row["boundaries"]] = _resample_scalar_variable(
                row, ds, preset.water_level_var, preset, native_minutes, output_minutes,
                logger=logger, cache=corner_cache,
            )
        progress_callback(0.1 + 0.7 * (idx + 1) / max(1, total_points))
        _log_progress(logger, idx, total_points)
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        _report_seam_consistency(
            "water level", output_dict, original_coords, point_sources,
            [s["label"] for s in sources], logger,
        )
    logger.info("Water level extraction complete")

    if model_type == "delft3d4":
        boundary_data = _IMPL._assemble_bct_data(
            {key: _IMPL.np.asarray(values, dtype=float) for key, values in output_dict.items()}
        )
        output_path = _default_output_path(output_path, output_dir, model_def_path, ".bct")
        _IMPL.writer_delft3d4.write_bct(
            output_path,
            boundary_data,
            timing["reference_date"],
            timing["tstart"],
            timing["tstop"],
            timestep,
        )
    else:
        boundary_data = {
            key: {"time": output_time_list, "values": values}
            for key, values in output_dict.items()
        }
        output_path = _default_output_path(output_path, output_dir, model_def_path, ".bc")
        _IMPL.writer_delft3dfm.write_bc(
            output_path,
            boundary_data,
            "waterlevelbnd",
            "m",
            str(timing["start_datetime"]),
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")
    result = {
        "output_files": [output_path],
        "n_boundaries": len(boundary_data),
        "n_timesteps": len(output_time_list),
        "repaired_points": qc_summary["repaired_points"],
        "rejected_boundaries": qc_summary["rejected_boundaries"],
    }
    if multi_source:
        result["source_assignments"] = qc_summary["source_assignments"]
        result["intertidal_points"] = qc_summary.get("intertidal_points", [])
    return result


def extract_currents(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: str | None = None,
    data_files: List[str] = None,
    preset_name: str = "TrilaWatt",
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = "N",
    output_path: str | None = None,
    output_dir: str = ".",
    source_groups: Optional[List[Dict[str, Any]]] = None,
    max_repair_distance_m: Optional[float] = None,
    logger=None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event=None,
) -> Dict[str, Any]:
    if logger is None:
        logger = _IMPL.logging.getLogger("repd3d.bc")
    multi_source = bool(source_groups)
    if multi_source:
        logger.info(f"Starting currents extraction — {len(source_groups)} source group(s)")
    else:
        preset = _IMPL.preset_module.get_bc_preset(preset_name)
        if not preset.current_u_var or not preset.current_v_var:
            raise ValueError(
                f"Preset '{preset.name}' has no current_u_var/current_v_var. "
                "Use a preset that supports currents (e.g. TrilaWatt, Copernicus SSH/Currents) or fill in via Custom mode."
            )
        logger.info(f"Starting currents extraction — preset: {preset.name}")
    _IMPL._check_cancel(cancel_event)
    timing, crs_type = _IMPL._load_model_definition(model_type, model_def_path, grid_path, logger)
    _IMPL._check_cancel(cancel_event)
    bnd_loc_geo, _native, _pairs = _IMPL._load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        sources = _build_boundary_sources(source_groups, "currents", model_type, timing, logger, cancel_event=cancel_event)
        original_coords = {
            row["boundaries"]: (float(row["lat"]), float(row["lon"]))
            for _, row in bnd_loc_geo.iterrows()
        }
        bnd_loc_geo, _unused_pairs, qc_summary, point_sources = _prepare_multi_source_sampling(
            model_type,
            "currents",
            sources,
            bnd_loc_geo,
            _pairs,
            logger,
            max_repair_distance_m=max_repair_distance_m,
            progress_callback=lambda p: progress_callback(0.1 * p),
            cancel_event=cancel_event,
        )
    else:
        ds = _IMPL._open_dataset(data_files, preset, timing, logger)
        _IMPL._check_cancel(cancel_event)
        qc_vars = [preset.current_u_var, preset.current_v_var]
        qc_context = _qc_context(ds, qc_vars, preset, logger, cancel_event)
        bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
            model_type,
            "currents",
            ds,
            bnd_loc_geo,
            _pairs,
            qc_vars,
            preset,
            qc_context,
            logger,
            progress_callback=lambda p: progress_callback(0.1 * p),
            cancel_event=cancel_event,
        )
    progress_callback(0.1)
    if not multi_source:
        native_minutes = _native_time_minutes_for_model(model_type, ds, preset, timing)
    output_minutes = _build_output_minutes(model_type, timing, timestep)
    output_time_list = output_minutes.tolist()
    total_points = len(bnd_loc_geo)
    logger.info(f"Extracting currents for {total_points} boundary points...")
    if multi_source:
        caches = _prefetch_source_caches(sources, bnd_loc_geo, point_sources, logger, cancel_event)
    else:
        corner_cache = _prefetch_corner_series(
            ds, _sampling_points(bnd_loc_geo), [preset.current_u_var, preset.current_v_var], preset, logger, cancel_event
        )
    u_dict: Dict[str, List[float]] = {}
    v_dict: Dict[str, List[float]] = {}
    for idx, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        _IMPL._check_cancel(cancel_event)
        if multi_source:
            src_idx = point_sources[row["boundaries"]]
            src = sources[src_idx]
            u_dict[row["boundaries"]] = _resample_scalar_variable(
                row, src["ds"], src["preset"].current_u_var, src["preset"],
                src["native_minutes"], output_minutes, logger=logger, cache=caches[src_idx],
            )
            v_dict[row["boundaries"]] = _resample_scalar_variable(
                row, src["ds"], src["preset"].current_v_var, src["preset"],
                src["native_minutes"], output_minutes, logger=logger, cache=caches[src_idx],
            )
        else:
            u_dict[row["boundaries"]] = _resample_scalar_variable(
                row, ds, preset.current_u_var, preset, native_minutes, output_minutes,
                logger=logger, cache=corner_cache,
            )
            v_dict[row["boundaries"]] = _resample_scalar_variable(
                row, ds, preset.current_v_var, preset, native_minutes, output_minutes,
                logger=logger, cache=corner_cache,
            )
        progress_callback(0.1 + 0.7 * (idx + 1) / max(1, total_points))
        _log_progress(logger, idx, total_points)
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        source_labels = [s["label"] for s in sources]
        _report_seam_consistency("currents (u)", u_dict, original_coords, point_sources, source_labels, logger)
        _report_seam_consistency("currents (v)", v_dict, original_coords, point_sources, source_labels, logger)
    logger.info("Currents extraction complete")

    if model_type == "delft3d4":
        boundary_data: Dict[str, Dict[str, List[float]]] = {}
        for key in u_dict:
            bnd_name = key[:-2]
            suffix = key[-1]
            boundary_data.setdefault(bnd_name, {})
            boundary_data[bnd_name][f"u_{suffix}"] = u_dict[key]
            boundary_data[bnd_name][f"v_{suffix}"] = v_dict[key]
        output_path = _default_output_path(output_path, output_dir, model_def_path, "_currents.bct")
        _IMPL.writer_delft3d4.write_bct_currents(
            output_path,
            boundary_data,
            timing["reference_date"],
            timing["tstart"],
            timing["tstop"],
            timestep,
        )
    else:
        boundary_data = {
            key: {"time": output_time_list, "values_x": u_dict[key], "values_y": v_dict[key]}
            for key in u_dict
        }
        output_path = _default_output_path(output_path, output_dir, model_def_path, "_currents.bc")
        # D-Flow FM full velocity vector: the .bc holds the x/y components as
        # columns named ux/uy (manual E.2.4.3), and the boundary is referenced
        # in the .ext with quantity=uxuyadvectionvelocitybnd (manual C.6.3).
        # Components are West->East / South->North, exactly as extracted; no
        # rotation to the boundary normal is required (manual section 22.6).
        _IMPL.writer_delft3dfm.write_bc(
            output_path,
            boundary_data,
            ("ux", "uy"),
            "m/s",
            str(timing["start_datetime"]),
            vector_name="uxuyadvectionvelocitybnd",
        )
        logger.info(
            "FM currents .bc written with a Vector=uxuyadvectionvelocitybnd:ux,uy "
            "declaration over ux/uy columns; reference it in the .ext [Boundary] "
            "block with quantity=uxuyadvectionvelocitybnd"
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")
    result = {
        "output_files": [output_path],
        "n_boundaries": len(boundary_data),
        "n_timesteps": len(output_time_list),
        "repaired_points": qc_summary["repaired_points"],
        "rejected_boundaries": qc_summary["rejected_boundaries"],
    }
    if multi_source:
        result["source_assignments"] = qc_summary["source_assignments"]
        result["intertidal_points"] = qc_summary.get("intertidal_points", [])
    return result


def extract_salinity(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: str | None = None,
    data_files: List[str] = None,
    preset_name: str = "Custom",
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = "N",
    output_path: str | None = None,
    output_dir: str = ".",
    source_groups: Optional[List[Dict[str, Any]]] = None,
    max_repair_distance_m: Optional[float] = None,
    logger=None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event=None,
) -> Dict[str, Any]:
    if logger is None:
        logger = _IMPL.logging.getLogger("repd3d.bc")
    multi_source = bool(source_groups)
    if multi_source:
        logger.info(f"Starting salinity extraction — {len(source_groups)} source group(s)")
    else:
        preset = _IMPL.preset_module.get_bc_preset(preset_name)
        if not preset.salinity_var:
            raise ValueError(
                f"Preset '{preset.name}' has no salinity_var. Use 'Custom' mode and select a salinity variable."
            )
        logger.info(f"Starting salinity extraction — preset: {preset.name}")
    _IMPL._check_cancel(cancel_event)
    timing, crs_type = _IMPL._load_model_definition(model_type, model_def_path, grid_path, logger)
    _IMPL._check_cancel(cancel_event)
    bnd_loc_geo, _native, _pairs = _IMPL._load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        sources = _build_boundary_sources(source_groups, "salinity", model_type, timing, logger, cancel_event=cancel_event)
        original_coords = {
            row["boundaries"]: (float(row["lat"]), float(row["lon"]))
            for _, row in bnd_loc_geo.iterrows()
        }
        bnd_loc_geo, _unused_pairs, qc_summary, point_sources = _prepare_multi_source_sampling(
            model_type,
            "salinity",
            sources,
            bnd_loc_geo,
            _pairs,
            logger,
            max_repair_distance_m=max_repair_distance_m,
            progress_callback=lambda p: progress_callback(0.1 * p),
            cancel_event=cancel_event,
        )
    else:
        ds = _IMPL._open_dataset(data_files, preset, timing, logger)
        _IMPL._check_cancel(cancel_event)
        qc_context = _qc_context(ds, [preset.salinity_var], preset, logger, cancel_event)
        bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
            model_type,
            "salinity",
            ds,
            bnd_loc_geo,
            _pairs,
            [preset.salinity_var],
            preset,
            qc_context,
            logger,
            progress_callback=lambda p: progress_callback(0.1 * p),
            cancel_event=cancel_event,
        )
    progress_callback(0.1)
    if not multi_source:
        native_minutes = _native_time_minutes_for_model(model_type, ds, preset, timing)
    output_minutes = _build_output_minutes(model_type, timing, timestep)
    output_time_list = output_minutes.tolist()
    total_points = len(bnd_loc_geo)
    raw: Dict[str, List[float]] = {}
    logger.info(f"Extracting salinity for {total_points} boundary points...")
    if multi_source:
        caches = _prefetch_source_caches(sources, bnd_loc_geo, point_sources, logger, cancel_event)
    else:
        corner_cache = _prefetch_corner_series(
            ds, _sampling_points(bnd_loc_geo), [preset.salinity_var], preset, logger, cancel_event
        )
    for idx, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        _IMPL._check_cancel(cancel_event)
        if multi_source:
            src_idx = point_sources[row["boundaries"]]
            src = sources[src_idx]
            raw[row["boundaries"]] = _resample_scalar_variable(
                row, src["ds"], src["preset"].salinity_var, src["preset"],
                src["native_minutes"], output_minutes, logger=logger, cache=caches[src_idx],
            )
        else:
            raw[row["boundaries"]] = _resample_scalar_variable(
                row, ds, preset.salinity_var, preset, native_minutes, output_minutes,
                logger=logger, cache=corner_cache,
            )
        progress_callback(0.1 + 0.7 * (idx + 1) / max(1, total_points))
        _log_progress(logger, idx, total_points)
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        _report_seam_consistency(
            "salinity", raw, original_coords, point_sources, [s["label"] for s in sources], logger,
        )
    logger.info("Salinity extraction complete")

    if model_type == "delft3d4":
        boundary_data: Dict[str, Dict[str, List[float]]] = {}
        for key, values in raw.items():
            bnd_name = key[:-2]
            suffix = key[-1]
            boundary_data.setdefault(bnd_name, {})
            boundary_data[bnd_name][f"sal_{suffix}"] = values
        output_path = _default_output_path(output_path, output_dir, model_def_path, "_salinity.bcc")
        _IMPL.writer_delft3d4.write_bcc_salinity(
            output_path,
            boundary_data,
            timing["reference_date"],
            timing["tstart"],
            timing["tstop"],
            timestep,
        )
    else:
        boundary_data = {
            key: {"time": output_time_list, "values": values}
            for key, values in raw.items()
        }
        output_path = _default_output_path(output_path, output_dir, model_def_path, "_salinity.bc")
        _IMPL.writer_delft3dfm.write_bc(
            output_path,
            boundary_data,
            "salinitybnd",
            "ppt",
            str(timing["start_datetime"]),
        )

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")
    result = {
        "output_files": [output_path],
        "n_boundaries": len(boundary_data),
        "n_timesteps": len(output_time_list),
        "repaired_points": qc_summary["repaired_points"],
        "rejected_boundaries": qc_summary["rejected_boundaries"],
    }
    if multi_source:
        result["source_assignments"] = qc_summary["source_assignments"]
        result["intertidal_points"] = qc_summary.get("intertidal_points", [])
    return result


def extract_sediments(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: str | None = None,
    data_files: List[str] = None,
    preset_name: str = "Copernicus BGC",
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = "N",
    output_path: str | None = None,
    output_dir: str = ".",
    source_groups: Optional[List[Dict[str, Any]]] = None,
    max_repair_distance_m: Optional[float] = None,
    logger=None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event=None,
) -> Dict[str, Any]:
    if logger is None:
        logger = _IMPL.logging.getLogger("repd3d.bc")
    multi_source = bool(source_groups)
    if multi_source:
        logger.info(f"Starting sediment extraction — {len(source_groups)} source group(s)")
    else:
        preset = _IMPL.preset_module.get_bc_preset(preset_name)
        if not preset.sediment_conc_vars:
            raise ValueError(
                f"Preset '{preset.name}' has no sediment_conc_vars. Use 'Copernicus BGC' or 'Custom' with sediment variables filled in."
            )
        n_frac = len(preset.sediment_conc_vars)
        fraction_names = list(preset.sediment_fraction_names)
        if len(fraction_names) != n_frac:
            fraction_names = [f"fraction_{i + 1}" for i in range(n_frac)]
        split_ratios = list(preset.sediment_split_ratios)
        if len(split_ratios) != n_frac:
            split_ratios = [1.0] * n_frac
        unique_vars = list(dict.fromkeys(preset.sediment_conc_vars))
        logger.info(f"Starting sediment extraction — preset: {preset.name} ({n_frac} fraction(s))")
    _IMPL._check_cancel(cancel_event)
    timing, crs_type = _IMPL._load_model_definition(model_type, model_def_path, grid_path, logger)
    _IMPL._check_cancel(cancel_event)
    bnd_loc_geo, _native, _pairs = _IMPL._load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        sources = _build_boundary_sources(source_groups, "sediments", model_type, timing, logger, cancel_event=cancel_event)
        # Fraction metadata comes from the first built source so per-group
        # Custom overrides are honoured; _build_boundary_sources already
        # enforced an equal fraction count across groups.
        meta_preset = sources[0]["preset"]
        n_frac = len(meta_preset.sediment_conc_vars)
        fraction_names = list(meta_preset.sediment_fraction_names)
        if len(fraction_names) != n_frac:
            fraction_names = [f"fraction_{i + 1}" for i in range(n_frac)]
        split_ratios = list(meta_preset.sediment_split_ratios)
        if len(split_ratios) != n_frac:
            split_ratios = [1.0] * n_frac
        logger.info(f"Sediment fractions: {n_frac} — {', '.join(fraction_names)}")
        original_coords = {
            row["boundaries"]: (float(row["lat"]), float(row["lon"]))
            for _, row in bnd_loc_geo.iterrows()
        }
        bnd_loc_geo, _unused_pairs, qc_summary, point_sources = _prepare_multi_source_sampling(
            model_type,
            "sediments",
            sources,
            bnd_loc_geo,
            _pairs,
            logger,
            max_repair_distance_m=max_repair_distance_m,
            progress_callback=lambda p: progress_callback(0.05 * p),
            cancel_event=cancel_event,
        )
    else:
        ds = _IMPL._open_dataset(data_files, preset, timing, logger)
        _IMPL._check_cancel(cancel_event)
        qc_context = _qc_context(ds, unique_vars, preset, logger, cancel_event)
        bnd_loc_geo, _unused_pairs, qc_summary = _prepare_boundary_sampling(
            model_type,
            "sediments",
            ds,
            bnd_loc_geo,
            _pairs,
            unique_vars,
            preset,
            qc_context,
            logger,
            progress_callback=lambda p: progress_callback(0.05 * p),
            cancel_event=cancel_event,
        )
    progress_callback(0.05)
    if not multi_source:
        native_minutes = _native_time_minutes_for_model(model_type, ds, preset, timing)
    output_minutes = _build_output_minutes(model_type, timing, timestep)
    output_time_list = output_minutes.tolist()
    raw: Dict[Any, Dict[str, List[float]]] = {}
    if multi_source:
        caches = _prefetch_source_caches(sources, bnd_loc_geo, point_sources, logger, cancel_event)
        # Fraction index keys: each source maps fraction f_idx to its own variable name.
        total_steps = max(1, n_frac * max(1, len(bnd_loc_geo)))
        step = 0
        for f_idx in range(n_frac):
            raw[f_idx] = {}
            for _, row in bnd_loc_geo.iterrows():
                _IMPL._check_cancel(cancel_event)
                src_idx = point_sources[row["boundaries"]]
                src = sources[src_idx]
                var_name = src["preset"].sediment_conc_vars[f_idx]
                raw[f_idx][row["boundaries"]] = _resample_scalar_variable(
                    row, src["ds"], var_name, src["preset"], src["native_minutes"], output_minutes,
                    logger=logger, cache=caches[src_idx],
                )
                step += 1
                progress_callback(0.05 + 0.75 * step / total_steps)
        _report_seam_consistency(
            "sediments (fraction 1)", raw[0], original_coords, point_sources,
            [s["label"] for s in sources], logger,
        )
    else:
        corner_cache = _prefetch_corner_series(
            ds, _sampling_points(bnd_loc_geo), unique_vars, preset, logger, cancel_event
        )
        total_steps = max(1, len(unique_vars) * max(1, len(bnd_loc_geo)))
        step = 0
        for var_name in unique_vars:
            raw[var_name] = {}
            for _, row in bnd_loc_geo.iterrows():
                _IMPL._check_cancel(cancel_event)
                raw[var_name][row["boundaries"]] = _resample_scalar_variable(
                    row, ds, var_name, preset, native_minutes, output_minutes, logger=logger,
                    cache=corner_cache,
                )
                step += 1
                progress_callback(0.05 + 0.75 * step / total_steps)
    _IMPL._check_cancel(cancel_event)
    logger.info("Sediment extraction complete")

    def _fraction_series(f_idx: int) -> Dict[str, List[float]]:
        if multi_source:
            return raw[f_idx]
        return raw[preset.sediment_conc_vars[f_idx]]

    if model_type == "delft3d4":
        boundary_data: Dict[str, Dict[str, Dict[str, List[float]]]] = {}
        for f_idx, fraction_name in enumerate(fraction_names):
            ratio = float(split_ratios[f_idx])
            boundary_data[fraction_name] = {}
            for key, values in _fraction_series(f_idx).items():
                bnd_name = key[:-2]
                suffix = key[-1]
                boundary_data[fraction_name].setdefault(bnd_name, {})
                boundary_data[fraction_name][bnd_name][f"conc_{suffix}"] = (
                    _IMPL.np.asarray(values, dtype=float) * ratio
                ).tolist()
        output_path = _default_output_path(output_path, output_dir, model_def_path, ".bcc")
        _IMPL.writer_delft3d4.write_bcc(
            output_path,
            boundary_data,
            fraction_names,
            timing["reference_date"],
            timing["tstart"],
            timing["tstop"],
            timestep,
        )
        n_boundaries = len(bnd_loc_geo) // 2
    else:
        # D-Flow FM sediment BC: each fraction is its own quantity
        # sedfracbnd<fractionname> (manual Table C.9), and every forcing block's
        # Name must equal the polyline support point so FM can bind it (manual
        # C.5.1.2). Multiple fractions at one point become multiple [forcing]
        # blocks that share the same Name but differ in Quantity. The dict key
        # only needs to be unique; the emitted Name comes from the 'name' field.
        boundary_data = {}
        for f_idx, fraction_name in enumerate(fraction_names):
            ratio = float(split_ratios[f_idx])
            quantity = f"sedfracbnd{fraction_name}"
            for key, values in _fraction_series(f_idx).items():
                boundary_data[f"{key}__{fraction_name}"] = {
                    "name": key,
                    "quantity": quantity,
                    "unit": "kg/m3",
                    "time": output_time_list,
                    "values": (_IMPL.np.asarray(values, dtype=float) * ratio).tolist(),
                }
        output_path = _default_output_path(output_path, output_dir, model_def_path, "_sediments.bc")
        _IMPL.writer_delft3dfm.write_bc(
            output_path,
            boundary_data,
            "sedfracbnd",
            "kg/m3",
            str(timing["start_datetime"]),
        )
        n_boundaries = len(bnd_loc_geo)

    progress_callback(1.0)
    logger.info(f"Output written: {output_path}")
    result = {
        "output_files": [output_path],
        "n_boundaries": n_boundaries,
        "n_fractions": n_frac,
        "n_timesteps": len(output_time_list),
        "repaired_points": qc_summary["repaired_points"],
        "rejected_boundaries": qc_summary["rejected_boundaries"],
    }
    if multi_source:
        result["source_assignments"] = qc_summary["source_assignments"]
        result["intertidal_points"] = qc_summary.get("intertidal_points", [])
    return result


def extract_waves(model_type: str, model_def_path: str, boundary_path: str, grid_path: str | None = None, mdw_template_path: str | None = None, data_files: List[str] | None = None, preset_name: str = "EasyGSH", timestep: float = 20.0, utm_zone: int = 32, hemisphere: str = "N", output_path: str | None = None, output_dir: str = ".", source_groups: Optional[List[Dict[str, Any]]] = None, max_repair_distance_m: Optional[float] = None, logger=None, progress_callback: Callable[[float], None] = lambda x: None, cancel_event=None, wave_timing_source: str = "Auto from .mdf", wave_reference_date: str | None = None, wave_start_datetime: str | None = None, wave_end_datetime: str | None = None) -> Dict[str, Any]:
    if logger is None:
        logger = _IMPL.logging.getLogger("repd3d.bc")
    if model_type != "delft3d4":
        raise NotImplementedError(
            "Delft3D-FM has no wave boundary-condition quantity in its external "
            "forcings (D-Flow FM User Manual Table C.9). Wave forcing for an "
            "FM-coupled model is prescribed through the Delft3D-WAVE module "
            "(.mdw + .bcw), not a .bc file. Re-run wave extraction with model "
            "type 'delft3d4' (optionally supplying an .mdw template) to produce "
            "the .bcw/.mdw wave boundary, which the FM-coupled D-Waves run uses."
        )
    multi_source = bool(source_groups)
    if multi_source:
        logger.info(f"Starting wave extraction - {len(source_groups)} source group(s)")
    else:
        preset = _IMPL.preset_module.get_bc_preset(preset_name)
        logger.info(f"Starting wave extraction - preset: {preset.name}")
    _IMPL._check_cancel(cancel_event)
    timing, crs_type = _resolve_wave_timing_and_crs(model_type, model_def_path, grid_path, logger, wave_timing_source, wave_reference_date, wave_start_datetime, wave_end_datetime)
    _IMPL._check_cancel(cancel_event)
    bnd_loc_geo, _bnd_native, bnd_pairs = _IMPL._load_boundary_locations(model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger)
    _IMPL._check_cancel(cancel_event)
    if multi_source:
        sources = _build_boundary_sources(source_groups, "waves", model_type, timing, logger, cancel_event=cancel_event)
        bnd_loc_geo, bnd_pairs, qc_summary, point_sources = _prepare_multi_source_sampling(
            model_type, "waves", sources, bnd_loc_geo, bnd_pairs, logger,
            waves_use_a_only=True, max_repair_distance_m=max_repair_distance_m,
            progress_callback=lambda frac: progress_callback(0.1 * frac), cancel_event=cancel_event,
        )
    else:
        ds = _IMPL._open_dataset(data_files, preset, timing, logger)
        _IMPL._check_cancel(cancel_event)
        wave_qc_vars = [preset.wave_hs_var, preset.wave_tp_var]
        if preset.wave_dir_format == "vector":
            wave_qc_vars.extend([preset.wave_dir_x_var, preset.wave_dir_y_var])
        else:
            wave_qc_vars.append(preset.wave_dir_var)
        if preset.wave_dir_spread_var and not preset.wave_dir_spread_use_default:
            wave_qc_vars.append(preset.wave_dir_spread_var)
        qc_context = _qc_context(ds, wave_qc_vars, preset, logger, cancel_event)
        qc_context["derived_vars"] = _wave_derived_variables(preset, wave_qc_vars)
        _log_derived_policy(logger, qc_context["derived_vars"])
        bnd_loc_geo, bnd_pairs, qc_summary = _prepare_boundary_sampling(model_type, "waves", ds, bnd_loc_geo, bnd_pairs, wave_qc_vars, preset, qc_context, logger, waves_use_a_only=True, progress_callback=lambda frac: progress_callback(0.1 * frac), cancel_event=cancel_event)
    progress_callback(0.1)
    if not multi_source:
        native_minutes = _native_time_minutes(ds, preset, timing)
    output_minutes = _build_wave_output_minutes(timing, timestep)
    bnd_loc_a = bnd_loc_geo[bnd_loc_geo["boundaries"].str.endswith("_a")].copy()
    if bnd_loc_a.empty:
        bnd_loc_a = bnd_loc_geo
    total_points = len(bnd_loc_a)
    logger.info(f"Extracting wave data for {total_points} boundary points...")
    if multi_source:
        caches = _prefetch_source_caches(sources, bnd_loc_a, point_sources, logger, cancel_event)
    else:
        corner_cache = _prefetch_corner_series(ds, _sampling_points(bnd_loc_a), wave_qc_vars, preset, logger, cancel_event)
    wave_data: Dict[str, Dict[str, List[float]]] = {}
    for idx, (_, row) in enumerate(bnd_loc_a.iterrows()):
        _IMPL._check_cancel(cancel_event)
        if multi_source:
            src_idx = point_sources[row["boundaries"]]
            src = sources[src_idx]
            point_preset = src["preset"]
            point_ds = src["ds"]
            point_native_minutes = src["native_minutes"]
            point_cache = caches[src_idx]
        else:
            point_preset = preset
            point_ds = ds
            point_native_minutes = native_minutes
            point_cache = corner_cache
        hs_stats: Dict[str, int] = {}
        tp_stats: Dict[str, int] = {}
        hs_native = _IMPL.interp.bilinear_interpolation(row["lat"], row["lon"], point_ds, variable=point_preset.wave_hs_var, lat_coord=point_preset.lat_coord, lon_coord=point_preset.lon_coord, stats=hs_stats, cache=point_cache)
        tp_native = _IMPL.interp.bilinear_interpolation(row["lat"], row["lon"], point_ds, variable=point_preset.wave_tp_var, lat_coord=point_preset.lat_coord, lon_coord=point_preset.lon_coord, stats=tp_stats, cache=point_cache)
        _check_plausible_range(str(row["boundaries"]), point_preset.wave_hs_var, hs_native, point_preset)
        _check_plausible_range(str(row["boundaries"]), point_preset.wave_tp_var, tp_native, point_preset)
        _log_fill_stats(logger, str(row["boundaries"]), point_preset.wave_hs_var, hs_stats)
        _log_fill_stats(logger, str(row["boundaries"]), point_preset.wave_tp_var, tp_stats)
        point_data = {
            "hs": _resample_linear_series(point_native_minutes, hs_native, output_minutes).tolist(),
            "tp": _resample_linear_series(point_native_minutes, tp_native, output_minutes).tolist(),
        }
        if point_preset.wave_dir_format == "vector":
            dir_x = _IMPL.interp.bilinear_interpolation(row["lat"], row["lon"], point_ds, variable=point_preset.wave_dir_x_var, lat_coord=point_preset.lat_coord, lon_coord=point_preset.lon_coord, cache=point_cache)
            dir_y = _IMPL.interp.bilinear_interpolation(row["lat"], row["lon"], point_ds, variable=point_preset.wave_dir_y_var, lat_coord=point_preset.lat_coord, lon_coord=point_preset.lon_coord, cache=point_cache)
            dir_x_resampled = _resample_linear_series(point_native_minutes, dir_x, output_minutes)
            dir_y_resampled = _resample_linear_series(point_native_minutes, dir_y, output_minutes)
            point_data["dir"] = _IMPL._wave_direction_from_vectors(dir_x_resampled, dir_y_resampled).tolist()
        else:
            dir_native = _IMPL.interp.bilinear_interpolation(row["lat"], row["lon"], point_ds, variable=point_preset.wave_dir_var, lat_coord=point_preset.lat_coord, lon_coord=point_preset.lon_coord, cache=point_cache)
            point_data["dir"] = _resample_circular_degrees(point_native_minutes, dir_native, output_minutes).tolist()
        if point_preset.wave_dir_spread_var and not point_preset.wave_dir_spread_use_default:
            spread_stats: Dict[str, int] = {}
            spread_native = _IMPL.interp.bilinear_interpolation(row["lat"], row["lon"], point_ds, variable=point_preset.wave_dir_spread_var, lat_coord=point_preset.lat_coord, lon_coord=point_preset.lon_coord, stats=spread_stats, cache=point_cache)
            _log_fill_stats(logger, str(row["boundaries"]), point_preset.wave_dir_spread_var, spread_stats)
            point_data["spread"] = _resample_linear_series(point_native_minutes, spread_native, output_minutes).tolist()
        else:
            point_data["spread"] = [point_preset.wave_dir_spread_default] * len(output_minutes)
            logger.warning(f"No directional spreading variable - using default {point_preset.wave_dir_spread_default} deg")
        point_data["time"] = output_minutes.tolist()
        clean_name = row["boundaries"][:-2] if row["boundaries"].endswith("_a") else row["boundaries"]
        wave_data[clean_name] = point_data
        progress_callback(0.1 + 0.7 * (idx + 1) / total_points)
        if (idx + 1) % max(1, total_points // 5) == 0:
            logger.info(f"  Processed {idx + 1}/{total_points} points")
    _IMPL._check_cancel(cancel_event)
    logger.info("Wave extraction complete")
    output_files: List[str] = []
    if model_type == "delft3d4":
        if output_path is None:
            output_path = str(Path(output_dir) / "wave.bcw")
        _IMPL.writer_delft3d4.write_bcw(output_path, wave_data, timing["reference_date"])
        output_files.append(output_path)
        if mdw_template_path:
            mdw_out = str(Path(output_dir) / f"{Path(mdw_template_path).stem}_updated.mdw")
            _IMPL.writer_mdw.write_mdw(mdw_template_path, bnd_pairs, mdw_out)
            output_files.append(mdw_out)
            logger.info(f"Updated .mdw written: {mdw_out}")
        else:
            logger.warning("No .mdw template provided - skipping .mdw update. The .bcw boundary names will not be linked to a wave model definition.")
    else:
        # Unreachable: the early guard rejects non-delft3d4 model types. Kept as
        # a defensive stop so no invalid FM wave .bc is ever written.
        raise NotImplementedError(
            "Delft3D-FM wave boundary conditions are not written as a .bc file; "
            "use model type 'delft3d4' to produce the .bcw/.mdw wave boundary."
        )
    progress_callback(1.0)
    logger.info(f"Output written: {', '.join(output_files)}")
    n_timesteps = len(next(iter(wave_data.values()))["hs"]) if wave_data else 0
    result = {"output_files": output_files, "n_boundaries": len(wave_data), "n_timesteps": n_timesteps, "repaired_points": qc_summary["repaired_points"], "rejected_boundaries": qc_summary["rejected_boundaries"]}
    if multi_source:
        result["source_assignments"] = qc_summary["source_assignments"]
        result["intertidal_points"] = qc_summary.get("intertidal_points", [])
    return result


_IMPL._fast_assess_variable_quality = _fast_assess_variable_quality
_IMPL._fast_evaluate_point_quality = _fast_evaluate_point_quality
_IMPL._evaluate_point_status = _evaluate_point_status
_IMPL._prepare_boundary_sampling = _prepare_boundary_sampling
_IMPL.extract_water_level = extract_water_level
_IMPL.extract_currents = extract_currents
_IMPL.extract_salinity = extract_salinity
_IMPL.extract_sediments = extract_sediments
_IMPL.extract_waves = extract_waves

for _name in dir(_IMPL):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_IMPL, _name)

globals()["_fast_assess_variable_quality"] = _fast_assess_variable_quality
globals()["_fast_evaluate_point_quality"] = _fast_evaluate_point_quality
globals()["_evaluate_point_status"] = _evaluate_point_status
globals()["_prepare_boundary_sampling"] = _prepare_boundary_sampling
globals()["extract_water_level"] = extract_water_level
globals()["extract_currents"] = extract_currents
globals()["extract_salinity"] = extract_salinity
globals()["extract_sediments"] = extract_sediments
globals()["extract_waves"] = extract_waves
globals()["set_intertidal_tolerance"] = set_intertidal_tolerance

__all__ = [name for name in globals() if not name.startswith("__")]


def __getattr__(name: str) -> Any:
    return getattr(_IMPL, name)


if __name__ == "__main__":
    import logging
    import tempfile
    import pandas as pd
    import xarray as xr
    from presets import TRILAWATT

    print("=== backend_bc_extractor.py self-test ===\n")
    all_pass = True

    def check(name: str, result: Any, expected: Any) -> None:
        nonlocal_all = globals()
        ok = result == expected
        print(f"  {name}: [{'PASS' if ok else 'FAIL'}]")
        if not ok:
            nonlocal_all["all_pass"] = False
            print(f"    got {result!r}, expected {expected!r}")

    def make_dataset(corner_modes: Dict[str, str]):
        time = pd.date_range("2020-01-01", periods=4, freq="h")
        lat = [0.0, 1.0]
        lon = [0.0, 1.0]
        values = _IMPL.np.ones((len(time), len(lat), len(lon)), dtype=float)
        corner_lookup = {"Q11": (0, 0), "Q12": (0, 1), "Q21": (1, 0), "Q22": (1, 1)}
        for key, mode in corner_modes.items():
            i, j = corner_lookup[key]
            if mode == "all_nan":
                values[:, i, j] = _IMPL.np.nan
            elif mode == "partial":
                values[1:-1, i, j] = _IMPL.np.nan
        return xr.Dataset({"sea_surface_height_2d": (("time", "lat", "lon"), values)}, coords={"time": time, "lat": lat, "lon": lon})

    qc_full = _IMPL._build_qc_context(make_dataset({}), ["sea_surface_height_2d"], TRILAWATT)
    fast_full = _fast_evaluate_point_quality(0.2, 0.2, make_dataset({}), ["sea_surface_height_2d"], TRILAWATT, qc_full)
    check("4 fully valid corners fast-pass", fast_full["decision"], "accept")

    # Point (0.2, 0.2): Q11 64 %, Q12 16 %, Q21 16 %, Q22 4 % weight.
    ds_three = make_dataset({"Q12": "all_nan"})
    qc_three = _IMPL._build_qc_context(ds_three, ["sea_surface_height_2d"], TRILAWATT)
    fast_three = _fast_evaluate_point_quality(0.2, 0.2, ds_three, ["sea_surface_height_2d"], TRILAWATT, qc_three)
    check("3 fully valid corners: 16 % all-NaN corner hard-fails fast QC", fast_three["decision"], "hard_fail")

    # The former 3-of-4 override: three corners fully valid, the fourth wet for
    # one timestep only (a storm). It carried 16 % weight here, 51 % in the
    # 2026-10-13 incident. It must now go through the strict assessment and be
    # rejected (valid fraction 25 % < 75 %, weight 16 % > 5 %).
    set_intertidal_tolerance(False)
    ds_storm = make_dataset({"Q12": "partial"})
    ds_storm["sea_surface_height_2d"].values[:, 0, 1] = _IMPL.np.nan
    ds_storm["sea_surface_height_2d"].values[2, 0, 1] = 2.0
    qc_storm = _IMPL._build_qc_context(ds_storm, ["sea_surface_height_2d"], TRILAWATT)
    fast_storm = _fast_evaluate_point_quality(0.2, 0.2, ds_storm, ["sea_surface_height_2d"], TRILAWATT, qc_storm)
    check("storm-only corner no longer fast-accepted (no 3-of-4 override)", fast_storm["decision"], "ambiguous")
    storm_status = _evaluate_point_status(0.2, 0.2, ds_storm, ["sea_surface_height_2d"], TRILAWATT, qc_storm)
    check("storm-only corner rejected by strict QC", storm_status["reliable"], False)
    check("storm-only corner verdict came from strict path", storm_status["path"], "strict")
    # The same storm corner at 0.04 % weight (point 0.98, 0.02) is harmless: strict QC passes.
    light_status = _evaluate_point_status(0.98, 0.02, ds_storm, ["sea_surface_height_2d"], TRILAWATT, qc_storm)
    check("storm-only corner at negligible weight still passes strict QC",
          (light_status["reliable"], light_status["path"]), (True, "strict"))

    # Plausibility guard: impossible sampled values abort before anything is written.
    try:
        _check_plausible_range("BND_a", TRILAWATT.water_level_var, [1.0, -83.4, 0.2], TRILAWATT)
        check("plausibility guard raises on -83 m water level", "no error", "ValueError")
    except ValueError:
        check("plausibility guard raises on -83 m water level", "ValueError", "ValueError")
    try:
        _check_plausible_range("BND_a", TRILAWATT.water_level_var, [1.0, -1.8, 4.5], TRILAWATT)
        check("plausibility guard accepts realistic water levels", "no error", "no error")
    except ValueError:
        check("plausibility guard accepts realistic water levels", "ValueError", "no error")
    check("plausibility guard skips unlisted variables", _plausible_range_for(TRILAWATT, "some_other_var"), None)

    ds_hard_fail = make_dataset({"Q21": "all_nan", "Q22": "all_nan"})
    qc_hard_fail = _IMPL._build_qc_context(ds_hard_fail, ["sea_surface_height_2d"], TRILAWATT)
    fast_hard_fail = _fast_evaluate_point_quality(
        0.25, 0.75, ds_hard_fail, ["sea_surface_height_2d"], TRILAWATT, qc_hard_fail
    )
    check("fast hard-fail detected", fast_hard_fail["decision"], "hard_fail")
    hard_fail_summary = _format_fast_assessment_summary(
        "sea_surface_height_2d",
        fast_hard_fail["assessments"]["sea_surface_height_2d"],
    )
    check("fast hard-fail summary includes bad weight", "bad weight" in hard_fail_summary, True)
    check("fast hard-fail summary includes all-NaN", "all-NaN" in hard_fail_summary, True)

    ds_tidal = make_dataset({"Q11": "partial", "Q21": "partial"})
    set_intertidal_tolerance(False)
    qc_tidal_off = _IMPL._build_qc_context(ds_tidal, ["sea_surface_height_2d"], TRILAWATT)
    check("intertidal off fails", _evaluate_point_status(0.5, 0.02, ds_tidal, ["sea_surface_height_2d"], TRILAWATT, qc_tidal_off)["reliable"], False)
    set_intertidal_tolerance(True)
    qc_tidal_on = _IMPL._build_qc_context(ds_tidal, ["sea_surface_height_2d"], TRILAWATT)
    check("intertidal on passes", _evaluate_point_status(0.5, 0.02, ds_tidal, ["sea_surface_height_2d"], TRILAWATT, qc_tidal_on)["path"], "intertidal")
    # allow_intertidal=False (Neumann probes) must not reuse the tolerant verdict cached above.
    strict_only = _evaluate_point_status(0.5, 0.02, ds_tidal, ["sea_surface_height_2d"], TRILAWATT, qc_tidal_on, allow_intertidal=False)
    check("allow_intertidal=False bypasses the cached intertidal verdict", strict_only["reliable"], False)

    print("\nDerived wave quantities (gaps allowed, a value somewhere required):")
    time_d = pd.date_range("2020-01-01", periods=6, freq="h")
    hs_d = _IMPL.np.ones((6, 2, 2))
    tp_d = _IMPL.np.ones((6, 2, 2))
    tp_d[1:5, :, :] = _IMPL.np.nan                      # Tp undefined 4 of 6 steps in every corner (calm weather)
    coords_d = {"time": time_d, "lat": [0.0, 1.0], "lon": [0.0, 1.0]}
    ds_d = xr.Dataset({"hs": (("time", "lat", "lon"), hs_d), "tp": (("time", "lat", "lon"), tp_d)}, coords=coords_d)
    ctx_plain = _IMPL._build_qc_context(ds_d, ["hs", "tp"], TRILAWATT)
    check("undeclared: Tp gaps fail strict QC", _evaluate_point_status(0.5, 0.5, ds_d, ["hs", "tp"], TRILAWATT, ctx_plain, allow_intertidal=False)["reliable"], False)
    ctx_d = _IMPL._build_qc_context(ds_d, ["hs", "tp"], TRILAWATT)
    ctx_d["derived_vars"] = ("tp",)
    status_d = _evaluate_point_status(0.5, 0.5, ds_d, ["hs", "tp"], TRILAWATT, ctx_d, allow_intertidal=False)
    check("derived: Tp gaps accepted on the fast path", (status_d["reliable"], status_d["path"]), (True, "fast"))
    hs_p = hs_d.copy()
    hs_p[0, 0, 0] = _IMPL.np.nan                        # one Hs gap -> ambiguous fast verdict -> strict path
    ds_p = ds_d.assign(hs=(("time", "lat", "lon"), hs_p))
    ctx_p = _IMPL._build_qc_context(ds_p, ["hs", "tp"], TRILAWATT)
    ctx_p["derived_vars"] = ("tp",)
    status_p = _evaluate_point_status(0.5, 0.5, ds_p, ["hs", "tp"], TRILAWATT, ctx_p, allow_intertidal=False)
    check("derived: Tp gaps accepted on the strict path", (status_p["reliable"], status_p["path"]), (True, "strict"))
    tp_none = _IMPL.np.full((6, 2, 2), _IMPL.np.nan)
    ds_n = ds_d.assign(tp=(("time", "lat", "lon"), tp_none))
    for tolerance_flag in (False, True):
        ctx_n = _IMPL._build_qc_context(ds_n, ["hs", "tp"], TRILAWATT)
        ctx_n["derived_vars"] = ("tp",)
        status_n = _evaluate_point_status(0.5, 0.5, ds_n, ["hs", "tp"], TRILAWATT, ctx_n, allow_intertidal=tolerance_flag)
        check(f"derived: Tp without any value rejected (intertidal={tolerance_flag})",
              (status_n["reliable"], any("no source value" in issue for issue in status_n["issues"])), (False, True))
    tp_corner = tp_none.copy()
    tp_corner[2, 1, 1] = 5.0                            # a single value in one corner is enough
    ds_c = ds_d.assign(tp=(("time", "lat", "lon"), tp_corner))
    ctx_c = _IMPL._build_qc_context(ds_c, ["hs", "tp"], TRILAWATT)
    ctx_c["derived_vars"] = ("tp",)
    check("derived: one value in one corner accepted", _evaluate_point_status(0.5, 0.5, ds_c, ["hs", "tp"], TRILAWATT, ctx_c, allow_intertidal=False)["reliable"], True)
    series_c = _IMPL.interp.bilinear_interpolation(0.5, 0.5, ds_c, variable="tp", lat_coord="lat", lon_coord="lon")
    check("derived: extracted Tp has no NaN (held/bridged)", (bool(_IMPL.np.isnan(series_c).any()), float(series_c.min()), float(series_c.max())), (False, 5.0, 5.0))
    check("wave derived variables of the TrilaWatt preset",
          _wave_derived_variables(TRILAWATT, [TRILAWATT.wave_hs_var, TRILAWATT.wave_tp_var, TRILAWATT.wave_dir_x_var, TRILAWATT.wave_dir_y_var, TRILAWATT.wave_dir_spread_var]),
          (TRILAWATT.wave_tp_var, TRILAWATT.wave_dir_spread_var))

    manual_timing = _build_manual_wave_timing("20180101", "2018-05-01 00:00", "2018-05-01 01:00")
    check("manual wave tstart", manual_timing["tstart"], 172800.0)
    check("manual wave first time value", _build_wave_time_values(manual_timing, 3, 20.0)[0], 172800.0)
    check("manual wave output minutes inclusive", _build_wave_output_minutes(manual_timing, 20.0).tolist(), [172800.0, 172820.0, 172840.0, 172860.0])

    native_minutes = _IMPL.np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0], dtype=float)
    requested_hourly = _IMPL.np.array([0.0, 60.0, 120.0], dtype=float)
    requested_45 = _IMPL.np.array([0.0, 45.0, 90.0], dtype=float)
    base_series = _IMPL.np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], dtype=float)
    check(
        "wave subsample exact multiple",
        _resample_linear_series(native_minutes, base_series, requested_hourly).tolist(),
        [1.0, 4.0, 7.0],
    )
    check(
        "wave linear interpolate non-multiple",
        _resample_linear_series(native_minutes, base_series, requested_45).tolist(),
        [1.0, 3.25, 5.5],
    )
    check(
        "wave same timestep unchanged",
        _resample_linear_series(native_minutes, base_series, native_minutes).tolist(),
        base_series.tolist(),
    )
    wrapped_deg = _resample_circular_degrees(
        _IMPL.np.array([0.0, 60.0], dtype=float),
        _IMPL.np.array([350.0, 10.0], dtype=float),
        _IMPL.np.array([30.0], dtype=float),
    )
    check("wave circular direction avoids 180 crossing", round(float(wrapped_deg[0]) % 360.0, 6), 0.0)
    dir_x_native = _IMPL.np.array([1.0, 0.0, -1.0], dtype=float)
    dir_y_native = _IMPL.np.array([0.0, 1.0, 0.0], dtype=float)
    dir_x_hourly = _resample_linear_series(
        _IMPL.np.array([0.0, 30.0, 60.0], dtype=float),
        dir_x_native,
        _IMPL.np.array([0.0, 60.0], dtype=float),
    )
    dir_y_hourly = _resample_linear_series(
        _IMPL.np.array([0.0, 30.0, 60.0], dtype=float),
        dir_y_native,
        _IMPL.np.array([0.0, 60.0], dtype=float),
    )
    check(
        "wave vector direction sampled after component resample",
        _IMPL._wave_direction_from_vectors(dir_x_hourly, dir_y_hourly).tolist(),
        _IMPL._wave_direction_from_vectors(_IMPL.np.array([1.0, -1.0]), _IMPL.np.array([0.0, 0.0])).tolist(),
    )

    # EasyGSH -12 h: _open_dataset hands over source stamps 12 h before the model
    # window; the native axis must land back on the model window so the sample
    # taken at source time t - 12 h is written at model time t (v1 behaviour).
    from presets import EASYGSH
    offset_timing = {"reference_date_iso": "2003-01-01", "tstart": 217440.0, "tstop": 217500.0,
                     "start_datetime": pd.Timestamp("2003-06-01 00:00"), "end_datetime": pd.Timestamp("2003-06-01 01:00")}
    offset_ds = xr.Dataset(coords={EASYGSH.time_coord: pd.date_range("2003-05-31 12:00", periods=4, freq="20min")})
    offset_out = _build_output_minutes("delft3d4", offset_timing, 20.0)
    check("EasyGSH offset: native axis on the model clock (D3D4)",
          _native_time_minutes_for_model("delft3d4", offset_ds, EASYGSH, offset_timing).tolist(), offset_out.tolist())
    check("EasyGSH offset: native axis on the model clock (FM)",
          _native_time_minutes_for_model("fm", offset_ds, EASYGSH, offset_timing).tolist(), [0.0, 20.0, 40.0, 60.0])
    check("EasyGSH offset: wave native axis on the model clock",
          _native_time_minutes(offset_ds, EASYGSH, offset_timing).tolist(), offset_out.tolist())
    check("EasyGSH offset: series written positionally, last hour not clamped",
          _resample_linear_series(_native_time_minutes(offset_ds, EASYGSH, offset_timing),
                                  _IMPL.np.array([1.0, 2.0, 3.0, 4.0]), offset_out).tolist(),
          [1.0, 2.0, 3.0, 4.0])

    stream = io.StringIO()
    logger = logging.getLogger("backend_bc_test")
    logger.handlers[:] = []
    handler = logging.StreamHandler(stream)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    bnd_loc = pd.DataFrame([{"boundaries": "South_1_a", "lat": 0.2, "lon": 0.2}, {"boundaries": "South_1_b", "lat": 0.8, "lon": 0.0}])
    bnd_pairs = [{"bnd_name": "South_1", "x_a": 0.0, "y_a": 0.0, "x_b": 1.0, "y_b": 1.0}]
    _prepare_boundary_sampling("delft3d4", "water level", ds_three, bnd_loc, bnd_pairs, ["sea_surface_height_2d"], TRILAWATT, qc_three, logger)
    check("partially covered end is repaired toward its partner, not override-accepted",
          "South_1_a moved" in stream.getvalue() and "override" not in stream.getvalue(), True)

    # extract_water_level re-reads the .bnd to find N-typed sections, so the
    # probe needs a real file even though _load_boundary_locations is stubbed.
    _probe_bnd_dir = tempfile.TemporaryDirectory()
    _probe_bnd_path = str(Path(_probe_bnd_dir.name) / "bnd.bnd")
    Path(_probe_bnd_path).write_text(
        "BND                  Z T     1     1     1    10  1.0000000e+002\n"
    )

    def run_non_wave_timestep_probe(extractor_name: str, model_type: str, timestep_value: float) -> Dict[str, Any]:
        class DummyPreset:
            name = "Dummy"
            lat_coord = "lat"
            lon_coord = "lon"
            time_coord = "time"
            water_level_var = "wl"
            current_u_var = "u"
            current_v_var = "v"
            salinity_var = "sal"
            sediment_conc_vars = ["sed"]
            sediment_fraction_names = ["mud"]
            sediment_split_ratios = [1.0]

        native_time = pd.date_range("2018-01-01 01:00:00", periods=7, freq="20min")
        probe_ds = xr.Dataset(coords={"time": native_time, "lat": [0.0, 1.0], "lon": [0.0, 1.0]})
        timing = {
            "reference_date": "20180101",
            "reference_date_iso": "2018-01-01",
            "tstart": 60.0,
            "tstop": 180.0,
            "start_datetime": pd.Timestamp("2018-01-01 01:00:00"),
            "end_datetime": pd.Timestamp("2018-01-01 03:00:00"),
        }
        boundary_df = pd.DataFrame(
            [
                {"boundaries": "BND_a", "lat": 0.0, "lon": 0.0},
                {"boundaries": "BND_b", "lat": 1.0, "lon": 1.0},
            ]
        )
        pair_data = [{"bnd_name": "BND"}]
        # Values sit inside PLAUSIBLE_RANGES (currents < 10 m/s, salinity < 50).
        series_lookup = {
            "wl": _IMPL.np.array([10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0], dtype=float) - 10.0,
            "u": _IMPL.np.array([2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6], dtype=float),
            "v": _IMPL.np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6], dtype=float),
            "sal": _IMPL.np.array([40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0], dtype=float),
            "sed": _IMPL.np.array([50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0], dtype=float),
        }

        orig_get_preset = _IMPL.preset_module.get_bc_preset
        orig_load_model = _IMPL._load_model_definition
        orig_load_bnd = _IMPL._load_boundary_locations
        orig_open_dataset = _IMPL._open_dataset
        orig_build_qc = _IMPL._build_qc_context
        orig_prepare = _prepare_boundary_sampling
        orig_check_cancel = _IMPL._check_cancel
        orig_interp = _IMPL.interp.bilinear_interpolation
        orig_write_bct = _IMPL.writer_delft3d4.write_bct
        orig_write_bct_currents = _IMPL.writer_delft3d4.write_bct_currents
        orig_write_bcc = _IMPL.writer_delft3d4.write_bcc
        orig_write_bcc_salinity = _IMPL.writer_delft3d4.write_bcc_salinity
        orig_write_bc = _IMPL.writer_delft3dfm.write_bc
        captured: Dict[str, Any] = {}

        try:
            _IMPL.preset_module.get_bc_preset = lambda _name: DummyPreset()
            _IMPL._load_model_definition = lambda *args, **kwargs: (timing, "Cartesian")
            _IMPL._load_boundary_locations = lambda *args, **kwargs: (boundary_df.copy(), None, list(pair_data))
            _IMPL._open_dataset = lambda *args, **kwargs: probe_ds
            _IMPL._build_qc_context = lambda *args, **kwargs: {}
            globals()["_prepare_boundary_sampling"] = lambda *args, **kwargs: (
                args[3],
                args[4],
                {"repaired_points": [], "rejected_boundaries": []},
            )
            _IMPL._check_cancel = lambda *_args, **_kwargs: None
            _IMPL.interp.bilinear_interpolation = (
                lambda lat, lon, ds_obj, variable=None, lat_coord=None, lon_coord=None, stats=None, cache=None: series_lookup[variable].copy()
            )
            _IMPL.writer_delft3d4.write_bct = lambda *args, **kwargs: captured.update(
                {"kind": "write_bct", "path": args[0], "boundary_data": args[1], "timing_args": args[2:6]}
            )
            _IMPL.writer_delft3d4.write_bct_currents = lambda *args, **kwargs: captured.update(
                {"kind": "write_bct_currents", "path": args[0], "boundary_data": args[1], "timing_args": args[2:6]}
            )
            _IMPL.writer_delft3d4.write_bcc = lambda *args, **kwargs: captured.update(
                {"kind": "write_bcc", "path": args[0], "boundary_data": args[1], "timing_args": args[2:7]}
            )
            _IMPL.writer_delft3d4.write_bcc_salinity = lambda *args, **kwargs: captured.update(
                {"kind": "write_bcc_salinity", "path": args[0], "boundary_data": args[1], "timing_args": args[2:6]}
            )
            _IMPL.writer_delft3dfm.write_bc = lambda *args, **kwargs: captured.update(
                {
                    "kind": "write_bc",
                    "path": args[0],
                    "boundary_data": args[1],
                    "quantity": args[2],
                    "unit": args[3],
                    "reference": args[4],
                }
            )
            captured["result"] = globals()[extractor_name](
                model_type,
                "flow.mdf",
                _probe_bnd_path,
                grid_path="grid.grd",
                data_files=["dummy.nc"],
                preset_name="Dummy",
                timestep=timestep_value,
                output_dir="out",
                logger=logger,
            )
        finally:
            _IMPL.preset_module.get_bc_preset = orig_get_preset
            _IMPL._load_model_definition = orig_load_model
            _IMPL._load_boundary_locations = orig_load_bnd
            _IMPL._open_dataset = orig_open_dataset
            _IMPL._build_qc_context = orig_build_qc
            globals()["_prepare_boundary_sampling"] = orig_prepare
            _IMPL._check_cancel = orig_check_cancel
            _IMPL.interp.bilinear_interpolation = orig_interp
            _IMPL.writer_delft3d4.write_bct = orig_write_bct
            _IMPL.writer_delft3d4.write_bct_currents = orig_write_bct_currents
            _IMPL.writer_delft3d4.write_bcc = orig_write_bcc
            _IMPL.writer_delft3d4.write_bcc_salinity = orig_write_bcc_salinity
            _IMPL.writer_delft3dfm.write_bc = orig_write_bc

        return captured

    wl_probe = run_non_wave_timestep_probe("extract_water_level", "delft3d4", 60.0)
    check("water level D3D4 resampled row count", wl_probe["result"]["n_timesteps"], 3)
    check("water level D3D4 sampled values", wl_probe["boundary_data"]["BND"]["wl_a"], [0.0, 3.0, 6.0])

    current_probe = run_non_wave_timestep_probe("extract_currents", "delft3dfm", 60.0)
    fm_current = current_probe["boundary_data"]["BND_a"]
    check("currents FM minute axis uses requested timestep", fm_current["time"], [0.0, 60.0, 120.0])
    check("currents FM sampled u values", [round(v, 9) for v in fm_current["values_x"]], [2.0, 2.3, 2.6])
    check("currents FM sampled v values", [round(v, 9) for v in fm_current["values_y"]], [3.0, 3.3, 3.6])

    sal_probe = run_non_wave_timestep_probe("extract_salinity", "delft3d4", 60.0)
    check("salinity D3D4 sampled values", sal_probe["boundary_data"]["BND"]["sal_a"], [40.0, 43.0, 46.0])

    sed_probe = run_non_wave_timestep_probe("extract_sediments", "delft3dfm", 60.0)
    # FM sediment blocks are keyed <point>__<fraction>; the emitted INI Name is
    # the bare support point and the quantity is sedfracbnd<fraction>.
    fm_sed = sed_probe["boundary_data"]["BND_a__mud"]
    check("sediments FM block Name is the support point", fm_sed["name"], "BND_a")
    check("sediments FM quantity is sedfracbnd<frac>", fm_sed["quantity"], "sedfracbndmud")
    check("sediments FM minute axis uses requested timestep", fm_sed["time"], [0.0, 60.0, 120.0])
    check("sediments FM sampled values", fm_sed["values"], [50.0, 53.0, 56.0])

    # ── Multi-source (cross-border) fallback ─────────────────────────────
    def make_region_dataset(lon_values, fill_value):
        time = pd.date_range("2020-01-01", periods=4, freq="h")
        lat = [0.0, 1.0]
        values = _IMPL.np.full((len(time), len(lat), len(lon_values)), fill_value, dtype=float)
        return xr.Dataset(
            {"sea_surface_height_2d": (("time", "lat", "lon"), values)},
            coords={"time": time, "lat": lat, "lon": list(lon_values)},
        )

    def make_source(label, ds):
        qc = _IMPL._build_qc_context(ds, ["sea_surface_height_2d"], TRILAWATT)
        qc["_intertidal_tolerance"] = False
        return {
            "label": label,
            "preset": TRILAWATT,
            "variables": ["sea_surface_height_2d"],
            "ds": ds,
            "qc_context": qc,
            "native_minutes": _IMPL.np.array([0.0, 60.0, 120.0, 180.0]),
        }

    ms_stream = io.StringIO()
    ms_logger = logging.getLogger("backend_bc_ms_test")
    ms_logger.handlers[:] = []
    ms_logger.addHandler(logging.StreamHandler(ms_stream))
    ms_logger.setLevel(logging.INFO)

    ms_sources = [
        make_source("west", make_region_dataset([0.0, 1.0], 1.0)),
        make_source("east", make_region_dataset([2.0, 3.0], 2.0)),
    ]
    ms_points = pd.DataFrame(
        [
            {"boundaries": "A_a", "lat": 0.5, "lon": 0.2},
            {"boundaries": "A_b", "lat": 0.5, "lon": 0.8},
            {"boundaries": "B_a", "lat": 0.5, "lon": 2.2},
            {"boundaries": "B_b", "lat": 0.5, "lon": 2.8},
            {"boundaries": "C_a", "lat": 0.5, "lon": 4.9},
            {"boundaries": "C_b", "lat": 0.5, "lon": 2.6},
            {"boundaries": "D_a", "lat": 0.5, "lon": 0.5},
            {"boundaries": "D_b", "lat": 0.5, "lon": 2.4},
        ]
    )
    ms_geo, _ms_pairs, ms_summary, ms_assign = _prepare_multi_source_sampling(
        "delft3d4", "water level", ms_sources, ms_points, [], ms_logger,
        max_repair_distance_m=1000.0,
    )
    check("multi-source west assignment", ms_assign["A_a"], 0)
    check("multi-source east assignment", ms_assign["B_b"], 1)
    check("straddling section split assignment", (ms_assign["D_a"], ms_assign["D_b"]), (0, 1))
    check("uncovered section rejected under repair cap", ms_summary["rejected_boundaries"], ["C"])
    check("rejected section removed from sampling", "C_a" in set(ms_geo["boundaries"]), False)
    ms_log_text = ms_stream.getvalue()
    check("coverage report logged", "Coverage report" in ms_log_text, True)
    check("no-source point flagged", "NO SOURCE" in ms_log_text, True)
    check("straddle logged", "section D straddles sources" in ms_log_text, True)
    check("source labels recorded", ms_summary["source_assignments"]["A_b"], "west")

    # Junction guard: two near-coincident corner endpoints must not switch
    # source when one source can serve both.
    pin_stream = io.StringIO()
    pin_logger = logging.getLogger("backend_bc_pin_test")
    pin_logger.handlers[:] = []
    pin_logger.addHandler(logging.StreamHandler(pin_stream))
    pin_logger.setLevel(logging.INFO)
    pin_sources = [
        make_source("narrow", make_region_dataset([0.0, 1.0], 1.0)),
        make_source("wide", make_region_dataset([0.0, 1.5, 3.0], 5.0)),
    ]
    pin_points = pd.DataFrame(
        [
            {"boundaries": "H_a", "lat": 0.5, "lon": 0.2},
            {"boundaries": "H_b", "lat": 0.5, "lon": 0.999},
            {"boundaries": "I_a", "lat": 0.5, "lon": 1.001},
            {"boundaries": "I_b", "lat": 0.5, "lon": 2.5},
        ]
    )
    _pin_geo, _pin_pairs, _pin_summary, pin_assign = _prepare_multi_source_sampling(
        "delft3d4", "water level", pin_sources, pin_points, [], pin_logger,
    )
    check("junction endpoints pinned to common source", (pin_assign["H_b"], pin_assign["I_a"]), (1, 1))
    check("junction pinning logged", "pinned to" in pin_stream.getvalue(), True)

    seam_flagged = _report_seam_consistency(
        "water level",
        {"H_b": [1.0, 1.0, 1.0], "I_a": [1.35, 1.35, 1.35]},
        {"H_b": (0.5, 0.999), "I_a": (0.5, 1.001)},
        {"H_b": 0, "I_a": 1},
        ["narrow", "wide"],
        pin_logger,
        threshold=0.10,
    )
    check("seam report flags 0.35 jump at corner", len(seam_flagged), 1)
    seam_ok = _report_seam_consistency(
        "water level",
        {"H_b": [1.0, 1.0, 1.0], "I_a": [1.02, 1.02, 1.02]},
        {"H_b": (0.5, 0.999), "I_a": (0.5, 1.001)},
        {"H_b": 0, "I_a": 0},
        ["narrow", "wide"],
        pin_logger,
        threshold=0.10,
    )
    check("seam report passes small differences", len(seam_ok), 0)

    check(
        "quantity variables for waves (vector preset)",
        _quantity_variables(TRILAWATT, "waves")[:2],
        [TRILAWATT.wave_hs_var, TRILAWATT.wave_tp_var],
    )

    # Per-group Custom overrides: each group gets its own cloned preset and
    # the shared Custom singleton must stay untouched.
    def make_named_var_dataset(var_name):
        time = pd.date_range("2018-01-01", periods=4, freq="h")
        values = _IMPL.np.ones((4, 2, 2), dtype=float)
        return xr.Dataset(
            {var_name: (("time", "lat", "lon"), values)},
            coords={"time": time, "lat": [0.0, 1.0], "lon": [0.0, 1.0]},
        )

    custom_singleton = _IMPL.preset_module.get_bc_preset("Custom")
    singleton_wl_before = custom_singleton.water_level_var
    override_timing = {
        "reference_date": "20180101",
        "reference_date_iso": "2018-01-01",
        "tstart": 0.0,
        "tstop": 180.0,
        "start_datetime": pd.Timestamp("2018-01-01 00:00:00"),
        "end_datetime": pd.Timestamp("2018-01-01 03:00:00"),
    }
    orig_open_dataset_ms = _IMPL._open_dataset
    try:
        _IMPL._open_dataset = (
            lambda files, preset, timing, log: make_named_var_dataset(preset.water_level_var)
        )
        override_sources = _build_boundary_sources(
            [
                {
                    "data_files": ["region_one.nc"],
                    "preset_name": "Custom",
                    "preset_overrides": {
                        "time_coord": "time", "lat_coord": "lat", "lon_coord": "lon",
                        "water_level_var": "wl_region_one",
                    },
                },
                {
                    "data_files": ["region_two.nc"],
                    "preset_name": "Custom",
                    "preset_overrides": {
                        "time_coord": "time", "lat_coord": "lat", "lon_coord": "lon",
                        "water_level_var": "wl_region_two",
                    },
                },
            ],
            "water_level", "delft3d4", override_timing, ms_logger,
        )
    finally:
        _IMPL._open_dataset = orig_open_dataset_ms
    check(
        "custom groups keep independent variable mappings",
        (override_sources[0]["variables"], override_sources[1]["variables"]),
        (["wl_region_one"], ["wl_region_two"]),
    )
    check(
        "custom preset singleton untouched by group overrides",
        _IMPL.preset_module.get_bc_preset("Custom").water_level_var,
        singleton_wl_before,
    )
    try:
        _build_boundary_sources(
            [{"data_files": ["x.nc"], "preset_name": "Custom", "preset_overrides": {"not_a_field": 1}}],
            "water_level", "delft3d4", override_timing, ms_logger,
        )
        check("unknown override field rejected", "no error", "ValueError")
    except ValueError:
        check("unknown override field rejected", "ValueError", "ValueError")

    # ── Neumann sections must never be sampled as water levels ───────────
    _probe_bnd_dir.cleanup()
    with tempfile.TemporaryDirectory() as tmp:
        bnd_mixed = Path(tmp) / "mixed.bnd"
        bnd_mixed.write_text(
            "North01              N T  1104   921  1104   514  0.0000000e+000\n"
            "West_1               Z T     2   922    31   922  1.0000000e+002\n"
            "West_2               Z T    32   922    61   922  1.0000000e+002\n"
            "South01              N T     1   550     1   921  0.0000000e+000\n"
        )
        n_names = neumann_section_names("delft3d4", str(bnd_mixed))
        check("N sections detected", n_names, {"North01", "South01"})
        check("FM (.pli) has no typed sections",
              neumann_section_names("delft3dfm", str(bnd_mixed)), set())

        guard_df = pd.DataFrame({
            "lat": [54.0] * 8, "lon": [8.0] * 8,
            "boundaries": ["North01_a", "North01_b", "West_1_a", "West_1_b",
                           "West_2_a", "West_2_b", "South01_a", "South01_b"],
        })
        guard_pairs = [{"name": nm} for nm in ("North01", "West_1", "West_2", "South01")]
        kept_df, kept_pairs = _drop_neumann_points(
            guard_df, guard_pairs, n_names, ms_logger
        )
        check("neumann points dropped", len(kept_df), 4)
        check("only Z sections survive",
              sorted(kept_df["boundaries"]),
              ["West_1_a", "West_1_b", "West_2_a", "West_2_b"])
        check("bnd_pairs filtered in step",
              [p["name"] for p in kept_pairs], ["West_1", "West_2"])
        check("index reset after drop", list(kept_df.index), [0, 1, 2, 3])

        untouched_df, untouched_pairs = _drop_neumann_points(
            guard_df, guard_pairs, set(), ms_logger
        )
        check("no N sections -> untouched", len(untouched_df), 8)
        check("no N sections -> pairs untouched", len(untouched_pairs), 4)

    print(f"\n{'=' * 40}")
    print("ALL TESTS PASSED" if all_pass else "SOME TESTS FAILED")
    print(f"{'=' * 40}")
