# -*- coding: utf-8 -*-
"""
backend_neumann.py — Generate Delft3D-4 Neumann (alongshore water-level
gradient) boundary forcing from a parent hydrodynamic model.

Physics (Delft3D-FLOW User Manual §9.4.1, Roelvink & Walstra 2004):
    Cross-shore lateral boundaries can be forced with the alongshore
    water-level gradient g(t) = dzeta/dx [-] instead of absolute water
    levels, in combination with ONE offshore water-level boundary
    (reflection alpha = 0) that anchors the absolute level. Both lateral
    boundaries take the gradient in the SAME fixed model direction (the
    manual's Table 9.4 applies identical phase/amplitude relations to both
    sides), so the sign convention here is the model grid's positive index
    direction normal to the section — not the outward normal.

File formats (manual A.2.10 and A.2.15):
    .bnd type letter 'N', no reflection coefficient on N records.
    .bct tables identical in structure to water-level tables, with
    parameter records "neumann         (n)  end A/B", unit '[-]'.
    The manual further recommends Cstbnd = #YES# in the .mdf to switch
    off normal-gradient advection terms at the boundaries (§9.4.1).

Method:
    For every unique boundary support point (section endpoints A/B, shared
    corners deduplicated so adjacent sections get bit-identical series), a
    small probe cloud is laid out around the point in UTM metres. Water
    levels are sampled at each probe from the parent model via the BC
    extraction engine (presets, QC, multi-year stitching and cross-border
    source_groups all apply). Per timestep a plane
        eta(x, y, t) = a(t) + b(t)*x + c(t)*y
    is fitted over the cloud (linear least squares, one precomputed
    pseudo-inverse applied to all timesteps at once) and the gradient
    (b, c) is projected onto the alongshore unit vector. The 'centered'
    method is the two-probe finite-difference special case.

    Probes are screened with strict QC only (never the intertidal tolerance):
    a probe on a drying flat cannot feed a gradient fit. When a coastal
    support point is left with too few probes, the point is slid along its
    section toward the other end until the cloud resolves, and the shift is
    reported loudly; Delft3D interpolates the gradient linearly between the
    section ends, so the end carries the nearest gradient that can be fitted.

Source data is the parent-model NetCDF only (preferred sources 1-2 of the
requirements). Datum/timezone consistency is the preset's responsibility
(time_offset_hours); all probes are sampled from the same source machinery
so they share datum and clock by construction.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

import backend_bc_extractor as bce
import coordinate_utils
import file_parsers
import writer_delft3d4

logger_default = logging.getLogger("repd3d.neumann")

PROBE_ALONG_SPAN_FACTOR = 4.0   # default alongshore probe span = 4 x parent grid spacing
PROBE_CROSS_SPAN_FACTOR = 2.0   # default cross-shore probe span = 2 x parent grid spacing
NEUMANN_SPIKE_THRESHOLD = 5.0e-4  # [-] typical tide/surge gradients are 1e-6..1e-5


# ─────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ─────────────────────────────────────────────────────────────────────────────

def _unit(vector: Tuple[float, float]) -> Tuple[float, float]:
    norm = math.hypot(vector[0], vector[1])
    if norm == 0.0:
        raise ValueError("Zero-length direction vector in grid geometry.")
    return (vector[0] / norm, vector[1] / norm)


def _grid_xy(grd: Dict[str, Any], mi: int, ni: int) -> Tuple[float, float]:
    """Read (x, y) at 1-based (m, n), mirroring file_parsers indexing."""
    return (
        float(grd["x_values"][ni - 1][mi]),
        float(grd["y_values"][ni - 1][mi]),
    )


def _is_missing(value: float) -> bool:
    return value == 0.0 or abs(value) > 1.0e8 or abs(value + 999.999) < 1.0


def _section_alongshore_direction(
    grd: Dict[str, Any], bnd: Dict[str, Any], flip_sign: bool
) -> Tuple[float, float]:
    """
    Alongshore unit vector for a cross-shore section: the model grid's
    positive index direction normal to the section (positive m for
    constant-m sections, positive n for constant-n sections).
    """
    m_a, n_a, m_b, n_b = bnd["m_a"], bnd["n_a"], bnd["m_b"], bnd["n_b"]
    m_max, n_max = grd["m"], grd["n"]
    mid_m = min(max((m_a + m_b) // 2, 1), m_max)
    mid_n = min(max((n_a + n_b) // 2, 1), n_max)

    if m_a == m_b:      # section runs along n → alongshore is +m
        lo, hi = max(mid_m - 1, 1), min(mid_m + 1, m_max)
        candidates = [(lo, hi)] + [(mm, mm + 1) for mm in range(max(lo - 3, 1), min(hi + 3, m_max))]
        axis = "m"
    elif n_a == n_b:    # section runs along m → alongshore is +n
        lo, hi = max(mid_n - 1, 1), min(mid_n + 1, n_max)
        candidates = [(lo, hi)] + [(nn, nn + 1) for nn in range(max(lo - 3, 1), min(hi + 3, n_max))]
        axis = "n"
    else:
        raise ValueError(
            f"Section '{bnd['bnd_name']}' is diagonal in grid indices "
            f"(m {m_a}->{m_b}, n {n_a}->{n_b}); Neumann generation supports "
            "constant-m or constant-n sections only."
        )

    for lo, hi in candidates:
        if lo >= hi:
            continue
        if axis == "m":
            x1, y1 = _grid_xy(grd, lo, mid_n)
            x2, y2 = _grid_xy(grd, hi, mid_n)
        else:
            x1, y1 = _grid_xy(grd, mid_m, lo)
            x2, y2 = _grid_xy(grd, mid_m, hi)
        if any(_is_missing(v) for v in (x1, y1, x2, y2)):
            continue
        direction = _unit((x2 - x1, y2 - y1))
        if flip_sign:
            direction = (-direction[0], -direction[1])
        return direction
    raise ValueError(
        f"Could not derive the grid direction for section '{bnd['bnd_name']}' "
        "(missing grid coordinates near the section)."
    )


def _azimuth_deg(direction: Tuple[float, float]) -> float:
    """Compass azimuth (deg from north, clockwise) of a UTM direction."""
    return (math.degrees(math.atan2(direction[0], direction[1])) + 360.0) % 360.0


# ─────────────────────────────────────────────────────────────────────────────
# Gradient estimation
# ─────────────────────────────────────────────────────────────────────────────

def _fit_gradient_series(
    probe_offsets_m: np.ndarray,   # (n_probes, 2) UTM offsets from support point
    probe_series: np.ndarray,      # (n_probes, n_times) water levels
    e_along: Tuple[float, float],
) -> np.ndarray:
    """
    Fit eta = a + b*dx + c*dy per timestep and project (b, c) onto e_along.
    Falls back to a 1-D alongshore fit when the probe cloud is collinear.
    Returns g(t) with shape (n_times,).
    """
    n_probes = probe_offsets_m.shape[0]
    if n_probes < 2:
        raise ValueError("Gradient fit needs at least 2 probes.")
    design = np.column_stack([np.ones(n_probes), probe_offsets_m])
    if np.linalg.matrix_rank(design) >= 3:
        coeffs = np.linalg.pinv(design) @ probe_series      # (3, n_times)
        return coeffs[1] * e_along[0] + coeffs[2] * e_along[1]
    # Collinear cloud: fit along the alongshore coordinate only.
    s = probe_offsets_m @ np.asarray(e_along)
    if len(np.unique(np.round(s, 3))) < 2:
        raise ValueError(
            "Probe cloud is collinear with no alongshore separation - "
            "cannot estimate the alongshore gradient."
        )
    design_1d = np.column_stack([np.ones(n_probes), s])
    coeffs = np.linalg.pinv(design_1d) @ probe_series       # (2, n_times)
    return coeffs[1]


def _probe_offsets(
    method: str,
    e_along: Tuple[float, float],
    along_span_m: float,
    cross_span_m: float,
) -> np.ndarray:
    """Probe cloud offsets in UTM metres around a support point."""
    e_cross = (-e_along[1], e_along[0])
    if method == "centered":
        deltas = [(-along_span_m / 2.0, 0.0), (along_span_m / 2.0, 0.0)]
    else:
        deltas = [
            (da, dc)
            for da in (-along_span_m / 2.0, 0.0, along_span_m / 2.0)
            for dc in (-cross_span_m / 2.0, 0.0, cross_span_m / 2.0)
        ]
    return np.array(
        [
            (
                da * e_along[0] + dc * e_cross[0],
                da * e_along[1] + dc * e_cross[1],
            )
            for da, dc in deltas
        ]
    )


NEUMANN_SLIDE_STEPS = 6   # bisection steps when sliding a support point along its section


def _find_usable_support_position(
    easting: float,
    northing: float,
    targets: List[Tuple[str, str, float, float]],
    evaluate: Callable[[float, float], List[Dict[str, Any]]],
    min_needed: int,
    max_distance_m: Optional[float] = None,
    n_steps: int = NEUMANN_SLIDE_STEPS,
) -> Optional[Dict[str, Any]]:
    """
    Slide a support point along its section toward the other end until the
    probe cloud there has at least ``min_needed`` usable probes.

    ``targets`` holds (section_name, target_name, easting, northing) per owning
    section; ``evaluate(easting, northing)`` returns the usable probes of the
    cloud centred there. Bisection over the fraction of the section length
    (the same scheme as the BC engine's point repair) finds the smallest
    passing displacement per target; the closest one across all targets wins.
    Returns None when no position within ``max_distance_m`` passes.
    """
    best: Optional[Dict[str, Any]] = None
    for section_name, target_name, t_east, t_north in targets:
        length = math.hypot(t_east - easting, t_north - northing)
        if length == 0.0:
            continue
        at_target = evaluate(t_east, t_north)
        best_here: Optional[Tuple[float, List[Dict[str, Any]]]] = (
            (1.0, at_target) if len(at_target) >= min_needed else None
        )
        lo, hi = 0.0, 1.0
        for _ in range(n_steps):
            frac = (lo + hi) / 2.0
            cand = evaluate(easting + frac * (t_east - easting), northing + frac * (t_north - northing))
            if len(cand) >= min_needed:
                best_here = (frac, cand)
                hi = frac
            else:
                lo = frac
        if best_here is None:
            continue
        frac, cand = best_here
        distance = frac * length
        if max_distance_m is not None and distance > max_distance_m:
            continue
        if best is None or distance < best["distance_m"]:
            best = {
                "easting": easting + frac * (t_east - easting),
                "northing": northing + frac * (t_north - northing),
                "distance_m": distance,
                "fraction": frac,
                "section": section_name,
                "target_name": target_name,
                "usable": cand,
            }
    return best


def _parent_grid_spacing_m(source: Dict[str, Any]) -> float:
    """Median parent grid spacing in metres (mean of lat/lon spacing)."""
    preset = source["preset"]
    ds = source["ds"]
    lats = np.asarray(ds[preset.lat_coord].values, dtype=float)
    lons = np.asarray(ds[preset.lon_coord].values, dtype=float)
    lat_mid = float(np.median(lats))
    dlat = float(np.median(np.abs(np.diff(lats)))) if len(lats) > 1 else 0.0
    dlon = float(np.median(np.abs(np.diff(lons)))) if len(lons) > 1 else 0.0
    spacing_lat = dlat * 111320.0
    spacing_lon = dlon * 111320.0 * math.cos(math.radians(lat_mid))
    spacing = np.mean([s for s in (spacing_lat, spacing_lon) if s > 0.0])
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("Could not derive parent grid spacing from the source dataset.")
    return float(spacing)


# ─────────────────────────────────────────────────────────────────────────────
# Main generator
# ─────────────────────────────────────────────────────────────────────────────

def generate_neumann_boundaries(
    model_type: str,
    model_def_path: str,
    boundary_path: str,
    grid_path: str,
    data_files: Optional[List[str]] = None,
    preset_name: str = "TrilaWatt",
    source_groups: Optional[List[Dict[str, Any]]] = None,
    section_names: Optional[List[str]] = None,
    method: str = "plane",
    along_span_m: Optional[float] = None,
    cross_span_m: Optional[float] = None,
    flip_sign: bool = False,
    spike_threshold: float = NEUMANN_SPIKE_THRESHOLD,
    timestep: float = 20.0,
    utm_zone: int = 32,
    hemisphere: str = "N",
    output_path: Optional[str] = None,
    output_dir: str = ".",
    qc_plot: bool = True,
    include_water_level_sections: bool = True,
    max_repair_distance_m: Optional[float] = None,
    logger=None,
    progress_callback: Callable[[float], None] = lambda x: None,
    cancel_event=None,
) -> Dict[str, Any]:
    """
    Generate g(s,t) Neumann forcing for the cross-shore lateral boundary
    sections of a Delft3D-4 model, from a parent-model NetCDF (single
    source or cross-border source_groups), and write an N-type .bct.

    Sections are taken from the .bnd: those typed 'N', or an explicit
    section_names list (with a loud warning if their .bnd type is not yet
    'N'). Shared corners of adjacent sections receive bit-identical series
    by construction.

    With include_water_level_sections=True (the default) the remaining
    sections of the .bnd — the ordinary type-Z water-level boundaries — are
    extracted too and written into the SAME .bct, in .bnd order. Delft3D
    binds tables to sections by name and needs one .bct covering every
    section, so a Neumann-only file leaves the model unrunnable. Set it to
    False to get the gradient sections alone.
    """
    if logger is None:
        logger = logger_default
    if model_type != "delft3d4":
        raise ValueError(
            "Neumann gradient generation currently supports Delft3D-4 only "
            "(.bnd type 'N' + .bct tables)."
        )
    if method not in ("plane", "centered"):
        raise ValueError(f"Unknown method '{method}' (use 'plane' or 'centered').")

    logger.info(f"Starting Neumann gradient generation — method: {method}")
    bce._check_cancel(cancel_event)

    timing, crs_type = bce._load_model_definition(model_type, model_def_path, grid_path, logger)
    bce._check_cancel(cancel_event)

    bnd_list = file_parsers.parse_bnd(boundary_path)
    bnd_list = file_parsers.extract_boundary_coordinates(grid_path, bnd_list)
    grd = file_parsers.parse_grd(grid_path)

    # ── Select target sections ──────────────────────────────────────────
    if section_names:
        wanted = {name.strip() for name in section_names}
        targets = [b for b in bnd_list if b["bnd_name"] in wanted]
        missing = wanted - {b["bnd_name"] for b in targets}
        if missing:
            raise ValueError(f"Sections not found in .bnd: {', '.join(sorted(missing))}")
    else:
        targets = [b for b in bnd_list if b["bnd_type"].upper() == "N"]
        if not targets:
            types_found = sorted({b["bnd_type"] for b in bnd_list})
            raise ValueError(
                "No sections typed 'N' (Neumann) in the .bnd and no section_names "
                f"given. Types found: {', '.join(types_found)}. Either set the lateral "
                "sections' type to N in the .bnd, or pass section_names explicitly."
            )
    for b in targets:
        if b["bnd_type"].upper() != "N":
            logger.warning(
                f"[neumann] Section {b['bnd_name']} is typed '{b['bnd_type']}' in the .bnd. "
                "Delft3D will only interpret the values as gradients if the type is 'N' "
                "(and N records take NO reflection coefficient — manual A.2.10)."
            )
    logger.info(f"[neumann] Target sections: {', '.join(b['bnd_name'] for b in targets)}")

    # ── Coordinates: UTM + geographic per endpoint ──────────────────────
    is_spherical = str(crs_type).lower().startswith("spher")
    supports: Dict[Tuple[float, float], Dict[str, Any]] = {}
    section_meta: List[Dict[str, Any]] = []
    for b in targets:
        e_along = _section_alongshore_direction(grd, b, flip_sign)
        entry = {"bnd": b, "e_along": e_along, "ends": {}}
        logger.info(
            f"[neumann] {b['bnd_name']}: alongshore direction {_azimuth_deg(e_along):.1f} degN"
            f"{' (sign flipped)' if flip_sign else ''}"
        )
        for end in ("a", "b"):
            x, y = b[f"x_{end}"], b[f"y_{end}"]
            if is_spherical:
                lon, lat = x, y
                easting, northing = coordinate_utils.latlon_to_utm(lat, lon, zone=utm_zone, hemisphere=hemisphere)
            else:
                easting, northing = x, y
                lat, lon = coordinate_utils.utm_to_latlon(easting, northing, zone=utm_zone, hemisphere=hemisphere)
            key = (round(easting, 1), round(northing, 1))
            if key not in supports:
                supports[key] = {
                    "easting": easting, "northing": northing,
                    "lat": lat, "lon": lon,
                    "e_along": e_along,
                    "names": [],
                }
            supports[key]["names"].append(f"{b['bnd_name']}_{end}")
            entry["ends"][end] = key
        section_meta.append(entry)
    logger.info(
        f"[neumann] {len(supports)} unique support point(s) for "
        f"{len(targets)} section(s) (shared corners deduplicated)."
    )

    # ── Sources (single or cross-border groups) ─────────────────────────
    groups = source_groups or [{"data_files": list(data_files or []), "preset_name": preset_name}]
    sources = bce._build_boundary_sources(groups, "water_level", model_type, timing, logger, cancel_event=cancel_event)
    labels = [s["label"] for s in sources]
    spacing_m = _parent_grid_spacing_m(sources[0])
    span_along = float(along_span_m or PROBE_ALONG_SPAN_FACTOR * spacing_m)
    span_cross = float(cross_span_m or PROBE_CROSS_SPAN_FACTOR * spacing_m)
    logger.info(
        f"[neumann] Parent grid spacing ~{spacing_m:.0f} m; probe spans: "
        f"alongshore {span_along:.0f} m, cross-shore {span_cross:.0f} m."
    )

    output_minutes = bce._build_output_minutes(model_type, timing, timestep)
    n_times = len(output_minutes)

    # ── Pass 1: QC every probe before sampling anything ─────────────────
    # _evaluate_point_status reads precomputed validity masks and is cheap;
    # _resample_scalar_variable reads a full time series per probe and is not.
    # Interleaving them meant an unusable support point was only discovered
    # after every earlier point had been sampled — a multi-hour wait to learn
    # something knowable up front. Screen all probes first, then report EVERY
    # unusable point at once.
    probes_per_point = len(_probe_offsets(method, (1.0, 0.0), span_along, span_cross))
    total_probes = max(1, len(supports) * probes_per_point)
    min_needed = 2 if method == "centered" else 3
    logger.info(
        f"[neumann] Screening {total_probes} probe(s) — {probes_per_point} per support "
        f"point x {len(supports)} point(s)."
    )

    def _screen_probes(easting: float, northing: float, e_along: Tuple[float, float]) -> List[Dict[str, Any]]:
        """Strict-QC the probe cloud around (easting, northing); return the usable probes."""
        usable: List[Dict[str, Any]] = []
        for d_east, d_north in _probe_offsets(method, e_along, span_along, span_cross):
            p_east, p_north = easting + d_east, northing + d_north
            p_lat, p_lon = coordinate_utils.utm_to_latlon(p_east, p_north, zone=utm_zone, hemisphere=hemisphere)
            for src_idx, src in enumerate(sources):
                # allow_intertidal=False: a probe whose stencil dries cannot feed a
                # gradient fit. Dry periods bridged at different times at different
                # probes read as metre-scale alongshore differences, i.e. spurious
                # gradients (South_1_a, 2026-10-13: +-1.9e-4 against +-1.4e-5 offshore).
                status = bce._evaluate_point_status(
                    p_lat, p_lon, src["ds"], src["variables"], src["preset"], src["qc_context"],
                    allow_intertidal=False,
                )
                if status["reliable"]:
                    usable.append({
                        "offset": (d_east, d_north), "source": src_idx,
                        "lat": p_lat, "lon": p_lon,
                    })
                    break
        return usable

    owners: Dict[Tuple[float, float], List[Tuple[Dict[str, Any], str]]] = {}
    for entry in section_meta:
        for end in ("a", "b"):
            owners.setdefault(entry["ends"][end], []).append((entry, end))

    probe_plan: Dict[Tuple[float, float], List[Dict[str, Any]]] = {}
    unusable: List[str] = []
    for idx, (key, sup) in enumerate(supports.items()):
        bce._check_cancel(cancel_event)
        usable = _screen_probes(sup["easting"], sup["northing"], sup["e_along"])
        point_label = "/".join(sup["names"])
        n_dropped = probes_per_point - len(usable)
        if n_dropped:
            logger.warning(
                f"[neumann] {point_label}: {n_dropped}/{probes_per_point} probe(s) failed strict QC "
                "in every source and were dropped (land, drying flat or missing data in the stencil)."
            )
        if len(usable) < min_needed:
            # Coastal end on a drying flat / land: slide the support point along
            # its section toward the other end until the cloud is resolvable.
            # Delft3D interpolates the gradient linearly between the two ends,
            # so the end receives the nearest gradient that can actually be fitted.
            slide_targets = []
            for entry, end in owners.get(key, []):
                other_end = "b" if end == "a" else "a"
                other = supports[entry["ends"][other_end]]
                slide_targets.append((
                    entry["bnd"]["bnd_name"], f"{entry['bnd']['bnd_name']}_{other_end}",
                    other["easting"], other["northing"],
                ))
            found = _find_usable_support_position(
                sup["easting"], sup["northing"], slide_targets,
                lambda e, n, _dir=sup["e_along"]: _screen_probes(e, n, _dir),
                min_needed, max_repair_distance_m,
            )
            if found is None:
                cap_note = (
                    f" within the {max_repair_distance_m:.0f} m repair cap"
                    if max_repair_distance_m is not None else ""
                )
                unusable.append(f"{point_label} ({len(usable)} usable; no resolvable position along the section{cap_note})")
            else:
                logger.warning(
                    f"[neumann] {point_label}: only {len(usable)}/{probes_per_point} probe(s) pass strict QC at the "
                    f"section end (coastal or intertidal stencil). Gradient evaluated {found['distance_m']:.0f} m "
                    f"inboard along {found['section']} toward {found['target_name']} "
                    f"({len(found['usable'])}/{probes_per_point} probes usable there). Delft3D interpolates linearly "
                    "between the section ends, so this end carries the nearest resolvable gradient."
                )
                if found["distance_m"] > bce.REPAIR_DISTANCE_ADVISORY_M:
                    logger.warning(
                        f"[neumann] {point_label}: the shift exceeds the {bce.REPAIR_DISTANCE_ADVISORY_M:.0f} m "
                        "advisory distance - check the QC plot and consider ending the section further offshore."
                    )
                usable = found["usable"]
                sup["easting"], sup["northing"] = found["easting"], found["northing"]
                sup["lat"], sup["lon"] = coordinate_utils.utm_to_latlon(
                    found["easting"], found["northing"], zone=utm_zone, hemisphere=hemisphere
                )
                sup["slid_m"] = found["distance_m"]
                sup["slid_toward"] = found["target_name"]
        probe_plan[key] = usable
        progress_callback(0.1 + 0.15 * (idx + 1) / max(1, len(supports)))

    if unusable:
        raise ValueError(
            f"{len(unusable)} support point(s) have fewer than {min_needed} usable "
            f"probe(s), so no gradient can be fitted there: {'; '.join(unusable)}. "
            "Nothing was sampled - fix this before rerunning. The generator already "
            "searched along each section toward the other end for a resolvable probe "
            "cloud. Options: enlarge along_span_m / shrink cross_span_m, use "
            "method='centered' (2 probes on the section line itself, no cross-shore "
            "offsets), raise max_repair_distance_m if one is set, extend the source "
            "coverage, or end the section further offshore."
        )

    # ── Pass 2: sample the surviving probes ─────────────────────────────
    gradient_by_key: Dict[Tuple[float, float], np.ndarray] = {}
    stats_rows: List[str] = []
    source_usage: Dict[str, int] = {}
    usable_total = max(1, sum(len(v) for v in probe_plan.values()))
    probes_done = 0
    logger.info(f"[neumann] Screening passed; extracting {usable_total} probe series...")
    # Prefetch every probe's stencil corners per source in one pass per
    # variable (see interpolation.CornerSeriesCache) instead of one full
    # read of the window per probe.
    probe_caches = [
        bce._prefetch_corner_series(
            src["ds"],
            [(p["lat"], p["lon"]) for plan in probe_plan.values() for p in plan if p["source"] == src_idx],
            src["variables"], src["preset"], logger, cancel_event,
        )
        for src_idx, src in enumerate(sources)
    ]

    for key, sup in supports.items():
        bce._check_cancel(cancel_event)
        point_label = "/".join(sup["names"])
        kept_offsets: List[Tuple[float, float]] = []
        kept_series: List[List[float]] = []
        probe_sources: List[str] = []
        for probe in probe_plan[key]:
            src = sources[probe["source"]]
            series = bce._resample_scalar_variable(
                {"lat": probe["lat"], "lon": probe["lon"], "boundaries": f"{point_label} probe"},
                src["ds"], src["preset"].water_level_var, src["preset"],
                src["native_minutes"], output_minutes, logger=logger,
                cache=probe_caches[probe["source"]],
            )
            kept_offsets.append(probe["offset"])
            kept_series.append(series)
            probe_sources.append(labels[probe["source"]])
            source_usage[labels[probe["source"]]] = (
                source_usage.get(labels[probe["source"]], 0) + 1
            )
            probes_done += 1
            progress_callback(0.25 + 0.35 * probes_done / usable_total)

        series_matrix = np.asarray(kept_series, dtype=float)
        if np.isnan(series_matrix).any():
            bad = int(np.isnan(series_matrix).sum())
            raise ValueError(
                f"Support point {point_label}: {bad} NaN sample(s) after QC-passed "
                "sampling — refusing to silently fill. Check the source time coverage."
            )
        g = _fit_gradient_series(np.asarray(kept_offsets), series_matrix, sup["e_along"])
        gradient_by_key[key] = g

        spikes = int((np.abs(g) > spike_threshold).sum())
        stats_rows.append(
            f"[neumann]   {point_label:<24} mean {np.mean(g):+.3e}  min {np.min(g):+.3e}  "
            f"max {np.max(g):+.3e}  RMS {np.sqrt(np.mean(g ** 2)):.3e}"
            + (f"  ({spikes} sample(s) above {spike_threshold:g}!)" if spikes else "")
        )
        if spikes:
            logger.warning(
                f"[neumann] {point_label}: {spikes} gradient sample(s) exceed "
                f"{spike_threshold:g} — inspect the QC plot before using this boundary."
            )
        if len(set(probe_sources)) > 1:
            logger.info(
                f"[neumann] {point_label}: probes span sources "
                f"{sorted(set(probe_sources))} (cross-border fit)."
            )

    probe_caches = None   # release the prefetched probe series before the section caches are built

    logger.info(f"[neumann] Gradient statistics [-] per support point:")
    for row in stats_rows:
        logger.info(row)

    # ── Water-level (type Z) sections of the same .bnd ───────────────────
    target_names = {b["bnd_name"] for b in targets}
    other_sections = [b for b in bnd_list if b["bnd_name"] not in target_names]
    wl_series: Dict[str, Dict[str, List[float]]] = {}
    wl_qc: Dict[str, Any] = {"repaired_points": [], "rejected_boundaries": []}

    if other_sections and include_water_level_sections:
        logger.info(
            f"[neumann] Extracting water levels for the {len(other_sections)} "
            "non-Neumann section(s) so the .bct covers the whole .bnd."
        )
        wl_series, wl_qc = _sample_water_level_sections(
            model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere,
            {b["bnd_name"] for b in other_sections}, sources, output_minutes,
            max_repair_distance_m, logger,
            lambda p: progress_callback(0.6 + 0.3 * p), cancel_event,
        )
    elif other_sections:
        logger.warning(
            f"[neumann] include_water_level_sections=False: the {len(other_sections)} "
            "non-Neumann section(s) are NOT in the output. Delft3D needs one .bct "
            "covering every .bnd section, so this file is incomplete on its own."
        )
    progress_callback(0.9)

    # ── Assemble in .bnd order and write ────────────────────────────────
    meta_by_name = {entry["bnd"]["bnd_name"]: entry for entry in section_meta}
    ordered_sections: List[Dict[str, Any]] = []
    incomplete: List[str] = []
    for b in bnd_list:
        name = b["bnd_name"]
        if name in meta_by_name:
            entry = meta_by_name[name]
            ordered_sections.append({
                "name": name,
                "kind": "neumann",
                "a": gradient_by_key[entry["ends"]["a"]].tolist(),
                "b": gradient_by_key[entry["ends"]["b"]].tolist(),
            })
        elif include_water_level_sections:
            ends = wl_series.get(name, {})
            if "a" in ends and "b" in ends:
                ordered_sections.append({
                    "name": name, "kind": "wl", "a": ends["a"], "b": ends["b"],
                })
            else:
                incomplete.append(name)

    if incomplete:
        logger.warning(
            f"[neumann] {len(incomplete)} water-level section(s) produced no usable "
            f"series and are missing from the .bct: {', '.join(incomplete)}. "
            "Delft3D will not run until every .bnd section has a table — extend the "
            "source coverage or relax the QC for these points."
        )

    if output_path is None:
        output_path = str(Path(output_dir) / f"{Path(model_def_path).stem}_neumann.bct")
    writer_delft3d4.write_bct_mixed(
        output_path, ordered_sections, timing["reference_date"],
        timing["tstart"], timing["tstop"], timestep,
    )
    n_wl = sum(1 for s in ordered_sections if s["kind"] == "wl")
    n_grad = len(ordered_sections) - n_wl
    logger.info(
        f"Output written: {output_path} "
        f"({n_grad} Neumann + {n_wl} water-level section(s))"
    )
    order_preview = [s["name"] for s in ordered_sections]
    if len(order_preview) > 6:
        order_preview = order_preview[:3] + ["..."] + order_preview[-2:]
    logger.info(
        f"[neumann] Table order mirrors the .bnd ({', '.join(order_preview)}). Delft3D binds "
        "'Boundary Section : k' positionally, so never re-order the .bct by hand - re-order the "
        ".bnd and regenerate. Parameter labels use the Delft3D-FLOW GUI spelling "
        "'neumann         (n)  end A' (9 + 2 spaces); FLOW rejects any other spacing."
    )

    # ── Mdf advice + validation checklist ───────────────────────────────
    try:
        mdf = file_parsers.parse_mdf(model_def_path)
        cstbnd = next((v for k, v in mdf.items() if k.lower() == "cstbnd"), "")
        if "yes" not in str(cstbnd).lower():
            logger.warning(
                "[neumann] The .mdf has no 'Cstbnd = #YES#'. The manual recommends it "
                "for water-level + Neumann coastal setups (switches off normal-gradient "
                "advection at the boundaries, §9.4.1)."
            )
    except Exception:
        pass
    logger.info(
        "[neumann] Validation checklist before enabling morphology: run a short "
        "fixed-bed simulation and verify (1) gradient sign — water should slope the "
        "same way as in the parent model, flip with flip_sign=True if not; (2) no "
        "velocity jets at the lateral boundaries; (3) smooth transitions across "
        "section corners; (4) stable boundary cells through low water; (5) offshore "
        "WL boundary uses reflection alpha = 0."
    )

    # ── QC plot ─────────────────────────────────────────────────────────
    output_files = [output_path]
    if qc_plot:
        try:
            plot_path = _write_qc_plot(
                Path(output_dir) / f"{Path(model_def_path).stem}_neumann_qc.png",
                supports, gradient_by_key, output_minutes, sources, span_along, span_cross,
            )
            output_files.append(plot_path)
            logger.info(f"[neumann] QC plot written: {plot_path}")
        except Exception as exc:  # plotting must never kill the run
            logger.warning(f"[neumann] QC plot skipped: {exc}")

    progress_callback(1.0)
    return {
        "output_files": output_files,
        "n_boundaries": len(ordered_sections),
        "n_timesteps": n_times,
        "repaired_points": wl_qc.get("repaired_points", []),
        "rejected_boundaries": wl_qc.get("rejected_boundaries", []),
        "source_assignments": {
            "/".join(sup["names"]): "mixed/see log" for sup in supports.values()
        } if len(sources) > 1 else {},
        "gradient_stats": stats_rows,
        "probe_source_usage": source_usage,
        "n_neumann_sections": n_grad,
        "n_water_level_sections": n_wl,
        "incomplete_sections": incomplete,
    }


def _sample_water_level_sections(
    model_type: str,
    boundary_path: str,
    grid_path: str,
    crs_type: str,
    utm_zone: int,
    hemisphere: str,
    wanted_names: set,
    sources: List[Dict[str, Any]],
    output_minutes,
    max_repair_distance_m: Optional[float],
    logger,
    progress_callback: Callable[[float], None],
    cancel_event,
) -> Tuple[Dict[str, Dict[str, List[float]]], Dict[str, Any]]:
    """
    Sample end A/B water levels for the ordinary (type Z) .bnd sections.

    Deliberately reuses the BC extraction engine's own sampling path rather
    than the Neumann probe loop: type-Z boundaries get the full QC treatment
    including point relocation, the coverage report and multi-source
    fallback, so a section here is served exactly as a standalone water-level
    run would serve it. `sources` is already priority-ordered, and the
    single-source case is just a one-entry list.

    Returns ({section_name: {'a': [...], 'b': [...]}}, qc_summary).
    """
    bnd_loc_geo, _native, bnd_pairs = bce._load_boundary_locations(
        model_type, boundary_path, grid_path, crs_type, utm_zone, hemisphere, logger
    )
    keep = bnd_loc_geo["boundaries"].map(lambda name: str(name)[:-2] in wanted_names)
    bnd_loc_geo = bnd_loc_geo[keep].reset_index(drop=True)
    bnd_pairs = [p for p in (bnd_pairs or []) if p["name"] in wanted_names]

    bnd_loc_geo, _pairs, qc_summary, point_sources = bce._prepare_multi_source_sampling(
        model_type, "water level", sources, bnd_loc_geo, bnd_pairs, logger,
        max_repair_distance_m=max_repair_distance_m,
        progress_callback=lambda p: progress_callback(0.5 * p),
        cancel_event=cancel_event,
    )

    series: Dict[str, Dict[str, List[float]]] = {}
    total = max(1, len(bnd_loc_geo))
    caches = bce._prefetch_source_caches(sources, bnd_loc_geo, point_sources, logger, cancel_event)
    for idx, (_, row) in enumerate(bnd_loc_geo.iterrows()):
        bce._check_cancel(cancel_event)
        name = str(row["boundaries"])
        src_idx = point_sources[row["boundaries"]]
        src = sources[src_idx]
        values = bce._resample_scalar_variable(
            row, src["ds"], src["preset"].water_level_var, src["preset"],
            src["native_minutes"], output_minutes, logger=logger, cache=caches[src_idx],
        )
        series.setdefault(name[:-2], {})[name[-1]] = (
            np.asarray(values, dtype=float).tolist()
        )
        progress_callback(0.5 + 0.5 * (idx + 1) / total)
    return series, qc_summary


def _write_qc_plot(path, supports, gradient_by_key, output_minutes, sources, span_along, span_cross) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_map, ax_ts) = plt.subplots(
        1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 1.6]}
    )
    for src in sources:
        preset = src["preset"]
        lats = np.asarray(src["ds"][preset.lat_coord].values, dtype=float)
        lons = np.asarray(src["ds"][preset.lon_coord].values, dtype=float)
        ax_map.add_patch(
            plt.Rectangle(
                (lons.min(), lats.min()), lons.max() - lons.min(), lats.max() - lats.min(),
                fill=False, linewidth=1.5, label=src["label"],
            )
        )
    hours = np.asarray(output_minutes, dtype=float) / 60.0
    for key, sup in supports.items():
        ax_map.plot(sup["lon"], sup["lat"], "o", ms=5)
        label = sup["names"][0]
        if sup.get("slid_m"):
            label += f" (shifted {sup['slid_m'] / 1000.0:.1f} km inboard)"
        ax_map.annotate(label, (sup["lon"], sup["lat"]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")
        ax_ts.plot(hours, gradient_by_key[key], linewidth=0.9, label="/".join(sup["names"]))
    ax_map.set_xlabel("lon"); ax_map.set_ylabel("lat")
    ax_map.set_title(f"Support points (probe span {span_along:.0f}x{span_cross:.0f} m)")
    ax_map.legend(loc="best", fontsize=7)
    ax_ts.set_xlabel("hours since reference date")
    ax_ts.set_ylabel("g = dWL/dl alongshore [-]")
    ax_ts.set_title("Neumann gradient time series")
    ax_ts.legend(loc="best", fontsize=7, ncol=2)
    ax_ts.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return str(path)


# ─────────────────────────────────────────────────────────────────────────────
# Self-tests
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=== backend_neumann.py self-test ===\n")
    all_pass = True

    def check(name, result, expected):
        ok = result == expected
        print(f"  {name}: [{'PASS' if ok else 'FAIL'}]")
        if not ok:
            globals()["all_pass"] = False
            print(f"    got {result!r}, expected {expected!r}")

    def check_close(name, result, expected, tol=1e-9):
        ok = np.allclose(result, expected, atol=tol)
        print(f"  {name}: [{'PASS' if ok else 'FAIL'}]")
        if not ok:
            globals()["all_pass"] = False
            print(f"    got {result!r}, expected {expected!r}")

    # Plane fit recovers a known linear field exactly
    e_along = (1.0, 0.0)
    offsets = _probe_offsets("plane", e_along, 2000.0, 1000.0)
    check("plane cloud has 9 probes", offsets.shape, (9, 2))
    b_true = np.array([1.0e-5, -2.0e-5, 3.0e-6])       # gradient east
    c_true = np.array([4.0e-6, 0.0, -1.0e-6])          # gradient north
    a_true = np.array([0.5, -0.2, 1.0])
    series = (
        a_true[None, :]
        + offsets[:, 0:1] * b_true[None, :]
        + offsets[:, 1:2] * c_true[None, :]
    )
    check_close("plane fit recovers alongshore gradient",
                _fit_gradient_series(offsets, series, e_along), b_true)
    check_close("plane fit projects onto rotated alongshore",
                _fit_gradient_series(offsets, series, (0.0, 1.0)), c_true)

    # Centered finite difference equals the plane fit on a linear field
    offsets_c = _probe_offsets("centered", e_along, 2000.0, 0.0)
    check("centered cloud has 2 probes", offsets_c.shape, (2, 2))
    series_c = a_true[None, :] + offsets_c[:, 0:1] * b_true[None, :]
    check_close("centered FD recovers gradient",
                _fit_gradient_series(offsets_c, series_c, e_along), b_true)

    # Collinear fallback: 1-D alongshore fit
    line_offsets = np.array([[-1000.0, 0.0], [0.0, 0.0], [1000.0, 0.0]])
    line_series = a_true[None, :] + line_offsets[:, 0:1] * b_true[None, :]
    check_close("collinear cloud falls back to 1-D fit",
                _fit_gradient_series(line_offsets, line_series, e_along), b_true)

    # Sign flip via direction vector
    check_close("flipped direction negates gradient",
                _fit_gradient_series(offsets, series, (-1.0, 0.0)), -b_true)

    # Along-section fallback: the cloud only resolves >= 1500 m from the coastal end.
    def fake_evaluate(easting, northing):
        n_ok = 9 if easting >= 1500.0 else 1
        return [{"offset": (0.0, 0.0), "source": 0, "lat": 0.0, "lon": 0.0}] * n_ok

    found = _find_usable_support_position(
        0.0, 0.0, [("North_1", "North_1_a", 10000.0, 0.0)], fake_evaluate, 3
    )
    check("slide finds the smallest resolvable displacement", found["distance_m"], 1562.5)
    check("slide reports the target end", found["target_name"], "North_1_a")
    check("slide returns the probes found there", len(found["usable"]), 9)
    check("slide respects the repair cap",
          _find_usable_support_position(0.0, 0.0, [("North_1", "North_1_a", 10000.0, 0.0)], fake_evaluate, 3, 1000.0),
          None)
    check("slide gives up when nothing resolves",
          _find_usable_support_position(0.0, 0.0, [("S", "S_b", 1000.0, 0.0)], fake_evaluate, 3), None)
    closer = _find_usable_support_position(
        0.0, 0.0,
        [("North_1", "North_1_a", 10000.0, 0.0), ("North_2", "North_2_b", 4000.0, 0.0)],
        fake_evaluate, 3,
    )
    # North_2 resolves at 0.375 x 4000 = 1500 m, North_1 at 0.15625 x 10000 = 1562.5 m.
    check("slide prefers the nearer resolvable position across owners",
          (closer["section"], closer["distance_m"]), ("North_2", 1500.0))

    # Alongshore direction from a synthetic north-up grid: constant-m section
    fake_grd = {
        "m": 3, "n": 3,
        "x_values": [[0.0, 0.0, 100.0, 200.0]] * 3,
        "y_values": [
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 100.0, 100.0],
            [0.0, 200.0, 200.0, 200.0],
        ],
    }
    fake_grd["x_values"] = [[0.0, 0.0, 100.0, 200.0],
                            [0.0, 0.0, 100.0, 200.0],
                            [0.0, 0.0, 100.0, 200.0]]
    fake_grd["y_values"] = [[0.0, 0.0, 0.0, 0.0],
                            [0.0, 100.0, 100.0, 100.0],
                            [0.0, 200.0, 200.0, 200.0]]
    sec_const_m = {"bnd_name": "West", "m_a": 1, "n_a": 1, "m_b": 1, "n_b": 3}
    check_close("constant-m section -> alongshore = +m (east)",
                _section_alongshore_direction(fake_grd, sec_const_m, False), (1.0, 0.0))
    sec_const_n = {"bnd_name": "North", "m_a": 1, "n_a": 3, "m_b": 3, "n_b": 3}
    check_close("constant-n section -> alongshore = +n (north)",
                _section_alongshore_direction(fake_grd, sec_const_n, False), (0.0, 1.0))
    check_close("flip_sign flips the direction",
                _section_alongshore_direction(fake_grd, sec_const_m, True), (-1.0, 0.0))
    try:
        _section_alongshore_direction(
            fake_grd, {"bnd_name": "Diag", "m_a": 1, "n_a": 1, "m_b": 3, "n_b": 3}, False
        )
        check("diagonal section rejected", "no error", "ValueError")
    except ValueError:
        check("diagonal section rejected", "ValueError", "ValueError")

    # Writer: N-type .bct structure
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        bct_path = str(Path(tmp) / "test_neumann.bct")
        writer_delft3d4.write_bct_neumann(
            bct_path,
            {"WestLat": {"g_a": [1.0e-5, 2.0e-5, 3.0e-5], "g_b": [0.0, -1.0e-5, -2.0e-5]}},
            "20180101", 0.0, 40.0, 20.0,
        )
        raw = open(bct_path, "rb").read()
        text = raw.decode("ascii")
        check("neumann parameter end A present", "'neumann         (n)  end A'" in text, True)
        check("neumann unit is dimensionless", "unit '[-]'" in text, True)
        check("records-in-table counted", "records-in-table     3" in text, True)
        check("CRLF line endings", raw.count(b"\r\n") > 0 and raw.count(b"\n") == raw.count(b"\r\n"), True)
        check("no reflection/water-elevation label", "water elevation" in text, False)

    print(f"\n{'=' * 40}")
    print("ALL TESTS PASSED" if all_pass else "SOME TESTS FAILED")
    print(f"{'=' * 40}")
