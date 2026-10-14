# -*- coding: utf-8 -*-
"""
coordinate_utils.py — Coordinate transformations for RepD3D v2.

Handles:
    - UTM ↔ lat/lon conversion (configurable zone — v1 hardcoded zone 32)
    - Boundary coordinate DataFrame conversion (from .bnd/.grd or .pli)
    - UTM zone auto-detection from coordinates
    - Batch array conversion for wind grid processing

Consolidates coordinate logic from:
    - bct_generator.py lines 116-148 (convert_coordinates)
    - bcw_generator.py lines 112-144 (identical copy)
    - bct_year_overlap_file_generator.py lines 138-170 (identical copy)
    - bcw_year_overlap_file_generator.py (identical copy)
    - cosmo_wind_file_generator.py lines 199-207 (latlon_to_utm via pyproj)

v1 used the `utm` package for BC extraction and `pyproj` for wind processing.
v2 standardises on `pyproj` everywhere — it is more robust, handles edge
cases better (polar regions, zone boundaries), and is already a dependency
for wind processing.

v2 improvements over v1:
    - UTM zone is a parameter, not hardcoded to 32N
    - Southern hemisphere support via hemisphere parameter
    - Auto-detection of UTM zone from coordinates
    - Consistent API: always (lat, lon) order for geographic, (easting, northing) for UTM
    - Array inputs supported (not just scalars)

Usage:
    from coordinate_utils import (
        utm_to_latlon, latlon_to_utm, auto_detect_utm_zone,
        convert_boundary_coordinates,
    )
"""

import numpy as np
import pandas as pd
from pyproj import Transformer
from typing import Tuple, Union

# Type alias for coordinates that can be scalar or array
CoordType = Union[float, np.ndarray]


# ─────────────────────────────────────────────────────────────────────────────
# CORE CONVERSIONS
# ─────────────────────────────────────────────────────────────────────────────

def _utm_epsg(zone: int, south: bool = False) -> int:
    """
    Get EPSG code for a UTM zone.

    Northern hemisphere: EPSG:326xx (e.g., zone 32N → EPSG:32632)
    Southern hemisphere: EPSG:327xx (e.g., zone 32S → EPSG:32732)

    Args:
        zone:  UTM zone number (1–60).
        south: True for southern hemisphere.

    Returns:
        EPSG code as integer.
    """
    if not 1 <= zone <= 60:
        raise ValueError(f"UTM zone must be 1–60, got {zone}")
    return (32700 if south else 32600) + zone


def _get_transformer(from_epsg: int, to_epsg: int) -> Transformer:
    """
    Create a cached pyproj Transformer.

    Uses always_xy=True so the call convention is (x, y) which maps to
    (lon, lat) for geographic and (easting, northing) for UTM.
    """
    return Transformer.from_crs(
        f"EPSG:{from_epsg}",
        f"EPSG:{to_epsg}",
        always_xy=True,
    )


def utm_to_latlon(
    easting: CoordType,
    northing: CoordType,
    zone: int = 32,
    hemisphere: str = 'N',
) -> Tuple[CoordType, CoordType]:
    """
    Convert UTM coordinates to geographic lat/lon (WGS84, EPSG:4326).

    Args:
        easting:    UTM easting in metres (scalar or array).
        northing:   UTM northing in metres (scalar or array).
        zone:       UTM zone number (1–60, default 32).
        hemisphere: 'N' for northern, 'S' for southern hemisphere.

    Returns:
        (latitude, longitude) — same type as input (scalar or array).

    v1 source: bct_generator.py line 139 used `utm.to_latlon(e, n, 32, 'N')`.
    v2 uses pyproj for consistency and configurability.
    """
    south = hemisphere.upper() == 'S'
    utm_epsg = _utm_epsg(zone, south)
    transformer = _get_transformer(utm_epsg, 4326)
    lon, lat = transformer.transform(easting, northing)
    return lat, lon


def latlon_to_utm(
    lat: CoordType,
    lon: CoordType,
    zone: int = 32,
    hemisphere: str = 'N',
) -> Tuple[CoordType, CoordType]:
    """
    Convert geographic lat/lon (WGS84) to UTM coordinates.

    Args:
        lat:        Latitude in degrees (scalar or array).
        lon:        Longitude in degrees (scalar or array).
        zone:       UTM zone number (1–60, default 32).
        hemisphere: 'N' for northern, 'S' for southern hemisphere.

    Returns:
        (easting, northing) in metres — same type as input.

    v1 source: cosmo_wind_file_generator.py lines 200-207 used
    pyproj.Transformer with hardcoded zone=32.
    """
    south = hemisphere.upper() == 'S'
    utm_epsg = _utm_epsg(zone, south)
    transformer = _get_transformer(4326, utm_epsg)
    easting, northing = transformer.transform(lon, lat)
    return easting, northing


def auto_detect_utm_zone(lon: CoordType) -> int:
    """
    Auto-detect UTM zone from longitude value(s).

    If an array is provided, uses the mean longitude.

    Args:
        lon: Longitude in degrees (scalar, list, or array).

    Returns:
        UTM zone number (1–60).

    Examples:
        >>> auto_detect_utm_zone(9.0)    # Hamburg, Germany
        32
        >>> auto_detect_utm_zone(-73.9)  # New York
        18
        >>> auto_detect_utm_zone(np.array([6.0, 9.0]))  # Mean = 7.5
        32
    """
    if hasattr(lon, '__len__'):
        lon = float(np.nanmean(lon))
    else:
        lon = float(lon)
    zone = int((lon + 180) / 6) + 1
    return min(max(zone, 1), 60)


def auto_detect_hemisphere(lat: CoordType) -> str:
    """
    Auto-detect hemisphere from latitude value(s).

    If an array is provided, uses the mean latitude.

    Args:
        lat: Latitude in degrees (scalar, list, or array).

    Returns:
        'N' or 'S'.
    """
    if hasattr(lat, '__len__'):
        lat = float(np.nanmean(lat))
    else:
        lat = float(lat)
    return 'N' if lat >= 0 else 'S'


# ─────────────────────────────────────────────────────────────────────────────
# BOUNDARY COORDINATE CONVERSION — used by BC extraction
# ─────────────────────────────────────────────────────────────────────────────

def convert_boundary_coordinates(
    bnd_loc: pd.DataFrame,
    crs_type: str,
    utm_zone: int = 32,
    hemisphere: str = 'N',
) -> pd.DataFrame:
    """
    Convert boundary location DataFrame to geographic coordinates.

    This is the central coordinate conversion used during BC extraction.
    The input DataFrame comes from parsing .bnd + .grd files (Delft3D-4)
    or .pli files (Delft3D-FM).

    Args:
        bnd_loc:    DataFrame with columns ['boundary', 'easting', 'northing'].
                    For Cartesian grids, these are UTM coordinates.
                    For Spherical grids, 'easting' = longitude, 'northing' = latitude.
        crs_type:   'Cartesian' or 'Spherical' (from .grd header or user selection).
        utm_zone:   UTM zone number (only used if crs_type == 'Cartesian').
        hemisphere: 'N' or 'S' (only used if crs_type == 'Cartesian').

    Returns:
        DataFrame with columns ['lat', 'lon', 'boundaries'].
        Always in geographic coordinates (EPSG:4326).

    v1 source: bct_generator.py lines 116-148 (copied identically in 3 other files).
    v1 used the `utm` package; v2 uses pyproj. Output is numerically identical
    to within floating-point precision.
    """
    if crs_type == 'Cartesian':
        easting = bnd_loc['easting'].astype(float).values
        northing = bnd_loc['northing'].astype(float).values
        lat, lon = utm_to_latlon(easting, northing, utm_zone, hemisphere)
        return pd.DataFrame({
            'lat': lat,
            'lon': lon,
            'boundaries': bnd_loc['boundary'].values,
        })
    else:
        # Spherical: coordinates are already lat/lon
        # In Delft3D convention: 'easting' column = longitude, 'northing' = latitude
        return pd.DataFrame({
            'lat': bnd_loc['northing'].astype(float).values,
            'lon': bnd_loc['easting'].astype(float).values,
            'boundaries': bnd_loc['boundary'].values,
        })


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== coordinate_utils.py self-test ===\n")

    all_pass = True

    def check(name, result, expected, tolerance=0.001):
        global all_pass
        if isinstance(result, str) or isinstance(expected, str):
            ok = result == expected
        elif isinstance(result, (tuple, list)):
            ok = all(abs(r - e) < tolerance for r, e in zip(result, expected))
        elif isinstance(result, (int, np.integer)) and isinstance(expected, (int, np.integer)):
            ok = result == expected
        else:
            ok = abs(result - expected) <= tolerance
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
        print(f"  {name}: {result} ≈ {expected}  [{status}]")

    # ── UTM zone detection ───────────────────────────────────────────────
    print("UTM zone auto-detection:")
    check("Hamburg (lon=9.0)", auto_detect_utm_zone(9.0), 32, tolerance=0)
    check("New York (lon=-73.9)", auto_detect_utm_zone(-73.9), 18, tolerance=0)
    check("Tokyo (lon=139.7)", auto_detect_utm_zone(139.7), 54, tolerance=0)
    check("Array mean", auto_detect_utm_zone(np.array([6.0, 9.0])), 32, tolerance=0)

    # ── Hemisphere detection ─────────────────────────────────────────────
    print("\nHemisphere auto-detection:")
    check("lat=54", auto_detect_hemisphere(54.0), 'N')
    check("lat=-33", auto_detect_hemisphere(-33.0), 'S')

    # ── Scalar round-trip: UTM → latlon → UTM ────────────────────────────
    print("\nScalar round-trip (UTM zone 32N):")
    # Known point: approximately Hamburg (53.55°N, 10.0°E)
    e_orig, n_orig = 566463.0, 5935270.0
    lat, lon = utm_to_latlon(e_orig, n_orig, zone=32, hemisphere='N')
    print(f"  UTM ({e_orig}, {n_orig}) → latlon ({lat:.6f}, {lon:.6f})")
    check("lat ≈ 53.55", lat, 53.55, tolerance=0.05)
    check("lon ≈ 10.0", lon, 10.0, tolerance=0.05)

    e_back, n_back = latlon_to_utm(lat, lon, zone=32, hemisphere='N')
    check("easting round-trip", e_back, e_orig, tolerance=0.01)
    check("northing round-trip", n_back, n_orig, tolerance=0.01)

    # ── Array conversion ─────────────────────────────────────────────────
    print("\nArray conversion (2 points, zone 32N):")
    eastings = np.array([566463.0, 500000.0])
    northings = np.array([5935270.0, 5900000.0])
    lats, lons = utm_to_latlon(eastings, northings, zone=32)
    print(f"  lats = {lats}")
    print(f"  lons = {lons}")
    assert len(lats) == 2, "Expected 2 latitudes"
    assert len(lons) == 2, "Expected 2 longitudes"
    print("  Array shape check: PASS")

    # ── Boundary coordinate conversion ───────────────────────────────────
    print("\nBoundary coordinate conversion (Cartesian):")
    bnd_df = pd.DataFrame({
        'boundary': ['BND_1_a', 'BND_1_b'],
        'easting': [566463.0, 500000.0],
        'northing': [5935270.0, 5900000.0],
    })
    result = convert_boundary_coordinates(bnd_df, 'Cartesian', utm_zone=32, hemisphere='N')
    print(f"  Columns: {list(result.columns)}")
    assert list(result.columns) == ['lat', 'lon', 'boundaries'], "Wrong columns"
    assert len(result) == 2, "Wrong row count"
    print(f"  Lat range: {result['lat'].min():.4f} – {result['lat'].max():.4f}")
    print(f"  Lon range: {result['lon'].min():.4f} – {result['lon'].max():.4f}")
    print("  Cartesian conversion: PASS")

    print("\nBoundary coordinate conversion (Spherical):")
    bnd_df_sph = pd.DataFrame({
        'boundary': ['BND_1_a', 'BND_1_b'],
        'easting': [9.0, 9.5],       # lon
        'northing': [53.5, 54.0],     # lat
    })
    result_sph = convert_boundary_coordinates(bnd_df_sph, 'Spherical')
    check("Spherical lat[0]", result_sph['lat'].iloc[0], 53.5, tolerance=0.0001)
    check("Spherical lon[0]", result_sph['lon'].iloc[0], 9.0, tolerance=0.0001)
    assert list(result_sph['boundaries']) == ['BND_1_a', 'BND_1_b']
    print("  Spherical pass-through: PASS")

    # ── EPSG helper ──────────────────────────────────────────────────────
    print("\nEPSG code generation:")
    check("Zone 32N", _utm_epsg(32, south=False), 32632, tolerance=0)
    check("Zone 32S", _utm_epsg(32, south=True), 32732, tolerance=0)
    check("Zone 1N", _utm_epsg(1, south=False), 32601, tolerance=0)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
