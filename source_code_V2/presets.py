# -*- coding: utf-8 -*-
"""
presets.py — Data source configuration presets for RepD3D v2.

Each preset encapsulates the variable names, coordinate names, dimension
order, time offsets, and other metadata for a specific oceanographic or
atmospheric data source. This makes the backend extraction code fully
data-source-agnostic — it reads the preset to know which variable to
look up, which coordinate names to use, etc.

Preset Categories:
    1. BC Presets (BCPreset): For boundary condition data sources
       - EasyGSH, TrilaWatt, Copernicus SSH, Copernicus Waves, Copernicus BGC, Custom
    2. Wind Presets (WindPreset): For atmospheric reanalysis data sources
       - COSMO REA6, ICON-DREAM-EU, ICON-DREAM-Global, ERA5

All the quirks and gotchas for each data source are documented in the
preset fields. This is the single source of truth for data source
configuration — if a variable name changes or a new data source is
added, only this file needs updating.

Usage:
    from presets import (
        EASYGSH, TRILAWATT, COPERNICUS_SSH, COPERNICUS_WAVES,
        COSMO_REA6, ICON_DREAM_EU, ERA5,
        get_bc_preset, get_wind_preset,
        BC_PRESET_NAMES, WIND_PRESET_NAMES,
    )
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict


# ─────────────────────────────────────────────────────────────────────────────
# BC PRESETS — Boundary condition data sources
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BCPreset:
    """Configuration for a boundary-condition data source."""

    # ── Identity ─────────────────────────────────────────────────────────
    name: str
    description: str

    # ── Coordinate names in the NetCDF ───────────────────────────────────
    lat_coord: str = 'lat'
    lon_coord: str = 'lon'
    time_coord: str = 'time'

    # ── Dimension order (for documentation and validation) ───────────────
    # Standard CF convention: (time, lat, lon)
    # EasyGSH exception: (lon, lat, time)
    dim_order: str = '(time, lat, lon)'

    # ── Water level / current variables ──────────────────────────────────
    water_level_var: Optional[str] = None
    current_u_var: Optional[str] = None
    current_v_var: Optional[str] = None

    # ── Wave variables ───────────────────────────────────────────────────
    wave_hs_var: Optional[str] = None       # Significant wave height
    wave_tp_var: Optional[str] = None       # Peak period
    wave_tm_var: Optional[str] = None       # Mean period (optional)
    wave_dir_var: Optional[str] = None      # Direction in degrees (Copernicus)
    wave_dir_x_var: Optional[str] = None    # Direction x-component (EasyGSH/TrilaWatt)
    wave_dir_y_var: Optional[str] = None    # Direction y-component (EasyGSH/TrilaWatt)
    wave_dir_spread_var: Optional[str] = None  # Directional spreading

    # ── Wave direction format ────────────────────────────────────────────
    # 'vector':  x/y components → requires atan2 conversion (EasyGSH, TrilaWatt)
    # 'degrees': direct nautical degrees (Copernicus)
    wave_dir_format: str = 'vector'

    # ── Wave directional spreading fallback ──────────────────────────────
    # If wave_dir_spread_var is None (e.g. Copernicus), the backend writes
    # this constant value for every timestep. Spec §2.3.
    wave_dir_spread_default: float = 30.0
    wave_dir_spread_use_default: bool = False

    # ── Sediment variables ───────────────────────────────────────────────
    # One entry per sediment fraction declared in the .sed file. If a single
    # NetCDF variable (e.g. SPM) feeds multiple fractions, repeat the same
    # variable name and use sediment_split_ratios to partition it.
    sediment_conc_vars: List[str] = field(default_factory=list)
    sediment_split_ratios: List[float] = field(default_factory=list)
    sediment_fraction_names: List[str] = field(default_factory=list)

    # ── Salinity variable ───────────────────────────────────────────────
    # Written to .bcc with 'Salinity end A/B uniform' parameter labels.
    # Spec §3.5. Typical names: 'so', 'salinity', 'vosaline', 'salt'.
    salinity_var: Optional[str] = None

    # ── Time handling ────────────────────────────────────────────────────
    # Offset applied to data timestamps (hours).
    # EasyGSH = -12 (data starts 12h before stated coverage_start)
    # All others = 0
    time_offset_hours: float = 0.0

    # Default timestep in the source data (minutes)
    default_timestep_minutes: float = 20.0

    # ── File organisation ────────────────────────────────────────────────
    # 'yearly':  one file per year (EasyGSH, TrilaWatt)
    # 'daily':   one file per day (Copernicus)
    # 'single':  user provides file(s), no assumption
    file_pattern: str = 'single'

    # ── CRS ──────────────────────────────────────────────────────────────
    native_crs: str = 'EPSG:4326'

    # ── Fill value ───────────────────────────────────────────────────────
    # Value used for land/missing cells in the source data
    fill_value: Optional[float] = None  # None = use NetCDF _FillValue


# ── EasyGSH (German Bight, BAW) ─────────────────────────────────────────
EASYGSH = BCPreset(
    name='EasyGSH',
    description='German Bight hindcast (BAW), ~1km, EPSG:4326',
    lat_coord='lat',
    lon_coord='lon',
    time_coord='nMesh2_data_time',       # Non-standard name
    dim_order='(lon, lat, time)',          # UNUSUAL — not CF standard
    water_level_var='Mesh2_face_Wasserstand_2d',
    wave_hs_var='Mesh2_face_signifikante_Wellenhoehe_2d',
    wave_tp_var='Mesh2_face_Peak_Wellenperiode_2d',
    wave_tm_var='Mesh2_face_mittlere_Wellenperiode_2d',
    wave_dir_x_var='Mesh2_face_Wellenrichtungsvektor_x_2d',
    wave_dir_y_var='Mesh2_face_Wellenrichtungsvektor_y_2d',
    wave_dir_spread_var='Mesh2_face_Richtungsaufweitung_der_Wellen_2d',
    wave_dir_format='vector',
    time_offset_hours=-12.0,  # CRITICAL — see handoff doc Section 5.1
    default_timestep_minutes=20.0,
    file_pattern='yearly',
    # Typical filename: 2015_1000m_waterlevel_2D.nc, 2015_1000m_wave_2D.nc
)

# ── TrilaWatt (Trilateral Wadden Sea, BAW) ──────────────────────────────
TRILAWATT = BCPreset(
    name='TrilaWatt',
    description='Trilateral Wadden Sea (BAW), ~500m, EPSG:4326',
    lat_coord='lat',
    lon_coord='lon',
    time_coord='time',
    dim_order='(time, lat, lon)',          # Standard CF
    water_level_var='sea_surface_height_2d',
    current_u_var='current_velocity_2d_x',
    current_v_var='current_velocity_2d_y',
    wave_hs_var='significant_wave_height',
    wave_tp_var='wave_peak_period',
    wave_tm_var='wave_mean_period',
    wave_dir_x_var='wave_to_direction_x',
    wave_dir_y_var='wave_to_direction_y',
    wave_dir_spread_var='wave_directional_spread',
    wave_dir_format='vector',
    time_offset_hours=0.0,  # No offset (unlike EasyGSH)
    default_timestep_minutes=20.0,
    file_pattern='yearly',
    fill_value=1e+31,
    # Regional files: tides_2017_ger.nc, tides_2017_dk.nc, tides_2017_nl.nc
    # Waves: waves_2017_ger.nc, etc.
)

# ── Copernicus Marine — SSH & Currents ──────────────────────────────────
COPERNICUS_SSH = BCPreset(
    name='Copernicus SSH/Currents',
    description='CMEMS NW Shelf (NWSHELF_ANALYSIS_FORECAST_PHY_004_013), ~1.5km',
    lat_coord='lat',
    lon_coord='lon',
    time_coord='time',
    dim_order='(time, lat, lon)',
    water_level_var='zos',
    current_u_var='uo',
    current_v_var='vo',
    time_offset_hours=0.0,
    default_timestep_minutes=15.0,
    file_pattern='daily',
    # Daily files: metoffice_foam1_amm15_NWS_SSH_*.nc (~50 MB/day)
    # Currents:    metoffice_foam1_amm15_NWS_SSC_*.nc (~151 MB/day)
)

# ── Copernicus Marine — Waves ───────────────────────────────────────────
COPERNICUS_WAVES = BCPreset(
    name='Copernicus Waves',
    description='CMEMS NW Shelf waves, ~1.5km, 3-hourly',
    lat_coord='latitude',     # NOTE: Inconsistent with SSH!
    lon_coord='longitude',    # NOTE: Inconsistent with SSH!
    time_coord='time',
    dim_order='(time, lat, lon)',
    wave_hs_var='VHM0',
    wave_tp_var='VTPK',
    wave_tm_var='VTM02',
    wave_dir_var='VMDR',      # Direction in DEGREES (not vectors!)
    wave_dir_format='degrees',
    wave_dir_spread_var=None,  # Not directly available in Copernicus
    wave_dir_spread_use_default=True,
    wave_dir_spread_default=30.0,
    time_offset_hours=0.0,
    default_timestep_minutes=180.0,  # 3 hours
    file_pattern='daily',
    # Files: metoffice_wave_amm15_NWS_WAV_*.nc (~59 MB/day)
)

# ── Copernicus Marine — BGC (Biogeochemistry, sediments) ───────────────
COPERNICUS_BGC = BCPreset(
    name='Copernicus BGC',
    description='CMEMS NW Shelf biogeochemistry — SPM only (single fraction)',
    lat_coord='latitude',
    lon_coord='longitude',
    time_coord='time',
    dim_order='(time, lat, lon)',
    sediment_conc_vars=['SPM'],
    sediment_fraction_names=['spm'],
    sediment_split_ratios=[1.0],
    time_offset_hours=0.0,
    default_timestep_minutes=1440.0,  # daily
    file_pattern='daily',
)

# ── Custom (user-defined) ──────────────────────────────────────────────
CUSTOM_BC = BCPreset(
    name='Custom',
    description='User-defined NetCDF dataset — fill in variable names in GUI',
)


# ─────────────────────────────────────────────────────────────────────────────
# WIND PRESETS — Atmospheric reanalysis data sources
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class WindPreset:
    """Configuration for an atmospheric reanalysis data source."""

    # ── Identity ─────────────────────────────────────────────────────────
    name: str
    description: str

    # ── Grid type ────────────────────────────────────────────────────────
    # 'rotated_ll':     Rotated lat/lon (COSMO)
    # 'unstructured':   Triangular mesh (ICON)
    # 'regular_ll':     Regular lat/lon (ERA5)
    grid_type: str

    # ── Variable names ───────────────────────────────────────────────────
    u_var: str = 'u10'
    v_var: str = 'v10'
    p_var: str = 'sp'

    # ── Temperature variable ─────────────────────────────────────────────
    # Name of the 2m air temperature variable (in Kelvin). Used for .amt
    # extraction. Typical values:
    #   - COSMO REA6:       't2m' (from T_2M files)
    #   - ERA5:             't2m'
    #   - ICON-DREAM-EU:    't2m' (from T_2M files)
    #   - ICON-DREAM-Global:'t'   (multi-level; lowest level is taken)
    t_var: Optional[str] = None

    # True when t_var is a multi-level field that needs level selection
    # (ICON-DREAM-Global: 'generalVerticalLayer' dim, lowest ≈ surface).
    t_multilevel: bool = False

    # ── Pressure unit ────────────────────────────────────────────────────
    # 'Pa' → divide by 100 for mbar (Delft3D expects mbar)
    p_unit: str = 'Pa'

    # ── Vector rotation ─────────────────────────────────────────────────
    # COSMO: MUST rotate despite GRIB flag saying otherwise (see handoff doc)
    # ICON, ERA5: no rotation needed
    needs_rotation: bool = False

    # ── Auxiliary file ───────────────────────────────────────────────────
    # COSMO: cosmo_latlon_rot.grb (true geographic coords + rotation info)
    # ICON:  ICON-DREAM-EU_grid.nc (cell centres, connectivity)
    # ERA5:  none needed
    needs_aux_file: bool = False

    # ── Coordinate handling ──────────────────────────────────────────────
    # ICON coordinates are in RADIANS (must convert to degrees)
    coord_unit: str = 'degrees'

    # ERA5 longitude: 0°–360° (must wrap to -180°–180° for European domains)
    lon_convention: str = 'geographic'  # 'geographic' (-180 to 180) or '0_360'

    # Source sampling CRS used during interpolation. Current atmospheric
    # presets all sample in geographic lat/lon, even for Cartesian target grids.
    sampling_crs: str = 'geographic'

    # ── Interpolation method ─────────────────────────────────────────────
    # 'griddata':    scipy.griddata — for structured 2D grids (COSMO, ERA5)
    # 'kdtree_idw':  KDTree + IDW — for unstructured grids (ICON)
    interpolation_method: str = 'griddata'

    # ── Time structure ───────────────────────────────────────────────────
    # 'flat':       standard time dimension (COSMO, ERA5)
    # 'time_step':  time(249) × step(3) needs flattening (ICON)
    time_structure: str = 'flat'

    # ── COSMO rotation pole (only used if needs_rotation=True) ───────────
    rotation_pole_lat: Optional[float] = None
    rotation_pole_lon: Optional[float] = None

    # ── File organisation ────────────────────────────────────────────────
    file_pattern: str = 'monthly'


# ── COSMO REA6 (DWD, ~6km, Europe) ─────────────────────────────────────
COSMO_REA6 = WindPreset(
    name='COSMO REA6',
    description='DWD Reanalysis, ~6km, Europe (rotated grid)',
    grid_type='rotated_ll',
    needs_rotation=True,
    needs_aux_file=True,
    interpolation_method='griddata',
    rotation_pole_lat=39.25,
    rotation_pole_lon=-162.0,
    t_var='t2m',
    # U files: UV/U_10M.2D.YYYYMM.grb
    # V files: UV/V_10M.2D.YYYYMM.grb
    # P files: PS/PS.2D.YYYYMM.grb
    # T files: T_2M.2D.YYYYMM.grb  (variable: t2m, unit: K)
    # Aux:     cosmo_latlon_rot.grb or COSMO_REA6_CONST_withOUTsponge.nc
)

# ── ICON-DREAM-EU (DWD, ~7km, Europe) ──────────────────────────────────
ICON_DREAM_EU = WindPreset(
    name='ICON-DREAM-EU',
    description='DWD ICON, ~7km, Europe (unstructured triangular)',
    grid_type='unstructured',
    needs_rotation=False,
    needs_aux_file=True,
    coord_unit='radians',               # Cell centres in RADIANS!
    interpolation_method='kdtree_idw',
    time_structure='time_step',          # time(249) × step(3) → flatten
    t_var='t2m',
    # U files: ICON-DREAM-EU_YYYYMM_U_10M_hourly.grb
    # V files: ICON-DREAM-EU_YYYYMM_V_10M_hourly.grb
    # P files: ICON-DREAM-EU_YYYYMM_PS_hourly.grb
    # T files: ICON-DREAM-EU_YYYYMM_T_2M_hourly.grb  (variable: t2m, unit: K)
    # Grid:    ICON-DREAM-EU_grid.nc
)

# ── ICON-DREAM-Global (DWD, ~13km, Global) ─────────────────────────────
ICON_DREAM_GLOBAL = WindPreset(
    name='ICON-DREAM-Global',
    description='DWD ICON Global, ~13km, Global (unstructured triangular)',
    grid_type='unstructured',
    needs_rotation=False,
    needs_aux_file=True,
    coord_unit='radians',
    interpolation_method='kdtree_idw',
    time_structure='time_step',
    p_var='sp',
    t_var='t',           # 3D 'Temperature' field on generalVerticalLayer
    t_multilevel=True,   # Pick lowest level (near-surface) during extraction
    # U files: ICON-DREAM-Global_YYYYMM_U_10M_hourly.grb
    # V files: ICON-DREAM-Global_YYYYMM_V_10M_hourly.grb
    # P files: ICON-DREAM-Global_YYYYMM_PS_hourly.grb  (variable: sp)
    # T files: ICON-DREAM-Global_YYYYMM_T_hourly.grb
    #          (variable: t, dims: time × step × generalVerticalLayer × values, K)
    # Grid:    ICON-DREAM-Global_grid.nc
)

# ── ERA5 (ECMWF, ~28km / 0.25°, Global) ────────────────────────────────
ERA5 = WindPreset(
    name='ERA5',
    description='ECMWF Reanalysis, ~28km (0.25°), Global',
    grid_type='regular_ll',
    needs_rotation=False,
    needs_aux_file=False,
    lon_convention='0_360',              # Must wrap for European domains
    interpolation_method='griddata',
    file_pattern='monthly',
    t_var='t2m',
    # Combined file: era5_YYYY_MM_DD.grib (~4 GB/month at global scale)
    # Temperature variable: t2m (2m air temperature, K)
)

# ─────────────────────────────────────────────────────────────────────────────
# PRESET LOOKUP
# ─────────────────────────────────────────────────────────────────────────────

_BC_PRESETS: Dict[str, BCPreset] = {
    'EasyGSH': EASYGSH,
    'TrilaWatt': TRILAWATT,
    'Copernicus SSH/Currents': COPERNICUS_SSH,
    'Copernicus Waves': COPERNICUS_WAVES,
    'Copernicus BGC': COPERNICUS_BGC,
    'Custom': CUSTOM_BC,
}

_WIND_PRESETS: Dict[str, WindPreset] = {
    'COSMO REA6': COSMO_REA6,
    'ICON-DREAM-EU': ICON_DREAM_EU,
    'ICON-DREAM-Global': ICON_DREAM_GLOBAL,
    'ERA5 (.nc)': ERA5,
}
_WIND_PRESET_ALIASES: Dict[str, str] = {
    'ERA5': 'ERA5 (.nc)',
}

# Lists for GUI dropdown menus
BC_PRESET_NAMES: List[str] = list(_BC_PRESETS.keys())
WIND_PRESET_NAMES: List[str] = list(_WIND_PRESETS.keys())


def get_bc_preset(name: str) -> BCPreset:
    """
    Look up a BC preset by name.

    Args:
        name: Preset name (must be in BC_PRESET_NAMES).

    Returns:
        BCPreset instance.

    Raises:
        KeyError if name not found.
    """
    if name not in _BC_PRESETS:
        raise KeyError(
            f"Unknown BC preset '{name}'. "
            f"Available: {BC_PRESET_NAMES}"
        )
    return _BC_PRESETS[name]


def get_wind_preset(name: str) -> WindPreset:
    """
    Look up a Wind preset by name.

    Args:
        name: Preset name (must be in WIND_PRESET_NAMES).

    Returns:
        WindPreset instance.

    Raises:
        KeyError if name not found.
    """
    preset_name = _WIND_PRESET_ALIASES.get(name, name)
    if preset_name not in _WIND_PRESETS:
        raise KeyError(
            f"Unknown Wind preset '{name}'. "
            f"Available: {WIND_PRESET_NAMES}"
        )
    return _WIND_PRESETS[preset_name]


# ─────────────────────────────────────────────────────────────────────────────
# VARIABLE NAME MAPPING SUMMARY (for reference / documentation)
# ─────────────────────────────────────────────────────────────────────────────
#
# | Quantity    | EasyGSH                                       | TrilaWatt                    | Copernicus        |
# |-------------|-----------------------------------------------|------------------------------|-------------------|
# | Water level | Mesh2_face_Wasserstand_2d                     | sea_surface_height_2d        | zos               |
# | Current U   | (not in WL file)                              | current_velocity_2d_x        | uo                |
# | Current V   | (not in WL file)                              | current_velocity_2d_y        | vo                |
# | Wave Hs     | Mesh2_face_signifikante_Wellenhoehe_2d        | significant_wave_height      | VHM0              |
# | Wave Tp     | Mesh2_face_Peak_Wellenperiode_2d              | wave_peak_period             | VTPK              |
# | Wave Dir    | x/y vectors (2 vars)                          | x/y vectors (2 vars)         | VMDR (degrees)    |
# | Dir Spread  | Mesh2_face_Richtungsaufweitung_der_Wellen_2d  | wave_directional_spread      | (not available)   |
# | Time dim    | nMesh2_data_time                              | time                         | time              |
# | Lat coord   | lat                                           | lat                          | lat / latitude    |
# | Lon coord   | lon                                           | lon                          | lon / longitude   |
# | Dim order   | (lon, lat, time)                              | (time, lat, lon)             | (time, lat, lon)  |
# | Time step   | 20 min                                        | 20 min                       | 15 min / 3h       |
# | Time offset | -12 hours                                     | None                         | None              |
# | Files/period| 1 per year                                    | 1 per year per region        | 1 per day         |
#
# | Wind Property     | COSMO REA6        | ICON-DREAM-EU           | ERA5              |
# |-------------------|-------------------|-------------------------|-------------------|
# | Grid type         | Rotated lat/lon   | Unstructured triangular | Regular lat/lon   |
# | Resolution        | ~6 km             | ~7 km                   | ~28 km (0.25°)    |
# | Rotation needed   | YES               | NO                      | NO                |
# | U variable        | u10               | u10                     | u10               |
# | V variable        | v10               | v10                     | v10               |
# | P variable        | sp (Pa)           | sp (Pa)                 | sp (Pa)           |
# | Coord units       | degrees (2D)      | RADIANS (1D)            | degrees (1D)      |
# | Interpolation     | griddata          | KDTree + IDW            | griddata          |
# | Lon convention    | geographic        | geographic              | 0°-360°           |
# | Time structure    | flat              | time × step             | flat              |


# ─────────────────────────────────────────────────────────────────────────────
# SELF-TEST
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=== presets.py self-test ===\n")
    all_pass = True

    def check(name, result, expected):
        global all_pass
        ok = result == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
            print(f"  {name}: [{status}]  got {result!r}, expected {expected!r}")
        else:
            print(f"  {name}: [{status}]")

    # ── BC preset lookups ────────────────────────────────────────────────
    print("BC preset lookups:")
    check("EasyGSH name", get_bc_preset('EasyGSH').name, 'EasyGSH')
    check("EasyGSH time_offset", get_bc_preset('EasyGSH').time_offset_hours, -12.0)
    check("EasyGSH time_coord", get_bc_preset('EasyGSH').time_coord, 'nMesh2_data_time')
    check("EasyGSH WL var", get_bc_preset('EasyGSH').water_level_var, 'Mesh2_face_Wasserstand_2d')
    check("EasyGSH dir format", get_bc_preset('EasyGSH').wave_dir_format, 'vector')

    check("TrilaWatt offset=0", get_bc_preset('TrilaWatt').time_offset_hours, 0.0)
    check("TrilaWatt WL var", get_bc_preset('TrilaWatt').water_level_var, 'sea_surface_height_2d')
    check("TrilaWatt has currents", get_bc_preset('TrilaWatt').current_u_var is not None, True)

    check("Copernicus SSH var", get_bc_preset('Copernicus SSH/Currents').water_level_var, 'zos')
    check("Copernicus SSH step", get_bc_preset('Copernicus SSH/Currents').default_timestep_minutes, 15.0)

    check("Copernicus Waves lat", get_bc_preset('Copernicus Waves').lat_coord, 'latitude')
    check("Copernicus Waves dir format", get_bc_preset('Copernicus Waves').wave_dir_format, 'degrees')
    check("Copernicus Waves no spread", get_bc_preset('Copernicus Waves').wave_dir_spread_var, None)

    # ── Wind preset lookups ──────────────────────────────────────────────
    print("\nWind preset lookups:")
    check("COSMO rotation", get_wind_preset('COSMO REA6').needs_rotation, True)
    check("COSMO pole lat", get_wind_preset('COSMO REA6').rotation_pole_lat, 39.25)
    check("COSMO pole lon", get_wind_preset('COSMO REA6').rotation_pole_lon, -162.0)
    check("COSMO interp", get_wind_preset('COSMO REA6').interpolation_method, 'griddata')

    check("ICON rotation", get_wind_preset('ICON-DREAM-EU').needs_rotation, False)
    check("ICON coords", get_wind_preset('ICON-DREAM-EU').coord_unit, 'radians')
    check("ICON interp", get_wind_preset('ICON-DREAM-EU').interpolation_method, 'kdtree_idw')
    check("ICON time", get_wind_preset('ICON-DREAM-EU').time_structure, 'time_step')

    check("ERA5 rotation", get_wind_preset('ERA5').needs_rotation, False)
    check("ERA5 lon", get_wind_preset('ERA5').lon_convention, '0_360')
    check("ERA5 no aux", get_wind_preset('ERA5').needs_aux_file, False)

    # ── Name lists for GUI ───────────────────────────────────────────────
    print("\nPreset name lists:")
    check("BC names count", len(BC_PRESET_NAMES), 6)
    check("Wind names count", len(WIND_PRESET_NAMES), 4)
    check("EasyGSH in BC list", 'EasyGSH' in BC_PRESET_NAMES, True)
    check("COSMO in Wind list", 'COSMO REA6' in WIND_PRESET_NAMES, True)
    check("Custom in BC list", 'Custom' in BC_PRESET_NAMES, True)

    # ── Error handling ───────────────────────────────────────────────────
    print("\nError handling:")
    try:
        get_bc_preset('NonExistent')
        check("bad BC preset raises", False, True)
    except KeyError:
        check("bad BC preset raises", True, True)

    try:
        get_wind_preset('NonExistent')
        check("bad Wind preset raises", False, True)
    except KeyError:
        check("bad Wind preset raises", True, True)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED — check output above")
    print(f"{'=' * 50}")
