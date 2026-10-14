# Boundary-condition extraction – inputs

The Delft3D model files are included here. The EasyGSH-DB source data is not (3.5 GB and 26 GB); download it and point the module at it.

## Included

| Folder | Files | Used for |
|--------|-------|----------|
| `flow/` | `flow.grd`, `flow.bnd`, `flow.mdf` | Delft3D-FLOW grid, open boundaries and master file. The `.mdf` defines the simulation period 2013-10-20 to 2013-12-25 with reference date 2013-01-01 (Cartesian, UTM 32N). |
| `wave/` | `wave_outer.grd`, `wave_outer.bnd`, `wave.mdw` | Delft3D-WAVE outer grid, its boundaries and the `.mdw` template whose boundary section the module fills in. The template is the example model's `wave.mdw` (reference date 2013-01-01, `TSeriesFile = wave.bcw`). |

## Download separately

EasyGSH-DB, synoptic hydrodynamics, 1000 m grid, 20 min, year 2013, from <https://mdi-de.baw.de/easygsh/>:

| File | Size | Module setting |
|------|------|----------------|
| `2013_1000m_waterlevel_2D.nc` | 3.5 GB | Data files, extraction type *Water Level* |
| `2013_1000m_wave_2D.nc` | 26 GB | Data files, extraction type *Waves* |

The files can be anywhere on disk; the module only needs their paths.

## Settings used for the expected outputs

| Setting | Water level | Waves |
|---------|-------------|-------|
| Model type | Delft3D-4 | Delft3D-4 |
| Data source | EasyGSH | EasyGSH |
| `.mdf` file | `flow/flow.mdf` | `flow/flow.mdf` |
| `.grd` / `.bnd` | `flow/flow.grd`, `flow/flow.bnd` | `wave/wave_outer.grd`, `wave/wave_outer.bnd` |
| `.mdw` template | – | `wave/wave.mdw` |
| Wave timing source | – | Auto from .mdf |
| Allow intertidal boundaries | off | off |
| UTM zone | 32 N | 32 N |
| Output timestep | 20 min | 20 min |

The results are in `../output/water_level/` and `../output/waves/`.
