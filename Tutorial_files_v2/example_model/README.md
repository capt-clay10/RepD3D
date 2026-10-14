# Example model (Delft3D-4 FLOW + WAVE, German Bight, 2013-10-20 to 2013-12-25)

| File(s) | Origin |
|---------|--------|
| `flow.grd`, `flow.enc`, `flow.dep`, `flow.bnd`, `flow.mdf`, `flow.obs`, `flow.thd`, `flow.fou` | FLOW set-up of the v1.0 tutorial, unchanged |
| `wave.grd`, `wave.enc`, `wave.dep`, `wave_outer.grd`, `wave_outer.enc`, `wave_outer.dep`, `wave_outer.bnd`, `wave.loc`, `wave.mdw` | WAVE set-up of the v1.0 tutorial, unchanged. `wave.mdw` already contains the boundary blocks and `TSeriesFile = wave.bcw`. |
| `flow.bct` | Generated with RepD3D v2.0.0 (`../02_boundary_condition_extraction/output/water_level/`) from EasyGSH-DB 2013 |
| `wave.bcw` | Generated with RepD3D v2.0.0 (`../02_boundary_condition_extraction/output/waves/`) from EasyGSH-DB 2013 |

## Before running

Copy the meteo files from `../03_meteo_file_generator/output/` into this folder. `flow.mdf` references `cosmo_2013.amu`, `cosmo_2013.amv` and `cosmo_2013.amp`; `wave.mdw` references `xwind_cosmo_2013.wnd` and `ywind_cosmo_2013.wnd`. The `.mdw` also lists `xwind2_cosmo_2013.wnd` and `ywind2_cosmo_2013.wnd`; if your Delft3D-WAVE version insists on every listed file, provide copies under those names.

The v2.0.0 version of the wave master file written by the extraction module is `../02_boundary_condition_extraction/output/waves/wave_updated.mdw`. It is the same as `wave.mdw` here except that the constant `WaveHeight`, `Period`, `Direction` and `DirSpreading` lines are left out of each boundary block, because the time series in `wave.bcw` supply those values.
