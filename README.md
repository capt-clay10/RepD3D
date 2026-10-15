# RepD3D – Representative-period Detection & Delft3D Boundary-data Toolbox

[![GitHub release](https://img.shields.io/github/v/release/capt-clay10/RepD3D?style=for-the-badge)](https://github.com/capt-clay10/RepD3D/releases)
[![License](https://img.shields.io/github/license/capt-clay10/RepD3D?style=for-the-badge)](LICENSE)
![Python](https://img.shields.io/badge/python-3.14-blue?style=for-the-badge)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23361961-blue?style=for-the-badge)](https://doi.org/10.5281/zenodo.23361961)


![RepD3D logo](https://github.com/user-attachments/assets/b1de47b1-be95-4b24-98d8-cc44c64f2e5a)

RepD3D is a prioritized graphical toolbox (also runnable as plain Python) that helps coastal and estuarine modellers

* **identify statistically representative simulation periods** from long-term wind records (Soares et al., 2024),
* **extract space- and time-varying boundary conditions** (water level, currents, waves, salinity, sediments and biogeochemistry, Neumann gradients, wind, pressure and temperature) from public hindcast and reanalysis data sets, and
* **write ready-to-use Delft3D input files** for both **Delft3D-4** and **Delft3D FM**.

## Versions in this repository

| Version | Where | Status |
|---------|-------|--------|
| **v2.0** (current) | `source_code_v2/`, `tutorial_v2/`, executable on the [Releases](https://github.com/capt-clay10/RepD3D/releases) page | Actively developed |
| **v1.0** | `source_code_v1/`, `Tutorial_files_v1/`, release [v1.0.0](https://github.com/capt-clay10/RepD3D/releases/tag/v1.0.0) | Archived; this is the version described in the MethodsX paper |

v2 is a rewrite with a new interface and new data sources. v1 remains available unchanged for anyone reproducing the published workflow. Every archived version has a DOI on Zenodo: <https://doi.org/10.5281/zenodo.23361961> always resolves to the latest one.

<p align="center">
  <img src="https://github.com/user-attachments/assets/e6c2d1b7-dcac-491e-adca-170dd708c083" alt="RepD3D v2.0 launcher with the five modules" width="420">
</p>

## What's new in v2.0

* **More data sources:** EasyGSH-DB, TrilaWatt, Copernicus Marine (NW Shelf physics, waves, biogeochemistry) and a *Custom* preset for other NetCDF files; wind, pressure and temperature from COSMO-REA6, ICON-DREAM-EU, ICON-DREAM-Global and ERA5.
* **Delft3D FM output** (`.bc`, `.pli`, `_bnd.ext`) next to the Delft3D-4 formats (`.bct`, `.bcw`, `.bnd`, `.mdw`); FM outputs still need testing.
* **Multi-year extraction in one run** across source files; the year-overlap tools of v1 are no longer needed.
* **Cross-border extraction** from several source files or cutouts (for example the German and Danish TrilaWatt grids) as prioritised source groups in one boundary set.
* **Quality control of every support point:** values come from wet cells only, dry periods are bridged without extrapolation, points on drying cells are moved along their boundary section, and implausible values stop the run with the point name. Everything is reported in the log.
* **Neumann (alongshore water-level gradient) boundaries** for Delft3D-4, derived from a parent model.
* **Synthetic boundary files** from a few hand-entered control points, and sea-level-change adjustment of existing files.
* **Temperature fields** (`.amt`) in addition to wind and pressure (`.amu`, `.amv`, `.amp`, `.wnd`); validation time series (CSV) at user-defined points.
* A data inspector lists the variables of a loaded file and picks the matching preset; the settings of each module can be saved and reloaded.
* Two corrections relative to v1 that change results: the COSMO-REA6 wind direction (about 2 degrees in the German Bight) and the meteo grid origin keyword (fields were placed half a cell to the north-east). Details in the [release notes](RELEASE_NOTES_v2.0.0.md).

## Modules

### Representative period

Scans a long-term wind record (CSV with `date` as `YYYYMMDDHH`, `speed` in m/s, `dir` in degrees) for windows that best reproduce the wind climate of a reference period. Window lengths are pandas offset strings (`12MS` = 12 calendar months, `6W`, `30D`), any subset of the 16 direction sectors and 12 Beaufort classes can be selected (all-class or selective-class correlation, Soares et al., 2024), and the result is a ranked list of candidate periods.

![Representative period window with the tutorial wind record](https://github.com/user-attachments/assets/2641c415-47c4-4624-b30f-30ac1974dd70)


The wind-rose viewer shows the reference period, each candidate and an overlay of the two, so a candidate can be judged by eye as well as by score.

![Wind-rose comparison of the reference period and a candidate](https://github.com/user-attachments/assets/71ba097f-5488-4e7c-ab34-f38dcbc4dcb8)

### Boundary-condition extraction

Reads a Delft3D-4 (`.mdf`, `.grd`, `.bnd`, `.mdw`) or Delft3D FM set-up, finds every open-boundary support point, and extracts water level, currents, waves, salinity, sediments or biogeochemistry for the simulation period from the chosen data set. Several files can be combined as source groups for cross-border domains. The log reports the quality check of every point, including points that were moved along their section because they sat on a drying cell.

![Boundary-condition extraction of EasyGSH water levels for the tutorial model](https://github.com/user-attachments/assets/9415b4ba-e78d-4386-802f-77be2788bf85)


For waves the same window takes the `.mdw` template and writes the `.bcw` with the time column aligned to the model reference date read from the `.mdf`.

![Wave extraction with the .mdw template and wave timing source](https://github.com/user-attachments/assets/916e42d4-e98a-45ba-8c48-ae7076f8082a)


### Meteo file generator

Samples wind, pressure or temperature from COSMO-REA6, ICON-DREAM-EU, ICON-DREAM-Global or ERA5 onto an equidistant Delft3D grid exported from QuickPlot (`.mat` v7) and writes the `.amu`, `.amv`, `.amp` or `.amt` files, optionally with the `xwind`/`ywind` `.wnd` variants and validation time series at named points. COSMO winds are rotated from the rotated-pole grid to true north. Enter `dx` and `dy` as the spacing of your grid; they are written to the file header.

![Meteo file generator with COSMO-REA6 input](https://github.com/user-attachments/assets/9a301c4e-7854-40f9-a79b-37e0a3ec90fa)


### Sea-level change

Adds a constant offset or a linear rise to every section of an existing `.bct`.

![Sea-level change applied to a .bct file](https://github.com/user-attachments/assets/998a0352-06ff-4f44-893b-7fea7f6f1e90)


### Synthetic boundary files

Builds a boundary file from a few control dates and values per boundary end (linear or step-hold), using the model files to get the boundary names. Useful for tests and scenario runs.

![Synthetic water-level boundary file from control points](https://github.com/user-attachments/assets/f194d292-29e9-44cc-858f-b11a037e0b64)


## Quick start (executable)

1. Download `RepD3D_v2.0.0.exe` from the [Releases](https://github.com/capt-clay10/RepD3D/releases) page.
2. Double-click to launch. The first start takes a moment while the single-file executable unpacks, and Windows SmartScreen may ask for confirmation because the file is unsigned.
3. Pick a module, point it at your data and Delft3D files, and follow the numbered steps. The log panel on the right shows progress, timing and quality-check messages.

All times in RepD3D input and output are UTC.

## Tutorial files

`tutorial_v2/` holds small inputs to try each module: the tutorial wind record for the representative period, the German Bight example model (`.grd`, `.bnd`, `.mdf`, `.mdw`) for boundary extraction, the 6 km target grid for the meteo generator, and a `.bct` for the sea-level module. The large source data (EasyGSH-DB, COSMO-REA6 and so on) is not included; download it from the providers listed below. See [`tutorial_v2/README.md`](tutorial_v2/README.md).

## Run from source

Everything the application needs lives in one flat folder, `source_code_v2/`: the Python modules, the icons, and the auxiliary grid files for the COSMO-REA6 and ICON-DREAM presets (`cosmo_latlon_rot.grb`, `ICON-DREAM-EU_grid.nc`, `ICON-DREAM-Global_grid.nc`). Keep that folder intact and start `main.py` from inside it.

```bash
git clone https://github.com/capt-clay10/RepD3D.git
cd RepD3D

# eccodes (needed by cfgrib for GRIB input) is easiest from conda-forge
conda create -n repd3d python=3.14 eccodes -c conda-forge
conda activate repd3d
pip install -r requirements.txt

cd source_code_v2
python main.py
```

`requirements.txt` pins the versions the v2.0.0 release was built and tested with (Windows, Python 3.14). Other platforms are untested, but nothing in the code is Windows-specific, except for the packaged executable.

## Supported data sets

| Data set | Used for | Source |
|----------|----------|--------|
| **EasyGSH-DB** (1000 m, German Bight) | water level, waves | <https://mdi-de.baw.de/easygsh/> |
| **TrilaWatt** (trilateral Wadden Sea) | water level, currents, waves | <https://dl.datenrepository.baw.de/7000/B3955.02.04.70237/> |
| **Copernicus Marine** NW Shelf (physics, waves, biogeochemistry) | SSH, currents, waves, salinity, SPM | <https://data.marine.copernicus.eu/> |
| **COSMO-REA6** (DWD, 6 km, Europe) | wind, pressure, temperature | <https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/> |
| **ICON-DREAM-EU / ICON-DREAM-Global** (DWD, 6.5 km / 13 km) | wind, pressure, temperature | <https://opendata.dwd.de/climate_environment/REA/> |
| **ERA5** (ECMWF, 0.25 degrees, global), NetCDF download | wind, pressure, temperature | <https://cds.climate.copernicus.eu/> |
| **Custom NetCDF** | any gridded variable on a lat/lon grid | – |

Please read the documentation of each dataset to understand how it was produced, and cite it when you use it (see below).

## Good to know

* ERA5 must be supplied as NetCDF, the format offered by the Climate Data Store; ERA5 GRIB is not read.
* Cartesian grids are assumed to be in the UTM zone entered in the form (default 32N).
* Boundary names in `.bnd` files are expected in the form `Name_1` (FLOW) and `Name1` (WAVE).
* Large source files are read lazily; the first access to a 25 GB EasyGSH wave file can take a minute.

## Build the Windows executable

The release executable is a single-file, windowed PyInstaller build made from `source_code_v2/`. From the repository root, in the environment above with `pyinstaller` installed:

```
pyinstaller --noconfirm --clean --onefile --windowed --name RepD3D-v2_win ^
  --icon source_code_v2/repd3d.ico ^
  --add-data "source_code_v2/*.py;." --add-data "source_code_v2/*.ico;." --add-data "source_code_v2/*.png;." ^
  --add-data "source_code_v2/*.jpg;." --add-data "source_code_v2/*.grb;." --add-data "source_code_v2/*.nc;." ^
  --add-binary "<env>/Library/bin/eccodes.dll;." --add-binary "<env>/Library/bin/eccodes_memfs.dll;." ^
  --collect-all customtkinter --collect-all pyproj --collect-all cfgrib --collect-all eccodes ^
  --collect-all netCDF4 --collect-all dask --collect-all matplotlib ^
  --hidden-import backend_bc_extractor_impl --hidden-import backend_neumann --hidden-import gui_windrose ^
  --hidden-import xarray.namedarray.daskmanager --hidden-import xarray.backends.netCDF4_ ^
  --hidden-import xarray.backends.scipy_ --hidden-import pyproj.datadir --hidden-import PIL._tkinter_finder ^
  source_code_v2/main.py
```

`auto_py_to_exe_settings.json` holds the same options for the auto-py-to-exe GUI; replace the placeholder path of the two eccodes DLLs with your own environment first. Notes:

* `--collect-all dask` is required: nothing imports dask statically, but multi-file extraction uses it through xarray and the frozen build fails at that point without it.
* The ICON grid files in `source_code_v2/` contain only the cell-centre coordinates (`clat`, `clon`) that RepD3D reads; the full DWD grid files are not needed.
* Keep `source_code_v2/` flat. The asset lookup and the build options assume every module and asset sits directly in that folder.

## Citation

If you use RepD3D in your research, please cite:

```bibtex
@article{Soares2025,
  title   = {RepD3D: A tool for representative period identification and associated boundary condition extraction},
  journal = {MethodsX},
  volume  = {14},
  pages   = {103109},
  year    = {2025},
  author  = {C. C. Soares and A. Knies and C. Winter},
  doi     = {10.1016/j.mex.2024.103109}
}
```

* Representative-period algorithm: **Soares, C.C., Galiforni-Silva, F., Winter, C., 2024. Representative residual transport pathways in a mixed-energy open tidal system. Journal of Sea Research 201, 102530. https://doi.org/10.1016/j.seares.2024.102530**
* Toolbox: **Soares, C. C., Knies, A., & Winter, C. (2025). RepD3D: A tool for representative period identification and associated boundary condition extraction. MethodsX, 14, 103109. https://doi.org/10.1016/j.mex.2024.103109**
* Software archive: **https://doi.org/10.5281/zenodo.23361961** (concept DOI; use the version DOI shown on Zenodo for the release you used)

### Please also acknowledge the underlying data sets

* **EasyGSH-DB:** Hagen, R., Plüß, A., Schrage, N., Dreier, N. (2020): EasyGSH-DB: Themengebiet - synoptische Hydrodynamik. Bundesanstalt für Wasserbau. https://doi.org/10.48437/02.2020.K2.7000.0004
* **COSMO-REA6:** Bollmeyer, C., Keller, J.D., Ohlwein, C., Wahl, S., Crewell, S., Friederichs, P., Hense, A., Keune, J., Kneifel, S., Pscheidt, I., Redl, S., Steinke, S., 2015. Towards a high-resolution regional reanalysis for the European CORDEX domain. Q.J.R. Meteorol. Soc. 141, 1–15. https://doi.org/10.1002/qj.2486
* **ICON-DREAM:** Deutscher Wetterdienst, ICON-DREAM: ICON Dual-resolution Reanalysis for Emulators, Applications and Monitoring, version 1. https://doi.org/10.5676/dwd/icon-dream_v1
* **ERA5:** Hersbach, H., et al., 2020. The ERA5 global reanalysis. Q.J.R. Meteorol. Soc. 146, 1999–2049. https://doi.org/10.1002/qj.3803
* **TrilaWatt:** Bundesanstalt für Wasserbau data repository, https://dl.datenrepository.baw.de/7000/B3955.02.04.70237/ (cite as given on the data page).
* **Copernicus Marine Service:** cite the product DOI shown on the product page of each data set you download.

## Contributing

Bug reports, feature requests and pull requests are welcome. Please open an issue first to discuss your idea.

## License

Released under the **MIT License** – see [`LICENSE`](LICENSE).

## Contact

Clayton C. Soares – `clayton.soares@ifg.uni-kiel.de`
Christian Winter
Arne Knies

> *Always validate numerical-model results yourself!* RepD3D provides no guarantees of suitability for any specific purpose.
