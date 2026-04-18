# RepD3D – Representative‐period Detection & Delft3D Boundary‑data Toolbox

![GitHub release](https://img.shields.io/github/v/release/capt-clay10/RepD3D?style=for-the-badge)
![License](https://img.shields.io/github/license/capt-clay10/RepD3D?style=for-the-badge)
![Python](https://img.shields.io/badge/python-3.9%2B-blue?style=for-the-badge)

![RepD3D logo](https://github.com/user-attachments/assets/b1de47b1-be95-4b24-98d8-cc44c64f2e5a)

RepD3D is a Windows‑first graphical user interface (GUI) and open‑source Python library that helps coastal and estuarine modelers

* **identify statistically representative simulation periods** from long‑term wind records using the algorithm of *Soares et al., 2024*.
* **extract space‑ and time‑varying boundary conditions** (water level, waves, wind & pressure) from public hindcast data sets, and
* **write fully‑formatted Delft3D input files** (`.bct`, `.bcw`, `.amu/.amv/.amp`, `.wnd`).

The toolbox is packaged as a portable **RepD3D.exe** file, as well as pure Python code (`main.py`), which runs on any platform with the required libraries.

## In development
* Inclusion of TrilaWatt dataset (https://dl.datenrepository.baw.de/7000/B3955.02.04.70237/)
* Inclusion of Icon-EU wind fields (https://huggingface.co/datasets/openclimatefix/dwd-icon-eu/tree/main)
* Inclusion of a universal boundary condition extractor (add source file -> select variable for extraction -> write to Delft3D compatible files) 

---

## Key features

| Module | Purpose | Output files |
|--------|---------|--------------|
| **RPI** | Identify all‑class *or* selective‑class representative periods with user‑defined window length. | Output rank list `.txt`, interactive wind‑rose viewer |
| **Extract all** | Generation of water‑level (`.bct`) *and* wave (`.bcw`) files, plus updated `.mdw` and CSV of boundary coordinates. | `.bct`, `.bcw`, updated `.mdw`, `*_locations.csv` |
| **Bct / Bcw year‑overlap** | Build multi‑year time‑series for water level or waves. | `.bct` / `.bcw` |
| **SLR adjuster** | Add gradual or constant sea‑level rise to existing `.bct`. | new `.bct` |
| **Wind‑field generator** | Convert COSMO‑REA6 (10 m U & V, surface P) to Delft3D wind/pressure grids on a user grid. | `.amu`, `.amv`, `.amp`, `xwind.wnd`, `ywind.wnd` |

---

## Quick start (GUI)

1. **Download** the latest `RepD3D.exe` from the [releases](https://github.com/⋯/releases) page and place it in your project folder.
2. **Double‑click** to launch.  A splash screen will appear followed by the main menu.
3. **Browse** to a working directory and select the desired sub‑module.
4. **Follow on‑screen prompts** to choose grids, boundaries, NetCDF files and parameters.
5. Output files are written to the working directory (all times in UTC).

> **Tip:** The GUI may appear unresponsive while heavy interpolation loops are running; a progress console at the bottom streams real‑time log output.

---

## Install from source (CLI)

```bash
# 1. Clone the repo
$ git clone https://github.com/capt-clay10/RepD3D.git
$ cd RepD3D

# 2. Create environment (Python ≥3.9)
$ conda create -n repd3d python=3.9
$ conda activate repd3d

# 3. Install dependencies
$ pip install -r requirements.txt  # see list below

# 4. Run command‑line interface
$ python main.py
```

### Python dependencies

```
ast cfgrib dask ecCodes h5py h5netcdf matplotlib numpy pandas pyproj scikit‑learn scipy tqdm utm xarray windrose ttkbootstrap pillow
```

> Some packages (`cfgrib`, `ecCodes`) require C‑libraries.  On Windows we recommend the conda‑forge channel.

---

## Required data sets

| Data set | Coverage | Used for | Source |
|----------|----------|----------|--------|
| **EasyGSH‑DB 1000 m** | German North Sea, 1996‑2016 | Water level & wave forcing | <https://mdi-de.baw.de/easygsh/> |
| **COSMO‑REA6 hourly 2D** | Europe, 1995‑2019 | Wind (U, V) & surface pressure | <https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/> |

Download files of interest into the working directory following the folder scheme below:

```
Working directory/
    YYYY_1000m_wave_2D.nc
    YYYY_1000m_waterlevel_2D
    COSMO_YYYY/
        ├─ PS/        # surface pressure GRIB     (PS.*.grb)
        ├─ UV/        # 10 m wind  GRIB           (U_*.grb, V_*.grb)
```

---

## Typical workflow

1. **Design grids & boundaries** in Delft3D‑RGFGRID / FLOW / WAVE; save `.grd`, `.bnd`, `.mdf`, `.mdw`.
2. **Run RPI** to choose an unfiltered, reduced representative period matching your study goals.
3. **Use Extract modules** to build boundary time‑series for that period.
4. **(Optional) Add SLR** or generate COSMO wind fields.
5. **Launch Delft3D**, reference the generated files, and simulate!

---

## Screenshots

| Screenshot |
|------------|
| **Main menu**<br>![Main menu](https://github.com/user-attachments/assets/9c4f226a-a729-42ea-b14a-40902fdca57f) |

| **Representative-period viewer**<br>![Representative-period viewer](https://github.com/user-attachments/assets/6ac10213-0ff1-4087-bf77-10bf52bcd61a) |

| **EasyGSH extractor**<br>![EasyGSH extractor](https://github.com/user-attachments/assets/088161cd-965c-4e2c-b4b1-4d395344e5cd) |

| **COSMO wind extractor**<br>![COSMO wind extractor](https://github.com/user-attachments/assets/ea3ef267-2420-42e5-8378-68f9cee96505) |
---


## Citation

If you use RepD3D in your research, please cite the following:

```
@article{Soares2025,
  title   = {RepD3D: A tool for representative period identification and associated boundary condition extraction},
  journal = {MethodsX},
  volume  = {14},
  pages   = {103109},
  year    = {2025},
  author  = {C. C. Soares and A. Knies and C. Winter},
  doi     = {10.1016/j.mex.2024.103109}
}

* Citation for using Representative period algorithm (source paper): **Soares, C.C., Galiforni-Silva, F., Winter, C., 2024. Representative residual transport pathways in a mixed-energy open tidal system. Journal of Sea Research 201, 102530. https://doi.org/10.1016/j.seares.2024.102530**
* Citation for using the RepD3D tool box (source paper): **Soares, C. C., Knies, A., & Winter, C. (2025). RepD3D: A tool for representative period identification and associated boundary condition extraction. MethodsX, 14, 103109. https://doi.org/10.1016/j.mex.2024.103109**


```

### Please also acknowledge the underlying data sets:

* EasyGSH‑DB (Hagen et al., 2020)  
* COSMO‑REA6 (Bollmeyer et al., 2015)
  
* Source of easygsh data  **https://mdi-de.baw.de/easygsh/Easy_Viewer_syn.html#home**
* Citations for using EasyGSH data : **Hagen, R., Plüß, A., Schrage, N., Dreier, N. (2020): EasyGSH-DB: Themengebiet - synoptische Hydrodynamik. Bundesanstalt für Wasserbau. https://doi.org/10.48437/02.2020.K2.7000.0004**

  * Please read the source document to understand how EasyGSH datasets are generated. Here are some quick points.
    * The data provided are the results of a numerical simulation gridded over 1km and provided every 20 minutes. 
    * The numerical modelling approach used to generate the data utilizes annually updated bathymetry, tidal dynamics simulated by the Untrim2 modelling system, using tidal constituents at the open boundaries (corrected for external surge), waves computed using a combination of the model UnK (Schneggenburger et al., 2000) and SWAN for near-shore physical processes. **This code does not extract SWAN-generated data**

* Source for COSMO-REA6 data **https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/**
* Citations for using COMSO : **Bollmeyer, C., Keller, J.D., Ohlwein, C., Wahl, S., Crewell, S., Friederichs, P., Hense, A., Keune, J., Kneifel, S., Pscheidt, I., Redl, S., Steinke, S., 2015. Towards a high‐resolution regional reanalysis for the European CORDEX domain. Q.J.R. Meteorol. Soc. 141, 1–15. https://doi.org/10.1002/qj.2486**
   
  * Please read the source document to understand how COSMO-REA6 datasets are generated. **https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/help_COSMO_REA6/**
    * COSMO data is provided on a rotated pole, which is currently corrected in the algorithm to true north
    * COSMO data U_10m, V_10m and PS is used in this toolbox
    * Additionally a Delft3D grid of 6X6 km, extracted as .mat (v7) is required
---

## Contributing

Bug reports, feature requests and pull requests are welcome!  Please open an issue first to discuss your ideas.

> **Road‑map:** Support additional hindcast products (TriwaWaTT, ERA5), automatic download helpers, and native Linux GUI.

---

## License

This project is released under the **MIT License** – see [`LICENSE`](LICENSE).

---

## Contact

Clayton C. Soares – `clayton.soares@ifg.uni-kiel.de`  
Christian Winter  
Arne Knies

> *Always validate numerical‑model results yourself!*  RepD3D provides no guarantees for suitability for any specific purpose.
