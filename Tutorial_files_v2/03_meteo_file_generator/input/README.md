# Meteo file generator – inputs

The target grid is included. The COSMO-REA6 source data is not (about 1 GB per file); download it and point the module at the folder.

## Included

| File | Used for |
|------|----------|
| `Delft3D_grid_6x6.mat` | Equidistant 6 km x 6 km target grid, 36 x 47 nodes, UTM 32N, exported from Delft3D QuickPlot as MATLAB v7 `.mat`. The values are sampled at its nodes. |

## Download separately

COSMO-REA6 hourly 2D fields from <https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/> for the months September to December 2013 (the example model runs 2013-10-20 to 2013-12-25). Put them in one folder with the layout the DWD server uses:

```
cosmo_wind_2013/
    UV/
        U_10M.2D.201309.grb   U_10M.2D.201310.grb   U_10M.2D.201311.grb   U_10M.2D.201312.grb
        V_10M.2D.201309.grb   V_10M.2D.201310.grb   V_10M.2D.201311.grb   V_10M.2D.201312.grb
    PS/
        PS.2D.201309.grb      PS.2D.201310.grb      PS.2D.201311.grb      PS.2D.201312.grb
```

On the server the files are under `U_10M/`, `V_10M/` and `PS/`. A flat folder with all twelve files also works. The COSMO coordinate file (`cosmo_latlon_rot.grb`) is packaged with RepD3D and does not need to be downloaded.

## Settings used for the expected outputs

| Setting | Value |
|---------|-------|
| Extraction file type | Wind fields (u, v, p) |
| Source | COSMO REA6 |
| Data directory | the folder above |
| `.mat` grid file | `Delft3D_grid_6x6.mat` |
| CRS, UTM zone | Cartesian, 32 |
| dx, dy | 6000, 6000 (the spacing of the `.mat` grid) |
| Reference time | 2013-01-01 00:00:00 (the example model's reference date) |
| Output filename | `cosmo_2013` |
| Also generate .xwind/.ywind files | on (`xwind_cosmo_2013.wnd` and `ywind_cosmo_2013.wnd` are byte-identical copies of `.amu` and `.amv`; the example wave model references them) |
| Validation point | `helgoland`, 54.1832 N, 7.8856 E |

The results are in `../output/`.
