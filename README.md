# ![repd3d_logo](https://github.com/user-attachments/assets/b1de47b1-be95-4b24-98d8-cc44c64f2e5a) RepD3D: A tool for representative period identification and associated boundary condition extraction

## Download RepD3D.exe or run the main.py file
  
### What does this code do?
1) Identify Representative period based on hindcast wind data (all-class correlation and selective-class correlation) - Universally valid
2) Extract time-series water level 2D information for any designed boundaries within the EasyGSH model domain  (data found under the synoptic simulation, UnTRIM2, 1000m grid section.) - German North Sea only
3) Extract time series wave/Sea-state data 2D (significant height, peak period, direction, directional spread) for any designed boundaries within the EasyGSH model domain (data found under the synoptic simulation, UnTRIM2, 1000m grid section.) - German Nort Sea only
4) Extract spatial wind and pressure fields and create Delft3D-compatible wind files. - Europe wide
#
### Data sources and citations
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
#
### Representative period identification
* Citation for using Representative period algorithm (source paper): **Soares, C.C., Galiforni-Silva, F., Winter, C., 2024. Representative residual transport pathways in a mixed-energy open tidal system. Journal of Sea Research 201, 102530. https://doi.org/10.1016/j.seares.2024.102530**
#
### Algorithm notes
* It uses bilinear interpolation to extract EasyGSH data
* The wave parameters are converted to nautical convention and wave-from orientation
* The COSMO data is corrected to true North
#
### Working on
* Adding variability in wave parameters for climate change scenarios.
* Including TrilaWatt dataset, an extension of EasyGSH dataset, coverage from Netherlands to Denmark
* Opportunity to create synthetic boundary conditions
#
* Snippet of the GUI startup

![Startup_GUI](https://github.com/user-attachments/assets/0e070b96-c9b1-4ff3-81c4-80b23c6271bc)

* Snippet of the Representative period module

![RPI_module_GUI](https://github.com/user-attachments/assets/f1b3c0b7-c6d9-4140-b3b6-cbdb8a8c5af4)

* Snippet of the Wind file generator- COSMO REA6

![COSMO_GUI](https://github.com/user-attachments/assets/978632dd-5e57-4648-af9e-244afd53cabc)

* Snippet of water and wave time series file generator module using EasyGSH
![Process_allfiles_GUI](https://github.com/user-attachments/assets/0a22c74f-dbb6-45f9-919b-3f39449277a4)

#
### Python environment requirements
* ast
* cfgrib
* csv
* dask
* datetime
* ecCodes
* h5py
* math
* numpy 
* os
* pandas
* pyproj
* re
* requests
* sci-kit-learn
* scipy
* statistics
* sys 
* time
* tqdm
* utm 
* xarray

### As a personal note from a fellow modeller ###
Always validate results yourself!
