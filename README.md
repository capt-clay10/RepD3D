![repd3d_logo](https://github.com/user-attachments/assets/b1de47b1-be95-4b24-98d8-cc44c64f2e5a)
# RepD3D: A tool for representative period identification and associated boundary condition extraction


*The tool is adapted for EasyGSH and COSMO extraction, but the source code is open and can be adapted to other datasets.*

Source of easygsh data  **https://mdi-de.baw.de/easygsh/Easy_Viewer_syn.html#home**
Citations for using EasyGSH data : **Hagen, R., Plüß, A., Schrage, N., Dreier, N. (2020): EasyGSH-DB: Themengebiet - synoptische Hydrodynamik. Bundesanstalt für Wasserbau. https://doi.org/10.48437/02.2020.K2.7000.0004**

Please read the source document to understand how EasyGSH datasets are generated. Here are some quick points.
* The data provided are the results of a numerical simulation gridded over 1km and provided every 20 minutes. 
* The numerical modelling approach used to generate the data utilizes annually updated bathymetry, tidal dynamics simulated by the Untrim2 modelling system, using tidal constituents at the open boundaries (corrected for external surge), waves computed using a combination of the model UnK (Schneggenburger et al., 2000) and SWAN for near-shore physical processes. **This code does not extract SWAN-generated data**

Source for COSMO data **https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/**

Citations for using COMSO : **Bollmeyer, C., Keller, J.D., Ohlwein, C., Wahl, S., Crewell, S., Friederichs, P., Hense, A., Keune, J., Kneifel, S., Pscheidt, I., Redl, S., Steinke, S., 2015. Towards a high‐resolution regional reanalysis for the European CORDEX domain. Q.J.R. Meteorol. Soc. 141, 1–15. https://doi.org/10.1002/qj.2486**

Citation for using Representative period algorithm (source paper): **Soares, C.C., Galiforni-Silva, F., Winter, C., 2024. Representative residual transport pathways in a mixed-energy open tidal system. Journal of Sea Research 201, 102530. https://doi.org/10.1016/j.seares.2024.102530**

**What does this code do?**
1) Identify Representative period based on wind data (all-class correlation and selective-class correlation) - Universally valid
2) Extract time-series water level 2D information for any designed boundaries within the EasyGSH model domain  (data found under the synoptic simulation, UnTRIM2, 1000m grid section.) - German North Sea only
3) Extract time series wave/Sea-state data 2D (significant height, peak period, direction, directional spread) for any designed boundaries within the EasyGSH model domain (data found under the synoptic simulation, UnTRIM2, 1000m grid section.) - German Nort Sea only
4) Extract spatial wind and pressure fields and create Delft3D-compatible wind files. - Europe wide

**Important notes on the code function**
* The code uses bilinear interpolation to extract the concerned point's water level and wave data.
* Please check the extent of your grid and the EasyGSH dataset limits before using the code, it has to fall within the limits of EasyGSH
* The wave direction is calculated according to Nautical convention and is in From-orientation
* The identification of a representative period is based on the algorithm explained in Soares et al.,(2024)
* The COSMO REA6 dataset is rotated to true north

### Packages used in this project

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
* requests
* sci-kit-learn
* scipy
* statistics
* sys 
* time
* tqdm
* utm 
* xarray


# Working on

* Adding variability in wave parameters for climate change scenarios.
* Including TrilaWatt dataset, an extension of EasyGSH dataset, coverage from Netherlands to Denmark
* Opportunity to create synthetic boundary conditions


### Two choices are presented, you can run the <ins>main.py</ins> script in your Python environment to extract your files or use the standalone GUI <ins>RepD3d.exe</ins>, alternatively, the gui_generator.py script is also provided in case you want to make changes to the source code and make a new GUI. Pyinstaller or Auto-py-to-exe can convert the <ins>gui_imporved.py</ins> to the executable GUI EasyD3d. 

* Snippet of the GUI startup

![Startup_GUI](https://github.com/user-attachments/assets/0e070b96-c9b1-4ff3-81c4-80b23c6271bc)

* Snippet of the Representative period module

![RPI_module_GUI](https://github.com/user-attachments/assets/f1b3c0b7-c6d9-4140-b3b6-cbdb8a8c5af4)

* Snippet of the Wind file generator- COSMO REA6

![COSMO_GUI](https://github.com/user-attachments/assets/978632dd-5e57-4648-af9e-244afd53cabc)

* Snippet of water and wave time series file generator module using EasyGSH

![Process_allfiles_GUI](https://github.com/user-attachments/assets/0a22c74f-dbb6-45f9-919b-3f39449277a4)

### Information about the BCW file

**An important thing to note is that this code replicates the file format for a SWAN wave model with uniform wave boundary conditions along one boundary line for multiple time points. so essentially if you make several boundaries one can create a space and time-varying effect. Disclaimer: read the manual**
**Inserting just the Cosmo Meteofiles also achieves good wave results**

**The direction is between 0-360, in the nautical convention**
**The directional spread is in degrees**

#### The creator of this script recommends the following steps to create and use the bcw file:

1) Use the Delft3d flow GUI to make boundaries on your wave grid (**boundary name should be the same as your wave boundary name**) and save the wave.bnd file. ***NOTE: Use this opportunity to make several boundaries this will create the effect of space-varying conditions. (at least 2km in length))***
2) Use the now-generated wave.bnd file and the wave.grd in the script to generate the boundary_location.csv
3) You can now run the script to generate the bcw file and the new mdw file.
4) Once completed, open the mdw file and add the keyword as in the manual with the appropriate format (**TSeriesFile= wave.bcw**)
5) Please follow the Delft3D manual instructions for proper implementation. 

* Note, first make a mdw file with a dummy boundary, since the code reads the .mdw file to look for the term 'Boundary' to store the boundary information
* Please use a boundary naming convention without an underscore so for example when making the wave.bnd file names the boundaries as North1 instead of North_1. 

### As a personal note from a fellow modeller ###
Always validate results yourself!
