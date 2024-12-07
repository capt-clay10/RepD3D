# -*- coding: utf-8 -*-
"""
Converting COSMO files to Delft3D format
"""


def create_wind_fields_cosmo(grid_ed_path, const_path, cosmo_files_path,
                             file_name, ref_time, crs_type):

    import os
    import numpy as np
    import scipy.io
    from scipy.interpolate import griddata
    import datetime as dt
    import cfgrib
    import xarray as xr
    from tqdm import tqdm
    from pyproj import Proj, Transformer
    import warnings
    # Suppress specific warning about ecCodes version
    warnings.filterwarnings(
        "ignore", message="ecCodes 2.31.0 or higher is recommended. You are running version 2.26.0")

    # %% Functions
    
    # Calculate rotation angle between cosmo and true north
    def calculate_theta(rlat, rlon):
        """
        Calculate the rotation angle theta for each grid point using rotated lat/lon
        and the rotated pole information.
        Inputs:
            rlat, rlon: Rotated latitude and longitude grids (2D arrays, in degrees)
            latsp, lonsp: Latitude and longitude of the rotated pole (scalars, in degrees)
        Output:
            theta: Rotation angle for each grid point (2D array, in radians)
        """
        # Pole information extracted from section 6
        # https://opendata.dwd.de/climate_environment/REA/COSMO_REA6/help_COSMO_REA6/COSMO_REA6_Starting_example.pdf
        lonsp = -162  # Rotated pole longitude
        latsp = 39.25  # Rotated pole latitude
        
        # Convert degrees to radians
        rlat_rad = np.radians(rlat)
        rlon_rad = np.radians(rlon)
        latsp_rad = np.radians(latsp)
        lonsp_rad = np.radians(lonsp)

        # Compute theta
        theta = np.arctan2(
            np.sin(rlon_rad - lonsp_rad),
            np.cos(latsp_rad) * np.tan(rlat_rad) - np.sin(latsp_rad) * np.cos(rlon_rad - lonsp_rad)
        )
        return theta
    
    # Function to process the COSMO data and write it to files
    def process_data(year, cosmo_dir, X, Y, xcos, ycos, theta, files, REF):
        """
        Extract and interpoalte COSMO data onto a user defined equidistant grid
        Inputs:
            year: starting year of cosmo data
            cosmo_dir = directory of cosmo files (follow recommended folder structure)
            X = Delft3D grid X
            Y = Delft3D grid Y
            xcos, ycos = Rotation corrected UTM converted COSMO
            theta = Angle of rotation for correcting U and V component
            files = Blank document for amu amv and amp files
            REF = reference time for wind files
        Output:
            theta: Rotation angle for each grid point (2D array, in radians)
        """
        # List all UV and PS files in the COSMO directory
        A = [f for f in os.listdir(os.path.join(
            cosmo_dir, 'UV')) if f.endswith('.grb')]
        B = [f for f in os.listdir(os.path.join(
            cosmo_dir, 'PS')) if f.endswith('.grb')]

        for i in tqdm(range(0, len(A) // 2), desc='Extracting wind and pressure fields', total=(len(A) // 2), leave=True, mininterval=0.1):
            # Define file paths for U, V, and PS components
            WU_file = os.path.join(cosmo_dir, 'UV/', 'U_' + A[i][2:])
            WV_file = os.path.join(
                cosmo_dir, 'UV/', 'V_' + A[i + len(A) // 2][2:])
            PS_file = os.path.join(cosmo_dir, 'PS/', B[i])

            # Open GRIB files using cfgrib
            WU = cfgrib.open_dataset(WU_file)
            WV = cfgrib.open_dataset(WV_file)
            PS = cfgrib.open_dataset(PS_file)

            for kk in range(len(WU.time)):
                # Process components
                DATA_U = WU.u10.values[kk] # vector
                DATA_V = WV.v10.values[kk] # vector
                DATA_P = PS.sp.values[kk] # scalar
                
                # Slice data with delft grid overlap
                wnd_cosu = DATA_U[mm1:mm2, nn1:nn2]
                wnd_cosv = DATA_V[mm1:mm2, nn1:nn2]

                # Rotate U and V components to align with the geographical grid
                U_corrected = wnd_cosu * np.cos(theta) - wnd_cosv * np.sin(theta)
                V_corrected = wnd_cosu * np.sin(theta) + wnd_cosv * np.cos(theta)
                
                # Interpolate the rotated U and V components onto the Delft3D grid
                U_interp = griddata(
                    (xcos.flatten(), ycos.flatten()), U_corrected.flatten(), (X, Y), method='linear'
                )
                V_interp = griddata(
                    (xcos.flatten(), ycos.flatten()), V_corrected.flatten(), (X, Y), method='linear'
                )
                P_interp = griddata(
                    (xcos.flatten(), ycos.flatten()), DATA_P[mm1:mm2, nn1:nn2].flatten(), (X, Y), method='linear'
                ) / 100  # Convert to mbar
    
                # Replace NaN values with placeholder
                U_interp[np.isnan(U_interp)] = 9999.00
                V_interp[np.isnan(V_interp)] = 9999.00
                P_interp[np.isnan(P_interp)] = 9999.00

                # Calculate time difference and format time string
                T = WU.time.values[kk].astype('datetime64[s]').tolist()
                diff = round((T - REF).total_seconds() / 3600)
                t = f'TIME = {diff} hours since {REF.strftime("%Y-%m-%d %H:%M:%S")} +00:00 '

                # Write time header to files
                for key in ['amu', 'amv', 'amp', 'xwind', 'ywind']:
                    files[key].write(f'{t}\n')

                # Write data to files
                write_data(files['amu'], U_interp, a, b)
                write_data(files['amv'], V_interp, a, b)
                write_data(files['amp'], P_interp, a, b, is_pressure=True)
                write_data(files['xwind'], U_interp, a, b)
                write_data(files['ywind'], V_interp, a, b)

    def write_data(fid, data, a, b, is_pressure=False):
        if is_pressure:
            Fa = '%6.2f'
            Fb = ' %7.2f'
            Fc = ' %7.2f\n'
        else:
            Fa = '%6.2f'
            Fb = ' %7.2f'
            Fc = ' %7.2f\n'

        Fb1 = ''.join([Fb] * (b - 2))
        Zeil = [Fa + Fb1 + Fc]
        Fd = '%6.2f'
        ZeilP = [Fd + Fb1 + Fc]

        for jj in range(a):
            if jj == 0:
                fid.write(Zeil[0] % tuple(data[jj, :]))
            else:
                fid.write(ZeilP[0] % tuple(data[jj, :]))

    # %% inputs

    # Set paths and reference datetime
    dir_cosmo = cosmo_files_path

    REF = dt.datetime.strptime(ref_time, '%Y-%m-%d %H:%M:%S')
    yr = REF.year

    
    # %% Start pre-processing
    if crs_type == 'Cartesian':
        # Step 1: Load the COSMO rotated and true lat and long
        const_file = const_path
    
        ds = xr.open_dataset(
            const_file,
            engine="cfgrib",
            filter_by_keys={'typeOfLevel': 'surface'}  
        )
    
        # Extract all coordinates
        lat = ds['RLAT'].values  # True latitudes
        lon = ds['RLON'].values  # True longitudes
        rlat  =  ds['latitude'].values  # Rotated latitudes
        rlon  =  ds['longitude'].values  # Rotated latitudes 
        
        # Step 2 Calculate rotation angle theta between rotated and true geographic
        theta = calculate_theta(rlat, rlon)
    
        # Step 3 : Convert True LAT and LON to UTM (32N) for finding overlap
        # Define the transformer for lat/lon to UTM 32N
        latlon_to_utm = Transformer.from_proj(
            Proj(proj="latlong", datum="WGS84"),  # Source: Latitude/Longitude (WGS84)
            Proj(proj="utm", zone=32, datum="WGS84")  # Target: UTM Zone 32N
        )
    
        x_cosmo, y_cosmo = latlon_to_utm.transform(lon.flatten(), lat.flatten())
        x_cosmo = x_cosmo.reshape(lon.shape)
        y_cosmo = y_cosmo.reshape(lat.shape)
    
        # Step 4: Load the Delft3D equidistant grid file in UTM 32N
        grid_ed = scipy.io.loadmat(grid_ed_path)
    
        # Extract Delft3D grid coordinates
        X = grid_ed['data']['X'][0][0][:-1, :-1]  # Longitude in UTM
        Y = grid_ed['data']['Y'][0][0][:-1, :-1]  # Latitude in UTM
    
        # Extract Delft3D grid information for file headers (amu,amv,amp)
        n_cols = X.shape[1] # number of columns
        n_rows = X.shape[0] # number of rows
        x_llcorner = X.min() # Lower left x ccordinate from delft morpho grid
        y_llcorner = Y.min() # Lower left y ccordinate 
        unit = 'm' # UTM
    
        # Step 5: Identify the min/max extent bounds of the Delft3D grid
        x_min, x_max = np.nanmin(X), np.nanmax(X)  # Lon
        y_min, y_max = np.nanmin(Y), np.nanmax(Y)  # Lat
    
        # Step 6: Identify the indices for the COSMO data
        # Create a mask to find the indices that fall within the Delft3D extent
        x_mask = (x_cosmo >= x_min) & (x_cosmo <= x_max) # Lon
        y_mask = (y_cosmo >= y_min) & (y_cosmo <= y_max) # Lat
        final_mask = x_mask & y_mask
    
        # Find the bounding box of the matching indices
        valid_indices = np.where(final_mask)
    
        # A buffer is added to make sure a complete overlap
        buffer = 2
        mm1 = max(0, np.min(valid_indices[0]) - buffer)
        mm2 = min(x_cosmo.shape[0], np.max(valid_indices[0]) + buffer)
        nn1 = max(0, np.min(valid_indices[1]) - buffer)
        nn2 = min(x_cosmo.shape[1], np.max(valid_indices[1]) + buffer)
    
        # Step 7: Slice theta and coordinates for overlapping region
        theta = theta[mm1:mm2, nn1:nn2]
        xcos = x_cosmo[mm1:mm2, nn1:nn2] # Rotation corrected UTM 
        ycos = y_cosmo[mm1:mm2, nn1:nn2] # Rotation corrected UTM 
        
    else: # Spherical
    
        # Step 1: Load the COSMO rotated and true lat and long
        const_file = const_path
    
        ds = xr.open_dataset(
            const_file,
            engine="cfgrib",
            filter_by_keys={'typeOfLevel': 'surface'}  
        )
    
        # Extract all coordinates
        lat = ds['RLAT'].values  # True latitudes
        lon = ds['RLON'].values  # True longitudes
        rlat  =  ds['latitude'].values  # Rotated latitudes
        rlon  =  ds['longitude'].values  # Rotated latitudes 
        
        # Step 2 Calculate rotation angle theta between rotated and true geographic
        theta = calculate_theta(rlat, rlon)
        
        x_cosmo = lon
        y_cosmo = lat
        
        # Step 3: Load the Delft3D equidistant grid file in spherical
        grid_ed = scipy.io.loadmat(grid_ed_path)
    
        # Extract Delft3D grid coordinates
        X = grid_ed['data']['X'][0][0][:-1, :-1]  # Longitude 
        Y = grid_ed['data']['Y'][0][0][:-1, :-1]  # Latitude 
        
        # Extract Delft3D grid information for file headers (amu,amv,amp)
        n_cols = X.shape[1] # number of columns
        n_rows = X.shape[0] # number of rows
        x_llcorner = X.min() # Lower left x ccordinate from delft morpho grid
        y_llcorner = Y.min() # Lower left y ccordinate 
        unit = 'deg' # Spherical
        
        # Step 4: Identify the min/max extent bounds of the Delft3D grid
        x_min, x_max = np.nanmin(X), np.nanmax(X)  # Lon
        y_min, y_max = np.nanmin(Y), np.nanmax(Y)  # Lat
    
        # Step 5: Identify the indices for the COSMO data
        # Create a mask to find the indices that fall within the Delft3D extent
        x_mask = (x_cosmo >= x_min) & (x_cosmo <= x_max) # Lon
        y_mask = (y_cosmo >= y_min) & (y_cosmo <= y_max) # Lat
        final_mask = x_mask & y_mask
    
        # Find the bounding box of the matching indices
        valid_indices = np.where(final_mask)
    
        # A buffer is added to make sure a complete overlap
        buffer = 2
        mm1 = max(0, np.min(valid_indices[0]) - buffer)
        mm2 = min(x_cosmo.shape[0], np.max(valid_indices[0]) + buffer)
        nn1 = max(0, np.min(valid_indices[1]) - buffer)
        nn2 = min(x_cosmo.shape[1], np.max(valid_indices[1]) + buffer)
    
        # Step 6: Slice theta and coordinates for overlapping region
        theta = theta[mm1:mm2, nn1:nn2]
        xcos = x_cosmo[mm1:mm2, nn1:nn2] # Rotation corrected 
        ycos = y_cosmo[mm1:mm2, nn1:nn2] # Rotation corrected 
        
    
    # %% start file generation

    # Open files for writing processed data
    files = {
        'amu': open(os.path.join(dir_cosmo, f'{file_name}.amu'), 'w+'),
        'amv': open(os.path.join(dir_cosmo, f'{file_name}.amv'), 'w+'),
        'amp': open(os.path.join(dir_cosmo, f'{file_name}.amp'), 'w+'),
        'xwind': open(os.path.join(dir_cosmo, f'xwind_{file_name}.wnd'), 'w+'),
        'ywind': open(os.path.join(dir_cosmo, f'ywind_{file_name}.wnd'), 'w+')
    }

    # Replace zero values in X and Y with NaN
    X[X == 0] = np.nan
    Y[Y == 0] = np.nan
    a, b = X.shape  # Get the shape of the grid

    # Define headers for different data files
    HEADER = {
        'amu': [
            'FileVersion = 1.03',
            'Filetype = meteo_on_equidistant_grid',
            'NODATA_value = 9999.00',
            f'n_cols = {n_cols}',
            f'n_rows = {n_rows}',
            f'grid_unit = {unit}',
            f'x_llcorner = {x_llcorner}',
            f'y_llcorner = {y_llcorner}',
            'dx = 6000',
            'dy = 6000',
            'n_quantity = 1',
            'quantity1 = x_wind',
            'unit1 = m s-1'
        ],
        'amv': [
            'FileVersion = 1.03',
            'Filetype = meteo_on_equidistant_grid',
            'NODATA_value = 9999.00',
            f'n_cols = {n_cols}',
            f'n_rows = {n_rows}',
            f'grid_unit = {unit}',
            f'x_llcorner = {x_llcorner}',
            f'y_llcorner = {y_llcorner}',
            'dx = 6000',
            'dy = 6000',
            'n_quantity = 1',
            'quantity1 = y_wind',
            'unit1 = m s-1'
        ],
        'amp': [
            'FileVersion = 1.03',
            'Filetype = meteo_on_equidistant_grid',
            'NODATA_value = 9999.00',
            f'n_cols = {n_cols}',
            f'n_rows = {n_rows}',
            f'grid_unit = {unit}',
            f'x_llcorner = {x_llcorner}',
            f'y_llcorner = {y_llcorner}',
            'dx = 6000',
            'dy = 6000',
            'n_quantity = 1',
            'quantity1 = air_pressure',
            'unit1 = mbar'
        ],
        'xwind': [
            'FileVersion = 1.03',
            'Filetype = meteo_on_equidistant_grid',
            'NODATA_value = 9999.00',
            f'n_cols = {n_cols}',
            f'n_rows = {n_rows}',
            f'grid_unit = {unit}',
            f'x_llcorner = {x_llcorner}',
            f'y_llcorner = {y_llcorner}',
            'dx = 6000',
            'dy = 6000',
            'n_quantity = 1',
            'quantity1 = x_wind',
            'unit1 = m s-1'
        ],
        'ywind': [
            'FileVersion = 1.03',
            'Filetype = meteo_on_equidistant_grid',
            'NODATA_value = 9999.00',
            f'n_cols = {n_cols}',
            f'n_rows = {n_rows}',
            f'grid_unit = {unit}',
            f'x_llcorner = {x_llcorner}',
            f'y_llcorner = {y_llcorner}',
            'dx = 6000',
            'dy = 6000',
            'n_quantity = 1',
            'quantity1 = y_wind',
            'unit1 = m s-1'
        ]
    }

    for key in HEADER.keys():
        for line in HEADER[key]:
            files[key].write(line + '\n')

    # %% Run the data processing for the year
    process_data(yr, dir_cosmo, X, Y, xcos, ycos, theta, files, REF)

    # Close files
    for f in files.values():
        f.close()

