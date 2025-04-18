# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 09:47:27 2023

@author: AG-Admin
"""

# -*- coding: utf-8 -*-
"""
Part two of the process
This script extracts boundary conditions - water level from the EasyGsh data
It requires the csv file generated from the part one of the process
It requires the netcdf file with the water level data
It requires the mdf file
The output file (.bct) will have the same name as the mdf file

"""

# %% Over head function


def bct_year_overlap_file_generator(boundaries, nc_files, mdf_file, step, bct_file_name, crs_type):
    # %% Import packages
    from datetime import timedelta
    from datetime import datetime
    import pandas as pd
    import numpy as np
    import os
    import csv
    from tqdm import tqdm
    from scipy import interpolate
    from scipy.interpolate import interp1d
    import re

    # %% input the files

    try:
        import utm
    except ModuleNotFoundError as err:
        # Error handling
        print(
            str(err) + (' Module utm doesnt exist please install it in your environment,',
                        'conda code: conda install - c conda-forge utm',
                        'pip code: pip install utm'))

    try:
        import xarray as xr
    except ModuleNotFoundError as err_4:
        # Error handling
        print(
            str(err_4) +
            ' This package also requires extra dependencies like netCDF4, h5netcdf and possibly scipy')

    # %% Create functions

    def value_from_txt_file(file, string_name):
        file1 = open(file, "r")
        for line in file1:
            # checking string is present in line or not
            if '=' in line:
                if string_name in line:
                    val = line.split('=')
                    string_val = val[1].strip()
                    break
                    file1.close()  # close file
                else:
                    print('{} is not in the file'.format(string_name))
        return string_val

    def convert_flt_to_sci_not(fltt, prec, exp_digits):
        s = "%.*e" % (prec, fltt)
        # print(f's: {s}')
        if s == 'nan':
            # TODO: is it a good idea to replace nan with 0?
            s = "%.*e" % (prec, 0)
        mantissa, exp = s.split('e')
        # print(f'mantissa: {mantissa}')
        # print(f'exp: {exp}')
        # add 1 to digits as 1 is taken by sign +/-
        return "%se%+0*d" % (mantissa, exp_digits + 1, int(exp))

    def convert_list_to_sci_not(input_list, prec, exp_digits):
        converted = []
        for flt in input_list:
            sci = convert_flt_to_sci_not(
                fltt=flt, prec=prec, exp_digits=exp_digits)
            converted.append(sci)

        return converted

    def add_blank_pos_val(input_str):
        """Add leading blank for positive value."""

        if input_str[:1] != '-':
            output_str = f' {input_str}'
        else:
            output_str = input_str
        return output_str

    def fill_zeros_with_interpolation(arr):
        # Identify the indices of the non-zero elements
        non_zero_indices = np.where(arr != 0)[0]
        zero_indices = np.where(arr == 0)[0]

        if non_zero_indices.size == 0:
            # If there are no non-zero values to interpolate from, return the original array
            return arr

        # Extract the non-zero elements and their indices
        non_zero_values = arr[non_zero_indices]

        # Create an interpolation function based on the non-zero elements
        interp_func = interp1d(non_zero_indices, non_zero_values,
                               kind='linear', fill_value="extrapolate")

        # Interpolate the zero elements
        arr[zero_indices] = interp_func(zero_indices)

        return arr

    # Ensure the nc_files list is ordered by year
    def extract_year(file_path):
        # Extract the file name from the path
        file_name = file_path.split('/')[-1].split('\\')[-1]  # Handle both Unix and Windows paths
        # Use regex to find the first numeric sequence (assumes it's the year)
        match = re.search(r'\d{4}', file_name)
        if match:
            return int(match.group(0))
        else:
            raise ValueError(f"Cannot extract year from file name: {file_name}")
        
    def fill_nans_with_interpolation(data_array):
        """Interpolate to fill NaNs in the time series data"""
        filled_data = data_array.copy()
        filled_data = np.nan_to_num(filled_data)  # Temporarily fill NaNs with 0
        return fill_zeros_with_interpolation(filled_data)
    
    def convert_coordinates(bnd_loc, crs_type, utm_zone=32, hemisphere='N'):
        """
        Convert boundary coordinates between Cartesian (UTM) and Spherical (Lat/Lon).
        
        Parameters:
        - bnd_loc: pd.DataFrame with columns ['boundary', 'easting', 'northing']
        - crs_type: str, either 'Cartesian' or 'Spherical'
        - utm_zone: int, UTM zone (default=32)
        - hemisphere: str, 'N' or 'S' for northern or southern hemisphere (default='N')
        
        Returns:
        - bnd_loc_geo: pd.DataFrame with columns ['lat', 'lon', 'boundaries']
        """
        # Initialize an empty DataFrame
        bnd_loc_geo = pd.DataFrame(columns=['lat', 'lon', 'boundaries'])
        
        if crs_type == 'Cartesian':
            # Extract UTM coordinates and boundary names
            easting = bnd_loc['easting'].astype(float)
            northing = bnd_loc['northing'].astype(float)
            bnd = bnd_loc['boundary']
            
            # Convert UTM to lat/lon
            geo_coords = [utm.to_latlon(e, n, utm_zone, hemisphere) for e, n in zip(easting, northing)]
            bnd_loc_geo = pd.DataFrame(geo_coords, columns=['lat', 'lon'])
            bnd_loc_geo['boundaries'] = bnd
        else:
            # For spherical coordinates, directly assign
            bnd_loc_geo['lat'] = bnd_loc['northing']
            bnd_loc_geo['lon'] = bnd_loc['easting']
            bnd_loc_geo['boundaries'] = bnd_loc['boundary']
        
        return bnd_loc_geo

    def bilinear_interpolation(lat, lon, ds, variable):
        lat_vals = ds.coords['lat'].values
        lon_vals = ds.coords['lon'].values

        lat_idx1 = np.searchsorted(lat_vals, lat) - 1
        lat_idx2 = lat_idx1 + 1
        lon_idx1 = np.searchsorted(lon_vals, lon) - 1
        lon_idx2 = lon_idx1 + 1

        lat_idx1 = np.clip(lat_idx1, 0, len(lat_vals) - 1)
        lat_idx2 = np.clip(lat_idx2, 0, len(lat_vals) - 1)
        lon_idx1 = np.clip(lon_idx1, 0, len(lon_vals) - 1)
        lon_idx2 = np.clip(lon_idx2, 0, len(lon_vals) - 1)

        lat1, lat2 = lat_vals[lat_idx1], lat_vals[lat_idx2]
        lon1, lon2 = lon_vals[lon_idx1], lon_vals[lon_idx2]

        Q11 = ds.sel(lat=lat1, lon=lon1)[variable].values

        Q12 = ds.sel(lat=lat1, lon=lon2)[variable].values
        Q21 = ds.sel(lat=lat2, lon=lon1)[variable].values
        Q22 = ds.sel(lat=lat2, lon=lon2)[variable].values

        # Check 1: Replace NaNs with 0 and interpolate
        if np.isnan(Q11).any():
            Q11 = fill_nans_with_interpolation(Q11)
        if np.isnan(Q12).any():
            Q12 = fill_nans_with_interpolation(Q12)
        if np.isnan(Q21).any():
            Q21 = fill_nans_with_interpolation(Q21)
        if np.isnan(Q22).any():
            Q22 = fill_nans_with_interpolation(Q22)

        values = [
            (Q11, (lat2 - lat) * (lon2 - lon)),
            (Q21, (lat - lat1) * (lon2 - lon)),
            (Q12, (lat2 - lat) * (lon - lon1)),
            (Q22, (lat - lat1) * (lon - lon1))
        ]

        if values:
            # Calculate the weighted sum
            interp_value = sum(val * weight for val, weight in values) / \
                sum(weight for _, weight in values)
        else:
            # All values are NaN, handle accordingly
            interp_value = np.nan  # or set to a default value, e.g., 0

        return interp_value

    print(".")
    # %% Extract information from mdf file

    string1 = 'Tstart'
    string2 = 'Tstop'
    string3 = 'Itdate'  # reference time

    file1 = open(mdf_file, "r")
    for line in file1:
        # checking string is present in line or not
        if string1 in line:
            values = line.split('=')
            tstart_val = values[1].strip()
            break
            file1.close()  # close file

    file1 = open(mdf_file, "r")
    for line in file1:
        # checking string is present in line or not
        if string2 in line:
            values_2 = line.split('=')
            tstop_val = values_2[1].strip()
            break
            file1.close()  # close file

    file1 = open(mdf_file, "r")
    for line in file1:
        # checking string is present in line or not
        if string3 in line:
            values_3 = line.split('=')
            ref_time_unedited = values_3[1].strip()
            break
            file1.close()  # close file

    # time series
    start = float(tstart_val)  # from mdf file
    stop = float(tstop_val)  # from mdf file
    ref_time = ref_time_unedited[1:11]
    # remove the hyphen for the bct file format
    reference_time = ref_time.replace('-', '')
    # step = 2.0000000e+001  # 20 minute step # adding here for understanding
    print(".")

    # %% extract start time and end time from mdf

    time_start = ref_time+" 00:00:00"  # Assuming it always starts at 00
    date_format_str = "%Y-%m-%d %H:%M:%S"

    # Calculate number of hours between ref time and sim time
    start_time_steps = int(start/60)  # to convert minutes to hours
    end_time_steps = int(stop/60)

    # create datetime object from timestamp string
    extracted_time = datetime.strptime(time_start, date_format_str)

    start_time = extracted_time + timedelta(hours=start_time_steps)
    # Convert datetime object to string in specific format
    start_time = start_time .strftime("%Y-%m-%d %H:%M:%S")
    start_time_og = start_time

    end_time = extracted_time + timedelta(hours=end_time_steps)
    end_time = end_time .strftime("%Y-%m-%d %H:%M:%S")
    end_time_og = end_time
    print(".")
    # %% correcting for 12 hour time difference in gsh

    time_start = start_time
    time_end = end_time
    date_format_str = "%Y-%m-%d %H:%M:%S"

    # create datetime object from timestamp string
    extracted_time = datetime.strptime(time_start, date_format_str)
    n = -12

    start_time = extracted_time + timedelta(hours=n)

    # Convert datetime object to string in specific format
    start_time = start_time .strftime("%Y-%m-%d %H:%M:%S")

    ##
    # create datetime object from timestamp string
    extracted_time = datetime.strptime(time_end, date_format_str)

    end_time = extracted_time + timedelta(hours=n)

    # Convert datetime object to string in specific format
    end_time = end_time .strftime("%Y-%m-%d %H:%M:%S")
    print(".")

    # %% Configuring time step to adhere to the coupling interval
    one_time_step = step
    min_data_time_step = int(20)

    if one_time_step <= 20:
        time_step_data = int(1)
    elif one_time_step > 20:
        time_step_data = int(one_time_step / min_data_time_step)
    print(".")
    # %% Open input files

    bnd_loc = pd.read_csv(
        boundaries, names=['boundary', 'easting', 'northing'], )

    nc_files = nc_files
    # Ensure the nc_files list is ordered by year
    nc_files.sort(key=extract_year)

    # Combine datasets from all .nc files
    datasets = []  # Store processed datasets
    end_time_previous = start_time  # Initialize the first dataset's start time
    
    for nc_file in nc_files:
        data = xr.open_dataset(nc_file, engine='h5netcdf')
        # start_data_time = data.nMesh2_data_time[0]
        end_data_time = data.nMesh2_data_time[-1]
    
        # Select the time range for the current dataset
        dataset = data.sel(nMesh2_data_time=slice(end_time_previous, end_data_time, time_step_data))
        datasets.append(dataset)  # Add to the list of datasets
    
        # Update the end_time_previous for the next dataset
        end_time_previous = end_data_time


    print(".")
    # %% Convert to geographic coordinates

    bnd_loc_geo = convert_coordinates(bnd_loc, crs_type)
    print(".")

    # %%  Converting time to scientific notations and check for records in table

    # create a range of the input time
    float_range = np.arange(start, stop + step, step).tolist()

    # for the bct file number of records in a tbale
    record_in_table = len(float_range)

    time_list = convert_list_to_sci_not(input_list=float_range,
                                        prec=7, exp_digits=3)
    print(".")

    # %%  Separating all boundary values into a dictionary with the boundary name as key
    output_dict = {}
    output_dict_2 = {}

    # Iterate through each boundary location
    for index, row in tqdm(bnd_loc_geo.iterrows(), desc='Extracting water level', total=len(bnd_loc_geo), leave=True, mininterval=0.1):
        
        water_level_combined = []  # List to store water levels from all datasets
    
        # Interpolate water levels from all datasets
        for dataset in datasets:  # `datasets` contains the preprocessed datasets from the .nc files
    
            water_level = bilinear_interpolation(
                row['lat'], row['lon'], dataset, variable='Mesh2_face_Wasserstand_2d'
            )
    
            water_level_combined.append(water_level)
    
        # Concatenate water levels across all datasets for the current boundary
        wl_combined = np.concatenate(water_level_combined)
    
        # Store the combined water levels in the output dictionary
        output_dict[row['boundaries']] = wl_combined
        

    for bnd_point_key, bnd_point_wl_list in output_dict.items():
        # Generate dict where keyword is bnd name and value is a list containing 3 list:
        # (1) time (2) wl at point a, (3) wl at point b

        bnd_point_wl_list_sci_not = convert_list_to_sci_not(
            input_list=bnd_point_wl_list,
            prec=7, exp_digits=3)  # convert to scientific notation

        # For every 'a' point: (1) create new list to write into dict (2) append array of time values
        if bnd_point_key[-1:] == 'a':
            list_in_dict = []
            list_in_dict.append(time_list)  # append array of time values
        list_in_dict.append(bnd_point_wl_list_sci_not)
        output_dict_2[bnd_point_key[:-2]] = list_in_dict

    print(".")
    print(".")
    print("Water level dataset extracted")
    print(".")
    print('Writing file')
    print(".")
    print(".")
    print(".")
    # %% write the bct file

    try:
        os.remove(bct_file_name)
    except FileNotFoundError:
        pass
    section_number = 1
    for key in output_dict_2:
        bn_name = str(key)
        header_lines = ["table-name           'Boundary Section : {}'".format(section_number),
                        "contents             'Uniform             '",
                        "location             '{}              '".format(
                            bn_name),
                        "time-function        'non-equidistant'",
                        "reference-time       {}".format(reference_time),
                        "time-unit            'minutes'",
                        "interpolation        'linear'",
                        "parameter            'time                '                     unit '[min]'",
                        "parameter            'water elevation (z)  end A'               unit '[m]'",
                        "parameter            'water elevation (z)  end B'               unit '[m]'",
                        "records-in-table     {}".format(record_in_table)]

        with open(bct_file_name, 'a', newline='') as f:
            for one_line in header_lines:
                f.write(one_line)
                f.write('\r\n')

            csv_writer = csv.writer(f)
            count = 0
            bnd_data_list = output_dict_2[key]
            for row in bnd_data_list[0]:
                # Set values to write, add leading blank for positive values
                time_val = add_blank_pos_val(input_str=bnd_data_list[0][count])
                wl_a = add_blank_pos_val(input_str=bnd_data_list[1][count])
                wl_b = add_blank_pos_val(input_str=bnd_data_list[2][count])

                # Generate row content as single string
                row_str = f'{time_val} {wl_a} {wl_b}'

                # Write row to file
                csv_writer.writerow([row_str])

                count += 1
        # Increment section_number for the next iteration
        section_number += 1
