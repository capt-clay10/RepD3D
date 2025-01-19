import warnings
from datetime import timedelta
from datetime import datetime
import sea_level_change
import mdw_writer
import bcw_year_overlap_file_generator
import bcw_generator
import bct_year_overlap_file_generator
import bct_generator
import time
import os
import extract_from_d3d_files
import output_methods
import rep_period
import plot_windroses
import cosmo_wind_file_generator
import ast
"""RUN THIS FILE"""
# %% import modules
warnings.filterwarnings("ignore", category=DeprecationWarning)


if __name__ == '__main__':

    t = time.time()  # start the time counter

    # %% Input data process 1
    print("""
          -----------------------------------------------------------------
          """)
          
    
    print("""
          RepD3D - toolbox for representative period identification and 
          associated boundary condition extraction
          
          Created by Clayton Soares
          Co-developer : Arne Knies
          Source: https://github.com/capt-clay10/RepD3D.git
          """)
    
    
    print("""
          -----------------------------------------------------------------
          """)
    
    
    
    path_req = input(
        '\n\nEnter the input/output path here (w/o quotation marks), eg E:/extract_bc : ')
    path = path_req  # 'F:/test'
    os.chdir(path)

    print('.')
    print('.')
    print('.')
    
    print("""
          --------------------------------NOTE-----------------------------
          """)
    print("Please read carefully the input criteria",
          "and choose which file you would like")
    print()
    print('Types of files offered:',
          'For all files, type 1',
          'For bct file overlapping over years, type 2',
          'For bcw file overlapping over years, type 3',
          'For boundary location csv file, type 4',
          'For boundary location and mdw file, type 5',
          'For adding sea level change to .bct files, type 6',
          'For identifying Representative period, type 7',
          'For generating COSMO wind field files, type 8\n', sep='\n')
    print('Important information : for choice 1 , boundary files cannot be generated with overlapping input years')
    
    print("""
          --------------------------------NOTE-----------------------------
          """)
    print()
    req = input("Enter number : ")
    print('.')
    print('.')
    print('.')
    choice = float(req)
    # %% CHOICE 1 extract all
    if choice == 1:
        # load time extracting function
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
                    # else:
                        # print('{} is not in the file'.format(string_name))
            return string_val
        
        print("""
              --------------------------------NOTE-----------------------------
              """)

        print("Please read carefully the input criteria,",
              "some requests are for the wave grid while some for flow\nPlease use file extensions for Delft3D files")
        
        print("""
              --------------------------------NOTE-----------------------------
              """)

        grid_req = input('Enter name of the Flow grid file : ')
        grid_input = grid_req

        bnd_req = input('Enter name of the Flow bnd file : ')
        bnd_input = bnd_req

        grid_wave_req = input('Enter name of the Wave grid file : ')
        grid_wave_input = grid_wave_req

        bnd_wave_req = input('Enter name of the Wave bnd file : ')
        bnd_wave_input = bnd_wave_req

        # Prompt user to enter years as a list
        years = input("Enter the flow year (e.g., 2002): ")
    
        # Generate the NetCDF file names using the naming pattern
        nc_file = f"{years}_1000m_waterlevel_2D.nc" 

        mdf_file_req = input('Enter the mdf file name : ')
        mdf_file = mdf_file_req  # 'test.mdf'

        # Prompt user to enter years as a list
        years_w = input("Enter the wave year (e.g., 2002): ")
   
        # Generate the NetCDF file names using the naming pattern
        nc_file_wave = f"{years_w}_1000m_wave_2D.nc" 

        mdw_file_req = input('Enter the mdw file name : ')
        mdw_file = mdw_file_req  # 'test.mdw'

        # request for time step
        step_req = input(
            'Enter time step to extract water level data (max resolution is 20 mins) format- 20 : ')
        # 2.0000000e+001  # 20 minute step # max resolution for gsh data
        step = float(step_req)

        step_wave_req = input(
            'Enter time step to extract WAVE data (max resolution is 20 mins, should be multiples of 20) format- 20 : ')
        step_wave = float(step_wave_req)

        # Extract start and end time from mdf file
        string1 = 'Tstart'
        tstart_val = value_from_txt_file(file=mdf_file, string_name=string1)
        string2 = 'Tstop'
        tstop_val = value_from_txt_file(file=mdf_file, string_name=string2)
        string3 = 'Itdate'  # reference time
        ref_time_unedited = value_from_txt_file(
            file=mdf_file, string_name=string3)
        start = float(tstart_val)  # from mdf file
        stop = float(tstop_val)  # from mdf file
        ref_time = ref_time_unedited[1:11]
        # remove the hyphen for the bct file format
        reference_time = ref_time.replace('-', '')
        # extract start time and end time from mdf
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
        end_time = extracted_time + timedelta(hours=end_time_steps)
        end_time = end_time .strftime("%Y-%m-%d %H:%M:%S")
        print('.')
        print('.')
        print('.')
        print('.')
        print("1 of 6")
        # output files
        # Use mdf file to extract bct file name and output file
        name_with_dot = mdf_file.partition('.')
        name_until_dot = name_with_dot[0]
        bct_file_name = '{}.bct'.format(name_until_dot)
        path_out_file = '{}.csv'.format(name_until_dot)

        # Use mdw file to extract bcw file name and output file
        wave_name_with_dot = mdw_file.partition('.')
        wave_name_until_dot = wave_name_with_dot[0]
        bcw_file = '{}.bcw'.format(wave_name_until_dot)
        wave_path_out_file = '{}.csv'.format(wave_name_until_dot)
        print('.')
        print("2 of 6")
        print('.')
        print('.')
        # Create the csv file for flow boundaries
        bnd_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
            path_bnd=bnd_input)

        coord_from_d3d_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
            path_grd=grid_input,
            request_list=bnd_grd_indices_output)

        crs_type_flow = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
        
        output_methods.write_bnd_coord_ascii(
            bnd_data_list=coord_from_d3d_grd_output, out_path=path_out_file)
        print('.')
        print('.')
        print('.')
        print(f'Grid type for flow is {crs_type_flow}')
        print('The process of creating',
              'the boundary location csv file for flow is completed - 3 of 6')
        print('.')
        print('.')
        print('.')

        # Create the csv file for wave boundaries
        bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
            path_bnd=bnd_wave_input)

        coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
            path_grd=grid_wave_input,
            request_list=bnd_wave_grd_indices_output)

        crs_type_wave = extract_from_d3d_files.extract_crs_type(path_grd=grid_wave_input)
        
        output_methods.write_bnd_coord_ascii(
            bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
        print('.')
        print('.')
        print('.')
        print(f'Grid type for wave is {crs_type_wave}')
        print('The process of creating',
              ' the boundary location csv file for Wave is completed - 4 of 6')

        # Create the bct file
        boundaries = path_out_file  # the csv file generated from process one
        bct = bct_generator.bct_file_generator(
            boundaries=boundaries, nc_file=nc_file, mdf_file=mdf_file, step=step,
            bct_file_name=bct_file_name, crs_type=crs_type_flow)

        t_2 = time.time()

        # end the time counter
        print('.')
        print('The process of extracting water level has now completed in : ')
        elapsed = time.time() - t
        print(str(elapsed) + " sec - 5 of 6")

        # Create the bcw file
        boundaries_wave = wave_path_out_file
        bcw = bcw_generator.bcw_file_generator(
            boundaries_wave=boundaries_wave, nc_file_wave=nc_file_wave, mdw_file=mdw_file, start_time=start_time,
            end_time=end_time, step_wave=step_wave,
            bcw_file_name=bcw_file, crs_type=crs_type_wave)

        print('.')
        print('The process of extracting wave boundary conditions has now completed in : ')
        elapsed = time.time() - t_2
        print(str(elapsed) + " sec - 6 of 6")
        print('.')
        print('.')
        print('.')
        # Write the new mdw file
        mdw_writer.write_mdw_file(
            mdw_file=mdw_file, boundaries_wave=boundaries_wave)
        print('New mdw file created')
        print('.')
        elapsed_final = time.time() - t
        print(f'Total time taken for both files is {elapsed_final/60} mins')

     # %% CHOICE 2 for bct overlapping years
    elif choice == 2:
        grid_req = input('Enter name of the Flow grid file : ')
        grid_input = grid_req

        bnd_req = input('Enter name of the Flow bnd file : ')
        bnd_input = bnd_req

        # Prompt user to enter years as a list
        years = input("Enter the list of years (e.g., ['2002', '2003']): ")
        years = eval(years)  # Convert string input to a Python list
    
        # Generate the NetCDF file names using the naming pattern
        nc_files = [f"{year}_1000m_waterlevel_2D.nc" for year in years]

        mdf_file_req = input('Enter the mdf file name : ')
        mdf_file = mdf_file_req  # 'test.mdf'

        step_req = input(
            'Enter time step to extract data (max resolution is 20mins) format 20 : ')
        step = float(step_req)  # 20 minute step # max resolution for gsh data

        # output files
        # Use mdf file to extract bct file and output file
        name_with_dot = mdf_file.partition('.')
        name_until_dot = name_with_dot[0]
        bct_file_name = '{}.bct'.format(name_until_dot)
        path_out_file = '{}.csv'.format(name_until_dot)
        print('.')
        print("1 of 3")

        # Create the csv file for flow boundaries
        bnd_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
            path_bnd=bnd_input)

        coord_from_d3d_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
            path_grd=grid_input,
            request_list=bnd_grd_indices_output)
        
        crs_type_flow = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)

        output_methods.write_bnd_coord_ascii(
            bnd_data_list=coord_from_d3d_grd_output, out_path=path_out_file)
        print('.')
        print('.')
        print('.')
        print(f'Grid type is {crs_type_flow}')
        print('The process of creating',
              'the boundary location csv file for flow is completed - 2 of 3')

        # Create the bct file
        boundaries = path_out_file  # the csv file generated from process one
        bct = bct_year_overlap_file_generator.bct_year_overlap_file_generator(
            boundaries=boundaries, nc_files=nc_files, mdf_file=mdf_file, step=step,
            bct_file_name=bct_file_name, crs_type=crs_type_flow)

        # end the time counter
        print('.')
        print('The process of extracting water level has now completed in : ')
        elapsed = time.time() - t
        print(str(elapsed) + " sec - 3 of 3")

    # %% Choice 3 for overlapping bcw file
    elif choice == 3:
        # BCW section input files
        grid_wave_req = input('Enter name of the Wave grid file : ')
        grid_wave_input = grid_wave_req

        bnd_wave_req = input('Enter name of the Wave bnd file : ')
        bnd_wave_input = bnd_wave_req

        # Prompt user to enter years as a list
        years = input("Enter the list of years (e.g., ['2002', '2003']): ")
        years = eval(years)  # Convert string input to a Python list
    
        # Generate the NetCDF file names using the naming pattern
        nc_files_wave = [f"{year}_1000m_wave_2D.nc" for year in years]

        mdw_file_req = input('Enter the mdw file name : ')
        mdw_file = mdw_file_req  # 'test.mdw'

        start_time_req = input(
            'Enter the simulation start time in the format YYYY-MM-DD hh:mm:ss : ')
        start_time = start_time_req  # '2015-02-01 00:00:00'

        end_time_req = input(
            'Enter the simulation end time in the format YYYY-MM-DD hh:mm:ss : ')
        end_time = end_time_req  # '2015-03-14 00:00:00'

        step_wave_req = input(
            'Enter time step to extract WAVE data (max resolution is 20mins, should be multiples of 20) format 20 : ')

        # 2.0000000e+001  # 20 minute step # max resolution for gsh data
        step_wave = float(step_wave_req)

        # bcw output file
        # Use mdw file to extract bcw file and output file
        wave_name_with_dot = mdw_file.partition('.')
        wave_name_until_dot = wave_name_with_dot[0]
        bcw_file = '{}.bcw'.format(wave_name_until_dot)
        wave_path_out_file = '{}.csv'.format(wave_name_until_dot)
        print('.')
        print("1 of 3")

        # Create the csv file for wave boundaries
        bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
            path_bnd=bnd_wave_input)

        coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
            path_grd=grid_wave_input,
            request_list=bnd_wave_grd_indices_output)

        crs_type_wave = extract_from_d3d_files.extract_crs_type(path_grd=grid_wave_input)
        
        output_methods.write_bnd_coord_ascii(
            bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
        print('.')
        print('.')
        print('.')
        print(f'Grid type is {crs_type_wave}')
        print('The process of creating',
              ' the boundary location csv file for wave is completed - 2 of 3')
        print('.')
        print('.')
        print('.')
        print('Initiating wave parameter extraction')

        # Create the bcw file
        boundaries_wave = wave_path_out_file
        bcw = bcw_year_overlap_file_generator.bcw_year_overlap_file_generator(
            boundaries_wave=boundaries_wave, nc_files_wave =nc_file_wave,
            mdw_file=mdw_file, start_time=start_time,
            end_time=end_time, step_wave=step_wave,
            bcw_file_name=bcw_file, crs_type=crs_type_wave)
        print('.')
        print(
            'The process of extracting wave boundary conditions has now completed in : ')
        elapsed = time.time() - t
        print(str(elapsed) + " sec - 3 of 3")
        print('.')
        print('.')
        # Write the new mdw file
        mdw_writer.write_mdw_file(
            mdw_file=mdw_file, boundaries_wave=boundaries_wave)
        print('New mdw file created')
    # %% CHOICE 4 make boundary csv file
    elif choice == 4:
        grid_req = input('Enter name of the grid file : ')
        grid_input = grid_req

        bnd_req = input('Enter name of the bnd file : ')
        bnd_input = bnd_req

        # Use grid file to extract bct file and output file
        wave_name_with_dot = grid_input.partition('.')
        wave_name_until_dot = wave_name_with_dot[0]
        wave_path_out_file = '{}.csv'.format(wave_name_until_dot)
        print('.')
        print("1 of 2")

        #  Create the csv file for wave boundaries
        bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
            path_bnd=bnd_input)

        coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
            path_grd=grid_input,
            request_list=bnd_wave_grd_indices_output)

        crs_type = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
        
        output_methods.write_bnd_coord_ascii(
            bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
        print('.')
        print('.')
        print('.')
        print(f'Grid type is {crs_type}')
        print('.')
        print('.')
        print('The process of creating',
              ' the boundary location csv file is completed - 2 of 2')

    #  Write boundary csv and mdw file
    elif choice == 5:
        grid_req = input('Enter name of the grid file : ')
        grid_input = grid_req

        bnd_req = input('Enter name of the bnd file : ')
        bnd_input = bnd_req

        mdw_file_req = input('Enter the mdw file name : ')
        mdw_file = mdw_file_req  # 'test.mdw'

        # Use grid file to extract bct file and output file
        wave_name_with_dot = grid_input.partition('.')
        wave_name_until_dot = wave_name_with_dot[0]
        wave_path_out_file = '{}.csv'.format(wave_name_until_dot)
        print('.')
        print("1 of 3")

        #  Create the csv file for wave boundaries
        bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
            path_bnd=bnd_input)

        coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
            path_grd=grid_input,
            request_list=bnd_wave_grd_indices_output)

        crs_type = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
        
        output_methods.write_bnd_coord_ascii(
            bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
        print('.')
        print('.')
        print('.')
        print(f'Grid type is {crs_type}')
        print('.')
        print('.')
        print('The process of creating',
              ' the boundary location csv file is completed - 2 of 3')

        # create mdw file
        boundaries_wave = wave_path_out_file
        mdw_writer.write_mdw_file(
            mdw_file=mdw_file, boundaries_wave=boundaries_wave)
        print('New mdw file created')
#%% Add sea level
    elif choice == 6:
        bct_req = input('Enter name of the bct file : ')
        bct_file_name = bct_req

        cons_grad = input('Enter 1 for gradual and 2 for constant increase :')
        cons_grad = int(cons_grad)
        if cons_grad == 1:
            type_inc = False
            type_add = 'gradual'
        elif cons_grad == 2:
            type_inc = True
            type_add = 'constant'
        else:
            print("Please only type 1 or 2")

        change_req = input('Enter sea level increase in m ( eg: 1) : ')
        change_amount = float(change_req)

        sea_level_change.add_wl(
            bct_file_name, sea_level_change=change_amount, constant=type_inc)

        print(f'sea level change added to water levels, type {type_add}')
# %% Representative period
    elif choice == 7:

        def parse_list_input(input_str):
            # Parse the input string as a Python expression (list)
            converted_list = ast.literal_eval(input_str)
            if isinstance(converted_list, list):
                return converted_list
        
        print("""
              --------------------------------NOTE-----------------------------
              """)

        print("Important note, the file should have only 3 columns, in the following order date, speed, dir")
        print("Important note, the file should be comma separated (,)")
        print("First row should be column headers")

        print("""
              --------------------------------NOTE-----------------------------
              """)
              
        file_input = input(
            'Enter file name of the wind data (w/o quotation marks), eg helgoland_wind.txt/csv : ')
        print('\n\n')

        out_file_name = input(
            'Enter desired output file name : ')
        
        print("""
              --------------------------------NOTE-----------------------------
              """)

        print("\n\nIn the next step you need to provide a list of relevant parameters, for how to prescribe them: ",
              "Directional orientations can be prescribed as ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE','SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW']",

              "Speed classes can be prescribed as [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] ",
              "Moving window duration can be prescribed as ['2MS', '3MS', '4MS', '5MS', '6MS']\n\n", sep='\n\n')
        
        print("""
              --------------------------------NOTE-----------------------------
              """)

        quad = input(
            "\nProvide list of directional orientations of relavance as mentioned above (minimum 2) :  ")
        quad = parse_list_input(quad)
        print()
        spd = input(
            "Provide list of speed classes of relavance as mentioned above (minimum 2) :  ")
        spd = parse_list_input(spd)

        start_time_total = input(
            "\nInsert start time of total reference period eg  : YYYY-MM-DD hh:mm:ss :")
        end_time_total = input(
            "Insert end time of total reference period eg  : YYYY-MM-DD hh:mm:ss :")

        frequency = input(
            "\nChoose moving window duration of relavance as mentioned above :  ")
        frequency = parse_list_input(frequency)

        print("\nNow insert details for the scanning period for prospect periods\n")
        start_time = input(
            "Insert start time of scanning window period eg  : YYYY-MM-DD hh:mm:ss :")

        print("\nEnd time must be (Total period end time - largest period frequency)\n")
        end_time = input(
            "Insert end time of scanning period eg  : YYYY-MM-DD hh:mm:ss :")

        out_path = input(
            "Insert output path for file saving eg : D:/extract_bc  :")

        output_file = rep_period.identify_rep_period(file_input, out_file_name,
                                                     quad, spd, start_time_total,
                                                     end_time_total,
                                                     frequency, start_time,
                                                     end_time, out_path)

        print('Representative wind file has been generated')

        plot_windroses.plot_windroses(file_input, start_time_total, end_time_total,
                                      output_file, parent=False)
#%% COSMO wind file generator
    elif choice == 8:
        
        print("""
              --------------------------------NOTE-----------------------------
              """)
              
        print("""
                \nThe main path containing the COSMO files should be structured as such\n\nTWO folders in the main path\n\nFolder 1 should be named UV and should have all the U and V monthly cosmo files you wish to extract from.
                \nFolder 2 should be named PS and should have all the PS monthly files.
                \nMake sure auxilary file cosmo_latlon_rot.grb is in the main path
                \nThe COSMO files can be found at:
                \nhttps://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/ 
                \nOn the webpage look for PS, U_10M and V_10M and download all monthly files necesssary and unzip them - use 7-Zip. Delete the zip files before generating the wind field files.
                \n In the current RepD3D version the components are rotated to true North
                """)

        print('\n\nPlease make sure the auxillary files, your equidistant 6X6 delft grid(.mat (V7)) and cosmo_latlon_rot.grb are in the main path as the cosmo folders UV and PS.')
        print('\nThey should be outside the UV and PS folders.')
        
        print("""
              --------------------------------NOTE-----------------------------
              """)

        cosmo_path_req = input(
            "\n\nMain path with COSMO (U,V,PS) data folders and aux mat files:")
        path = cosmo_path_req
        os.chdir(path)
        ref_time = input(
            """
            \nReference time: \nShould be same as your model reference date\n
            \nHowever it is better to have the real COSMO data start a little after the reference date, preferable download data from a month later
            \nso if your reference date is 2001-01-01 00:00:00 then cosmo data should start from 2001-02-01 00:00:00
            \n\neg format YYYY-MM-DD hh:mm:ss  :
            """)

        output_file_name = input(
            "\n\nOutput file name \neg: cosmo_2011 :")

        db_file = input(
            "\n\nEquidistant Delft3D grid (6km resolution) as .mat file \neg: DB_6km.mat :")
        
        crs_type = input(
            "\n\nInsert Coordinate system of Delft3D grid \neg: Cartesian or Spherical :")
        
        # cosmo correction file
        const_path = f'{cosmo_path_req}/cosmo_latlon_rot.grb'

        # Delete .idx files before starting the process
        for root, dirs, files in os.walk(cosmo_path_req):
            for file in files:
                if file.endswith('.idx'):
                    os.remove(os.path.join(root, file))

        cosmo_wind_file_generator.create_wind_fields_cosmo(db_file,
                                                           const_path,
                                                           cosmo_path_req,
                                                           output_file_name,
                                                           ref_time,crs_type)

    else:
        print("""
              --------------------------------NOTE-----------------------------
              """)
        print("You probably din't insert the number right, Please run again! ")
        
        print("""
              --------------------------------NOTE-----------------------------
              """)
