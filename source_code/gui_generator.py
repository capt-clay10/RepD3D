import tkinter as tk
from tkinter import filedialog, messagebox, Toplevel
from ttkbootstrap import Style
from PIL import Image, ImageTk
import sys
import time
import threading
import queue
import numpy as np
import output_methods
import extract_from_d3d_files
import bct_generator
import bct_year_overlap_file_generator
import bcw_generator
import bcw_year_overlap_file_generator
import mdw_writer
import sea_level_change
import rep_period
import plot_windroses
import cosmo_wind_file_generator
import ast
import os
from datetime import datetime, timedelta


def show_splash(app, resource_path):
    splash = tk.Toplevel()
    splash.overrideredirect(True)

    # Load splash screen image
    splash_image = tk.PhotoImage(file=resource_path('repd3d_logo.png'))
    splash_label = tk.Label(splash, image=splash_image)
    splash_label.pack()

    # Center splash screen
    window_width = 360
    window_height = 360
    screen_width = splash.winfo_screenwidth()
    screen_height = splash.winfo_screenheight()
    position_top = int(screen_height / 2 - window_height / 2)
    position_right = int(screen_width / 2 - window_width / 2)
    splash.geometry(
        f'{window_width}x{window_height}+{position_right}+{position_top}')

    def close_splash():
        splash.destroy()
        app.deiconify()  # Show the main window

    # Close splash screen after 3 seconds
    splash.after(3500, close_splash)
    splash.mainloop()

# %% Define application
class Application(tk.Tk):
    def toggle_fullscreen(event=None):
        app.attributes("-fullscreen", True)

    def end_fullscreen(event=None):
        app.attributes("-fullscreen", False)

    def __init__(self):
        super().__init__()
        self.withdraw()  # Hide the main window initially
        self.title("RepD3D")
        self.geometry('500x575')
        self.resizable(0, 0)

        style = Style(theme="darkly")

        self.task_queue = queue.Queue()

        # Load logo
        self.img = tk.PhotoImage(
            file=self.resource_path('repd3d_logo.png'))
        self.iconphoto(True, self.img)

        self.setup_ui()
        self.after(0, lambda: show_splash(
            self, self.resource_path))  # Show splash screen

    def resource_path(self, relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        try:
            base_path = sys._MEIPASS
        except Exception:
            base_path = os.path.abspath(".")

        return os.path.join(base_path, relative_path)

    def setup_ui(self):
        # First Frame for browsing directory
        frame1 = tk.Frame(self)
        frame1.pack(padx=10, pady=10, fill=tk.X)

        path_label = tk.Label(frame1, text="Select the input/output path:")
        path_label.pack(side=tk.LEFT, padx=10)

        self.path_entry = tk.Entry(frame1, width=30)
        self.path_entry.pack(side=tk.LEFT, padx=10)

        browse_button = tk.Button(
            frame1, text="Browse", command=self.browse_path)
        browse_button.pack(side=tk.LEFT, pady=10)

        # Second Frame for choosing file type using radio buttons
        frame2 = tk.Frame(self, borderwidth=1, relief='solid')
        frame2.pack(padx=10, pady=10)

        file_types = [
            "Identify representative period",
            "Extract all files",
            "Bct file overlapping over years",
            "Bcw file overlapping over years",
            "Boundary location CSV file",
            "Boundary location and mdw file",
            "Add sea level change to .bct files",
            "Generate ST-varying wind field files- COSMO"
        ]

        self.choice_var = tk.IntVar()

        for idx, file_type in enumerate(file_types, start=1):
            radio_button = tk.Radiobutton(frame2, text=file_type,
                                          variable=self.choice_var, value=idx, font=('Times', 14))
            radio_button.pack(anchor=tk.W, pady=5)

        # Submit button
        submit_button = tk.Button(
            self, text="Submit", command=self.submit_choice, width=20)
        submit_button.pack(pady=10)

        # Information frame
        frame3 = tk.Frame(self)
        frame3.pack()


        frame4 = tk.Frame(self, borderwidth=0.2, relief='solid')
        frame4.pack(padx=10, pady=10, side='left')

        text = """
        \nGUI created by : Clayton Soares\ncontact: clayton.soares@ifg.uni-kiel.de\nCo-developer : Arne Knies
        \nSource: https://github.com/capt-clay10/RepD3D.git

        
        """
        permanent_text_label = tk.Label(
            frame4, text=text, font=('Helvetic', 10), justify=tk.LEFT, wraplength=450)
        permanent_text_label.pack(anchor='sw')

        # Start the queue processing
        self.after(100, self.process_queue)

    def browse_path(self):
        folder_selected = filedialog.askdirectory()
        self.path_entry.delete(0, tk.END)
        self.path_entry.insert(0, folder_selected)

    def run_task_in_thread(self, task, *args):
        def wrapper():
            result = task(*args)
            self.task_queue.put((task.__name__, result))
        thread = threading.Thread(target=wrapper)
        thread.start()

    def process_queue(self):
        try:
            while True:
                task_name, result = self.task_queue.get_nowait()
                if task_name == "long_running_task":
                    messagebox.showinfo("Info", result)
                elif task_name == "submit_choice":
                    self.handle_task_result(result)
        except queue.Empty:
            pass
        self.after(100, self.process_queue)

    def handle_task_result(self, result):
        # Process the result of the task here
        pass

    def submit_choice(self):
        selected_choice = self.choice_var.get()
        if not self.path_entry.get():
            messagebox.showwarning(
                "Warning", "Please browse for a directory before submitting.")
        elif selected_choice == 0:
            messagebox.showwarning(
                "Warning", "Please select a file type before submitting.")
        elif selected_choice in [1, 2, 3, 4, 5, 6, 7, 8]:
            t = time.time()  # start the time counter

            new_window = tk.Toplevel(self)
            new_window.title("Processing Window")
            # new_window.geometry('1000x1000')
            new_window.grab_set()
            new_window.resizable(1, 0)
            # Bind the F11 key to toggle full screen for the internal window
            new_window.bind("<F11>", self.toggle_fullscreen)
            # Bind the Escape key to exit full screen for the internal window
            new_window.bind("<Escape>", self.end_fullscreen)

            # Function to browse for files
            def browse_file(entry_widget):
                file_path = filedialog.askopenfilename()
                entry_widget.delete(0, "end")
                entry_widget.insert(0, file_path)

            main_frame = tk.Frame(new_window)
            main_frame.pack(fill='both', expand=True, padx=10, pady=10)

            # Console frame
            frame_console = tk.Frame(
                new_window, width=200, height=5, borderwidth=1, relief='solid')
            frame_console.pack(side='bottom', fill='both', padx=10, pady=10)

            # Create a Text widget to display console output
            console_output = tk.Text(
                frame_console, wrap=tk.WORD, width=100, height=5)
            console_output.pack(padx=10, pady=10)
            
            # Disclaimer label
            disclaimer_label = tk.Label(
                frame_console,
                text="When running it might look unresponsive but it's working, please be patient...",
                font=('Times', 10),
                fg='red',
            )
            disclaimer_label.pack(pady=5)

            # Redirect sys.stdout to the Text widget
            class ConsoleRedirector:
                def __init__(self, text_widget):
                    self.text_space = text_widget
                    self._stdout = sys.stdout
                    self._stderr = sys.stderr

                def write(self, message):
                    if not message.endswith('\n'):
                        message += '\n'
                    self.text_space.insert(tk.END, message)
                    # Automatically scroll to the end
                    self.text_space.see(tk.END)
                    self.text_space.update_idletasks()  # Update the widget

                def flush(self):
                    self.text_space.update_idletasks()  # Ensure the widget updates

            # console_redirector = ConsoleRedirector(console_output)
            # Redirect stdout to the Text widget
            sys.stdout = ConsoleRedirector(console_output)
            sys.stderr = ConsoleRedirector(console_output)

            def process_choice():
                if selected_choice == 1:
                    self.identify_representative_period(main_frame)

                elif selected_choice == 2:
                    self.process_all_files(main_frame)

                elif selected_choice == 3:
                    self.process_bct_overlap_files(main_frame)

                elif selected_choice == 4:
                    self.process_bcw_overlap_files(main_frame)

                elif selected_choice == 5:
                    self.generate_boundary_csv(main_frame)

                elif selected_choice == 6:
                    self.generate_boundary_mdw(main_frame)

                elif selected_choice == 7:
                    self.add_sea_level(main_frame)
                    
                elif selected_choice == 8:
                    self.generate_wind_files(main_frame)

            # Run the processing function in a separate thread to avoid blocking the GUI
            threading.Thread(target=process_choice).start()
#%% Extract all files
    def process_all_files(self, main_frame):

        # Frame for Flow Module
        frameup = tk.Frame(main_frame, width=250, borderwidth=1, relief='solid')
        frameup.pack(side='left', fill='both', expand=True, padx=10)
        left_label = tk.Label(frameup, text="Flow Module", font=("Times", 16))
        left_label.pack(side='top', padx=10, pady=10)

        # Frame for Wave Module
        framedown = tk.Frame(main_frame, width=250,
                             borderwidth=1, relief='solid')
        framedown.pack(side='right', fill='both', expand=True, padx=10)
        right_label = tk.Label(
            framedown, text="Wave Module", font=("Times", 16))
        right_label.pack(side='top', padx=10, pady=10)
        t = time.time()
        # Function to browse for files

        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)

        # Input buttons for Flow Module
        mdf_label = tk.Label(frameup, text="MDF file:",
                             font='Times').pack(pady=5)
        mdf_entry = tk.Entry(frameup, width=50)
        mdf_entry.pack()
        mdf_button = tk.Button(frameup, text="Browse",
                               command=lambda: browse_file(mdf_entry)).pack(pady=20)

        grd_label = tk.Label(frameup, text=".grd file (flow):",
                             font='Times').pack(pady=5)
        grd_entry = tk.Entry(frameup, width=50)
        grd_entry.pack()
        grd_button = tk.Button(frameup, text="Browse",
                               command=lambda: browse_file(grd_entry)).pack(pady=20)

        bnd_label = tk.Label(frameup, text=".bnd file (flow)\n boundary naming format\n Name_1 :",
                             font='Times').pack(pady=5)
        bnd_entry = tk.Entry(frameup, width=50)
        bnd_entry.pack()
        bnd_button = tk.Button(frameup, text="Browse",
                               command=lambda: browse_file(bnd_entry)).pack(pady=20)

        nc_label = tk.Label(
            frameup, text=".nc waterlevel file (EasyGSH):", font='Times').pack(pady=5)
        nc_entry = tk.Entry(frameup, width=50)
        nc_entry.pack()
        nc_button = tk.Button(frameup, text="Browse",
                              command=lambda: browse_file(nc_entry)).pack(pady=20)

        step_types = ["20 mins", "40 mins", "60 mins", "80 mins"]

        stepf_label = tk.Label(
            frameup, text="Time step for waterlevel extraction:", font='Times').pack(pady=5)
        selected_step_f = tk.IntVar()

        for idx, step_type in enumerate(step_types, start=1):
            radio_button_f = tk.Radiobutton(
                frameup, text=step_type, variable=selected_step_f, value=idx)
            radio_button_f.pack(anchor='c')

        # Input buttons for Wave Module
        mdw_label = tk.Label(framedown, text="MDW file:",
                             font='Times').pack(pady=5)
        mdw_entry = tk.Entry(framedown, width=50)
        mdw_entry.pack()
        mdw_button = tk.Button(framedown, text="Browse",
                               command=lambda: browse_file(mdw_entry)).pack(pady=20)

        grdw_label = tk.Label(
            framedown, text=".grd file (wave):", font='Times').pack(pady=5)
        grdw_entry = tk.Entry(framedown, width=50)
        grdw_entry.pack()
        grdw_button = tk.Button(framedown, text="Browse",
                                command=lambda: browse_file(grdw_entry)).pack(pady=20)

        bndw_label = tk.Label(
            framedown, text=".bnd file (wave)\n boundary naming format\n Name1 :", font='Times').pack(pady=5)
        bndw_entry = tk.Entry(framedown, width=50)
        bndw_entry.pack()
        bndw_button = tk.Button(framedown, text="Browse",
                                command=lambda: browse_file(bndw_entry)).pack(pady=20)

        ncw_label = tk.Label(
            framedown, text=".nc wave file (EasyGSH):", font='Times').pack(pady=5)
        ncw_entry = tk.Entry(framedown, width=50)
        ncw_entry.pack()
        ncw_button = tk.Button(framedown, text="Browse",
                               command=lambda: browse_file(ncw_entry)).pack(pady=20)

        stepw_label = tk.Label(
            framedown, text="Time step for wave extraction:", font='Times').pack(pady=5)
        step_types_w = ["20 mins", "40 mins", "60 mins", "80 mins", "120 mins"]
        selected_step_w = tk.IntVar()

        for idx_w, step_type_w in enumerate(step_types_w, start=1):
            radio_button_w = tk.Radiobutton(
                framedown, variable=selected_step_w, text=step_type_w, value=idx_w)
            radio_button_w.pack(anchor='c')

        def check_submit():
            t = time.time()
            if not mdf_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .mdf file")
            elif not grd_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .grd (flow) file")
            elif not bnd_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .bnd (flow) file")
            elif not nc_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .nc waterlevel file")
            elif not mdw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .mdw file")
            elif not grdw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .grd (wave) file.")
            elif not bndw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .bnd (wave) file.")
            elif not ncw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .nc wave file.")
            else:
                # Input files from browsed
                step_w = selected_step_w.get()
                step_f = selected_step_f.get()

                def step_number(step):
                    if step == 1:
                        step_c = 20
                    elif step == 2:
                        step_c = 40
                    elif step == 3:
                        step_c = 60
                    else:
                        step_c = 80
                    return step_c

                def step_number_w(step):
                    if step == 1:
                        step_c = 20
                    elif step == 2:
                        step_c = 40
                    elif step == 3:
                        step_c = 60
                    elif step == 4:
                        step_c = 80
                    else:
                        step_c = 120
                    return step_c

                step = step_number(step_f)
                step_wave = step_number_w(step_w)

                def value_from_txt_file(file, string_name):
                    file1 = open(file, "r")
                    for line in file1:
                        if '=' in line:
                            if string_name in line:
                                val = line.split('=')
                                string_val = val[1].strip()
                                break
                    file1.close()
                    return string_val

                # Files
                mdf_file = mdf_entry.get()
                mdw_file = mdw_entry.get()
                grid_input = grd_entry.get()
                bnd_input = bnd_entry.get()
                grid_wave_input = grdw_entry.get()
                bnd_wave_input = bndw_entry.get()
                nc_file = nc_entry.get()
                nc_file_wave = ncw_entry.get()

                # Extract start and end time from mdf file
                string1 = 'Tstart'
                tstart_val = value_from_txt_file(
                    file=mdf_file, string_name=string1)
                string2 = 'Tstop'
                tstop_val = value_from_txt_file(
                    file=mdf_file, string_name=string2)
                string3 = 'Itdate'  # reference time
                ref_time_unedited = value_from_txt_file(
                    file=mdf_file, string_name=string3)
                start = float(tstart_val)  # from mdf file
                stop = float(tstop_val)  # from mdf file
                ref_time = ref_time_unedited[1:11]
                # remove the hyphen for the bct file format
                reference_time = ref_time.replace('-', '')
                time_start = ref_time+" 00:00:00"  # Assuming it always starts at 00
                date_format_str = "%Y-%m-%d %H:%M:%S"
                start_time_steps = int(start/60)  # to convert minutes to hours
                end_time_steps = int(stop/60)
                extracted_time = datetime.strptime(time_start, date_format_str)
                start_time = extracted_time + timedelta(hours=start_time_steps)
                start_time = start_time .strftime("%Y-%m-%d %H:%M:%S")
                end_time = extracted_time + timedelta(hours=end_time_steps)
                end_time = end_time .strftime("%Y-%m-%d %H:%M:%S")

                # Output files
                name_with_dot = mdf_file.partition('.')
                name_until_dot = name_with_dot[0]
                bct_file_name = '{}.bct'.format(name_until_dot)
                path_out_file = '{}.csv'.format(name_until_dot)
                wave_name_with_dot = mdw_file.partition('.')
                wave_name_until_dot = wave_name_with_dot[0]
                bcw_file = '{}.bcw'.format(wave_name_until_dot)
                wave_path_out_file = '{}.csv'.format(wave_name_until_dot)

                # Create the csv file for flow boundaries
                bnd_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
                    path_bnd=bnd_input)
                coord_from_d3d_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
                    path_grd=grid_input, request_list=bnd_grd_indices_output)
                crs_type_flow = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
                output_methods.write_bnd_coord_ascii(
                    bnd_data_list=coord_from_d3d_grd_output, out_path=path_out_file)
                print(f'Grid type for flow is {crs_type_flow}')
                print('.')
                print('.')
                # Create the csv file for wave boundaries
                bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
                    path_bnd=bnd_wave_input)
                coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
                    path_grd=grid_wave_input, request_list=bnd_wave_grd_indices_output)
                crs_type_wave = extract_from_d3d_files.extract_crs_type(path_grd=grid_wave_input)
                output_methods.write_bnd_coord_ascii(
                    bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
                print(f'Grid type for wave is {crs_type_wave}')
                print('.')
                print('.')
                # Create the bct file
                boundaries = path_out_file  # the csv file generated from process one
                bct = bct_generator.bct_file_generator(
                    boundaries=boundaries, nc_file=nc_file, mdf_file=mdf_file, step=step,
                    bct_file_name=bct_file_name, crs_type=crs_type_flow)

                # Create the bcw file
                boundaries_wave = wave_path_out_file
                bcw = bcw_generator.bcw_file_generator(
                    boundaries_wave=boundaries_wave, nc_file_wave=nc_file_wave, mdw_file=mdw_file, start_time=start_time,
                    end_time=end_time, step_wave=step_wave,
                    bcw_file_name=bcw_file, crs_type=crs_type_wave)

                # Write the new mdw file
                mdw_writer.write_mdw_file(
                    mdw_file=mdw_file, boundaries_wave=boundaries_wave)

                elapsed_final = (time.time() - t) / 60
                print(
                    f'Total time taken for both files is {np.round(elapsed_final, 2)} mins')

        # Submit code for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(side='bottom', fill='both', padx=10, pady=10)
        
        # Load and display the image
        image_path = self.resource_path('easygsh_extent.jpg')
        img = Image.open(image_path)
        img = img.resize((400, 400), Image.LANCZOS)  # Resize the image if necessary
        img = ImageTk.PhotoImage(img)
        image_label = tk.Label(frame_submit, image=img)
        image_label.image = img  # Keep a reference to avoid garbage collection
        image_label.pack(padx=10, pady=10)
        
        # Citation text 
        text_2 = """Source for EasyGSH data: https://mdi-de.baw.de/easygsh/Easy_Viewer_syn.html#home\n\nCitations for using EasyGSH data:\n Hagen, R., Plüß, A., Schrage, N., Dreier, N. (2020): EasyGSH-DB: Themengebiet - synoptische Hydrodynamik. Bundesanstalt für Wasserbau. https://doi.org/10.48437/02.2020.K2.7000.0004\n"""

        permanent_text_label_2 = tk.Label(
            frame_submit, text=text_2, justify=tk.LEFT, wraplength=400, font=('Times', 11))
        permanent_text_label_2.pack()

        
        submit_button = tk.Button(
            frame_submit, text="Extract Boundary conditions", command=check_submit, width=30)
        submit_button.pack(pady=10, anchor='s')
    # def process_bct_files(self, main_frame):

    #     # Frame for Flow Module
    #     frameup = tk.Frame(main_frame, width=250, borderwidth=1, relief='solid')
    #     frameup.pack(fill='both', expand=True, padx=10)
    #     left_label = tk.Label(frameup, text="Flow Module", font=("Times", 16))
    #     left_label.pack(side='top', padx=10, pady=10)

    #     # Function to browse for files
    #     def browse_file(entry_widget):
    #         file_path = filedialog.askopenfilename()
    #         entry_widget.delete(0, "end")
    #         entry_widget.insert(0, file_path)

    #     # Input buttons for Flow Module
    #     mdf_label = tk.Label(frameup, text="MDF file:",
    #                          font='Times').pack(pady=5)
    #     mdf_entry = tk.Entry(frameup, width=50)
    #     mdf_entry.pack()
    #     mdf_button = tk.Button(
    #         frameup, text="Browse", command=lambda: browse_file(mdf_entry)).pack(pady=20)

    #     grd_label = tk.Label(
    #         frameup, text=".grd file (flow):", font='Times').pack(pady=5)
    #     grd_entry = tk.Entry(frameup, width=50)
    #     grd_entry.pack()
    #     grd_button = tk.Button(
    #         frameup, text="Browse", command=lambda: browse_file(grd_entry)).pack(pady=20)

    #     bnd_label = tk.Label(
    #         frameup, text=".bnd file (flow)\n boundary naming format\n Name_1 :", font='Times').pack(pady=5)
    #     bnd_entry = tk.Entry(frameup, width=50)
    #     bnd_entry.pack()
    #     bnd_button = tk.Button(
    #         frameup, text="Browse", command=lambda: browse_file(bnd_entry)).pack(pady=20)

    #     nc_label = tk.Label(
    #         frameup, text=".nc waterlevel file (EasyGSH):", font='Times').pack(pady=5)
    #     nc_entry = tk.Entry(frameup, width=50)
    #     nc_entry.pack()
    #     nc_button = tk.Button(frameup, text="Browse",
    #                           command=lambda: browse_file(nc_entry)).pack(pady=20)

    #     step_types = ["20 mins", "40 mins", "60 mins", "80 mins"]

    #     stepf_label = tk.Label(
    #         frameup, text="Time step for waterlevel extraction:", font='Times').pack(pady=5)
    #     selected_step_f = tk.IntVar()

    #     for idx, step_type in enumerate(step_types, start=1):
    #         radio_button_f = tk.Radiobutton(
    #             frameup, text=step_type, variable=selected_step_f, value=idx)
    #         radio_button_f.pack(anchor='c')

    #     def check_submit():
    #         t = time.time()
    #         if not mdf_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .mdf file")
    #         elif not grd_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .grd (flow) file")
    #         elif not bnd_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .bnd (flow) file")
    #         elif not nc_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .nc waterlevel file")
    #         else:
    #             # Input files from browsed
    #             step_f = selected_step_f.get()

    #             def step_number(step):
    #                 if step == 1:
    #                     step_c = 20
    #                 elif step == 2:
    #                     step_c = 40
    #                 elif step == 3:
    #                     step_c = 60
    #                 else:
    #                     step_c = 80
    #                 return step_c

    #             step = step_number(step_f)

    #             def value_from_txt_file(file, string_name):
    #                 file1 = open(file, "r")
    #                 for line in file1:
    #                     if '=' in line:
    #                         if string_name in line:
    #                             val = line.split('=')
    #                             string_val = val[1].strip()
    #                             break
    #                 file1.close()
    #                 return string_val

    #             # Files
    #             mdf_file = mdf_entry.get()
    #             grid_input = grd_entry.get()
    #             bnd_input = bnd_entry.get()
    #             nc_file = nc_entry.get()

    #             # Extract start and end time from mdf file
    #             string1 = 'Tstart'
    #             tstart_val = value_from_txt_file(
    #                 file=mdf_file, string_name=string1)
    #             string2 = 'Tstop'
    #             tstop_val = value_from_txt_file(
    #                 file=mdf_file, string_name=string2)
    #             string3 = 'Itdate'  # reference time
    #             ref_time_unedited = value_from_txt_file(
    #                 file=mdf_file, string_name=string3)
    #             start = float(tstart_val)  # from mdf file
    #             stop = float(tstop_val)  # from mdf file
    #             ref_time = ref_time_unedited[1:11]
    #             # remove the hyphen for the bct file format
    #             reference_time = ref_time.replace('-', '')
    #             time_start = ref_time+" 00:00:00"  # Assuming it always starts at 00
    #             date_format_str = "%Y-%m-%d %H:%M:%S"
    #             start_time_steps = int(start/60)  # to convert minutes to hours
    #             end_time_steps = int(stop/60)
    #             extracted_time = datetime.strptime(time_start, date_format_str)
    #             start_time = extracted_time + timedelta(hours=start_time_steps)
    #             start_time = start_time .strftime("%Y-%m-%d %H:%M:%S")
    #             end_time = extracted_time + timedelta(hours=end_time_steps)
    #             end_time = end_time .strftime("%Y-%m-%d %H:%M:%S")

    #             # Output files
    #             name_with_dot = mdf_file.partition('.')
    #             name_until_dot = name_with_dot[0]
    #             bct_file_name = '{}.bct'.format(name_until_dot)
    #             path_out_file = '{}.csv'.format(name_until_dot)

    #             # Create the csv file for flow boundaries
    #             bnd_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
    #                 path_bnd=bnd_input)
    #             coord_from_d3d_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
    #                 path_grd=grid_input, request_list=bnd_grd_indices_output)
    #             crs_type_flow = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
    #             output_methods.write_bnd_coord_ascii(
    #                 bnd_data_list=coord_from_d3d_grd_output, out_path=path_out_file)
    #             print(f'Grid type for flow is {crs_type_flow}')
    #             print('.')
    #             print('.')
    #             # Create the bct file
    #             boundaries = path_out_file  # the csv file generated from process one
    #             bct = bct_generator.bct_file_generator(
    #                 boundaries=boundaries, nc_file=nc_file, mdf_file=mdf_file, step=step,
    #                 bct_file_name=bct_file_name, crs_type=crs_type_flow)

    #             elapsed_final = (time.time() - t) / 60
    #             print(
    #                 f'Total time taken for bct file is {np.round(elapsed_final, 2)} mins')

    #     # Submit code for execution
    #     frame_submit = tk.Frame(main_frame, width=200,
    #                             height=20, borderwidth=1, relief='solid')
    #     frame_submit.pack(fill='both', padx=10, pady=10)

    #     submit_button = tk.Button(
    #         frame_submit, text="Extract waterlevel timeseries (bct)", command=check_submit, width=30)
    #     submit_button.pack(pady=10, anchor='s')
# %% BCT overlap file
    def process_bct_overlap_files(self, main_frame):
        # Frame for Flow Module
        frameup = tk.Frame(main_frame, borderwidth=1, relief='solid')
        frameup.pack(fill='both', expand=True, padx=10)
    
        # Left section for MDF, GRD, and BND files
        left_frame = tk.Frame(frameup, borderwidth=1, relief='solid', padx=10, pady=10)
        left_frame.pack(side='left', fill='both', expand=True)
    
        left_label = tk.Label(left_frame, text="Flow Module", font=("Times", 16))
        left_label.pack(side='top', padx=10, pady=10)
    
        # Function to browse for files
        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)
    
        # Input buttons for Flow Module
        mdf_label = tk.Label(left_frame, text="MDF file:", font='Times')
        mdf_label.pack(pady=5)
        mdf_entry = tk.Entry(left_frame, width=50)
        mdf_entry.pack()
        mdf_button = tk.Button(left_frame, text="Browse", command=lambda: browse_file(mdf_entry))
        mdf_button.pack(pady=10)
    
        grd_label = tk.Label(left_frame, text=".grd file (flow):", font='Times')
        grd_label.pack(pady=5)
        grd_entry = tk.Entry(left_frame, width=50)
        grd_entry.pack()
        grd_button = tk.Button(left_frame, text="Browse", command=lambda: browse_file(grd_entry))
        grd_button.pack(pady=10)
    
        bnd_label = tk.Label(left_frame, text=".bnd file (flow)\n boundary naming format\n Name_1 :", font='Times')
        bnd_label.pack(pady=5)
        bnd_entry = tk.Entry(left_frame, width=50)
        bnd_entry.pack()
        bnd_button = tk.Button(left_frame, text="Browse", command=lambda: browse_file(bnd_entry))
        bnd_button.pack(pady=10)
    
        # Right section for .nc files and Time Step options
        right_frame = tk.Frame(frameup, borderwidth=1, relief='solid', padx=10, pady=10)
        right_frame.pack(side='right', fill='both', expand=True)
    
        # Frame for displaying selected .nc files
        nc_files = []  # List to store selected .nc files
    
        def browse_nc_files():
            """Browse and select multiple .nc files."""
            file_paths = filedialog.askopenfilenames(filetypes=[("NetCDF Files", "*.nc")])
            if file_paths:
                nc_files.extend(file_paths)
                update_nc_files_display()
    
        def update_nc_files_display():
            """Update the display of selected .nc files."""
            for widget in nc_files_frame.winfo_children():
                widget.destroy()
    
            for idx, file in enumerate(nc_files, start=1):
                label = tk.Label(nc_files_frame, text=f"{idx}. {file}", font='Times', anchor='w', justify='left', wraplength=300)
                label.pack(fill='x', padx=5, pady=2)
    
        # Title for .nc files section
        nc_title = tk.Label(right_frame, text="Select all water level .nc Files", font=("Times", 14))
        nc_title.pack(anchor='n', pady=10)
    
        nc_files_frame = tk.Frame(right_frame, borderwidth=1, relief='sunken', width=300, height=200)
        nc_files_frame.pack(fill='both', padx=5, pady=5, expand=True)
    
        browse_nc_button = tk.Button(right_frame, text="Browse .nc Files", command=browse_nc_files)
        browse_nc_button.pack(pady=10)
    
        # Time Step Options
        stepf_label = tk.Label(right_frame, text="Time step for waterlevel extraction:", font='Times')
        stepf_label.pack(pady=10, anchor='w')
    
        selected_step_f = tk.IntVar(value=1)
        step_types = ["20 mins", "40 mins", "60 mins", "80 mins"]
        for idx, step_type in enumerate(step_types, start=1):
            radio_button_f = tk.Radiobutton(right_frame, text=step_type, variable=selected_step_f, value=idx)
            radio_button_f.pack(anchor='w')

        def check_submit():
            t = time.time()
            if not mdf_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .mdf file")
            elif not grd_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .grd (flow) file")
            elif not bnd_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .bnd (flow) file")
            elif not nc_files:
                messagebox.showwarning(
                    "Warning", "Please browse for the .nc waterlevel file")
            else:
                # Input files from browsed
                step_f = selected_step_f.get()

                def step_number(step):
                    if step == 1:
                        step_c = 20
                    elif step == 2:
                        step_c = 40
                    elif step == 3:
                        step_c = 60
                    else:
                        step_c = 80
                    return step_c

                step = step_number(step_f)

                def value_from_txt_file(file, string_name):
                    file1 = open(file, "r")
                    for line in file1:
                        if '=' in line:
                            if string_name in line:
                                val = line.split('=')
                                string_val = val[1].strip()
                                break
                    file1.close()
                    return string_val

                # Files
                mdf_file = mdf_entry.get()
                grid_input = grd_entry.get()
                bnd_input = bnd_entry.get()


                # Extract start and end time from mdf file
                string1 = 'Tstart'
                tstart_val = value_from_txt_file(
                    file=mdf_file, string_name=string1)
                string2 = 'Tstop'
                tstop_val = value_from_txt_file(
                    file=mdf_file, string_name=string2)
                string3 = 'Itdate'  # reference time
                ref_time_unedited = value_from_txt_file(
                    file=mdf_file, string_name=string3)
                start = float(tstart_val)  # from mdf file
                stop = float(tstop_val)  # from mdf file
                ref_time = ref_time_unedited[1:11]
                # remove the hyphen for the bct file format
                reference_time = ref_time.replace('-', '')
                time_start = ref_time+" 00:00:00"  # Assuming it always starts at 00
                date_format_str = "%Y-%m-%d %H:%M:%S"
                start_time_steps = int(start/60)  # to convert minutes to hours
                end_time_steps = int(stop/60)
                extracted_time = datetime.strptime(time_start, date_format_str)
                start_time = extracted_time + timedelta(hours=start_time_steps)
                start_time = start_time .strftime("%Y-%m-%d %H:%M:%S")
                end_time = extracted_time + timedelta(hours=end_time_steps)
                end_time = end_time .strftime("%Y-%m-%d %H:%M:%S")

                # Output files
                name_with_dot = mdf_file.partition('.')
                name_until_dot = name_with_dot[0]
                bct_file_name = '{}.bct'.format(name_until_dot)
                path_out_file = '{}.csv'.format(name_until_dot)

                # Create the csv file for flow boundaries
                bnd_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
                    path_bnd=bnd_input)
                coord_from_d3d_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
                    path_grd=grid_input, request_list=bnd_grd_indices_output)
                crs_type_flow = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
                output_methods.write_bnd_coord_ascii(
                    bnd_data_list=coord_from_d3d_grd_output, out_path=path_out_file)
                print(f'Grid type for flow is {crs_type_flow}')
                print('.')
                print('.')
                # Create the bct file
                boundaries = path_out_file  # the csv file generated from process one
                bct = bct_year_overlap_file_generator.bct_year_overlap_file_generator(
                    boundaries=boundaries, nc_files=nc_files,
                    mdf_file=mdf_file, step=step,
                    bct_file_name=bct_file_name, crs_type=crs_type_flow)

                elapsed_final = (time.time() - t) / 60
                print(
                    f'Total time taken for the bct overlap file is {np.round(elapsed_final, 2)} mins')

        # Submit button for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(side='bottom', fill='both', padx=10, pady=10)

        submit_button = tk.Button(
            frame_submit, text="Extract waterlevel timeseries (bct)", command=check_submit, width=30)
        submit_button.pack(pady=10, anchor='s')
    # def process_bcw_files(self, main_frame):
    #     # Frame for Wave Module
    #     framedown = tk.Frame(main_frame, width=250,
    #                          borderwidth=1, relief='solid')
    #     framedown.pack(side='left', fill='both', expand=True, padx=10)
    #     right_label = tk.Label(
    #         framedown, text="Wave Module", font=("Times", 16))
    #     right_label.pack(side='top', padx=10, pady=10)

    #     # Frame for right Module
    #     frameright = tk.Frame(main_frame, width=250,
    #                           borderwidth=1, relief='solid')
    #     frameright.pack(side='right', fill='both', expand=True, padx=10)

    #     t = time.time()
    #     # Function to browse for files

    #     def browse_file(entry_widget):
    #         file_path = filedialog.askopenfilename()
    #         entry_widget.delete(0, "end")
    #         entry_widget.insert(0, file_path)

    #     # Input buttons for Wave Module
    #     mdw_label = tk.Label(framedown, text="MDW file:",
    #                          font='Times').pack(pady=5)
    #     mdw_entry = tk.Entry(framedown, width=50)
    #     mdw_entry.pack()
    #     mdw_button = tk.Button(
    #         framedown, text="Browse", command=lambda: browse_file(mdw_entry)).pack(pady=20)

    #     grdw_label = tk.Label(
    #         framedown, text=".grd file (wave):", font='Times').pack(pady=5)
    #     grdw_entry = tk.Entry(framedown, width=50)
    #     grdw_entry.pack()
    #     grdw_button = tk.Button(
    #         framedown, text="Browse", command=lambda: browse_file(grdw_entry)).pack(pady=20)

    #     bndw_label = tk.Label(
    #         framedown, text=".bnd file (wave)\n boundary naming format\n Name1 :", font='Times').pack(pady=5)
    #     bndw_entry = tk.Entry(framedown, width=50)
    #     bndw_entry.pack()
    #     bndw_button = tk.Button(
    #         framedown, text="Browse", command=lambda: browse_file(bndw_entry)).pack(pady=20)

    #     ncw_label = tk.Label(
    #         framedown, text=".nc wave file (EasyGSH):", font='Times').pack(pady=5)
    #     ncw_entry = tk.Entry(framedown, width=50)
    #     ncw_entry.pack()
    #     ncw_button = tk.Button(
    #         framedown, text="Browse", command=lambda: browse_file(ncw_entry)).pack(pady=20)

    #     stepw_label = tk.Label(
    #         frameright, text="Time step for wave extraction:", font='Times').pack(pady=5)
    #     step_types_w = ["20 mins", "40 mins", "60 mins", "80 mins", "120 mins"]
    #     selected_step_w = tk.IntVar()

    #     for idx_w, step_type_w in enumerate(step_types_w, start=1):
    #         radio_button_w = tk.Radiobutton(
    #             frameright, variable=selected_step_w, text=step_type_w, value=idx_w)
    #         radio_button_w.pack(anchor='c')

    #     # Start date
    #     s_date_label = tk.Label(
    #         frameright, text="\nStart date of simulation\n (YYYY-MM-DD hh:mm:ss):", font='Times').pack(pady=5)
    #     s_date_entry = tk.Entry(frameright, width=50)
    #     s_date_entry.pack(pady=5)

    #     # End date
    #     e_date_label = tk.Label(
    #         frameright, text="\nEnd date of simulation\n (YYYY-MM-DD hh:mm:ss):", font='Times').pack(pady=5)
    #     e_date_entry = tk.Entry(frameright, width=50)
    #     e_date_entry.pack(pady=5)

    #     def check_submit():
    #         t = time.time()
    #         if not mdw_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .mdw file")
    #         elif not grdw_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .grd (wave) file.")
    #         elif not bndw_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .bnd (wave) file.")
    #         elif not ncw_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please browse for the .nc wave file.")
    #         elif not s_date_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please write the start date for simulation.")
    #         elif not e_date_entry.get():
    #             messagebox.showwarning(
    #                 "Warning", "Please write the end date for simulation.")
    #         else:
    #             # Input files from browsed
    #             step_w = selected_step_w.get()

    #             def step_number(step):
    #                 if step == 1:
    #                     step_c = 20
    #                 elif step == 2:
    #                     step_c = 40
    #                 elif step == 3:
    #                     step_c = 60
    #                 elif step == 4:
    #                     step_c = 80
    #                 else:
    #                     step_c = 120
    #                 return step_c

    #             step_wave = step_number(step_w)

    #             def value_from_txt_file(file, string_name):
    #                 file1 = open(file, "r")
    #                 for line in file1:
    #                     if '=' in line:
    #                         if string_name in line:
    #                             val = line.split('=')
    #                             string_val = val[1].strip()
    #                             break
    #                 file1.close()
    #                 return string_val

    #             # Files
    #             mdw_file = mdw_entry.get()
    #             grid_wave_input = grdw_entry.get()
    #             bnd_wave_input = bndw_entry.get()
    #             nc_file_wave = ncw_entry.get()
    #             start_time = s_date_entry.get()
    #             end_time = e_date_entry.get()

    #             # Output file
    #             wave_name_with_dot = mdw_file.partition('.')
    #             wave_name_until_dot = wave_name_with_dot[0]
    #             bcw_file = '{}.bcw'.format(wave_name_until_dot)
    #             wave_path_out_file = '{}.csv'.format(wave_name_until_dot)

    #             # Create the csv file for wave boundaries
    #             bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
    #                 path_bnd=bnd_wave_input)
    #             coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
    #                 path_grd=grid_wave_input, request_list=bnd_wave_grd_indices_output)
    #             crs_type_wave = extract_from_d3d_files.extract_crs_type(path_grd=grid_wave_input)
    #             output_methods.write_bnd_coord_ascii(
    #                 bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
    #             print(f'Grid type for wave is {crs_type_wave}')
    #             print('.')
    #             print('.')
    #             # Create the bcw file
    #             boundaries_wave = wave_path_out_file
    #             bcw = bcw_generator.bcw_file_generator(
    #                 boundaries_wave=boundaries_wave, nc_file_wave=nc_file_wave, mdw_file=mdw_file,
    #                 start_time=start_time, end_time=end_time, step_wave=step_wave,
    #                 bcw_file_name=bcw_file, crs_type=crs_type_wave)

    #             print(
    #                 'The process of extracting wave boundary conditions has now completed.')
    #             elapsed = (time.time() - t) / 60
    #             print(f'Total time taken: {np.round(elapsed, 2)} mins')

    #             # Write the new mdw file
    #             mdw_writer.write_mdw_file(
    #                 mdw_file=mdw_file, boundaries_wave=boundaries_wave)
    #             print('New mdw file created')

    #     # Submit button for execution
    #     frame_submit = tk.Frame(main_frame, width=200,
    #                             height=20, borderwidth=1, relief='solid')
    #     frame_submit.pack(side='bottom', fill='both', padx=10, pady=10)

    #     submit_button = tk.Button(
    #         frame_submit, text="Extract wave timeseries (bcw)", command=check_submit, width=30)
    #     submit_button.pack(pady=10, anchor='s')
# %% BCW overlap file generator
    def process_bcw_overlap_files(self, main_frame):

        # Frame for left
        framedown = tk.Frame(main_frame, borderwidth=1, relief='solid')
        framedown.pack( fill='both', expand=True, padx=10, pady=10)
        
        # Left section for MDF, GRD, and BND files
        left_frame = tk.Frame(framedown, borderwidth=1, relief='solid', padx=10, pady=10)
        left_frame.pack(side='left', fill='both', expand=True)
    
        # Title for the left frame
        left_label = tk.Label(left_frame, text="Wave Module", font=("Times", 16))
        left_label.pack(side='top', padx=10, pady=10)
    
        # Function to browse for files
        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)
    
        # Input buttons for Wave Module
        mdw_label = tk.Label(left_frame , text="MDW file:", font='Times')
        mdw_label.pack(pady=5)
        mdw_entry = tk.Entry(left_frame , width=50)
        mdw_entry.pack()
        mdw_button = tk.Button(left_frame , text="Browse", command=lambda: browse_file(mdw_entry))
        mdw_button.pack(pady=10)
    
        grdw_label = tk.Label(left_frame , text=".grd file (wave):", font='Times')
        grdw_label.pack(pady=5)
        grdw_entry = tk.Entry(left_frame , width=50)
        grdw_entry.pack()
        grdw_button = tk.Button(left_frame , text="Browse", command=lambda: browse_file(grdw_entry))
        grdw_button.pack(pady=10)
    
        bndw_label = tk.Label(left_frame , text=".bnd file (wave)\n boundary naming format\n Name1 :", font='Times')
        bndw_label.pack(pady=5)
        bndw_entry = tk.Entry(left_frame , width=50)
        bndw_entry.pack()
        bndw_button = tk.Button(left_frame , text="Browse", command=lambda: browse_file(bndw_entry))
        bndw_button.pack(pady=10)
    
        # Frame for right Module
        frameright = tk.Frame(framedown, borderwidth=1, relief='solid')
        frameright.pack(side='right', fill='both', expand=True, padx=10, pady=10)
    
        # List to store selected .nc wave files
        nc_files_wave = []
    
        def browse_nc_wave_files():
            """Browse and select multiple .nc wave files."""
            file_paths = filedialog.askopenfilenames(filetypes=[("NetCDF Files", "*.nc")])
            if file_paths:
                nc_files_wave.extend(file_paths)
                update_nc_wave_files_display()
    
        def update_nc_wave_files_display():
            """Update the display of selected .nc wave files."""
            for widget in ncw_files_frame.winfo_children():
                widget.destroy()
    
            for idx, file in enumerate(nc_files_wave, start=1):
                label = tk.Label(ncw_files_frame, text=f"{idx}. {file}", font='Times', anchor='w', justify='left', wraplength=300)
                label.pack(fill='x', padx=5, pady=2)
    
        # Title for .nc wave files section
        ncw_title = tk.Label(left_frame, text="Select all wave .nc files", font=("Times", 14))
        ncw_title.pack(anchor='n', pady=10)
    
        # Frame for displaying selected .nc wave files
        ncw_files_frame = tk.Frame(left_frame, borderwidth=1, relief='sunken', width=300, height=200)
        ncw_files_frame.pack(fill='both', padx=5, pady=5, expand=True)
    
        # Button to browse .nc wave files
        browse_ncw_button = tk.Button(left_frame, text="Browse .nc Wave Files", command=browse_nc_wave_files)
        browse_ncw_button.pack(pady=10)
    
        # Time Step Options
        stepw_label = tk.Label(frameright, text="Time step for wave extraction:", font='Times')
        stepw_label.pack(pady=10, anchor='w')
    
        selected_step_w = tk.IntVar(value=1)
        step_types_w = ["20 mins", "40 mins", "60 mins", "80 mins", "120 mins"]
        for idx_w, step_type_w in enumerate(step_types_w, start=1):
            radio_button_w = tk.Radiobutton(frameright, variable=selected_step_w, text=step_type_w, value=idx_w)
            radio_button_w.pack(anchor='w')
    
        # Start date
        s_date_label = tk.Label(frameright, text="\nStart date of simulation\n (YYYY-mm-dd HH:MM:SS):", font='Times')
        s_date_label.pack(pady=5)
        s_date_entry = tk.Entry(frameright, width=50)
        s_date_entry.pack(pady=5)
    
        # End date
        e_date_label = tk.Label(frameright, text="\nEnd date of simulation\n (YYYY-mm-dd HH:MM:SS):", font='Times')
        e_date_label.pack(pady=5)
        e_date_entry = tk.Entry(frameright, width=50)
        e_date_entry.pack(pady=5)
    
        def check_submit():
            t = time.time()
            if not mdw_entry.get():
                messagebox.showwarning("Warning", "Please browse for the .mdw file")
            elif not grdw_entry.get():
                messagebox.showwarning("Warning", "Please browse for the .grd (wave) file.")
            elif not bndw_entry.get():
                messagebox.showwarning("Warning", "Please browse for the .bnd (wave) file.")
            elif not nc_files_wave:
                messagebox.showwarning("Warning", "Please browse for the .nc wave files.")
            elif not s_date_entry.get():
                messagebox.showwarning("Warning", "Please write the start date for simulation.")
            elif not e_date_entry.get():
                messagebox.showwarning("Warning", "Please write the end date for simulation.")
            else:
                # Input files from browsed
                step_wave = [20, 40, 60, 80, 120][selected_step_w.get() - 1]
                mdw_file = mdw_entry.get()
                grid_wave_input = grdw_entry.get()
                bnd_wave_input = bndw_entry.get()
                start_time = s_date_entry.get()
                end_time = e_date_entry.get()
    
                # Output file
                wave_name_with_dot = mdw_file.partition('.')
                wave_name_until_dot = wave_name_with_dot[0]
                bcw_file = f"{wave_name_until_dot}.bcw"
                wave_path_out_file = f"{wave_name_until_dot}.csv"

                # Create the csv file for wave boundaries
                bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
                    path_bnd=bnd_wave_input)
                coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
                    path_grd=grid_wave_input, request_list=bnd_wave_grd_indices_output)
                crs_type_wave = extract_from_d3d_files.extract_crs_type(path_grd=grid_wave_input)
                output_methods.write_bnd_coord_ascii(
                    bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
                print(f'Grid type for wave is {crs_type_wave}')
                print('.')
                print('.')
                # Create the bcw file
                boundaries_wave = wave_path_out_file
                bcw = bcw_year_overlap_file_generator.bcw_year_overlap_file_generator(
                    boundaries_wave=boundaries_wave, nc_files_wave=nc_files_wave,
                    mdw_file=mdw_file, start_time=start_time,
                    end_time=end_time, step_wave=step_wave,
                    bcw_file_name=bcw_file, crs_type=crs_type_wave)
                
                print(
                    'The process of extracting wave boundary conditions has now completed.')
                elapsed = (time.time() - t) / 60
                print(f'Total time taken: {np.round(elapsed, 2)} mins')

                # Write the new mdw file
                mdw_writer.write_mdw_file(
                    mdw_file=mdw_file, boundaries_wave=boundaries_wave)
                print('New mdw file created')

        # Submit button for execution
        frame_submit = tk.Frame(frameright, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(side='bottom', fill='both', padx=10, pady=10)

        submit_button = tk.Button(
            frame_submit, text="Extract wave timeseries (bcw)", command=check_submit, width=30)
        submit_button.pack(pady=10, anchor='s')
        
# %% Create boundary files
    def generate_boundary_csv(self, main_frame):
        # Frame for Wave/Flow Module
        framedown = tk.Frame(main_frame, width=250,
                             borderwidth=1, relief='solid')
        framedown.pack(side='top', fill='both', expand=True, padx=10)
        right_label = tk.Label(
            framedown, text="Wave/Flow Module", font=("Times", 16))
        right_label.pack(side='top', padx=10, pady=10)
        t = time.time()
        # Function to browse for files

        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)

        # Input buttons for .grd and .bnd files
        grdw_label = tk.Label(framedown, text=".grd file:",
                              font='Times').pack(pady=5)
        grdw_entry = tk.Entry(framedown, width=50)
        grdw_entry.pack()
        grdw_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(grdw_entry)).pack(pady=20)

        bndw_label = tk.Label(framedown, text=".bnd file",
                              font='Times').pack(pady=5)
        bndw_entry = tk.Entry(framedown, width=50)
        bndw_entry.pack()
        bndw_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(bndw_entry)).pack(pady=20)

        def check_submit():
            t = time.time()
            if not grdw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .grd file.")
            elif not bndw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .bnd file.")
            else:
                # Input files
                grid_input = grdw_entry.get()
                bnd_input = bndw_entry.get()

                # Use grid file to extract bct file and output file
                wave_name_with_dot = grid_input.partition('.')
                wave_name_until_dot = wave_name_with_dot[0]
                wave_path_out_file = '{}.csv'.format(wave_name_until_dot)

                # Create the csv file for wave boundaries
                bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
                    path_bnd=bnd_input)
                coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
                    path_grd=grid_input, request_list=bnd_wave_grd_indices_output)
                crs_type = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
                output_methods.write_bnd_coord_ascii(
                    bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
                print(f'Grid type is {crs_type}')
                print('.')
                print('.')
                print(
                    'The process of creating the boundary location csv file is completed.')

        # Submit button for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(fill='both', padx=10, pady=10)

        submit_button = tk.Button(
            frame_submit, text="Create boundary locations CSV file", command=check_submit, width=30)
        submit_button.pack(pady=10, anchor='s')
# %% Generate boundary files and mdw files
    def generate_boundary_mdw(self, main_frame):
        # Frame for Wave Module
        framedown = tk.Frame(main_frame, width=250,
                             borderwidth=1, relief='solid')
        framedown.pack(side='top', fill='both', expand=True, padx=10)
        right_label = tk.Label(
            framedown, text="Wave Module", font=("Times", 16))
        right_label.pack(side='top', padx=10, pady=10)
        t = time.time()
        # Function to browse for files

        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)

        # Input buttons for .grd, .bnd, and .mdw files
        grdw_label = tk.Label(framedown, text=".grd file:",
                              font='Times').pack(pady=5)
        grdw_entry = tk.Entry(framedown, width=50)
        grdw_entry.pack()
        grdw_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(grdw_entry)).pack(pady=20)

        bndw_label = tk.Label(framedown, text=".bnd file:",
                              font='Times').pack(pady=5)
        bndw_entry = tk.Entry(framedown, width=50)
        bndw_entry.pack()
        bndw_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(bndw_entry)).pack(pady=20)

        mdw_label = tk.Label(framedown, text="MDW file:",
                             font='Times').pack(pady=5)
        mdw_entry = tk.Entry(framedown, width=50)
        mdw_entry.pack()
        mdw_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(mdw_entry)).pack(pady=20)

        def check_submit():
            t = time.time()
            if not grdw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .grd file.")
            elif not bndw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .bnd file.")
            elif not mdw_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .mdw file.")
            else:
                # Input files
                grid_input = grdw_entry.get()
                bnd_input = bndw_entry.get()
                mdw_file = mdw_entry.get()

                # Use grid file to extract bct file and output file
                wave_name_with_dot = grid_input.partition('.')
                wave_name_until_dot = wave_name_with_dot[0]
                wave_path_out_file = '{}.csv'.format(wave_name_until_dot)

                # Create the csv file for wave boundaries
                bnd_wave_grd_indices_output = extract_from_d3d_files.extract_bnd_grd_indices(
                    path_bnd=bnd_input)
                coord_from_d3d_wave_grd_output = extract_from_d3d_files.extract_coord_from_d3d_grd(
                    path_grd=grid_input, request_list=bnd_wave_grd_indices_output)
                crs_type = extract_from_d3d_files.extract_crs_type(path_grd=grid_input)
                output_methods.write_bnd_coord_ascii(
                    bnd_data_list=coord_from_d3d_wave_grd_output, out_path=wave_path_out_file)
                print(f'Grid type is {crs_type}')
                print('.')
                print('.')
                print('Boundary CSV file generated')
                # Create the new mdw file
                mdw_writer.write_mdw_file(
                    mdw_file=mdw_file, boundaries_wave=wave_path_out_file)
                print('New mdw file created.')
                elapsed_final = (time.time() - t) / 60
                print(
                    f'\nTotal time taken {np.round(elapsed_final, 2)} mins')

        # Submit button for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(fill='both', padx=10, pady=10)

        submit_button = tk.Button(
            frame_submit, text="Generate boundary csv and mdw files", command=check_submit, width=30)
        submit_button.pack(pady=10, anchor='s')
# %% Add Sea level
    def add_sea_level(self, main_frame):
        # Frame for Flow Module
        framedown = tk.Frame(main_frame, width=250,
                             borderwidth=1, relief='solid')
        framedown.pack(side='top', fill='both', expand=True, padx=10)
        right_label = tk.Label(
            framedown, text="Flow Module", font=("Times", 16))
        right_label.pack(side='top', padx=10, pady=10)
        t = time.time()
        # Function to browse for files

        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)

        # Input buttons for .bct file
        bct_label = tk.Label(framedown, text=".bct file:",
                             font='Times').pack(pady=5)
        bct_entry = tk.Entry(framedown, width=50)
        bct_entry.pack()
        bct_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(bct_entry)).pack(pady=20)

        # Radio buttons for type of sea level change
        type_choice = ["Gradual", "Constant"]
        selected_step_c = tk.IntVar()

        for idx_c, step_type_c in enumerate(type_choice, start=1):
            radio_button_c = tk.Radiobutton(
                framedown, variable=selected_step_c, text=step_type_c, value=idx_c)
            radio_button_c.pack(anchor='c')

        # Entry for amount of change
        change_label = tk.Label(
            framedown, text="Amount of sea level rise (m):", font='Times').pack(pady=5)
        change_entry = tk.Entry(framedown, width=50)
        change_entry.pack(pady=5)

        def check_submit_bct():
            if not bct_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the .bct file.")
            elif selected_step_c.get() == 0:
                messagebox.showwarning(
                    "Warning", "Please select type of sea level change.")
            else:
                bct_file_name = bct_entry.get()
                change_amount = float(change_entry.get())
                type_inc = selected_step_c.get() == 2  # True if constant, False if gradual
                sea_level_change.add_wl(
                    bct_file_name, sea_level_change=change_amount, constant=type_inc)
                type_add = "constant" if type_inc else "gradual"
                print(
                    f'Sea level change added to water levels, type {type_add}')

        # Submit button for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(fill='both', padx=10, pady=10)

        submit_button = tk.Button(
            frame_submit, text="Add sea level", command=check_submit_bct, width=20)
        submit_button.pack(pady=10, anchor='s')
#%% RPI module
    def identify_representative_period(self, main_frame):
        # Frame for Representative Period
        framedown = tk.Frame(main_frame, width=250,
                             borderwidth=1, relief='solid')
        framedown.pack(side='left', fill='both', expand=True, padx=10)
        left_label = tk.Label(
            framedown, text="Representative period", font=("Times", 16))
        left_label.pack(side='top', padx=10, pady=10)

        frameup = tk.Frame(main_frame, width=250, borderwidth=1, relief='solid')
        frameup.pack(side='right', fill='both', expand=True, padx=10)
        # Function to browse for files

        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)
        
        text_file = """
        WIND file: \n\nThe file should have only 3 columns, in the following order
        \ndate, speed, dir \nThe file should be comma separated (,)\n Please remove nans and inconsistencies before hand
        """
        
        # Input button for WIND file
        wind_label = tk.Label(framedown, text=text_file,
                              font=('Times',12)).pack(pady=5)
        wind_entry = tk.Entry(framedown, width=50)
        wind_entry.pack()
        wind_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(wind_entry)).pack(pady=20)

        # Output name
        output_label = tk.Label(
            framedown, text="Output filename:", font=('Times', 12)).pack(pady=5)
        output_entry = tk.Entry(framedown, width=50)
        output_entry.pack(pady=5)

        text = """
        \n            Format to insert parameters:
        \n\nDirectional orientations:\n['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE','SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW']\n
        \nSpeed classes:\n[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]\n
        \nMoving window duration:\n['2MS', '3MS', '4MS', '5MS', '6MS']\n
        """

        permanent_text_label = tk.Label(
            framedown, text=text, justify=tk.LEFT, wraplength=400, font=('Times', 12))
        permanent_text_label.pack()

        # Quad selection
        quad_label = tk.Label(frameup, text="\nDirectional orientations (minimum 2 or all):", font=(
            'Times', 12)).pack(pady=5)

        # Adding ALL toggle for quad_entry
        def toggle_quad():
            if all_quad_var.get():
                quad_entry.config(state='disabled')
            else:
                quad_entry.config(state='normal')

        all_quad_var = tk.BooleanVar()
        all_quad_check = tk.Checkbutton(
            frameup, text="ALL", variable=all_quad_var, command=toggle_quad)
        all_quad_check.pack()

        quad_entry = tk.Entry(frameup, width=50)
        quad_entry.pack(pady=5)

        # Speed selection
        spd_label = tk.Label(
            frameup, text="Speed classes (minimum 2 or all):", font=('Times', 12))
        spd_label.pack(pady=5)

        # Adding ALL toggle for spd_entry
        def toggle_spd():
            if all_spd_var.get():
                spd_entry.config(state='disabled')
            else:
                spd_entry.config(state='normal')

        all_spd_var = tk.BooleanVar()
        all_spd_check = tk.Checkbutton(
            frameup, text="ALL", variable=all_spd_var, command=toggle_spd)
        all_spd_check.pack()

        spd_entry = tk.Entry(frameup, width=50)
        spd_entry.pack(pady=5)

        # Frequency selection
        frq_label = tk.Label(frameup, text="Moving window duration (recommended upto '24MS'):",
                             font=('Times', 12)).pack(pady=5)
        frq_entry = tk.Entry(frameup, width=50)
        frq_entry.pack(pady=5)

        # Total reference period
        t_start_label = tk.Label(
            frameup, text="\n\nStart time of reference period:\n eg  : YYYY-MM-DD hh:mm:ss ", font=('Times', 12)).pack()
        t_start_entry = tk.Entry(frameup, width=50)
        t_start_entry.pack(pady=5)

        t_end_label = tk.Label(
            frameup, text="End time of reference period:", font=('Times', 12)).pack()
        t_end_entry = tk.Entry(frameup, width=50)
        t_end_entry.pack(pady=5)

        # Compare periods
        start_label = tk.Label(
            frameup, text="\n\nStart time of scanning window period", font=('Times', 12)).pack()
        start_entry = tk.Entry(frameup, width=50)
        start_entry.pack(pady=5)

        end_label = tk.Label(
            frameup, text=" End time of scanning window period\nmust be at most\nReference period end time - largest moving window duration ", font=('Times', 12)).pack()
        end_entry = tk.Entry(frameup, width=50)
        end_entry.pack(pady=5)

        def browse_path(entry_widget):
            file_path = filedialog.askdirectory()
            entry_widget.config(state='normal')
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)
            entry_widget.config(state='disabled')

        # Out path
        output_wind_label = tk.Label(frameup, text="\nSelect output path",
                                     font=('Times', 12)).pack(pady=5)
        output_wind_entry = tk.Entry(frameup, width=50)
        output_wind_entry.pack()
        output_wind_button = tk.Button(
            frameup, text="Browse", command=lambda: browse_path(output_wind_entry)).pack(pady=20)


        def parse_list_input(input_str):
            # Parse the input string as a Python expression (list)
            converted_list = ast.literal_eval(input_str)
            if isinstance(converted_list, list):
                return converted_list

        def check_submit_rr():
            t = time.time()
            if not wind_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the wind file.")
            else:
                file_input = wind_entry.get()
                output_name = output_entry.get()
                # Use preset values if "ALL" is checked
                quad = (['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW']
                        if all_quad_var.get() else parse_list_input(quad_entry.get()))
                spd = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
                       if all_spd_var.get() else parse_list_input(spd_entry.get()))
                frequency = parse_list_input(frq_entry.get())
                start_time_total = t_start_entry.get()
                end_time_total = t_end_entry.get()
                start_time = start_entry.get()
                end_time = end_entry.get()
                out_path = output_wind_entry.get()

                output_file = rep_period.identify_rep_period(
                    file_input, output_name, quad, spd, start_time_total,
                    end_time_total, frequency, start_time, end_time, out_path)
                print('Representative period file has been generated')
                elapsed_final = (time.time() - t) / 60
                print(
                    f'\nTotal time taken {np.round(elapsed_final, 2)} mins')

                # Create a new window to display the plot
                plot_window = Toplevel(main_frame)
                plot_window.title("Windroses Plot")

                plot_windroses.plot_windroses(file_input, start_time_total, end_time_total,
                                              output_file, parent=plot_window)

        # Submit button for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(side='bottom', fill='both', padx=10, pady=10)
        
        # Load and display the image
        image_path = self.resource_path('classes.jpg')
        img = Image.open(image_path)
        img = img.resize((425, 500), Image.LANCZOS)  # Resize the image if necessary
        img = ImageTk.PhotoImage(img)
        image_label = tk.Label(frame_submit, image=img)
        image_label.image = img  # Keep a reference to avoid garbage collection
        image_label.pack(padx=10, pady=10)

        # Citation text 
        text_2 = "\nCitation for use of Rep period algorithm:\nSoares, C.C., Galiforni-Silva, F., Winter, C., 2024. Representative residual transport pathways in a mixed-energy open tidal system. Journal of Sea Research 201, 102530. https://doi.org/10.1016/j.seares.2024.102530\n"

        permanent_text_label_2 = tk.Label(
            frame_submit, text=text_2, justify=tk.LEFT, wraplength=400, font=('Times', 11))
        permanent_text_label_2.pack()
                
        submit_button = tk.Button(
            frame_submit, text="Identify Rep period", command=check_submit_rr, width=20)
        submit_button.pack(pady=15, anchor='s')
#%% COSMO Wind generator
    def generate_wind_files(self, main_frame):
        # Frame for Wave Module
        framedown = tk.Frame(main_frame, width=250,
                             borderwidth=1, relief='solid')
        framedown.pack(side='left', fill='both', expand=True, padx=10)
        right_label = tk.Label(
            framedown, text="COSMO- Wind files", font=("Times", 16))
        right_label.pack(side='top', padx=10, pady=10)

        frameright = tk.Frame(main_frame, width=250,
                              borderwidth=1, relief='solid')
        frameright.pack(side='right', fill='both', expand=True, padx=10)
        label = tk.Label(
            frameright, text="Important notes", font=("Times", 16))
        label.pack(side='top', padx=10, pady=10)

        t = time.time()

        # Function to browse for files

        def browse_path(entry_widget):
            file_path = filedialog.askdirectory()
            entry_widget.config(state='normal')
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)
            entry_widget.config(state='disabled')

        # Input buttons for .grd, .bnd, and .mdw files
        cosmo_file_path_label = tk.Label(framedown, text="Main path with COSMO (U,V,PS) data folders:",
                                         font='Times').pack(pady=5)
        
        cosmo_file_path_entry = tk.Entry(framedown, width=50, state='readonly')
        cosmo_file_path_entry.pack()
        cosmo_file_path_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_path(cosmo_file_path_entry))
        cosmo_file_path_button.pack(pady=20)

        ref_time_label = tk.Label(framedown, text="Reference time: \nRead notes\n eg YYYY-MM-DD hh:mm:ss",
                                  font='Times').pack(pady=5)
        ref_time_entry = tk.Entry(framedown, width=50)
        ref_time_entry.pack()

        output_filename_label = tk.Label(framedown, text="Output file name:\neg: cosmo_2011",
                                         font='Times').pack(pady=5)
        output_filename_entry = tk.Entry(framedown, width=50)
        output_filename_entry.pack()
        
        def browse_file(entry_widget):
            file_path = filedialog.askopenfilename()
            entry_widget.delete(0, "end")
            entry_widget.insert(0, file_path)
            
        grid_ed_path = tk.Label(framedown, text="\nDelft3d Equidistant grid (6km resolution) as .mat:\neg: GNS_6km.mat",
                                         font='Times').pack(pady=5)
        grid_ed_path_entry = tk.Entry(framedown, width=50)
        grid_ed_path_entry.pack()
        grid_ed_button = tk.Button(
            framedown, text="Browse", command=lambda: browse_file(grid_ed_path_entry)).pack(pady=20)
        
        # Adding radio buttons for coordinate system selection
        coord_system_label = tk.Label(framedown, text="Select Coordinate System of Delft3D Grid:",
                                       font='Times').pack(pady=10)
        coord_system = tk.StringVar(value="Cartesian")  # Default selection
        
        cartesian_button = tk.Radiobutton(framedown, text="Cartesian", variable=coord_system, value="Cartesian",
                                          font='Times')
        cartesian_button.pack(anchor='w', padx=20)
        
        spherical_button = tk.Radiobutton(framedown, text="Spherical", variable=coord_system, value="Spherical",
                                          font='Times')
        spherical_button.pack(anchor='w', padx=20)

        text = """
        \nThe main path containing the COSMO files should be structured as follows:\n\nFolder 1 should be named UV and should have all the U and V monthly cosmo grib files you wish to extract from.
        \nFolder 2 should be named PS and should have all the PS monthly grib files.
        \nThe reference time should be the same as your model reference date. However it is better to have the real COSMO data start a little after the reference date, preferable download data from a month later so if your reference date is 2001-01-01 00:00:00 then cosmo data should start from 2001-02-01 00:00:00
        \nProvide an equidistant 6X6 km grid, extracted from Delft3D in .mat format
        \nIn the current RepD3D version the components are rotated to true North
        \n Output files are .amu, .amv, .amp, .xwind and .ywind
        """

        permanent_text_label = tk.Label(
            frameright, text=text, justify=tk.LEFT, wraplength=400, font=('Times', 12))
        permanent_text_label.pack()

        def check_submit():
            t = time.time()
            if not cosmo_file_path_entry.get():
                messagebox.showwarning(
                    "Warning", "Please browse for the path with the U,V,PS files")
            elif not output_filename_entry.get():
                messagebox.showwarning(
                    "Warning", "Please provide the output file name.")
            elif not ref_time_entry.get():
                messagebox.showwarning(
                    "Warning", "Please provide the reference time.")
            else:
                # Input files
                script_dir = os.path.dirname(os.path.abspath(__file__))
                const_path = os.path.join(script_dir, 'cosmo_latlon_rot.grb')
                grid_ed_path = grid_ed_path_entry.get()
                cosmo_files_path = cosmo_file_path_entry.get()
                file_name = output_filename_entry.get()
                ref_time = ref_time_entry.get()
                crs_type = coord_system.get()

                # Delete .idx files before starting the process
                for root, dirs, files in os.walk(cosmo_files_path):
                    for file in files:
                        if file.endswith('.idx'):
                            os.remove(os.path.join(root, file))

                cosmo_wind_file_generator.create_wind_fields_cosmo(grid_ed_path,
                                                                   const_path,
                                                                   cosmo_files_path,
                                                                   file_name,
                                                                   ref_time, crs_type=crs_type)

                print('S-T varying wind field file have been generated')
                elapsed_final = (time.time() - t) / 60
                print(
                    f'\nTotal time taken {np.round(elapsed_final, 2)} mins')

        # Submit button for execution
        frame_submit = tk.Frame(main_frame, width=200,
                                height=20, borderwidth=1, relief='solid')
        frame_submit.pack(side='bottom', fill='both', padx=10, pady=10)
        
        # Citation text 
        text_2 = """
        \nThe COSMO files can be found at:
        \nhttps://opendata.dwd.de/climate_environment/REA/COSMO_REA6/hourly/2D/ 
        \nOn the webpage look for PS, U_10M and V_10M and download all monthly files necesssary and unzip them - use 7-Zip. Delete the zip files before generating the wind field files.
        \nCitation for using COSMO data: \nBollmeyer, C., Keller, J.D., Ohlwein, C., Wahl, S., Crewell, S., Friederichs, P., Hense, A., Keune, J., Kneifel, S., Pscheidt, I., Redl, S., Steinke, S., 2015. Towards a high‐resolution regional reanalysis for the European CORDEX domain. Q.J.R. Meteorol. Soc. 141, 1–15. https://doi.org/10.1002/qj.2486\n
        """

        permanent_text_label_2 = tk.Label(
            frame_submit, text=text_2, justify=tk.LEFT, wraplength=400, font=('Times', 11))
        permanent_text_label_2.pack()

        # Load and display the image
        image_path = self.resource_path('Cosmo_extent.jpg')
        img = Image.open(image_path)
        img = img.resize((400, 200), Image.LANCZOS)  # Resize the image if necessary
        img = ImageTk.PhotoImage(img)
        image_label = tk.Label(frame_submit, image=img)
        image_label.image = img  # Keep a reference to avoid garbage collection
        image_label.pack(padx=10, pady=10)
                
        submit_button = tk.Button(
            frame_submit, text="Generate wind files", command=check_submit, width=30)
        submit_button.pack(pady=10, anchor='s')


if __name__ == "__main__":

    app = Application()
    app.mainloop()
