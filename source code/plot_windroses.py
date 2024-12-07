def plot_windroses(overall_file, total_start, total_end, rep_file, parent=None, initial_index=0):
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from windrose import WindroseAxes
    import pandas as pd
    import seaborn as sns
    import numpy as np
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    import tkinter as tk
    import matplotlib.gridspec as gridspec
    import matplotlib.colorbar as colorbar

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 10

    spo_wind = pd.read_csv(overall_file, delimiter=',', skiprows=1, names=['date', 'speed', 'dir'])

    spo_wind_2 = spo_wind[(spo_wind['dir'] <= 360) & (spo_wind['dir'] >= 0) &
                          (spo_wind['speed'] >= 0) & (spo_wind['speed'] <= 35)]

    date = pd.to_datetime(spo_wind_2['date'].astype(str), format='%Y%m%d%H')

    data = [date, spo_wind_2['speed'], spo_wind_2['dir']]
    headers = ["datetime", "speed", "dir"]
    obs_data = pd.concat(data, axis=1, keys=headers)
    obs_data_all = obs_data.set_index(['datetime'])

    start_time = total_start
    end_time = total_end

    rep_file = pd.read_csv(rep_file)
    
    sns.set(style='white')
    font_prop_1 = FontProperties(family='Times New Roman', size=14)
    font_prop_2 = FontProperties(family='Times New Roman', size=12)

    fig_width_cm = 20
    fig_height_cm = 10  
    fig_width_in = fig_width_cm / 2.54
    fig_height_in = fig_height_cm / 2.54

    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    
    gs = gridspec.GridSpec(2, 3, height_ratios=[4, 0.2])

    ax1 = fig.add_subplot(gs[0, 0], projection='windrose')
    ax2 = fig.add_subplot(gs[0, 1], projection='windrose')
    ax3 = fig.add_subplot(gs[0, 2], projection='windrose')

    cbar_ax = fig.add_subplot(gs[1, :])
    
    canvas = None # Initialize canvas
    
    def plot_for_index(index):
        cp_start = rep_file['start_point'][index]
        cp_end = rep_file['end_point'][index]
        freq = rep_file['period_freq'][index]

        obs_data_97_16 = obs_data_all.loc[start_time:end_time]
        prospect_data = obs_data_all.loc[cp_start:cp_end]

        ax1.clear()
        ax2.clear()
        ax3.clear()
        cbar_ax.clear()

        def plot_windrose(ax, data_dir, data_speed, title):
            ax.contourf(data_dir, data_speed, bins=[
                0, 0.5, 1.5, 3.3, 5.5, 7.9, 10.7, 13.8, 17.1, 20.7, 24.4, 28.4], nsector=36, cmap=plt.cm.Reds, normed=True)
            ax.set_thetagrids(range(0, 360, 90), labels=[90, 0, 270, 180])
            ax.set_yticks(np.arange(0, 8, step=1))
            ax.set_yticklabels(np.arange(0, 8, step=1))
            ax.tick_params(axis='x', labelsize=12, pad=2)
            ax.tick_params(axis='y', labelsize=12, color='white')
            ax.grid(True, linestyle='--', alpha=0.8)
            ax.set_theta_zero_location('W', offset=-180)
            ax.set_title(title, fontsize=14, fontproperties=font_prop_2)

        plot_windrose(ax1, obs_data_97_16['dir'], obs_data_97_16['speed'], f'Total reference period\n{start_time} to {end_time}')
        plot_windrose(ax2, prospect_data['dir'], prospect_data['speed'], f'Top rep period\n{cp_start} to {cp_end}\n{freq} months')

        ax3.contour(obs_data_97_16['dir'], obs_data_97_16['speed'], bins=[
            0, 0.5, 1.5, 3.3, 5.5, 7.9, 10.7, 13.8, 17.1, 20.7, 24.4, 28.4], nsector=36, colors='k', normed=True)

        ax3.contourf(prospect_data['dir'], prospect_data['speed'], bins=[
            0, 0.5, 1.5, 3.3, 5.5, 7.9, 10.7, 13.8, 17.1, 20.7, 24.4, 28.4], nsector=36, cmap=plt.cm.Reds, normed=True, alpha=0.6)

        ax3.set_thetagrids(range(0, 360, 90), labels=[90, 0, 270, 180])
        ax3.set_yticks(np.arange(0, 8, step=1))
        ax3.set_yticklabels(np.arange(0, 8, step=1))
        ax3.tick_params(axis='x', labelsize=12, pad=2)
        ax3.tick_params(axis='y', labelsize=12, color='white')
        ax3.grid(True, linestyle='--', alpha=0.8)
        ax3.set_theta_zero_location('W', offset=-180)
        ax3.set_title('Overlap Heatmap', fontsize=14, fontproperties=font_prop_2)

        norm = plt.Normalize(vmin=0, vmax=28.4)
        cbar = colorbar.ColorbarBase(cbar_ax, cmap=plt.cm.Reds, norm=norm, orientation='horizontal')
        cbar.set_label('Wind speed (m/s)', fontproperties=font_prop_2, fontsize=12)

        if parent:
            canvas.draw()
        else:
            plt.draw()
           
    if parent:
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    else:
        plot_for_index(initial_index)
        plt.show()
    
    # Function to handle next button click
    def next_plot():
        nonlocal initial_index
        if initial_index < len(rep_file) - 1:
            initial_index += 1
        plot_for_index(initial_index)
    
    # Function to handle previous button click
    def prev_plot():
        nonlocal initial_index
        if initial_index > 0:
            initial_index -= 1
        plot_for_index(initial_index)
    
        # Function to handle first button click
    def first_plot():
        nonlocal initial_index
        initial_index = 0
        plot_for_index(initial_index)
    
    # Navigation buttons
    if parent:
        nav_frame = tk.Frame(parent)
        nav_frame.pack(side=tk.BOTTOM, fill=tk.X)
        prev_button = tk.Button(nav_frame, text="Previous", command=prev_plot)
        prev_button.pack(side=tk.LEFT)
        first_button = tk.Button(nav_frame, text="First", command=first_plot)
        first_button.pack(side =tk.LEFT, expand=True, anchor=tk.CENTER)
        next_button = tk.Button(nav_frame, text="Next", command=next_plot)
        next_button.pack(side=tk.RIGHT)
