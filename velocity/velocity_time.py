import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)
from analysis_functions_tracks import *
from well_plate_dictionary import *
import numba
import seaborn as sns
import matplotlib.lines as mlines
from scipy import stats



## Velocity vs time

plate = 1738
av_type = 'median'
# all_wells = ["C2", "C4", "C6", "C8", "C10", "G2"]
# all_wells = ['C3', 'C5', 'C7', 'C9', 'C11', 'G3']
# all_wells = ['E2', 'E4', 'E6', 'E8', 'E10', 'G4']

stim = "No stimulation"
t_cell = "Specific CD4"

for stim in all_stim:
    for t_cell in all_t_cell:
        
        all_wells = []

        for group in all_groups:
            all_wells += get_wells(group, stim, t_cell)
        print(all_wells)

        fig, ax = plt.subplots()

        # next(ax._get_lines.prop_cycler)
        for well in all_wells:
            try:
                filepath = get_xml_file_vel(plate, well)
                tracks = read_xml(filepath)
                # tracks = np.array(read_xml(filepath), dtype=object)
                print(well, len(tracks))
                time_plot, vel_plot = get_velocities(tracks)
                time_bins, vel_bins_mean, vel_bins_median, vel_bins_mode, vel_sd = bin_velocities(time_plot, vel_plot, max_frames=216, n_bins=12)
                # ax.scatter(np.array(time_plot), vel_plot, marker='o', color="tab:blue", s=10, alpha=0.5)
                group, t_cell, stim = get_well_info(well)
                if av_type == 'mean':
                    ax.plot(time_bins, vel_bins_mean, "-o", label=group + ", " + well)
                elif av_type == 'median':
                    ax.plot(time_bins, vel_bins_median, "-o", label=group + ", " + well)
                # ax.plot(np.array(time_bins), vel_bins_mode, "-o", label=group + ", " + well)
                # ax.errorbar(np.array(time_bins), vel_bins, yerr=vel_sd, fmt='none', capsize=2, ecolor="tab:blue", alpha=0.5)
            except:
                print("No tracks for " + well)
                next(ax._get_lines.prop_cycler)
        ax.set_xlim(0,72)
        # ax.set_ylim(8,32)
        ax.set_xlabel("Time (hrs)")
        if av_type == 'mean':
            ax.set_ylabel("Mean Instantaneous speed (pixels/hr)")
        elif av_type == 'median':
            ax.set_ylabel("Median Instantaneous speed (pixels/hr)")
        ax.legend()
        ax.set_title("Plate " + str(plate))

        folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/velocity_time/" + av_type + "/"
        # file_name = str(plate) + "_CD4_stim_thresh5_b12.png"
        file_name = str(plate) + "_" + stim + "_" + t_cell + "_b12.png"
        # file_name = str(plate) + "_all_thresh5_b12.png"
        # file_name = str(plate) + "_" + well + "_errorbars.png"
        if not os.path.exists(folder):
            os.makedirs(folder)
        plt.savefig(folder + file_name)
        plt.close()
        # plt.show()


