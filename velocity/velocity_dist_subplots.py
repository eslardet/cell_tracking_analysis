import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)
from analysis_functions_tracks import *
from well_plate_dictionary import *
import numba
import seaborn as sns
import matplotlib.lines as mlines
from scipy import stats



plate = 1737

stim = "With stimulation"
t_cell = "Specific CD4"

all_wells = []
for group in all_groups:
    all_wells += get_wells(group, stim, t_cell)

plot_wells = all_wells

# plot_wells = []
# plot_groups = ["Uninfected healthy control", "HAM"]
# for group in plot_groups:
#     plot_wells += get_wells(group, stim, t_cell)

## Subplots over time
time_ints = [[i,i+12] for i in range(0,72,12)]

fig, axes = plt.subplots(2, 3, figsize=(15,10), sharey=True, sharex=True)

for i, ax in enumerate(axes.flat):
    for well in all_wells:
        if well not in plot_wells:
            next(ax._get_lines.prop_cycler)
        else:
            ax = plot_velocity_distribution(plate, well, time_ints[i], kde=False, bins=20, max_hist=60, ax=ax)
    ax.set_title(str(time_ints[i][0]) + "-" + str(time_ints[i][1]) + "hrs")
    ax.set_xlim(0,60)
    ax.set_ylim(0,0.1)

handles = []
# colours = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]
# colours = ["tab:red", "tab:purple", "tab:brown"]
cmap = plt.get_cmap("tab10")

for j in range(len(plot_wells)):
    well = plot_wells[j]
    group, t_cell, stim = get_well_info(well)
    handles.append(mlines.Line2D([], [], color=cmap[j], label=group + ", " + str(plate) + ", " + well))
fig.legend(handles=handles, loc='upper center')
folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/velocity_distribution/subplots/"
# filename = str(plate) + "_" + well + "_hist_thresh4.png"
# filename = str(plate) + "_" + str(time_int[0]) + "-" + str(time_int[1]) + "hrs_ham_vs_control.png"
filename = str(plate) + "_ham_vs_control_int12_bin20.png"
if not os.path.exists(folder):
    os.makedirs(folder)
# plt.savefig(folder + filename, bbox_inches='tight')
# plt.close()
plt.show()




