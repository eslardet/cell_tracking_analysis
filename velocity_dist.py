from analysis_functions_tracks import *
from well_plate_dictionary import *
import numba
import seaborn as sns
import matplotlib.lines as mlines
from scipy import stats



# print(tracks.shape)
# print(len(tracks))
# print(tracks[0][0])
plate = 1738
well = "G2"
all_wells = ["C2", "C4", "C6", "C8", "C10", "G2"]
# plot_wells = ["C8", "C10", "G2"]
plot_wells = all_wells
# time_interval = [24,30]

## Velocity distributions

fig, axes = plt.subplots(2, 3, figsize=(15,10), sharey=True, sharex=True)
time_ints = [[i,i+12] for i in range(0,72,12)]


for i, ax in enumerate(axes.flat):
    for well in all_wells:
        if well in ["C2", "C4", "C6"]:
            next(ax._get_lines.prop_cycler)
        else:
            ax = plot_velocity_distribution(plate, well, time_ints[i], kde=False, bins=20, max_hist=60, ax=ax)
    ax.set_title(str(time_ints[i][0]) + "-" + str(time_ints[i][1]) + "hrs")
    ax.set_xlim(0,60)
    ax.set_ylim(0,0.1)

handles = []
colours = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]
# colours = ["tab:red", "tab:purple", "tab:brown"]

for j in range(len(plot_wells)):
    well = plot_wells[j]
    group, t_cell, stim = get_well_info(well)
    handles.append(mlines.Line2D([], [], color=colours[j], label=group + ", " + str(plate) + ", " + well))
fig.legend(handles=handles, loc='upper center')

    
# fig, ax = plt.subplots()

# # for time_int in [[i,i+12] for i in range(0,60,12)]:
# #     ax = plot_velocity_distribution(plate, well, time_int, kde=False, bins=20, ax=ax)

# for well in all_wells:
#     ax = plot_velocity_distribution(plate, well, time_int, kde=False, bins=40, max_hist=60, ax=ax)

# ax.set_xlim(0,60)
# ax.set_ylim(0,0.06)

folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/velocity_distribution/subplots/"
# filename = str(plate) + "_" + well + "_hist_thresh4.png"
# filename = str(plate) + "_" + str(time_int[0]) + "-" + str(time_int[1]) + "hrs_ham_vs_control.png"
filename = str(plate) + "_ham_vs_control_int12_bin20.png"
if not os.path.exists(folder):
    os.makedirs(folder)
plt.savefig(folder + filename, bbox_inches='tight')
plt.close()
# plt.show()



