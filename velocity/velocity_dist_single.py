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

well = "G2"

stim = "With stimulation"
t_cell = "Specific CD4"
all_wells = []
for group in all_groups:
    all_wells += get_wells(group, stim, t_cell)
time_interval = [24,30]

fig, ax = plt.subplots()


## Compare time intervals
# for time_int in [[i,i+12] for i in range(0,60,12)]:
#     ax = plot_velocity_distribution(plate, well, time_int, kde=False, bins=20, max_hist=60, ax=ax)
# filename = str(plate) + "_ham_vs_control_int12_bin20.png"

## Compare groups
for well in all_wells:
    ax = plot_velocity_distribution(plate, well, time_interval, kde=False, bins=20, ax=ax)
filename = str(plate) + "_" + well + "_hist_thresh4.png"

ax.set_xlim(0,60)
ax.set_ylim(0,0.06)

folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/velocity_distribution/subplots/"
if not os.path.exists(folder):
    os.makedirs(folder)
# plt.savefig(folder + filename, bbox_inches='tight')
# plt.close()
plt.show()




