import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

from analysis_functions_image import *
from analysis_functions_tracks import *
import numpy as np
import csv
from scipy.stats import wilcoxon, ranksums, ttest_ind
import seaborn as sns
import time

def read_coloc_increase(cell_type, stim, t_cell, slice_compare):
    data_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plot_data/colocalization_increase/"
    file_name = cell_type + "_" + stim + "_" + t_cell + "_" + str(slice_compare[1]) + ".txt"
    with open(data_folder + file_name, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        data = np.array(list(reader)[4:])

    group = np.array(data[:, 0], dtype=int)
    manders_increase = np.array(data[:, 3], dtype=float)

    return group, manders_increase

def plot_coloc(cell_type, stim, t_cell, slice_compare, ax=None):

    group, manders_increase = read_coloc_increase(cell_type, stim, t_cell, slice_compare)

    control_group = manders_increase[np.where(group == 1)[0]]
    AC_lPVL_group = manders_increase[np.where(group == 2)[0]]
    AC_hPVL_group = manders_increase[np.where(group == 3)[0]]
    HAM_group = manders_increase[np.where(group == 4)[0]]

    if ax == None:
        fig, ax = plt.subplots()

    # ax.boxplot([control_group, AC_lPVL_group, AC_hPVL_group, HAM_group], patch_artist=True, boxprops=dict(facecolor='tab:blue', color='black', alpha=0.5), medianprops=dict(color='black'))
    # ax.scatter(group, manders_increase, marker='^', s=12, color='k')

    sns.boxplot(x=group, y=manders_increase, ax=ax)
    sns.swarmplot(x=group, y=manders_increase, color=".25", ax=ax)

    # ax.set_xticks([0,1,2,3], ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"], fontsize=8)
    # ax.set_ylabel("Manders colocalization coefficient increase (%)")
    # ax.set_title(cell_type + ", " + stim + ", " + t_cell + ", t=" + str(slice_compare[1]//3) + "hrs")
    # ax.set_ylim([-10,300])

    return ax

    # folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization_increase/Microglia_CD4_stim/box/"
    # file_name = cell_type + "_" + stim + "_" + t_cell + "_" + str(slice_compare[1]) + ".png"
    # plt.savefig(folder + file_name)
    # plt.close()
    # plt.show()



cell_type = 'Microglia & T-cells'
stim = "With stimulation"
t_cell = "Specific CD4"
slice_compare = [3, 108]

t0 = time.time()
# plot_coloc(cell_type, stim, t_cell, slice_compare)

fig, axes = plt.subplots(2, 3, figsize=(15,10), sharey=True, sharex=True)
slice_compare_list = [36, 72, 108, 144, 180, 216]

for i, ax in enumerate(axes.flat):
    ax = plot_coloc(cell_type, stim, t_cell, [3, slice_compare_list[i]], ax=ax)
    ax.set_xticks([0,1,2,3], ["Uninfected\n healthy control", "AC lPVL", "AC hPVL", "HAM"], fontsize=8)
    # ax.get_xaxis().set_visible(False)
    ax.set_ylabel("Manders colocalization coefficient increase (%)")
    ax.set_title("t=" + str(slice_compare_list[i]//3) + "hrs")
    ax.set_ylim([-10,300])

fig.legend(["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"], loc='upper center', ncol=4)

folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization_increase/Microglia_CD4_stim/box/"
if not os.path.exists(folder):
    os.makedirs(folder)
file_name = "all.png"
# plt.savefig(folder + file_name)
# plt.close()

print(time.time()-t0)

plt.show()


## Stats ##
# group, manders_increase = read_coloc_increase(cell_type, stim, t_cell, slice_compare)

# control_group = manders_increase[np.where(group == 1)[0]]
# AC_lPVL_group = manders_increase[np.where(group == 2)[0]]
# AC_hPVL_group = manders_increase[np.where(group == 3)[0]]
# HAM_group = manders_increase[np.where(group == 4)[0]]


# # control_group = np.concatenate((control_group, AC_lPVL_group, AC_hPVL_group))
# # print(control_group, HAM_group)

# # print(wilcoxon(control_group, HAM_group))

# print(ranksums(control_group, HAM_group))
# print(ranksums(control_group, HAM_group, alternative='less'))

# print(ttest_ind(control_group, HAM_group, equal_var=False))

# def t_test(x,y,alternative='both-sided'):
#     _, double_p = ttest_ind(x,y,equal_var = False)
#     if alternative == 'both-sided':
#         pval = double_p
#     elif alternative == 'greater':
#         if np.mean(x) > np.mean(y):
#             pval = double_p/2.
#         else:
#             pval = 1.0 - double_p/2.
#     elif alternative == 'less':
#         if np.mean(x) < np.mean(y):
#             pval = double_p/2.
#         else:
#             pval = 1.0 - double_p/2.
#     return pval

# print(t_test(control_group, HAM_group, alternative='less'))