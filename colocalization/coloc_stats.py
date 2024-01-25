import os, sys, inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0, parentdir)

from analysis_functions_image import *
from analysis_functions_tracks import *
import numpy as np
import csv
from scipy.stats import wilcoxon, ranksums, ttest_ind
import seaborn as sns

def read_coloc_increase(cell_type, stim, t_cell, slice_compare):
    data_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plot_data/colocalization_increase/"
    file_name = cell_type + "_" + stim + "_" + t_cell + "_" + str(slice_compare[1]) + ".txt"
    with open(data_folder + file_name, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        data = np.array(list(reader)[4:])

    group = np.array(data[:, 0], dtype=int)
    manders_increase = np.array(data[:, 3], dtype=float)

    return group, manders_increase

def t_test(x,y,alternative='both-sided'):
    _, double_p = ttest_ind(x,y,equal_var = False)
    if alternative == 'both-sided':
        pval = double_p
    elif alternative == 'greater':
        if np.mean(x) > np.mean(y):
            pval = double_p/2.
        else:
            pval = 1.0 - double_p/2.
    elif alternative == 'less':
        if np.mean(x) < np.mean(y):
            pval = double_p/2.
        else:
            pval = 1.0 - double_p/2.
    return pval


cell_type = 'Microglia & T-cells'
stim = "With stimulation"
t_cell = "Specific CD4"
slice_compare = [3, 108]


group, manders_increase = read_coloc_increase(cell_type, stim, t_cell, slice_compare)

control_group = manders_increase[np.where(group == 1)[0]]
AC_lPVL_group = manders_increase[np.where(group == 2)[0]]
AC_hPVL_group = manders_increase[np.where(group == 3)[0]]
HAM_group = manders_increase[np.where(group == 4)[0]]


# control_group = np.concatenate((control_group, AC_lPVL_group, AC_hPVL_group))
# print(control_group, HAM_group)

# print(wilcoxon(control_group, HAM_group))

print(ranksums(control_group, HAM_group))
print(ranksums(control_group, HAM_group, alternative='less'))

print(ttest_ind(control_group, HAM_group, equal_var=False))



print(t_test(control_group, HAM_group, alternative='less'))