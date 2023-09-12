import os
import csv
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/coloc/"

file_list = ["AutoCorrelation on VID1737_red_F4_1_t0", "AutoCorrelation on VID1737_red_F4_1_end"]
# file_name = "AutoCorrelation on VID1737_red_F4_1.tif kept stack full"

def func(r, alpha, coeff):
    return coeff*np.exp(-r/alpha)

def fit_curve(r_plot, corr):
    alpha, coeff = curve_fit(func, r_plot, corr)[0]
    return alpha, coeff

fig, ax = plt.subplots()
for file_name in file_list:
    file_path = folder + file_name
    with open(file_path) as f:
        reader = csv.reader(f, delimiter="\t")
        r = list(reader)[1:200]


    radius = []
    for item in r:
        radius.append(float(item[0]))

    corr = []
    for item in r:
        corr.append(np.abs(float(item[1])))

    ax.plot(radius, corr, label=file_name)
    alpha, coeff = fit_curve(radius, corr)
    print(alpha)
    r_plot = np.linspace(0, np.max(radius), 100)
    ax.plot(r_plot, func(r_plot, alpha, coeff), '--', label=r"$\xi=$" + str(round(alpha,2)))

ax.set_xlabel("r (pixels)")
ax.set_ylabel("C(r)")
# ax.set_ylim([0,1])
ax.legend()


plate = "1737"
well = "F4"
plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation/"
plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')

# plt.show()