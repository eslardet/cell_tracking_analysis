import os
import csv
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit




def func(r, xi, coeff):
    return np.exp(coeff)*np.exp(-r/xi)

def fit_curve(r_plot, corr):
    xi, coeff = curve_fit(func, r_plot, corr)[0]
    return xi, coeff

def fit_log_line(r_plot, corr, r_min, r_max):
    r_plot = np.array(r_plot)
    corr = np.array(corr)
    idx1 = np.where(r_plot<r_max)[0]
    idx2 = np.where(r_plot>r_min)[0]
    idx = list(set(idx1) & set(idx2))
    corr = corr[idx]
    r_plot = r_plot[idx]

    alpha, coeff = np.polyfit(r_plot, np.log(corr), deg=1)
    xi = -1/alpha
    return xi, coeff


def plot_fit(file_list, save_plot=False, show_plot=True):
    folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/coloc/"
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
        xi, coeff = fit_curve(radius, corr)
        print(xi)
        r_plot = np.linspace(0, np.max(radius), 100)
        ax.plot(r_plot, func(r_plot, xi, coeff), '--', label=r"$\xi=$" + str(round(xi,2)))

    ax.set_xlabel("r (pixels)")
    ax.set_ylabel("C(r)")
    # ax.set_ylim([0,1])
    ax.legend()

    if save_plot == True:
        plate = "1737"
        well = "F4"
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation/"
        plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()

def plot_fit_log(file_list, r_min, r_max, save_plot=False, show_plot=True, x_lim=True):
    folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/coloc/"
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
        xi, coeff = fit_log_line(radius, corr, r_min, r_max)
        print(xi, coeff)
        r_plot = np.linspace(r_min, r_max, 100)
        # r_plot = np.logspace(np.log(10), np.log(50), 100, base=np.exp(1))
        ax.plot(r_plot, func(r_plot, xi, coeff), '--', label=r"$\xi=$" + str(round(xi,2)))

        ax.set_yscale('log')

    ax.set_xlabel("r (pixels)")
    ax.set_ylabel("C(r)")
    # ax.set_ylim([0,1])
    if x_lim == True:
        ax.set_xlim(0,r_max*2)
    ax.legend()

    if save_plot == True:
        plate = "1737"
        well = "F4"
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/correlation_functions/"
        plt.savefig(plot_folder + str(plate) + '_' + well + '_log.png', bbox_inches='tight')

    if show_plot == True:
        plt.show()


file_list = ["AutoCorrelation on VID1737_red_F4_1_t0", "AutoCorrelation on VID1737_red_F4_1_end"]
# file_name = "AutoCorrelation on VID1737_red_F4_1.tif kept stack full"

r_min = 15
r_max = 50
save_plot = True
show_plot = True
x_lim = True

plot_fit_log(file_list, r_min, r_max, save_plot, show_plot, x_lim)