import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/coloc/"

filename = "Correlation of images_end.csv"
file_list = ["Correlation of images_t0.csv", "Correlation of images_end.csv"]


def get_exponent(r_plot, corr, min_r=5, max_r=50):
    idx1 = np.where(r_plot<max_r)[0]
    idx2 = np.where(r_plot>min_r)[0]
    idx = list(set(idx1) & set(idx2))
    # print(corr_bin_av[idx])

    exponent = np.polyfit(x=np.log10(r_plot[idx[0]:idx[-1]]), y=np.log10(np.abs(corr[idx[0]:idx[-1]])), deg=1)[0]
    return exponent

def func(r, alpha, coeff, C_inf):
    return r**(-alpha)*np.exp(coeff) + C_inf

def fit_curve(r_plot, corr):
    alpha, coeff, p_inf = curve_fit(func, r_plot, corr)[0]
    return alpha, coeff, p_inf

fig, ax = plt.subplots()
for filename in file_list:
    with open(folder + filename) as f:
        reader = csv.reader(f, delimiter=",")
        r = np.array(list(reader))

    # print(type(r[6,1]))
    pixels = r[1:,1].astype(float)
    corr = r[1:,2].astype(float)
    corr = corr/corr[0]


    ## Binning
    corr_r_max = np.max(pixels)
    # corr_r_max = 100
    corr_r_min = 0
    r_bin_num = 200

    bin_size = (corr_r_max-corr_r_min) / r_bin_num
    r_plot = np.linspace(corr_r_min, corr_r_max, num=r_bin_num, endpoint=False) + bin_size/2

    corr_plot = []
    for i in range(r_bin_num):
        lower = r_plot[i]
        try:
            upper = r_plot[i+1]
        except:
            upper = corr_r_max+1
        idx = np.where((pixels>lower)&(pixels<upper))[0]
        if len(idx) != 0:
            # c = np.mean(corr[idx])
            corr_plot.append(np.mean(corr[idx]))

    ax.plot(r_plot, corr_plot, label=filename)
    alpha, coeff, C_inf = fit_curve(r_plot, corr_plot)
    ax.plot(r_plot, func(r_plot, alpha, coeff, C_inf), '--', label=filename + " fit")
    print(1/alpha)


# ax.set_xlim()
ax.legend()
plt.show()