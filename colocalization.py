import numpy as np
# from analysis_functions_xml import *
from analysis_functions_image import *
from scipy import ndimage as ndi
from skimage import measure, segmentation, filters
from skimage.measure import manders_coloc_coeff, intersection_coeff, pearson_corr_coeff

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)

def get_manders_coeff(plate, well, slice, green_thresh=0.25, red_thresh=0.15):
    im_red = get_image(plate, well, "red")[slice]
    im_green = get_image(plate, well, "green")[slice]

    # ## First segment with thresholding
    im_red_thresh = im_red/np.max(im_red) > red_thresh
    im_green_thresh = im_green/np.max(im_green) > green_thresh

    coeff = measure.manders_coloc_coeff(im_green_thresh, im_red_thresh)
    return coeff

def get_manders_increase(plate, well, compare, green_thresh=0.25, red_thresh=0.15):
    coeff0 = get_manders_coeff(plate, well, compare[0], green_thresh, red_thresh)
    coeff1 = get_manders_coeff(plate, well, compare[1], green_thresh, red_thresh)
    
    percentage_increase = (coeff1-coeff0)/coeff0*100

    return percentage_increase

def plot_manders_vs_time(plate, well, slice_range, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True):
    green_path = get_tif_file(plate, well, "green")
    red_path = get_tif_file(plate, well, "red")
    if not os.path.exists(green_path) or not os.path.exists(red_path):
        print("File does not exist: " + str(plate) + ", " + str(well))
        return
    coeff = []
    for i in slice_range:
        if time_av == True:
            av = 0
            for j in range(10):
                av += get_manders_coeff(plate, well, i+j, green_thresh, red_thresh)
            coeff.append(av/10)
        else:
            coeff.append(get_manders_coeff(plate, well, i, green_thresh, red_thresh))

    fig,ax = plt.subplots()

    ax.plot(slice_range/3, coeff, 'o-')
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Manders' colocalization coefficient")
    ax.set_ylim([0,1])

    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization/" + str(plate) + "/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')
    if show_plot == True:
        plt.show()
    plt.close()

def plot_manders_increase(cell_type, group_list, stim_list, t_cell_list, slice_compare, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True):
    plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization_increase/"
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)
    plot_name = cell_type + "_" + stim_list[-1] + "_" + t_cell_list[-1] + "_" + str(slice_compare[1]) + ".png"
    if os.path.exists(os.path.join(plot_folder, plot_name)):
        print("Already plotted!")
        print(os.path.join(plot_folder, plot_name))
        return
    
    fig, ax = plt.subplots(figsize=(12,6))

    x_ticks = []

    plate_list = get_plates(cell_type)
    pos = 1
    for group in group_list:
        for stim in stim_list:
            for t_cell in t_cell_list:
                well_list = get_wells(group, stim, t_cell)
                for plate in plate_list:
                    for well in well_list:
                        green_path = get_tif_file(plate, well, "green")
                        red_path = get_tif_file(plate, well, "red")
                        if os.path.exists(green_path) and os.path.exists(red_path):
                    
                            percentage_increase = get_manders_increase(plate, well, slice_compare, green_thresh, red_thresh)
                            if percentage_increase > 0:
                                ax.scatter(pos, percentage_increase, color='k')
                                ax.annotate(str(plate) + ", " + well, (pos, percentage_increase))
                x_ticks.append(group + "\n " + stim + "\n " + t_cell)
                pos += 1
    ax.set_ylabel("Percentage increase in Manders' Coefficient (%)")
    ax.set_ylim(0,800)
    set_axis_style(ax,x_ticks)

    ax.set_title(cell_type + ", compare slices " + str(slice_compare[0]) + " and " + str(slice_compare[1]))

    if save_plot == True:
        plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/colocalization_increase/"
        if not os.path.exists(plot_folder):
            os.makedirs(plot_folder)
        plot_name = cell_type + "_" + stim + "_" + t_cell + "_" + str(slice_compare[1]) + ".png"
        plt.savefig(os.path.join(plot_folder, plot_name))

    if show_plot == True:
        plt.show()

    plt.close()


plate = 1728
well = "B2"

slice_range = np.arange(1, 202, 10)

# t0 = time.time()
# plot_manders_vs_time(plate, well, slice_range)
# print(time.time()-t0)

# for well in ["F2", "F3", "F4", "F5", "G2", "G3", "G4", "G5"]:
# for well in all_wells:
#     t0 = time.time()
#     plot_manders_vs_time(plate, well, slice_range)
#     # print(get_manders_coeff(plate, well, -1))
#     print(time.time()-t0)


cell_type = 'Microglia & T-cells'
group_list = ["Uninfected healthy control", "AC lPVL", "AC hPVL", "HAM"]
# group_list = ["Uninfected healthy control"]
stim_list = ["With stimulation"]
t_cell_list = ["Specific CD4"]

# slice_compare = [1, 72]
# slice_compare = [1, 215]
slice_compare = [1, 36]

# t0 = time.time()
# plot_manders_increase(cell_type, group_list, stim_list, t_cell_list, slice_compare, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True)
# print(time.time()-t0)

for cell_type in ["Astrocytes & T-cells", "Microglia & T-cells"]:
    t0 = time.time()
    for stim in all_stim:
        for t_cell in all_t_cell:
            stim_list = [stim]
            t_cell_list = [t_cell]
            for slice_compare in [[1, 36]]:
                t0 = time.time()
                plot_manders_increase(cell_type, group_list, stim_list, t_cell_list, slice_compare, green_thresh=0.25, red_thresh=0.15, time_av=False, show_plot=False, save_plot=True)
                print(time.time()-t0)

# plate = 1728
# well = "B2"
# slice = 1
# im_red = get_image(plate, well, "red")[slice]
# im_green = get_image(plate, well, "green")[slice]
# red_thresh = 0.15
# # ## First segment with thresholding
# im_red_thresh = im_red/np.max(im_red) > red_thresh

# io.imshow(im_red_thresh)
# plt.show()