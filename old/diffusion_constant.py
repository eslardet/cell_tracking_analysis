import os
import analysis_functions_tracks as fun
import numpy as np
import time
from matplotlib import pyplot as plt

save_plot=True
plate = 'VID1713'
plot_name = plate + "_short.png"
show_plot=False
max_t = 3

# plate = "1712"
# well = "B2"
# pic = "1"

# file_path = "tracks/VID" + plate + "_day_red_" + well + "_" + pic + "_Tracks_bright.xml"

plate_dict = {'VID1711': "Astrocytes & T-cells", 'VID1712': "Microglia & T-cells", 'VID1713': "T-cells only"}
number_dict = {2:'low PVL, not stim. \n', 3: 'low PVL, stim. \n', 4: 'high PVL1, not stim. \n', 5: 'high PVL1, stim. \n', 6: 'high PVL2, not stim. \n', 7: 'high PVL2, stim. \n', 8: 'HAM1, not stim. \n', 9: 'HAM1, stim. \n', 10: 'HAM2, not stim. \n', 11: 'HAM2, stim. \n'}
cell_dict = {'B': 'non-specific CD4', 'C': 'specific CD4', 'D': 'not-specific CD8', 'E': 'specific CD8'}
cell_dict_tcell = {'B': 'non-specific CD4', 'C': 'S CD4 + non-specific CD8', 'D': 'specific CD8'}

fig, ax = plt.subplots(figsize=(15,8))


diff_const = []
treatment = []

for cell in ['B', 'D']:
    for number in np.arange(2,12):
        well = cell + str(number)

        xml_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml_short"
        phase = 'red'
        # well = 'B10'
        pic = '1'

        filename = plate + '_' + phase + '_' + well + '_' + pic
        file_path = os.path.join(xml_folder, filename + '.xml')

        # tracks = fun.read_xml(file_path)

        msd = fun.get_msd(file_path, plot=False)
        diff_const.append(np.polyfit(np.arange(max_t), msd[:max_t], 1)[0]/4)
        if plate == "VID1713":
            treatment.append(number_dict[number] + cell_dict_tcell[cell])
        else:
            treatment.append(number_dict[number] + cell_dict[cell])

ax.scatter(treatment, diff_const)
ax.set_ylabel("Diffusion constant (pixels/ hour)")
plt.xticks(fontsize=7,rotation=45)
ax.set_title("T-cells only")

if save_plot == True:
    plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/diffusion_const"
    plt.savefig(os.path.join(plot_folder, plot_name))

if show_plot == True:
    plt.show()