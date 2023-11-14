import numpy as np
# from analysis_functions_xml import *
from analysis_functions_image import *


plate = 1737
well = "E5"
# frames = np.arange(1,108,5)
frames = np.arange(1,216,10)

fig, ax = plt.subplots()

phase = "red"
im = get_image(plate, well, phase)
coverage = []
for i in frames:
    a = im[i]/np.max(im[i]) > 0.15
    coverage.append(a.sum()/a.size * 100)
times = frames / 3
coverage = np.array(coverage)
ax.plot(times, coverage, 'o-', label="Microglia", color="tab:red")

phase = "green"
im = get_image(plate, well, phase)
coverage = []
for i in frames:
    a = scaled_threshold(im, threshold=0.75)
    coverage.append(a.sum()/a.size * 100)
times = frames / 3
coverage = np.array(coverage)
ax.plot(times, coverage, 'o-', label="T-cells", color="tab:green")

ax.set_xlabel("Time (hours)")
ax.set_ylabel("Coverage (%)")
ax.legend()
# ax.set_ylim([0,100])

plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/coverage/"
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')

plt.show()
