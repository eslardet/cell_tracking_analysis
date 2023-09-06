import numpy as np
from analysis_functions_xml import *

# file_path = get_xml_file(plate='1727', well='B2')

# tracks = read_xml(file_path)
# msd = []
# count = 0


# for track in tracks:
#     t0 = track[0][0]
#     if t0 <= 2:
#         count += 1
#         for spot in track:
#             tau = spot[0]
#             # print(tau)
# print(count)
# print(len(tracks))

# sdt = []
# for t in range(200):
#     sdt_t = []
#     for track in tracks:
#         t0 = track[0][0]
#         if t0 <= 2:
#             n_frames = track[-1][0]-track[0][0]
#             if t<n_frames:
#                 diff = np.array(track[t][1:]) - np.array(track[0][1:])
#                 sdt_t.append(diff[0]**2 + diff[1]**2)
#     sdt.append(np.mean(sdt_t))

# plt.plot(np.arange(len(sdt)), sdt)
# plt.show()

