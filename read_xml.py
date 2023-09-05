from bs4 import BeautifulSoup as bs
import lxml
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import os

fileName = "FakeTracks2.xml"

with open(fileName, 'r') as f:
    data = f.readlines()
data = "".join(data)
bs_data = bs(data, features="xml")


## Numpy array format
n_frames = 50

n_tracks = int(bs_data.find("Tracks")["nTracks"])

tracks = np.zeros((n_frames, n_tracks, 2))

track_num = 0
for track in bs_data.find_all("particle"):
    start_frame = int(track.find_all("detection")[0]['t'])
    last_frame = int(track.find_all("detection")[-1]['t'])
    tracks[:start_frame, track_num, :] = None
    tracks[last_frame:, track_num, :] = None
    for frame in track.find_all("detection"):
        t = int(frame['t'])
        # print(np.c_[float(frame['x']), float(frame['y'])])
        tracks[t, track_num, :] = np.c_[float(frame['x']), float(frame['y'])]
    track_num += 1

print(tracks[44,:,:])

# ## List format
# tracks = []
# for track in bs_data.find_all("particle"):
#     pos = []
#     for frame in track.find_all("detection"):
#     # for t in range(num_frame):
#         # frame = track.find_all("detection")
#         pos.append([float(frame['t']), float(frame['x']), float(frame['y'])])
#     tracks.append(pos)



## Animation
# track_0 = tracks[0]

# plt.rcParams["animation.html"] = "jshtml"
# plt.ioff()
# plt.rcParams['animation.embed_limit'] = 2**128

# fig, ax = plt.subplots(figsize=(5,5))

# points, = plt.plot([], [], 'o', zorder=1)
# # ax.set_xlim(0, 128)
# # ax.set_ylim(0, 128)
# # plt.show()


# def init():
#     ax.set_xlim(0, 128)
#     ax.set_ylim(0, 128)
#     return points,

# def update(n):
#     points.set_data(tracks[n,:,0], tracks[n,:,1])
#     ax.set_title("t = " + str(n), fontsize=10, loc='left')
    
#     return points,


# ani = FuncAnimation(fig, update, init_func=init, frames=len(tracks), interval=10, blit=True)
# plt.show()
# # folder = os.path.abspath('../animations_cell')
# # filename = 'test' + '.mp4'
# # if not os.path.exists(folder):
# #     os.makedirs(folder)
# # ani.save(os.path.join(folder, filename))

# # ani.save('test.mp4')
# # plt.show()