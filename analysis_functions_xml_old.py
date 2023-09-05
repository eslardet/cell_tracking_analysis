import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup as bs
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import os


def read_xml(file_path):
    with open(file_path, 'r') as f:
        data = f.readlines()
    data = "".join(data)
    bs_data = bs(data, features="xml")

    ## Numpy array format
    n_tracks = int(bs_data.find("Tracks")["nTracks"])
    n_frames = int(bs_data.find("particle")["nSpots"])

    tracks = np.zeros((n_frames, n_tracks, 2))

    track_num = 0
    for track in bs_data.find_all("particle"):
        start_frame = int(track.find_all("detection")[0]['t'])
        last_frame = int(track.find_all("detection")[-1]['t'])
        tracks[:start_frame, track_num, :] = None
        tracks[last_frame:, track_num, :] = None
        for frame in track.find_all("detection"):
            t = int(frame['t'])
            tracks[t, track_num, :] = np.c_[float(frame['x']), float(frame['y'])]
        track_num += 1
    
    return tracks

def animate_xml(tracks, x_max=1408, y_max=1040):

    plt.rcParams["animation.html"] = "jshtml"
    plt.ioff()
    plt.rcParams['animation.embed_limit'] = 2**128

    fig, ax = plt.subplots(figsize=(5,5))

    points, = plt.plot([], [], 'o', zorder=1)
    lines, = plt.plot([], [], '-')

    def init():
        ax.set_xlim(0, x_max)
        ax.set_ylim(0, y_max)
        return points, lines,

    def update(n):
        points.set_data(tracks[n,:,0], y_max-tracks[n,:,1])
        for i in range(np.shape(tracks)[1]):
            # lines.set_data(tracks[:n,i,0], y_max-tracks[:n,i,1])
            plt.plot(tracks[:n+1,i,0], y_max-tracks[:n+1,i,1], '-', color="gray", lw=0.5)
        ax.set_title("t = " + str(n), fontsize=10, loc='left')
        
        return points, lines,


    ani = FuncAnimation(fig, update, init_func=init, frames=len(tracks), interval=10, blit=True)
    # plt.show()
    folder = os.path.abspath('../animations_cell')
    filename = 'test' + '.mp4'
    if not os.path.exists(folder):
        os.makedirs(folder)
    ani.save(os.path.join(folder, filename))


def plot_msd(tracks):
    msd = []
    n_frames = len(tracks)
    n_tracks = np.shape(tracks)[1]

    for t in range(n_frames):
        diff = tracks[t,:] - tracks[0,:]
        msd.append(np.mean(diff[:,0]**2 + diff[:,1]**2))

    # print(msd[:5])

    plt.loglog(range(n_frames), msd)
    plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
    plt.legend()
    plt.xlabel("Time (hours)")
    plt.ylabel("MSD (pixels)")
    plt.show()