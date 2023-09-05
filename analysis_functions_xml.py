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

    # ## Numpy array format
    # n_tracks = int(bs_data.find("Tracks")["nTracks"])
    # n_frames = int(bs_data.find("particle")["nSpots"])

    tracks = []
    for track in bs_data.find_all("particle"):
        # n_tracks = track["nSpots"]
        # start_frame = int(track.find_all("detection")[0]['t'])
        # last_frame = int(track.find_all("detection")[-1]['t'])

        tracks.append([])
        for frame in track.find_all("detection"):
            t = int(frame['t'])
            tracks[-1].append([t, float(frame['x']), float(frame['y'])])
    
    return tracks


## Squared distance travelled (MSD with no average over t_0 (always set to 0))
def get_sdt(file_path, max_frames=25, plot=False):
    tracks = read_xml(file_path)
    msd = []

    for t in range(max_frames):
        msd_t = []
        for track in tracks:
            n_frames = track[-1][0]-track[0][0]
            if t<n_frames:
                diff = np.array(track[t][1:]) - np.array(track[0][1:])
                msd_t.append(diff[0]**2 + diff[1]**2)
        msd.append(np.mean(msd_t))

    if plot == True:
        plt.loglog(range(max_frames), msd)
        # plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
        # plt.legend()
        plt.xlabel("Time (hours)")
        plt.ylabel("MSD (pixels)")
        plt.show()
    
    return msd

def get_msd(file_path, max_frames=25, plot=False):
    tracks = read_xml(file_path)
    msd = []

    for t in range(max_frames):
        msd_t = []
        for track in tracks:
            n_frames = track[-1][0]-track[0][0]
            if t<n_frames:
                for t0 in range(n_frames-t):
                    diff = np.array(track[t+t0][1:]) - np.array(track[t0][1:])
                    msd_t.append(diff[0]**2 + diff[1]**2)
        msd.append(np.mean(msd_t))

    if plot == True:
        plt.loglog(range(max_frames), msd)
        # plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
        # plt.legend()
        plt.xlabel("Time (hours)")
        plt.ylabel("MSD (pixels)")
        plt.show()
    
    return msd

def plot_sdt_vs_t0(file_path, tau_range, plot_name, max_frames=25, max_plot=15, log=False):

    fig, ax = plt.subplots()
    tracks = read_xml(file_path)

    for tau in tau_range:
        msd = []
        for t0 in range(max_frames):
            msd_t = []
            for track in tracks:
                n_frames = track[-1][0]-track[0][0]
                if tau+t0<n_frames:
                    diff = np.array(track[tau+t0][1:]) - np.array(track[t0][1:])
                    msd_t.append(diff[0]**2 + diff[1]**2)
            msd.append(np.mean(msd_t))

        if log == True:
            ax.loglog(range(max_frames), msd, label=r"$\tau=$" + str(tau))
        else:
            ax.plot(range(max_frames), msd, label=r"$\tau=$" + str(tau))
        # plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
    ax.set_xlabel(r"$t_0$ (hours)")
    ax.set_ylabel("SDT (pixels)")
    ax.set_xlim(right=max_plot)
    ax.legend()
    # plt.show()

    plot_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/plots/sdt_vs_t0"
    plt.savefig(os.path.join(plot_folder, plot_name))
    
    # return msd

## Squared distance travelled (MSD with no average over t_0 (always set to 0))
def get_sdt_vs_t0(file_path, tau, max_frames=25):
    tracks = read_xml(file_path)
    sdt = []

    for t0 in range(max_frames):
        msd_t = []
        for track in tracks:
            n_frames = track[-1][0]-track[0][0]
            if t0+tau<n_frames:
                diff = np.array(track[t0+tau][1:]) - np.array(track[t0][1:])
                msd_t.append(diff[0]**2 + diff[1]**2)
        sdt.append(np.mean(msd_t))
    
    return sdt