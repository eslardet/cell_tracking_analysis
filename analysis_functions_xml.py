import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup as bs
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import os
from scipy import stats

def get_xml_file(plate, well, phase="green"):
    plate = str(plate)
    xml_folder = "/Users/el2021/OneDrive - Imperial College London/PhD/Incucyte/xml/" + plate
    filename = 'VID' + plate + '_' + phase + '_' + well + '_1'
    file_path = os.path.join(xml_folder, filename + '.xml')
    return file_path

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



def get_msd(file_path, max_frames=200, plot=False, remove_outliers=False, return_sd=False):
    tracks = read_xml(file_path)
    msd = []
    sd = []
    for t in range(max_frames):
        msd_t = []
        for track in tracks:
            n_frames = track[-1][0]-track[0][0]
            if t<n_frames:
                for t0 in range(n_frames-t):
                    diff = np.array(track[t+t0][1:]) - np.array(track[t0][1:])
                    msd_t.append(diff[0]**2 + diff[1]**2)
        if remove_outliers == True:
            msd_t_filtered = sorted(msd_t)[:-10]
            mean_t = np.mean(msd_t_filtered)
            msd.append(mean_t)
        else:
            mean_t = np.mean(msd_t)
            sd_t = np.std(msd_t)
            msd.append(mean_t)
            sd.append(sd_t)
        
    if plot == True:
        plt.loglog(range(max_frames), msd)
        # plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
        # plt.legend()
        plt.xlabel("Time (hours)")
        plt.ylabel("MSD (pixels)")
        plt.show()
    
    if return_sd == True:
        return msd, sd
    else:
        return msd

def get_msd_individual_tracks(file_path, max_frames=200):
    tracks = read_xml(file_path)
    msd_all = []
    for track in tracks:
        msd = []
        n_frames = track[-1][0]-track[0][0]
        for t in range(max_frames):
            msd_t = []
            if t<n_frames:
                for t0 in range(n_frames-t):
                    diff = np.array(track[t+t0][1:]) - np.array(track[t0][1:])
                    msd_t.append(diff[0]**2 + diff[1]**2)
            if len(msd_t) > 0:
                mean_t = np.mean(msd_t)
                msd.append(mean_t)
        msd_all.append(msd)
    return msd_all


def get_msd_filter(file_path, max_frames=200, d_threshold=400, d_negative=True, plot=False):
    msd_all = get_msd_individual_tracks(file_path, max_frames=max_frames)
    all_index = []
    for i in range(len(msd_all)):
        msd = msd_all[i]
        n_frames = len(msd)
        if n_frames > 20:
            # max_t = n_frames // 3
            min_frame = 15
            if n_frames < 72:
                max_frame = n_frames
            else:
                max_frame = 72
            d = np.polyfit(np.arange(max_frame-min_frame)/3, msd[min_frame:max_frame], 1)[0]/4
            if d < d_threshold:
                if d_negative == False:
                    if d > 0:
                        all_index.append(i)
                else:
                    all_index.append(i)

    tracks = read_xml(file_path)
    msd = []
    for t in range(max_frames):
        msd_t = []
        # for track in tracks:
        for i in all_index:
            track = tracks[i]
            n_frames = track[-1][0]-track[0][0]
            if t<n_frames:
                for t0 in range(n_frames-t):
                    diff = np.array(track[t+t0][1:]) - np.array(track[t0][1:])
                    msd_t.append(diff[0]**2 + diff[1]**2)
        msd_t_filtered = sorted(msd_t)[:-10]
        mean_t = np.mean(msd_t_filtered)
        msd.append(mean_t)
        
    if plot == True:
        plt.loglog(range(max_frames), msd)
        # plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
        # plt.legend()
        plt.xlabel("Time (hours)")
        plt.ylabel("MSD (pixels)")
        plt.show()
    
    return msd


## Squared distance travelled 
# MSD with no average over t_0 (always set to 0)
# Only use tracks with t_0 <= 2
def get_sdt(file_path, max_frames=200, plot=False):
    tracks = read_xml(file_path)
    sdt = []
    for t in range(max_frames):
        sdt_t = []
        for track in tracks:
            t0 = track[0][0]
            if t0 <= 2:
                n_frames = track[-1][0]-track[0][0]
                if t<n_frames:
                    diff = np.array(track[t][1:]) - np.array(track[0][1:])
                    sdt_t.append(diff[0]**2 + diff[1]**2)
        sdt.append(np.mean(sdt_t))

    if plot == True:
        plt.plot(range(max_frames), sdt)
        # plt.loglog(np.arange(1, 3, 1), 500*np.arange(1, 3, 1), label=r"$\propto t$")
        # plt.legend()
        plt.xlabel("Time (hours)")
        plt.ylabel("MSD (pixels)")
        plt.show()
    
    return sdt

def get_sdt_individual_tracks(file_path, max_frames=200):
    tracks = read_xml(file_path)
    sdt_all = []
    for track in tracks:
        t0 = track[0][0]
        if t0 <= 2:
            sdt = []
            n_frames = track[-1][0]-track[0][0]
            for t in range(max_frames):
                sdt_t = []
                if t<n_frames:
                    diff = np.array(track[t][1:]) - np.array(track[0][1:])
                    sdt_t.append(diff[0]**2 + diff[1]**2)
                if len(sdt_t) > 0:
                    mean_t = np.mean(sdt_t)
                    sdt.append(mean_t)
            sdt_all.append(sdt)
    return sdt_all

def plot_sdt_vs_t0(file_path, tau_range, plot_name, max_frames=200, max_plot=15, log=False):

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