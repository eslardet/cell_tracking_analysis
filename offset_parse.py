import numpy as np
import os, csv
import matplotlib.pyplot as plt


def get_wells_time(plate, axis):
    offset_folder = "/Volumes/T7 Shield/Incucyte_data/offset/"
    file = offset_folder + str(plate) + "_" + axis + '_offset.txt'

    with open(file) as f:
        reader = csv.reader(f, delimiter="\n")
        r = list(reader)

    all_wells = r[6][0].split('\t')[2:]
    frames = len(r[7:])
    r = r[7:]
    t = []
    for f in range(frames):
        t.append(float(r[f][0].split('\t')[1]))

    return all_wells, t


def get_offsets(plate, axis, well):
    offset_folder = "/Volumes/T7 Shield/Incucyte_data/offset/"

    file = offset_folder + str(plate) + "_" + axis + '_offset.txt'

    with open(file) as f:
        reader = csv.reader(f, delimiter="\n")
        r = list(reader)

    all_wells = r[6][0].split('\t')[2:]
    well_index = all_wells.index(well)
    frames = len(r[7:])
    r = r[7:]

    t = []
    offsets = []
    for f in range(frames):
        t.append(float(r[f][0].split('\t')[1]))
        offsets.append([float(i) for i in r[f][0].split('\t')[2:]])

    # converts offset to lists for each well instead
    offsets = list(map(list, zip(*offsets)))

    return t, offsets[well_index]

def plot_offsets(plate, well):
    fig, ax = plt.subplots(2,1)
    time, offset = get_offsets(plate=plate, axis='x', well=well)
    ax[0].plot(time, offset, '-x')
    ax[0].set_title('x offset')
    time, offset = get_offsets(plate=plate, axis='y', well=well)
    ax[1].plot(time, offset, '-x')
    ax[1].set_title('y offset')

    for a in ax:
        a.set_xlabel('time (hours)')
        a.set_ylabel('offset')
    fig.tight_layout()

    plot_folder = "/Volumes/T7 Shield/Incucyte_data/offset_plots/"
    plt.savefig(plot_folder + str(plate) + '_' + well + '.png', bbox_inches='tight')
    plt.close()
    # plt.show()

# plot_offsets(plate=1727, well='B3')
for plate in [1727]:
    for axis in ['x', 'y']:
        wells, time = get_wells_time(plate, axis)
        for well in wells:
            plot_offsets(plate, well)


## To pre-process offsets
# def save_offsets

# plate = 1722
# well = 'B2'

# folder = "/Volumes/T7 Shield/Incucyte_data/offset_processed/"
# save_offset = open(folder + str(plate) + '_' + well + '.txt', 'w')

# t, x_offset = get_offsets(plate=plate, axis='x', well=well)
# t, y_offset = get_offsets(plate=plate, axis='y', well=well)

# ref_x = x_offset[1]
# ref_y = y_offset[1]

# for frame in range(20):
#     print(frame, x_offset[frame]-ref_x, y_offset[frame]-ref_y)