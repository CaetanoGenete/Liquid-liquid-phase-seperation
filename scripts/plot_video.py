import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

import sys
import struct

from _plot_common import *

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Unexpected argument count! (Expected 2)."

    with open(sys.argv[1], 'rb') as file:
        title   = read_string(file)
        x_label = read_string(file)
        y_label = read_string(file)
        x_scale = read_string(file)
        y_scale = read_string(file)

        fig = plt.figure()
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xscale(x_scale)
        plt.yscale(y_scale)

        #unused
        plot_counts = unpack_ull(file)
        assert plot_counts == 1, "Only one video is allowed to be displayed!"

        value_type_size, space_type_size = struct.unpack("<BB", file.read(2))
        value_type = fp_type_from_size(value_type_size)
        #time_type  = fp_type_from_size(time_type_size)
        space_type = fp_type_from_size(space_type_size)

        min_max = np.fromfile(file, value_type, 2)
        vmin = min_max[0]
        vmax = min_max[1]
        if((min_max == 0).all()):
            vmin = None
            vmax = None
            
        scale = np.fromfile(file, space_type, 2)
        shape = np.fromfile(file, np.int64, 2)

        interval = struct.unpack("<I", file.read(4))[0]

        frames_count = unpack_ull(file)

        frames = []
        for i in range(frames_count):
            x = np.fromfile(file, value_type, np.prod(shape)).reshape(shape)
            frames.append([plt.imshow(x, vmin=vmin, vmax=vmax)])

        anim = ArtistAnimation(fig, frames, interval=interval, blit=True)

        plt.colorbar()
        plt.show()

