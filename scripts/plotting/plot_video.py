import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
import struct
import argparse

from _plot_common import *

if __name__ == "__main__":
    argParser = argparse.ArgumentParser(prog="PlotVideo")
    argParser.add_argument("filepath")
    argParser.add_argument("--out-dir", dest="outputdir")
    argParser.add_argument("-o", action='store_true')
    argParser.add_argument("-s", action="store_true")

    args = vars(argParser.parse_args())

    # pixel in inches
    px = 1/plt.rcParams['figure.dpi']

    with open(args["filepath"], "rb") as file:
        plot_counts = unpack_ull(file)

        title   = read_string(file)
        x_label = read_string(file)
        y_label = read_string(file)
        x_scale = read_string(file)
        y_scale = read_string(file)

        fig, axs = plt.subplots(ncols=plot_counts, figsize=(1920*px, 1080*px))
        fig.suptitle(title)

        meta_args = {}
        meta_data_counts = unpack_ull(file)

        for i in range(meta_data_counts):
            name = read_string(file)
            dtype = fp_type_from_size(struct.unpack("<B", file.read(1))[0])
            meta_args[name] = np.fromfile(file, dtype=dtype, count=1)[0]

        frames_counts = []
        shapes = []
        frames = []

        for ax in axs.flatten():
            value_type_size, space_type_size = struct.unpack("<BB", file.read(2))
            value_type = fp_type_from_size(value_type_size)
            #time_type = fp_type_from_size(time_type_size)
            space_type = fp_type_from_size(space_type_size)

            sub_title = read_string(file)
            ax.set_title(sub_title)
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_xscale(x_scale)
            ax.set_yscale(y_scale)

            scale = np.fromfile(file, space_type, 2)
            shape = np.fromfile(file, np.int64, 2)
            shapes.append(shape)

            interval = struct.unpack("<I", file.read(4))[0]

            frame_counts = unpack_ull(file)
            frames_counts.append(frame_counts)

            frames.append(np.fromfile(file, value_type, frame_counts * np.prod(shape)).reshape(-1, *shape))

        times = np.fromfile(file, dtype=np.float64, count=np.max(frames_counts))

        time_txt = fig.text(.5, 0.85, f"Time = {0:.3f}", fontsize=14, transform=fig.transFigure, ha="center")

        images = []
        for ax, shape in zip(axs, shapes):
            images.append(ax.imshow(np.zeros(shape), cmap="seismic", interpolation=None, origin="lower"))

        for ax, image in zip(axs, images):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(image, cax=cax, orientation='vertical')

        def animate(i):
            vmin = np.inf
            vmax = -np.inf

            for image_index, image in enumerate(images):
                if(i < frames_counts[image_index]):
                    frame = frames[image_index][i]
                    vmin = np.min((vmin, np.min(frame)))
                    vmax = np.max((vmax, np.max(frame)))
                    image.set_data(frame)

            for image in images:
                image.set_clim(vmin, vmax)

            time_txt.set_text(f"Time = {times[i]:.3f}")

        anim = FuncAnimation(fig, animate, interval=16, frames=np.max(frames_counts), repeat_delay=2000)        

        if args["outputdir"] is not None:
            anim.save(args["outputdir"])
        elif args["o"] == True:
            filepath = args["filepath"]
            ext_pos = filepath.rfind('.')

            anim.save(filepath[:ext_pos] + ".mp4")

        if args["s"] == True:
            plt.show()