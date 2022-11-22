import numpy as np
import matplotlib.pyplot as plt
import struct

import sys

def unpack_ull(file):
    return struct.unpack("<Q", file.read(8))[0]


def read_string(file):
    str_len = unpack_ull(file)
    str = file.read(str_len).decode('UTF-8')

    return str    


def fp_type_from_size(size):
    if size == 8:
        return np.float64
    elif size == 4:
        return np.float32
    else:
        raise ValueError(f"{size}-byte floating point values are not currently supported!")


if __name__ == "__main__":
    if(len(sys.argv) != 2):
        print("Unexpected argument count! (Expected 2).")
        sys.exit()

    with open(sys.argv[1], 'rb') as file:
        title   = read_string(file)
        x_label = read_string(file)
        y_label = read_string(file)
        x_scale = read_string(file)
        y_scale = read_string(file)

        line_counts = unpack_ull(file)

        plt.figure()
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xscale(x_scale)
        plt.yscale(y_scale)

        for line_index in range(line_counts):
            label = read_string(file)
            colour = read_string(file)
            linestyle = read_string(file)

            x_type_size, y_type_size = struct.unpack("<BB", file.read(2))
            x_type = fp_type_from_size(x_type_size)
            y_type = fp_type_from_size(y_type_size)

            samples_count = unpack_ull(file)

            x = np.fromfile(file, dtype=x_type, count=samples_count)
            y = np.fromfile(file, dtype=y_type, count=samples_count)

            plt.plot(x, y, label=label, c=colour, ls=linestyle)

        plt.grid(True, which="both", ls=":", color='0.65')
        plt.legend()
        plt.show()