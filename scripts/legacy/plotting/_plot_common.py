import struct
import numpy as np

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