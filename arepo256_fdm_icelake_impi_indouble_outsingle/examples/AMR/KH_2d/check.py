#!/usr/bin/env python3
""" Idealized test, AMR Kelvin-Helmholtz (non-static) in 2D. """
import os.path
import h5py
import numpy as np
import sys


def verify_result(path):
    """ Unit test function, validate results. """
    snap_path = '%s/output/snap_004.hdf5' % path

    result = {}
    if not os.path.isfile(snap_path):
        return False, ['Expected snapshot file apparently not produced.']

    with h5py.File(snap_path, 'r') as f:
        for k in f['PartType0'].keys():
            result[k] = f['PartType0'][k][()]

    solution_dens = [0.952269494, 2.09526524]  # min max

    if np.abs(result['Density'].min() - solution_dens[0]) > 1e-5:
        return False, [
            'Density field mismatch with recorded solution (min) in the past.'
        ]
    if np.abs(result['Density'].max() - solution_dens[1]) > 1e-5:
        return False, [
            'Density field mismatch with recorded solution (max) in the past.'
        ]

    # return success
    return True, []


if __name__ == '__main__':
    path = sys.argv[1]
    print("examples/AMR/KH_2d/check.py: checking output in directory " + path)

    status_ok, info = verify_result(path=path)

    if status_ok:
        print("check.py: success!")
        sys.exit(0)  # normal exit

    for msg in info:
        print('check.py: failed! ' + msg)
    sys.exit(1)


def visualize_result(path, Lx, Ly):
    pass
