#!/usr/bin/env python3
""" Idealized test, AMR sedov in 2D. """
import os.path
import h5py
import numpy as np
import sys


def verify_result(path):
    """ Unit test function, validate results. """
    theta = np.deg2rad(57)

    def analytic_grad_rho(x, y):
        return (np.cos(theta) * np.ones_like(x),
                np.sin(theta) * np.ones_like(x), np.zeros_like(x))

    # load
    snap_path = '%s/output/snap_001.hdf5' % path

    result = {}
    if not os.path.isfile(snap_path):
        return False, ['Expected snapshot file apparently not produced.']

    with h5py.File(snap_path, 'r') as f:
        for k in f['PartType0'].keys():
            result[k] = f['PartType0'][k][()]

    # discard border cells
    x, y, z = result['Coordinates'].T

    margin = 0.1
    inside = (x > margin) & (x < 1 - margin) & (y > margin) & (y < 1 - margin)

    gx, gy, gz = result['DensityGradient'][inside, :].T

    # compare to reference gradient (analytic)
    grefx, grefy, grefz = analytic_grad_rho(x[inside], y[inside])

    err = np.sqrt((gx - grefx)**2 + (gy - grefy)**2 + (gz - grefz)**2)
    err_max = np.max(err)
    err_L2 = np.sqrt(np.mean(err**2))

    error_msgs = []

    if err_max > 1e-12:
        error_msgs.append('Absolute error too large, maximum of %g.' % err_max)
    if err_L2 > 1e-13:
        error_msgs.append('L2 error too large, mean of %g.' % err_L2)

    if len(error_msgs):
        return False, error_msgs

    # return success
    return True, []


if __name__ == '__main__':
    path = sys.argv[1]
    print("examples/AMR/gradients_2d/check.py: checking output in directory " +
          path)

    status_ok, info = verify_result(path=path)

    if status_ok:
        print("check.py: success!")
        sys.exit(0)  # normal exit

    for msg in info:
        print('check.py: failed! ' + msg)
    sys.exit(1)


def visualize_result(path, Lx, Ly):
    """ Needs cleanup. """
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 12))
    plt.axvline(0.5)
    plt.quiver(x[inside],
               y[inside],
               gx - np.cos(theta),
               gy - np.sin(theta),
               scale=5)

    plt.tight_layout()
    plt.savefig(path + 'vis/gradients.pdf')
