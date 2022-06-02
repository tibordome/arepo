#!/usr/bin/env python3
""" Idealized test, AMR sedov in 3D. """
import os.path
import h5py
import numpy as np
import sys


def verify_result(path):
    """ Unit test function, validate results. """
    snap_path = '%s/output/snap_030.hdf5' % path

    result = {}
    if not os.path.isfile(snap_path):
        return False, ['Expected snapshot file apparently not produced.']

    with h5py.File(snap_path, 'r') as f:
        for k in f['PartType0'].keys():
            result[k] = f['PartType0'][k][()]

    solution_dens = [0.01842026, 15.53200192]  # min max

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
    status_ok, info = verify_result(path=sys.argv[1])

    if status_ok:
        print("check.py: success!")
        sys.exit(0)  # normal exit

    for msg in info:
        print('check.py: failed! ' + msg)
    sys.exit(1)


def visualize_result(path, Lx, Ly, Lz=1.0):
    """ Derive and plot a series of 1D radial density profiles. """
    import matplotlib.pyplot as plt
    from scipy.stats import binned_statistic

    # radial profile at a few times
    snaps = [5, 10, 15, 20, 25, 30]
    ptName = 'PartType0'

    fig = plt.figure(figsize=[14, 10])
    ax = fig.add_subplot(111)
    ax.set_xlabel('radius')
    ax.set_ylabel('density')

    for i, snap in enumerate(snaps):
        # load
        data = {}
        with h5py.File(path + 'output/snap_' + '%03d.hdf5' % snap, 'r') as f:
            for k in f[ptName].keys():
                data[k] = f[ptName][k][()]
            data['Time'] = f['Header'].attrs['Time']

        # plot
        rad = np.sqrt( (data['Coordinates'][:,0] - Lx/2)**2 + \
                       (data['Coordinates'][:,1] - Ly/2)**2 + \
                       (data['Coordinates'][:,2] - Lz/2)**2 )

        nBins = int(np.sqrt(data['Density'].size)) * 2
        dens, r_pts, _ = binned_statistic(rad,
                                          data['Density'],
                                          statistic='mean',
                                          bins=nBins,
                                          range=[0, Lx])

        ax.plot(r_pts[:-1] + Lx / nBins / 2,
                dens,
                linestyle='-',
                label='t = %.2f' % data['Time'])

        #plt.plot(t[:,1],t[:,2], label="analytical solution",zorder=0) # todo

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(os.path.join(path, "vis/sedov_profile_3d.pdf"))
