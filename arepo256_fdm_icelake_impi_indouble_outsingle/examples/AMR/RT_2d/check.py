#!/usr/bin/env python3
""" Idealized test, AMR sedov in 2D. """
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

    # check
    return False, ['verification of result needed (vs analytical solution)']


if __name__ == '__main__':
    path = sys.argv[1]
    print("examples/AMR/sedov_2d/check.py: checking output in directory " +
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
    import matplotlib.pyplot as p
    import arepo

    fig = p.figure(figsize=[14, 1.4 * 14])

    ax1 = p.subplot2grid((17, 5), (0, 0), rowspan=8)
    ax2 = p.subplot2grid((17, 5), (0, 1), rowspan=8)
    ax3 = p.subplot2grid((17, 5), (0, 2), rowspan=8)
    ax4 = p.subplot2grid((17, 5), (0, 3), rowspan=8)
    ax5 = p.subplot2grid((17, 5), (0, 4), rowspan=8)

    ax6 = p.subplot2grid((17, 5), (8, 0), rowspan=8)
    ax7 = p.subplot2grid((17, 5), (8, 1), rowspan=8)
    ax8 = p.subplot2grid((17, 5), (8, 2), rowspan=8)
    ax9 = p.subplot2grid((17, 5), (8, 3), rowspan=8)
    ax10 = p.subplot2grid((17, 5), (8, 4), rowspan=8)

    cbar = p.subplot2grid((17, 5), (16, 0), colspan=4)
    cbar2 = p.subplot2grid((17, 5), (16, 4), colspan=4)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]
    nums = [0, 5, 10, 15, 20, 25, 30, 35, 40]

    cmin = 0.9
    cmax = 2.1

    for j in np.arange(len(nums)):
        sn = arepo.Simulation(path + "output/snap_%.3d.hdf5" % nums[j])

        sn.plot_AMRslice(sn.rho,
                         axes=axes[j],
                         colorbar=False,
                         cmap='parula',
                         res=2048,
                         vmin=cmin,
                         vmax=cmax,
                         gradient=sn.grar)
        axes[j].xaxis.set_visible(False)
        axes[j].yaxis.set_visible(False)
        axes[j].set_frame_on(False)
        axes[j].annotate("t = %.2g" % sn.time, [0.02, 1.87], color='black')

    arr = np.arange(0, 256)
    arr = np.array([
        arr,
    ])
    cbar.imshow(arr,
                aspect='auto',
                extent=(cmin, cmax, 0, 1),
                cmap='parula',
                rasterized=True,
                origin="lower")
    cbar.set_frame_on(True)
    cbar.set_xlabel("$\\rho$")
    cbar.set_yticks([])

    j = 9
    sn = arepo.Simulation(path + "output/snap_%.3d.hdf5" % 40)

    sn.plot_AMRslice(sn.amrlevel,
                     axes=axes[j],
                     colorbar=False,
                     cmap='RdBu',
                     res=2048)
    axes[j].xaxis.set_visible(False)
    axes[j].yaxis.set_visible(False)
    axes[j].set_frame_on(False)
    axes[j].annotate("t = %.2g" % sn.time, [0.02, 1.87], color='black')

    minlevel = sn.amrlevel.min()
    maxlevel = sn.amrlevel.max()

    arr = np.arange(0, 256)
    arr[:128] = 0.
    arr[128:] = 1.
    arr = np.array([
        arr,
    ])
    cbar2.imshow(arr,
                 aspect='auto',
                 extent=(minlevel, maxlevel, 0, 1),
                 cmap='RdBu',
                 rasterized=True,
                 origin="lower")
    cbar2.set_frame_on(True)
    cbar2.set_xlabel("AMR level")
    cbar2.set_yticks([])
    cbar2.set_xticks([minlevel + 0.25, maxlevel - 0.25])
    cbar2.set_xticklabels(["%d" % minlevel, "%d" % maxlevel])

    p.tight_layout(h_pad=1., w_pad=1.)
    p.savefig(path + "vis/RT_time.pdf")
