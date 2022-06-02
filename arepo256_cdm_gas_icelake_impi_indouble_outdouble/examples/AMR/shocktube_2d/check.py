#!/usr/bin/env python3
""" Idealized test, AMR shocktube in 2D. """
import h5py
import numpy as np
import sys


def verify_result(path):
    """ Unit test function, validate results. """
    snap_path = '%s/output/snap_010.hdf5' % path

    result = {}

    with h5py.File(snap_path, 'r') as f:
        for k in f['PartType0'].keys():
            result[k] = f['PartType0'][k][()]

    solution = {
        'Density': [0.11534932, 1.03735481],  # min max of each at t=0.1
        'InternalEnergy': [1.52313356, 3.56966442],
        'Pressure': [0.08976760, 1.05274379]
    }

    for key in solution:
        if np.abs(result[key].min() - solution[key][0]) > 1e-5:
            return False, [
                '%s field mismatch with recorded solution (min) in the past.' %
                key
            ]
        if np.abs(result[key].max() - solution[key][1]) > 1e-5:
            return False, [
                '%s field mismatch with recorded solution (max) in the past.' %
                key
            ]

    # return success
    return True, []


if __name__ == '__main__':
    path = sys.argv[1]
    print('examples/AMR/shocktube_2d/check.py: checking output in directory ' +
          path)

    status_ok, info = verify_result(path=path)

    if status_ok:
        print('check.py: success!')
        sys.exit(0)  # normal exit

    for msg in info:
        print('check.py: failed! ' + msg)
    sys.exit(1)


def visualize_result(path, Lx, Ly):
    """ Needs cleanup. """
    import matplotlib.pyplot as p
    import arepo

    toinch = 0.393700787
    p.figure(figsize=np.array([14.7, 16]) * toinch, dpi=300)

    ax1 = p.subplot2grid((2, 2), (0, 0))
    ax2 = p.subplot2grid((2, 2), (0, 1))
    ax3 = p.subplot2grid((2, 2), (1, 0))
    ax4 = p.subplot2grid((2, 2), (1, 1))

    exact = np.loadtxt('e1rpex.out')  # NEEDED

    axes = [ax1, ax2, ax3, ax4]
    sims = [
        '../output/snap_010.hdf5', '../output_arepo/snap_010.hdf5',
        '../output_albada/snap_010.hdf5', '../output_minbee/snap_010.hdf5',
        '../output_superbee/snap_010.hdf5', '../output_vanleer/snap_010.hdf5'
    ]
    label = ['None', 'Arepo', 'Albada', 'Minbee', 'Superbee', 'Vanleer']

    sims = [
        '../output_superbee/snap_020.hdf5', '../output_vanleer/snap_020.hdf5',
        '../output_albada/snap_020.hdf5', '../output_minbee/snap_020.hdf5',
        '../output_arepo/snap_020.hdf5'
    ]
    label = ['Superbee', 'van Leer', 'van Albada', 'Minbee', 'Springel 2010']

    colors = [
        '#0072bd', '#d95319', '#edb120', '#7e2f8e', '#77ac30', '#4dbeee',
        '#a2142f'
    ]
    zorders = [30, 40, 50, 60, 70]

    lines = []

    l, = ax1.plot(exact[:, 0], exact[:, 1], ls='--', lw=0.5, c='k', zorder=10)
    lines.append(l)
    label.insert(0, 'analytical solution')

    ax2.plot(exact[:, 0], exact[:, 3], ls='--', lw=0.5, zorder=10)
    ax3.plot(exact[:, 0], exact[:, 2], ls='--', lw=0.5, zorder=10)
    ax4.plot(exact[:, 0], exact[:, 4], ls='--', lw=0.5, zorder=10)

    for j in np.arange(len(sims)):
        sn = arepo.Simulation(sims[j])
        sn.set_center([0.25, 1.2])
        g = sn.get_AMRline(sn.rho, gradient=sn.grar, axis=1, box=1.)
        l, = ax1.plot(g['x2'] - 0.7,
                      g['grid'],
                      lw=0.5,
                      color=colors[j],
                      zorder=zorders[j])
        lines.append(l)
        g = sn.get_AMRline(sn.pres, gradient=sn.grap, axis=1, box=1.)
        ax2.plot(g['x2'] - 0.7,
                 g['grid'],
                 lw=0.5,
                 color=colors[j],
                 zorder=zorders[j])
        g = sn.get_AMRline(sn.vel_y, gradient=sn.grav[:, 3:6], axis=1, box=1.)
        ax3.plot(g['x2'] - 0.7,
                 g['grid'],
                 lw=0.5,
                 color=colors[j],
                 zorder=zorders[j])
        g = sn.get_AMRline(sn.u, axis=1, box=1.)
        ax4.plot(g['x2'] - 0.7,
                 g['grid'],
                 lw=0.5,
                 color=colors[j],
                 zorder=zorders[j])

    ax1.set_ylim(0., 1.1)
    ax2.set_ylim(0., 1.1)
    ax3.set_ylim(-0.1, 1.5)
    ax4.set_ylim(1.9, 3.7)

    ax1.set_ylabel(r'$\rho$')
    ax2.set_ylabel(r'$v$')
    ax3.set_ylabel(r'$P$')
    ax4.set_ylabel(r'$u$')
    ax4.set_xlabel(r'$x$')
    ax3.set_xlabel(r'$x$')

    p.figlegend(lines,
                label,
                'lower center',
                ncol=3,
                handlelength=3,
                frameon=False)
    p.tight_layout()
    p.subplots_adjust(bottom=0.18)
    p.show()
    p.savefig('shocktube_compare.pdf', dpi=300)
