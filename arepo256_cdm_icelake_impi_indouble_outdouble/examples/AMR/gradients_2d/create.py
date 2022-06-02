#!/usr/bin/env python3
""" Idealized test, 2D test for density gradient computation with AMR. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 1.0
Ly = 1.0


def create_ics(path, filename='ics.hdf5'):
    # config
    gamma = 1.4
    pressure = 1.0

    np.random.seed(42)
    theta = np.deg2rad(57)

    # helper functions
    def analytic_rho(x, y):
        return np.cos(theta) * x + np.sin(theta) * y

    def make_half_grid(N, half):
        ndim = 2
        Ntot = int(N**ndim / 2)
        dx = Lx / float(N)
        x, y = np.meshgrid(0.5 * half + (np.arange(N / 2) + 0.5) * dx,
                           (np.arange(N) + 0.5) * dx,
                           indexing="ij")
        x = x.flatten()
        y = y.flatten()
        #z = np.zeros_like(x)
        vel = np.zeros([Ntot, 3])
        rho = analytic_rho(x, y)
        vol = (dx**ndim) * np.ones_like(x)
        mass = rho * vol
        ids = np.arange(Ntot, dtype="i")
        u = pressure / ((gamma - 1.0) * rho)

        r = {
            'ntot': Ntot,
            'ids': ids,
            'x': x,
            'y': y,
            'vel': vel,
            'rho': rho,
            'mass': mass,
            'u': u,
            'vol': vol
        }
        return r

    def combine(field):
        return np.concatenate([g1[field], g2[field]], axis=0)

    # create
    level = 6
    g1 = make_half_grid(1 << (level + 0), 0)
    g2 = make_half_grid(1 << (level + 1), 1)

    N = g1['ntot'] + g2['ntot']

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    u = np.zeros(N, dtype='float64')
    mass = np.zeros(N, dtype='float64')

    # assign gas cell properties
    pos[:, 0] = combine('x')
    pos[:, 1] = combine('y')

    vel[:] = combine('vel')
    mass[:] = combine('mass')
    u[:] = combine('u')
    # vol
    # rho

    ids = np.arange(N) + 1

    # write
    pt0 = {
        'Coordinates': pos,
        'Velocities': vel,
        'Masses': mass,
        'InternalEnergy': u,
        'ParticleIDs': ids
    }
    h = {'Flag_DoublePrecision': 1}

    utils.write_ic_file(os.path.join(path, filename), {'PartType0': pt0},
                        boxSize=Lx,
                        headerAttrs=h)


if __name__ == '__main__':
    path = sys.argv[1]
    print("examples/AMR/gradients_2d/create.py: creating ICs in directory " +
          path)
    create_ics(path)
