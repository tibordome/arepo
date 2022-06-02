#!/usr/bin/env python3
""" Idealized test, AMR shocktube in 2D. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 0.5
Ly = 2.0


def create_ics(path, filename='ics.hdf5'):
    # config
    N_x = 128 / 4
    N_y = 128
    GAMMA = 1.4

    # derived
    dx = Lx / N_x
    dy = Ly / N_y
    N = int(N_x * N_y)

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    rho = np.zeros(N, dtype='float64')
    u = np.zeros(N, dtype='float64')

    # assign gas cell properties
    mesh = np.meshgrid(np.arange(N_x), np.arange(N_y))

    pos[:, 0] = (mesh[0] * dx).reshape(N) + 0.5 * dx
    pos[:, 1] = (mesh[1] * dy).reshape(N) + 0.5 * dy

    ind = int(N_y / 2 * N_x)

    rho[:ind] = 1.0
    rho[ind:] = 0.125

    mass = rho * dx * dy

    u[:ind] = 1.0 / ((GAMMA - 1) * 1.0)
    u[ind:] = 0.1 / ((GAMMA - 1) * 0.125)

    vel[:ind, 1] = 0.75
    vel[ind:, 1] = 0.0

    id = np.arange(N) + 1

    # write
    pt0 = {
        'Coordinates': pos,
        'Velocities': vel,
        'Masses': mass,
        'InternalEnergy': u,
        'ParticleIDs': id
    }
    h = {'Flag_DoublePrecision': 1}

    utils.write_ic_file(os.path.join(path, filename), {'PartType0': pt0},
                        boxSize=Ly,
                        headerAttrs=h)


if __name__ == '__main__':
    path = sys.argv[1]
    print("examples/AMR/shocktube_2d/create.py: creating ICs in directory " +
          path)
    create_ics(path)
