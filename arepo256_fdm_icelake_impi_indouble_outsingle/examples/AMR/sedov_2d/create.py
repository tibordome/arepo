#!/usr/bin/env python3
""" Idealized test, AMR sedov in 2D. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 1.0
Ly = 1.0

numTasksMPI = 4


def create_ics(path, filename='ics.hdf5'):
    # config
    N_x = 64
    N_y = 64

    P = 1e-4
    GAMMA = 5.0 / 3.0

    # derived
    dx = Lx / N_x
    dy = Ly / N_y
    N = N_x * N_y

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    u = np.zeros(N, dtype='float64')
    mass = np.zeros(N, dtype='float64')

    # assign gas cell properties
    mesh = np.meshgrid(np.arange(N_x), np.arange(N_y))

    pos[:, 0] = (mesh[0] * dx).reshape(N_x * N_y) + 0.5 * dx
    pos[:, 1] = (mesh[1] * dy).reshape(N_x * N_y) + 0.5 * dy

    mass[:] = 1.0 * dx * dy
    u[:] = P / ((GAMMA - 1) * 1.0)

    E0 = 1.0 / mass[0] / 4.0

    u[int((N_x / 2 - 1) * N_y + (N_y / 2 - 1))] = E0
    u[int((N_x / 2 - 1) * N_y + (N_y / 2))] = E0
    u[int((N_x / 2) * N_y + (N_y / 2 - 1))] = E0
    u[int((N_x / 2) * N_y + (N_y / 2))] = E0

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
    print("examples/AMR/sedov_2d/create.py: creating ICs in directory " + path)
    create_ics(path)
