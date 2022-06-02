#!/usr/bin/env python3
""" Idealized test, AMR sedov in 3D. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 1.0
Ly = 1.0
Lz = 1.0

numTasksMPI = 4


def create_ics(path, filename='ics.hdf5'):
    # config
    N_x = 16
    N_y = 16
    N_z = 16

    P = 1e-4
    GAMMA = 5.0 / 3.0

    # derived
    dx = Lx / N_x
    dy = Ly / N_y
    dz = Lz / N_z
    N = N_x * N_y * N_z

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    u = np.zeros(N, dtype='float64')
    mass = np.zeros(N, dtype='float64')

    # assign gas cell properties
    mesh = np.mgrid[0:N_x, 0:N_y, 0:N_z]

    pos[:, 0] = (mesh[0] * dx).reshape(N_x * N_y * N_z) + 0.5 * dx
    pos[:, 1] = (mesh[1] * dy).reshape(N_x * N_y * N_z) + 0.5 * dy
    pos[:, 2] = (mesh[2] * dz).reshape(N_x * N_y * N_z) + 0.5 * dz

    mass[:] = 1.0 * dx * dy * dz
    u[:] = P / ((GAMMA - 1) * 1.0)

    E0 = 1.0 / mass[0] / 8.0

    u[int(N_x / 2 - 1) * N_y * N_z + int(N_y / 2 - 1) * N_z +
      int(N_z / 2 - 1)] = E0
    u[int(N_x / 2 - 1) * N_y * N_z + int(N_y / 2 - 1) * N_z +
      int(N_z / 2)] = E0
    u[int(N_x / 2 - 1) * N_y * N_z + int(N_y / 2) * N_z +
      int(N_z / 2 - 1)] = E0
    u[int(N_x / 2 - 1) * N_y * N_z + int(N_y / 2) * N_z + int(N_z / 2)] = E0
    u[int(N_x / 2) * N_y * N_z + int(N_y / 2 - 1) * N_z +
      int(N_z / 2 - 1)] = E0
    u[int(N_x / 2) * N_y * N_z + int(N_y / 2 - 1) * N_z + int(N_z / 2)] = E0
    u[int(N_x / 2) * N_y * N_z + int(N_y / 2) * N_z + int(N_z / 2 - 1)] = E0
    u[int(N_x / 2) * N_y * N_z + int(N_y / 2) * N_z + int(N_z / 2)] = E0

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
    print("examples/AMR/sedov_3d/create.py: creating ICs in directory " + path)
    create_ics(path)
