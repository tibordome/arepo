#!/usr/bin/env python3
""" Idealized test, AMR Rayleigh-Taylor in 2D. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 0.5
Ly = 2.0

numTasksMPI = 4  # ngbtree error with 1, amr_ngb error with >=2


def create_ics(path, filename='ics.hdf5'):
    # config
    N_x = 64
    N_y = 4 * N_x

    P = 2.5
    g = -0.1
    rho1 = 1.0
    rho2 = 2.0
    GAMMA = 1.4

    # derived
    dx = Lx / N_x
    dy = Ly / N_y
    N = N_x * N_y

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    u = np.zeros(N, dtype='float64')
    rho = np.zeros(N, dtype='float64')

    # assign gas cell properties
    mesh = np.meshgrid(np.arange(N_x), np.arange(N_y))

    pos[:, 0] = (mesh[0] * dx).reshape(N) + 0.5 * dx
    pos[:, 1] = (mesh[1] * dy).reshape(N) + 0.5 * dy

    vel[:, 1] = 0.0025 * (1.0 - np.cos(4.0 * np.pi * pos[:, 0])) * (
        1.0 - np.cos(np.pi * pos[:, 1]))

    rho[:int(N_y / 2 * N_x)] = rho1
    rho[int(N_y / 2 * N_x):] = rho2
    mass = rho * dx * dy

    u[:] = (P + g * (pos[:, 1] - 1.0) * rho) / ((GAMMA - 1) * rho)

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
    create_ics(sys.argv[1])
