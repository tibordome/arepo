#!/usr/bin/env python3
""" Idealized test, AMR Kelvin-Helmholtz (non-static) in 2D. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('./examples/')  # test.sh
import utils

Lx = 1.0
Ly = 1.0

numTasksMPI = 4  # ngbtree error with 1


def create_ics(path, filename='ics.hdf5'):
    # config
    N_x = 32
    N_y = 32

    P = 2.5
    omega0 = 0.1
    sigma = 0.05 / np.sqrt(2)
    GAMMA = 5.0 / 3.0

    # derived
    dx = Lx / N_x
    dy = Ly / N_y

    N = int(N_x * N_y / 2 + N_x * N_y * 4 / 2)
    ind4 = int(N_x * N_y / 4)

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    u = np.zeros(N, dtype='float64')
    mass = np.zeros(N, dtype='float64')

    # set gas cell propeties
    mesh = np.meshgrid(np.arange(N_x), np.arange(N_y))

    posx = (mesh[0] * dx).reshape(N_x * N_y) + 0.5 * dx
    posy = (mesh[1] * dy).reshape(N_x * N_y) + 0.5 * dy

    mesh = np.meshgrid(np.arange(N_x * 2), np.arange(N_y * 2))
    posx2 = (mesh[0] * dx / 2.).reshape(N_x * N_y * 4) + 0.5 * dx / 2.0
    posy2 = (mesh[1] * dy / 2.).reshape(N_x * N_y * 4) + 0.5 * dy / 2.0

    pos[:ind4, 0] = posx[:ind4]
    pos[:ind4, 1] = posy[:ind4]
    mass[:ind4] = 1. * dx * dy
    u[:ind4] = P / ((GAMMA - 1) * 1.)
    vel[:ind4, 0] = -0.5

    pos[ind4:ind4 + N_x * N_y * 2, 0] = posx2[N_x * N_y:N_x * N_y * 3]
    pos[ind4:ind4 + N_x * N_y * 2, 1] = posy2[N_x * N_y:N_x * N_y * 3]
    mass[ind4:ind4 + N_x * N_y * 2] = 2. * dx * dy / 4.
    u[ind4:ind4 + N_x * N_y * 2] = P / ((GAMMA - 1) * 2.)
    vel[ind4:ind4 + N_x * N_y * 2, 0] = +0.5

    pos[-ind4:, 0] = posx[-ind4:]
    pos[-ind4:, 1] = posy[-ind4:]
    mass[-ind4:] = 1. * dx * dy
    u[-ind4:] = P / ((GAMMA - 1) * 1.)
    vel[-ind4:, 0] = -0.5

    # perturbation
    vel[:, 1] = omega0 * np.sin(4 * np.pi * pos[:, 0]) * (
        np.exp(-(pos[:, 1] - 0.25)**2 * 0.5 /
               (sigma**2)) + np.exp(-(pos[:, 1] - 0.75)**2 * 0.5 / (sigma**2)))

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
    print("examples/AMR/KH_2d/create.py: creating ICs in directory " + path)
    create_ics(path)
