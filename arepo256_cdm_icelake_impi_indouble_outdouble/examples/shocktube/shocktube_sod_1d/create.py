#!/usr/bin/env python3
""" Idealized test, Shocktube (Sod) in true 1D, e.g. Springel+ (2010) Section 8.2. """
import os.path
import sys
import numpy as np

sys.path.insert(0, './')  # test.py
sys.path.insert(0, 'examples/')  # test.sh
import utils

# 1D code can only be used serially
numTasksMPI = 1
Lx = 20.0
Ly = 2.0
numPart = 100


def create_ics(path, filename='ics.hdf5'):
    # config
    gamma = 1.4  # adiabatic index
    P1 = 1.0  # left pressure
    P2 = 0.1795  # right pressure
    rho1 = 1.0  # left density
    rho2 = 0.25  # right density

    # derived properties
    Nx = numPart
    dx = Lx / Nx

    # allocate
    pos = np.zeros((Nx, 3), dtype='float32')
    vel = np.zeros((Nx, 3), dtype='float32')
    mass = np.zeros(Nx, dtype='float32')
    u = np.zeros(Nx, dtype='float32')
    id = np.arange(Nx, dtype='int32') + 1

    # assign gas cell properties
    for i in range(Nx):
        pos[i, 0] = i * dx + dx / 2.0
        pos[i, 1] = 0.0
        pos[i, 2] = 0.0

        # left
        if pos[i, 0] < Lx / 2.0:
            mass[i] = rho1 * dx
            u[i] = P1 / rho1 / (gamma - 1.0)

        # right
        if pos[i, 0] >= Lx / 2.0:
            mass[i] = rho2 * dx
            u[i] = P2 / rho2 / (gamma - 1.0)

    # write
    pt0 = {
        'Coordinates': pos,
        'Velocities': vel,
        'Masses': mass,
        'InternalEnergy': u,
        'ParticleIDs': id
    }

    utils.write_ic_file(os.path.join(path, filename), {'PartType0': pt0},
                        boxSize=Lx)


if __name__ == '__main__':
    path = sys.argv[1]
    print(
        'examples/shocktube/shocktube_sod_1d/create.py: creating ICs in directory '
        + path)

    try:
        create_ics(path, filename='ics.hdf5')
    except:
        sys.exit(1)
