#!/usr/bin/env python3
""" Idealized test, check force accuracy. Common ICs for all subtests. """
import os.path
import sys
import numpy as np

sys.path.append('./')  # test.py
sys.path.append('examples/')  # test.sh
import utils

Lx = Ly = Lz = 100.0
N = 100

# Result seems to depend on compiler/MPI library and number of MPI tasks;
# seems to work with ICC/Intel MPI and 1, 2, 4, or 8 tasks (?)
numTasksMPI = 8


def create_ics(path, filename='ics.hdf5'):
    # random coordinates in [49.5, 50.5]
    np.random.seed(424242)
    pos = np.random.uniform(low=49.5, high=50.5, size=(N, 3)).astype('float32')

    # zero velocities
    vel = np.zeros((N, 3), dtype='float32')

    # mass of first particle 1.0, otherwise zero
    mass = np.zeros(N, dtype='float32')
    mass[0] = 1.0

    # sequential IDs, starting at 1
    ids = np.arange(N, dtype='int32') + 1

    # write
    pt2 = {
        'Coordinates': pos,
        'Velocities': vel,
        'Masses': mass,
        'ParticleIDs': ids
    }

    utils.write_ic_file(os.path.join(path, filename), {'PartType2': pt2},
                        boxSize=Lx)


if __name__ == '__main__':
    path = sys.argv[1]
    print('examples/ForceLawTests/create.py: creating ICs in directory ' +
          path)

    try:
        create_ics(path, filename='ics.hdf5')
    except:
        sys.exit(1)
