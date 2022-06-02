#!/usr/bin/env python3
""" @package examples/collisionless_galaxy_3d/create.py
Code creates the output list; ICs need to be present already
collisionless_galaxy_3d uses GalIC initial condition code, Model M4

created by Rainer Weinberger, last modified 29.04.2018
"""
# load libraries
import sys  # system calls
import numpy as np  # scientific computing package
import h5py  # hdf5 format
import os  # operating system interface
import os.path
from subprocess import call  # execute bash commands

createICs = False
# input
simulation_directory = str(sys.argv[1])
print('examples/collisionless_galaxy_3d/create.py ' + simulation_directory)
# creating ICs
if createICs:
    status = call([
        'git', 'clone', 'https://github.com/djurin/GALIC.git',
        os.path.join(simulation_directory, 'galic')
    ])
    if status != 0:
        print('CREATE: ERROR: git clone failed!')
        sys.exit(status)
    cwd = os.getcwd()
    os.chdir(os.path.join(simulation_directory, 'galic'))
    ## copy Makefile and config file
    str = 'SYSTYPE="Darwin"'
    file = open('Makefile.systype', 'w')
    file.write(str + '\n')
    file.close()
    call([
        'cp',
        os.path.join(simulation_directory, 'Makefile.galic'), 'Makefile'
    ])
    str = 'DOUBLEPRECISION=1 \nHAVE_HDF5 \nVER_1_1 \n'
    file = open('Config.sh', 'w')
    file.write(str + '\n')
    file.close()
    status = call(['make'])
    if status != 0:
        print('CREATE: ERROR: make failed!')
        sys.exit(status)
    status = call([
        'mpiexec', '-np', '1',
        os.path.join('.', 'galic'),
        os.path.join(cwd, 'examples', 'isolated_galaxy_collisionless_3d',
                     'param_galic.txt')
    ])
    if status != 0:
        print('CREATE: ERROR: execution failed!')
        sys.exit(status)
    os.chdir(cwd)
else:
    call([
        'cp',
        os.path.join(simulation_directory, 'ICs', 'ICs'),
        os.path.join(simulation_directory, 'snap_010')
    ])
# set output times
outputTimes = np.linspace(0.0, 1.0, 10, dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=int)
# write output list file
data = np.array([outputTimes, ones]).T
np.savetxt(os.path.join(simulation_directory, 'output_list.txt'),
           data,
           fmt='%g %1.f')
