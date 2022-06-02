#!/usr/bin/env python3
""" @package examples/galaxy_merger_star_formation_3d/create.py
Code creates the output list; ICs need to be present already
galaxy_merger_star_formation_3d uses MakeNewDisk and CombineGalaxies codes
to create the initial conditions

created by Rainer Weinberger, last modified 09.03.2019
"""
# load libraries
import sys  # system calls
import os  # operating system calls
import os.path
import numpy as np  # scientific computing package
import h5py  # hdf5 format
from subprocess import call  # execute shell commands
# input
simulation_directory = str(sys.argv[1])
print('create.py ' + simulation_directory)
# set output times
outputTimes = np.linspace(0.0, 3.0, 32, dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=int)
# write output list file
data = np.array([outputTimes, ones]).T
np.savetxt(os.path.join(simulation_directory, 'output_list.txt'),
           data,
           fmt='%g %1.f')

if len(sys.argv) > 2:
    arepopath = sys.argv[2]
else:
    # TODO: find a better way to determine AREPO path
    arepopath = os.path.join(simulation_directory, '..', '..', '..')

simdir_abs = os.path.abspath(simulation_directory)
# create backgroundgrid ICs from SPH ICs
## compile Arepo with ADDBACKGROUNDGRID
os.chdir(arepopath)
res = call([
    'make',
    'CONFIG=' + os.path.join(simdir_abs, 'Config_ADDBACKGROUNDGRID.sh'),
    'BUILD_DIR=' + os.path.join(simdir_abs, 'build_ADDBACKGROUNDGRID'),
    'EXEC=' + os.path.join(simdir_abs, 'Arepo_ADDBACKGROUNDGRID')
])
if res != 0:
    sys.exit(int(res))

## execute Arepo with ADDBACKGROUNDGRID from run directory
os.chdir(simdir_abs)
res = call([
    'mpiexec', '-n', '1', './Arepo_ADDBACKGROUNDGRID',
    'param_ADDBACKGROUNDGRID.txt'
])
sys.exit(int(res))
