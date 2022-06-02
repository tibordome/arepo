#!/usr/bin/env python3
"""
Code that creates 3d cosmological box ICs;
some ICs provided, but can be created also at startup if
ic_creation is set to 'music' or 'ngenic'.

For ic_creation == 'copy':
Code creates the output list; ICs need to be present already.

created by Rainer Weinberger, last modified 28.02.2019
"""
# load libraries
import sys  # system calls
import numpy as np  # scientific computing package
import h5py  # hdf5 format
import os  # operating system interface
import os.path
from subprocess import call  # execute bash commands

# make sure that this file's directory is in sys.path
# (see also https://bugs.python.org/issue17639)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import create_parameters

## create new ics with 'music' or 'ngenic' or just 'copy' existing ones to the run directory
ic_creation = 'copy'
# input
simulation_directory = str(sys.argv[1])
print('create.py ' + simulation_directory)
# initial conditions: either copy or create with code
if ic_creation == 'copy':
    ## copy and use provided initial conditions
    call([
        'cp',
        os.path.join(simulation_directory, create_parameters.ic_directory),
        os.path.join(simulation_directory, 'ics')
    ])

elif ic_creation == 'music':
    ## create new initial conditions with the MUSIC code
    status = call([
        'git', 'clone', 'https://bitbucket.org/ohahn/music',
        os.path.join(simulation_directory, 'music')
    ])
    if status != 0:
        print('CREATE: ERROR: hg clone failed!')
        sys.exit(status)
    cwd = os.getcwd()
    os.chdir(os.path.join(simulation_directory, 'music'))
    status = call(['make'])
    if status != 0:
        print('CREATE: ERROR: make failed!')
        sys.exit(status)
    status = call(['./MUSIC', os.path.join(cwd, 'param_music.txt')])
    if status != 0:
        print('CREATE: ERROR: execution failed!')
        sys.exit(status)
    os.chdir(cwd)

elif ic_creation == 'ngenic':
    ## create new initial conditions with the N-GenIC code
    status = call([
        'git', 'clone',
        'https://gitlab.mpcdf.mpg.de/ext-c2c74fbfcdff/ngenic.git',
        os.path.join(simulation_directory, 'ngenic')
    ])
    if status != 0:
        print('CREATE: ERROR: git clone failed!')
        sys.exit(status)
    cwd = os.getcwd()
    os.chdir(os.path.join(simulation_directory, 'ngenic'))
    status = call(['make'])
    if status != 0:
        print('CREATE: ERROR: make failed!')
        sys.exit(status)
    status = call([
        'mpiexec', '-np', '1', './N-GenIC',
        os.path.join(cwd, 'param_ngenic.txt')
    ])
    if status != 0:
        print('CREATE: ERROR: execution failed!')
        sys.exit(status)
    os.chdir(cwd)
else:
    print(
        'CREATE: ERROR: no valid option for ic creation! choose "copy", "music" or "ngenic"'
    )
    sys.exit(1)
# set output times
outputTimes = np.array(create_parameters.outputTimes, dtype=np.float64)
ones = np.ones(outputTimes.shape, dtype=int)

# copy treecool file to run directory
if len(sys.argv) > 2:
    arepopath = sys.argv[2]
else:
    arepopath = '.'
call([
    'cp',
    os.path.join(arepopath, 'data', 'TREECOOL_ep'),
    os.path.join(simulation_directory, 'TREECOOL_ep')
])

# write output list file
data = np.array([outputTimes, ones]).T
np.savetxt(os.path.join(simulation_directory, 'output_list.txt'),
           data,
           fmt='%g %1.f')
