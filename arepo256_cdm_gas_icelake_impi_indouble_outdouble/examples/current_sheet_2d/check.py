#!/usr/bin/env python3
""" @package ./examples/current_sheet_2d/check.py
Code that checks results of 2d current sheet problem

created by Rainer Weinberger, last modified 13.03.2019
"""
# load libraries
import sys  ## system calls
import numpy as np  ## load numpy
import h5py  ## load h5py; needed to read snapshots
import os  # file specific calls
import os.path
import matplotlib
import matplotlib.pyplot as plt  # needs to be active for plotting!
matplotlib.rc_file_defaults()

makeplots = True
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print('current_sheet_2d: checking simulation output in directory ' +
      simulation_directory)

FloatType = np.float64  # double precision: np.float64, for single use np.float32

## open initial conditiions to get parameters
try:
    data = h5py.File(os.path.join(simulation_directory, 'IC.hdf5'), 'r')
except:
    print('Could not open initial conditions!')
    sys.exit(-1)
Boxsize = FloatType(data['Header'].attrs['BoxSize'])
NumberOfCells = np.int32(data['Header'].attrs['NumPart_Total'][0])
CellsPerDimension = np.sqrt(NumberOfCells)  ## 2d sim

## parameters for initial state
density_0 = 1.0
velocity_radial_0 = -1.0  ## radial inflow velocity
pressure_0 = 1.0e-4
gamma = 5. / 3.  ## note: this has to be consistent with the parameter settings for Arepo!
utherm_0 = pressure_0 / (gamma - 1.0) / density_0

## maximum L1 error after one propagation; based on tests
DeltaMaxAllowed = 0.1 * (FloatType(CellsPerDimension) / 150.0)**-1

Time = []
MagneticEnergy = []
# loop over all output files
i_file = 0
while True:
    # try to read in snapshot
    directory = os.path.join(simulation_directory, 'output')
    filename = 'snap_%03d.hdf5' % i_file
    try:
        data = h5py.File(os.path.join(directory, filename), 'r')
    except (OSError, IOError):
        # should have at least 9 snapshots
        if i_file <= 8:
            print('Could not find snapshot ' + filename + '!')
            sys.exit(1)
        break
    # get simulation data

    B = np.array(data['PartType0']['MagneticField'], dtype=FloatType)

    B2 = B[:, 0]**2 + B[:, 1]**2 + B[:, 2]**2

    Time.append(FloatType(data['Header'].attrs['Time']))
    MagneticEnergy.append(np.sum(B2) / 8.0 / np.pi)

    i_file += 1

Time = np.array(Time)
MagneticEnergy = np.array(MagneticEnergy)

if makeplots:
    if not os.path.exists(os.path.join(simulation_directory, 'plots')):
        os.mkdir(os.path.join(simulation_directory, 'plots'))

    fig, ax = plt.subplots(1)
    ax.plot(Time, MagneticEnergy / MagneticEnergy[0])
    fig.savefig(
        os.path.join(simulation_directory, 'plots', 'MagneticEnergy.pdf'))

    for i_file in [0, 2, 4, 8]:
        filename = 'snap_%03d.hdf5' % i_file
        try:
            data = h5py.File(os.path.join(directory, filename), 'r')
        except (OSError, IOError):
            break

        VoronoiPos = np.array(data['PartType0']['Coordinates'],
                              dtype=np.float64)
        MagneticField = np.array(data['PartType0']['MagneticField'],
                                 dtype=np.float64)
        Time = data['Header'].attrs['Time']
        Boxsize = data['Header'].attrs['BoxSize']
        NumberOfCells = np.int32(data['Header'].attrs['NumPart_Total'][0])

        Nplot = 256
        from scipy import spatial  # needed for KDTree that we use for nearest neighbour search
        Edges1d = np.linspace(0.,
                              Boxsize,
                              Nplot + 1,
                              endpoint=True,
                              dtype=np.float64)
        Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
        xx, yy = np.meshgrid(Grid1d, Grid1d)
        Grid2D = np.array([xx.reshape(Nplot**2), yy.reshape(Nplot**2)]).T
        dist, cells = spatial.KDTree(VoronoiPos[:, :2]).query(Grid2D, k=1)

        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = plt.axes([0, 0, 1, 1])
        ax.streamplot(Grid1d,
                      Grid1d,
                      MagneticField[cells, 0].reshape((Nplot, Nplot)),
                      MagneticField[cells, 1].reshape((Nplot, Nplot)),
                      density=2.,
                      linewidth=0.5,
                      arrowsize=1e-6,
                      integration_direction='both',
                      minlength=0.3)
        ax.set_xlim(0, Boxsize)
        ax.set_ylim(0, Boxsize)
        plt.text(0.02,
                 0.95,
                 r'$N=%d^2,\ t=%3.1f$' %
                 (np.int32(np.sqrt(NumberOfCells)), Time),
                 color='k',
                 transform=fig.transFigure,
                 bbox={
                     'boxstyle': 'round',
                     'facecolor': 'w'
                 })

        fig.savefig(os.path.join(simulation_directory, 'plots',
                                 'streamlines_%03d.pdf' % i_file),
                    dpi=300)
        plt.close(fig)

## if everything is ok
if not (np.min(MagneticEnergy / MagneticEnergy[0]) > 0.92):
    sys.exit(1)
