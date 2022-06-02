#!/usr/bin/env python3
""" @package ./examples/OrszagTang/create_common.py
Code that creates 2d Orszag Tang  initial conditions
parameters are identical to Pakmor & Springel (2013), MNRAS 432, 176

created by Ruediger Pakmor, last modified 10.12.2021 -- comments welcome
"""

#### load libraries
import sys  # system specific calls
import os.path
import numpy as np  ## load numpy
import h5py  ## load h5py; needed to write initial conditions in hdf5 format

# make sure that this file's directory is in sys.path
# (see also https://bugs.python.org/issue17639)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import create_parameters

if hasattr(create_parameters, 'numTasksMPI'):
    numTasksMPI = create_parameters.numTasksMPI


FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

# parameters
Boxsize = FloatType(1.0)
CellsPerDimension = IntType(256)
gamma  = 5. / 3
fourpi = 4. * np.pi


def create_ics(path, filename='ics.hdf5'):
    NumberOfCells = CellsPerDimension**2

    # set up grid: cartesian 2d grid
    # spacing
    dx = Boxsize / FloatType(CellsPerDimension)
    # position of first and last cell
    pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx

    # set up grid
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy = np.meshgrid(Grid1d, Grid1d)
    Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
    Pos[:, 0] = xx.reshape(NumberOfCells)
    Pos[:, 1] = yy.reshape(NumberOfCells)
    # calculate distance from center
    xPosFromCenter = (Pos[:, 0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:, 1] - 0.5 * Boxsize)
    Radius = np.sqrt(xPosFromCenter**2 + yPosFromCenter**2)

    Density = np.ones(NumberOfCells, dtype=FloatType) * gamma**2 / fourpi

    Vel       =   np.zeros([NumberOfCells, 3], dtype=FloatType)
    Vel[:, 0] = - np.sin( 2. * np.pi * Pos[:, 1] )
    Vel[:, 1] = + np.sin( 2. * np.pi * Pos[:, 0] )

    Vel[:, 0] += create_parameters.vx_boost

    Bfld       =   np.zeros([NumberOfCells, 3], dtype=FloatType)
    Bfld[:, 0] = - np.sin( 2. * np.pi * Pos[:, 1] )
    Bfld[:, 1] = + np.sin( 4. * np.pi * Pos[:, 0] )

    Uthermal = 1.5 / Density * gamma / fourpi

    # write *.hdf5 file; minimum number of fields required by Arepo
    with h5py.File(os.path.join(path, filename), 'w') as IC:
        # create HDF5 groups
        header = IC.create_group('Header')
        part0 = IC.create_group('PartType0')

        # header entries
        NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
        header.attrs.create('NumPart_ThisFile', NumPart)
        header.attrs.create('NumPart_Total', NumPart)
        header.attrs.create('NumPart_Total_HighWord', np.zeros(6, dtype=IntType))
        header.attrs.create('MassTable', np.zeros(6, dtype=IntType))
        header.attrs.create('Time', 0.0)
        header.attrs.create('Redshift', 0.0)
        header.attrs.create('BoxSize', Boxsize)
        header.attrs.create('NumFilesPerSnapshot', 1)
        header.attrs.create('Omega0', 0.0)
        header.attrs.create('OmegaB', 0.0)
        header.attrs.create('OmegaLambda', 0.0)
        header.attrs.create('HubbleParam', 1.0)
        header.attrs.create('Flag_Sfr', 0)
        header.attrs.create('Flag_Cooling', 0)
        header.attrs.create('Flag_StellarAge', 0)
        header.attrs.create('Flag_Metals', 0)
        header.attrs.create('Flag_Feedback', 0)
        if Pos.dtype == np.float64:
            header.attrs.create('Flag_DoublePrecision', 1)
        else:
            header.attrs.create('Flag_DoublePrecision', 0)

        # copy datasets
        part0.create_dataset('ParticleIDs', data=np.arange(1, NumberOfCells + 1))
        part0.create_dataset('Coordinates', data=Pos)
        part0.create_dataset('Masses', data=Density)
        part0.create_dataset('Velocities', data=Vel)
        part0.create_dataset('InternalEnergy', data=Uthermal)
        part0.create_dataset('MagneticField', data=Bfld)


if __name__ == '__main__':
    simulation_directory = str(sys.argv[1])
    if len(sys.argv) > 3:
        CellsPerDimension = IntType(sys.argv[3])
    print(__file__ + ': creating ICs in directory ' + simulation_directory)
    create_ics(simulation_directory, filename='ics.hdf5')
