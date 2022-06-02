#!/usr/bin/env python3
""" @package ./examples/Yee_Cartesian_2d/create.py
Code that creates 2d Yee vortex on initially Cartesian mesh
parameters are identical to Pakmor et al. (2016), MNRAS 455, 1134
"""

#### load libraries
import sys  # system specific calls
import numpy as np  ## load numpy
import h5py  ## load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print('examples/Yee_Cartesian_2d/create.py: creating ICs in directory ' +
      simulation_directory)
# initial condition parameters
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(10.0)
if len(sys.argv) > 3:
    CellsPerDimension = IntType(sys.argv[3])
else:
    CellsPerDimension = IntType(50)
NumberOfCells = CellsPerDimension**2

#### parameters
Tinf = FloatType(1.0)
beta = FloatType(5.0)
gamma = FloatType(
    1.4
)  ## note: this has to be consistent with the parameter settings for Arepo!


# set up grid: cartesian 2d grid
## spacing
dx = Boxsize / FloatType(CellsPerDimension)
## position of first and last cell
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx

## set up grid
Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
xx, yy = np.meshgrid(Grid1d, Grid1d)
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:, 0] = xx.reshape(NumberOfCells)
Pos[:, 1] = yy.reshape(NumberOfCells)
## calculate distance from center
xPosFromCenter = (Pos[:, 0] - 0.5 * Boxsize)
yPosFromCenter = (Pos[:, 1] - 0.5 * Boxsize)
Radius = np.sqrt(xPosFromCenter**2 + yPosFromCenter**2)


# set up hydrodynamical quantitites
Temperature = Tinf - (
    (gamma - 1.0) * beta * beta /
    (8.0 * np.pi * np.pi * gamma) * np.exp(1.0 - Radius * Radius))
exponent = 1.0 / (gamma - 1.0)

Density = Temperature**exponent
Vel = np.zeros([NumberOfCells, 3], dtype=FloatType)
Vel[:, 0] = -0.5 * yPosFromCenter * beta / np.pi * np.exp(
    0.5 * (1.0 - Radius * Radius))
Vel[:, 1] = 0.5 * xPosFromCenter * beta / np.pi * np.exp(
    0.5 * (1.0 - Radius * Radius))
Uthermal = Temperature / (gamma - 1.0)
# write *.hdf5 file; minimum number of fields required by Arepo
IC = h5py.File(simulation_directory + '/IC.hdf5', 'w')

## create hdf5 groups
header = IC.create_group('Header')
part0 = IC.create_group('PartType0')

## header entries
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

## copy datasets
part0.create_dataset('ParticleIDs', data=np.arange(1, NumberOfCells + 1))
part0.create_dataset('Coordinates', data=Pos)
part0.create_dataset('Masses', data=Density)
part0.create_dataset('Velocities', data=Vel)
part0.create_dataset('InternalEnergy', data=Uthermal)

## close file
IC.close()
