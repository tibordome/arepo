#!/usr/bin/env python3
""" @package ./examples/khi_mhd_2d/create.py
Create initial conditions for 2D magnetohydrodynamic
Kelvin-Helmholtz instability.

The test is described in the paper:

Berlok & Pfrommer (2019), On the Kelvin-Helmholtz instability with smooth
initial conditions - linear theory and simulations, MNRAS, 485, 908-923
https://doi.org/10.1093/mnras/stz379

On top of this we add a magnetic field in the x-direction which has a strength
large enough to inhibit the growth rate of the KHI. Here we have beta=10, the
paper has beta=5.

The KHI is excited by initializing the fastest growing linear eigenmode
for this equilibrium. The linear solution is loaded from a .p file.
More information on the equilibrium and the method for initializing the
problem will be provided in the near future.

One note on irregular vs Cartesian mesh:
The KHI needs much higher resolution and/or a higher initial amplitude of the
perturbation when using an irregular mesh than when using a Cartesian mesh.
This is because the background velocity profile does not have div v = 0 to
machine precision on a irregular mesh.

created by Thomas Berlok, last modified 19.09.2018 -- comments welcome
"""

# load libraries
import os.path
import sys  # load sys; needed for exit codes
import numpy as np  # load numpy
import h5py  # load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print('examples/khi_mhd_2d/create.py: creating ICs in directory ' +
      simulation_directory)
# initial condition parameters
FilePath = os.path.join(simulation_directory, 'IC.hdf5')

# Load linear solution
dat1 = np.loadtxt(
    os.path.join(simulation_directory, 'khi-mhd-uniform-real.txt'))
dat2 = np.loadtxt(
    os.path.join(simulation_directory, 'khi-mhd-uniform-imag.txt'))
f = {
    'drho': dat1[:, 0] + 1j * dat2[:, 0],
    'dvx': dat1[:, 1] + 1j * dat2[:, 1],
    'dvy': dat1[:, 2] + 1j * dat2[:, 2],
    'dT': dat1[:, 3] + 1j * dat2[:, 3],
    'dbx': dat1[:, 4] + 1j * dat2[:, 4],
    'dby': dat1[:, 5] + 1j * dat2[:, 5],
    'dA': dat1[:, 6] + 1j * dat2[:, 6]
}

kx = 5.4234088072021045

Lx = 2 * np.pi / kx
Ly = 2.0

gamma = 5. / 3.
uflow = 1
a = 0.05
p0 = 1.0
rho0 = 1.0
beta = 10.0
B0 = np.sqrt(2 * p0 / beta) * np.sqrt(4 * np.pi)
y1 = 0.5
y2 = 1.5

# Grid configuration
hexagonal = False
Ny = 128
Nx = int(np.ceil(Lx / Ly * Ny))
dx = Lx / Nx
dy = Ly / Ny

# Amplitude of perturbation
ampl = 1e-2


def to_real(x, compl):
    return 2 * (np.cos(kx * x) * compl.real - np.sin(kx * x) * compl.imag)


def interpolate(z, f):
    """Interpolate component f onto the grid z"""
    import numpy as np
    N = len(f)
    Lz = 2
    dz = Lz / N

    ak = 2 * np.fft.rfft(f) / N
    n = np.arange(N // 2 + 1)

    def to_grid(z):
        cos = np.sum(ak[1:].real * np.cos(2.0 * np.pi * n[1:] *
                                          (z - dz / 2) / Lz))
        sin = -np.sum(ak[1:].imag * np.sin(2.0 * np.pi * n[1:] *
                                           (z - dz / 2) / Lz))
        y = ak[0].real / 2.0 + cos + sin
        return y

    to_grid_v = np.vectorize(to_grid)

    return to_grid_v(z)


def perturb(x, y, key):
    """Interpolate linear solution onto x and y"""
    if key == 'dp':
        dr = interpolate(y, (f['drho'] + f['dT']).real)
        di = interpolate(y, (f['drho'] + f['dT']).imag)
    else:
        dr = interpolate(y, f[key].real)
        di = interpolate(y, f[key].imag)
    d = dr + 1j * di
    return ampl * to_real(x, d)


def rho_an(x, y):
    return rho0 + rho0 * perturb(x, y, 'drho')


def vx_an(x, y):
    vx = uflow*(np.tanh((y-y1)/a) - np.tanh((y-y2)/a) - 1.0) \
        + perturb(x, y, 'dvx')
    return vx


def vy_an(x, y):
    return perturb(x, y, 'dvy')


def bx_an(x, y):
    return B0 * (1 + perturb(x, y, 'dbx'))


def by_an(x, y):
    return B0 * perturb(x, y, 'dby')


def p_an(x, y):
    return p0 + p0 * perturb(x, y, 'dp')


def Az_an(x, y):
    return -B0 * perturb(x, y, 'dA')


def make_ics(filename='ics.hdf5', writefunc=None):
    # allocate
    pos = np.zeros((Nx * Ny, 3), dtype='float64')
    vel = np.zeros((Nx * Ny, 3), dtype='float64')
    bfld = np.zeros((Nx * Ny, 3), dtype='float64')
    A = np.zeros((Nx * Ny, 3), dtype='float64')
    dens = np.zeros(Nx * Ny, dtype='float64')
    u = np.zeros(Nx * Ny, dtype='float64')
    id = np.arange(Nx * Ny, dtype='int64') + 1

    x = np.arange(Nx) * dx + dx / 2
    y = np.arange(Ny) * dy + dy / 2
    xx, yy = np.meshgrid(x, y)
    if hexagonal:
        for j in range(Ny):
            xx[j, :] += (j % 2) * 0.45 * dx

    # IPython.embed()
    x = xx.flatten()
    y = yy.flatten()
    index = range(Nx * Ny)

    rho = rho_an(x, y)
    vx = vx_an(x, y)
    vy = vy_an(x, y)
    P = p_an(x, y)

    bx = bx_an(x, y)
    by = by_an(x, y)

    Az = Az_an(x, y)

    pos[index, 0] = x
    pos[index, 1] = y
    pos[index, 2] = 0.0

    dens[index] = rho
    u[index] = P / rho / (gamma - 1.0)
    vel[index, 0] = vx
    vel[index, 1] = vy

    bfld[index, 0] = bx
    bfld[index, 1] = by

    A[index, 2] = Az

    pt0 = {
        'Coordinates': pos,
        'Velocities': vel,
        'Masses': dens,
        'MagneticField': bfld,
        'MagneticVectorPotential': A / np.sqrt(4 * np.pi),
        'InternalEnergy': u,
        'ParticleIDs': id
    }

    writefunc(filename, {'PartType0': pt0}, boxSize=Ly)
    return True


def write_ic_file(fileName, partTypes, boxSize, massTable=None):
    """ Helper to write a HDF5 IC file. partTypes is a dictionary with keys of
    the form PartTypeX, each of which is its own dictionary of particle fields
    and ndarrays. boxSize is a scalar float, and massTable a 6-element float
    array, if specified. """

    nPartTypes = 6

    with h5py.File(fileName, 'w') as f:
        # write each PartTypeX group and datasets
        for ptName in partTypes.keys():
            g = f.create_group(ptName)
            for field in partTypes[ptName]:
                g[field] = partTypes[ptName][field]

        # set particle counts
        NumPart = np.zeros(nPartTypes, dtype='int64')
        for ptName in partTypes.keys():
            ptNum = int(ptName[-1])
            NumPart[ptNum] = partTypes[ptName]['ParticleIDs'].size

        # create standard header
        h = f.create_group('Header')
        h.attrs['BoxSize'] = boxSize
        h.attrs['NumFilesPerSnapshot'] = 1
        h.attrs['NumPart_ThisFile'] = NumPart
        h.attrs['NumPart_Total'] = NumPart
        h.attrs['NumPart_Total_HighWord'] = np.zeros(nPartTypes, dtype='int64')

        for k in ['Time', 'Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
            h.attrs[k] = 0.0
        for k in ['Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback']:
            h.attrs['Flag_%s' % k] = 0
        h.attrs['Flag_DoublePrecision'] = 1

        if massTable is not None:
            h.attrs['MassTable'] = massTable
        else:
            h.attrs['MassTable'] = np.zeros(nPartTypes, dtype='float64')


make_ics(filename=FilePath, writefunc=write_ic_file)
