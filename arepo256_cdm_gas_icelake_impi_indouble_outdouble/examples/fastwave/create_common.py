#!/usr/bin/env python3
""" @package ./examples/fastwave/create_common.py
Create initial conditions for 2D fast magnetosonic wave
with Braginskii viscosity.

The tests differ from fastwave_2d mainly in the type of timestepping employed.
- fastwave_local_timestepping_2d:
  Here we use local time steps.
- fastwave_subcycling_2d:
  Here we use simple subcycling with local time steps.
- fastwave_RKL2_timestepping_2d:
  Here we use RKL2 timestepping with global time steps.

Standard MHD can also be tested by setting nu = 0
in this script and commenting out the Braginskii flags in Config.sh
and the Braginskii inputs in param.txt

The test is described in the code paper (see Figure 3 therein)

Braginskii viscosity on an unstructured, moving mesh -- accelerated with
super-time-stepping by T. Berlok, R. Pakmor & and C. Pfrommer (2020)

https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2919B/abstract


created by Thomas Berlok, last modified 24.11.2020 -- comments welcome
"""

# load libraries
import os.path
import sys
import numpy as np
import h5py    # load h5py; needed to write initial conditions in hdf5 format

# make sure that this file's directory is in sys.path
# (see also https://bugs.python.org/issue17639)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import create_check_parameters


simulation_directory = str(sys.argv[1])
print(__file__ + ": creating ICs in directory " +
      simulation_directory)

# initial condition parameters
FilePath = simulation_directory + '/ics.hdf5'

# Viscosity coefficient
nu = create_check_parameters.nu

# For high nu solutions, choose between 1 and -1
# to test the two different solutions
sign = -1

# Standing or traveling wave
traveling = False

# Grid configuration
hexagonal = True
Nx = 32
Ny = 32

# Grid size in x- and y-direction
Lx = 1.0
Ly = Lx*Ny/Nx

dx = Lx / Nx
dy = Ly / Ny

# Dimensionless amplitude of perturbation
A = 1e-3

# Wavenumbers
ikx = 1
iky = create_check_parameters.iky

gamma = 5.0/3.0

rho0 = 1.0
T0 = 1.0
P0 = rho0*T0

# Sound speed
cs = np.sqrt(gamma*P0/rho0)

# Wave vector and its modulus
kx = 2*np.pi*ikx/Lx
ky = 2*np.pi*iky/Ly
k = np.sqrt(kx*kx + ky*ky)
kperp = k

nx = kx/k
ny = ky/k

(Bx0, By0, Bz0) = (0, 0, 1)
B20 = Bx0**2 + By0**2 + Bz0**2
va = np.sqrt(B20/(4*np.pi*rho0))


vph = np.sqrt(cs**2 + va**2)

high_nu_limit = cs**2 + va**2 < (k*nu/6)**2

omega = k*np.sqrt(vph**2 - (k*nu/6)**2 + 0j) - 1j*k**2.0*nu/6
if high_nu_limit:
    omega *= sign

# Real and imaginary part
omega0 = omega.real
gam = -omega.imag


if traveling or high_nu_limit:
    def rho_an(s, t):
        """Analytical density"""
        return rho0*(1 + A*np.cos(kperp*s)*np.exp(-gam*t))

    def Bz_an(s, t):
        """Analytical z-component of B-field"""
        return Bz0*(1 + A*np.cos(kperp*s)*np.exp(-gam*t))

    def vperp_an(s, t):
        """Velocity in direction of kperp"""
        fun = A*np.sin(kperp*s)*gam/kperp*np.exp(-gam*t)
        return fun

    # def rho_an(s, t):
    #     """Analytical density"""
    #     return rho0*(1 + A*np.exp(1j*kperp*s - 1j*omega*t).real)

    # def Bz_an(s, t):
    #     """Analytical z-component of B-field"""
    #     return Bz0*(1 + A*np.exp(1j*kperp*s - 1j*omega*t).real)

    # def vperp_an(s, t):
    #     """Velocity in direction of kperp"""
    #     fun = (A*omega/kperp*np.exp(1j*kperp*s - 1j*omega*t)).real
    #     return fun

# Standing wave
else:

    def rho_an(s, t):
        """Analytical density"""
        return rho0*(1 + A*np.cos(kperp*s)*np.sin(omega0*t)*np.exp(-gam*t))

    def Bz_an(s, t):
        """Analytical z-component of B-field"""
        return Bz0*(1 + A*np.cos(kperp*s)*np.sin(omega0*t)*np.exp(-gam*t))

    def vperp_an(s, t):
        """Velocity in direction of kperp"""
        fun = -A*np.sin(kperp*s)*(omega0*np.cos(omega0*t) -
                                  gam*np.sin(omega0*t))
        fun *= np.exp(-gam*t)/kperp
        return fun


def vx_an(s, t):
    """Analytic x-velocity"""
    return vperp_an(s, t)*nx


def vy_an(s, t):
    """Analytic y-velocity"""
    return vperp_an(s, t)*ny


def create_ics(filename='ics.hdf5', writefunc=None):
    if create_check_parameters.points_file:
        # Load grid points
        dat = np.loadtxt(os.path.join(
            simulation_directory, create_check_parameters.points_file
        ))
        N = dat.shape[0]
        x = dat[:, 0]
        y = dat[:, 1]
    else:
        N = Nx * Ny
        x = np.arange(Nx)*dx + dx/2
        y = np.arange(Ny)*dy + dy/2
        xx, yy = np.meshgrid(x, y)
        if hexagonal:
            for j in range(Ny):
                xx[j, :] += (j % 2)*0.45*dx

        x = xx.flatten()
        y = yy.flatten()

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    bfld = np.zeros((N, 3), dtype='float64')
    dens = np.zeros(N, dtype='float64')
    u = np.zeros(N, dtype='float64')
    id = np.arange(N, dtype='int64') + 1

    s = nx*x + ny*y

    t = 0.0

    pos[:, 0] = x
    pos[:, 1] = y

    vel[:, 0] = vx_an(s, t)
    vel[:, 1] = vy_an(s, t)

    bfld[:, 0] = 0.0
    bfld[:, 1] = 0.0
    bfld[:, 2] = Bz_an(s, t)

    dens[:] = rho_an(s, t)
    P = P0 + (dens - rho0)/rho0*cs**2
    u[:] = P / dens / (gamma-1)

    pt0 = {'Coordinates': pos, 'Velocities': vel, 'Masses': dens,
           'MagneticField': bfld, 'InternalEnergy': u, 'ParticleIDs': id}

    writefunc(filename, {'PartType0': pt0}, boxSize=Lx)
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
        h.attrs['NumPart_Total_HighWord'] = np.zeros(nPartTypes,
                                                     dtype='int64')

        for k in ['Time', 'Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
            h.attrs[k] = 0.0
        for k in ['Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback']:
            h.attrs['Flag_%s' % k] = 0
        h.attrs['Flag_DoublePrecision'] = 1

        if massTable is not None:
            h.attrs['MassTable'] = massTable
        else:
            h.attrs['MassTable'] = np.zeros(nPartTypes, dtype='float64')


create_ics(filename=FilePath, writefunc=write_ic_file)
