""" @package ./examples/brag_2d/create.py
Create initial conditions for 2D anisotropic diffusion
with Braginskii viscosity.

The test is described in the code paper (see Figure 2 therein)

Braginskii viscosity on an unstructured, moving mesh -- accelerated with
super-time-stepping by T. Berlok, R. Pakmor & and C. Pfrommer (2020)

https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2919B/abstract


created by Thomas Berlok, last modified 24.11.2020 -- comments welcome
"""

# load libraries
import sys    # load sys; needed for exit codes
import numpy as np    # load numpy
import h5py    # load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/brag_2d/create.py: creating ICs in directory " +
      simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/ics.hdf5'

hexagonal = True

Lx = 1.0
gamma = 5/3.0  # adiabatic index

Nx = 128
Ny = 4
Ly = Lx*Ny/Nx
dx = Lx / Nx
dy = Ly / Ny
print(Ly)

ampl = 1e-1

rho0 = 1.0
p0 = 1.0

beta = 1e4

# β = 2 cs²/va²
# cs² = p/ρ
# va² = B²/4πρ
# B² = 2 cs²/β 4πρ = 8 π p/β

B0 = np.sqrt(8*np.pi*p0/beta)
kx = 2*np.pi/Lx
va = B0/np.sqrt(4*np.pi*rho0)
omega = kx*va
tau_a = 1/omega


# Velocity
def q(x):
    from scipy.special import erf
    # from math import erf
    return 1.5 - 0.5*(erf((x-0.75)/0.05) - erf((x-0.25)/0.05))


def rho_an(x, y):
    return rho0*np.ones_like(x)


def bx_an(x, y):
    return B0/np.sqrt(2)*np.ones_like(x)


def by_an(x, y):
    return B0/np.sqrt(2)*np.ones_like(x)


def p_an(x, y):
    return p0*np.ones_like(x)


def Az_an(x, y):
    return np.zeros_like(x)


def create_ics(filename='ics.hdf5', writefunc=None):

    # allocate
    pos = np.zeros((Nx*Ny, 3), dtype='float64')
    vel = np.zeros((Nx*Ny, 3), dtype='float64')
    bfld = np.zeros((Nx*Ny, 3), dtype='float64')
    A = np.zeros((Nx*Ny, 3), dtype='float64')
    dens = np.zeros(Nx*Ny, dtype='float64')
    u = np.zeros(Nx*Ny, dtype='float64')
    id = np.arange(Nx*Ny, dtype='int64') + 1

    # assign gas cell properties

    x = np.arange(Nx)*dx + dx/2
    y = np.arange(Ny)*dy + dy/2
    xx, yy = np.meshgrid(x, y)
    if hexagonal:
        for j in range(Ny):
            xx[j, :] += (j % 2)*0.45*dx

    # IPython.embed()
    x = xx.flatten()
    y = yy.flatten()
    index = range(Nx*Ny)

    rho = rho_an(x, y)
    P = p_an(x, y)

    bx = bx_an(x, y)
    by = by_an(x, y)

    Az = Az_an(x, y)

    pos[index, 0] = x
    pos[index, 1] = y
    pos[index, 2] = 0.0

    vel[index, 1] = q(x)

    dens[index] = rho
    u[index] = P/rho/(gamma-1.0)

    bfld[index, 0] = bx
    bfld[index, 1] = by

    A[index, 2] = Az

    pt0 = {'Coordinates': pos, 'Velocities': vel, 'Masses': dens,
           'MagneticField': bfld, 'MagneticVectorPotential':
           A/np.sqrt(4*np.pi),
           'InternalEnergy': u, 'ParticleIDs': id}

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
sys.exit(0)
