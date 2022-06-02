""" @package ./examples/brag_decay_3d/create.py
Create initial conditions for 3D decay of velocity profile with Braginskii
viscosity but *without* hydro update.

The test is described in the code paper (see Figure 1 therein)

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
print("examples/brag_decay_3d/create.py: creating ICs in directory " +
      simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/ics.hdf5'

# Grid size in x- and y-direction
Lx = 1.0
Ly = 1.0
Lz = 1.0
hexagonal = True
Nx = 16
Ny = 16
Nz = 16
dx = Lx / Nx
dy = Ly / Ny
dz = Lz / Nz

# Dimensionless amplitude of perturbation
A = 1

# Wavenumbers
ikx = 0
iky = 0
ikz = 1

rho = 1.0
T = 1.0
P = rho*T

# Wave vector and its modulus
kx = 2*np.pi*ikx/Lx
ky = 2*np.pi*iky/Ly
kz = 2*np.pi*ikz/Lz
k = np.sqrt(kx*kx + ky*ky + kz*kz)

nx = kx/k
ny = ky/k
nz = kz/k

nu = 1e-3
u0 = 1
rho = 1.0

gamma = 4/3*nu*k**2.


def v_an(s, t):
    return A*np.sin(k*s)*np.exp(-gamma*t)


def u_an(s, t):
    return u0 + 0.5*rho*A**2*np.cos(k*s)**2*(1.0 - np.exp(-2*gamma*t))


def create_ics(filename='ics.hdf5', writefunc=None):

    # allocate
    pos = np.zeros((Nx*Ny*Nz*2, 3), dtype='float64')
    vel = np.zeros((Nx*Ny*Nz*2, 3), dtype='float64')
    bfld = np.zeros((Nx*Ny*Nz*2, 3), dtype='float64')
    dens = np.zeros(Nx*Ny*Nz*2, dtype='float64')
    u = np.zeros(Nx*Ny*Nz*2, dtype='float64')
    id = np.arange(Nx*Ny*Nz*2, dtype='int64') + 1
    index = range(Nx*Ny*Nz*2)

    x = np.arange(Nx)*dx + dx*0.01
    y = np.arange(Ny)*dy + dy*0.01
    z = np.arange(Nz)*dz + dz*0.01
    xx, yy, zz = np.meshgrid(x, y, z)

    x = xx.flatten()
    y = yy.flatten()
    z = zz.flatten()

    x = np.hstack([x, x + dx*0.51])
    y = np.hstack([y, y + dy*0.51])
    z = np.hstack([z, z + dz*0.51])

    np.random.seed(7)
    x += 2*(np.random.random(Nx*Ny*Nz*2) - 0.5)*0.*dx
    y += 2*(np.random.random(Nx*Ny*Nz*2) - 0.5)*0.*dx
    z += 2*(np.random.random(Nx*Ny*Nz*2) - 0.5)*0.*dx

    pos[index, 0] = x
    pos[index, 1] = y
    pos[index, 2] = z

    dens[index] = rho
    u[index] = 1.0
    vel[index, 0] = nx*np.sin(kx*x + ky*y + kz*z)
    vel[index, 1] = ny*np.sin(kx*x + ky*y + kz*z)
    vel[index, 2] = nz*np.sin(kx*x + ky*y + kz*z)
    bfld[index, 0] = nx
    bfld[index, 1] = ny
    bfld[index, 2] = nz

    pt0 = {'Coordinates': pos, 'Velocities': vel, 'Masses': dens,
           'MagneticField': bfld,
           'InternalEnergy': u, 'ParticleIDs': id}

    writefunc(filename, {'PartType0': pt0}, boxSize=Ly)
    return True
# def create_ics(filename='ics.hdf5', writefunc=None):

#     # allocate
#     pos = np.zeros((Nx*Ny*Nz, 3), dtype='float64')
#     vel = np.zeros((Nx*Ny*Nz, 3), dtype='float64')
#     bfld = np.zeros((Nx*Ny*Nz, 3), dtype='float64')
#     dens = np.zeros(Nx*Ny*Nz, dtype='float64')
#     u = np.zeros(Nx*Ny*Nz, dtype='float64')
#     id = np.arange(Nx*Ny*Nz, dtype='int64') + 1

#     # derived properties
#     dx = Lx / Nx
#     dy = Ly / Ny
#     dz = Lz / Nz

#     index = 0
#     # assign gas cell properties
#     for ix in range(Nx):
#         x = (0.5+ix) * dx
#         for iy in range(Ny):
#             y = (0.5+iy) * dy + (ix % 2) * 0.45 * dy
#             for iz in range(Nz):
#                 z = (0.5+iz) * dz + (iy % 2) * 0.45 * dz

#                 pos[index, 0] = x
#                 pos[index, 1] = y
#                 pos[index, 2] = z

#                 s = nx*x + ny*y + nz*z

#                 vel[index, 0] = nx*v_an(s, t=0)
#                 vel[index, 1] = ny*v_an(s, t=0)
#                 vel[index, 2] = nz*v_an(s, t=0)

#                 bfld[index, 0] = nx
#                 bfld[index, 1] = ny
#                 bfld[index, 2] = nz

#                 dens[index] = rho
#                 u[index] = u0

#                 index += 1

#     pt0 = {'Coordinates': pos, 'Velocities': vel, 'Masses': dens,
#            'MagneticField': bfld, 'InternalEnergy': u, 'ParticleIDs': id}

#     writefunc(filename, {'PartType0': pt0}, boxSize=Lx)
#     return True


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
