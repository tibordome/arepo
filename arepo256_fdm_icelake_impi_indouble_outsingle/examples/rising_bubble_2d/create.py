""" @package ./examples/rising_bubble_2d/create.py
Create initial conditions for 2D bubble problem.

This example considers the buoyant rise of a bubble in an isothermal
intracluster medium.

The bubble dynamics in 2D differ substantially from the 3D case.
Additionally, bubbles that are inflated by a jet (instead of simply
inserted as an under-dense region in the ICs) differ substantially
from this simple setup. This test should therefore be regarded as a
minimal working example (MWE) of the influence of Braginskii
viscosity on bubble dynamics. More details, references, and
3D simulations of this setup can be found in Oliver Franke's
MSc thesis.

The primary motivation is to optimize and stress test the Braginskii
implementation with RKL2 and local time steps. We use a very
large viscosity coefficient in order to make the Braginskii simulation
expensive and the corresponding ideal MHD simulation cheap.

Due to the high ν_∥, the morphology of the bubble differs a lot
between ideal MHD and Braginskii viscosity. This is clearly seen by
eye. We compare the dye entropy of the passive scalar which follows the
bubble in the simulation. This is defined as (Lecoanet et al 2016)

    S = -∫ c ln(c) ρ dV

The ideal MHD simulation gives for the 10 snapshots,

S =

0    0.03429394377850212
1    0.04188833004392935
2    0.049859951889823326
3    0.05695533924961305
4    0.06512466961158118
5    0.07315008496007883
6    0.08113528332583417
7    0.0888704506401326
8    0.09672110204793481
9    0.10399294550392374
10   0.1102453905712064

A Braginskii run performed on 28/11/2020 gives S

0   0.03429394377850212
1   0.03859151895217329
2   0.04313806998245085
3   0.04686904841837425
4   0.050088266867879205
5   0.05263918205165713
6   0.055010710757032545
7   0.05766591940753624
8   0.06021175334789205
9   0.0626810758929754
10  0.06497914061457936

The Braginskii run has S < 0.07 at all times while ideal MHD has
S > 0.07 already at snapshot 5. This is a quantitative estimate of
the amount of mixing in the simulation.

The Braginskii test passes if S < 0.075 in the final snapshot.
The amount of mixing is not converged at this resolution,
so the criterion for a passed test is not as rigorous for this test as
for the other Braginskii tests.


created by Thomas Berlok, last modified 30.11.2020 -- comments welcome
"""

# load libraries
import sys    # load sys; needed for exit codes
import numpy as np    # load numpy
import h5py    # load h5py; needed to write initial conditions in hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/rising_bubble_2d/create.py: creating ICs in directory " +
      simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/ics.hdf5'

res = 5000 # Number of cells.

# For now, IC's with 10000 generating points not uploaded
# (it is 0.3 MB...)
assert res in [2000, 5000]

# Make sure it is the same in Arepo!
BoxSize = 6.

# 2D or 3D initial conditions?
dims = 2

# MHD simulation or not?
mhd = True
(bx, by, bz) = (1, 0, 0)  # Magnetic field direction
assert np.allclose(np.sqrt(bx**2 + by**2 + bz**2), 1.0), 'b not unit vector!'

# Whether to show mesh
show = False

# Whether to enable outer boundaries or not
use_boundaries = True
# Fraction of radius considered active region when above is true
fraction = 0.92

# Physical parameters
gamma = 5/3.0
rho0 = 1
cs = 1
r0 = 1
beta0 = 1e8  # plasma-beta at r=0
va = cs*np.sqrt(2/beta0)  # Alfven speed at r=0
B = va*np.sqrt(4*np.pi*rho0)  # Magnetic field strength

# Bubble parameters
set_bubble = True
xb1, yb1, zb1 = 0, 0.25, 0  # Center of bubble 1
xb2, yb2, zb2 = 0, -0.25, 0  # Center of bubble 2
rb = 0.25  # bubble radius
rhob = 0.01*rho0  # bubble density
a = 0.1*rb  # Smoothing parameter


def rho_an(r):
    """density as a function of radius"""
    return rho0*(1 + (r/r0)**2)**(-3/4)


def drhodr_an(r):
    """density derivative wrt r as a function of radius"""
    return -3/2*r/r0**2*rho0*(1 + (r/r0)**2)**(-7/4)


def p_an(r):
    """pressure as a function of radius"""
    return cs**2*rho_an(r)


def dpdr_an(r):
    """pressure derivative wrt r as a function of radius"""
    return cs**2*drhodr_an(r)


def g_an(r):
    """Gravitational acceleration as a function of radius.
       Derived from dp/dr = - rho*g
    """
    # return -cs**2*drhodr(r)/rho(r)
    return -3/2*cs**2*r/(r**2 + r0**2)


def phi_an(r):
    """Gravitational potential as a function of radius.
    """
    return 3/4*cs**2*np.log(r**2 + r0**2)


def create_ics(filename='ics.hdf5', writefunc=None):

    # Create grid
    dat = np.loadtxt(simulation_directory + '/bubble_points_{:}.txt'.format(res))
    x = dat[:, 0] - BoxSize/2
    y = dat[:, 1] - BoxSize/2
    z = np.zeros_like(x)

    # Radius
    r = np.sqrt(x**2 + y**2 + z**2)

    if use_boundaries:
        # Remove points outside the box
        # ind = r <= BoxSize/2
        # x = x[ind]
        # y = y[ind]
        # z = z[ind]
        # r = r[ind]

        # Find boolean arrays for inside and outside active computational domain
        # Here active domain hardcoded to be a fraction of the radius of the box
        ind = r < BoxSize/2*fraction
        not_ind = np.logical_not(ind)

        # Count number of boundary points
        Nbound = np.sum(not_ind)

        # Sort arrays with boundary points first, followed by active points
        x = np.hstack([x[not_ind], x[ind]])
        y = np.hstack([y[not_ind], y[ind]])
        z = np.hstack([z[not_ind], z[ind]])
        r = np.hstack([r[not_ind], r[ind]])

    # Total number of grid points
    N = len(x)

    # allocate
    pos = np.zeros((N, 3), dtype='float64')
    vel = np.zeros((N, 3), dtype='float64')
    dens = np.zeros(N, dtype='float64')
    u = np.zeros(N, dtype='float64')
    scalar = np.zeros(N, dtype='float64')
    index = np.arange(N, dtype='int64')

    # Plot mesh
    if show:
        import matplotlib.pyplot as plt
        fig = plt.figure(1)
        plt.clf()
        plt.plot(x, y, '.')
        if use_boundaries:
            plt.plot(x[:Nbound], y[:Nbound], '.')

        plt.show()

    # Set x, y, and z positions, density, and internal energy
    rho = rho_an(r)
    P = p_an(r)

    # Reduce density inside bubbles
    if set_bubble:
        dist1 = np.sqrt((x-xb1)**2 + (y-yb1)**2 + (z-zb1)**2)
        dist2 = np.sqrt((x-xb2)**2 + (y-yb2)**2 + (z-zb2)**2)

        def density_profile(r1, r2):
            """
            Density profile which smoothly changes from analytical profile
            to reduced values inside the bubbles.
            r1 is distance to the closest bubble center
            r2 is distance from cluster center
            """
            # Analytical expression for the density
            rho = rhob + (1 + np.tanh((r1 - rb)/a))*(rho_an(r2)-rhob)/2
            return rho

        def scalar_profile(r1):
            """
            Profile which smoothly changes from 1 to 0 outside the bubble.
            Used to track the bubble material and how it mixes with the
            surroundings.
            r1 is distance to the closest bubble center
            """
            scalar = (1 + np.tanh(-(r1-rb)/a))/2
            return scalar

        dist = np.min(np.vstack([dist1, dist2]), axis=0)
        rho = density_profile(dist, r)
        scalar = scalar_profile(dist)

    pos[index, 0] = x + BoxSize/2
    pos[index, 1] = y + BoxSize/2
    pos[index, 2] = z + BoxSize/2

    dens[index] = rho
    u[index] = P/rho/(gamma-1.0)

    if use_boundaries:
        # Indices of boundary points will be between 0 and 10000000.
        # Active physical points will have index larger than 10000000
        index[Nbound:] += 10000000 - Nbound + 1

    pt0 = {'Coordinates': pos, 'Velocities': vel, 'Masses': dens, 'rho': dens,
           'InternalEnergy': u, 'PassiveScalars': scalar, 'ParticleIDs': index}

    # Set magnetic field
    if mhd:
        bfld = pos = np.zeros((N, 3), dtype='float64')
        bfld[:, :] = np.array([bx, by, bz])*B
        pt0.update({'MagneticField': bfld})

    writefunc(filename, {'PartType0': pt0}, boxSize=BoxSize)
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
