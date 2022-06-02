#!/usr/bin/env python3
""" @package ./examples/khi_hydro_2d/check.py
Code that checks results of 2d KHI problem

created by Thomas Berlok, last modified 12.09.2018 -- comments welcome
"""

# load libraries
import sys  # load sys; needed for exit codes
import numpy as np  # load numpy
import h5py  # load h5py; needed to read snapshots
import os  # file specific calls
import os.path
import matplotlib
import matplotlib.pyplot as plt  ## needs to be active for plotting!
matplotlib.rc_file_defaults()

# Needed for fitting exponential growth
from scipy.optimize import curve_fit

# Maximum relative error allowed on growth rate
tol = 1e-2

# Save a plot (needs matplotlib)
makeplots = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print('khi_hydro_2d: checking simulation output in directory ' +
      simulation_directory)

directory = os.path.join(simulation_directory, 'output')

# # Load linear solution
# filename = os.path.join('examples', 'khi_hydro_2d', 'kh-with-delta.p')
# system = pickle.load(open(filename, 'rb'))
# grid = system.grid
# system.grid.make_grid()

# # Get background equilibrium functions
# bg = system.make_background(True)

# Maximally unstable wavenumber
# kx = system.kx

# Frequency of the KHI
omega = 1.4035132749787171 - 0.5422062481760126j
# Growth rate of the KHI
sigma = omega.real

# Box size
kx = 3.5128286
Lx = 2 * np.pi / kx
Ly = 2.0

# Get background equilibrium functions
u0 = 1.0
rho0 = 1.0
T0 = 1.0
p0 = rho0 * T0
delta = 1.0

y1 = 0.5
y2 = 1.5
a = 0.05

# adiabatic index
gamma = 5.0 / 3.0


def v0_an(y):
    """Background velocity shear"""
    return u0 * (np.tanh((y - y1) / a) - np.tanh((y - y2) / a) - 1.0)


def rho0_an(y):
    """Background density profile"""
    return rho0 * (1 + delta / 2 * (np.tanh((y - y1) / a) - np.tanh(
        (y - y2) / a)))


def load_snap(filename):
    """Helper function for loading snapshots"""
    ptName = 'PartType0'
    fields = {}
    with h5py.File(os.path.join(directory, filename), 'r') as f:
        for k in f[ptName].keys():
            fields[k] = f[ptName][k][()]

        time = f['Header'].attrs['Time']
    return fields, time


rms_vx = []
rms_vy = []
rms_rho = []
rms_p = []
time = []

for i_file in range(11):
    filename = 'snap_%03d.hdf5' % i_file
    f, t = load_snap(filename)
    pos = f['CenterOfMass']
    vel = f['Velocities']
    f.update({'vx': vel[:, 0]})
    f.update({'vy': vel[:, 1]})
    f.update({'vx0': v0_an(pos[:, 1])})
    f.update({'rho0': rho0_an(pos[:, 1])})
    f.update({'vol': f['Masses'] / f['Density']})
    vol = np.sum(f['vol'])

    rms_vx.append(np.sum(np.abs(f['vx'] - f['vx0']) * f['vol']) / vol)
    rms_vy.append(np.sum(np.abs(f['vy']) * f['vol']) / vol)
    rms_rho.append(
        np.sum(np.abs(f['Density'] / f['rho0'] - 1) * f['vol']) / vol)
    rms_p.append(
        np.sum(
            np.abs(f['InternalEnergy'] *
                   (gamma - 1) * f['Density'] / p0 - 1) * f['vol']) / vol)
    time.append(t)

time = np.array(time)
rms_vx = np.array(rms_vx)
rms_vy = np.array(rms_vy)
rms_rho = np.array(rms_rho)
rms_p = np.array(rms_p)

first = 0
last = int(10)


# Exponential growth function
def func(x, a, b):
    return a * np.exp(b * x)


def lin_func(x, a, b):
    return a + b * x


def fitting(vals):
    popt, pcov = curve_fit(lin_func, time[first:last],
                           np.log(vals[first:last]))
    sigma_fit = popt[1]
    err = np.abs((sigma_fit - sigma) / sigma)
    return sigma_fit, err


sigma_vx, err_vx = fitting(rms_vx)
sigma_vy, err_vy = fitting(rms_vy)
sigma_rho, err_rho = fitting(rms_rho)
sigma_p, err_p = fitting(rms_p)
# printing results
print('examples/khi_hydro_2d/check.py: relative error on growth fit:')
print('\t vx: %g' % err_vx)
print('\t vy: %g' % err_vy)
print('\t rho: %g' % err_rho)
print('\t p: %g' % err_p)

if makeplots:
    if not os.path.exists(os.path.join(simulation_directory, 'plots')):
        os.mkdir(os.path.join(simulation_directory, 'plots'))

    plt.figure(2)
    plt.clf()
    plt.semilogy(time, rms_vx, label=r'$\delta v_x$')
    plt.semilogy(time, rms_vy, label=r'$\delta v_y$')
    plt.semilogy(time, rms_rho, label=r'$\delta \rho/\rho$')
    plt.semilogy(time, rms_p, label=r'$\delta p/p$')
    plt.semilogy(time[:last],
                 rms_vy[0] * np.exp(sigma_vy * time[:last]),
                 '--',
                 label='Theory')

    plt.legend(frameon=False)
    plt.xlabel(r'$t$')
    plt.savefig(
        os.path.join(simulation_directory, 'plots', 'khi_hydro_2d_growth.pdf'))

    for i_file in [0, 10]:
        filename = 'snap_%03d.hdf5' % i_file
        try:
            data = h5py.File(os.path.join(directory, filename), 'r')
        except (OSError, IOError):
            break

        print(filename)
        VoronoiPos = np.array(data['PartType0']['Coordinates'],
                              dtype=np.float64)
        Velocity = np.array(data['PartType0']['Velocities'], dtype=np.float64)
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
        pc = ax.pcolormesh(Edges1d,
                           Edges1d,
                           Velocity[cells, 1].reshape((Nplot, Nplot)),
                           rasterized=True)
        ax.set_xlim(0, Boxsize)
        ax.set_ylim(0, Boxsize)
        plt.text(0.02,
                 0.95,
                 r'$N=%d,\ t=%3.1f, v_y=\left\{%g,%g\right\}$' %
                 (np.int32(NumberOfCells), Time, Velocity[:, 1].min(),
                  Velocity[:, 1].max()),
                 color='w',
                 transform=fig.transFigure)

        fig.savefig(os.path.join(simulation_directory, 'plots',
                                 'velocity_%03d.pdf' % i_file),
                    dpi=300)
        plt.close(fig)
# criteria for failing the test
if np.max([err_vx, err_vy, err_rho, err_p]) <= tol:
    print('normal exit')
else:
    sys.exit(1)
