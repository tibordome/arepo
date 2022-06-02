#!/usr/bin/env python3
""" @package ./examples/khi_mhd_2d/check.py
Code that checks results of 2d KHI problem

created by Thomas Berlok, last modified 19.09.2018 -- comments welcome
"""

# load libraries
import sys  # load sys; needed for exit codes
import numpy as np  # load numpy
import h5py  # load h5py; needed to read snapshots
import os  # file specific calls
import os.path
import matplotlib
import matplotlib.pyplot as plt  # needs to be active for plotting!
matplotlib.rc_file_defaults()

# Needed for fitting exponential growth
from scipy.optimize import curve_fit

# Maximum relative error allowed on growth rate
tol = 2.5e-2

# Save a plot (needs matplotlib)
makeplots = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print('examples/khi_mhd_2d/check.py: checking simulation output in directory ' +
      simulation_directory)

directory = os.path.join(simulation_directory, 'output')

# Some parameters
gamma = 5. / 3.
uflow = 1.0
a = 0.05
p0 = 1.0
rho0 = 1.0
beta = 10.0
B0 = np.sqrt(2 * p0 / beta) * np.sqrt(4 * np.pi)
y1 = 0.5
y2 = 1.5

# Maximally unstable wavenumber
kx = 5.4234088072021045

# Growth rate of the KHI
sigma = 1.649689877861002

# Box size
Lx = 2 * np.pi / kx
Ly = 2.0


def vx_an(y):
    """Background shear profile"""
    vx = uflow * (np.tanh((y - y1) / a) - np.tanh((y - y2) / a) - 1.0)
    return vx


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
rms_bx = []
rms_by = []
rms_rho = []
rms_T = []
time = []

for i_file in range(11):
    filename = 'snap_%03d.hdf5' % i_file
    f, t = load_snap(filename)
    pos = f['Coordinates']
    vel = f['Velocities']
    bfld = f['MagneticField']
    f.update({'vx': vel[:, 0]})
    f.update({'vy': vel[:, 1]})
    f.update({'vx0': vx_an(pos[:, 1])})
    f.update({'vol': f['Masses'] / f['Density']})
    vol = np.sum(f['vol'])

    rms_vx.append(np.sum(np.abs(f['vx'] - f['vx0']) * f['vol']) / vol)
    rms_vy.append(np.sum(np.abs(f['vy']) * f['vol']) / vol)
    rms_bx.append(np.sum(np.abs(bfld[:, 0] - B0) * f['vol']) / vol)
    rms_by.append(np.sum(np.abs(bfld[:, 1]) * f['vol']) / vol)
    rms_rho.append(np.sum(np.abs(f['Density'] / rho0 - 1) * f['vol']) / vol)
    rms_T.append(
        np.sum(np.abs(f['InternalEnergy'] * (gamma - 1) - 1) * f['vol']) / vol)
    time.append(t)

time = np.array(time)
rms_vx = np.array(rms_vx)
rms_vy = np.array(rms_vy)
rms_bx = np.array(rms_bx)
rms_by = np.array(rms_by)
rms_rho = np.array(rms_rho)
rms_T = np.array(rms_T)

first = 1
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
sigma_bx, err_bx = fitting(rms_bx)
sigma_by, err_by = fitting(rms_by)
sigma_rho, err_rho = fitting(rms_rho)
sigma_p, err_T = fitting(rms_T)
# printing results
print('examples/khi_mhd_2d/check.py: relative error on growth fit:')
print('\t vx: %g' % err_vx)
print('\t vy: %g' % err_vy)
print('\t bx: %g' % err_bx)
print('\t by: %g' % err_by)
print('\t rho: %g' % err_rho)
print('\t T: %g' % err_T)

if makeplots:
    if not os.path.exists(os.path.join(simulation_directory, 'plots')):
        os.mkdir(os.path.join(simulation_directory, 'plots'))

    plt.figure(2)
    plt.clf()
    plt.semilogy(time, rms_vx, label=r'$\delta v_x$')
    plt.semilogy(time, rms_vy, label=r'$\delta v_y$')
    plt.semilogy(time, rms_bx, label=r'$\delta b_x$')
    plt.semilogy(time, rms_by, label=r'$\delta b_y$')
    plt.semilogy(time, rms_rho, label=r'$\delta \rho/\rho$')
    plt.semilogy(time, rms_T, label=r'$\delta T/T$')
    plt.semilogy(time[:last],
                 rms_vy[0] * np.exp(sigma_vy * time[:last]),
                 '--',
                 label='Theory')

    plt.legend(frameon=False)
    plt.xlabel(r'$t$')
    plt.savefig(
        os.path.join(simulation_directory, 'plots', 'khi_mhd_2d_growth.pdf'))

    for i_file in [0, 10]:
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
        Edges1dx = np.linspace(0.,
                               Lx,
                               Nplot + 1,
                               endpoint=True,
                               dtype=np.float64)
        Grid1dx = 0.5 * (Edges1dx[1:] + Edges1dx[:-1])
        Edges1dy = np.linspace(0.,
                               Ly,
                               Nplot + 1,
                               endpoint=True,
                               dtype=np.float64)
        Grid1dy = 0.5 * (Edges1dy[1:] + Edges1dy[:-1])
        xx, yy = np.meshgrid(Grid1dx, Grid1dy)
        Grid2D = np.array([xx.reshape(Nplot**2), yy.reshape(Nplot**2)]).T
        dist, cells = spatial.KDTree(VoronoiPos[:, :2]).query(Grid2D, k=1)

        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = plt.axes([0, 0, 1, 1])
        pc = ax.pcolormesh(Edges1dx,
                           Edges1dy,
                           MagneticField[cells, 1].reshape((Nplot, Nplot)),
                           rasterized=True)
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly)
        plt.text(0.02,
                 0.95,
                 r'$N=%d,\ t=%3.1f, b_y=\left\{%g,%g\right\}$' %
                 (np.int32(NumberOfCells), Time, MagneticField[:, 1].min(),
                  MagneticField[:, 1].max()),
                 color='w',
                 transform=fig.transFigure)

        fig.savefig(os.path.join(simulation_directory, 'plots',
                                 'magneticfield_%03d.pdf' % i_file),
                    dpi=300)
        plt.close(fig)
# criteria for failing the test
if np.max([err_vx, err_vy, err_bx, err_by, err_rho, err_T]) <= tol:
    print('normal exit')
else:
    sys.exit(1)
