""" @package ./examples/brag_decay_3d/check.py
Code that checks results of 3D decay of velocity profile with Braginskii
viscosity

created by Thomas Berlok, last modified 18.07.2019 -- comments welcome
"""

# load libraries
import sys    # load sys; needed for exit codes
import numpy as np    # load numpy
import h5py    # load h5py; needed to read snapshots
import glob

# Save a plot (needs matplotlib)
makeplots = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print("examples/brag_decay_3d/check.py: checking simulation output in directory "
      + simulation_directory)

directory = simulation_directory+"/output/"

""" initial condition parameters """
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


def load_snap(filename):
    """Helper function for loading snapshots"""
    ptName = 'PartType0'
    fields = {}
    with h5py.File(directory + filename, 'r') as f:
        for k in f[ptName].keys():
            fields[k] = f[ptName][k][()]

        time = f['Header'].attrs['Time']
    return fields, time


def rms_diff(x, y):
    """rms-difference between x and y"""
    return np.sqrt(np.mean((x - y)**2.0))


tol0 = 0.027
tol1 = 0.035
msg0 = 'Velocity error ({}) too large'
msg1 = 'Internal energy error ({}) too large'

fnames = glob.glob(directory + '/snap_*.hdf5')

t = []
ak = []
failed = False
for snap in range(1, len(fnames)):
    f, time = load_snap('snap_{:03d}.hdf5'.format(snap))

    x = f['Coordinates'][:, 0]
    y = f['Coordinates'][:, 1]
    z = f['Coordinates'][:, 2]
    vx = f['Velocities'][:, 0]
    vy = f['Velocities'][:, 1]
    vz = f['Velocities'][:, 2]
    Bx = f['MagneticField'][:, 0]
    By = f['MagneticField'][:, 1]
    Bz = f['MagneticField'][:, 2]
    u = f['InternalEnergy']

    s = nx*x + ny*y + nz*z
    v = nx*vx + ny*vy + nz*vz

    index = np.argsort(s)
    s = s[index]
    u = u[index]
    v = v[index]

    err0 = rms_diff(v, v_an(s, time))
    err1 = rms_diff(u, u_an(s, time))

    # print(err0, err1)

    if not err0 <= tol0:
        failed = True
        print(msg0.format(err0))
    if not err1 <= tol1:
        failed = True
        print(msg1.format(err1))

    S = np.sin(k*s)/len(s)
    C = np.cos(k*s)/len(s)

    # Fourier amplitude of perturbation
    ak.append(2*np.sqrt((S*v).sum()**2 + (C*v).sum()**2))

    t.append(time)


if makeplots:
    import os

    vis_path = simulation_directory + '/vis'
    try:
        os.mkdir(vis_path)
    except:
        pass

    from matplotlib import rc
    columnwidth = 336*2.5
    dpi = 72*2.5
    aspect = 1
    figwidth = columnwidth/dpi
    fontsize = 12
    rc('figure', figsize=(figwidth, aspect*figwidth), dpi=dpi)
    rc('savefig', dpi=dpi)
    rc('xtick', direction='out')
    rc('ytick', direction='out')
    rc('axes', labelsize=fontsize, titlesize=fontsize)
    rc('font', size=fontsize)
    rc('legend', fontsize=fontsize)
    rc('xtick', labelsize=fontsize)
    rc('ytick', labelsize=fontsize)
    rc('font', family='serif', serif='cm')
    rc('text', usetex=True)
    rc('text.latex', preamble=[
        r'\usepackage[T1]{fontenc}',
        r'\usepackage[english]{babel}',
        r'\usepackage[utf8]{inputenc}',
        r'\usepackage{lmodern}',
        r'\usepackage{microtype}',
        r'\usepackage{amsmath}',
        r'\usepackage{bm}'])

    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.clf()
    fig, axes = plt.subplots(num=1, sharex=True, nrows=2)

    f, t = load_snap('snap_000.hdf5')

    t = []
    ak = []
    for j in range(0, len(fnames), 1):
        f, time = load_snap('snap_{:03d}.hdf5'.format(j))

        x = f['Coordinates'][:, 0]
        y = f['Coordinates'][:, 1]
        z = f['Coordinates'][:, 2]
        vx = f['Velocities'][:, 0]
        vy = f['Velocities'][:, 1]
        vz = f['Velocities'][:, 2]
        Bx = f['MagneticField'][:, 0]
        By = f['MagneticField'][:, 1]
        Bz = f['MagneticField'][:, 2]
        u = f['InternalEnergy']

        s = nx*x + ny*y + nz*z

        index = np.argsort(s)
        s = s[index]

        v = nx*vx[index] + ny*vy[index] + nz*vz[index]
        u = u[index]

        S = np.sin(k*s)/len(s)
        C = np.cos(k*s)/len(s)

        # Fourier amplitude of perturbation
        ak.append(2*np.sqrt((S*v).sum()**2 + (C*v).sum()**2))
        t.append(time)

        axes[0].plot(s, v, label=r"$c t = {:1.2f}$".format(time))
        axes[0].plot(s, v_an(s, time), 'k--')
        axes[1].plot(s, u, label=r"$c t = {:1.2f}$".format(time))
        axes[1].plot(s, u_an(s, time), 'k--')

    axes[0].set_ylim(-1, 1.8)
    axes[1].set_ylabel(r"$\varepsilon/(\rho c^2)$")
    axes[0].legend(frameon=False, ncol=3)
    axes[1].set_xlabel(r"$x_\parallel/L$")
    axes[0].set_ylabel(r"$v_\parallel/c$")
    # axes[0].set_title("3D Braginskii ({}x{}x{})".format(Nx, Ny, Nz))
    fig.savefig(vis_path+"/3D_Brag_test_profiles.pdf")

    t = np.array(t)

    plt.figure(2)
    plt.clf()
    plt.semilogy(t, ak, label='Simulation')
    plt.semilogy(t, np.exp(-gamma*t), 'k--', label='Theory')
    plt.ylabel(r"$v_k$")
    plt.xlabel(r"$t$")
    plt.legend(frameon=False)
    plt.savefig(vis_path+"/3D_Brag_test_decay.pdf")

# criteria for failing the test
if failed:
    sys.exit(1)

print("normal exit")
