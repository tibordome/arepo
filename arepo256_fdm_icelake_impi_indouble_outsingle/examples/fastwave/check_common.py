#!/usr/bin/env python3
""" @package ./examples/fastwave/check_common.py
Code that checks results of 2D fast magnetosonic wave

created by Thomas Berlok, last modified 16.07.2019 -- comments welcome
"""

# load libraries
import glob
import os.path
import sys
import numpy as np    # load numpy
import h5py    # load h5py; needed to read snapshots

# make sure that this file's directory is in sys.path
# (see also https://bugs.python.org/issue17639)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import create_check_parameters


# Maximum relative error allowed on decay rate
tol = 2.5e-2

# Save a plot (needs matplotlib)
makeplots = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print(__file__ + ": checking simulation output in directory "
      + simulation_directory)

directory = simulation_directory+"/output/"


# Viscosity coefficient
nu = create_check_parameters.nu

# For high nu solutions, choose between 1 and -1
# to test the two different solutions
sign = -1

# Standing or traveling wave
traveling = False

# Grid configuration
hexagonal = True
Nx = 64
Ny = 64

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
B0 = np.sqrt(B20)
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


def load_snap(filename):
    """Helper function for loading snapshots"""
    ptName = 'PartType0'
    fields = {}
    with h5py.File(directory + filename, 'r') as f:
        for k in f[ptName].keys():
            fields[k] = f[ptName][k][()]

        time = f['Header'].attrs['Time']
    return fields, time


def get_ak():
    fnames = glob.glob(directory + '/snap_*.hdf5')

    t = []
    ak = []
    for snap in range(1, len(fnames)):
        f, time = load_snap('snap_{:03d}.hdf5'.format(snap))

        x = f['Coordinates'][:, 0]
        y = f['Coordinates'][:, 1]
        vx = f['Velocities'][:, 0]
        vy = f['Velocities'][:, 1]
        Bz = f['MagneticField'][:, 2]
        rho = f['Density']

        s = nx*x + ny*y
        v = nx*vx + ny*vy

        index = np.argsort(s)
        s = s[index]
        Bz = Bz[index]
        v = v[index]
        rho = rho[index] - 1

        S = np.sin(k*s)/len(s)
        C = np.cos(k*s)/len(s)

        # Fourier amplitude of perturbation
        ak.append(2*np.sqrt((S*rho).sum()**2 + (C*rho).sum()**2))

        t.append(time)
    return (np.array(t), np.array(ak))


def rms_diff(x, y):
    """rms-difference between x and y"""
    return np.sqrt(np.mean((x - y)**2.0))


def fit_gamma_and_omega0(t, ak):
    from scipy.optimize import curve_fit

    if traveling or high_nu_limit:
        omega0 = 0
    else:
        from scipy.signal import argrelextrema
        # Find first local maximum
        index = argrelextrema(ak, np.greater)
        t = t[index]
        ak = ak[index]
        omega0 = np.mean(2*np.pi/(2*np.diff(t)))

    # Exponential growth function
    def func(x, a, b):
        return a*np.exp(b*x)

    def lin_func(x, a, b):
        return a + b*x

    # popt, pcov = curve_fit(lin_func, t, np.log(ak))
    popt, pcov = curve_fit(func, t, ak)

    return -popt[1], omega0


tol0 = 3e-5
tol1 = 3e-5
tol2 = create_check_parameters.tol2
msg0 = 'Density error ({}) too large'
msg1 = 'Magnetic field error ({}) too large'
msg2 = 'Velocity error ({}) too large'

failed = False

fnames = glob.glob(directory + '/snap_*.hdf5')
ak_v = []
ak_rho = []
ak_B = []
time = []
for i in range(1, len(fnames)):
    f, t = load_snap('snap_{:03d}.hdf5'.format(i))

    x = f['Coordinates'][:, 0]
    y = f['Coordinates'][:, 1]

    vx = f['Velocities'][:, 0]
    vy = f['Velocities'][:, 1]
    Bz = f['MagneticField'][:, 2]
    rho = f['Density']

    s = nx*x + ny*y

    ind = np.argsort(s)

    s = nx*x + ny*y
    v = nx*vx + ny*vy

    index = np.argsort(s)
    s = s[index]
    Bz = Bz[index]
    v = v[index]
    rho = rho[index]

    S = np.sin(k*s)/len(s)
    C = np.cos(k*s)/len(s)

    # Fourier amplitude of perturbation
    ak_rho.append(2*np.sqrt((S*(rho-1)).sum()**2 + (C*(rho-1)).sum()**2))
    ak_v.append(2*np.sqrt((S*v).sum()**2 + (C*v).sum()**2))
    ak_B.append(2*np.sqrt((S*(Bz-B0)).sum()**2 + (C*(Bz-B0)).sum()**2))

    time.append(t)

    err0 = rms_diff(rho, rho_an(s, t))
    err1 = rms_diff(Bz, Bz_an(s, t))
    err2 = rms_diff(vx, vx_an(s, t))

    if err0 > tol0:
        failed = True
        print(msg0.format(err0))
    if err1 > tol1:
        failed = True
        print(msg1.format(err1))
    if err2 > tol2:
        print(msg2.format(err2))
        failed = True

    if np.abs(f['MagneticFieldDivergence']).max() > 1e-16:
        failed = True

ak_v = np.array(ak_v)
ak_rho = np.array(ak_rho)
ak_B = np.array(ak_B)

time = np.array(time)
gam_rho, omega0_rho = fit_gamma_and_omega0(time, ak_rho)
gam_v, omega0_v = fit_gamma_and_omega0(time, ak_v)
gam_B, omega0_B = fit_gamma_and_omega0(time, ak_B)

print('Theory', gam, omega0)
print('fits give')
print(gam_rho, gam_v, gam_B)
print(omega0_rho, omega0_v, omega0_B)


if makeplots:
    import os
    import subprocess

    vis_path = simulation_directory + '/vis'
    try:
        os.mkdir(vis_path)
    except:
        pass

    from matplotlib import rc
    columnwidth = 336*2.5
    dpi = 72*2.5
    aspect = 2
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
    fig, axes = plt.subplots(num=1, sharex=True, nrows=4)

    f, t = load_snap('snap_000.hdf5')

    x = f['Coordinates'][:, 0]
    y = f['Coordinates'][:, 1]
    vx = f['Velocities'][:, 0]
    Bz = f['MagneticField'][:, 2]
    rho = f['Density']
    u = f['InternalEnergy']

    s = nx*x + ny*y

    ind = np.argsort(s)

    im0 = axes[0].plot(s[ind], rho[ind], s[ind], rho_an(s[ind], t), '--')
    im1 = axes[1].plot(s[ind], Bz[ind], s[ind], Bz_an(s[ind], t), '--')
    im2 = axes[2].plot(s[ind], vx[ind], s[ind], vx_an(s[ind], t), '--')
    im3, = axes[3].plot(s[ind], u[ind])

    A = 1e-3
    axes[0].set_ylim(-A, A)
    axes[1].set_ylim(-A, A)

    for ax in axes:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axes[2].set_xlabel(r"$s$")

    axes[0].set_ylabel(r"$\rho/\rho$")
    axes[1].set_ylabel(r"$\delta B_z/B_z$")
    axes[2].set_ylabel(r"$v_s$")

    axes[0].legend([im0[0], im0[1]], ["Simulation", "Theory"], frameon=False,
                   loc="upper right")

    for i in range(1, len(fnames)):
        f, t = load_snap('snap_{:03d}.hdf5'.format(i))

        x = f['Coordinates'][:, 0]
        y = f['Coordinates'][:, 1]

        vx = f['Velocities'][:, 0]
        Bz = f['MagneticField'][:, 2]
        rho = f['Density']
        u = f['InternalEnergy']

        s = nx*x + ny*y

        ind = np.argsort(s)

        im0[0].set_data(s[ind], rho[ind] - rho0)
        im1[0].set_data(s[ind], Bz[ind] - Bz0)
        im2[0].set_data(s[ind], vx[ind])
        im3.set_data(s[ind], u[ind])
        axes[3].set_ylim(u.min(), u.max())

        im0[1].set_data(s[ind], rho_an(s[ind], t) - rho0)
        im1[1].set_data(s[ind], Bz_an(s[ind], t) - Bz0)
        im2[1].set_data(s[ind], vx_an(s[ind], t))

        axes[0].set_title(r"$t={:1.2f}$".format(t))

        fig.savefig(vis_path + "/frame_{}.png".format(i))

    # if ffmpeg exists, make a movie
    try:
        cmd = ['ffmpeg',
               '-framerate',
               '10',
               '-i',
               'frame_%d.png',
               '-pix_fmt',
               'yuv420p',
               'movie_nu_{}_sign{}.mp4'.format(nu, sign),
               '-y']
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT,
                                         cwd=vis_path)
    except:
        pass

    plt.figure(2)
    plt.clf()
    fig, axes = plt.subplots(num=2, sharex=True)
    axes.semilogy(time, np.abs(ak_rho), label=r'$\rho_k$')
    axes.semilogy(time, np.abs(ak_v), label=r'$v_k$')
    axes.semilogy(time, np.abs(ak_B), label=r'$B_k$')
    axes.semilogy(time, A*np.exp(-gam*np.array(time)), 'k--', label='Theory')
    axes.legend()
    axes.set_title('fits {:1.4f}, {:1.4f}, {:1.4f}, theory {:1.4f}'.format(
                   gam_rho, gam_v, gam_B, gam))
    fig.savefig(vis_path + "/decay_nu_{}_sign{}.pdf".format(nu, sign))

# criteria for failing the test
if failed:
    sys.exit(1)

print("normal exit")
