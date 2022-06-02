""" @package ./examples/brag_2d/check.py
Code that checks results of 2D anisotropic Braginskii

created by Thomas Berlok, last modified 9.10.2019 -- comments welcome
"""

# load libraries
import sys    # load sys; needed for exit codes
import numpy as np    # load numpy
import h5py    # load h5py; needed to read snapshots
import glob
from numpy import sqrt, sin, cos, exp, pi

# Save a plot (needs matplotlib)
makeplots = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print("examples/brag_2d/check.py: checking simulation output in directory "
      + simulation_directory)

directory = simulation_directory+"/output/"

nu = 1e-2

rho = 1

nx = 256

Lx = 1
x0 = -0.5

u0 = 3/2

dx = Lx/nx

# Grid
xg = (np.arange(nx) + 0.5)*dx + x0

N = 25
n = np.arange(N+1)

kn = n*2*np.pi/Lx
gamma_n = 5/6*nu*kn**2

bx = by = 1/sqrt(2)
nx = bx
ny = -by

cn = np.zeros(N+1)
cn[0] = 2
cn[1:] = -2*sin(3*n[1:]*pi/2)/(n[1:]*pi)*exp(-(n[1:]*pi)**2/400)


def vx_an(x, t):
    def temp(x):
        return -3/10*np.sum(cn*cos(kn*x)*(1 - exp(-gamma_n*t)))
    temp_v = np.vectorize(temp)
    return temp_v(x)


def vy_an(x, t):
    def temp(x):
        return 1/10*np.sum(cn*cos(kn*x)*(1 + 9*exp(-gamma_n*t)))
    temp_v = np.vectorize(temp)
    return temp_v(x)


def u_an(x, t):

    def temp(x):
        heating = 0
        for ii in range(1, N):
            for jj in range(1, N):
                heating += cn[ii]*cn[jj]*sqrt(gamma_n[ii]*gamma_n[jj]) / \
                          (gamma_n[ii] + gamma_n[jj]) * \
                          sin(kn[ii]*x)*sin(kn[jj]*x) * \
                          (1 - exp(-(gamma_n[ii] + gamma_n[jj])*t))
        return heating

    temp_v = np.vectorize(temp)
    return u0 + 9/10*rho*temp_v(x)


def pa_an(x, t):
    def temp(x):
        return -3/2*rho*nu*np.sum(cn*kn*sin(kn*x)*exp(-gamma_n*t))
    temp_v = np.vectorize(temp)
    return temp_v(x)


def load_snap(filename):
    """Helper function for loading snapshots"""
    ptName = 'PartType0'
    fields = {}
    with h5py.File(filename, 'r') as f:
        for k in f[ptName].keys():
            fields[k] = f[ptName][k][()]

        time = f['Header'].attrs['Time']
    return fields, time


def rms_diff(x, y):
    """rms-difference between x and y"""
    return np.sqrt(np.mean((x - y)**2.0))


tol0 = 5e-3
tol1 = 5e-3
tol2 = 5e-3
tol3 = 5e-3
msg0 = 'vx error ({}) too large'
msg1 = 'vy field error ({}) too large'
msg2 = 'InternalEnergy error ({}) too large'
msg3 = 'PresAnis error ({}) too large'

failed = False

fnames = glob.glob(directory + '/snapshot_*.hdf5')
fnames = np.sort(fnames)
for i in range(len(fnames)):
    f, t = load_snap(fnames[i])

    s = f['Coordinates'][:, 0]

    vx = f['Velocities'][:, 0]
    vy = f['Velocities'][:, 1]

    u = f['InternalEnergy']
    pa = f['PressureAnisotropy']

    index = np.argsort(s)

    err0 = rms_diff(vx, vx_an(s-0.5, t))
    err1 = rms_diff(vy, vy_an(s-0.5, t))
    err2 = rms_diff(u, u_an(s-0.5, t))
    err3 = rms_diff(pa, pa_an(s-0.5, t))

    # print(err0, err1, err2, err3)

    if not err0 <= tol0:
        failed = True
        print(msg0.format(err0))
    if not err1 <= tol1:
        failed = True
        print(msg1.format(err1))
    if not err2 <= tol2:
        print(msg2.format(err2))
        failed = True
    if not err3 <= tol3:
        print(msg3.format(err3))
        failed = True


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

    for i in range(len(fnames)):
        plt.figure(1)
        plt.clf()
        fig, axes = plt.subplots(num=1, nrows=2, ncols=2, sharex=True)

        axes[0, 0].set_ylabel(r'$v_x/c$')
        axes[0, 1].set_ylabel(r'$v_y/c$')
        axes[1, 1].set_ylabel(r'$\varepsilon/(\rho c^2)$')
        axes[1, 0].set_ylabel(r'$\Delta p/p$')
        axes[1, 1].set_xlim(-0.5, 0.5)
        axes[1, 1].set_xlim(-0.5, 0.5)
        axes[1, 1].set_xlabel('$x/L$')
        axes[1, 0].set_xlabel('$x/L$')

        line0, = axes[0, 0].plot(xg, vx_an(xg, 0), 'k--', label=r'$t=0$')
        axes[0, 1].plot(xg, vy_an(xg, 0), 'k--')
        axes[1, 1].plot(xg, u_an(xg, 0), 'k--')
        axes[1, 0].plot(xg, pa_an(xg, 0), 'k--')

        line1, = axes[0, 0].plot(xg, vx_an(xg, 1000), 'k-', label=r'$t\rightarrow \infty$')
        axes[0, 1].plot(xg, vy_an(xg, 1000), 'k-')
        axes[1, 1].plot(xg, u_an(xg, 1000), 'k-')
        axes[1, 0].plot(xg, pa_an(xg, 1000), 'k-')

        f, t = load_snap(fnames[i])

        s = f['Coordinates'][:, 0]

        vx = f['Velocities'][:, 0]
        vy = f['Velocities'][:, 1]
        u = f['InternalEnergy']
        pa = f['PressureAnisotropy']

        ind = np.argsort(s)

        line2, = axes[0, 0].plot(s[ind]-0.5, vx[ind], 'C0-', label='Arepo')
        axes[0, 1].plot(s[ind]-0.5, vy[ind], 'C0-', label='Arepo')
        axes[1, 1].plot(s[ind]-0.5, u[ind], 'C0-', label='Arepo')
        axes[1, 0].plot(s[ind]-0.5, pa[ind], 'C0-', label='Arepo')

        line3, = axes[0, 0].plot(xg, vx_an(xg, t), 'C1--', label='Theory')
        axes[0, 1].plot(xg, vy_an(xg, t), 'C1--')
        axes[1, 1].plot(xg, u_an(xg, t), 'C1--')
        axes[1, 0].plot(xg, pa_an(xg, t), 'C1--')

        axes[0, 0].legend(frameon=False)
        axes[0, 0].legend((line2, line3, line0, line1),
                       ('Arepo', 'Theory', r'$t = 0$', r'$t\rightarrow \infty$'),
                       frameon=False, ncol=1)

        axes[0, 0].set_title(r"$t={:1.2f}$".format(t))

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
               'brag_2d.mp4',
               '-y']
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT,
                                         cwd=vis_path)
    except:
        pass


# criteria for failing the test
if failed:
    sys.exit(1)

print("normal exit")
