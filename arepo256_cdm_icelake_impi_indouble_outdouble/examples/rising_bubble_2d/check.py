#!/usr/bin/env python3
""" @package ./examples/rising_bubble_2d/check.py
Code that checks results of 2d bubble problem

created by Thomas Berlok, last modified 30.11.2020 -- comments welcome
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

# Save a plot (needs matplotlib)
makeplots = False
if len(sys.argv) > 2:
    if sys.argv[2] == 'True':
        makeplots = True
    else:
        makeplots = False

simulation_directory = str(sys.argv[1])
print('examples/rising_bubble_2d/check.py: checking simulation output in directory ' +
      simulation_directory)

directory = os.path.join(simulation_directory, 'output')

def load_snap(filename):
    """Helper function for loading snapshots"""
    ptName = 'PartType0'
    fields = {}
    with h5py.File(os.path.join(filename), 'r') as f:
        for k in f[ptName].keys():
            fields[k] = f[ptName][k][()]

        time = f['Header'].attrs['Time']
    return fields, time


import glob
fnames = glob.glob(directory + '/snapshot_*.hdf5')
fnames = np.sort(fnames)

for ii in range(len(fnames)):
    f, t = load_snap(fnames[ii])
    c = f['PassiveScalars']
    index = c>0.
    index
    c = c[index]
    M = f['Masses'][index]
    S = np.sum(-c*np.log(c)*M)
    print('snap {:2d},  S = {:1.3f}'.format(ii, S))


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


    import matplotlib as mpl
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from scipy.spatial import Voronoi, voronoi_plot_2d

    import glob
    fnames = glob.glob(directory + '/snapshot_*.hdf5')
    fnames = np.sort(fnames)

    # for ii in range(len(fnames)):
    for ii in [0, 5, 10]:
        f, t = load_snap(fnames[ii])

        x = f['Coordinates'][:, 0]
        y = f['Coordinates'][:, 1]

        points = np.vstack([x, y]).T

        # generate Voronoi tessellation
        vor = Voronoi(points)

        for key in ['Density', 'PassiveScalars', 'BragViscositySubsteps', 'TimebinHydro']:

            plt.figure(1)
            plt.clf()
            fig, axes = plt.subplots(num=1, sharex=True)

            vmin = np.min(f[key])
            vmax = np.max(f[key])

            # normalize chosen colormap
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
            if key == 'Density':
                mapper = cm.ScalarMappable(norm=norm, cmap=cm.RdBu)
            elif key == 'PassiveScalars':
                vmin = 0
                vmax = 1
                norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
                mapper = cm.ScalarMappable(norm=norm, cmap=cm.RdBu_r)
            else:
                mapper = cm.ScalarMappable(norm=norm, cmap=cm.Accent)

            axes.set_xlim(x.min(), x.max())
            axes.set_ylim(y.min(), y.max())

            # plot Voronoi diagram, and fill finite regions with color
            # voronoi_plot_2d(vor, ax=axes, show_points=False, show_vertices=False, s=1)
            for r in range(len(vor.point_region)):
                region = vor.regions[vor.point_region[r]]
                if not -1 in region:
                    polygon = [vor.vertices[i] for i in region]
                    axes.fill(*zip(*polygon), color=mapper.to_rgba(f[key][r]))
            fig.colorbar(mapper)
            axes.set_title(r"$t={:1.2f}$".format(t))

            fig.savefig(vis_path + "/" + key + "frame_{}.pdf".format(ii))
            print(key + ' image {} done'.format(ii))

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

# criteria for failing the test
if S > 0.075:
    print('Intregrated dye entropy larger than expected for Braginskii')
    print('S = ', S)
    sys.exit(1)

print('normal exit')
