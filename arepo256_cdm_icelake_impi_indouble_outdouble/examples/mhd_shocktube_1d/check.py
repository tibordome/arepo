#!/usr/bin/env python3
""" @package ./examples/mhd_shocktube_1d/check.py
Code that checks results of 1d mhd shocktube problem

created by Rainer Weinberger, last modified 12.03.2019
"""
# load libraries
import sys  ## load sys; needed for exit codes
import numpy as np  ## load numpy
import h5py  ## load h5py; needed to read snapshots
import os  # file specific calls
import os.path
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc_file_defaults()

IC_FILENAME = 'ics.hdf5'
CreateReferenceSolution = False

FloatType = np.float64  # double precision: np.float64, for single use np.float32


def verify_result(path):
    ## open initial conditions to get parameters
    try:
        data = h5py.File(os.path.join(path, IC_FILENAME), 'r')
    except (OSError, IOError):
        return False, ['Could not open initial conditions!']
    Boxsize = FloatType(data['Header'].attrs['BoxSize'])
    NumberOfCells = np.int32(data['Header'].attrs['NumPart_Total'][0])
    CellsPerDimension = np.sqrt(NumberOfCells)  ## 2d sim
    data.close()
    # loop over all output files
    status = True
    i_file = 1
    info = []
    while True:
        # get simulation data
        directory = os.path.join(path, 'output')
        filename = 'snap_%03d.hdf5' % i_file
        try:
            x, vel, rho, u, B = read_sim_data(os.path.join(
                directory, filename))
        except (OSError, IOError):
            # should have at least three snapshots
            if i_file <= 2:
                info.append('Could not find snapshot ' + filename + '!')
                return False, info
            break

        absB = np.sqrt(B[:, 0]**2 + B[:, 1]**2 + B[:, 2]**2)
        alphaB = (B[:, 1] / B[:, 2])

        verificationData = np.array([x[:, 0], rho, vel[:, 0], u, absB]).T

        if CreateReferenceSolution:
            checkData = verificationData[::40, :]
            np.savetxt(os.path.join(path, 'data_%03d.txt' % i_file), checkData)
        else:
            checkData, rho_ref, vel_ref, u_ref, absB_ref = read_check_data(
                os.path.join(path, 'data_%03d.txt' % i_file), x)

            delta_rho = rho - rho_ref
            delta_vel = vel[:, 0] - vel_ref
            delta_u = u - u_ref
            delta_B = absB - absB_ref

            res_scaling = 200. / np.float64(len(rho))
            tolerance_rho = 0.02 * res_scaling
            tolerance_vel = 0.04 * res_scaling
            tolerance_u = 0.09 * res_scaling
            tolerance_B = 0.05 * res_scaling

            success = (
                np.std(delta_rho) <= tolerance_rho and
                np.std(delta_vel) <= tolerance_vel and
                np.std(delta_u) <= tolerance_u and
                np.std(delta_B) <= tolerance_B
            )
            status = status and success

            info.append(
                'standard deviations of absolute error and tolerance (density, velocity, int. energy, magnetic field:'
            )
            info.append('{} {}'.format(np.std(delta_rho), tolerance_rho))
            info.append('{} {}'.format(np.std(delta_vel), tolerance_vel))
            info.append('{} {}'.format(np.std(delta_u), tolerance_u))
            info.append('{} {}'.format(np.std(delta_B), tolerance_B))

        i_file += 1
    return status, info


def visualize_result(path, Lx, Ly):
    # loop over all output files
    status = True
    i_file = 1
    info = []
    while True:
        # get simulation data
        directory = os.path.join(path, 'output')
        filename = 'snap_%03d.hdf5' % i_file
        try:
            # try to read in snapshot
            x, vel, rho, u, B = read_sim_data(os.path.join(
                directory, filename))
        except (OSError, IOError):
            break

        absB = np.sqrt(B[:, 0]**2 + B[:, 1]**2 + B[:, 2]**2)
        alphaB = (B[:, 1] / B[:, 2])

        verificationData = np.array([x[:, 0], rho, vel[:, 0], u, absB]).T

        checkData, rho_ref, vel_ref, u_ref, absB_ref = read_check_data(
            os.path.join(path, 'data_%03d.txt' % i_file), x)

        fig, ax = plt.subplots(4, sharex=True, figsize=np.array([6.9, 6.0]))
        fig.subplots_adjust(left=0.09, bottom=0.09, right=0.98, top=0.98)

        ax[0].plot(x[:, 0], rho, 'b', zorder=3, label='Arepo')
        ax[0].plot(x[:, 0], rho, 'r+', zorder=2, label='Arepo cells')
        ax[0].plot(checkData[:, 0],
                   checkData[:, 1],
                   c='k',
                   zorder=1,
                   lw=0.7,
                   label='Exact solution')
        ax[0].set_ylabel(r'density')

        ax[1].plot(x[:, 0], vel[:, 0], 'b', zorder=3)
        ax[1].plot(x[:, 0], vel[:, 0], 'r+', zorder=2)
        ax[1].plot(checkData[:, 0], checkData[:, 2], c='k', zorder=1, lw=0.7)
        ax[1].set_ylabel(r'vel')

        ax[2].plot(x[:, 0], u, 'b', zorder=3)
        ax[2].plot(x[:, 0], u, 'r+', zorder=2)
        ax[2].plot(checkData[:, 0], checkData[:, 3], c='k', zorder=1, lw=0.7)
        ax[2].set_ylabel(r'u$_{th}$')

        ax[3].plot(x[:, 0], absB, 'b', zorder=3)
        ax[3].plot(x[:, 0], absB, 'r+', zorder=2)
        ax[3].plot(checkData[:, 0], checkData[:, 4], c='k', zorder=1, lw=0.7)
        ax[3].set_ylabel(r'abs(B)')
        ax[3].set_xlabel(r'position')
        ax[3].set_xlim([-0.1, 2.6])

        fig.align_ylabels(ax[:])
        ax[0].legend(loc='upper right', frameon=False, fontsize=8)

        if not os.path.exists(os.path.join(path, 'plots')):
            os.mkdir(os.path.join(path, 'plots'))

        fig.savefig(os.path.join(path, 'plots', 'figure_%03d.pdf' % i_file))

        i_file += 1


def read_sim_data(path):
    data = h5py.File(path, 'r')
    # get simulation data
    x = np.array(data['PartType0']['Coordinates'], dtype=FloatType)
    vel = np.array(data['PartType0']['Velocities'], dtype=FloatType)
    rho = np.array(data['PartType0']['Density'], dtype=FloatType)
    u = np.array(data['PartType0']['InternalEnergy'], dtype=FloatType)
    B = np.array(data['PartType0']['MagneticField'], dtype=FloatType)
    data.close()
    return x, vel, rho, u, B


def read_check_data(path, x):
    checkData = np.loadtxt(path)

    rho_ref = np.interp(x[:, 0], checkData[:, 0], checkData[:, 1])
    vel_ref = np.interp(x[:, 0], checkData[:, 0], checkData[:, 2])
    u_ref = np.interp(x[:, 0], checkData[:, 0], checkData[:, 3])
    absB_ref = np.interp(x[:, 0], checkData[:, 0], checkData[:, 4])

    return checkData, rho_ref, vel_ref, u_ref, absB_ref


if __name__ == '__main__':
    makeplots = False
    if len(sys.argv) > 2:
        makeplots = sys.argv[2] == 'True'

    simulation_directory = str(sys.argv[1])
    print('mhd_shocktube_1d: checking simulation output in directory ' +
          simulation_directory)

    # perform checks
    status_ok, info = verify_result(simulation_directory)
    for msg in info:
        print(msg)
    # optionally plot
    if makeplots:
        visualize_result(simulation_directory, None, None)

    sys.exit(int(not status_ok))
