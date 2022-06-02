#!/usr/bin/env python3
""" @package ./examples/polytrope_1d_spherical/check.py
Code that checks results of 1d polytrope test problem

created by Rainer Weinberger, last modified 04.03.2019
"""
# load libraries
import os  # file specific calls
import os.path
import sys  # needed for exit codes
import numpy as np  # scientific computing package
import h5py  # hdf5 format
import matplotlib
import matplotlib.pyplot as plt  ## needs to be active for plotting!
matplotlib.rc_file_defaults()

IC_FILENAME = 'ics.hdf5'
FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32


def verify_result(path):
    ## open initial conditions to get parameters
    try:
        data = h5py.File(os.path.join(path, IC_FILENAME), 'r')
    except (OSError, IOError):
        return False, ['Could not open initial conditions!']

    Boxsize = FloatType(data['Header'].attrs['BoxSize'])
    NumberOfCells = IntType(data['Header'].attrs['NumPart_Total'][0])

    ## maximum L1 error after one propagation; based on experience
    DeltaMaxAllowed = 0.001 * FloatType(NumberOfCells / 256.0)**-2
    # Initial conditions as reference
    i_snap = 0
    directory = os.path.join(path, 'output')
    filename = 'snap_%03d.hdf5' % i_snap
    data = h5py.File(os.path.join(directory, filename), 'r')
    Pos_ref = np.array(data['PartType0']['Coordinates'], dtype=FloatType)[:, 0]
    Density_ref = np.array(data['PartType0']['Density'], dtype=FloatType)
    Velocity_ref = np.array(data['PartType0']['Velocities'],
                            dtype=FloatType)[:, 0]
    Uthermal_ref = np.array(data['PartType0']['InternalEnergy'],
                            dtype=FloatType)

    Accel_ref = np.array(data['PartType0']['Acceleration'],
                         dtype=FloatType)[:, 0]
    GradPress_ref = np.array(data['PartType0']['PressureGradient'],
                             dtype=FloatType)[:, 0]

    msg = []
    while True:
        i_snap += 1
        # compare to snapshots
        filename = 'snap_%03d.hdf5' % i_snap
        try:
            data = h5py.File(os.path.join(directory, filename), 'r')
        except (OSError, IOError):
            # should have at least 9 snapshots
            if i_snap <= 8:
                msg.append('Could not find snapshot ' + filename + '!')
                return False, msg
            break
        Pos = np.array(data['PartType0']['Coordinates'], dtype=FloatType)[:, 0]
        Density = np.array(data['PartType0']['Density'], dtype=FloatType)
        Velocity = np.array(data['PartType0']['Velocities'],
                            dtype=FloatType)[:, 0]
        Uthermal = np.array(data['PartType0']['InternalEnergy'],
                            dtype=FloatType)

        Accel = np.array(data['PartType0']['Acceleration'],
                         dtype=FloatType)[:, 0]
        GradPress = np.array(data['PartType0']['PressureGradient'],
                             dtype=FloatType)[:, 0]

        # compare to ICs by interpolating to IC positions
        delta_dens = np.interp(Pos_ref, Pos, Density) - Density_ref
        delta_vel = np.interp(Pos_ref, Pos, Velocity) - Velocity_ref
        delta_uthermal = np.interp(Pos_ref, Pos, Uthermal) - Uthermal_ref

        ## L1 norm
        i_check, = np.where(Pos_ref < 0.7)

        L1_dens = np.average(np.abs(delta_dens[i_check]))
        L1_vel = np.average(np.abs(delta_vel[i_check]))
        L1_uthermal = np.average(np.abs(delta_uthermal[i_check]))
        # printing results
        msg += [
            'examples/polytrope_spherical_1d/check.py: L1 error of ' +
            filename + ':',
            '\t density: %g' % L1_dens,
            '\t velocity: %g' % L1_vel,
            '\t specific internal energy: %g' % L1_uthermal,
            '\t tolerance: %g for %d cells' % (DeltaMaxAllowed, NumberOfCells)
        ]
        # criteria for failing the test
        success = (
            L1_dens <= DeltaMaxAllowed and
            L1_vel <= DeltaMaxAllowed and
            L1_uthermal <= DeltaMaxAllowed
        )
        if not success:
            return False, msg

    return True, msg


def visualize_result(path, Lx, Ly):
    ## open initial conditions as reference
    data = h5py.File(os.path.join(path, IC_FILENAME), 'r')

    i_snap = 0
    directory = os.path.join(path, 'output')
    filename = 'snap_%03d.hdf5' % i_snap
    data = h5py.File(os.path.join(directory, filename), 'r')
    Pos_ref = np.array(data['PartType0']['Coordinates'], dtype=FloatType)[:, 0]
    Density_ref = np.array(data['PartType0']['Density'], dtype=FloatType)

    Accel_ref = np.array(data['PartType0']['Acceleration'],
                         dtype=FloatType)[:, 0]
    GradPress_ref = np.array(data['PartType0']['PressureGradient'],
                             dtype=FloatType)[:, 0]

    ## plots
    while True:
        i_snap += 1
        # compare to snapshots
        filename = 'snap_%03d.hdf5' % i_snap
        try:
            data = h5py.File(os.path.join(directory, filename), 'r')
        except (OSError, IOError):
            break
        Pos = np.array(data['PartType0']['Coordinates'], dtype=FloatType)[:, 0]
        Density = np.array(data['PartType0']['Density'], dtype=FloatType)
        Velocity = np.array(data['PartType0']['Velocities'],
                            dtype=FloatType)[:, 0]

        Accel = np.array(data['PartType0']['Acceleration'],
                         dtype=FloatType)[:, 0]
        GradPress = np.array(data['PartType0']['PressureGradient'],
                             dtype=FloatType)[:, 0]

        fig, ax = plt.subplots(ncols=1,
                               nrows=4,
                               sharex=True,
                               figsize=np.array([6.9, 6.0]))
        fig.subplots_adjust(left=0.13, bottom=0.09, right=0.98, top=0.98)

        ax[0].plot(Pos, Density, 'b', label='evolved profile')
        ax[0].plot(Pos, Density, 'r+', label='Arepo cells')
        ax[0].plot(Pos_ref,
                   Density_ref,
                   'k--',
                   label='initial profile',
                   lw=0.7)
        ax[0].set_ylabel(r'density')
        ax[0].fill_between([0.1, 0.7], [0.0, 0.0], [1.01, 1.01],
                           color='k',
                           alpha=0.2)
        ax[0].set_ylim(0., 1.)

        ax[1].plot(Pos, Density / Density_ref, 'b')
        ax[1].plot(Pos_ref, [1.0] * Pos_ref.shape[0], 'k--', lw=0.7)
        ax[1].set_ylabel(r'rel. density')
        ax[1].set_ylim([0.99, 1.01])
        ax[1].fill_between([0.1, 0.7], [0.99, 0.99], [1.01, 1.01],
                           color='k',
                           alpha=0.2)

        ax[2].plot(Pos, Velocity, 'b')
        ax[2].plot(Pos_ref, [0.0] * Pos_ref.shape[0], 'k--', lw=0.7)
        ax[2].set_ylabel(r'velocity')
        ax[2].set_ylim([-0.01, 0.01])
        ax[2].fill_between([0.1, 0.7], [-0.01, -0.01], [0.01, 0.01],
                           color='k',
                           alpha=0.2)

        ax[3].plot(Pos[:-1], Accel[:-1] - GradPress[:-1] / Density[:-1], 'b')
        ax[3].plot(Pos_ref[:-1],
                   Accel_ref[:-1] - GradPress_ref[:-1] / Density_ref[:-1],
                   'k--',
                   lw=0.7)
        ax[3].set_ylim([-0.5, 0.1])

        ax[3].set_ylabel(r'net accel.')
        ax[3].set_xlabel(r'radius')
        ax[3].set_xlim([0.0, 0.8])
        ax[3].fill_between([0.1, 0.7], [-0.5, -0.5], [0.1, 0.1],
                           color='k',
                           alpha=0.2)

        ax[0].legend(loc=1, frameon=False)

        fig.align_ylabels(ax[:])

        if not os.path.exists(os.path.join(path, 'plots')):
            os.mkdir(os.path.join(path, 'plots'))
        fig.savefig(os.path.join(path, 'plots', 'profiles_%03d.pdf' % i_snap))


if __name__ == '__main__':
    makeplots = False
    if len(sys.argv) > 2:
        makeplots = sys.argv[2] == 'True'

    simulation_directory = str(sys.argv[1])
    print(
        'examples/polytrope_1d_spherical/check.py: checking simulation output in directory '
        + simulation_directory)

    # perform checks
    status_ok, info = verify_result(simulation_directory)
    for msg in info:
        print(msg)
    # optionally plot
    if makeplots:
        visualize_result(simulation_directory, None, None)

    sys.exit(int(not status_ok))
