#!/usr/bin/env python3
"""
Code that checks results of 1d shocktube problem
"""
# load libraries
import sys  # needed for exit codes
import numpy as np  # scientific computing package
import h5py  # hdf5 format
import os  # file specific calls
import os.path
import matplotlib
import matplotlib.pyplot as plt  ## needs to be active for plotting!
matplotlib.rc_file_defaults()

# use relative import for local module if not executed directly in order to
# prevent importing several modules of the same name
if __name__ == '__main__':
    # make sure that this file's directory is in sys.path
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from Riemann import *  ## Riemann solver
    import check_parameters
else:
    from .Riemann import *  ## Riemann solver
    from . import check_parameters

IC_FILENAME = 'ics.hdf5'
gamma = 1.4  ## needs to be identical to Config.sh (= 5/3 if not specified there)!
Dtype = np.float64  # double precision: np.float64, for single use np.float32
## open initial conditions to get parameters
forceExitOnError = False  ## exits immediately when tolerance is exceeded

# check functions


def verify_result(path):
    try:
        W_L, W_R, position_0 = read_ics(path)
    except (OSError, IOError):
        return False, [
            'Could not find file {}'.format(os.path.join(path, IC_FILENAME))
        ]

    # loop over all output files
    ReturnFlag = 0
    info = []
    i_file = -1
    while True:
        i_file += 1

        ## read in data
        directory = os.path.join(path, 'output')
        filename = 'snap_%03d.hdf5' % i_file
        info.append('Trying to open ' + os.path.join(directory, filename))
        try:
            W, time, position = read_data(os.path.join(directory, filename))
        except (OSError, IOError):
            # should have at least 6 snapshots
            if i_file <= 5:
                info.append('Could not find snapshot ' + filename + '!')
                return False, info
            break
        info.append('Analyzing ' + filename)

        # perform checks
        flag, msg = CheckL1Error(position[:, 0], W, W_L, W_R, gamma,
                                 position_0, time)
        ReturnFlag += flag
        info += msg
        if ReturnFlag > 0 and forceExitOnError:
            info.append('ERROR: exceeding tolerance!')
            return False, info

        flag, msg = CheckTotalVariation(position[:, 0], W, W_L, W_R, gamma,
                                        position_0, time)
        ReturnFlag += flag
        info += msg
        if ReturnFlag > 0 and forceExitOnError:
            info.append('ERROR: exceeding tolerance!')
            return False, info

        flag, msg = CheckWidthOfDiscontinuities(position[:, 0], W, W_L, W_R,
                                                gamma, position_0, time)
        ReturnFlag += flag
        info += msg
        if ReturnFlag > 0 and forceExitOnError:
            info.append('ERROR: exceeding tolerance!')
            return False, info

    return ReturnFlag == 0, info


def visualize_result(path, Lx, Ly):
    try:
        W_L, W_R, position_0 = read_ics(path)
    except (OSError, IOError):
        return 1
    # loop over all output files
    ReturnFlag = 0
    i_file = -1
    while True:
        i_file += 1
        ## read in data
        directory = os.path.join(path, 'output')
        filename = 'snap_%03d.hdf5' % i_file
        try:
            W, time, position = read_data(os.path.join(directory, filename))
        except (OSError, IOError):
            break
        ReturnFlag += PlotSimulationData(position[:, 0], W, W_L, W_R, gamma,
                                         position_0, time, path, i_file)
    return ReturnFlag


def read_ics(path):
    data = h5py.File(os.path.join(path, IC_FILENAME), 'r')
    IC_position = np.array(data['PartType0']['Coordinates'], dtype=np.float64)
    IC_mass = np.array(data['PartType0']['Masses'], dtype=np.float64)
    IC_velocity = np.array(data['PartType0']['Velocities'], dtype=np.float64)
    IC_internalEnergy = np.array(data['PartType0']['InternalEnergy'],
                                 dtype=np.float64)
    NumberOfCells = np.int32(data['Header'].attrs['NumPart_Total'][0])
    DeltaMaxAllowed = 2.0 / np.float(NumberOfCells)

    ## get parameters of Riemann problem from ICs
    dx = IC_position[1, 0] - IC_position[0, 0]  ## assuming a uniform grid
    W_L = np.array([
        IC_mass[0] / dx, IC_velocity[0, 0],
        (gamma - 1.0) * IC_internalEnergy[0] * IC_mass[0] / dx
    ])
    W_R = np.array([
        IC_mass[-1] / dx, IC_velocity[-1, 0],
        (gamma - 1.0) * IC_internalEnergy[-1] * IC_mass[-1] / dx
    ])
    i_0, = np.where((IC_mass[:] / dx == W_R[0]) & (IC_velocity[:, 0] == W_R[1])
                    & ((gamma - 1.0) * IC_internalEnergy[-1] * IC_mass[-1] /
                       dx == W_R[2]))
    position_0 = 0.5 * (IC_position[i_0[0] - 1, 0] + IC_position[i_0[0], 0]
                        )  ## discontinuity at interface position

    return W_L, W_R, position_0


def read_data(path):
    ## open hdf5 file
    data = h5py.File(path, 'r')
    ## get data from snapshot
    time = np.float(data['Header'].attrs['Time'])
    position = np.array(data['PartType0']['Coordinates'], dtype=np.float64)
    density = np.array(data['PartType0']['Density'], dtype=np.float64)
    vel = np.array(data['PartType0']['Velocities'], dtype=np.float64)
    internalEnergy = np.array(data['PartType0']['InternalEnergy'],
                              dtype=np.float64)
    ## convert to more useful data structure
    W = np.array(
        [density, vel[:, 0], (gamma - 1.0) * internalEnergy * density],
        dtype=np.float64).T  ## shape: (n, 3)

    return W, time, position


def CheckL1Error(Pos, W, W_L, W_R, gamma, position_0, time):
    r"""
    Compare the hydrodynamical quantities of the simulation with the exact
    solution of the Riemann problem, calculate the L1 error and check whether
    avarge L1 error is acceptable

    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data

    \return status: zero if everything is within tolerances, otherwise 1
    """
    xx, W_exact, PosOfCharacteristics = RiemannProblem(Pos, position_0, W_L,
                                                       W_R, gamma, time)

    norm = np.abs(W_exact)
    norm[:, 1] = 1  ## use absolute error in velocity component

    ## calculate L1 norm
    delta = np.abs(W_exact - W) / norm
    L1 = np.average(delta, axis=0)

    ## tolerance value; found empirically, fist order convergence!
    L1MaxAllowed = 2.0 / np.float(Pos.shape[0])

    if np.all(L1 <= L1MaxAllowed):
        return 0, [
            'CheckL1Error: L1 error fine: %g %g %g; tolerance %g' %
            (L1[0], L1[1], L1[2], L1MaxAllowed)
        ]
    else:
        return 1, [
            'CheckL1Error: ERROR: L1 error too large: %g %g %g; tolerance %g' %
            (L1[0], L1[1], L1[2], L1MaxAllowed)
        ]


def CheckTotalVariation(Pos, W, W_L, W_R, gamma, position_0, time):
    r"""
    Compare the total variation in simulation quantities with the total
    variation in the analytic solution of the Riemann problem

    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data

    \return status: zero if everything is within tolerances, otherwise 1
    """
    xx, W_exact, PosOfCharacteristics = RiemannProblem(Pos, position_0, W_L,
                                                       W_R, gamma, time)

    TotalVariationSim = np.zeros(3)
    TotalVariationExact = np.zeros(3)
    ratio = np.ones(3)

    i_sorted = np.argsort(Pos)  ## sorted by position
    dW = W[i_sorted[1:], :] - W[
        i_sorted[:-1], :]  ## difference of neighbouring cells
    dW_exact = W_exact[i_sorted[1:], :] - W_exact[i_sorted[:-1], :]

    for i in np.arange(3):
        i1_pos, = np.where(dW[:, i] >= 0)
        i1_neg, = np.where(dW[:, i] < 0)
        TotalVariationSim[i] = np.sum(dW[i1_pos, i], axis=0) - np.sum(
            dW[i1_neg, i])
        TotalVariationExact[i] = np.sum(dW_exact[i1_pos, i], axis=0) - np.sum(
            dW_exact[i1_neg, i])
        # handle case where TotalVariationExact is zero, only allowing it when
        # TotalVariationSim is also zero
        if TotalVariationExact[i] == 0:
            if TotalVariationSim[i] == TotalVariationExact[i]:
                ratio[i] = 1
            else:
                ratio[i] = np.inf
        else:
            ratio[i] = TotalVariationSim[i] / TotalVariationExact[i]

    MaxRatioAllowed = check_parameters.MaxRatioAllowed
    if np.all(ratio <= MaxRatioAllowed):
        return 0, [
            'CheckTotalVariation: TotalVariation Sim/Exact fine: %g %g %g, tolerance: %g'
            % (ratio[0], ratio[1], ratio[2], MaxRatioAllowed)
        ]
    else:
        return 1, [
            'CheckTotalVariation: ERROR: TotalVariation Sim/Exact: %g %g %g, tolerance: %g'
            % (ratio[0], ratio[1], ratio[2], MaxRatioAllowed)
        ]


def CheckWidthOfDiscontinuities(Pos, W, W_L, W_R, gamma, position_0, time):
    r"""
    Measure the width of the fluid discontinuities in simulation quantities
    to assess the numerical diffusivity of the scheme

    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data

    \return status: zero if everything is within tolerances, otherwise 1
    """
    xx, W_exact, PosOfCharacteristics = RiemannProblem(Pos, position_0, W_L,
                                                       W_R, gamma, time)

    for i, pos_char in enumerate(PosOfCharacteristics):
        ReturnFlag = 0
        ## PosOfCharacteristics will give different values in index 0 and 1 (3 and 4)
        ## if it is a rarefaction; in this case, don't check, otherwise do check
        if i == 1 or i == 4:
            continue
        if i == 0 or i == 3:
            if PosOfCharacteristics[i] != PosOfCharacteristics[i + 1]:
                continue

        ## get hydrodynamic states left and right of jump
        xx = np.array([10.0 + pos_char - 1e-3, 10.0 + pos_char + 1e-3])
        xx, W_exact, dummy = RiemannProblem(xx, 10.0, W_L, W_R, gamma, time)

        ## get threshold values to measure how many cells are within this interval
        jump = W_exact[1, :] - W_exact[0, :]
        percentiles = np.array(
            [W_exact[0, :] + 0.05 * jump, W_exact[0, :] + 0.95 * jump])
        percentile_05 = np.min(percentiles, axis=0)
        percentile_95 = np.max(percentiles, axis=0)
        i_sorted = np.argsort(np.abs(Pos - pos_char - position_0))

        i_low = np.full(3, i_sorted[0], dtype=np.int32)
        i_high = np.full(3, i_sorted[0], dtype=np.int32)

        ## check by how many cells 5th to 95th percentile of jump are sampled
        for j in np.arange(3):
            while W[i_low[j], j] > percentile_05[j] and W[
                    i_low[j], j] < percentile_95[j]:
                i_low[j] -= 1
            while W[i_low[j], j] > percentile_05[j] and W[
                    i_low[j], j] < percentile_95[j]:
                i_high[j] += 1

        ## sufficient for exact Riemann solver
        MaxNumerOfCells = 4

        info = []
        if i == 2:
            info.append(
                'CheckWidthOfDiscontinuities: density jump at contact discontinuity resolved by %d cells (5th to 95th precentile), tolerance: %d'
                % (i_high[0] - i_low[0], MaxNumerOfCells))
        else:
            info.append(
                'CheckWidthOfDiscontinuities: density, velocity and pressure jump at shock resolved by %d, %d and %d cells (5th to 95th precentile), tolerance: %d'
                % (i_high[0] - i_low[0], i_high[1] - i_low[1],
                   i_high[2] - i_low[2], MaxNumerOfCells))
        if np.any(i_high - i_low > MaxNumerOfCells):
            info.append(
                'CheckWidthOfDiscontinuities: ERROR: discontinuity too wide!')
            ReturnFlag += 1

    return ReturnFlag, info


def PlotSimulationData(Pos, W, W_L, W_R, gamma, position_0, time,
                       simulation_directory, i_file):
    r"""
    Plot density, velocity, specific internal energy and pressure of
    Simulation and exact solution on top

    \param[in] Pos: shape (n,) array with 1d positions
    \param[in] W: shape (n,3) array with density, velocity and pressure of each cell
    \param[in] W_L: shape (3,) array with left initial condition state (density, velocity and pressure)
    \param[in] W_R: shape (3,) array with right initial condition state (density, velocity and pressure)
    \param[in] gamma: adiabatic index of gas
    \param[in] position_0: initial position of discontinuity
    \param[in] time: time of snapshot data
    \param[in] simulation_directory: path to simulation

    \return status: zero if everything went fine
    """
    xx = np.linspace(Pos.min(), Pos.max(), 1000)
    xx, W_exact, dummy = RiemannProblem(xx, position_0, W_L, W_R, gamma, time)

    min = np.zeros(4)
    max = np.zeros(4)
    min[0] = np.min(W_exact[:, 0]) - 0.1
    max[0] = np.max(W_exact[:, 0]) + 0.1
    min[1] = np.min(W_exact[:, 1]) - 0.1
    max[1] = np.max(W_exact[:, 1]) + 0.1
    min[2] = np.min(W_exact[:, 2] / W_exact[:, 0] / (gamma - 1.0)) - 0.1
    max[2] = np.max(W_exact[:, 2] / W_exact[:, 0] / (gamma - 1.0)) + 0.1
    min[3] = np.min(W_exact[:, 2]) - 0.1
    max[3] = np.max(W_exact[:, 2]) + 0.1

    ## plot
    fig, ax = plt.subplots(4, sharex=True, figsize=np.array([6.9, 6.0]))
    fig.subplots_adjust(left=0.13, bottom=0.09, right=0.98, top=0.98)

    ax[0].plot(xx, W_exact[:, 0], color='k', lw=0.7, label='Exact solution')
    ax[1].plot(xx, W_exact[:, 1], color='k', lw=0.7)
    ax[2].plot(xx,
               W_exact[:, 2] / W_exact[:, 0] / (gamma - 1.0),
               color='k',
               lw=0.7)
    ax[3].plot(xx, W_exact[:, 2], color='k', lw=0.7)

    ax[0].plot(Pos, W[:, 0], '+r', mec='r', mfc='None', label='Arepo cells')
    ax[1].plot(Pos, W[:, 1], '+r', mec='r', mfc='None')
    ax[2].plot(Pos,
               W[:, 2] / W[:, 0] / (gamma - 1.0),
               '+r',
               mec='r',
               mfc='None')
    ax[3].plot(Pos, W[:, 2], '+r', mec='r', mfc='None')

    ax[0].plot(Pos, W[:, 0], 'b', label='Arepo')
    ax[1].plot(Pos, W[:, 1], 'b')
    ax[2].plot(Pos, W[:, 2] / W[:, 0] / (gamma - 1.0), 'b')
    ax[3].plot(Pos, W[:, 2], 'b')

    ax[3].set_xlim([0, 20])
    for i_plot in np.arange(4):
        ax[i_plot].set_ylim([min[i_plot], max[i_plot]])

    ax[0].legend(loc='upper right', frameon=False, fontsize=8)

    ## set labels
    ax[3].set_xlabel(r'pos')
    ax[0].set_ylabel(r'density')
    ax[1].set_ylabel(r'velocity')
    ax[2].set_ylabel(r'spec. int. energy')
    ax[3].set_ylabel(r'pressure')
    fig.align_ylabels(ax[:])

    if not os.path.exists(os.path.join(simulation_directory, 'plots')):
        os.mkdir(os.path.join(simulation_directory, 'plots'))

    print(
        os.path.join(simulation_directory, 'plots',
                     'figure_%03d.pdf' % i_file))
    fig.savefig(
        os.path.join(simulation_directory, 'plots',
                     'figure_%03d.pdf' % i_file))
    plt.close(fig)

    return 0


if __name__ == '__main__':
    makeplots = False
    if len(sys.argv) > 2:
        makeplots = sys.argv[2] == 'True'

    simulation_directory = str(sys.argv[1])
    print(
        os.path.dirname(__file__) +
        ': checking simulation output in directory ' + simulation_directory)

    ReturnFlag = 0

    # plot data if you want
    if makeplots:
        ReturnFlag += visualize_result(simulation_directory, None, None)
        if ReturnFlag > 0 and forceExitOnError:
            print('ERROR: something went wrong in plot')
            sys.exit(ReturnFlag)
    # perform checks
    status_ok, info = verify_result(simulation_directory)
    ReturnFlag += not status_ok
    for msg in info:
        print(msg)

    if ReturnFlag == 0:
        print('check.py: success!')
    else:
        print('check.py: failed!')

    sys.exit(ReturnFlag)
