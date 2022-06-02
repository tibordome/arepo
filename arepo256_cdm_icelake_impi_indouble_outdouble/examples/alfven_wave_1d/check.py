#!/usr/bin/env python3
""" @package ./examples/alfven_wave_1d/check.py
Code that checks results of 1d Alfven wave propagation problem

created by Alessandro Stenghel and Federico Marinacci, 
last modified 13.7.2020 -- comments welcome
"""

# load libraries
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os     # file specific calls
import os.path

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

# initial state -- copied from create.py
density_0 = FloatType(1.0)
velocity_0 = FloatType(0.0)
pressure_0 = FloatType(1.0)
gamma = FloatType(5.0) / FloatType(3.0)
gamma_minus_one = gamma - FloatType(1.0)
delta = FloatType(1e-6)    # relative velocity perturbation
uthermal_0 = pressure_0 / density_0 / gamma_minus_one
bfield_0= FloatType(1.0)
k_z = FloatType(2.0*np.pi)
omega = bfield_0*k_z/np.sqrt(density_0) 


def verify_result(path):
    # open initial conditions to get parameters
    try:
        with h5py.File(os.path.join(path, 'ics.hdf5'), 'r') as data:
            Boxsize = FloatType(data["Header"].attrs["BoxSize"])
            NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0])
    except (OSError, IOError):
        return False, ['Could not open initial conditions!']
    # maximum L1 error after two propagations
    DeltaMaxAllowed = 1e-4 * FloatType(NumberOfCells)**-2
    #    loop over all output files; need to be at times when analytic
    #    solution equals the initial conditions
    i_file = 0
    status = True
    info = []
    error_data = []
    while True:
        # get simulation data
        directory = os.path.join(path, 'output')
        filename = 'snap_%03d.hdf5' % i_file
        try:
            # try to read in snapshot
            time, Pos, Density, Mass, Velocity, Uthermal, Bfield, Volume, Pressure = \
                read_sim_data(os.path.join(directory, filename))
        except (OSError, IOError):
            # should have at least 9 snapshots
            if i_file <= 8:
                status = False
                info.append('Could not find snapshot ' + filename + '!')
            break

        # calculate analytic solution at cell positions
        Density_ref, Velocity_ref, Bfield_ref, Pressure_ref, L1_dens, L1_vel_y, L1_vel_z, L1_bfield_y, L1_bfield_z, L1_pressure = \
            calc_ref_data(time, Pos, Density, Velocity, Bfield, Volume, Pressure)

        # printing results
        info.append("alfven_wave_1d: L1 error of " + filename +":")
        info.append("\t density: %g" % L1_dens)
        info.append("\t velocity y: %g" % L1_vel_y)
        info.append("\t velocity z: %g" % L1_vel_z)
        info.append("\t magnetic field y: %g" % L1_bfield_y)
        info.append("\t magnetic field z: %g" % L1_bfield_z)
        info.append("\t pressure: %g" % L1_pressure)
        info.append("\t tolerance: %g for %d cells" % (DeltaMaxAllowed, NumberOfCells))

        error_data.append(np.array([time, L1_dens, L1_vel_y, L1_vel_z,
                          L1_bfield_y, L1_bfield_z, L1_pressure], dtype=FloatType))

        # criteria for failing the test
        success = (
            L1_dens <= DeltaMaxAllowed and
            L1_vel_y <= DeltaMaxAllowed and
            L1_vel_z <= DeltaMaxAllowed and
            L1_bfield_y <= DeltaMaxAllowed and
            L1_bfield_z <= DeltaMaxAllowed and
            L1_pressure <= DeltaMaxAllowed
        )
        status = status and success

        i_file += 1

    #save L1 errors
    np.savetxt(os.path.join(path, 'error_%d.txt' % NumberOfCells), np.array(error_data))
    return status, info


def visualize_result(path, Lx, Ly):
    if not os.path.exists(os.path.join(path, 'plots')):
        os.mkdir(os.path.join(path, 'plots'))
    
    # open initial conditions to get parameters
    with h5py.File(os.path.join(path, 'ics.hdf5'), 'r') as data:
        Boxsize = FloatType(data["Header"].attrs["BoxSize"])
        NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0])

    # only import matplotlib if needed
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rc_file_defaults()
    plt.rcParams['text.usetex'] = True

    i_file = 1
    while True:
        # try to read in snapshot
        directory = os.path.join(path, 'output')
        filename = 'snap_%03d.hdf5' % i_file
        try:
            time, Pos, Density, Mass, Velocity, Uthermal, Bfield, Volume, Pressure = \
                read_sim_data(os.path.join(directory, filename))
        except (OSError, IOError):
            break
        # calculate analytic solution at cell positions
        Density_ref, Velocity_ref, Bfield_ref, Pressure_ref, L1_dens, L1_vel_y, L1_vel_z, L1_bfield_y, L1_bfield_z, L1_pressure = \
            calc_ref_data(time, Pos, Density, Velocity, Bfield, Volume, Pressure)
        # plot density
        f = plt.figure( figsize=(3.5,3.5) )
        ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
        dx = Boxsize / FloatType(Pos.shape[0])
        ax.plot( Pos[:,0], Density_ref , 'k', lw=0.7, label="Analytical solution" )
        ax.plot( Pos[:,0], Density, 'o-r', mec='r', mfc="None", label="Arepo" )
        ax.set_xlim( 0, Boxsize )
        ax.set_xlabel( "x" )
        ax.set_ylabel( "Density" )
        ax.legend( loc='upper right', frameon=False, fontsize=8 )
        ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_dens), loc='right', size=8 )
        plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        f.savefig(os.path.join(path, 'plots', 'density_%02d.pdf' % i_file))
        plt.close(f)

        # plot pressure
        f = plt.figure( figsize=(3.5,3.5) )
        ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
        ax.plot( Pos[:,0], Pressure_ref , 'k', lw=0.7, label="Analytical solution" )
        ax.plot( Pos[:,0], Pressure, 'o-r', mec='r', mfc="None", label="Arepo" )
        ax.set_xlim( 0, Boxsize )
        ax.set_xlabel( "x" )
        ax.set_ylabel( "Pressure" )
        ax.legend( loc='upper right', frameon=False, fontsize=8 )
        ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_pressure), loc='right', size=8 )
        plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        f.savefig(os.path.join(path, 'plots', 'pressure_%02d.pdf' % i_file))
        plt.close(f)

        # plot velocities
        f = plt.figure( figsize=(3.5,3.5) )
        ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
        ax.plot( Pos[:,0], Velocity_ref[:,1] , ':k', lw=0.7, label="Analytical solution y" )
        ax.plot( Pos[:,0], Velocity[:,1] , 'o-m', mec='m', mfc="None", label="Arepo v y" )
        ax.set_xlim( 0, Boxsize )
        ax.set_xlabel( "x" )
        ax.set_ylabel( "Velocity y" )
        ax.legend( loc='upper right', frameon=False, fontsize=8 )
        ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_vel_y), loc='right', size=8 )
        plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        f.savefig(os.path.join(path, 'plots', 'velocityy_%02d.pdf' % i_file))
        plt.close(f)

        f = plt.figure( figsize=(3.5,3.5) )
        ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
        ax.plot( Pos[:,0], Velocity_ref[:,2] , ':k', lw=0.7, label="Analytical solution z" )
        ax.plot( Pos[:,0], Velocity[:,2], 'o-c', mec='c', mfc="None", label="Arepo v z" )
        ax.set_xlim( 0, Boxsize )
        ax.set_xlabel( "x" )
        ax.set_ylabel( "Velocity z" )
        ax.legend( loc='upper right', frameon=False, fontsize=8 )
        ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_vel_z), loc='right', size=8 )
        plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        f.savefig(os.path.join(path, 'plots', 'velocityz_%02d.pdf' % i_file))
        plt.close(f)

        #plot Bfields
        f = plt.figure( figsize=(3.5,3.5) )
        ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
        ax.plot( Pos[:,0], Bfield_ref[:,1] , ':k', lw=0.7, label="Analytical solution y" )
        ax.plot( Pos[:,0], Bfield[:,1], 'o-m', mec='m', mfc="None", label="Arepo B y" )
        ax.set_xlim( 0, Boxsize )
        ax.set_xlabel( "x" )
        ax.set_ylabel( "Magnetic Field y" )
        ax.legend( loc='upper right', frameon=False, fontsize=8 )
        ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_bfield_y), loc='right', size=8 )
        plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        f.savefig(os.path.join(path, 'plots', 'bfieldy_%02d.pdf' % i_file))
        plt.close(f)

        f = plt.figure( figsize=(3.5,3.5) )
        ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
        ax.plot( Pos[:,0], Bfield_ref[:,2] , ':k', lw=0.7, label="Analytical solution z" )
        ax.plot( Pos[:,0], Bfield[:,2], 'o-c', mec='c', mfc="None", label="Arepo B z" )
        ax.set_xlim( 0, Boxsize )
        ax.set_xlabel( "x" )
        ax.set_ylabel( "Magnetic Field z" )
        ax.legend( loc='upper right', frameon=False, fontsize=8 )
        ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_bfield_z), loc='right', size=8 )
        plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
        f.savefig(os.path.join(path, 'plots', 'bfieldz_%02d.pdf' % i_file))
        plt.close(f)

        i_file += 1


def read_sim_data(path):
    data = h5py.File(path, 'r')
    # get simulation data
    ## simulation data
    time = FloatType(data['Header'].attrs['Time'])
    Pos = np.array(data["PartType0"]["CenterOfMass"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    Bfield = np.array(data["PartType0"]["MagneticField"], dtype = FloatType)/FloatType(np.sqrt(4.0*np.pi))
    data.close()
    Volume = Mass / Density
    Pressure = gamma_minus_one*Density*Uthermal
    return time, Pos, Density, Mass, Velocity, Uthermal, Bfield, Volume, Pressure


def calc_ref_data(time, Pos, Density, Velocity, Bfield, Volume, Pressure):
    # calculate analytic solution at cell positions
    Density_ref = np.full(Pos.shape[0], density_0, dtype=FloatType)
    Velocity_ref = np.zeros(Pos.shape, dtype=FloatType)
    Pressure_ref = np.full(Pos.shape[0], pressure_0, dtype=FloatType)
    Bfield_ref = np.zeros(Pos.shape, dtype=FloatType)

    ## perturbations
    Velocity_ref[:,0] = velocity_0
    Velocity_ref[:,1] = delta*np.sin(k_z*Pos[:,0]-omega*time)
    Velocity_ref[:,2] = delta*np.cos(k_z*Pos[:,0]-omega*time)
    Bfield_ref[:,0] = bfield_0
    Bfield_ref[:,1] = -k_z*bfield_0/omega*Velocity_ref[:,1]
    Bfield_ref[:,2] = -k_z*bfield_0/omega*Velocity_ref[:,2]

    # compare data
    ## density
    abs_delta_dens = np.abs(Density - Density_ref)
    L1_dens = np.average(abs_delta_dens, weights=Volume)

    ## velocity 
    abs_delta_vel_y = np.abs(Velocity - Velocity_ref)[:,1]
    L1_vel_y = np.average(abs_delta_vel_y, weights=Volume)

    abs_delta_vel_z = np.abs(Velocity - Velocity_ref)[:,2]
    L1_vel_z = np.average(abs_delta_vel_z, weights=Volume)

    ## magnetic field 
    abs_delta_bfield_y = np.abs(Bfield -Bfield_ref)[:,1]
    L1_bfield_y = np.average(abs_delta_bfield_y, weights=Volume)

    abs_delta_bfield_z = np.abs(Bfield -Bfield_ref)[:,2]
    L1_bfield_z = np.average(abs_delta_bfield_z, weights=Volume)

    ## pressure
    abs_delta_pressure = np.abs(Pressure-Pressure_ref)
    L1_pressure = np.average(abs_delta_pressure, weights=Volume)

    return Density_ref, Velocity_ref, Bfield_ref, Pressure_ref, L1_dens, L1_vel_y, L1_vel_z, L1_bfield_y, L1_bfield_z, L1_pressure


if __name__ == '__main__':
    simulation_directory = str(sys.argv[1])
    print("alfven_wave_1d: checking simulation output in directory " + simulation_directory)
    makeplots = True
    if len(sys.argv) > 2:
        if sys.argv[2] == "True":
            makeplots = True
        else:
            makeplots = False
    # perform checks
    status_ok, info = verify_result(simulation_directory)
    for msg in info:
        print(msg)
    # optionally plot
    if makeplots:
        visualize_result(simulation_directory, None, None)
    sys.exit(int(not status_ok))
