#!/usr/bin/env python3
""" @package ./examples/OrszagTang_2d/check.py
Code that checks results for 2d Orszag Tang problem

created by Ruediger Pakmor, last modified 10.12.2021 -- comments welcome
"""
#### load libraries
import sys  ## load sys; needed for exit codes
import numpy as np  ## load numpy
import h5py  ## load h5py; needed to read snapshots
import os  # file specific calls
import os.path


FloatType = np.float64  # double precision: np.float64, for single use np.float32

DensityMinMaxRange = { 1  : { "Min" : [0.060,0.065], "Max" : [0.493,0.497] },
                       2  : { "Min" : [0.084,0.087], "Max" : [0.493,0.497] },
                       10 : { "Min" : [0.100,0.130], "Max" : [0.350,0.450] },
                     }


def verify_result(path):
    directory = os.path.join(path, 'output')
    info = []
    for snap, minmax in DensityMinMaxRange.items():
        filename = 'snap_%03d.hdf5' % snap
        snap_path = os.path.join(directory, filename)
        try:
            with h5py.File(snap_path, 'r') as data:
                Time    = FloatType(data['Header'].attrs['Time'])
                Density = np.array(data['PartType0']['Density'], dtype=FloatType)
        except (OSError, IOError):
            return False, ['Could not open file "%s"!' % snap_path]
        
        # check minimum
        if Density.min() < minmax["Min"][0] or Density.min() > minmax["Min"][1]:
            return False, [
                "Minimum density at time %g is %g, outside range [%g,%g]." %
                (Time, Density.min(), minmax["Min"][0], minmax["Min"][1])
            ]
        
        # check maximum
        if Density.max() < minmax["Max"][0] or Density.max() > minmax["Max"][1]:
            return False, [
                "Maximum density at time %g is %g, outside range [%g,%g]." %
                (Time, Density.max(), minmax["Max"][0], minmax["Max"][1])
            ]

        info.append(
            "t=%g, minimum density=%g, allowed range=[%g,%g]" %
            (Time, Density.min(), minmax["Min"][0], minmax["Min"][1])
        )
        info.append(
            "t=%g, maximum density=%g, allowed range=[%g,%g]" %
            (Time, Density.max(), minmax["Max"][0], minmax["Max"][1])
        )
    return True, info


if __name__ == '__main__':
    simulation_directory = str(sys.argv[1])
    print(__file__ + ': checking simulation output in directory ' +
          simulation_directory)
    # perform checks
    status_ok, info = verify_result(simulation_directory)
    for msg in info:
        print(msg)
    sys.exit(int(not status_ok))
