#!/usr/bin/env python3
""" Idealized test, check force accuracy: TreeEwald. """
import os.path
from os.path import isfile, getsize
import numpy as np
import sys


def verify_result(path):
    """ Unit test function, validate results. """
    file_path = os.path.join(path, 'output', 'forcetest.txt')

    # load
    if not isfile(file_path) or getsize(file_path) == 0:
        return False, ['Expected test file apparently not produced.']

    # columns (no PMGRID): type, task, id, ti, x, y, z, r, fx, fy, fz, ftx, fty, ftz, p, pt
    data = np.loadtxt(file_path)

    # check
    ptype, task, pid, ti = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    x, y, z, r = data[:, 4], data[:, 5], data[:, 6], data[:, 7]
    fx, fy, fz = data[:, 8], data[:, 9], data[:, 10]  # direct summation force
    ftx, fty, ftz = data[:, 11], data[:, 12], data[:, 13]  # tree force
    p = data[:, 14]  # direct summation potential
    pt = data[:, 15]  # direct tree potential

    f = np.sqrt(fx**2 + fy**2 + fz**2)
    ft = np.sqrt(ftx**2 + fty**2 + ftz**2)
    df = np.sqrt((fx - ftx)**2 + (fy - fty)**2 + (fz - ftz)**2)
    dp = np.abs(p - pt)

    # force is zero once per cycle
    w = np.where(f == 0.0)
    if len(w[0]) != 40:
        error_msgs.append(
            'Unexpected number of zero forces %d (should be 40).' % len(w[0]))

    w_nonzero = np.where(f > 0.0)
    df_rel = df[w_nonzero] / f[w_nonzero]
    dp_rel = np.abs(dp / p)

    # validate
    error_msgs = []

    if not df.max() <= 1.1e-5:
        error_msgs.append('Force error too large, maximum of %g.' % df.max())
    if not df_rel.max() <= 0.0004:
        error_msgs.append('Relative force error too large, maximum of %g.' %
                          df_rel.max())
    if not dp.max() <= 1.3:
        error_msgs.append('Potential error too large, maximum of %g.' %
                          dp.max())
    if not dp_rel.max() <= 0.2:
        error_msgs.append(
            'Relative potential error too large, maximum of %g.' %
            dp_rel.max())

    info = []
    info.append('force error max: %g' % df.max())
    info.append('rel force error max: %g' % df_rel.max())
    info.append('potential error max: %g' % dp.max())
    info.append('rel pot error max: %g' % dp_rel.max())

    if df.max() == 0 or dp.max() == 0:
        error_msgs.append(
            'Error is absolutely zero! Not likely, needs to be checked.')

    if error_msgs:
        return False, info + error_msgs

    # return success
    return True, info


def visualize_result(path, Lx, Ly):
    pass


if __name__ == '__main__':
    path = sys.argv[1]
    print(
        'examples/ForceLawTests/TreeEwald/check.py: checking output in directory '
        + path)

    status_ok, info = verify_result(path)

    if status_ok:
        print('check.py: success!')
        sys.exit(0)  # normal exit

    for msg in info:
        print('check.py: failed! ' + msg)
    sys.exit(1)
