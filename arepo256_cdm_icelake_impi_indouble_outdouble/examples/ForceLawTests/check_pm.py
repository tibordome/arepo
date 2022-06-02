#!/usr/bin/env python3
""" Idealized test, check force accuracy: TreePM. """
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

    data = np.loadtxt(file_path)

    # columns (w/ PMGRID): type, task, id, ti, x, y, z, r, fx, fy, fz,
    #   fsx, fsy, fsz, flx, fly, flz, ftx, fty, ftz, fpmx, fpmy, fpmz, p, ps, pl, pt, ppm
    ptype, task, pid, ti = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    x, y, z, r = data[:, 4], data[:, 5], data[:, 6], data[:, 7]
    fx, fy, fz = data[:, 8], data[:, 9], data[:, 10]  # direct summation force
    fsx, fsy, fsz = data[:, 11], data[:,
                                      12], data[:,
                                                13]  # direct short-range force
    flx, fly, flz = data[:, 14], data[:,
                                      15], data[:,
                                                16]  # direct long-range force
    ftx, fty, ftz = data[:, 17], data[:, 18], data[:, 19]  # tree force
    fpx, fpy, fpz = data[:, 20], data[:, 21], data[:, 22]  # pm force
    p, ps, pl, pt, pm = data[:,
                             23], data[:,
                                       24], data[:,
                                                 25], data[:,
                                                           26], data[:,
                                                                     27]  # potentials

    # calculate force and potential errors
    glx = fx - fsx  # difference of direct and short-range forces
    gly = fy - fsy
    glz = fz - fsz
    ffx = ftx + fpx  # TreePM total force
    ffy = fty + fpy
    ffz = ftz + fpz
    ptm = pt + pm  # TreePM total potential

    f = np.sqrt(fx**2 + fy**2 + fz**2)
    fs = np.sqrt(fsx**2 + fsy**2 + fsz**2)
    fl = np.sqrt(flx**2 + fly**2 + flz**2)
    gl = np.sqrt(glx**2 + gly**2 + glz**2)
    ft = np.sqrt(ftx**2 + fty**2 + ftz**2)
    fp = np.sqrt(fpx**2 + fpy**2 + fpz**2)
    ff = np.sqrt(ffx**2 + ffy**2 + ffz**2)

    df_t = np.sqrt((fx - ffx)**2 + (fy - ffy)**2 + (fz - ffz)**2)
    df_s = np.sqrt((fsx - ftx)**2 + (fsy - fty)**2 + (fsz - ftz)**2)
    df_l = np.sqrt((flx - fpx)**2 + (fly - fpy)**2 + (flz - fpz)**2)

    dp_t_rel = np.abs((p - ptm) / p)
    dp_s_rel = np.abs((ps - pt) / ps)
    dp_l_rel = np.abs((pl - pm) / pl)

    w = np.where(f > 0.0)
    df_t_rel = df_t[w] / f[w]

    w = np.where(fs > 0.0)
    df_s_rel = df_s[w] / fs[w]

    w = np.where(fp > 1e-10)
    df_l_rel = df_l[w] / fp[w]

    # force is zero once per cycle
    error_msgs = []

    w = np.where(f == 0.0)
    if len(w[0]) != 40:
        error_msgs.append(
            'Unexpected number of zero forces %d (should be 40).' % len(w[0]))

    w = np.where(fs <= 1e-40)
    if len(
            w[0]
    ) > 700:  # todo, seems strange (distribution of fs is continuous down to small values, then underflow)
        error_msgs.append(
            'Unexpected number of zero direct short-range forces %d (should be 40).'
            % len(w[0]))

    w = np.where(fp < 1e-10)  # nonzero but small
    if len(w[0]) != 40:
        error_msgs.append(
            'Unexpected number of near-zero PM forces %d (should be 40).' %
            len(w[0]))

    # validate
    if not df_t_rel.max() <= 0.02:
        error_msgs.append(
            'Relative total force error too large, maximum of %g.' %
            df_t_rel.max())
    if not df_s_rel.max(
    ) <= 3e10:  # todo, seems strange (~1e10 in TreePM, ~1e2 in TreePMWithHighResRegion)
        error_msgs.append(
            'Relative short-range force error too large, maximum of %g.' %
            df_s_rel.max())
    if not df_l_rel.max() <= 0.04:
        error_msgs.append(
            'Relative long-range force error too large, maximum of %g.' %
            df_l_rel.max())
    if not dp_t_rel.max() <= 0.5:
        error_msgs.append(
            'Relative total potential error too large, maximum of %g.' %
            dp_t_rel.max())
    if not dp_s_rel.max(
    ) <= 1e3:  # todo, seems strange (~1e3), TreePMWithHighResRegion
        error_msgs.append(
            'Relative short-range potential error too large, maximum of %g.' %
            dp_s_rel.max())
    if not dp_l_rel.max(
    ) <= 3e5:  # todo, seems strange (~1e5), TreePM and TreePMWithHighResRegion
        error_msgs.append(
            'Relative long-range potential error too large, maximum of %g.' %
            dp_l_rel.max())

    if df_t.max() == 0 or df_s.max() == 0 or df_l.max() == 0:
        error_msgs.append(
            'Error is absolutely zero! Not likely, needs to be checked.')

    if error_msgs:
        return False, error_msgs

    # return success
    return True, []


def visualize_result(path, Lx, Ly):
    pass


if __name__ == '__main__':
    path = sys.argv[1]
    print(__file__ + ': checking output in directory ' + path)

    status_ok, info = verify_result(path)

    if status_ok:
        print('check.py: success!')
        sys.exit(0)  # normal exit

    for msg in info:
        print('check.py: failed! ' + msg)
    sys.exit(1)
