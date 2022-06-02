#!/usr/bin/env python3
'''
This module is executed as part of the “verify results” step of running the test
sanity_checks/compilation_basic. It does not, however, verify any results,
because this test only checks if the code compiles correctly, i.e. the compiled
executables are not actually run and, correspondingly, there are no simulation
results to verify.
The test itself, including any checks for compilation errors, is performed in
create.py.
'''


def verify_result(path):
    # this test does not verify any results because it only checks if the code
    # compiles correctly
    return True, []


def visualize_result(path, Lx, Ly):
    pass
