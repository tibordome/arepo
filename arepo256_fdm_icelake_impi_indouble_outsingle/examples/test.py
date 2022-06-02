#!/usr/bin/env python3
"""
  Execute some (or all) simulation tests, creating ICs if necessary.
  After completion, verify results (vs. analytical or previous numerical solutions), optional visualization.
"""
import argparse
import ast
import collections
import enum
import importlib
import multiprocessing
import os
import os.path
import shutil
import subprocess
import sys
import traceback
from pathlib import Path
import numpy as np
import utils

# do not make .pyc files
sys.dont_write_bytecode = True
# do not make __pycache__ folders
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'

SCRIPT_PATH = Path(__file__).resolve()
PYTHON = 'python3'
CREATE_FILE_NAME = 'create.py'
CHECK_FILE_NAME = 'check.py'
CREATE_FUNCTION_NAME = 'create_ics'
CHECK_FUNCTION_NAME = 'verify_result'
IC_FILE_NAME = 'ics.hdf5'
CONFIG_FILE_NAME = 'Config.sh'
PARAM_FILE_NAME = 'param.txt'

# tests which match BLACKLIST are not included when running all tests
BLACKLIST = [
    # TODO: fix failing AMR tests
    'AMR/',
    # TODO: fix TreePM, TreePMWithHighResRegion
    'ForceLawTests/TreePM',
    'ForceLawTests/TreePMWithHighResRegion',
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--withvis',
                    action='store_true',
                    help='create visualization for the test(s) as well')
parser.add_argument(
    '--number-of-tasks',
    type=int,
    default=10,
    help='the number of MPI tasks to use for tests by default (default: '
    '%(default)d)')
parser.add_argument(
    '--number-of-compilers',
    type=int,
    default=8,
    help='the number of processes to use for compilation via make -j (default: '
    '%(default)d)')
parser.add_argument('--no-cleanup',
                    dest='cleanup',
                    action='store_false',
                    help='do not remove run files after successful tests')
parser.add_argument(
    '--print-timestep',
    type=float,
    default=0.1,
    help='the time step interval for printing occasional test progress '
    '(default: %(default)g)')
parser.add_argument('--no-print-output',
                    dest='print_output',
                    action='store_false',
                    help='do not print any text output generated by the tests')
parser.add_argument(
    '--print-all-output',
    action='store_true',
    help='print all text output generated by the test simulations')
parser.add_argument(
    'test_name',
    nargs='+',
    help='The name of a test to run (relative to the examples/ directory). '
    'To run all tests in a row, specify “all”')
args = parser.parse_args()


class TestResult(enum.Enum):
    FAIL = -1
    SKIP = 0
    SUCCESS = 1

CreateAttributes = collections.namedtuple(
    'CreateAttributes', 'compile_only, numTasksMPI, vis_Lx, vis_Ly'
)


def path_in_blacklist(path):
    path_str = str(path)
    for entry in BLACKLIST:
        if entry in path_str:
            return True
    return False


def import_test_module(package_name):
    module = None
    try:
        module = importlib.import_module(package_name)
    except:
        utils.fail('FAILED: Error on importing test modules.')
        utils.fail(' ' + traceback.format_exc())
    return module


def run_create_script(name, run_path, conn):
    pkg_str = '.'.join(name.parts)
    create = import_test_module(pkg_str + '.' + CREATE_FILE_NAME[:-len('.py')])
    if create is None:
        conn.send((TestResult.FAIL, None))
        return
    try:
        getattr(create, CREATE_FUNCTION_NAME)(path=str(run_path),
                                              filename=IC_FILE_NAME)
    except:
        utils.fail('FAILED: Could not generate ICs.')
        utils.fail(' ' + traceback.format_exc())
        conn.send((TestResult.FAIL, None))
        return
    conn.send((TestResult.SUCCESS, CreateAttributes(
        # some tests (optionally) specify that the code should only be compiled
        compile_only=getattr(create, 'COMPILE_ONLY', False),
        # allow test to (optionally) specify the number of MPI tasks
        numTasksMPI=getattr(create, 'numTasksMPI', None),
        # allow test to (optionally) specify visualization options
        vis_Lx=getattr(create, 'Lx', None),
        vis_Ly=getattr(create, 'Ly', None),
    )))


def run_check_script(name, run_path, vis, vis_Lx, vis_Ly, conn):
    pkg_str = '.'.join(name.parts)
    check = import_test_module(pkg_str + '.' + CHECK_FILE_NAME[:-len('.py')])
    if check is None:
        conn.send(TestResult.FAIL)
        return
    try:
        status_ok, info = check.verify_result(str(run_path))
        if not status_ok:
            utils.fail('FAILED.')
            for msg in info:
                utils.fail(' ' + msg)
            conn.send(TestResult.FAIL)
            return
    except:
        utils.fail('FAILED.')
        utils.fail(' ' + traceback.format_exc())
        conn.send(TestResult.FAIL)
        return
    else:
        utils.success('SUCCESS.')
        for msg in info:
            utils.success(' ' + msg)
    # visualization requested? call visualize_result() of the test, leave
    # output
    if vis:
        print(' - Creating visualization(s)...')
        if not (run_path / 'vis').is_dir():
            (run_path / 'vis').mkdir()
        try:
            check.visualize_result(str(run_path), vis_Lx, vis_Ly)
        except:
            print(' - An error occured during creation of visualization:')
            traceback.print_exc()
    conn.send(TestResult.SUCCESS)


def run_test(name, numTasksMPI=args.number_of_tasks, vis=False):
    """Run a single test."""
    name = Path(name)
    test_path = SCRIPT_PATH.parent / name
    if (not (test_path / CREATE_FILE_NAME).is_file()
            or not (test_path / CHECK_FILE_NAME).is_file()):
        print('[%s] test not found, skipping.' % name)
        return TestResult.SKIP

    print('[%s]' % name)

    numTasksMPI_default = numTasksMPI

    arepo_root = SCRIPT_PATH.parent.parent
    run_path_abs = arepo_root / 'run' / 'examples' / name
    # relative to the current working directory
    rel_path = Path(os.path.relpath(run_path_abs))
    # relative to AREPO root
    run_path = run_path_abs.relative_to(arepo_root)

    # prepare directory for run
    if rel_path.is_dir():
        print(' - Warning: Run path [%s] already exists, skipping.' % rel_path)
        return TestResult.SKIP

    shutil.copytree(test_path, rel_path)
    (rel_path / 'output').mkdir()
    (rel_path / 'build').mkdir()

    # check if the functions
    #   - create_ics() in {test}/create.py, and
    #   - verify_result() in {test}/check.py
    # exist
    has_func = False
    ast_module_create = None
    ast_module_check = None
    compile_only = False
    try:
        ast_module_create = ast.parse(
            (test_path / CREATE_FILE_NAME).read_text())
        ast_module_check = ast.parse((test_path / CHECK_FILE_NAME).read_text())
    except SyntaxError:
        pass
    if ast_module_create is not None and ast_module_check is not None:
        has_func_create = False
        for node in ast_module_create.body:
            if isinstance(
                    node,
                    ast.FunctionDef) and node.name == CREATE_FUNCTION_NAME:
                has_func_create = True
                break
        has_func_check = False
        for node in ast_module_check.body:
            if isinstance(
                    node,
                    ast.FunctionDef) and node.name == CHECK_FUNCTION_NAME:
                has_func_check = True
                break
        has_func = has_func_create and has_func_check

    print(' - Generating initial conditions...')

    if has_func:
        # new test format: execute create_ics() in {test}/create.py
        # run in a separate process to avoid side effects from code in test
        # scripts
        spawn = multiprocessing.get_context('spawn')
        conn_recv, conn_send = spawn.Pipe(duplex=False)
        process = spawn.Process(
            target=run_create_script, args=(name, rel_path, conn_send)
        )
        process.start()
        status, create_attrs = conn_recv.recv()
        process.join()

        if status == TestResult.FAIL:
            return TestResult.FAIL
        compile_only = create_attrs.compile_only
        if create_attrs.numTasksMPI is not None:
            numTasksMPI = create_attrs.numTasksMPI
        vis_Lx = create_attrs.vis_Lx
        vis_Ly = create_attrs.vis_Ly

        # vis? update Config.sh and parameter files (except for AMR runs)
        if vis and 'AMR/' not in str(name) and (
            vis_Lx is not None and vis_Ly is not None
        ):
            with open(rel_path / CONFIG_FILE_NAME, 'a') as f:
                f.write('\nVORONOI_FREQUENT_IMAGES\n')
            with open(rel_path / PARAM_FILE_NAME, 'a') as f:
                paramLines = utils.parameterFileVisOptions
                paramLines = paramLines.replace('MAX_X', '%.1f' % vis_Lx)
                paramLines = paramLines.replace('MAX_Y', '%.1f' % vis_Ly)
                f.write(paramLines)
    else:
        # old test format: execute Python file directly
        try:
            cmd = [PYTHON, test_path / CREATE_FILE_NAME, rel_path]
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            if args.print_output and output:
                print(output.decode())
        except subprocess.CalledProcessError as e:
            utils.fail(
                'FAILED: Could not generate ICs (direct call), exit code [%d].'
                % e.returncode)
            if args.print_output:
                utils.fail(' Output:')
                print(str(e.output, 'utf-8'))
            return TestResult.FAIL

    # compile AREPO
    cmd = [
        'make',
        'CONFIG=%s' % (run_path / CONFIG_FILE_NAME),
        'BUILD_DIR=%s' % (run_path / 'build'),
        'EXEC=%s' % (run_path / 'Arepo'), '-j',
        str(args.number_of_compilers)
    ]

    try:
        print(' - Compiling executable...')
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, cwd=arepo_root)
    except subprocess.CalledProcessError as e:
        utils.fail('FAILED: Compilation failed with code [%d].' % e.returncode)
        if args.print_output and e.output:
            utils.fail(' Output:')
            print(str(e.output, 'utf-8'))
        return TestResult.FAIL

    if not (rel_path / 'Arepo').is_file():
        utils.fail(
            'FAILED: Compilation failed to produce an Arepo executable.')
        return TestResult.FAIL

    # execute
    if compile_only:
        print(' - Test is defined with COMPILE_ONLY, not starting run.')
    else:
        print(' - Compilation done, starting run...')

        cmd = ['mpiexec', '-n', str(numTasksMPI), './Arepo', PARAM_FILE_NAME]

        prev_time = -args.print_timestep
        output = ''

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=rel_path)

        while True:
            # poll as simulation runs, print occasional progress
            line = p.stdout.readline().decode('utf-8')
            if not line: break
            if prev_time < 0.01: output += line  # store early output

            if args.print_output:
                if args.print_all_output:
                    print('     ' + line.rstrip())
                else:
                    if 'Time:' in line:
                        cur_time = float(line.split('Time:')[1].split(',')[0])
                        if cur_time - prev_time >= args.print_timestep:
                            print('     ' + line.rstrip())
                            prev_time = cur_time
                            continue
                    if 'TEST-FORCE-LAW' in line:
                        print('     ' + line.rstrip())
                        continue

        if 'TERMINATE:' in output:
            utils.fail('FAILED: Run failed with early termination.')
            if args.print_output:
                utils.fail(' Output:')
                print(output)
            return TestResult.FAIL

        try:
            p.wait(timeout=30)
            # check for run failure
            if p.returncode != 0:
                utils.fail('FAILED: Run failed with exit code [%d].' %
                           p.returncode)
                if args.print_output:
                    utils.fail(' Output:')
                    print(output)
                return TestResult.FAIL
        except subprocess.TimeoutExpired:
            utils.fail('FAILED: Run process did not terminate')
            if args.print_output:
                utils.fail(' Output:')
                print(output)
            return TestResult.FAIL

        # check for failed startup due to parameter file
        if output.count(
                '\n'
        ) < 10000 and 'Error in file %s:' % PARAM_FILE_NAME in output:
            utils.fail('FAILED:')
            for line in output.split('\n'):
                if 'Error' in line: utils.fail(line)
            return TestResult.FAIL

    print(' - Run done, verifying results...')

    if has_func:
        # new test format: execute verify_result() function in {test}/check.py
        # run in a separate process to avoid side effects from code in test
        # scripts (see above)
        spawn = multiprocessing.get_context('spawn')
        conn_recv, conn_send = spawn.Pipe(duplex=False)
        process = spawn.Process(
            target=run_check_script,
            args=(name, rel_path, vis, vis_Lx, vis_Ly, conn_send)
        )
        process.start()
        status = conn_recv.recv()
        process.join()
        if status == TestResult.FAIL:
            return TestResult.FAIL
    else:
        # old test format: execute Python file directly
        try:
            cmd = [PYTHON, test_path / CHECK_FILE_NAME, rel_path, str(vis)]
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            if args.print_output and output:
                for line in output.decode().split('\n'):
                    utils.success(line)
        except subprocess.CalledProcessError as e:
            utils.fail(
                'FAILED: Could not verify results or verification failed (direct call), exit code [%d].'
                % e.returncode)
            if args.print_output and e.output:
                utils.fail(' Output:')
                print(str(e.output, 'utf-8'))
            return TestResult.FAIL
        utils.success('SUCCESS.')

    if args.cleanup and not vis:
        # cleanup
        shutil.rmtree(rel_path)

    return TestResult.SUCCESS


def main():
    # set environment variables
    os.environ['NUMBER_OF_TASKS'] = str(args.number_of_tasks)
    os.environ['NUMBER_OF_COMPILERS'] = str(args.number_of_compilers)
    # which run, or all?
    examples_dir = SCRIPT_PATH.parent
    testNames = []
    if 'all' in args.test_name:
        testNames += [
            path.parent.relative_to(examples_dir)
            for path in sorted(examples_dir.glob('**/' + CREATE_FILE_NAME),
                               key=lambda p: str(p).casefold())
            if not path_in_blacklist(path)
        ]
    testNames += [Path(name) for name in args.test_name if name != 'all']
    # run tests sequentially, store success/failure status of each
    results = []
    for testName in testNames:
        cwd = Path.cwd()
        sys_path = sys.path.copy()
        results.append(run_test(testName, vis=args.withvis))
        # restore potentially modified CWD and sys.path
        if cwd != Path.cwd():
            os.chdir(cwd)
        sys.path = sys_path
        print()
    assert len(testNames) == len(results)
    # print summary if we ran more than one test
    if len(testNames) > 1:
        status = utils.success if results.count(
            TestResult.FAIL) == 0 else utils.fail
        status(
            'All tests done. [%d] skipped, [%d] succeeded, [%d] failed.' %
            (results.count(TestResult.SKIP), results.count(
                TestResult.SUCCESS), results.count(TestResult.FAIL)))
        if status == utils.fail:
            utils.fail('Failed tests:')
        for test_idx in range(len(testNames)):
            if results[test_idx] == TestResult.FAIL:
                utils.fail('\t[%s]' % testNames[test_idx])


if __name__ == '__main__':
    main()
