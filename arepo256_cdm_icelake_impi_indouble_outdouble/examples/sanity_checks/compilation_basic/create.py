#!/usr/bin/env python3
"""
This module is executed as part of the “create ICs” step of running the test
sanity_checks/compilation_basic. It does not, however, create any real ICs
because it only checks if the code compiles correctly. Instead, it creates and
compiles Config.sh files as specified by the file
sanity_checks/compilation_basic/compilation_test_definitions.yml, reporting
any errors encountered during compilation. That is, the actual compilation test
is performed in this file, while the compiled executables are not actually run
and, correspondingly, there are no simulation results to verify in check.py.
"""
import collections
import collections.abc
import os
import subprocess
import sys
from pathlib import Path
import yaml

COMPILE_ONLY = True

DEFAULT_NUMBER_OF_COMPILERS = 4
DEFINITION_FILE_NAME = 'compilation_test_definitions.yaml'
CONFIG_FILE_NAME = 'Config.sh'
EXEC_FILE_NAME = 'Arepo'
BUILD_DIR_NAME = 'build'


class CompilationError(RuntimeError):
    pass


# Define an OrderedSet class because this is (annoyingly) not provided in the
# Python standard library
class OrderedSet(collections.OrderedDict, collections.abc.MutableSet):
    def add(self, elem):
        self[elem] = None

    def discard(self, elem):
        self.pop(elem, None)

    def update(self, *args):
        for s in args:
            for e in s:
                self.add(e)

    def __le__(self, other):
        return all(e in other for e in self)

    def __lt__(self, other):
        return self <= other and self != other

    def __ge__(self, other):
        return all(e in self for e in other)

    def __gt__(self, other):
        return self >= other and self != other

    def __repr__(self):
        return 'OrderedSet([%s])' % (', '.join(map(repr, self.keys())))

    def __str__(self):
        return '{%s}' % (', '.join(map(repr, self.keys())))


def process_definitions(config_definitions):
    result = {}
    for test_name, options in config_definitions.items():
        result[test_name] = OrderedSet()
        for option in options:
            # entries of the form <NAME> mean that all the options defined for
            # the test NAME are copied into this test’s configuration
            if option.startswith('<') and option.endswith('>'):
                result[test_name] |= result[option[1:-1]]
            # entries of the form ~NAME mean that NAME is excluded from this
            # test’s configuration
            elif option.startswith('~'):
                option = option[1:]
                # if a whole test is being excluded, remove all its options
                if option.startswith('<') and option.endswith('>'):
                    result[test_name] -= result[option[1:-1]]
                # otherwise remove the single option itself
                else:
                    result[test_name].discard(option)
            # otherwise, the entry is simply added as a configuration option
            else:
                result[test_name].add(option)
    return result


def write_config(test_dir, options):
    config_path = test_dir / CONFIG_FILE_NAME
    config_path.write_text('\n'.join(options))
    return config_path


def create_ics(path, filename='ics.hdf5', tests_to_run=None):
    # This test does not create any real ICs because it only checks if the code
    # compiles correctly. Instead, it creates and compiles Config.sh files as
    # specified by the file compilation_test_definitions.yml.
    base_dir = Path(path)
    number_of_compilers = (os.environ['NUMBER_OF_COMPILERS']
                           if 'NUMBER_OF_COMPILERS' in os.environ else
                           DEFAULT_NUMBER_OF_COMPILERS)
    # read YAML file
    yaml_path = base_dir / DEFINITION_FILE_NAME
    with open(yaml_path) as yaml_file:
        config_definitions = yaml.load(yaml_file, yaml.BaseLoader)
    # parse configuration options
    configs = process_definitions(config_definitions)
    print(
        f'\tFound [{len(configs)}] compilation test definitions in {yaml_path}'
    )
    # process each test definition
    num_tests = 0
    for test_name, options in configs.items():
        if tests_to_run and test_name not in tests_to_run:
            continue
        print(f'\tRunning compilation test [{test_name}]...')
        # write configuration to Config.sh file
        test_dir = base_dir / test_name
        test_dir.mkdir()
        config_path = write_config(test_dir, options)
        # compile code with given options
        make_result = subprocess.run([
            'make',
            f'-j{number_of_compilers}',
            f'CONFIG={config_path}',
            f'BUILD_DIR={test_dir / BUILD_DIR_NAME}',
            f'EXEC={test_dir / EXEC_FILE_NAME}',
            # skip checks for Template-Config.sh and defines_extra, which are
            # always the same independent of Config.sh
            'build_with_config_check',
        ],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     universal_newlines=True)
        # check for errors
        try:
            make_result.check_returncode()
        except subprocess.CalledProcessError as e:
            config_text = config_path.read_text()
            error_msg = (
                f'\nCompilation test [{test_name}] failed with code '
                f'[{make_result.returncode}].\n\n'
                'Test definition:\n'
                f'{yaml.safe_dump(config_definitions[test_name])}\n'
                f'Contents of {config_path}:\n{config_text}\n\n'
                f'Compilation output:\n{make_result.stdout}')
            if __name__ == '__main__':
                print(f'FAILED: {error_msg}')
                sys.exit(make_result.returncode)
            raise CompilationError(error_msg) from None
        num_tests += 1
    print(f'All [{num_tests}] compilation tests have passed successfully')


if __name__ == '__main__':
    create_ics(sys.argv[1], tests_to_run=sys.argv[2:])
