#!/usr/bin/env python3
# This file processes the configurations options in Config.sh, producing
# three files:
#
# arepoconfig.h             to be included in each source file (via allvars.h)
# compile_time_info.c       code to be compiled in, which will print the
#                           configuration
# compile_time_info_hdf5.c  code to be compiled in, which will write the
#                           configuration to every HDF5 snapshot
#
from __future__ import print_function
import os.path
import sys
import warnings

COMMENT_CHAR = '#'
ASSIGNMENT_CHAR = '='
INCLUDE_DIRECTIVE = '!include'

if len(sys.argv) != 3:
    print('usage: python prepare-config.py <Config.sh> <build dir>\n')
    sys.exit(1)


def main():
    config_path, build_path = sys.argv[1:]
    outfile_path = os.path.join(build_path, 'arepoconfig.h')
    coutf_path = os.path.join(build_path, 'compile_time_info.c')
    coutf2_path = os.path.join(build_path, 'compile_time_info_hdf5.c')

    with \
     open(outfile_path, 'w') as OUTFILE, \
     open(coutf_path, 'w') as COUTF, \
     open(coutf2_path, 'w') as COUTF2 \
    :
        COUTF.write('#include <stdio.h>\n')
        COUTF.write('void output_compile_time_options(void)\n{\n')
        COUTF.write('printf(\n')

        COUTF2.write('#include <stdio.h>\n')
        COUTF2.write('#include "arepoconfig.h"\n')

        COUTF2.write('#ifdef HAVE_HDF5\n')
        COUTF2.write('#include <hdf5.h>\n')
        COUTF2.write('#include "hdf5_util.h"\n\n')

        COUTF2.write(
            'void write_compile_time_options_in_hdf5(hid_t handle)\n{\n')
        COUTF2.write('hid_t hdf5_dataspace, hdf5_attribute;\n')
        COUTF2.write('double val;\n')
        COUTF2.write('hid_t atype = H5Tcopy(H5T_C_S1);\n')
        COUTF2.write('H5Tset_size(atype, 1);\n')

        read_config(config_path, OUTFILE, COUTF, COUTF2)


def read_config(config_path, OUTFILE, COUTF, COUTF2):
    # Note that this method of parsing implies that there must always be
    # whitespace between any configuration option and a potential comment
    # following it!
    # Further, comments after a configuration option (on the same line) do
    # not technically need to be introduced by a special character;
    # everything after the first whitespace is simply thrown away. By
    # corollary, there must not be any whitespace between a configuration
    # option, the equals sign, and/or its value; configuration options
    # cannot contain any whitespace in general.
    with open(config_path) as config:
        for line_num, line in enumerate(config, start=1):
            fields = line.split()
            # ignore whitespace-only lines and comment lines
            if not fields or fields[0].startswith(COMMENT_CHAR):
                continue

            option = fields[0]
            if option == INCLUDE_DIRECTIVE:
                if len(fields) < 2:
                    warnings.warn(
                        'Found !include directive without value in <%s>, '
                        'line %d' % (config_path, line_num))
                    continue
                read_config(fields[1], OUTFILE, COUTF, COUTF2)
                continue

            subfields = option.split(ASSIGNMENT_CHAR)
            if len(subfields) < 2:
                subfields.append('')
            sub_key, sub_val = subfields[:2]

            OUTFILE.write('#define %s %s\n' % (sub_key, sub_val))
            COUTF.write('"        %s\\n"\n' % option)
            COUTF2.write('hdf5_dataspace = my_H5Screate(H5S_SCALAR);\n')

            # Code is written to COUTF2 directly, so no configuration option
            # may be defined more than once! (Otherwise, the code would attempt
            # to create the same HDF5 attribute more than once.)
            # However, if TOLERATE_WRITE_ERROR is set, it is possible to use
            # duplicate options, since the HDF5 errors will be ignored. (In
            # that case, the option/value that was specified first "wins".)
            if sub_val:
                COUTF2.write(
                    'hdf5_attribute = my_H5Acreate(handle, "%s", '
                    'H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);\n' %
                    sub_key)
                COUTF2.write('val = %s;\n' % sub_val)
                COUTF2.write(
                    'my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, '
                    '"%s");\n' % sub_key)
            else:
                COUTF2.write(
                    'hdf5_attribute = my_H5Acreate(handle, "%s", atype, '
                    'hdf5_dataspace, H5P_DEFAULT);\n' % sub_key)
                COUTF2.write(
                    'my_H5Awrite(hdf5_attribute, atype, "", "%s");\n' %
                    sub_key)

            COUTF2.write('my_H5Aclose(hdf5_attribute, "%s");\n' % sub_key)
            COUTF2.write('my_H5Sclose(hdf5_dataspace, H5S_SCALAR);\n\n')

        COUTF.write('"\\n");\n')
        COUTF.write('}\n')

        COUTF2.write('my_H5Tclose(atype);\n')
        COUTF2.write('}\n')
        COUTF2.write('#endif\n')


main()
