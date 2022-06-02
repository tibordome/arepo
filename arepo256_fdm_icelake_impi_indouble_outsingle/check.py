#!/usr/bin/env python3
''' check.py: helper to maintain code option consistency '''
import glob
import re
import sys
import os

# workaround for BlockingIOError on print() with make -j
import fcntl
fcntl.fcntl(1, fcntl.F_SETFL, 0)


def parseIf(string, defines, fin):
    ''' Parse pre-processor line from a source file. '''
    s = string
    while s:
        if s.startswith('defined'):
            s = s[7:]
            continue
        if s.startswith('//'):
            return
        m = re.search(r'[a-zA-Z_][a-zA-Z_0-9]*', s)
        if m is not None and m.start() == 0:
            defines.add(m.group())
            s = s[m.span()[1]:]
            continue
        if s.startswith('/*'):
            m = re.search(r'\*/', s)
            if m is not None:
                s = s[m.span()[1]:]
                continue
            else:
                return
        if s[0] == '\n':
            return
        if s == '\\\n':
            string = fin.readline()
            s = string
            continue
        if s[0] in '<>+-/=!&|() 0123456789\t':
            s = s[1:]
            continue
        print("Strange character in '%s' detected: '%s', skipping it." %
              (string, s[0]))
        s = s[1:]


def load(fin):
    ''' Wrapper for loading a list of strings from file. '''
    return set(i.strip() for i in fin if not i.startswith('#'))


def write(s, fout):
    ''' Wrapper for writing out a sorted list of strings to file. '''
    with open(fout, 'w') as fout:
        for i in sorted(s):
            fout.write(i + os.linesep)


# ----------------------------


def filter_code_ifdef(fin):
    ''' Find all macros used in #if/#ifdef/etc. directives in .c or .h files. '''
    defines = set()
    line = fin.readline()
    while line:
        s = line.lstrip()
        if s.startswith('#if ') or s.startswith('#if('):
            parseIf(s[4:], defines, fin)
        elif s.startswith('#elseif ') or s.startswith('#elseif('):
            parseIf(s[8:], defines, fin)
        elif s.startswith('#elif ') or s.startswith('#elif('):
            parseIf(s[6:], defines, fin)
        elif s.startswith('#ifdef ') or s.startswith('#ifndef '):
            if s.startswith('#ifdef '):
                s = s[7:]
            else:
                s = s[8:]
            s = s.lstrip()
            m = re.search(r'[a-zA-Z_][a-zA-Z_0-9]*', s)
            if m is not None and m.start() == 0:
                defines.add(m.group())
            else:
                print("Strange #ifdef/#ifndef: '%s'. Skipping it.", s)
        line = fin.readline()
    return defines


def filter_code_params():
    ''' Find all parameter names that can be read from a parameter file. '''
    with open('src/parameters.c', 'r') as fin:
        s = fin.read()
    d = re.findall(r'strcpy\(tag\[nt\],\s?"[a-zA-Z_0-9]+"\);', s)
    params = set(dd.split('"')[1] for dd in d)
    return params


def filter_template_config(fin):
    ''' Find all items of Template-Config.sh file. '''
    defines = set()
    for line in fin:
        s = line.split()
        if s:
            d = re.findall(r'^#*([a-zA-Z_][a-zA-Z_0-9]*)', s[0])
            defines.update(d)
    return defines


def filter_config(fin):
    ''' Find all active items on Config.sh file. '''
    defines = set()
    for line in fin:
        s = line.split()
        if s:
            d = re.findall(r'^([a-zA-Z_][a-zA-Z_0-9]*)', s[0])
            defines.update(d)
    return defines


def filter_makefile(fin):
    ''' Find all macros used in the given Makefile. '''
    defines = set()
    r = re.compile(
            r'^ifn?eq\s*\(\s*((([a-zA-Z_][a-zA-Z_0-9]*)\s*,)|,?)\s*'
            r'\$\(findstring\s+([a-zA-Z_][a-zA-Z_0-9]*)\s*,\s*'
            r'\$\(CONFIGVARS\)\s*\)\s*,?\s*\)')
    for line in fin:
        d = r.findall(line.strip())
        for groups in d:
            for match in groups:
                if match and ',' not in match:
                    defines.add(match)
    return defines


def filter_config_documentation():
    ''' Parse for all documented Config.sh options (in documentation/*). These can appear either
        as base options in core_config_options.md or as specialized options in any individual
        modules_*.md file. '''
    defines = set()
    files = ['documentation/core_config_options.md']
    files += glob.glob('documentation/modules_*.md')
    for filename in files:
        with open(filename, 'r') as fin:
            s = fin.read()
        d = re.findall(r'\n\n[a-zA-Z_0-9]+\n  ', s)
        defines.update(dd.strip() for dd in d if dd.strip() != 'OPTION')
    return defines


def filter_params_documentation():
    ''' Parse for all documented Config.sh options (in documentation/*). These can appear either
        as base options in core_param_options.md or as specialized options in any individual
        modules_*.md file. '''
    defines = set()
    # syntax in core_param_options.md file: same as Config.sh options above
    with open('documentation/core_param_options.md', 'r') as fin:
        s = fin.read()
    d = re.findall(r'\n\n[a-zA-Z_0-9]+\n  ', s)
    defines.update(dd.strip() for dd in d)
    # syntax in individual module files: * ``ParamName`` description here.
    files = glob.glob('documentation/modules_*.md')
    for filename in files:
        with open(filename, 'r') as fin:
            s = fin.read()
        d = re.findall(r'\n\* ``[a-zA-Z_0-9]+``', s)
        defines.update(dd.strip()[4:-2] for dd in d if dd.strip()[4:-2] != 'Any')
    return defines


# ----------------------------


def check_code(fin, fout, template, extra):
    ''' Check source files for illegal macros. '''
    allowed = filter_template_config(template)
    allowed.update(filter_config(extra))
    used = filter_code_ifdef(fin)
    diff = sorted(used - allowed)
    if diff:
        print()
        print('Illegal macros/options detected in file %s.' % fin.name)
        print(
            "Check for potential typos and add them either to the file 'Template-Config.sh' or 'defines_extra'"
        )
        print(
            "('defines_extra' is for macros which are either internally defined or should not appear in Template-Config.sh)."
        )
        print(
            "In case you want to suppress this check, build with 'make build' instead of 'make'."
        )
        print()
        for i in diff:
            print(i)
        print()
        return 1
    write(used, fout)


def check_makefile(fin, fout, template, extra):
    ''' Check Makefile for illegal options. '''
    allowed = filter_template_config(template)
    allowed.update(filter_config(extra))
    used = filter_makefile(fin)
    diff = sorted(used - allowed)
    if diff:
        print()
        print('Illegal macros/options detected in file %s.' % fin.name)
        print(
            "Check for potential typos and add them either to the file 'Template-Config.sh' or 'defines_extra'"
        )
        print(
            "('defines_extra' is for macros which are either internally defined or should not appear in Template-Config.sh)."
        )
        print(
            "In case you want to suppress this check, build with 'make build' instead of 'make'."
        )
        print()
        for i in diff:
            print(i)
        print()
        return 1
    write(used, fout)


def check_config(fin, fout, args):
    ''' Check Config.sh for illegal options (those which don't appear in any source files). '''
    allowed = set()
    for file in args:
        with open(file, 'r') as f:
            allowed.update(load(f))
    used = filter_config(fin)
    diff = sorted(used - allowed)
    if diff:
        print()
        print(
            'The following options are active in %s, but are not used in any of the'
            % fin.name)
        print('source code files being compiled into the final executable.')
        print('Please check for typos and deactivate the options.')
        print(
            "In case you want to suppress this check, build with 'make build' instead of 'make'."
        )
        print()
        print()
        for i in diff:
            print(i)
        print()
        return 1
    write(used, fout)


def check_documentation(fin, fout):
    ''' Check whether all Template-Config.sh and parameter file options are documented. '''
    ex = False
    # Template-Config.sh
    documented = filter_config_documentation()
    used = filter_template_config(fin)
    diff = sorted(used - documented)
    if diff:
        print()
        print('WARNING: The following options are undocumented, but appear in',
              fin.name)
        print('Please add a proper documentation description:')
        print()
        for i in diff:
            print(i)
        print()
        ex = True
    diff = sorted(documented - used)
    if diff:
        print()
        print(
            'WARNING: The following options are documented, but are (not used/deleted/misspelled) in %s:'
            % fin.name)
        print()
        for i in diff:
            print(i)
        print()
        ex = True
    # parameter file
    documented = filter_params_documentation()
    allowed = filter_code_params()
    diff = sorted(allowed - documented)
    if diff:
        print()
        print(
            'WARNING: The following parameters are undocumented, but appear in src/parameter.c'
        )
        print('Please add a proper documentation description:')
        print()
        for i in diff:
            print(i)
        print()
        ex = True
    diff = sorted(documented - allowed)
    if diff:
        print()
        print(
            'WARNING: The following parameters are documented, but are (not used/deleted/misspelled) in src/parameter.c:'
        )
        print()
        for i in diff:
            print(i)
        print()
        ex = True
    if ex:
        return 1
    write(used, fout)


def check_template_unused(fin, fout, args):
    ''' Check Template-Config.sh and defines_extra for unused options (those
        which don't appear in any source files or Makefiles). '''
    allowed = set()
    for file in args:
        with open(file, 'r') as f:
            allowed.update(load(f))
    used = filter_template_config(fin)
    diff = sorted(used - allowed)
    if diff:
        print()
        print(
            'The following options are defined in %s, but are not used in any of the'
            % fin.name)
        print('source code files or Makefiles.')
        print('Please check for typos and remove unused options.')
        print(
            "In case you want to suppress this check, build with 'make build' instead of 'make'."
        )
        print()
        print()
        for i in diff:
            print(i)
        print()
        return 1
    write(used, fout)


def check_template_duplicates(template, extra, fout):
    ''' Check for duplicates between Template-Config.sh and defines_extra. '''
    template_opts = load(template)
    extra_opts = load(extra)
    intersection = sorted(template_opts & extra_opts)
    if intersection:
        print()
        print(
            'The following options are defined in both %s and %s.'
            % (template.name, extra.name))
        print('Please check for typos and remove duplicated options.')
        print(
            "In case you want to suppress this check, build with 'make build' instead of 'make'."
        )
        print()
        print()
        print(repr(intersection))
        for i in intersection:
            print(i)
        print()
        return 1
    write(template_opts | extra_opts, fout)


if __name__ == '__main__':
    # check command line arguments
    if len(sys.argv) < 3:
        print('Too few command line arguments')
        sys.exit(1)
    mode = sys.argv[1]
    if mode not in (
        '1', '2', '3', '4', '5', '6',
        'code', 'config', 'makefile', 'documentation', 'template-unused',
         'template-duplicate',
    ):
        print('Unknown mode')
        sys.exit(1)
    # process
    with open(sys.argv[2], 'r') as fin:
        fout = sys.argv[3]
        rest = sys.argv[4:]
        if mode in ('1', 'code'):
            print('Checking %s for illegal define macros' % fin.name)
            with open(rest[0], 'r') as template, open(rest[1], 'r') as extra:
                result = check_code(fin, fout, template, extra)
        elif mode in ('2', 'config'):
            print('Checking active options of %s' % fin.name)
            result = check_config(fin, fout, rest)
        elif mode in ('3', 'makefile'):
            print('Checking %s for illegal define macros' % fin.name)
            with open(rest[0], 'r') as template, open(rest[1], 'r') as extra:
                result = check_makefile(fin, fout, template, extra)
        elif mode in ('4', 'documentation'):
            print('Checking %s for undocumented options' % fin.name)
            result = check_documentation(fin, fout)
        elif mode in ('5', 'template-unused'):
            print('Checking %s for unused options' % fin.name)
            result = check_template_unused(fin, fout, rest)
        elif mode in ('6', 'template-duplicate'):
            template, extra, fout = sys.argv[2:]
            print('Checking %s and %s for duplicate options' %
                (template, extra))
            with open(extra, 'r') as fextra:
                result = check_template_duplicates(fin, fextra, fout)
    # potentially exit with error status
    if result:
        sys.exit(result)
