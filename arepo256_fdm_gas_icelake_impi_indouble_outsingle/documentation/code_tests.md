Code Tests
==========

A number of code (i.e. integration/functional) tests are available in the ``examples/`` directory.
They can be used to verify that core functionality, accuracy, and/or
correctness of the code has not changed. A small helper script is available
to run one or more tests on the local code base, verify the results, and
optionally visualize the output.

Although such full tests are the ideal case, simpler compilation tests, which
only make sure that the code compiles for a given configuration, are also
present to provide some rudimentary coverage even in the absence of proper
tests. Such tests can be added in the
``examples/sanity_checks/compilation_examples/`` directory by dropping a
``Config.sh`` file in a sub-directory (see also
``examples/sanity_checks/compilation_examples/README.txt``).


Required libraries
------------------

.. |numpy| replace:: ``numpy``
.. _numpy: https://numpy.org/
.. |scipy| replace:: ``scipy``
.. _scipy: https://www.scipy.org/
.. |h5py| replace:: ``h5py``
.. _h5py: https://www.h5py.org/
.. |matplotlib| replace:: ``matplotlib``
.. _matplotlib: https://matplotlib.org/
.. |pip| replace:: ``pip``
.. _pip: https://pip.pypa.io/

The test scripts are written in Python and generally require the libraries
|numpy|_, |scipy|_, and |h5py|_ to read and verify simulation results, and
|matplotlib|_ to create plots. The compilation tests
``examples/sanity_checks/compilation_basic`` additionally use the
`PyYAML <https://pyyaml.org/>`_ library.

These libraries can be installed using the |pip|_ Python package manager::

    pip install numpy
    pip install scipy
    pip install h5py
    pip install matplotlib
    pip install pyyaml

Alternatively, the Anaconda Python distribution also contains all of these
libraries.


Running the test suite
----------------------

Execution of the ``test.sh`` script in the root directory will run a series of tests as
specified (or commented out) therein::

    $ ./test.sh

Alternatively, execution of the more advanced ``test.py`` script::

    $ ./test.py
    usage: test.py [-h] [--withvis] [--number-of-tasks NUMBER_OF_TASKS]
                   [--number-of-compilers NUMBER_OF_COMPILERS] [--no-cleanup]
                   [--print-timestep PRINT_TIMESTEP] [--no-print-output]
                   [--print-all-output]
                   test_name [test_name ...]


For example, to execute all tests sequentially::

    ./test.py all

Or, to run only the Kelvin–Helmholtz instability test in 2D, with
generation of visualization frames (and optionally a movie if ffmpeg is present)::

    ./test.py --withvis khi_hydro_2d

If visualization is requested, a ``vis/`` subdirectory within the test directory will
be made and left for inspection.


Adding a new test
-----------------

Each test resides in a subdirectory with a ``_1d``, ``_2d``, or ``_3d`` suffix, as appropriate.
To add a new test, create a subdirectory with the following required files::

    Config.sh
    param.txt
    create.py
    check.py

If the test does not generate its initial conditions on the fly, and instead requires
pre-made ICs, then the following file should also be present::

    ics.hdf5

The configuration and parameter files should set the code execution behavior. For
simplicity, disabled configuration options can be excluded (i.e. don't keep the whole
``Template-Config.sh`` file). The ``create.py`` file should be a simple Python script, with
a mandatory function: ``create_ics()``. This should be invoked if the script is executed
from the command line, following the existing examples.

If possible, ``create_ics()`` should generate the initial conditions of the test, such
that binary files do not have to be added to the repo. In this case, a file with the
name ``ics.hdf5`` should be created, and the function should return ``True``, to indicate
this.
For this purpose, a simple helper function ``write_ic_file()`` is available in the file ``examples/utils.py`` which will create a HDF5 format ICs file given a dictionary of particle properties (see example below).
Otherwise, if using pre-made ICs, the function can simply return immediately.

The Arepo executable is called from the respective run directory, which contains a copy of the
example directory with all its content and the output created by the ``create.py`` script. Therefore,
the relative path to the output directory (containing the initial condition file) during the Arepo run is ``./``.

In ``check.py``, two functions should be present: ``verify_result()`` and ``visualize_result()``.
The first should be invoked if the script is executed from the command line, and the script should
finish with an exit code ``0`` if successful, and nonzero otherwise.

The ``verify_result()`` function should assess the validity of the output of the test run.
It can do so in any way desired, most likely by comparison to an analytical/exact result, or
by comparison to a previously saved snapshot of the final state of the test, run on the code
base when the test was introduced. It should return a 2-tuple, the first entry being ``False``
to indicate test failure, or ``True`` to indicate success. The second entry is a list, optionally
empty, of string messages giving further information.

The optional ``visualize_result()`` function should do any visualization on the output of the run,
if possible, putting any results into a ``vis/`` subdirectory of the test, which is always kept
for inspection. Two simple helper functions are available, ``visualize_result_2d()`` and
``visualize_result_1d()``, which can be called if the test is of a standard type.
``visualize_result()`` does not need to return any value.

More precisely, here are the function signatures which the functions required by the test framework must be compatible with::

    # create.py:
    # return value is currently not checked
    def create_ics(path, filename='ics.hdf5') -> Optional[bool]:
        ...

    # check.py:
    def verify_result(path) -> Tuple[bool, list]:
        ...

    def visualize_result(path, Lx, Ly) -> Any:
        ...

where ``path`` is the path to the run directory.

Tests can also optionally define the following attributes in ``create.py`` to influence how the test is run::

    # create.py:
    # enforce a fixed number of MPI tasks to run the test with
    # (e.g. 1 for serial 1D tests)
    numTasksMPI = ...
    # plot options which will be used in the parameter file and
    # passed to visualize_result()
    Lx = ...
    Ly = ..

For a full working example of a 1D test, which verifies its result based on an exact
solution, see ``examples/shocktube/shocktube_sod_1d/``. For an example of a 2D test, which verifies its
result based on a previously saved output, see ``examples/khi_hydro_2d``. The easiest
way to add a new test is to start from an existing example.


Test considerations
-------------------

* If possible, generate the initial conditions for the test on the fly. This avoids binary
  files in the repo, enables transparency and reproducibility, and makes it possible to
  modify the test setup or use it in a different context.

* Tests should be short, as the entire test suite should be run on each commit. Ideally,
  stick to configurations which can complete on a single core in a few minutes or less,
  ideally in a few seconds. Memory usage should be minimized.

* Let us add the other common tests from the various technical Arepo papers as time permits,
  i.e. 1D acoustic wave, 1D interacting waves/compression, 2D and 3D Sedov-Taylor explosion,
  2D Gresho vortex, 2D Yee vortex, 2D and 3D Noh, 2D RT, 2D implosion, 1D and multi-D Evrard
  collapse, 1D and multi-D Zeldovich pancake, 2D Orszag-Tang, 1D MHD shocktube, magnetic
  rotor and other MHD tests, 3D turbulence driven boxes. We can also add 3D low resolution
  versions of canonical simulations: isolated disks, spherical cooling gas clouds, cosmological
  boxes, etc.

* Any updates to the framework itself are welcome.

