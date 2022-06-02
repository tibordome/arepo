Getting Started
===============

Prerequisites
-------------

AREPO needs the following non-standard libraries for compilation.

**Required in all cases:**

* **MPI** - the Message Passing Interface (version 1.0 or higher). Many vendor-supplied versions
  exist, in addition to excellent open source implementations, e.g.
  `OpenMPI <https://www.open-mpi.org/>`_ or
  `MPICH <https://www.mpich.org/>`_.

* **GSL** - the `GNU scientific library <https://www.gnu.org/software/gsl/>`_  is a required
  open-source package. AREPO needs this library for a few simple cosmological integrations at
  start-up, and for random number generation.

* **GMP** - the `GNU multiple precision arithmetic library <https://gmplib.org/>`_.
  is a required open-source package. AREPO needs this library for big number integer arithmetic
  in order to calculate exact geometric predicates if needed.

**Optional, required for specific features:**

* **FFTW** - the `Fastest Fourier Transform in the West <http://www.fftw.org>`_ is a required
  open-source library, needed for simulations that use the TreePM algorithm. Version **3.x** of
  FFTW is required. AREPO implements its own FFT parallelization, so it is not required to compile
  FFTW with parallel support.

* **HDF5** - the `Hierarchical Data Format <https://www.hdfgroup.org/solutions/hdf5>`_. This library has
  been developed by NCSA, is open-source, and optional. AREPO can be compiled without this library,
  but then input/output in HDF5 format is not supported. The old binary file formats are no longer
  updated, so HDF5 is highly recommended.

* **hwloc** : The `Portable Hardware Locality library <https://www.open-mpi.org/projects/hwloc/>`_
  is useful for allowing the code to probe the processor topology it is running on and enable a
  pinning of MPI threads to individual cores. This library is optional and only required if the
  ``IMPOSE_PINNING`` flag is set.

Note that if any of the above libraries is not installed in standard locations on your system, the
``Makefile`` provided with the code may need slight adjustments (see information on the ``SYSTYPE`` option below).
Similarly, compiler options, particularly with respect to optimizations, may need adjustment to the C compiler that is used.

The provided ``Makefile`` is compatible with GNU make, so typing ``make`` will then build the
executable assuming the default input configuration filename of ``Config.sh`` and the default
output executable name ``Arepo``. You customize these with the ``CONFIG=`` and ``EXEC=`` options
to make, for example::

  make CONFIG=Config_setup2.sh EXEC=Arepo_setup2


Configuration
-------------

There are two important ways in which the AREPO code is configured and controlled:

1. The ``Config.sh`` file contains a substantial number of compile-time options that need to be set
appropriately for the type of simulation that is to be carried out. The code needs to be
recompiled whenever one of these options is changed.

* See :doc:`core_config_options`.

2. There is also a text parameter file (e.g. ``param.txt``) that is passed as an argument at run-time
to the code. This file specifies a number of variables through simple keyword-value pairs.

* See :doc:`core_param_options`.


Compiling
---------

To build the code, do the following:

#. Copy the file ``Template-Config.sh`` to ``Config.sh``.
#. Edit ``Config.sh`` as needed for your application.
#. Specify the target system using the ``SYSTYPE`` option (see below).
#. Run ``make`` to produce the executable.

The ``Template-Config.sh`` contains a list of all available compile-time options, with most
of them commented out by default. Before compiling for a run, this template should be copied to a
new file called ``Config.sh`` for customziation. To activate a certain feature, the corresponding
parameter can then be commented in, and given the desired value, where appropriate.

.. warning::

  Whenever one of the compile-time options described below is modified, a full
  recompilation of the code may be necessary. To guarantee that this is done when a simple ``make`` is
  specified, all source files have been specified dependent on the Config.sh options. Alternatively,
  one can also issue the command ``make clean``, which will erase all object files, followed by ``make``.

This technique has the disadvantage that different simulations may require different
binaries of AREPO. If several simulations are run concurrently, there is hence the danger that a
simulation is started/resumed with the *wrong* binary. Note that while AREPO checks the plausibility
of some of the most important code options, this is not done for all of them. To minimize the risk of
using the wrong executable for a simulation, it is recommended to produce a separate executable for
each simulation that is run. For example, a good strategy is to make a copy of the whole code together
with its ``Makefile`` in the output directory of each simulation run, and then to use this copy to compile
the code and to run the simulation.

The ``SYSTYPE`` can be set in two ways, with an environment variable or with a file in the
build directory. If present, the file has priority over the environment variable.

* To set with an environment variable e.g. ``export SYSTYPE=Magny`` either at the command
  line or permanently in your ``.bashrc`` (or equivalent) file.

* Or, to set with a file:

  #. copy the ``Template-Makefile.systype`` file to ``Makefile.systype``
  #. uncomment your system in ``Makefile.systype``

If none of the systems in ``Template-Makefile.systype`` are appropriate for your machine, new systems can be added in the file ``makefiles/systypes.make``.

.. note::

  The ``Config.sh`` file should not be checked in to the repository. During development,
  new compile-time options should be added to the ``Template-Config.sh`` file only. Usually, they
  should be added there in the disabled/default version.

Compile-time checks of flags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The normal ``make`` command will first execute a check using the
python script ``check.py``. This ensures that all flags used in the
code via ``#if defined(MyFlag)`` or ``#ifdef MyFlag`` occur either in
``Template-Config.sh`` or ``defines_extra``. The main reason for doing
this is to prevent typos in the code such as ``#if defined(MyFlg)``,
which are generally very hard to find because the code might compile
normally, but the code followed by ``#if defined(MyFlg)`` just does
not appear in the executable.

New flags that can be activated should be included to
``Template-Config.sh``, while internally defined flags (such as header
file guards) are listed in ``defines_extra`` in order to exempt them
from the check.


Running
-------

After compiling AREPO, the executable ``Arepo`` should be present in
the main directory.
In order to start the simulation code, a “parameter file” needs to be specified, usually called ``param.txt``. An additional
optional numerical parameter can be used to signal whether a continuation from a set of restart
files, or from a snapshot file, is desired. A typical command to start the code looks like
the following::

  mpiexec -n 8 ./Arepo param.txt [RestartFlag]

This would start the code using 8 processors, assuming that the parallel environment uses the
``mpiexec`` command to start MPI applications. Depending on the operating system, other commands
may be required for this task, e.g. ``mpirun``. Note that the code can in principle be started
using an arbitrary number of processors, but the communication algorithms will be most efficient for
powers of 2. It is also possible to use a single processor only, in which case the code behaves
like a serial code.

The optional ``RestartFlag`` takes a numeric value, which invokes a partial mode of code operation.
The primary start-up options are just 0, 1, or 2.

* 0: (the default, if omitted) requests the code start from initial conditions.
* 1: signals a continuation from restart files (see `Restarting a run`_).
* 2: can be used to restart from a snapshot file produced by the code (see `Restarting a run`_).

Higher values for the start-up options will execute the special postprocessing options built into
the code. These are currently:

* 3: Execution of FOF/SUBFIND on a snapshot dump.
* 4: Making of an image slice.
* 5: Making of a projected image.
* 6: Conversion of snapshots from one file format to another.
* 7: Calculate a velocity power spectrum for the gas cells.
* 8: Make a grid/orthographic projection (using ray-tracing through the Voronoi mesh).
* 9: Make a projection along an arbitrary axis.
* 10: Make a perspective camera projection.
* 11: Calculate power spectra of various quantities for TRACER_PARTICLEs.
* 12: Calculate two-point correlation function (for a given particle type).
* 13: Calculate power spectrum (for a given particle type).
* 14: Write out the Voronoi mesh.
* 15: Run the post-processing shock finder.
* 16: Write out a two-dimensional slice of the Voronoi mesh.
* 17: Write out snapshot dump with measured gradients.
* 18: Recalculate gravitational potential values.
* 19: Calculate additional quantities from a snapshot dump.
* 20: Render Auriga movie frame from a snapshot dump.
* 21: Run SimpleX RT in post-processing mode on a snapshot.

More information about the additional parameters required for these options can be obtained by
starting AREPO without any options.

Interrupting a run
^^^^^^^^^^^^^^^^^^

AREPO supports running a simulation in several consecutive
installments, and resuming an already started simulation from
so-called restart files (a mechanism sometimes called checkpointing),
which are essentially memory dumps written to the output
directory. The code is written in such a way that the results are
completely independent of the number of restarts performed, which is a
crucial aspect for debugging the code if necessary, as this ensures
deterministic reproducibility.

A (regular) code termination always results in the writing of such
restart files, and is normally triggered when more than 85% of the
specified maximum runtime has been reached, or when the final
simulation time has been reached. It is important, in particular in
larger simulations, that the code is able to complete the output of
the restart file before, e.g. a job time limit is reached. Therefore,
the specified maximum runtime should never exceed the runtime limits
of the machine.

It is also possible to trigger a (regular) code interruption by
introducing a file called ``stop`` in the output directory.::

    touch stop

Furthermore, AREPO usually produces restart files in regular, user-specified time intervals to ensure runs can be resumed without a major loss
of computation time even if a run was terminated abruptly by the
operating system.

Restarting a run
^^^^^^^^^^^^^^^^

A simulation can be restarted in two ways in AREPO, which is indicated with the
additional flags 1 or 2 in the execution call. The flag 1 is the restarting
from a restart file and generally recommended if possible. This, however,
has the restriction that a run needs to be resumed with **exactly as many MPI
ranks** as it was originally started from. The restart command then looks
as follows::

    mpiexec -n 8 ./Arepo param.txt 1

Alternatively, a flag value of 2 signals that one wants to resume a
simulation from a regular snapshot output. This is possible but
discouraged because here the simulation results will not be binary
identical to a run without this restart. In particular in runs where
some fields (especially positions) are written in single precision,
this can also have effects on the overall quality of the simulation. A
command to restart AREPO from a snapshot looks like the following::

    mpiexec -n 8 ./Arepo param.txt 2 7

Here AREPO would read in the snapshot number 7 and resume the
simulation from its corresponding time. This second method is only
recommened if the regular start from restart files is not possible for
some reason.

Run status (“passive files”)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Apart from the ``stop`` file mentioned above, AREPO uses a few other empty files in the output directory as a simple form of communication while running.
They are the following:

* ``stop``: If this file is present, AREPO interrupts the run in a regular fashion, including creation of restart files (see `Interrupting a run`_).
  Afterwards, AREPO removes this file.
* ``restart``: If this file is present, AREPO immediately creates restart files at the end of the current time step *without* interrupting the run.
  Afterwards, AREPO removes this file.
* ``end``: The existence of this file indicates that the simulation has run up to its final time and ended without error.
  AREPO removes this file at the beginning of a normal run.
* ``cont``: The existence of this file indicates that the simulation had to be interrupted because it exceeded the CPU time limit (``TimeLimitCPU`` in ``param.txt``).
  If ``ResubmitOn`` is set, the ``ResubmitCommand`` was executed (see :doc:`core_param_options`).
  AREPO removes this file at the beginning of a normal run.

