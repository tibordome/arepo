
EOS_OPAL
========

OPAL equation of state including ionization to be used for stellar applications.
Reimplementation of the freely available Fortran code in C.


Usage
-----

Switch on `EOS_OPAL` in Config.sh and set `EOS_NSPECIES` to the number of
species that is provided. At least one species is required. The first species is
interpreted as the mass fraction of hydrogen which is used for the EOS call.


Additional Parameters
---------------------

* ``EosOpalTable`` is the filename of the Opal EOS data file. The path is
  expected to be relative to the current working directory. If this is not
  given, the code breaks down.


Additional Config.sh Options
----------------------------

EOS\_OPAL
  Enable for using the Opal equation of state.

EOS\_NSPECIES
  Number of species in xnuc.


Authors
-------

  * Sebastian Ohlmann (sebastian.ohlmann@mpcdf.mpg.de)


Usage Policy and Citation
-------------------------

Please contact the author before using this code for a new project. 

Papers to cite:

  * Ohlmann et al. (2017), Astronomy & Astrophysics, 599, A5
