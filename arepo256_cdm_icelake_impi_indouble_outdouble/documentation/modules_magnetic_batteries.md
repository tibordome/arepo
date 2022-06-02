
MAGNETIC_BATTERIES
==================

The module for magnetic batteries includes methods for generating magnetic fields via 
the Biermann (1950) battery, the Durrive and Langer (2015) battery, and injection
from SN explosions.


Usage
-----

The module provides independent switches for the three magnetic seeding mechanisms, 
plus a fourth one to include new gradients in the output. The Biermann and Durrive 
batteries introduce a new source term, while the SN-injection scheme builds upon
the ``GFM`` module and includes magnetic fields in the enrichment scheme from SN.


Additional Parameters
---------------------

* ``SupernovaInjectedMagneticEnergyInErgs`` Energy  injected in a dipolar magnetic field from each SN event.


Additional Config.sh Options
----------------------------

BIERMANN_BATTERY
  activates a source term implementing the Biermann (1950) battery.

DURRIVE_BATTERY
  activates a source term implementing the Durrive & Langer (2015) battery.

GFM_INJECT_B_FROM_SN
  activates the injection of a dipolar magnetic field during each SN event. The dipole is centered on the SN position, 
  has a randomly-oriented magnetic dipople, and its strength is set by the total energy injected (read in as a new simulation parameter)

MAGNETIC_BATTERIES_OUTPUT_GRADIENTS
  includes the electron number density, pressure, and momentum transfer rate gradients in the output


Authors
-------

* Enrico Garaldi


Usage Policy and Citation
-------------------------

Users are free to use this module  in their own work. We ask that papers based on this work cite

* Garaldi E., Pakmor R., Springel V., MNRAS, in press (2021)
