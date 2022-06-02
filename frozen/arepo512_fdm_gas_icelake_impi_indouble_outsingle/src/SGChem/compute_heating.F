#include "../../build/arepoconfig.h"
      subroutine compute_heating(flux, temp, rho, heating_rate)
#ifdef SGCHEM
      implicit none
#include "cool.h"
      REAL flux, temp, rho, heating_rate
      REAL kappa, specific_heating_rate
c
c Compute planck-mean kappa, given T, n   [cm^2 g^-1]
c
      call compute_md_kappa(temp, rho, kappa)
c
c Compute specific heating rate
c
      specific_heating_rate = flux * kappa
c
c Convert to volumetric heating rate
c
      heating_rate = rho * specific_heating_rate
c
#endif /* SGCHEM */
      return
      end

