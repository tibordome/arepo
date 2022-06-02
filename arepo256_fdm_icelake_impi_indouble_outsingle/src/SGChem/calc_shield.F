#include "arepoconfig.h"
#ifdef SGCHEM
      subroutine calc_shielding(yn, dl, temp, abh2, abco, abcI,
     $                          column_density_projection,
     $                          column_density_projection_H2,
     $                          column_density_projection_co,
     $                          column_density_projection_c,
     $                          fshield_H2, fshield_CO, fshield_C,
     $                          AV_mean, chi_mean)
      implicit NONE
#include "cool.h"
      REAL yn, dl, abh2, abco, abcI, temp
      DOUBLE PRECISION column_density_projection(NPIX)
      DOUBLE PRECISION column_density_projection_H2(NPIX)
      DOUBLE PRECISION column_density_projection_co(NPIX)
      DOUBLE PRECISION column_density_projection_c(NPIX)
      REAL fshield_H2, fshield_CO, fshield_C
      REAL AV_mean, chi_mean, AV, chi
      REAL local_NH, AV_local, local_NH2, local_NCO, local_NC, f_dust
      external chi
      integer I, nray, nray_start, nray_end
c
c Optically thin approx:
c
      if (iphoto .eq. 0) then
        fshield_H2 = 1.0
        fshield_CO = 1.0
        fshield_C  = 1.0
        AV_mean    = 0.0
        chi_mean   = 1.0
        return
      endif
c
      if (iphoto .ne. 6) then
        print*, "Unknown photochemistry option: ", iphoto
        ABORT(3)
      endif

      if (isnan(AV)) print*,"calc_shield.F AV NaN",AV,yn
c
c Total column density of HI nuclei within characteristic scale-length dl
c
      local_NH = dl * yn
      AV_local = AV_conversion_factor * dust_to_gas_ratio 
     $         * local_NH
c
c H2 self-shielding
c
      local_NH2 = local_NH * abh2
      if (local_NH2 .lt. 0.0) then
         print*,"local_NH2 very wrong ",local_NH2,local_NH,dl,yn,
     &        abh2
         local_NH2=0.0
      endif
      call calc_H2_shielding_factor(temp, local_NH2, 
     $                              column_density_projection_h2, 
     $                              fshield_H2)

      if (isnan(fshield_H2)) print*,"calc_shield.F fshield_H2 is NaN:",
     $         fshield_H2,yn

c
c CO self-shielding, shielding of CO by H2
c
#if CHEMISTRYNETWORK != 1
      local_NCO = local_NH * abco
      call calc_CO_shielding_factor(column_density_projection_h2, 
     $                              column_density_projection_co, 
     $                              local_NCO, local_NH2, fshield_CO)
#else
      fshield_CO = 1d0
#endif
c
c C self-shielding, shielding of C by H2
c
#if CHEMISTRYNETWORK != 1
      local_NC = local_NH * abcI
      call calc_C_shielding_factor(column_density_projection_h2,
     $                             column_density_projection_c,
     $                             local_NC, local_NH2, fshield_C)
#else
      fshield_C = 1d0
#endif
c
c Dust shielding -- we compute a mean AV, assuming a weighting factor exp(-2.5AV); this
c is appropriate for photoelectric heating and CO photodissociation, shouldn't be too
c much in error for other things
c
#if CHEMISTRYNETWORK != 1
      nray     = NPIX
      f_dust   = 0d0
      chi_mean = 0d0
      do I = 1, nray
        AV = AV_conversion_factor * column_density_projection(I) 
     $     * dust_to_gas_ratio + AV_local + AV_ext
        if (isnan(AV)) print*,"calc_shield.F AV loop NaN",AV,I,'/',nray,
     $       column_density_projection(I)
        f_dust   = f_dust + (1d0 / nray) * dexp(-2.5d0 * AV)
        chi_mean = chi_mean + (1d0 / nray) * chi(AV)
      enddo
      AV_mean  = -log(f_dust) / 2.5d0
      if (AV_mean .lt. 0d0) then
         chi_mean = 1d0
         AV_mean  = 0d0
      endif
#else
      chi_mean = 1d0
      AV_mean  = 0d0
#endif
      return
      end
#endif /* SGCHEM */
c=======================================================================
c
c    \\\\\\\\\\          E N D     F U N C T I O N        //////////
c    //////////         C A L C _ S H I E L D I N G       \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////                C H I                      \\\\\\\\\\
c
c=======================================================================
c
      REAL function chi(AV)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "isrf.h"
      REAL AV, dchi, dAV
      integer index

      if (AV .ne. AV) then
        print*, "Error: AV is NaN"
        stop
        
      endif

      if (AV .eq. 0) then
        chi = 1d0
        return
      elseif (AV .gt. 2d2) then
c Assume chi remains constant for AV > 200
        if (isrf_option .eq. 1) then
          chi = chi_tab_mathis(200)
        else
          chi = chi_tab_black(200)
        endif
        return
      else
        index = int(AV)
        if ( index.eq. 200 ) then
           print *, "Index = 200... Whoops."
           stop
        end if 
        if (isrf_option .eq. 1) then
          if (index .eq. 0) then
            dchi = chi_tab_mathis(1) - 1d0
            dAV  = AV - index
            chi  = 1d0 + dAV * dchi
          else
            dchi = chi_tab_mathis(index+1) - chi_tab_mathis(index)
            dAV  = AV - index
            chi  = chi_tab_mathis(index) + dAV * dchi
          endif
        else
          if (index .eq. 0) then
            dchi = chi_tab_black(1) - 1d0
            dAV  = AV - index
            chi  = 1d0 + dAV * dchi
          else
            dchi = chi_tab_black(index+1) - chi_tab_black(index)
            dAV  = AV - index
            chi  = chi_tab_black(index) + dAV * dchi
          endif
        endif
        return
      endif
#else
      chi = 1.0
#endif /* SGCHEM */
      return
      end
c
c=======================================================================
c
c
c    \\\\\\\\\\          E N D   F U N C T I O N          //////////
c    //////////                  C H I                    \\\\\\\\\\
c
c=======================================================================
c
      subroutine calc_H2_shielding_factor(temp, local_NH2,
     $                                    column_density_projection_h2, 
     $                                    fshield_H2)

#ifdef SGCHEM
      implicit none
#include "cool.h"
      REAL fshield_H2, temp, local_NH2
      REAL column_density_projection_H2(NPIX)
      REAL b, b5, h2dis0, h2dis1, alpha, term1, term2
      integer I, nray
c
c alpha depends on choice of H2 self-shielding function
c In the original DB96 function, for rotationally cold H2, alpha = 2.0
c In the revised WGH11 function, for rotationally hot H2,  alpha = 1.1
c
      alpha = 2.0
c
      nray = NPIX
c
c To compute H2 self-shielding, we assume that the local b param. is 
c representative of the whole flow. This isn't true in general, but 
c is still a minor approximation compared to the others we are making
c
      b  = dsqrt(kboltz * temp / mh)
      b5 = b / 1d5
#ifdef SGCHEM_GONG17_TEST
      b5 = 3d0
#endif
      fshield_H2 = 0d0
      do I = 1, NPIX
        h2dis0  = (local_NH2 + column_density_projection_h2(I)) / 5d14
        h2dis1  = dsqrt(1d0 + h2dis0)
        term1  = 0.965 / (1d0 + h2dis0/b5)**alpha
        term2  = 0.035 * dexp(-8.5d-4 * h2dis1) / h2dis1
        fshield_H2 = fshield_H2 + (term1 + term2) / nray
        if (isnan(fshield_H2)) print*, "Calc H2 shielding", I, local_NH2
     &       ,column_density_projection_h2(I),b5,h2dis0
      enddo
#endif /* SGCHEM */
      return
      end
c
c
c
      subroutine calc_CO_shielding_factor(column_density_projection_h2, 
     $                                    column_density_projection_co, 
     $                                    local_NCO, local_NH2, 
     $                                    fshield_CO)
#ifdef SGCHEM
      implicit none
#include "cool.h"
#include "shield_data.h"
#include "shield_data_old.h"
      REAL column_density_projection_h2(NPIX)
      REAL column_density_projection_co(NPIX)
      REAL fshield_CO, fshield_CO_ray, NCO, NH2
      REAL local_NCO, local_NH2
      integer nray, I
c
c CO self-shielding, shielding by H2.
c
      fshield_CO = 0d0
      nray       = NPIX

      do I = 1, nray
        NCO = column_density_projection_CO(I) + local_NCO
        NH2 = column_density_projection_h2(I) + local_NH2
#ifdef OLD_CO_SELF_SHIELDING
        call CO_shielding_Lee(NCO, NH2, fshield_CO_ray)
#else
        call CO_shielding_Visser(NCO, NH2, fshield_CO_ray)
#endif
        fshield_CO = fshield_CO + (1d0 / nray) * fshield_CO_ray
      enddo
 
#endif /* SGCHEM */
      return
      end
c
c CO self-shielding, shielding of CO by H2, based on Lee et al (1996)
c
      subroutine CO_shielding_Lee(NCO, NH2, fshield_CO)
#ifdef SGCHEM
#ifdef OLD_CO_SELF_SHIELDING
      implicit none
#include "cool.h"
#include "shield_data_old.h"
      REAL NCO, NH2, fshield_CO
      REAL theta_CO, theta_H2, dth, dNact
      integer ipos

      if (NCO .le. NCO_shield_min) then
        theta_CO = 1d0
      elseif (NCO .ge. NCO_shield_max) then
        theta_CO = CO_self_shielding(imax_COss)
      else
        ipos  = int(dlog10(NCO) / dNshield) - ioff_co
        if (ipos .eq. 0) then
          theta_CO = 1d0
        else
          dth = CO_self_shielding(ipos+1) -
     $          CO_self_shielding(ipos)
          dNact = dlog10(NCO) - (ipos + ioff_co) * dNshield
          theta_CO = (CO_self_shielding(ipos)
     $               + dth * (dNact / dNshield))
        endif
      endif

      if (NH2 .le. NH2_shield_min) then
        theta_H2 = 1d0
      elseif (NH2 .ge. NH2_shield_max) then
        theta_H2 = CO_H2_shielding(imax_COH2)
      else
        ipos  = int(dlog10(NH2) / dNshield) - ioff_h2
        if (ipos .eq. 0) then
          theta_H2 = 1d0
        else
          dth   = CO_H2_shielding(ipos+1) -
     $            CO_H2_shielding(ipos)
          dNact = dlog10(NH2) - (ipos + ioff_h2) * dNshield
          theta_H2 = (CO_H2_shielding(ipos)
     $             + dth * (dNact / dNshield))
        endif
      endif
c
      fshield_CO = theta_H2 * theta_CO
c
#endif /* OLD_CO_SELF_SHIELDING */
#endif /* SGCHEM */
      return
      end
c
c CO self-shielding, shielding of CO by H2, based on Visser et al (2009)
c
      subroutine CO_shielding_Visser(NCO, NH2, fshield_CO)
#ifdef SGCHEM
#ifndef OLD_CO_SELF_SHIELDING
      implicit none
#include "cool.h"
#include "shield_data.h"
      REAL NCO, NH2, fshield_CO
      REAL factor_CO, factor_H2, fsh_1, fsh_2
      integer CO_index, H2_index, ipos, ipos2
c
      if (NCO .le. NCO_shield_min) then
        CO_index = 1
        if (NH2 .lt. NH2_shield_min) then
          fshield_CO = 1d0
        elseif (NH2 .eq. NH2_shield_min) then
          fshield_CO = CO_shield_factors(num_CO_self_shield+1)
        elseif (NH2 .ge. NH2_shield_max) then
          ipos = num_CO_self_shield * (num_CO_H2_shield - 1) + 1
          fshield_CO = CO_shield_factors(ipos)
        else
          H2_index = 2 + int((log10(NH2) - 15d0) / 0.2d0)
          factor_H2 = (log10(NH2) - log10(H2_shield_columns(H2_index)))
     $              / 0.2d0
          ipos  = num_CO_self_shield * (H2_index - 1) + 1
          ipos2 = num_CO_self_shield * (H2_index - 0) + 1
          fshield_CO = (1d0-factor_H2) * log10(CO_shield_factors(ipos))
     $               +      factor_H2  * log10(CO_shield_factors(ipos2))
          fshield_CO = 1d1**(fshield_CO)
        endif
        return
      elseif (NCO .ge. NCO_shield_max) then
        CO_index = num_CO_self_shield
        if (NH2 .le. NH2_shield_min) then
          print*, "Warning: high CO with no H2!"
          fshield_CO = CO_shield_factors(num_CO_self_shield)
        elseif (NH2 .ge. NH2_shield_max) then
          fshield_CO = CO_shield_factors(num_factors)
        else
          H2_index = 2 + int((log10(NH2) - 15d0) / 0.2d0)
          factor_H2 = (log10(NH2) - log10(H2_shield_columns(H2_index)))
     $              / 0.2d0
          ipos  = num_CO_self_shield * (H2_index - 0)
          ipos2 = num_CO_self_shield * (H2_index + 1)
          fshield_CO = (1d0-factor_H2) * log10(CO_shield_factors(ipos))
     $               +      factor_H2  * log10(CO_shield_factors(ipos2))
          fshield_CO = 1d1**(fshield_CO)
        endif
        return
      else
        CO_index = 2 + int((log10(NCO) - 10d0) / 0.2d0)
        factor_CO = (log10(NCO) - log10(CO_shield_columns(CO_index)))
     $            / 0.2d0
c
        if (NH2 .le. NH2_shield_min) then
          ipos = CO_index
          fshield_CO = (1d0-factor_CO)* log10(CO_shield_factors(ipos))
     $               +      factor_CO * log10(CO_shield_factors(ipos+1))
          fshield_CO = 1d1**(fshield_CO)
        elseif (NH2 .ge. NH2_shield_max) then
          ipos = num_CO_self_shield * (num_CO_H2_shield - 1) + CO_index
          fshield_CO = (1d0-factor_CO)* log10(CO_shield_factors(ipos))
     $               +      factor_CO * log10(CO_shield_factors(ipos+1))
          fshield_CO = 1d1**(fshield_CO)
        else
          H2_index = 2 + int((log10(NH2) - 15d0) / 0.2d0)
          factor_H2 = (log10(NH2) - log10(H2_shield_columns(H2_index)))
     $              / 0.2d0
c
          ipos  = num_CO_self_shield * (H2_index - 1) + CO_index
          ipos2 = num_CO_self_shield * (H2_index - 0) + CO_index
          fsh_1 = (1d0 - factor_H2) * log10(CO_shield_factors(ipos))
     $          +        factor_H2  * log10(CO_shield_factors(ipos2))
c
          ipos  = num_CO_self_shield * (H2_index - 1) + CO_index + 1
          ipos2 = num_CO_self_shield * (H2_index - 0) + CO_index + 1
          fsh_2 = (1d0 - factor_H2) * log10(CO_shield_factors(ipos))
     $          +        factor_H2  * log10(CO_shield_factors(ipos2))
c
          fshield_CO = (1d0 - factor_CO) * fsh_1 + factor_CO * fsh_2
          fshield_CO = 1d1**(fshield_CO)
        endif
        return
      endif
#endif /* OLD_CO_SELF_SHIELDING */
#endif /* SGCHEM */
      return
      end
c
c
c
      subroutine calc_C_shielding_factor(column_density_projection_h2,
     $                                   column_density_projection_c,
     $                                   local_NC, local_NH2, fshield_C)
#ifdef SGCHEM
      implicit none
#include "cool.h"
      REAL column_density_projection_h2(NPIX)
      REAL column_density_projection_c(NPIX)
      REAL local_NC, local_NH2, fshield_C
      REAL NC, NH2, CH2, fshield_C_ray
      integer nray, I
c
c C self-shielding, shielding by H2
c
      fshield_C = 0d0
      nray      = NPIX

      do I = 1, nray
        NC  = column_density_projection_c(I)  + local_NC
        NH2 = column_density_projection_h2(I) + local_NH2
        CH2 = 2.8d-22 * NH2
#ifdef SGCHEM_GONG17_TEST
c Incorrect version of shielding factor used in Gong et al tests
        fshield_C_ray = dexp(-1.6d-17 * NC) * exp(-CH2) * (1d0 + CH2)
#else
c Correct version per Tielens & Hollenbach (1985)
        fshield_C_ray = dexp(-1.6d-17 * NC) * exp(-CH2) / (1d0 + CH2)
#endif
        fshield_C = fshield_C + fshield_C_ray / nray
      enddo

#endif /* SGCHEM */
      return
      end
