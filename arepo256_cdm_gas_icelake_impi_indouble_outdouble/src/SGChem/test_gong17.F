#include "arepoconfig.h"
      program test
#ifdef SGCHEM
#define NPIX  1
#include "sgchem_def.h"
      implicit NONE
      real*8 abundD, abundc, abundo, abundM, Zatom, tdust, G0,
     $       LWBGz, D2G, crion, initz, h2ex, h2kin, AVeff, time,
     $       dl, yn, divv, energy, current_redshift, sink_flux
      real*8 non_eq_abundances(SGCHEM_NUM_SPECIES)
      real*8 column_density_projection(NPIX)
      real*8 column_density_projection_H2(NPIX)
      real*8 column_density_projection_CO(NPIX)
      real*8 column_density_projection_C(NPIX)
      real*8 thermal_rates_out(SGCHEM_NUM_THERMAL_RATES)
      real*8 photorates_direct(NDIRECT)
      real*8 column(2000), columnH2(2000), columnCO(2000)
      real*8 DN(2000), columnC(2000)
      real*8 ntot, abe, temp, abundCI
      real*8 column_factor, NC
 
      integer LWBG, iphoto, isrf, iatom, id, I, iflag_h2_opacity, J, K

      abundD = 2.6d-5
      abundc = 1.6d-4
      abundo = 3.2d-4
      abundM = 1.7d-6
      Zatom  = 1.0
      tdust  = 2d1
      G0     = 1.7d0 / 2d0
      LWBG   = 1
      LWBGz  = 0d0
      D2G    = 1.0
      crion  = 2d-16
      initz  = 0d0
      AVeff  = 0d0
      h2ex   = 0.84d0
      h2kin  = 0.12d0
      iphoto = 6
      isrf   = 1
      iatom  = 3
      iflag_h2_opacity = 0

      column(1) = 1d17
      column_factor = 1d1**(2.5d-3)

      columnH2(1) = 0d0
      columnCO(1) = 0d0
      columnC(1)  = 0d0

      do J = 2, 2000
        column(J) = column(J-1) * column_factor
        columnH2(J) = 0d0
        columnCO(J) = 0d0
        columnC(J)  = 0d0
      enddo

      do J  = 1, 1999
        DN(J) = column(J+1) - column(J)
      enddo
      DN(2000) = 0.0

      call init_chemistry_parameters(abundD,
#ifndef SGCHEM_VARIABLE_Z
     $ abundc, abundo, abundM, Zatom,
#endif
     $ tdust, G0, LWBG, LWBGz,
#ifndef SGCHEM_VARIABLE_Z
     $ D2G,
#endif
#ifndef SGCHEM_VARIABLE_CRION
     $ crion,
#endif
     $ initz, AVeff, h2ex, h2kin, iphoto, isrf, iatom, iflag_h2_opacity)

      call coolinmo
      call cheminmo
      call init_tolerances

#ifdef SGCHEM_VARIABLE_Z
      call set_local_dust_abundance(D2G)
      call set_local_element_abundances(abundc, abundo, abundM,
     $                                  Zatom)
#endif

      divv = 1d-15

      yn   = 98.679
      current_redshift = 0d0

c Loop over column densities
      do J = 1, 2000

        dl = 0.0
c        print*, dl, DN(J)

c Set starting abundances to zero, except for H+, C+, He+
        do I = 1, SGCHEM_NUM_SPECIES
          non_eq_abundances(I) = 0.0
        enddo
c Set H+ and He+ abundances to small but non-zero values
        non_eq_abundances(IHP+1)  = 1d-4
        non_eq_abundances(IHEP+1) = 1d-5
c
c Assume that carbon starts as C+
c 
        non_eq_abundances(ICP+1) = abundc

        abe = non_eq_abundances(IHP+1) + non_eq_abundances(IHEP+1) 
     $      + non_eq_abundances(ICP+1)

        ntot = (1.1 + abe) * yn
c Fix temperature to 20 K, as in Gong et al. At this temperature, can safely assume gamma=5/3
        temp = 2d1
        energy = 1.5d0 * ntot * 1.38066e-16 * temp

c Columns from Treecol
        do I = 1, NPIX
          column_density_projection(I)    = 2.0 * column(J)
          column_density_projection_H2(I) = 2.0 * columnH2(J)
          column_density_projection_CO(I) = 2.0 * columnCO(J)
          column_density_projection_C(I)  = 2.0 * columnC(J)
        enddo
c Photorates from SPRAI
        do I = 1, NDIRECT
          photorates_direct(I) = 0d0
        enddo
        sink_flux = 0d0
        id        = 1

        abundCI = 0d0

        do K = 1, 300
          time = 3.16d13

          call evolve_abundances(time, dl, yn, divv, energy,
     $                           current_redshift, non_eq_abundances,
     $                           thermal_rates_out,
     $                           column_density_projection,
     $                           column_density_projection_H2,
     $                           column_density_projection_co,
     $                           column_density_projection_c,
     $                           photorates_direct,
     $                           tdust,sink_flux,id)
        enddo

        abundCI = abundc - non_eq_abundances(ICO+1) 
     $          - non_eq_abundances(ICP+1) 
     $          - non_eq_abundances(ICHX+1)
     $          - non_eq_abundances(IHCOP+1)


        print*, column(J), columnH2(J), non_eq_abundances 

        if (J .lt. 2000) then
        columnH2(J+1) = columnH2(J) + DN(J) 
     $                * non_eq_abundances(IH2+1) 
        columnCO(J+1) = columnCO(J) + DN(J) 
     $                * non_eq_abundances(ICO+1)
        columnC(J+1)  = columnC(J)  + DN(J) 
     $                * abundCI
        endif
      enddo


#endif
      end
