#include "arepoconfig.h"
      program test
#ifdef SGCHEM
#define NSIDE 2
#define NPIX  12*NSIDE*NSIDE
#include "sgchem_def.h"
      implicit NONE
      real*8 abundD, abundc, abundo, abundM, Zatom, tdust, G0,
     $       LWBGz, D2G, crion, initz, h2ex, h2kin, AVeff, time,
     $       dl, yn, divv, energy, current_redshift, sink_flux
      real*8 non_eq_abundances(SGCHEM_NUM_SPECIES)
      real*8 column_density_projection(NPIX)
      real*8 column_density_projection_H2(NPIX)
      real*8 column_density_projection_CO(NPIX)
      real*8 thermal_rates_out(SGCHEM_NUM_THERMAL_RATES)
      real*8 photorates_direct(NDIRECT)

      integer LWBG, iphoto, isrf, iatom, id, I, iflag_h2_opacity

      abundD = 0.0
      abundc = 1.4d-4
      abundo = 3.2d-4
      abundM = 1.0d-7
      Zatom  = 1.0
      tdust  = 2d0
      G0     = 1.7d0
      LWBG   = 0
      LWBGz  = 0d0
      D2G    = 1.0
      crion  = 3d-17
      initz  = 0d0
      AVeff  = 0d0
      h2ex   = 0.84d0
      h2kin  = 0.12d0
      iphoto = 6
      isrf   = 1
      iatom  = 3
      iflag_h2_opacity = 0

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

      do I = 1, 8
        photorates_direct(I) = 0d0
      enddo

      open(141, file='column_data.txt', status='old')
      do I = 1, 48
        read(141,*) column_density_projection(I)
      enddo
      do I = 1, 48
        read(141,*) column_density_projection_H2(I)
      enddo
      do I = 1, 48
        read(141,*) column_density_projection_CO(I)
      enddo
      close(141, status='keep')

      time = 1d11
      dl   = 1.307343520685416d17
      yn   = 17977877.9196187
      divv = -4.600292570881704d-13
      energy   = 5.729631140329736d-8
      current_redshift = 0d0
c Current abundances
      non_eq_abundances(IH2+1) = 0.499999762024257
      non_eq_abundances(IHP+1) = 6.252429774686641d-13
      non_eq_abundances(ICO+1) = 1.399999999999999d-4
c Columns from Treecol
c      do I = 1, NPIX
c        column_density_projection(I) = 0d0
c        column_density_projection_H2(I) = 0d0
c        column_density_projection_CO(I) = 0d0
c      enddo
      sink_flux = 0d0
      id   = 154209690

      call evolve_abundances(time, dl, yn, divv, energy,
     $                       current_redshift, non_eq_abundances,
     $                       thermal_rates_out,
     $                       column_density_projection,
     $                       column_density_projection_H2,
     $                       column_density_projection_co,
     $                       photorates_direct,
     $                       tdust,sink_flux,id)

      print*, non_eq_abundances, energy
      print*, 0.666 * energy / (0.6 * yn * 1.38e-16), tdust

#endif
      end
