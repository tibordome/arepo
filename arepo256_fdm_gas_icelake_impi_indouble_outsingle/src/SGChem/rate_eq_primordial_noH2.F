#include "arepoconfig.h"
#ifdef SGCHEM_NO_MOLECULES
#if CHEMISTRYNETWORK == 1
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
c    Written by S. Glover (2004-2006)
c
c  Simplified version of the primordial chemistry network that doesn't include
c  any molecular ion or molecule chemistry. (Species such as H2 are still advected
c  but retain their initial abundances). Primarily designed for testing purposes
c
      subroutine rate_eq(neq, time, y, ydot, rpar, ipar)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "non_eq.h"
      integer neq
c
      REAL time, y(nspec), ydot(nspec), rpar(nrpar)
      integer ipar(nipar)
c
      integer niter
      parameter (niter = 25)
      REAL eqtol
      parameter (eqtol = 1d-3)
c
      REAL    dtemp
      integer itemp, I, J, K
c
      REAL abundances(nabn), rates(nrates),
     $     rates_chem(nrates_chem), ylam, ylam_chem,
     $     CMB_temp, ylam_cmb, rates_cmb(nrates)
      REAL temp, yn, dl, divv, energy, gamma, yntot,
     $     ekn, mp
      parameter (mp = 1.6726d-24)

      REAL abhp, abh2, abhd, abdp, abhep, abhepp
c
      REAL total_rad_heating, heating_rate, flux

      REAL abe, abHI, abDI, abheI

      REAL abM

      REAL ph4,  ph5, ph8, ph9, ph35, ph38, ph39,
     $     ph40, ph41, ph42, ph43, ph44, ph45
 
      REAL cr1, cr1a, cr2, cr6

      REAL phih, hvphih, uvphih, uvgeom,
     $     hephih, hevphih

      REAL ynsq, rho, cs2, L_jeans

      REAL  ch1,  ch2,  ch3,  ch4,  ch5,  ch6,  ch8,  ch9,
     $     ch10, ch11, ch13, ch14, ch15, ch17, ch18, ch20,
     $     ch21, ch22, ch23, ch24, ch25, ch26, ch27, ch28, 
     $     ch29, ch30, ch31, ch32, ch33, ch98, ch99, 
     $     ch100, ch101, ch102, ch103, ch108, ch109, ch110,
     $     ch126, ch127, ch137, ch138, ch139, ch144, ch148,
     $     ch153, ch154, ch155, ch156, ch157, ch158, 
     $     ch159, ch160, ch161, ch162, ch163, ch164, ch165, 
     $     ch166, ch167, ch168, ch169, ch170, ch171, ch172, 
     $     ch173

      REAL  dtch1,  dtch2, dtch3,   dtch4,  dtch5,  dtch6,
     $      dtch8,  dtch9, dtch10, dtch11, dtch13, dtch14, 
     $     dtch15, dtch17, dtch18, dtch20, dtch21, dtch22, 
     $     dtch23, dtch24, dtch25, dtch26, dtch27, dtch28, 
     $     dtch29, dtch30, dtch31, dtch32, dtch33,
     $     dtch98, dtch99,  dtch100, dtch101, dtch102, 
     $     dtch103, dtch108, dtch109, dtch110,
     $     dtch126, dtch127, dtch137, dtch138, dtch139,
     $     dtch144, dtch148, dtch153, dtch154, dtch155,
     $     dtch156, dtch157, dtch158, dtch159, 
     $     dtch160, dtch161, dtch162, dtch163, dtch164, 
     $     dtch165, dtch166, dtch167, dtch168, dtch169, 
     $     dtch170, dtch171, dtch172, dtch173
c
      no_dchem = 0
#ifdef SGCHEM_NO_HIGHN_DCHEM
      if (yn .gt. 1d8) then
        no_dchem = 1
      endif
#endif
c
c Unpack parameters
c
      yn         = rpar(1)
      dl         = rpar(2)
      divv       = rpar(3)
c
      flux       = rpar(11)     ! Accretion luminosity from sinks
      phih       = rpar(12)     ! HI ionization rate
      hvphih     = rpar(13)     ! Heating due to HI photoionization
      uvgeom     = rpar(17)     ! uvgeom: Used for photoelectric heating
      hephih     = rpar(18)     ! HeI ionization rate
      hevphih    = rpar(19)     ! Heating due to HeI photoionization
c
      ynsq  = yn**2
c
      call validate_input(y, ipar)
      if (ipar(1) .eq. 1) then
        do I = 1, nspec
          ydot(I) = 0d0
        enddo
c We simply return here because we'll be redoing the step anyway
c (making any further work pointless) and because we may be doing 
c this on account of having negative chemical abundances, which 
c would otherwise cause serious problems later on.
        return
      endif
c
c Here, we call lwbg to calculate the current LWBG rate G0_curr:
c    
      call lwbg
c
c From this point on, we know that our tracked abundances cannot be negative
c
      abhp   = y(ihp)
      abh2   = y(ih2)
      abhd   = y(ihd)
      abdp   = y(idp)
      abhep  = y(ihep)
      abhepp = y(ihepp)
c
      abe  = abhp + abdp + abhep + 2d0 * abhepp
c
      abHI = 1d0 - 2d0 * abH2 - abhp - abhd
      abHI = max(0d0, abHI)
c
      abDI  = max(0d0, abundD - abdp - abhd)
      if (yn .gt. 9d7 .and. no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        abDI = max(abDI, atol(ihd))
      endif
c
      abheI = max(0d0, abhe - abhep - abhepp)
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute current adiabatic index and total number density
c
      yntot = (1d0 + abhe + abe - abh2) * yn
      ekn = energy / (kboltz * yntot)

#ifdef VARIABLE_GAMMA
      call lookup_temp_gamma(abh2, ekn, temp, gamma)
#else
c If H2 fraction is small (as will always be the case in primordial gas
c at n << 1e8 cm^-3), then gamma = 5/3 to within a small fraction of
c a percent
      gamma = GAMMA
      temp = (gamma - 1d0) * ekn
#endif
c     
c Compute local Jeans length (for use in cool_func)
c
      call compute_Jeans_length(energy, yn, gamma, L_jeans)
c
c Compute CMB temperature
c
      CMB_temp = 2.726d0 * (1d0 + redshift)
c
c Bounds checking
c
      if (temp .le. 0d0) then
        itemp   = 1
        dtemp   = 0d0
        temp    = tmin
        ipar(1) = 1
      elseif (temp .le. tmin) then
        itemp = 1
        dtemp = 0d0
      elseif (temp .ge. tmax) then
        itemp = nmd
        dtemp = 0d0
      else
        itemp = int(dlog10(temp) / dtlog) + 1
        if (itemp .le. 0 .or. itemp .gt. nmd) then
          print*, 'Fatal error in rate_eq.F', itemp, temp
          ABORT(1)
        endif
        dtemp = temp - temptab(itemp)
      endif
c
c Rates:
c
       dtch1  = dtchtab(1, itemp)
       dtch2  = dtchtab(2, itemp)
       dtch3  = dtchtab(3, itemp)
       dtch4  = dtchtab(4, itemp)
       dtch5  = dtchtab(5, itemp)
       dtch6  = dtchtab(6, itemp)
c
       dtch8  = dtchtab(8, itemp)
       dtch9  = dtchtab(9, itemp)
c 
       dtch10 = dtchtab(10,itemp)
       dtch11 = dtchtab(11,itemp) 
c
       dtch13 = dtchtab(13,itemp)
       dtch14 = dtchtab(14,itemp)
       dtch15 = dtchtab(15,itemp)
c
       dtch17 = dtchtab(17,itemp)
       dtch18 = dtchtab(18,itemp)
c
       dtch20 = dtchtab(20,itemp)
       dtch21 = dtchtab(21,itemp)
       dtch22 = dtchtab(22,itemp)
       dtch23 = dtchtab(23,itemp)
       dtch24 = dtchtab(24,itemp)
       dtch25 = dtchtab(25,itemp)
       dtch26 = dtchtab(26,itemp)
       dtch27 = dtchtab(27,itemp)
       dtch28 = dtchtab(28,itemp)
       dtch29 = dtchtab(29,itemp)
       dtch30 = dtchtab(30,itemp)
       dtch31 = dtchtab(31,itemp)
       dtch32 = dtchtab(32,itemp)
       dtch33 = dtchtab(33,itemp)
c
       dtch98 = dtchtab(98,itemp)
       dtch99 = dtchtab(99,itemp)
       dtch100 = dtchtab(100,itemp)
       dtch101 = dtchtab(101,itemp)
       dtch102 = dtchtab(102,itemp)
       dtch103 = dtchtab(103,itemp)
c
       dtch108 = dtchtab(108,itemp)
       dtch109 = dtchtab(109,itemp)
       dtch110 = dtchtab(110,itemp)
c
       dtch126 = dtchtab(126,itemp)
       dtch127 = dtchtab(127,itemp)
c
       dtch137 = dtchtab(137,itemp)
       dtch138 = dtchtab(138,itemp)
       dtch139 = dtchtab(139,itemp)
c
       dtch144 = dtchtab(144,itemp)
c
       dtch148 = dtchtab(148,itemp)
c
       dtch153 = dtchtab(153,itemp)
       dtch154 = dtchtab(154,itemp)
       dtch155 = dtchtab(155,itemp)
       dtch156 = dtchtab(156,itemp)
       dtch157 = dtchtab(157,itemp)
       dtch158 = dtchtab(158,itemp)
       dtch159 = dtchtab(159,itemp)
       dtch160 = dtchtab(160,itemp)
       dtch161 = dtchtab(161,itemp)
       dtch162 = dtchtab(162,itemp)
       dtch163 = dtchtab(163,itemp)
       dtch164 = dtchtab(164,itemp)
       dtch165 = dtchtab(165,itemp)
       dtch166 = dtchtab(166,itemp)
       dtch167 = dtchtab(167,itemp)
       dtch168 = dtchtab(168,itemp)
       dtch169 = dtchtab(169,itemp)
       dtch170 = dtchtab(170,itemp)
       dtch171 = dtchtab(171,itemp)
       dtch172 = dtchtab(172,itemp)
       dtch173 = dtchtab(173,itemp)
c
       ch1  = chtab(1, itemp) + dtemp * dtch1
       ch2  = chtab(2, itemp) + dtemp * dtch2
       ch3  = chtab(3, itemp) + dtemp * dtch3
       ch4  = chtab(4, itemp) + dtemp * dtch4
       ch5  = chtab(5, itemp) + dtemp * dtch5
       ch6  = chtab(6, itemp) + dtemp * dtch6
c       
       ch8  = chtab(8, itemp) + dtemp * dtch8
       ch9  = chtab(9, itemp) + dtemp * dtch9
c
       ch10 = chtab(10,itemp) + dtemp * dtch10 
       ch11 = chtab(11,itemp) + dtemp * dtch11
c
       ch13 = chtab(13,itemp) + dtemp * dtch13
       ch14 = chtab(14,itemp) + dtemp * dtch14
       ch15 = chtab(15,itemp) + dtemp * dtch15
c
       ch17 = chtab(17,itemp) + dtemp * dtch17
       ch18 = chtab(18,itemp) + dtemp * dtch18
c
       ch20 = chtab(20,itemp) + dtemp * dtch20
       ch21 = chtab(21,itemp) + dtemp * dtch21
       ch22 = chtab(22,itemp) + dtemp * dtch22
       ch23 = chtab(23,itemp) + dtemp * dtch23
       ch24 = chtab(24,itemp) + dtemp * dtch24
       ch25 = chtab(25,itemp) + dtemp * dtch25
       ch26 = chtab(26,itemp) + dtemp * dtch26
       ch27 = chtab(27,itemp) + dtemp * dtch27
       ch28 = chtab(28,itemp) + dtemp * dtch28
       ch29 = chtab(29,itemp) + dtemp * dtch29
       ch30 = chtab(30,itemp) + dtemp * dtch30
       ch31 = chtab(31,itemp) + dtemp * dtch31
       ch32 = chtab(32,itemp) + dtemp * dtch32
       ch33 = chtab(33,itemp) + dtemp * dtch33
c
       ch98 = chtab(98,itemp) + dtemp * dtch98
       ch99 = chtab(99,itemp) + dtemp * dtch99
       ch100 = chtab(100,itemp) + dtemp * dtch100
       ch101 = chtab(101,itemp) + dtemp * dtch101
       ch102 = chtab(102,itemp) + dtemp * dtch102
       ch103 = chtab(103,itemp) + dtemp * dtch103
c
       ch108 = chtab(108,itemp) + dtemp * dtch108
       ch109 = chtab(109,itemp) + dtemp * dtch109
       ch110 = chtab(110,itemp) + dtemp * dtch110
c
       ch126 = chtab(126,itemp) + dtemp * dtch126
       ch127 = chtab(127,itemp) + dtemp * dtch127
c
       ch137 = chtab(137,itemp) + dtemp * dtch137
       ch138 = chtab(138,itemp) + dtemp * dtch138
       ch139 = chtab(139,itemp) + dtemp * dtch139
c
       ch144 = chtab(144,itemp) + dtemp * dtch144
c
       ch148 = chtab(148,itemp) + dtemp * dtch148      
c
       ch153 = chtab(153,itemp) + dtemp * dtch153
       ch154 = chtab(154,itemp) + dtemp * dtch154
       ch155 = chtab(155,itemp) + dtemp * dtch155
       ch156 = chtab(156,itemp) + dtemp * dtch156
       ch157 = chtab(157,itemp) + dtemp * dtch157
       ch158 = chtab(158,itemp) + dtemp * dtch158
       ch159 = chtab(159,itemp) + dtemp * dtch159
       ch160 = chtab(160,itemp) + dtemp * dtch160
       ch161 = chtab(161,itemp) + dtemp * dtch161
       ch162 = chtab(162,itemp) + dtemp * dtch162
       ch163 = chtab(163,itemp) + dtemp * dtch163
       ch164 = chtab(164,itemp) + dtemp * dtch164
       ch165 = chtab(165,itemp) + dtemp * dtch165
       ch166 = chtab(166,itemp) + dtemp * dtch166
       ch167 = chtab(167,itemp) + dtemp * dtch167
       ch168 = chtab(168,itemp) + dtemp * dtch168
       ch169 = chtab(169,itemp) + dtemp * dtch169
       ch170 = chtab(170,itemp) + dtemp * dtch170
       ch171 = chtab(171,itemp) + dtemp * dtch171
       ch172 = chtab(172,itemp) + dtemp * dtch172
       ch173 = chtab(173,itemp) + dtemp * dtch173
c
c H2 and HD: not updated, so set rates of change to zero
c
      ydot(ih2) = 0d0
      ydot(ihd) = 0d0
c
c  H ionization due to HeII recombination.
c  First, estimate whether we're optically thick to ionizing radiation; if not, then
c  we can ignore this correction
      
      NH = 0.5d0 * dl * yn
c
      if (2d0 * NH * abhI .lt. 1.6d17) then
         ch18 = 0d0
      else
c Linearly interpolate between low density, high density limits
        ch18 = ch18 + (ch148 - ch18) * ((yn / 4d3) / (1d0 + yn / 4d3))
      endif
c
c Protons
c
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          + ch144 * abHI  * abHI          ! Collisional ion: H + H
     $          + ch158 * abheI * abHI          ! Collisional ion: H + He
     $           - ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $           + ch18 * abhep * abe           ! Ion. from He+ rec: He+ + e
     $           + ch32 * abhep * abHI          ! Charge transfer: He+ + H
     $           - ch33 * abheI * abhp) * yn    ! Charge transfer: He  + H+
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          + ph38  * abhI                  ! Photoion: H + PH
     $          + phih  * abHI                  ! HI photoionization
c
      if (no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        ydot(ihp) = ydot(ihp) 
     $            + (ch99 * abdp  * abHI          ! Charge transfer: H + D+
     $            -  ch98 * abhp  * abDI) * yn    ! Charge transfer: D + H+)
      endif
c
c  Ionized deuterium (D+)
c
      if (abundD .le. 0d0 .or. no_dchem .eq. 1) then
        ydot(idp) = 0d0
      else
        ydot(idp) = (ch98   * abhp  * abDI          ! Charge transfer: D + H+
     $             - ch99   * abdp  * abHI          ! Charge transfer: H + D+
     $             - ch13   * abdp  * abe           ! Case B recomb.
     $             + ch11   * abDI  * abe           ! Collisional ionization
     $             + ch32   * abhep * abDI          ! Charge transfer: He+ + D
     $             - ch33   * abheI * abdp) * yn    ! Charge transfer: He  + D+
     $             + cr1a   * abDI                  ! Cosmic ray ion: D + CR
     $             + ph38   * abDI                  ! Photoionization
      endif
c
c Singly-ionized helium (He+)
c
      if (abhe .eq. 0d0) then
        ydot(ihep) = 0d0
      else
        ydot(ihep) = (ch14    * abheI * abe
     $              - ch15    * abhep * abe
     $              - ch17    * abhep * abe
     $              + ch20    * abhepp * abe
     $              - ch32    * abhep * abHI
     $              + ch33    * abheI * abhp) * yn
     $              + cr6     * abheI
     $              + ph39    * abheI
     $              - ph40    * abhep
     $              + hephih  * abheI  ! He photoionization
      endif
c
      if (no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        ydot(ihep) = ydot(ihep)
     $             + (ch33    * abheI * abdp
     $             -  ch32    * abhep * abDI) * yn
      endif
c
c Doubly-ionized helium (He++)
c
      if (abhe .eq. 0d0) then
        ydot(ihepp)  = 0d0
      else
        ydot(ihepp) = (ch15 * abhep  * abe
     $              -  ch20 * abhepp * abe) * yn
     $              +  ph40 * abhep
      endif
c
c Chemical heating & cooling
c
      do I = 1, nrates_chem
        rates_chem(I) = 0d0
      enddo
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(8) = 13.6d0 * eV * (ch11 * abe + ch144 * abHI
     $              + ch158 * abheI) * abHI * ynsq
c
c HeI collisional ionization cooling -- 24.6 eV per ionization
c
      rates_chem(9) = 24.6d0 * eV * ch14 * abe * abheI * ynsq
c
c HeII collisional ionization cooling -- 54.4 eV per ionization
c
      rates_chem(10) = 54.4d0 * eV * ch15 * abe * abhep * ynsq
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(11) = kboltz * temp * abhp * ch13 * abe * 
     $                 ynsq
c
c HeII recombination cooling. Dielectronic recombination results
c in the loss of ~ 40.8eV per recombination, while for standard 
c radiative recombination, we make the same assumption as for HII
c
c Dielectronic:
      rates_chem(12) = 40.8d0 * eV * ch10 * abhep * abe * ynsq
c
c Radiative:
      rates_chem(13) = kboltz * temp * (ch17 - ch10) * abhep * abe * 
     $                 ynsq
c
c HeIII recombination cooling -- see note for HII
c
      rates_chem(14) = kboltz * temp * abhepp * ch20 * abe * 
     $                 ynsq
c
c HI photoionization heating
c 
      rates_chem(17) = - ph42 * abHI * yn
c
c HeI photoionization heating
c 
      rates_chem(18) = - ph43 * abheI * yn
c
c HeII photoionization heating
c 
      rates_chem(19) = - ph44 * abhepp * yn
c
      ylam_chem = 0d0
      do I = 1, nrates_chem
        ylam_chem = ylam_chem + rates_chem(I)
      enddo
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0d0
      enddo
c
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
      abundances(14) = abheI
      abundances(15) = abhep
      abundances(16) = abhepp
c
#ifdef ADIABATIC_DENSITY_THRESHOLD
      if (yn .ge. yn_adiabatic) then
        ydot(itmp) = 0d0
#ifdef DEBUG_RATE_EQ
        call validate_output(y,ydot,rpar)
#endif
        return
      endif
#endif
c
      call cool_func(temp, yn, L_jeans, divv, G0_curr, abundances, ylam,
     $               rates, 0, energy)
c
      if(ylam .ne. ylam) then
       print*,'RATE_EQ2: ylam NaN:', temp, yn, L_jeans, divv,
     $        'abund:', abundances, ylam,
     $        'rates:', rates,'energy:', energy
      endif
c
c Now correct for effects of CMB:
c
      if (redshift .gt. 0d0) then
        call cool_func(CMB_temp, yn, L_jeans, divv, G0_curr, abundances,
     $                 ylam_cmb, rates_cmb, 1, energy)
        if (ylam_cmb .ne. ylam_cmb) then
          print*,'RATE_EQ3: ylam_cmb NaN:', CMB_temp, yn, L_jeans, divv,
     $           'abund:', abundances,'ylam_cmb:', ylam_cmb,
     $           'rates:', rates_cmb,'energy:', energy
        endif
      else
        ylam_cmb = 0d0
      endif
c
c Finally, account for accretion luminosity heating. Note that since the
c effectiveness depends on the local temperature, we need to do this here
c rather than in sgchem.c or evolve_abundances.F
c
      total_rad_heating = 0.0
#ifdef SGCHEM_ACCRETION_LUMINOSITY
      rho   = (1d0 + 4d0 * abhe) * mp * yn
      call compute_heating(flux, temp, rho, heating_rate)
      total_rad_heating = heating_rate
#endif


c Include HI photoionization heating here; other photoheating terms already accounted for abov
      ydot(itmp) = -(ylam - ylam_cmb + ylam_chem) + total_rad_heating
     $             + hvphih * abHI * yn  ! HI photoionization heating
     $             + hevphih * abHeI * yn ! HeI photoionization heating

#ifdef THERMAL_INFO_DUMP
      dtcool_nopdv = rates(29) !hitchhiked H2_opacity_correction
      cool_h2_line = rates(2) - rates_cmb(2)
      cool_h2_cie  = rates(20) - rates_cmb(20)
      cool_h2_diss = rates_chem(1)
      heat_3b      = rates_chem(5) + rates_chem(6) + rates_chem(7)
#endif
c
#ifdef DEBUG_RATE_EQ
      if (ydot(itmp) .ne. ydot(itmp)) then
        print*, "Cooling rate is NaN: y=",y,' rpar=',rpar,' ylam=',ylam,
     $         ' ylam_cmb=',ylam_cmb,' ylam_chem=',ylam_chem,' rates='
     $        ,rates,' rates_chem=',rates_chem
#if SGCHEM_ACCRETION_LUMINOSITY
     $        ,total_rad_heating
#endif
      endif
      call validate_output(y,ydot,rpar)
#endif
c
#endif /* SGCHEM */
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              R A T E _ E Q                \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c    \\\\\\\\\\        B E G I N   F U N C T I O N        //////////
c    //////////          C O M P _ E Q U I L I B          \\\\\\\\\\
c
c=======================================================================
c
#ifdef SGCHEM
      REAL function comp_equilib(num, denom, max_abundance, id, yn, 
     $                           time_to_equilibrium)
      implicit NONE
#include "cool.h"
#include "non_eq.h"
       REAL num, denom, max_abundance, yn
       REAL time_to_equilibrium(num_eqb_species)
       integer id
c
       if (num .eq. 0d0) then
         comp_equilib = 0d0
       elseif (denom .eq. 0d0) then
c If this happens, then our assumption that this species is in
c equilibrium has fatally broken down...
         print*, "Creation without destruction in comp_equilib"
         print*, "Species = ", id
c         ABORT(9751)
         comp_equilib = -1d0
      else
#ifdef COLLECT_EQB_DATA
         time_to_equilibrium(id) = 1d0 / (denom * yn)
#endif
        comp_equilib = num / denom
        if (comp_equilib .gt. max_abundance) then
          print*, "COMP_EQUI Error: impossible abundance!"
          print*, "Species ", id
          print*, num, denom, num/denom
c          ABORT(9752)
          comp_equilib = -1d0
        endif
      endif
c
      return
      end
#endif /* SGCHEM */
c=======================================================================
c
c    \\\\\\\\\\          E N D     F U N C T I O N        //////////
c    //////////          C O M P _ E Q U I L I B          \\\\\\\\\\
c
c=======================================================================
#endif /* CHEMISTRYNETWORK == 1 */
#endif /* SGCHEM_NO_MOLECULES */
