/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Chemistry that calculates rages for the SGChem 
 * \details     
 */

#include <stdio.h>

#include "../allvars.h"
#include "../proto.h"

#if(SX_CHEMISTRY == 3)

static int debugID = 1016795866;

static double uvPumpFact = 6.94;   // UV Pumping factor (Draine & Bertoldi 1996)
static double minMassFrac = 1e-20; // minimum of the mass fraction to prevent NaN 
static double recFactThr = 0.0;    // threshold value that triggers the recombination correction
static double recFactCutout = 0.1; // minimum recombination factor that is worth to calcualte
                                   // 0.1 is an empirical value for which it seems to work

//======================================================================================================//
//                           Functions used by other parts of the code                                  //
//======================================================================================================//
/*!
 * \brief functions that returns Arepo image property value
 */
double sx_chem_image1(int index)
{
  return SphP[index].sxPhotonRates[RIH];
}
double sx_chem_image2(int index)
{
  return SphP[index].sxPhotonRates[HRIH];
}
double sx_chem_image_all(int index, int prop)
{
  return SphP[index].sxPhotonRates[prop];
}
void sx_chem_cell_props( int index, double flux[SX_NRATES], double *dens, double *utherm )
{
  *dens = (double)SphP[index].Density;
  *utherm = (double)SphP[index].Utherm;
  for (int r=0; r<SX_NRATES; r++ )
    flux[r] = (double)SphP[index].sxPhotonRates[r];

}
//======================================================================================================//
//                                    Chemistry initialization                                          //
//======================================================================================================//

/*!
 * \brief chemistry initialization at the beginning of the simulation
 */
void sx_chem_initialize(void)
{
  // initialize frequency bins [Hz]
  All.sxFreqBins[0] = TFHDUST;  // 5.6  eV - photo electric heating on dust
  All.sxFreqBins[1] = TFDH2;    // 11.2 eV - H2 dissociation by UV
  All.sxFreqBins[2] = TFIH;     // 13.6 eV - H ionisation + H2 dissociation by UV
  All.sxFreqBins[3] = TFIH2;    // 15.2 eV - H2, H ionisation
  All.sxFreqBins[4] = TFIHE;    // 24.6 eV - He, H, H2 ionisation
  All.sxFreqBins[5] = 10.*TFIH; // 136  eV - infinity

#if (SX_SOURCES>10)
  // Include source settings and set the values
#include "sx_test_sources.h"
  for(int f=0; f<SX_NSIGMA; f++)
    All.sxSigma[f] = sxTS_sigma[f];
  for(int f=0; f<SX_NENERGY; f++)
   All.sxEnergy[f] = sxTS_energy[f];
#else
  // DEBUG: cross sections and photon energies should be in long runs calculated at the beginning of each run
  // DEBUG: for now we will assume a black body with temperature set by All.sxTeff
  // initialize cross sections
  sx_chem_compute_cross_sections(All.sxTeff);
  // initialize mean photon energies
  sx_chem_compute_mean_photon_energies(All.sxTeff);
#endif

#if(SX_SOURCES == 5)
  // save the initial mass of all sink particles
  for(int f = 0; f < NSinksAllTasks; f++)
    {
      SinkP[f].MassOld = SinkP[f].Mass;
      //printf("SX: (%d) sink %d Mass %.03e \n",ThisTask,f,SinkP[f].Mass);
    }
#endif
}

/*!
 * \brief initialize chemistry at the beginning of the run
 */
void sx_chem_start_run(void)
{
  // loop over active particles
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      // set variables to zero
      memset(SphP[i].sxPhotonRates, 0, SX_NRATES * sizeof(double));
#ifdef SX_RADIATION_PRESSURE
      memset(sxCell[i].MomentumChange, 0, SX_NDIM * sizeof(double));
#endif
      
      // make an image of chemical abundances that will be used during the RT step
      for(int s = 0; s < SGCHEM_NUM_ADVECTED_SPECIES; s++)
	{
	  sxCell[i].TracAbund[s] = SphP[i].TracAbund[s];
	}      
      double meanweight = (SphP[i].TracAbund[IHATOM] + 2 * SphP[i].TracAbund[IH2] + SphP[i].TracAbund[IHP] 
			   + 4 * SphP[i].TracAbund[IHEATOM] + 4 *  SphP[i].TracAbund[IHEP]) /
	                  ( SphP[i].TracAbund[IHATOM] +  SphP[i].TracAbund[IH2] + 2 *  SphP[i].TracAbund[IHP] 
			    +  SphP[i].TracAbund[IHEATOM] + 2 *  SphP[i].TracAbund[IHEP]);
      double temp = SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * 2.0/3.0 * meanweight * PROTONMASS / BOLTZMANN; // K
#ifdef SX_CMB_TEMP_LIMIT
      // restricting the minimum temperature to the CMB threshold
      if (temp < sxRunData->tempCMB)
	{
	  //printf("SX: (%d:%d:%d) Particle with temperature %.1f hit the CMB temperature floor of %.1f \n",
	  //	 ThisTask, i, P[i].ID, temp, sxRunData->tempCMB);
	  temp = sxRunData->tempCMB;
	  SphP[i].Utherm = temp * ( All.UnitMass_in_g / All.UnitEnergy_in_cgs ) 
	    * BOLTZMANN / ( 2.0/3.0 * meanweight * PROTONMASS );
	}
#endif
#ifdef SX_RECOMBINE
      // Assumed GAMMA = 5/3
      if (temp < 6000.0)
        temp = 6000.0;
      double tinv = 1.0 / temp;
      double log_T = log10(temp);
      double density =  SphP[i].Density * sxRunData->density2cgs;    // density [g/cm^3]
      // NOTE: The 'numdens' assumes everything is fully ionized and so the number density is electron number density
      double numdens = density / ((1. + 4. * ABHE) * PROTONMASS);    // nucleon number density [1/cm^3]

      // Recombination rate coefficients from SGChem/cheminmo.F
      // The hydrogen case B rate is from Ferland et al (1992)
      double hII_rec_B = 2.753e-14 * pow(315614e0 * tinv, 1.500e0) /     // cm^3/s
        (pow(1e0 + pow(115188e0 * tinv, 0.407e0),2.242e0));

      // The helium case A and case B rates are from Hummer & Storey (1998)
      double heII_rec_A  = 1e-11 * (12.72e0 - 1.615e0 * log_T            // cm^3/s == 1/(s * 1/cm^3)
          - 0.3162e0 * log_T * log_T
          + 4.93e-2 * log_T * log_T * log_T) / sqrt(temp);

      double heII_rec_B  = 1e-11 * (11.19e0 - 1.676e0 * log_T            // cm^3/s == 1/(s * 1/cm^3)
          - 0.2852e0 * log_T * log_T
          + 4.433e-2 * log_T * log_T * log_T) / sqrt(temp);

      // Recombination rate to the He ground state, which by definition is just the difference between the
      // case A and case B rates
      double heII_rec_gs = heII_rec_A - heII_rec_B;

      // Fraction of recombinations to excited states of He that produce photons able to ionize H,
      // computed following Osterbrock & Ferland (2006). We assume a fixed critical density of 4e3 cm^-3
      // for the He 2^3S triplet state, corresponding to the value at T ~ 1e4 K. The real critical density
      // varies from ~6e3 cm^-3 at 6000 K to ~3.3e3 cm^-3 at T = 2.5e4 K, so in the worst case we have an
      // error of ~5% in the value of p

      double pfac = (1.0 / 6.0) + 0.56 / 12.0 + (0.597 + 0.403 / (1 + numdens / 4e3)) * 0.75;

      // Fraction of photons produced by recombination to He ground state that are absorbed by neutral H.
      // The remainder are absorbed by neutral He. Note that this calculation assumes that the ratio of
      // neutral He to neutral H is the same as the cosmological He/H ratio. This is only approximately
      // true in HII regions: the real ratio will depend on the gas temperature, the source spectrum and
      // also the gas density (which controls the effectiveness of charge transfer between H and He).

      double yfac = 0.677 - (0.054 * (ABHE - 0.079) / 0.021);

      double heII_rec_actual = yfac * heII_rec_gs + heII_rec_B;
      double hII_rec_actual  = hII_rec_B - (pfac * heii_rec_B + yfac * heii_rec_gs) * ABHE;

      sxCell[i].H_rec_fac  = sxRunData->runTime * numdens * hII_rec_actual;
      sxCell[i].He_rec_fac = sxRunData->runTime * numdens * heII_rec_actual;
      if (sxCell[i].H_rec_fac < recFactCutout)
        sxCell[i].H_rec_fac = 0.0;
      if (sxCell[i].He_rec_fac < recFactCutout)
        sxCell[i].He_rec_fac = 0.0;
#endif /* SX_RECOMBINE */
    }
}

/*!
 * \brief Initialize chemistry variables in the site structure
 */
void sx_chem_initialize_site(struct sxSite *site)
{
  // set all rates to zero
  memset(site->rates, 0, SX_NRATES * sizeof(double));
  // initial abundances from SGChem [number of nucleons]
  site->initxH   = SphP[site->index].TracAbund[IHATOM];
  site->initxH2  = SphP[site->index].TracAbund[IH2];
  site->initxHp  = SphP[site->index].TracAbund[IHP];
  site->initxHe  = SphP[site->index].TracAbund[IHEATOM];
  site->initxHep = SphP[site->index].TracAbund[IHEP];
  // need something to calculate otherwise enjoy some NaNs
  if (site->initxH == 0.0 || site->initxHp == 1.0)  site->initxH = minMassFrac;
  if (site->initxH <= minMassFrac)  site->initxH = minMassFrac;
  if (site->initxH2 <= minMassFrac) site->initxH2 = minMassFrac;
  if (site->initxHe <= minMassFrac) site->initxHe = minMassFrac;

  // abundances evolving during the RT [number of nucleons]
  site->xH = sxCell[site->index].TracAbund[IHATOM];
  site->xH2 = sxCell[site->index].TracAbund[IH2];
  site->xHp = sxCell[site->index].TracAbund[IHP];
  site->xHe = sxCell[site->index].TracAbund[IHEATOM];
  site->xHep = sxCell[site->index].TracAbund[IHEP];
  // need something to calculate otherwise enjoy some NaNs
  if (site->xH <= minMassFrac) site->xH = minMassFrac;
  if (site->xH2 <= minMassFrac) site->xH2 = minMassFrac;
  if (site->xHe <= minMassFrac) site->xHe = minMassFrac;

#ifdef SX_RECOMBINE
  site->H_rec_fac  = sxCell[site->index].H_rec_fac;
  site->He_rec_fac = sxCell[site->index].He_rec_fac;  
#endif
#ifdef SX_RADIATION_PRESSURE
  site->momentum = 0.0;
#endif

  site->density =  SphP[site->index].Density * sxRunData->density2cgs;               // density [g/cm^3]
  site->volume = SphP[site->index].Volume * sxRunData->volume2cgs;                   // volume [cm^3]
  site->meanNgbDist = sxCell[site->index].MeanNeighDistance * sxRunData->length2cgs; // distance [cm] 
  site->numdens   = site->density / ((1. + 4. * ABHE) * PROTONMASS);  // nucleon number density [1/cm^3]
  site->nNucleons = site->numdens * site->volume;                     // initial number of nucleon species

  // nucleon number density
  site->numdr = site->numdens * site->meanNgbDist;

  if (site->isSource) 
    {
#ifdef SX_RADIATION_PRESSURE
      // We ban the source cell from stealing all of the radiation pressure.
      // Anyway, it is usually smaller than the sink radius
      site->numdr = 0.0;  
#else
      // in case of sources we count just half the distance
      // DEBUG: notice that the number of photons in the cell is the combination of source and incoming photons
      //        in this case the path of the incoming photons is also only half!!! we neglect it... a bit...      
      site->numdr *= 0.5;
#endif
    }

  //if (site->isSource)
  //  printf("SX: (%d:%d:%d) Rot %d Hrf %.03e Herf %.03e xHp %.03e xHep %.03e \n",ThisTask,site->index,P[site->index].ID,
  //	   sxRunData->currentRot+1, site->H_rec_fac, site->He_rec_fac, site->initxHp, site->initxHep);
}

/*!
 * \brief Finalize chemistry variables in the site structure
 */
void sx_chem_finalize_site( struct sxSite *site )
{

  // update rates
  for (int i=0; i<SX_NRATES; i++ )
      SphP[site->index].sxPhotonRates[i] += site->rates[i];

  // update abundances
  sxCell[site->index].TracAbund[IHATOM] = site->xH;
  sxCell[site->index].TracAbund[IHEATOM] = site->xHe;
  sxCell[site->index].TracAbund[IH2] = site->xH2;
#ifdef SX_RECOMBINE
  sxCell[site->index].H_rec_fac = site->H_rec_fac;
  sxCell[site->index].He_rec_fac = site->He_rec_fac;
#endif

#ifdef SX_RADIATION_PRESSURE
  // calculate the momentum vector
  double momVector[SX_NDIM] = {0.,0.,0.};
  for(int d=0; d<SX_NDIR; d++ )
    {
      int app = sxSPDtoAPP[ site->index * SX_NDIR + d ];
      if ( app >= 0 )
        {
          double nPhotons = 0.;
          for(int f=0; f<SX_NFREQ; f++ )
            {
              nPhotons += (double)sxAPP[app].nPhotons[f] * (1.-site->nPhotonsRatio[f]);
            }
          momVector[0] += sxRunData->orthoBase[d][2][0] * nPhotons;
          momVector[1] += sxRunData->orthoBase[d][2][1] * nPhotons;
          momVector[2] += sxRunData->orthoBase[d][2][2] * nPhotons;
        }
    }
  // add momentum to the cell
  double unitNorm = sqrt( momVector[0]*momVector[0] + momVector[1]*momVector[1] + momVector[2]*momVector[2] );
  if (unitNorm>0.0)
    {
      sxCell[site->index].MomentumChange[0] += momVector[0]/unitNorm * site->momentum;
      sxCell[site->index].MomentumChange[1] += momVector[1]/unitNorm * site->momentum;
      sxCell[site->index].MomentumChange[2] += momVector[2]/unitNorm * site->momentum;
    }
#endif
}

/*!
 * \brief finalize chemistry at the end of the run
 */
void sx_chem_end_run(void)
{
  // loop over active particles
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef SX_RADIATION_PRESSURE
      // normalize the momentum vector and add momentum to the cell
      if ( sxCell[i].MomentumChange[0]!=0. || sxCell[i].MomentumChange[1]!=0. || sxCell[i].MomentumChange[2]!=0.)
	{
	  double momFactor = 1. / ( sxRunData->mass2cgs * sxRunData->velocity2cgs ); // * sxRunData->runTimeInv;
          double abs_mom2_old = pow(SphP[i].Momentum[0],2) + pow(SphP[i].Momentum[1],2) + pow(SphP[i].Momentum[2],2);
	  SphP[i].Momentum[0] += sxCell[i].MomentumChange[0] * momFactor;
	  SphP[i].Momentum[1] += sxCell[i].MomentumChange[1] * momFactor;
	  SphP[i].Momentum[2] += sxCell[i].MomentumChange[2] * momFactor;
	  P[i].Vel[0] = SphP[i].Momentum[0] / P[i].Mass;
	  P[i].Vel[1] = SphP[i].Momentum[1] / P[i].Mass;
	  P[i].Vel[2] = SphP[i].Momentum[2] / P[i].Mass;
	  double abs_mom2 = pow(SphP[i].Momentum[0],2) + pow(SphP[i].Momentum[1],2) + pow(SphP[i].Momentum[2],2);
	  // XXX: does this work in comoving coordinates?
	  SphP[i].Energy += 0.5 * (abs_mom2 - abs_mom2_old) / P[i].Mass;
	}
#endif

    }
}
//======================================================================================================//
//                                     Rate equation                                                    //
//======================================================================================================//

/*
 * \brief Solve H, He and H2 ionization at frequencies above 24.6 eV
 */
void sx_chem_solve_F246( struct sxSite *site )
{
  // calculate number of available species
  double nSpecH = site->nNucleons * site->xH;
  double nSpecH2 = site->nNucleons * site->xH2;
  double nSpecHe = site->nNucleons * site->xHe;
#ifdef SX_RECOMBINE
  // correcting number of absorbed photons for recombination
  if (site->He_rec_fac > recFactThr)
    nSpecHe += site->nNucleons * ABHE * site->He_rec_fac;
  if (site->H_rec_fac > recFactThr)
    nSpecH += site->nNucleons * site->H_rec_fac;
#endif
  double nSpec = nSpecH + nSpecH2 +nSpecHe;
  if (nSpec==0) return;  // return if there are no available species for a reaction
  
  // some initial constants
  double sigmaH = All.sxSigma[S246IH] * site->xH;
  double sigmaH2 = All.sxSigma[S246IH2] * site->xH2;
  double sigmaHe = All.sxSigma[S246IHE] * site->xHe;
  double sigma = sigmaH + sigmaH2 + sigmaHe;
  if (sigma==0) return;  // skip this frequency if there is not a cross-section for the reaction

  // caclulate number of attenuated photons
  double nPhot = site->nPhotons[F246] * (1.0-exp( -1.0 * site->numdr * sigma ));
  if (nPhot==0) return;
  double nPhotH = nPhot * sigmaH / sigma;
  double nPhotH2 = nPhot * sigmaH2 / sigma;
  double nPhotHe = nPhot * sigmaHe / sigma;
  if (nPhotH>nSpecH) nPhotH=nSpecH;
  if (nPhotH2>nSpecH2) nPhotH2=nSpecH2;
  if (nPhotHe>nSpecHe) nPhotHe=nSpecHe;
  nPhot = nPhotH + nPhotH2 + nPhotHe;

#ifdef SX_RADIATION_PRESSURE
  site->momentum += nPhotH  * ( All.sxFreqBins[F246] * PLANCK + All.sxEnergy[E246IH]  ) * sxRunData->cInv; // erg/c
  site->momentum += nPhotH2 * ( All.sxFreqBins[F246] * PLANCK + All.sxEnergy[E246IH2] ) * sxRunData->cInv; // erg/c
  site->momentum += nPhotHe * ( All.sxFreqBins[F246] * PLANCK + All.sxEnergy[E246IHE] ) * sxRunData->cInv; // erg/c
#endif
    
  double attRatio;
  
 // update rates and mass fractions of H
  attRatio = nPhotH / (site->nNucleons * site->initxH);
  double maxAttRatio = site->nPhotons[F246] * All.sxSigma[S246IH] * site->meanNgbDist / site -> volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RIH] += attRatio * sxRunData->runTimeInv;
  site->rates[HRIH] += attRatio * sxRunData->runTimeInv * All.sxEnergy[E246IH];
#ifdef SX_RECOMBINE
  double red_factor = (1.0 - nPhotH / nSpecH);
  if (!(red_factor > 1.0e-6)) red_factor = 0.0;
  site->H_rec_fac *= red_factor; 
  site->xH *= red_factor;
#else
  site->xH -= attRatio * site->initxH;
  if (site->xH < 0.0) site->xH = minMassFrac; // preventing negative values from rounding
#endif
  
  // update rates and mass fractions of H2
  attRatio = nPhotH2 / (site->nNucleons * site->initxH2);
  maxAttRatio = site->nPhotons[F246] * All.sxSigma[S246IH2] * site->meanNgbDist / site->volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RIH2] += attRatio * sxRunData->runTimeInv;
  site->rates[HRIH2] += attRatio * sxRunData->runTimeInv * All.sxEnergy[E246IH2];
  site->xH2 -= attRatio * site->initxH2;
  if (site->xH2 < 0.0) site->xH2 = minMassFrac; // preventing negative values from rounding
  // update rates and mass fractions of He
  attRatio = nPhotHe / (site->nNucleons * site->initxHe);
  maxAttRatio = site->nPhotons[F246] * All.sxSigma[S246IHE] * site->meanNgbDist / site -> volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RIHE] += attRatio * sxRunData->runTimeInv;
  site->rates[HRIHE] += attRatio * sxRunData->runTimeInv * All.sxEnergy[E246IHE];
#ifdef SX_RECOMBINE
  red_factor = (1.0 - nPhotHe / nSpecHe);
  if (!(red_factor > 1.0e-6)) red_factor = 0.0;
  site->He_rec_fac *= red_factor; 
  site->xHe *= red_factor;
#else
  site->xHe -= attRatio * site->initxHe;
  if (site->xHe < 0.0) site->xHe = minMassFrac; // preventing negative values from rounding
#endif  
  // update number of photons (DEBUG: the following line prevents negative values caused by rounding errors)
  site->nPhotons[F246] = (site->nPhotons[F246] - nPhot < 0) ? 0 : site->nPhotons[F246] - nPhot;
#ifdef SX_DISPLAY_STATS
  sxRunData->nAbsPhotons[F246] += nPhot;
#endif

}

/*
 * \brief Solve H and H2 ionization at frequencies above 15.2 eV
 */
void sx_chem_solve_F152( struct sxSite *site )
{
  // calculate number of available species                                                                                          
  double nSpecH = site->nNucleons * site->xH;
  double nSpecH2 = site->nNucleons * site->xH2;
#ifdef SX_RECOMBINE
  // correcting number of absorbed photons for recombination
  if (site->H_rec_fac > recFactThr)
    nSpecH += site->nNucleons * site->H_rec_fac;
#endif
  double nSpec = nSpecH + nSpecH2;
  if (nSpec==0) return;  // return if there are no available species for a reaction

  // some initial constants
  double sigmaH = All.sxSigma[S152IH] * site->xH;
  double sigmaH2 = All.sxSigma[S152IH2] * site->xH2;
  double sigma = sigmaH + sigmaH2;
  if (sigma==0) return;  // skip this frequency if there is not a cross-section for the reaction

  // caclulate number of attenuated photons
  double nPhot = site->nPhotons[F152] * (1.0-exp( -1.0 * site->numdr * sigma ));
  if (nPhot==0) return;
  double nPhotH = nPhot * sigmaH / sigma;
  double nPhotH2 = nPhot * sigmaH2 / sigma;
  if (nPhotH>nSpecH) nPhotH=nSpecH;
  if (nPhotH2>nSpecH2) nPhotH2=nSpecH2;
  nPhot = nPhotH + nPhotH2;

#ifdef SX_RADIATION_PRESSURE
  site->momentum += nPhotH  * ( All.sxFreqBins[F152] * PLANCK + All.sxEnergy[E152IH]  ) * sxRunData->cInv; // erg/c
  site->momentum += nPhotH2 * ( All.sxFreqBins[F152] * PLANCK + All.sxEnergy[E152IH2] ) * sxRunData->cInv; // erg/c
#endif

  double attRatio;

  // update rates and mass fractions of H
  attRatio = nPhotH / (site->nNucleons * site->initxH); 
  double maxAttRatio = site->nPhotons[F152] * All.sxSigma[S152IH] * site->meanNgbDist / site->volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RIH] += attRatio * sxRunData->runTimeInv;
  site->rates[HRIH] += attRatio * sxRunData->runTimeInv * All.sxEnergy[E152IH];
#ifdef SX_RECOMBINE
  double red_factor = (1.0 - nPhotH / nSpecH);
  if (!(red_factor > 1.0e-6)) red_factor = 0.0;
  site->H_rec_fac *= red_factor; 
  site->xH *= red_factor;
#else
  site->xH -= attRatio * site->initxH;
  if (site->xH < 0.0) site->xH = minMassFrac; // preventing negative values from rounding
#endif

  // update rates and mass fractions of H2
  attRatio = nPhotH2 / (site->nNucleons * site->initxH2);
  maxAttRatio = site->nPhotons[F152] * All.sxSigma[S152IH2] * site->meanNgbDist / site->volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RIH2] += attRatio * sxRunData->runTimeInv;
  site->rates[HRIH2] += attRatio * sxRunData->runTimeInv * All.sxEnergy[E152IH2];
  site->xH2 -= attRatio * site->initxH2;
  if (site->xH2 < 0.0) site->xH2 = minMassFrac; // preventing negative values from rounding

  // update number of photons (DEBUG: the following line prevents negative values caused by rounding errors)
  site->nPhotons[F152] = (site->nPhotons[F152] - nPhot < 0) ? 0 : site->nPhotons[F152] - nPhot;
#ifdef SX_DISPLAY_STATS
  sxRunData->nAbsPhotons[F152] += nPhot;
#endif

}

/*
 * \brief Solve frequencies between 13.6 eV - 15.2 eV
 */
void sx_chem_solve_F136( struct sxSite *site )
 {
  // calculate number of available species                                                                                          
  double nSpecH =  site->nNucleons * site->xH;
  double nSpecH2 = site->nNucleons * site->xH2;
#ifdef SX_RECOMBINE
  // correcting number of absorbed photons for recombination
  if (site->H_rec_fac > recFactThr)
    nSpecH += site->nNucleons * site->H_rec_fac;
#endif
  double nSpec = nSpecH + nSpecH2;
  if (nSpec==0) return;  // return if there are no available species for a reaction

  // some initial constants                                                                        
  double sigmaH =  All.sxSigma[S136IH] * site->xH;
  double sigmaH2 = All.sxSigma[S136DH2] * site->xH2;
  double sigma = sigmaH + sigmaH2;
  if (sigma==0) return;  // skip this frequency if there is not a cross-section for the reaction

  // caclulate number of attenuated photons                                                  
  double nPhot= site->nPhotons[F136] * (1.0-exp( -1.0 * site->numdr * sigma ));  
  if (nPhot==0) return;
  double nPhotH =  nPhot * sigmaH / sigma;
  double nPhotH2 = nPhot * sigmaH2 / sigma;
  if (nPhotH>nSpecH) nPhotH=nSpecH;
  if (nPhotH2>nSpecH2) nPhotH2=nSpecH2;
  nPhot = nPhotH + nPhotH2;

#ifdef SX_RADIATION_PRESSURE
  site->momentum += nPhot * ( All.sxFreqBins[F136] * PLANCK + All.sxEnergy[E136IH] ) * sxRunData->cInv; // erg/c
#endif

  /*
#ifdef SX_RECOMBINE
  if (site->isSource)
    printf("SX: (%d:%d:%d) Rot %d Hrf %.03e xHp %.03e nspecH %.03e %.03e nphot %.03e \n",ThisTask,site->index,P[site->index].ID,
	   sxRunData->currentRot+1, site->H_rec_fac, site->initxHp, site->nNucleons * site->xH, nSpecH, nPhotH);
#endif
  */

  double attRatio;

  // update rates and mass fractions of H                                                              
  attRatio = nPhotH / (site->nNucleons * site->initxH);
  double maxAttRatio = site->nPhotons[F136] * All.sxSigma[S136IH] * site->meanNgbDist / site->volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RIH] += attRatio * sxRunData->runTimeInv;  // [1/s]
  site->rates[HRIH] += attRatio * sxRunData->runTimeInv * All.sxEnergy[E136IH]; // [erg/s]
#ifdef SX_RECOMBINE
  double red_factor = (1.0 - nPhotH / nSpecH);
  if (!(red_factor > 1.0e-6)) red_factor = 0.0;
  site->H_rec_fac *= red_factor; 
  site->xH *= red_factor;
#else
  site->xH -= attRatio * site->initxH;
  if (site->xH < 0.0) site->xH = minMassFrac; // preventing negative values from rounding
#endif

  // update rates and mass fractions of LW
  attRatio = nPhotH2 / (site->nNucleons * site->initxH2);
  maxAttRatio = site->nPhotons[F136] * All.sxSigma[S136DH2] * site->meanNgbDist / site->volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RDH2] += attRatio * sxRunData->runTimeInv;
  site->xH2 -= attRatio * site->initxH2;
  if (site->xH2 < 0.0) site->xH2 = minMassFrac; // preventing negative values from rounding

  // update number of photons  (DEBUG: the following line prevents negative values caused by rounding errors)
  site->nPhotons[F136] = (site->nPhotons[F136] - nPhot < 0) ? 0 : site->nPhotons[F136] - nPhot;
#ifdef SX_DISPLAY_STATS
  sxRunData->nAbsPhotons[F136] += nPhot;
#endif

}

/*
 * \brief Solve LW radiation between 11.2 eV - 13.6 eV
 */
void sx_chem_solve_F112( struct sxSite *site )
{
  // calculate number of available species                                                                                          
  double nSpecH2 = site->nNucleons * site->xH2;
  double nSpecH  = site->nNucleons * site->xH;
  double nSpec = nSpecH + nSpecH2;

  if (nSpec==0) return;  // return if there are no available species for a reaction

  // some initial constants                                                                                 
  double sigmaH2 = All.sxSigma[S112DH2] * site->xH2 * (1.0 + uvPumpFact);
  double sigmaH = All.sxSigma[S112H] * site->xH;
  double sigma = sigmaH2 + sigmaH;
  if (sigma==0) return;  // skip this frequency if there is not a cross-section for the reaction

  // Calculate number of photons absorbed in cell
  // NB. Treatment of absorption by H2 is a crude approximation that doesn't properly account for H2 self-shielding and thus
  // will lead to PDRs that are too narrow in the cases where they are resolved
  double nPhot  = site->nPhotons[F112] * (1.0-exp( -1.0 * site->numdr * sigma ));
  if (nPhot==0) return;
  double nPhotH  = nPhot * sigmaH  / sigma;
  double nPhotH2 = nPhot * sigmaH2 / sigma;
  double nPhotdis = nPhotH2 / (1.0 + uvPumpFact);

  if (nPhotdis > nSpecH2)
    {
      nPhotdis = nSpecH2;
      nPhotH2  = (1.0 + uvPumpFact) * nPhotdis;
    }
  // No corresponding check for H because absorption of LW photons by H atoms doesn't destroy the atoms
  nPhot = nPhotH + nPhotH2;


#ifdef SX_RADIATION_PRESSURE
  // DEBUG: this should be discussed first
  //site->momentum += nPhotH2 * ( All.sxFreqBins[F112] * PLANCK + All.sxEnergy[E112DH2] ) * sxRunData->cInv; // erg/c
#endif

  double attRatio;
  // update rate and H2 fraction
  attRatio = nPhotdis / (site->nNucleons * site->initxH2);
  double maxAttRatio = site->nPhotons[F112] * All.sxSigma[S112DH2] * (1.0 + uvPumpFact) * site->meanNgbDist / site->volume;
  if (attRatio > maxAttRatio) attRatio = maxAttRatio;
  site->rates[RDH2] += attRatio * sxRunData->runTimeInv;
  site->xH2 -= attRatio * site->initxH2;

  // update number of photons  (DEBUG: the following line prevents negative values caused by rounding errors)
  site->nPhotons[F112] = (site->nPhotons[F112] - nPhot < 0) ? 0 : site->nPhotons[F112] - nPhot;
#ifdef SX_DISPLAY_STATS
  sxRunData->nAbsPhotons[F112] += nPhot;
#endif

}

/*!
 * \brief solve the rate equation for a site
 */
void sx_chem_solve_site(struct sxSite *site)
{
  /*
  if (site->isSource)
    {
      printf("SX: (%d:%d:%d) source phot [ ",ThisTask,P[site->index].ID,site->index);
      for(int i=0; i<SX_NFREQ; i++)
	{
	  printf("%.02e ",site->nPhotons[i]);
	}
      printf("] \n");
    }
  */

#if defined(SX_SKIP_RADIUS)
  if (sxCell[site->index].skip==0)
    {
#endif

  // solve frequencies 24.6+
  if( site->nPhotons[F246] > 0.0 )
    sx_chem_solve_F246( site );
  
  // solve frequencies 15.2 - 24.6 eV
  if( site->nPhotons[F152] > 0.0 )
    sx_chem_solve_F152( site );
  
  // solve frequencies btw. 13.6 - 15.2 eV
  if( site->nPhotons[F136] > 0.0 )
    sx_chem_solve_F136( site );
  
  // solve frequencies btw. 11.2 - 13.6 eV
  if( site->nPhotons[F112] > 0.0 )
    sx_chem_solve_F112( site );

#if defined(SX_SKIP_RADIUS)
    }
#endif

  for (int f=0; f<SX_NFREQ; f++)
    {
      // dump small number of photons
      if ( site->nPhotons[f] < site->nNucleons * All.sxMinNumPhotons )
	{
#ifdef SX_DISPLAY_STATS
	  sxRunData->nLostPhotons[f] += site->nPhotons[f];
#endif
	  site->nPhotons[f] = 0.0;
	}
      // precalculate photon attenuation ratio for each frequency
      site->nPhotonsRatio[f] = (site->nPhotonsOriginal[f]>0.0) ? site->nPhotons[f]/site->nPhotonsOriginal[f] : 0.0;
    }

}

#endif
