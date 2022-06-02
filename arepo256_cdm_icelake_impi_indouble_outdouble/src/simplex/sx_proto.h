/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Definitions of variables and structures in the SPRAI
 * \details     
 */

#ifndef SX_PROTO_H
#define SX_PROTO_H

#include <time.h>
#include <gsl/gsl_rng.h>
#include "../allvars.h"
#include "../dtypes.h"
#include "../SGChem/sgchem_def.h"
#include "sx_def.h"

//======================================================================================================//
//                                     Run variables structure                                          //
//======================================================================================================//

// Information about the RT run and various statistics
extern struct sxRunData_struct
{
  int currentStep;                     // current number of steps
  int currentRot;                      // current rotation step
  int nSteps;                          // total number of steps in the current run

  int nSourcesExport;                  // number of source particles on the local task
  int nSourcesLocal;                   // number of source sites on the local task
  int nSources;                        // number of all sources

  clock_t RunClock;                    // starting time of the run
#if SX_NUM_ROT > 1
  clock_t RotClock;                    // starting time of the rotation
#endif

  // set of directional bins
  // orthoBase is an array with a form [dir][vector x,y,z][coordinate] where z-vector is the directional vector
  double orthoBase[SX_NDIR][SX_NDIM][SX_NDIM]; // set of orthogonal vectors
  double solidAngles[SX_NDIR];                 // solid angles corresponding to each direction bin

#ifdef SX_DISPLAY_STATS
  double nInitPhotons[SX_NFREQ];       // initial number of photons in the system
  double nEmitPhotons[SX_NFREQ];       // number of photons emited by the sources
  double nAbsPhotons[SX_NFREQ];        // number of photons absorbed by gas
  double nLostPhotons[SX_NFREQ];       // number of lost photons due to numerical cutouts
  double nFinPhotons[SX_NFREQ];        // number of remaining photons after the run

  int nDifTransp;                      // number of diffusion transports during the RT step per task
  int nBalTransp;                      // number of ballistic transports during the RT step per task
  int nRecPackets;                     // number of photon packets received by tasks
  int nSendPackets;                    // number of photon packets send to other tasks
#endif

  // time settings
  double time;                         // current physical time of the simulation [s]
#if SX_RUNNING_MODE == 2 || SX_RUNNING_MODE == 1
  double runTime;                      // time of the radiation RT run [s]
  double runTimeInv;                   // 1/runTime
  double timeStep;                     // runTime/SX_NUM_STEPS = duration of the one RT step [s]
#endif

  // conversion factors from comoving-code-units to physical-cgs-units
  double mass2cgs;                     // mass [g]
  double density2cgs;                  // density [g/cm^3]
  double volume2cgs;                   // volume [cm^3]
  double length2cgs;                   // length [cm]
  double velocity2cgs;                 // velocity [cm/s]

  // Cosmological parameters
  double expansionFactor;              // current expansion factor
  double redshift;                     // current redshift
#ifdef SX_CMB_TEMP_LIMIT
  double tempCMB;                      // temperature of the CMB at the current redshift
#endif

  // Other settings
#ifdef SX_RADIATION_PRESSURE
  double cInv;                         // 1/c inverted speed of light
#endif

} *sxRunData;

//======================================================================================================//
//                                     Site structure                                                   //
//======================================================================================================//

// Used for solving of a rate equation on sites
struct sxSite
{
  int index;                           // index of a SPH particle (site)
  _Bool isSource;                      // 1 = cell with a star, 0 = normal cell
  double dt;                           // time step
  
  double nPhotons[SX_NFREQ];           // total number of photons
  double nPhotonsOriginal[SX_NFREQ];   // original total number of photons (before solving)
  double nPhotonsSource[SX_NFREQ];     // number of photons emitted by a source in the site
  double nPhotonsRatio[SX_NFREQ];      // ratio nPhotons/nPhotonsOriginal

#if SX_CHEMISTRY == 3
  double fluxes[SX_NFREQ];             // photon fluxes in the cell (defined in "sgchem_def.h")
  double rates[SX_NRATES];             // rates in the cell (defined in "sgchem_def.h")

#ifdef SX_RADIATION_PRESSURE
  double momentum;                     // momentum from the radiation pressure [erg/c]
#endif

  double xH, initxH;                   // nucleon fraction of H atoms
  double xH2, initxH2;                 // nucleon fraction of H2 molecules
  double xHp, initxHp;                 // nucleon fraction of H+ atoms
  double xHe, initxHe;                 // nucleon fraction of He atoms
  double xHep, initxHep;               // nucleon fraction of He+ atoms

  double nNucleons;                    // number of nucleons H, H2 and H+

  double density;                      // density in the cell in physical [g/cm^3]
  double volume;                       // volume of the cell in physical [cm^3]
  double meanNgbDist;                  // mean neighbor distance [cm]
  double numdens;                      // nucleon number density
  double numdr;                        // column nucleon density

#elif SX_CHEMISTRY == 4
  double nFract[SX_NMASSFRACT];        // number fraction
  double nNucleons;                    // number of nucleons H
  double density;                      // density in the cell in physical [g/cm^3]
  double volume;                       // volume of the cell in physical [cm^3]
  double numdens;                      // nucleon number density
  double numdr;                        // column nucleon density
  double temp;                         // temperature of the cell
  double heating;                      // cumulated heating rate
  double cooling;                      // cumulated cooling rate
  double scFactor;                     // sub-cycling factor in sites with too many photons
#endif
#ifdef SX_RECOMBINE
  double H_rec_fac;                    // to account for additional absorption from
  double He_rec_fac;                   // recombined H, He atoms
#endif
};

//======================================================================================================//
//                                     Source structure variables                                       //
//======================================================================================================//

// List of ionisation rates from different sources
extern struct sxSource_struct
{
  int index;                           // index of the gas particle in the SphP
  double ionRate[SX_NFREQ];            // ionisation rate [ph/s/UnitPhotons_per_s]
} *sxSources;

//======================================================================================================//
//                                     Additional SphP cell variables                                   //
//======================================================================================================//

// Array of SimpleX cells (parallel to SphP) that stores information needed only during the RT run
extern struct sxCell_struct
{
#ifdef SX_SKIP_RADIUS
  int skip;                            // skip calculations of radiation in this cell
#endif
#ifdef SX_RECOMBINE
  double H_rec_fac;                    // additional number of species from recombination
  double He_rec_fac;                   // additional number of species from recombination
#endif
  double MeanNeighDistance;                       // mean distance to the Delaunay neighbour [code units]
#if SX_CHEMISTRY == 3
  double TracAbund[SGCHEM_NUM_ADVECTED_SPECIES];  // image of the TracAbund from the SGChem
#endif
  int Source;                                     // index of a source in sxSources (no source = -1) 
#ifdef SX_RADIATION_PRESSURE
  double MomentumChange[SX_NDIM];      // cummulative momentum vector change [code units]
#endif
} *sxCell;

//======================================================================================================//
//                                     Photon Packet Transport variables                                //
//======================================================================================================//

extern int sxNumSPD;                   // Number of SphP Directions ( NumGas * SX_NDIRS )

extern struct sxABPP_struct            // Photon Packages structure
{
  MySxDouble nPhotons[SX_NFREQ];       // Photons for exchange
} *sxAPP, *sxBPP,                      // fake Photon Packages
  *sxRealAPP, *sxRealBPP;              // real Photon Packages
extern int sxMaxNumABPP;               // maximum number of PP in the APP/BPP arrays
extern int sxSwapABPP;                 // indication if the ABPP arrays are swapped
extern int sxNumReallocABPP;           // counter that counts number of memory reallocations
 
extern int sxNumAPP, sxNumBPP;                // current number of PP in the APP/BPP arrays
extern int *sxSPDtoAPP, *sxSPDtoBPP;          // fake connections between SPD and APP/BPP
extern int *sxRealSPDtoAPP, *sxRealSPDtoBPP;  // real connections between SPD and APP/BPP

//======================================================================================================//
//                                     Photon Packet Exchange variables                                 //
//======================================================================================================//

extern int sxNumTotEDP;                // number of all External DP
extern int sxNumTotEDPD;               // number of all External DP Directions ( sxNumTotEDP * SX_NDIRS )

extern struct sxQPP_struct             // Photons Packets added to the exchange Queue
{
  MySxDouble nPhotons[SX_NFREQ];       // photons for exchange
  int index;                           // SphP index of the particle on a corresponding task
} *sxQPP; 
extern int sxMaxNumQPP;                // number of photon packets
extern int *sxEDPDtoQPP;               // indexQPP[task][indexEDP][dir] ( sxNumTotEDPD )
extern int sxNumReallocQPP;            // counter that counts number of memory reallocations

extern struct sxEPP_struct             // Photon Packet Exchange structures
{
  MySxDouble nPhotons[SX_NFREQ];       // photons for exchange
  int index;                           // SphP index of the particle on a corresponding task
  int dir;                             // directional bin
} *sxExportPP;
extern int sxNumExportPP, sxNumImportPP;      // number of photon packets
extern int *sxNumEDP, *sxOffEDP;              // numEDP[task] and offsetsEDP[task] (NTask)
extern int *sxDPtoEDP;                        // indexEDP[indexDP] (Mesh.Ndp)

//======================================================================================================//
//                                     Other stuff variables                                            //
//======================================================================================================//

extern int SXPID;                      //DEBUG: ID of the debugged particle
extern int SXID;                       //DEBUG: current index and process of the sph particle in SphP and P arrays
                                       //       this index can be changed when the particle mash is recalculated!!!

extern gsl_rng *sxRand;   // random number used during the calculation of a directional vector

//======================================================================================================//
//                                     Commont function definitions                                     //
//======================================================================================================//

// Debugging functions
void sx_debug_start_run(void);
void sx_debug_end_run(void);
void sx_debug_site( struct sxSite *site, int );
void sx_debug_particle( int, int );
void sx_debug_sources(void);
void sx_debug_display_stats(void);
void sx_debug_display_memory(void);
void   sx_debug_timer_start(int);
double sx_debug_timer_end(int);
double sx_debug_timer_switch( int, int );

// Mathematical and numerical functions
void sx_math_compute_rotational_matrix( double M[3][3], double d );
double sx_math_compute_inner_product( double[SX_NDIM], double[SX_NDIM] );
double sx_math_integrate_discrete_function( double, double, double[], double[], int );

// Main functions
void sx_restart(void);
void sx_initialize(void);
void sx_evolve( double a, double t, double dt );
void sx_evolve_hd(void);
void sx_evolve_pp(void);

// General Photon Transport functions
void sx_phot_initialize(void);
void sx_phot_start_run(void);
void sx_phot_start_step(void);
void sx_phot_end_step(void);
#if SX_NUM_ROT > 1
void sx_phot_start_rot(void);
void sx_phot_end_rot(void);
#endif
void sx_phot_end_run(void);
_Bool sx_phot_initialize_site( struct sxSite *site , int );
_Bool sx_phot_check_n_photons( struct sxSite *site );

// Photon Sources functions
void sx_phot_sources_start_run(void);
void sx_phot_sources_end_run(void);

// Photon Exchange functions
void sx_phot_exchange_start_run(void);
void sx_phot_copy_PP_to_QPP( int dp, int d, MySxDouble nPhotons[SX_NFREQ] );
void sx_phot_copy_PP_to_APP( int i, int d, MySxDouble nPhotons[SX_NFREQ] );
void sx_phot_exchange_photons(void);
void sx_phot_exchange_end_run(void);

// Photon transport
void sx_phot_transport_start_run(void);
void sx_phot_diffuse_transport( struct sxSite *site );
void sx_phot_ballistic_transport( struct sxSite *site );
void sx_phot_transport_copy(void);
void sx_phot_transport_end_run(void);

// Chemistry specific functions
#if SX_CHEMISTRY == 3 // SGChem chemistry
void sx_chem_compute_cross_sections( double );
void sx_chem_compute_mean_photon_energies( double );
#endif

// Common chemistry functions
double sx_chem_image1(int index);
double sx_chem_image2(int index);
double sx_chem_image_all(int index, int prop);

void sx_chem_initialize(void);
void sx_chem_compute_ionisation_rates( double[SX_NFREQ], int );
#ifdef POPIII_SNE
void sx_chem_compute_popIII_ionisation_rates(double ionisationRate[SX_NFREQ], int index);
#endif
void sx_chem_start_run(void);
void sx_chem_initialize_site( struct sxSite *site );
void sx_chem_solve_site( struct sxSite *site );
void sx_chem_finalize_site( struct sxSite *site );
void sx_chem_end_run(void);

// PopIII Lookup tables
void setup_popiii_tables_(void);
void popiii_lookup_(double *accretion, double *mass, double *radius, double *temperature);

#endif
