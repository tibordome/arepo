/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne.c
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel
 * \brief       Enables Supernova feedback. This file includes the main functions
 *              which perform the energy injection and control the workflow of the
 *              module.
 * \details     Please contact the authors at robin.tress@uni-heidelberg.de
 *              and rowan.smith@manchester.ac.uk before using this to avoid overlapping
 *              of projects. And please report any issue you may encounter by using this
 *              routine.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

/*! \brief Main workflow of the SNe feedback routine
 *
 *  -the function is_it_time_for_a_new_sn() will take care about when it is time for a SN injection
 *  -the function determine_position_and_injection_region_of_sn() decides where to position the SN and
 *   which cells are affected by the energy injection (the indices of the interested cells on the local
 *   processor are stored in the variable indices of size local_n)
 *  -we then find the gas propreties of that region with find_gas_properties_of_injection_region() and
 *   smooth the region with distribute_mass()
 *  -Finally we inject the energy with do_sn_energy_injection(). The injection scheme used is returned
 *   (injection_scheme_used = 0 if thermal energy is injected, injection_scheme_used = 1 if the SN is not resolved and
 *   momentum is injected instead).
 *  -At last we inject tracer particles (if the flag is active) and log the specifics of the just created SN into the SN
 *   log file
 *
 *  TODO: -allow for mass injection with each SN (mass return from the stellar population)
 *        -allow for variable energy injection so this can be recycled for wind feedback
 */
void sne_feedback(void)
{
  double sne_pos[3], radius;
  int local_n, total_n, *indices; /*local_n is the number of injection cells on the local processor and indices stores their indices*/

  enum SNE_type timeForANewSN;

  mpi_printf("SNe Feedback %f\n", All.Time);

  double t0 = second();

#if defined(TRACER_MC) && defined(SINK_PARTICLES_FEEDBACK_RETURN_MASS)
  start_MC_tracer(N_tracer);
#endif

  while(1)
    {
      /*When*/
      timeForANewSN = is_it_time_for_a_new_sn();

      if(!timeForANewSN)
        {
          mpi_printf("SNe: no further supernovae for now\n");
          MPI_Barrier(MPI_COMM_WORLD);
          break;
        }
      mpi_printf("SNe: it is time for a new supernova!\n");
#ifdef POPIII_SNE
      /* Get mass, explosion energies */
      double sn_mass, sn_E51;
      sne_sink_feedback_PROPERTIES(&sn_mass, &sn_E51);
#else
      double sn_E51        = 1.0;
      double sn_mass       = 0.0;
#endif

      int sn_size;
      if(sn_E51 > 1.0)
        sn_size = sn_E51;
      else
        sn_size = 1;
      /*Where*/
      determine_position_and_injection_region_of_sn(sne_pos, timeForANewSN, &radius, &local_n, &total_n, &indices, sn_size);

      mpi_printf("SNe: Position %g, %g, %g\n", sne_pos[0], sne_pos[1], sne_pos[2]);

#ifdef POPIII_SNE
      if(sn_E51 == 0.0) /* fallback! */
        {
          mpi_printf("SNe: Fallback, no SN, mass= %g M_sun\n", sn_mass);
          MPI_Barrier(MPI_COMM_WORLD);
          myfree(indices);
          continue;
        }
#endif

#ifdef SNE_solarring
      double Rinner_cut, xgal, ygal, rgal;

      Rinner_cut = 50.;
      xgal       = sne_pos[0] - boxHalf_X;
      ygal       = sne_pos[1] - boxHalf_Y;
      rgal       = sqrt(xgal * xgal + ygal * ygal);

      if(rgal < Rinner_cut)
        {
          myfree(indices);
          continue;
        }
#endif

      /*Find gas properties of the injection region*/
      double total_mass, total_volume, mean_density, mean_utherm, mean_temp, vel_CM[3];
      find_gas_properties_of_injection_region(local_n, indices, &total_mass, &total_volume, &mean_density, &mean_utherm, &mean_temp,
                                              vel_CM);

      if(mean_density < 1e-29 / All.UnitDensity_in_cgs)
        {
          if(ThisTask == 0)
            warn("We lost a SN at position (x, y, z) = (%g, %g, %g)\n", sne_pos[0], sne_pos[1], sne_pos[2]);
          myfree(indices);
          continue;
        }

      /*distribute mass*/
      // distribute_mass(local_n, indices, mean_density, mean_utherm, vel_CM);

      /*inject mass from the SN event*/
      inject_mass(timeForANewSN, local_n, indices, total_volume, &total_mass, &mean_density, sn_mass);

      /*inject the energy*/
      int injection_scheme_used = do_sn_energy_injection(sne_pos, local_n, indices, radius, mean_density, total_mass, sn_E51);

#ifdef INJECT_TRACER_INTO_SN
      /*inject tracer particles*/
      int tracers_injected = inject_tracer_particles(local_n, indices, total_n, timeForANewSN);
#else
      int tracers_injected = 0;
#endif

      /*update_timesteps_around_injection_region(sne_pos);*/

      /*Log the SN just injected*/
      sne_log(sne_pos, total_n, injection_scheme_used, radius, mean_density, mean_temp, tracers_injected, sn_mass, sn_E51);

      myfree(indices);
#if defined(TREE_BASED_TIMESTEPS) && (defined(POPIII_SNE))
      /*update the timesteps of injection region and surroundings*/
      double centerBox[3] = {boxHalf_X, boxHalf_Y, boxHalf_Z};
      update_timesteps_around_injection_region(centerBox, boxHalf_X);
#endif
    }

#if defined(TRACER_MC) && defined(SINK_PARTICLES_FEEDBACK_RETURN_MASS)
  finish_MC_tracer();
#endif

  All.LastSNTime = All.Time;

  /*double t1 = second();
    mpi_printf("SNe: feedback done, took %g sec.\n", timediff(t0, t1));*/

  double tt0 = second();
#if defined(TREE_BASED_TIMESTEPS) && (!defined(POPIII_SNE))
  /*update the timesteps of injection region and surroundings*/
  double centerBox[3] = {boxHalf_X, boxHalf_Y, boxHalf_Z};
  update_timesteps_around_injection_region(centerBox, boxHalf_X);
#endif
  double tt1 = second();
  mpi_printf("SNe: check done 5, took %g sec.\n", timediff(tt0, tt1));

  double t1 = second();
  mpi_printf("SNe: feedback done, took %g sec.\n", timediff(t0, t1));

  return;
}

/*! \brief wrapper function to do the SN energy injection
 *
 *  everything is in internal units
 *
 *  \param sne_pos[3] position of SN
 *  \param local_n numer of injection particles on the local processor
 *  \param indices[] indices of the local injection cells
 *  \param radius radius of the injection region
 *  \param mean_density mean density of the injection region
 *  \param total_mass total mass of the cells in the injection region
 *
 *  returns the injection scheme used
 *
 */
int do_sn_energy_injection(double sne_pos[3], int local_n, int indices[], double radius, double mean_density, double total_mass,
                           double sn_E51)
{
  /*Decide wether we can do thermal energy injection or we have to do momentum injection*/

  int injection_scheme = decide_injection_scheme(mean_density, radius); /*1 = energy, 0 = momentum*/

  /*Do the injection*/
  if(injection_scheme)
    {
      mpi_printf("SNe: doing energy injection\n");
      energy_injection(local_n, indices, total_mass, sn_E51);
    }
  else
    {
      mpi_printf("SNe: doing momentum injection\n");
      momentum_injection(sne_pos, local_n, indices, mean_density, total_mass);
    }

  update_primitive_variables();
  mpi_printf("SNe: a star just exploded!!!\n");

  return injection_scheme;
}

/*! \brief decides if the SN would be resolved and we can use thermal energy injection or
 *         if we have to do momentum injection instead
 *
 *  \param mean_density mean density of the injection region
 *  \param radius radius of the injection region
 *
 *  returns 1 if we can use thermal energy injection
 *  returns 0 if we have to use momentum injection
 *
 */
int decide_injection_scheme(double mean_density, double radius)
{
  /* Given the mean density of the injection region we calculate the radius at the end
     of the Sedov-Taylor phase using formula from Draine 2011 p 4333 for a supernova with 1 E51 erg */

  double yn = (mean_density * All.UnitDensity_in_cgs) / ((1.0 + 4.0 * ABHE) * PROTONMASS);

  double R_sedov = 5.894e19 * pow(yn, -7. / 17) / All.UnitLength_in_cm; /*This is to match Gatto et al. 2014*/

#ifdef POPIII_SNE
  /* always do energy injection because momentum injection is not implemented in com. units */
  yn      = yn * All.cf_a3inv * All.HubbleParam * All.HubbleParam;
  R_sedov = 5.894e19 * pow(yn, -7. / 17) / All.UnitLength_in_cm;
  if(R_sedov < radius * All.cf_atime / All.HubbleParam)
    {
      warn("SNE WARNING: unresolved energy injection, R=%g R_sedov=%g n=%g \n", R_sedov, radius * All.cf_atime / All.HubbleParam, yn);
    }
  return 1;
#endif

  /* if then the injection radius is greater we are unresolved and have to do the momentum injection */
  if(radius <= R_sedov)
    return 1;
  else
    return 0;
}

/*! \brief injects 10^51 erg of THERMAL energy distributed over the injection cells
 *         with a fraction based on the mass fraction of each cell
 *
 *  \param local_n numer of injection particles on the local processor
 *  \param indices[] indices of the local injection cells
 *  \param total_mass total mass of the cells in the injection region
 *
 */
void energy_injection(int local_n, int indices[], double total_mass, double sn_E51)
{
  double total_energy = sn_E51 * 1E51 * All.HubbleParam / All.UnitEnergy_in_cgs;
  double gamma;

  for(int i = 0; i < local_n; i++)
    {
      int idx = indices[i];

      /* Set energies */
      double energy_fraction = P[idx].Mass / total_mass;
      SphP[idx].Utherm += total_energy * energy_fraction / P[idx].Mass;
      SphP[idx].Energy += All.cf_atime * All.cf_atime * total_energy * energy_fraction;
#ifdef MHD_THERMAL_ENERGY_SWITCH
      SphP[idx].Etherm = SphP[idx].Utherm * P[idx].Mass;
#endif

      /* Update other variables*/
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      // Update Entropy, Taken from SGChem
      gamma             = 5.0 / 3.0;  // Gas in injection region will be fully ionized, so gamma = 5/3
      SphP[idx].A       = (gamma - 1.0) * SphP[idx].Utherm / pow(SphP[idx].Density * All.cf_a3inv, GAMMA_MINUS1);
      SphP[idx].Entropy = log(SphP[idx].A) * P[idx].Mass;
#endif

#ifdef SGCHEM
      // chemistry (ionize everything)
#if CHEMISTRYNETWORK == 1
      SphP[idx].TracAbund[IH2]     = 0.;
      SphP[idx].TracAbund[IHP]     = 1.;
      SphP[idx].TracAbund[IDP]     = All.DeutAbund;
      SphP[idx].TracAbund[IHD]     = 0.;
      SphP[idx].TracAbund[IHEP]    = 0.;
      SphP[idx].TracAbund[IHEPP]   = ABHE;
      SphP[idx].TracAbund[IHATOM]  = 0.;
      SphP[idx].TracAbund[IHEATOM] = 0.;
      SphP[idx].TracAbund[IDATOM]  = 0.;
#endif
#if CHEMISTRYNETWORK == 4
      SphP[idx].TracAbund[IH2]    = 0.;
      SphP[idx].TracAbund[IHP]    = 1.;
      SphP[idx].TracAbund[IHATOM] = 0.;
#endif
#if CHEMISTRYNETWORK == 5
      SphP[idx].TracAbund[IH2]    = 0.;
      SphP[idx].TracAbund[IHP]    = 1.;
      SphP[idx].TracAbund[IHATOM] = 0.;
#ifdef SGCHEM_VARIABLE_Z
      SphP[idx].TracAbund[ICP] = SphP[idx].CarbAbund;
      SphP[idx].TracAbund[ICO] = 0.;
#else
      SphP[idx].TracAbund[ICP]    = All.CarbAbund;
      SphP[idx].TracAbund[ICO]    = 0.;
#endif
#endif /* CHEMISTRYNETWORK == 5 */
#if CHEMISTRYNETWORK == 15
      SphP[idx].TracAbund[IH2]     = 0.;
      SphP[idx].TracAbund[IHP]     = 1.;
      SphP[idx].TracAbund[ICHX]    = 0.;
      SphP[idx].TracAbund[IOHX]    = 0.;
      SphP[idx].TracAbund[ICO]     = 0.;
      SphP[idx].TracAbund[IHCOP]   = 0.;
      SphP[idx].TracAbund[IHATOM]  = 0.;
      SphP[idx].TracAbund[IHEP]    = ABHE;
      SphP[idx].TracAbund[IHEATOM] = 0.;
#ifdef SGCHEM_VARIABLE_Z
      SphP[idx].TracAbund[ICP]    = SphP[idx].CarbAbund;
      SphP[idx].TracAbund[ICATOM] = 0.;
      SphP[idx].TracAbund[IOATOM] = SphP[idx].OxyAbund;
      SphP[idx].TracAbund[IMP]    = SphP[idx].MAbund;
      SphP[idx].TracAbund[IMATOM] = 0.;
#else
      SphP[idx].TracAbund[ICP]    = All.CarbAbund;
      SphP[idx].TracAbund[ICATOM] = 0.;
      SphP[idx].TracAbund[IOATOM] = All.OxyAbund;
      SphP[idx].TracAbund[IMP]    = All.MAbund;
      SphP[idx].TracAbund[IMATOM] = 0.;
#endif
#endif /* CHEMISTRYNETWORK == 15 */

      for(int j = 0; j < SGCHEM_NUM_SPECIES; j++)
        SphP[idx].MassTracAbund[j] = P[idx].Mass * SphP[idx].TracAbund[j];
#endif /*SGCHEM*/
    }
}

/*! \brief Injects momentum energy to each affected cell.
 *
 *         The amount of momentum injected is calculated based on how much momentum a SN
 *         remnant should have at the end of the Sedov-Taylor phase given the ambient density.
 *         The interior of the bubble is then heated to at least 10^4 K.
 *
 *  \param sne_pos[] position of the SN
 *  \param local_n numer of injection particles on the local processor
 *  \param indices[] indices of the local injection cells
 *  \param mean_density mean density of the injection region
 *  \param total_mass total mass of the cells in the injection region
 *
 *
 */
void momentum_injection(double sne_pos[3], int local_n, int indices[], double mean_density, double total_mass)
{
  /*Use the same total momentum as in Gatto et al. 2014 ultimately from Blondin et al. 1998*/
  /* XXX: currently not compatible with comoving coordinates */
  /* No longer adding additional thermal energy as we could overshoot and are at the
     terminal momentum anyway*/

  double temp_set  = 1.e4; /*temperature of the interior of the bubble*/
  double energy_SN = 1e51 / All.UnitEnergy_in_cgs;

  double yn = (mean_density * All.UnitDensity_in_cgs) / ((1.0 + 4.0 * ABHE) * PROTONMASS);

  double total_momentum = 2.6e+05 * pow(yn, -2.0 / 17.0); /*Units of solar mass km/s*/
  total_momentum        = total_momentum * (SOLAR_MASS / All.UnitMass_in_g) * (1.e5 / All.UnitVelocity_in_cm_per_s);

  double total_energy = total_momentum * total_momentum / 2. / total_mass;

  mpi_printf("DEBUG_SNE: total_energy = %g erg\n", total_energy * All.UnitEnergy_in_cgs);
  total_momentum = (total_energy > energy_SN) ? sqrt(2. * total_mass * energy_SN) : total_momentum;

  double inj_momentum[3];

  for(int i = 0; i < local_n; i++)
    {
      int idx = indices[i];

      /* Set momentum (radially from sne location)*/
      double dx = P[idx].Pos[0] - sne_pos[0];
      double dy = P[idx].Pos[1] - sne_pos[1];
      double dz = P[idx].Pos[2] - sne_pos[2];
      double r  = sqrt(dx * dx + dy * dy + dz * dz);

      double energy_fraction = P[idx].Mass / total_mass;
      inj_momentum[0]        = total_momentum * energy_fraction * dx / r;
      inj_momentum[1]        = total_momentum * energy_fraction * dy / r;
      inj_momentum[2]        = total_momentum * energy_fraction * dz / r;

      /* Inject Momentum*/
      SphP[idx].Momentum[0] += inj_momentum[0];
      SphP[idx].Momentum[1] += inj_momentum[1];
      SphP[idx].Momentum[2] += inj_momentum[2];

      /* Update other variables*/
      // Save old kinetic energy
      double Ekin_old =
          0.5 * P[idx].Mass * (P[idx].Vel[0] * P[idx].Vel[0] + P[idx].Vel[1] * P[idx].Vel[1] + P[idx].Vel[2] * P[idx].Vel[2]);

      // Velocity (accroding to new momentum)
      P[idx].Vel[0] = SphP[idx].Momentum[0] / P[idx].Mass;
      P[idx].Vel[1] = SphP[idx].Momentum[1] / P[idx].Mass;
      P[idx].Vel[2] = SphP[idx].Momentum[2] / P[idx].Mass;

      // Compute new kinetic energy and change from old to new
      double Ekin_new =
          0.5 * P[idx].Mass * (P[idx].Vel[0] * P[idx].Vel[0] + P[idx].Vel[1] * P[idx].Vel[1] + P[idx].Vel[2] * P[idx].Vel[2]);
      double dEkin = Ekin_new - Ekin_old;

#ifdef SGCHEM
      // chemistry (ionize everything but He and metals)
#if CHEMISTRYNETWORK == 1
      SphP[i].TracAbund[IH2]    = 0.;
      SphP[i].TracAbund[IHP]    = 1.;
      SphP[i].TracAbund[IDP]    = All.DeutAbund;
      SphP[i].TracAbund[IHD]    = 0.;
      SphP[i].TracAbund[IHATOM] = 0.;
      SphP[i].TracAbund[IDATOM] = 0.;
#endif
#if CHEMISTRYNETWORK == 4
      SphP[i].TracAbund[IH2]    = 0.;
      SphP[i].TracAbund[IHP]    = 1.;
      SphP[i].TracAbund[IHATOM] = 0.;
#endif
#if CHEMISTRYNETWORK == 5
      SphP[i].TracAbund[IH2]    = 0.;
      SphP[i].TracAbund[IHP]    = 1.;
      SphP[i].TracAbund[IHATOM] = 0.;
#ifdef SGCHEM_VARIABLE_Z
      SphP[i].TracAbund[ICP] = SphP[i].CarbAbund;
      SphP[i].TracAbund[ICO] = 0.;
#else
      SphP[i].TracAbund[ICP]      = All.CarbAbund;
      SphP[i].TracAbund[ICO]      = 0.;
#endif
#endif /* CHEMISTRYNETWORK == 5 */
#if CHEMISTRYNETWORK == 15
      SphP[i].TracAbund[IH2]    = 0.;
      SphP[i].TracAbund[IHP]    = 1.;
      SphP[i].TracAbund[ICHX]   = 0.;
      SphP[i].TracAbund[IOHX]   = 0.;
      SphP[i].TracAbund[ICO]    = 0.;
      SphP[i].TracAbund[IHCOP]  = 0.;
      SphP[i].TracAbund[IHATOM] = 0.;
#ifdef SGCHEM_VARIABLE_Z
      SphP[i].TracAbund[ICP]    = SphP[i].CarbAbund;
      SphP[i].TracAbund[ICATOM] = 0.;
      SphP[i].TracAbund[IOATOM] = SphP[i].OxyAbund;
#else
      SphP[i].TracAbund[ICP]      = All.CarbAbund;
      SphP[i].TracAbund[ICATOM]   = 0.;
      SphP[i].TracAbund[IOATOM]   = All.OxyAbund;
#endif
#endif /* CHEMISTRYNETWORK == 15 */

      // temperature (bring everything to at least 10^4 K)
      double gas_energy_density = SphP[idx].Utherm * mean_density * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);

#if CHEMISTRYNETWORK == 1
      double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP] + SphP[idx].TracAbund[IHEP] +
                      2. * SphP[idx].TracAbund[IHEPP]) *
                     yn;
#elif CHEMISTRYNETWORK == 15
      double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP] + SphP[idx].TracAbund[IHEP]) * yn;
#else
      double yntot = (1.0 + ABHE - SphP[idx].TracAbund[IH2] + SphP[idx].TracAbund[IHP]) * yn;
#endif

      double gamma_minus1            = 2.0 / 3.0;  // No H2, so gamma=5/3
      double gas_temp                = gamma_minus1 * gas_energy_density / (yntot * BOLTZMANN);
      double totalThermalEnergyAdded = 0.;

      if(gas_temp < temp_set)
        {
          gas_energy_density = (1. / gamma_minus1) * yntot * BOLTZMANN * temp_set;
          double new_utherm  = gas_energy_density * pow(All.UnitLength_in_cm, 3) / (mean_density * All.UnitEnergy_in_cgs);

          totalThermalEnergyAdded = (new_utherm - SphP[idx].Utherm) * total_mass;
          if(totalThermalEnergyAdded > 10. * 1E51 / All.UnitEnergy_in_cgs)
            warn(
                "SNe: during an unresolved SN event we just injected considerably more thermal energy than expected: E_therm = %g "
                "erg\n",
                totalThermalEnergyAdded * All.UnitEnergy_in_cgs);

          SphP[idx].Utherm = new_utherm;
#ifdef MHD_THERMAL_ENERGY_SWITCH
          SphP[idx].Etherm = SphP[idex].Utherm * P[idx].Mass;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          // Update Entropy, Taken from SGChem
          SphP[idx].A       = gamma_minus1 * SphP[idx].Utherm / pow(SphP[idx].Density * All.cf_a3inv, GAMMA_MINUS1);
          SphP[idx].Entropy = log(SphP[idx].A) * P[idx].Mass;
#endif
        }
#endif /*SGCHEM*/  // TODO:without SGCHEM momentum injection will be wrong or not work at all
      // Total energy (account for changes in thermal, kinetic energy)
      SphP[idx].Energy += totalThermalEnergyAdded + dEkin;
    }
}
