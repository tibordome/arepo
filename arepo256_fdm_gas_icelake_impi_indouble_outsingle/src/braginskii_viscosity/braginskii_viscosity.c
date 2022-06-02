/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/braginskii_viscosity/braginskii_viscosity.c
 * \date        11/2020
 * \author      Thomas Berlok
 * \brief       Implementation of Braginskii viscosity.
 * \details     This implementation of Braginskii viscosity uses some of the
 *              infrastructure developed for cosmic ray diffusion by Ruediger
 *              Pakmor, see src/cosmic_rays_diffusion.c
 *
 *              A description of this module can be found in the paper
 *
 *              Braginskii viscosity on an unstructured, moving mesh -- accelerated
 *              with super-time-stepping, T. Berlok, R. Pakmor & and C. Pfrommer.
 *
 *              https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2919B/abstract
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"
#include "sys/stat.h"

#ifdef BRAGINSKII_VISCOSITY

#include <gsl/gsl_linalg.h>
#include "braginskii_viscosity.h"

/* structs first ... */
static struct diff_face_data *diff_face_data;

static struct corner_list *corner_list;

// Struct for storing values at corners which only need to be calculated once.
static struct corner_data
{
  int active;
  // The M matrix in Pakmor et al. 2016 (eq. 8)
  double matrix[NUMDIMS + 1][NUMDIMS + 1];

  // Magnetic field vector
  double bfld[3];

  // mass density
  // double density;

  // // pressure
  // double pressure;

  // eta_aniso
  double eta_aniso;

  int tetra;
  int fail;
} * corner_data;

static struct corner_vel *corner_vel;

// Struct for storing the velocity at the center of cells
#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
static struct center_vel *center_vel;
#endif

// Struct for communicating viscous fluxes between tasks
static struct viscflux_list_data *ViscFluxList;

// Struct for storing a copy of the velocities.
static struct sphp_copy *state0, *state1Exch;
// For super time stepping, we need to have copies of the velocity at several
// substeps between each super step.
#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
static struct sphp_copy *state_j, *state_jm1;  //, *state_jm2;
#endif

//
static struct grad_v *grad_v;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
static double *refCenter;
static int *pIsReflective;
#endif

static int cornerCount;

static int get_substeps_of_point(int point);

/*
  Initialize the viscosity coefficient and perform some consistency checks on the Config
  options
*/
void init_braginskii_viscosity(void)
{
  double L0 = All.UnitLength_in_cm;
  double m0 = All.UnitMass_in_g;
  double t0 = All.UnitTime_in_s;

  // TODO: Check Hubble constant scaling here
  if(All.HubbleParam > 0.)
    {
      L0 /= All.HubbleParam;
      m0 /= All.HubbleParam;
      t0 /= All.HubbleParam;
    }

#ifdef BRAGINSKII_SPITZER
  double coulomb_log = 30.0;

  double meanweight = PROTONMASS * 4.0 / (8.0 - 5.0 * (1 - HYDROGEN_MASSFRAC)); /* assuming full ionization */
  double nu0        = 0.96 * 0.75 / pow(ELECTRONCHARGE, 4.0) / coulomb_log * sqrt(PROTONMASS / M_PI) * pow(meanweight, 2.5);

  // Make dimensionless
  nu0 *= pow(L0, 6.0) / m0 / pow(t0, 4.0);

  /* Store. Note that nu_0 is just a convenient constant, not a viscosity
    coefficient. The viscosity cofficient is:
    \nu_\para = nu_0 \cdot (P/\rho)^{5/2}/\rho
    Here P and rho are in code units, i.e., All.UnitPressure_in_cgs and
    All.UnitPressure_in_cgs, which are given by
    All.UnitDensity_in_cgs = m0 / pow(L0, 3);
    All.UnitPressure_in_cgs = m0 / L0 / pow(t0, 2);
  */
  All.BragViscosityCoefficient = nu0;

  // Convert maximum coefficient to internal units
  All.BragViscosityMaximumCoefficient *= t0 / (L0 * L0);

  mpi_printf("Braginskii: BragViscosityMaximumCoefficient set to %g (int units)\n", All.BragViscosityMaximumCoefficient);

#else

  // Convert to internal units
  All.BragViscosityCoefficient *= t0 / (L0 * L0);

#endif
  mpi_printf("Braginskii: nu_0 set to %g (int units)\n", All.BragViscosityCoefficient);
}

/*!
 * Main function for the implementation of Braginskii viscosity via operator splitting.
 * Called in src/run.c
 */
void do_braginskii_viscosity(void)
{
  TIMER_START(CPU_BRAGINSKII);

#ifdef BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
  if(All.brag_diffusion_Ti_endstep != All.Ti_Current && All.brag_diffusion_Ti_begstep != All.Ti_Current)
    {
      mpi_printf("BRAGINSKII_VISCOSITY: Nothing to do on this timestep, begin=%d, end=%d, current=%d\n", All.brag_diffusion_Ti_begstep,
                 All.brag_diffusion_Ti_endstep, All.Ti_Current);
      return;
    }
#endif

  int CountAll;
  MPI_Allreduce(&TimeBinsHydro.NActiveParticles, &CountAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(CountAll == 0)
    {
      mpi_printf("BRAGINSKII_VISCOSITY: No gas cells active, skipping.\n");
      TIMER_STOP(CPU_BRAGINSKII);
      return;
    }
  else
    mpi_printf("BRAGINSKII_VISCOSITY: Doing diffusion for %d active cells.\n", CountAll);

  // Update and exchange primitive variables
  update_primitive_variables();

  TIMER_START(CPU_BRAGINSKII_EXCHANGE);
  exchange_primitive_variables();
  TIMER_STOP(CPU_BRAGINSKII_EXCHANGE);

  int substeps = prepare_stuff();
  set_diffusion_coefficients(diff_face_data);

// RKL2 super timestepping
#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
  diffusion_explicit_sts(substeps);
#else
  // Standard explicit update (perhaps with subcycling)
  diffusion_explicit(substeps);
#endif

  free_stuff();

  update_primitive_variables();
  TIMER_STOP(CPU_BRAGINSKII);
}

/* Explicit version of the diffusion */
void diffusion_explicit(int substeps)
{
  mpi_printf("BRAGINSKI_VISCOSITY: Doing explicit diffusion.\n");

  // Allocate memory
  state0 = (struct sphp_copy *)mymalloc_movable(&state0, "state0", NumGas * sizeof(struct sphp_copy));

  // Copy necessary information from SphP and P.
  copy_SphP(state0);

  // Update the state

#ifndef BRAGINSKII_VISCOSITY_SUBCYCLE
  calculate_and_add_fluxes(state0, state0, 1.0, 1.0);
#endif

#ifdef BRAGINSKII_VISCOSITY_SUBCYCLE
// Determine the number of subcycles
#ifdef BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
// int subcycles = All.BragViscositySubcycles;
#else
  MPI_Allreduce(&substeps, &All.BragViscositySubcycles, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
// int subcycles = substeps;
#endif

  int subcycles = All.BragViscositySubcycles;

  // Subcycle explicit update
  mpi_printf("BRAGINSKI_VISCOSITY: Subcycling with %d.\n", subcycles);
  for(int j = 0; j < subcycles; j++)
    {
      calculate_and_add_fluxes(state0, state0, 1.0 / subcycles, 1.0 / subcycles);
    }
#endif

  // Copy result back to main array
  set_SphP(state0);

  // Free memory
  myfree(state0);
}

/*
 Super time stepping (STS) version of Braginskii viscosity.
 The implementation is based on "A second-order accurate Super TimeStepping
 formulation for anisotropic thermal conduction" by Meyer, Balsara & Aslam 2012
 (see doi:10.1111/j.1365-2966.2012.20744.x). More specifically, this function
 updates the velocities using the second-order accurate Runge-Kutta Legendre
 method (RKL2). See especially equations 14-17 in Meyer et al. (2012).
 See section 3.2 in the Braginskii code paper.
*/
#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
void diffusion_explicit_sts(int substeps)
{
  mpi_printf("BRAGINSKI_VISCOSITY: Doing diffusion with RKL2 STS.\n");

  // Determine the number of stages in the super timestepping
#ifndef BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
  int stages = imax(substeps, 3);
  MPI_Allreduce(&stages, &All.BragViscosityRKL2Stages, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  int s = All.BragViscosityRKL2Stages;  // shorthand

  mpi_printf("BRAGINSKI_VISCOSITY: We have %d RKL2 stages.\n", s);

  // Allocate memory
  // state0 will contain τ L(Y⁰) which we calculate only once.
  state0 = (struct sphp_copy *)mymalloc_movable(&state0, "state0", NumGas * sizeof(struct sphp_copy));
  // System states at various steps (jm1 is j minus 1 etc)
  // state_jm2 = mymalloc_movable(&state_jm2, "state_jm2", NumGas * sizeof(struct sphp_copy));
  state_jm1 = (struct sphp_copy *)mymalloc_movable(&state_jm1, "state_jm1", NumGas * sizeof(struct sphp_copy));
  state_j   = (struct sphp_copy *)mymalloc_movable(&state_j, "state_j", NumGas * sizeof(struct sphp_copy));

  // STS weights and parameters
  double mu[s + 1];
  double mu_til[s + 1];
  double nu[s + 1];
  double gam_til[s + 1];
  double w1 = 4.0 / (s * s + s - 2.0);

  find_rkl_weights(mu, mu_til, nu, gam_til, s);

  // tau = dt*(s*s + s - 2.0)/4.0 is the super time step
  // tau/dt set to one, as dt is set to tau inside timestep.c
  double tau_over_dt = 1.0;

  /*
  The implementation commented out below is a conceptually simpler but
  more memory intensive version of the one actually used.
  The uncommented version has eliminated state_jm2

  // Copy necessary information from SphP and P, initialize some energies to zero
  // TODO: Rename state.Energy to state.delta_Energy for clarity?
  // Would then need to modify non-RKL2 as well.
  for(int part = 0; part < NumGas; part++)
    {
       // Initial state at j-2
      state_jm2[part].Energy      = SphP[part].Energy;
      state_jm2[part].Momentum[0] = SphP[part].Momentum[0];
      state_jm2[part].Momentum[1] = SphP[part].Momentum[1];
      state_jm2[part].Momentum[2] = SphP[part].Momentum[2];

      // Initial state at j-1
      state_jm1[part].Momentum[0] = SphP[part].Momentum[0];
      state_jm1[part].Momentum[1] = SphP[part].Momentum[1];
      state_jm1[part].Momentum[2] = SphP[part].Momentum[2];

      // Set energies to zero. We collect the total change in each cell
      // by adding up the changes during the RKL2 stages
      state_jm1[part].Energy = 0.0;
    }

  // Solution at j-1 (that is, at dt/τ = w_1/3.0) now contained in state_jm1
  calculate_and_add_fluxes(state_jm1, state_jm1, mu_til[1] * tau_over_dt, w1/3.0);


  for(int part = 0; part < NumGas; part++)
    {
      // Overwrite state0 with τ L(Y⁰)
      state0[part].Momentum[0] = (state_jm1[part].Momentum[0] - SphP[part].Momentum[0])/mu_til[1];
      state0[part].Momentum[1] = (state_jm1[part].Momentum[1] - SphP[part].Momentum[1])/mu_til[1];
      state0[part].Momentum[2] = (state_jm1[part].Momentum[2] - SphP[part].Momentum[2])/mu_til[1];

      // Copy energy change from first step to state_j
      state_j[part].Energy = state_jm1[part].Energy;
    }

  for(int j = 2; j < s + 1; j++)
    {
      // All terms on the RHS of equation 59 in the Braginskii code paper
      // except for the one proportional to \tilde(μ)_j are
      // stored in state_j
      for(int part = 0; part < NumGas; part++)
        {
          state_j[part].Momentum[0] = nu[j] * state_jm2[part].Momentum[0] +
                                      mu[j] * state_jm1[part].Momentum[0] +
                                 gam_til[j] * state0   [part].Momentum[0] +
                          (1.0 - mu[j] - nu[j]) *  SphP[part].Momentum[0];

          state_j[part].Momentum[1] = nu[j] * state_jm2[part].Momentum[1] +
                                      mu[j] * state_jm1[part].Momentum[1] +
                                 gam_til[j] * state0   [part].Momentum[1] +
                          (1.0 - mu[j] - nu[j]) *  SphP[part].Momentum[1];

          state_j[part].Momentum[2] = nu[j] * state_jm2[part].Momentum[2] +
                                      mu[j] * state_jm1[part].Momentum[2] +
                                 gam_til[j] * state0   [part].Momentum[2] +
                          (1.0 - mu[j] - nu[j]) *  SphP[part].Momentum[2];

        }

      // factorE is the time-advance during this substep (in units of dt/τ, hardcoded to 1)
      double factorE;
      if (j == 2)
        factorE = w1*2.0/3.0;
      else
        factorE = w1*j/2.0;

      // the last term in eq 59, \tilde(μ)_j τ L(Y_(j-1)), is added to state_j
      calculate_and_add_fluxes(state_jm1, state_j, mu_til[j] * tau_over_dt, factorE);

      // Copy from left to right (velocities only)
      for(int part = 0; part < NumGas; part++)
      {
        state_jm2[part].Momentum[0] = state_jm1[part].Momentum[0];
        state_jm2[part].Momentum[1] = state_jm1[part].Momentum[1];
        state_jm2[part].Momentum[2] = state_jm1[part].Momentum[2];

        state_jm1[part].Momentum[0] = state_j[part].Momentum[0];
        state_jm1[part].Momentum[1] = state_j[part].Momentum[1];
        state_jm1[part].Momentum[2] = state_j[part].Momentum[2];

      }

    }
  // state_j now has the solution for the velocities at t + tau and the
  // accumulated energy change.

  end of conceptually simpler but more mememory intensive version
  */

  // Copy necessary information from SphP and P, initialize some energies to zero
  // TODO: Rename state.Energy to state.delta_Energy for clarity?
  // Would then need to modify non-RKL2 as well.
  for(int part = 0; part < NumGas; part++)
    {
      // Initial state at j-2
      state_jm1[part].Momentum[0] = SphP[part].Momentum[0];
      state_jm1[part].Momentum[1] = SphP[part].Momentum[1];
      state_jm1[part].Momentum[2] = SphP[part].Momentum[2];

      // Set energies to zero. We collect the total change in each cell
      // by adding up the changes during the RKL2 stages
      state_jm1[part].Energy = 0.0;

      // Add contribution from j-2 to j (eq 59)
      state_j[part].Momentum[0] = nu[2] * state_jm1[part].Momentum[0];
      state_j[part].Momentum[1] = nu[2] * state_jm1[part].Momentum[1];
      state_j[part].Momentum[2] = nu[2] * state_jm1[part].Momentum[2];
    }

  // Solution at j-1 (that is, at dt/τ = w_1/3.0) now contained in state_jm1
  calculate_and_add_fluxes(state_jm1, state_jm1, mu_til[1] * tau_over_dt, w1 / 3.0);

  for(int part = 0; part < NumGas; part++)
    {
      // Overwrite state0 with τ L(Y⁰)
      state0[part].Momentum[0] = (state_jm1[part].Momentum[0] - SphP[part].Momentum[0]) / mu_til[1];
      state0[part].Momentum[1] = (state_jm1[part].Momentum[1] - SphP[part].Momentum[1]) / mu_til[1];
      state0[part].Momentum[2] = (state_jm1[part].Momentum[2] - SphP[part].Momentum[2]) / mu_til[1];

      // Copy energy change from first step to state_j
      state_j[part].Energy = state_jm1[part].Energy;

      // Add contribution from j-1 to j (eq 59)
      state_j[part].Momentum[0] += mu[2] * state_jm1[part].Momentum[0];
      state_j[part].Momentum[1] += mu[2] * state_jm1[part].Momentum[1];
      state_j[part].Momentum[2] += mu[2] * state_jm1[part].Momentum[2];

      // Add contribution from initial condition and L(Y⁰) to state j (eq 59)
      double zero_contr = 1.0 - mu[2] - nu[2];
      state_j[part].Momentum[0] += gam_til[2] * state0[part].Momentum[0] + zero_contr * SphP[part].Momentum[0];

      state_j[part].Momentum[1] += gam_til[2] * state0[part].Momentum[1] + zero_contr * SphP[part].Momentum[1];

      state_j[part].Momentum[2] += gam_til[2] * state0[part].Momentum[2] + zero_contr * SphP[part].Momentum[2];
    }

  for(int j = 2; j < s + 1; j++)
    {
      // factorE is the time-advance during this substep (in units of dt/τ, hardcoded to 1)
      double factorE;
      if(j == 2)
        factorE = w1 * 2.0 / 3.0;
      else
        factorE = w1 * j / 2.0;

      // the last term in eq 59, \tilde(μ)_j τ L(Y_(j-1)), is added to state_j
      calculate_and_add_fluxes(state_jm1, state_j, mu_til[j] * tau_over_dt, factorE);

      // Preparation for next step. Break out if this was the final step.
      if(j >= s)
        break;
      for(int part = 0; part < NumGas; part++)
        {
          // Store values in state_j[part] (this allowed us to get rid of state_jm2,
          // compare with commented code above)
          // these will be the new state_jm1
          double momentum[3];
          momentum[0] = state_j[part].Momentum[0];
          momentum[1] = state_j[part].Momentum[1];
          momentum[2] = state_j[part].Momentum[2];

          // Add contribution from j-2 to j (eq 59)
          state_j[part].Momentum[0] = nu[j + 1] * state_jm1[part].Momentum[0];
          state_j[part].Momentum[1] = nu[j + 1] * state_jm1[part].Momentum[1];
          state_j[part].Momentum[2] = nu[j + 1] * state_jm1[part].Momentum[2];

          // Add contribution from j-1 to j (eq 59)
          state_j[part].Momentum[0] += mu[j + 1] * momentum[0];
          state_j[part].Momentum[1] += mu[j + 1] * momentum[1];
          state_j[part].Momentum[2] += mu[j + 1] * momentum[2];

          // Add contribution from initial condition and L(Y⁰) to state j (eq 59)
          double zero_contr = 1.0 - mu[j + 1] - nu[j + 1];
          state_j[part].Momentum[0] += gam_til[j + 1] * state0[part].Momentum[0] + zero_contr * SphP[part].Momentum[0];

          state_j[part].Momentum[1] += gam_til[j + 1] * state0[part].Momentum[1] + zero_contr * SphP[part].Momentum[1];

          state_j[part].Momentum[2] += gam_til[j + 1] * state0[part].Momentum[2] + zero_contr * SphP[part].Momentum[2];

          // Set state_jm1 for next step
          state_jm1[part].Momentum[0] = momentum[0];
          state_jm1[part].Momentum[1] = momentum[1];
          state_jm1[part].Momentum[2] = momentum[2];
        }
    }
  // state_j now has the solution for the velocities at t + tau and the
  // accumulated energy change.

  // Add/subtract accumulated energy change and copy the velocities back to the
  // main array
  for(int part = 0; part < NumGas; part++)
    {
      SphP[part].Energy += state_j[part].Energy;
      SphP[part].Momentum[0] = state_j[part].Momentum[0];
      SphP[part].Momentum[1] = state_j[part].Momentum[1];
      SphP[part].Momentum[2] = state_j[part].Momentum[2];
    }

  // Free memory
  myfree(state_j);
  myfree(state_jm1);
  // myfree(state_jm2);
  myfree(state0);

  mpi_printf("BRAGINSKI_VISCOSITY: Finished diffusion with RKL2 STS.\n");
}
#endif

/*
 Calculates the flux of state1 and applies this flux (multiplied by the variable
 factor) to state2. For an explicit update, factor = 1 and state1=state2, but
 this is not the case when using the STS method.
 */
void calculate_and_add_fluxes(struct sphp_copy *state1, struct sphp_copy *state2, double factor, double factorE)
{
  TIMER_START(CPU_BRAGINSKII_CALCULATE_AND_ADD_FLUXES);
  // Compute vx, vy, vz and their gradients at corners
  state1Exch = (struct sphp_copy *)mymalloc("state1Exch", Mesh_nimport * sizeof(struct sphp_copy));
  exchange_state_variables(state1, state1Exch);
  compute_quantities_at_corners(state1, state1Exch);

#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
  center_vel = (struct center_vel *)mymalloc("center_vel", Mesh.Ndp * sizeof(struct center_vel));
  compute_velocity_at_centers(state1, state1Exch);
#else
  myfree(state1Exch);
#endif

  int Nflux    = 0;
  int MaxNflux = Mesh.Indi.AllocFacNflux;
  ViscFluxList =
      (struct viscflux_list_data *)mymalloc_movable(&ViscFluxList, "FluxList", MaxNflux * sizeof(struct viscflux_list_data));
  grad_v = (struct grad_v *)mymalloc("grad_v", sizeof(struct grad_v));
  double flux_E, flux_Mx, flux_My, flux_Mz;

  if(All.ComovingIntegrationOn)
    factor *= All.cf_atime * All.cf_atime;
  factorE *= All.cf_atime * All.cf_atime * All.cf_atime;

  // Loop over Voronoi interfaces
  for(int iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      if(!fd->active || fd->dt == 0.)
        continue;

      // Estimate velocity gradients at interface
      calculate_grad_v_in_local_coordinates(fd, grad_v);
#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
      calculate_simple_normal_gradients(fd, iface, grad_v);
#endif
      // Calculate flux across interface
      flux_E  = 0.0;
      flux_Mx = 0.0;
      flux_My = 0.0;
      flux_Mz = 0.0;
      if(fd->cornerCount - fd->failCount >= 2)
        {
#ifdef BRAGINSKII_ISOTROPIC
          calculate_isotropic_flux_across_face(fd, &flux_E, &flux_Mx, &flux_My, &flux_Mz, factor, factorE);
#else
          calculate_flux_across_face(fd, grad_v, &flux_E, &flux_Mx, &flux_My, &flux_Mz, factor, factorE);
#endif
        }

      // Directly apply flux or add to fluxlist for later communication
      Nflux = apply_flux_or_calculate_fluxlist(ViscFluxList, state2, fd, flux_E, flux_Mx, flux_My, flux_Mz, iface, &MaxNflux, Nflux);
    }

  // Communicate fluxlist between tasks.
  exchange_and_apply_fluxes(ViscFluxList, Nflux, state2);

  // Free memory
  myfree(grad_v);
  myfree(ViscFluxList);
#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
  myfree(center_vel);
  myfree(state1Exch);
#endif
  TIMER_STOP(CPU_BRAGINSKII_CALCULATE_AND_ADD_FLUXES);
}

/*
 This function calculates the flux of velocity and energy across the face
 fd by using precomuted gradients of the velocity (stored in grad_v).
 The resulting fluxes are stored in flux_E, flux_Mx, flux_My and flux_Mz.
 */
void calculate_flux_across_face(struct diff_face_data *fd, struct grad_v *grad_v, double *flux_E, double *flux_Mx, double *flux_My,
                                double *flux_Mz, double factor, double factorE)
{
#ifdef BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
  double dt = (All.brag_diffusion_Ti_endstep - All.brag_diffusion_Ti_begstep) * All.Timebase_interval;
#else
  double dt = fd->dt;
#endif

  // Density
  // double rho = fd->density;

  // Magnetic field direction vector in local coordinate system
  double bn = fd->bfld[0];
  double bm = fd->bfld[1];

  // Components of the unit vectors of the local coordinate system at the
  // interface. n is normal to the interface and m is parallel to the
  // the magnetic field component which lies in the plane of the interface.
  double nx = fd->nx;
  double mx = fd->mx;
  double ny = fd->ny;
  double my = fd->my;
  double nz = fd->nz;
  double mz = fd->mz;

  // Components of grad v in the local coordinate system
  double vn = grad_v->vn;
  double vm = grad_v->vm;

  // Components of v in the local coordinate system
  double dvndn = grad_v->dvndn;
  double dvmdn = grad_v->dvmdn;
  double dvndm = grad_v->dvndm;
  double dvmdm = grad_v->dvmdm;
  double dvpdp = grad_v->dvpdp;

  // Trace of the dyadic of b with the gradient of v (\b\b : \nabla \vec{v})
  double bbdV = bn * bn * dvndn + bm * bm * dvmdm + bn * bm * (dvmdn + dvndm);

  // Divergence of the velocity, v. (\nabla \cdot \vec{v})
  double divV = dvndn + dvmdm + dvpdp;

  // Pressure anisotropy
  double presaniso = fd->eta_aniso * (3.0 * bbdV - divV);

  // We can limit pressure anisotropy according to stability criterion for
  // firehose and mirror instability.
#ifdef BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY
  if(presaniso < -fd->B2)
    {
      presaniso = -fd->B2;
    }
  if(presaniso > 0.5 * fd->B2)
    {
      presaniso = 0.5 * fd->B2;
    }
#endif

  // Calculate normal and tangential velocity flux
  double flux_Mn = 0.5 * dt * factor * presaniso * (bn * bn - 1.0 / 3.0);
  double flux_Mm = 0.5 * dt * factor * presaniso * bn * bm;

  // Convert fluxes to Cartesian system
  *flux_Mx = flux_Mn * nx + flux_Mm * mx;
  *flux_My = flux_Mn * ny + flux_Mm * my;
  *flux_Mz = flux_Mn * nz + flux_Mm * mz;

  // Calculate the energy flux
  *flux_E = 0.5 * dt * factorE * presaniso * (bn * bn * vn + bn * bm * vm - vn / 3.0);
}

/*
 Isotropic Navier-Stokes viscosity. Sometimes useful for comparison with
 anisotropic viscosity.
 This function calculates the flux of velocity and energy across the face
 fd. The resulting fluxes are added to flux_E, flux_Mx, flux_My and flux_Mz.
 */
#ifdef BRAGINSKII_ISOTROPIC
void calculate_isotropic_flux_across_face(struct diff_face_data *fd, double *flux_E, double *flux_Mx, double *flux_My, double *flux_Mz,
                                          double factor, double factorE)
{
#ifdef BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP
  double dt = (All.brag_diffusion_Ti_endstep - All.brag_diffusion_Ti_begstep) * All.Timebase_interval;
#else
  double dt                   = fd->dt;
#endif

  // Components of the unit vector normal to the interface interface.
  double nx = fd->nx;
  double ny = fd->ny;
  double nz = fd->nz;
  double vx = 0.0;
  double vy = 0.0;
  double vz = 0.0;

  // Loop over all corners
  for(int icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_vel *cv  = &corner_vel[cl->index];

      // Velocity in Cartesian system
      vx += cl->weight * cv->vx;
      vy += cl->weight * cv->vy;
      vz += cl->weight * cv->vz;
    }

  double dvxdx = limited_derivative_iso(fd, 1);
  double dvxdy = limited_derivative_iso(fd, 2);
  double dvxdz = limited_derivative_iso(fd, 3);
  double dvydx = limited_derivative_iso(fd, 4);
  double dvydy = limited_derivative_iso(fd, 5);
  double dvydz = limited_derivative_iso(fd, 6);
  double dvzdx = limited_derivative_iso(fd, 7);
  double dvzdy = limited_derivative_iso(fd, 8);
  double dvzdz = limited_derivative_iso(fd, 9);

  // Divergence of velocity
  double divV = dvxdx + dvydy + dvzdz;

  // TODO: Set coefficients in a better way and make more tests.
  // Also check that the density dependence is correct.
  // No bulk viscosity for the moment
  double zeta = 0.0;

  // Hardcode shear viscosity for the moment
  double eta = fd->eta_aniso;

  // Flux in x-direction
  *flux_Mx += dt * factor * 0.5 *
              (nx * ((zeta - 2.0 / 3.0 * eta) * divV + 2.0 * eta * dvxdx) + ny * eta * (dvxdy + dvydx) + nz * eta * (dvxdz + dvzdx));

  // Flux in y-direction
  *flux_My += dt * factor * 0.5 *
              (nx * eta * (dvxdy + dvydx) + ny * ((zeta - 2.0 / 3.0 * eta) * divV + 2.0 * eta * dvydy) + nz * eta * (dvydz + dvzdy));

  // Flux in z-direction
  *flux_Mz += dt * factor * 0.5 *
              (nx * eta * (dvxdz + dvzdx) + ny * eta * (dvydz + dvzdy) + nz * ((zeta - 2.0 / 3.0 * eta) * divV + 2.0 * eta * dvzdz));

  // Calculate the energy flux
  *flux_E += dt * factorE * 0.5 *
             ((vx * nx + vy * ny + vz * nz) * (zeta - 2.0 / 3.0 * eta) * divV +
              eta * (2.0 * dvxdx * nx * vx + (dvxdy + dvydx) * ny * vx + (dvxdz + dvzdx) * nz * vx + (dvxdy + dvydx) * nx * vy +
                     2.0 * dvydy * ny * vy + (dvydz + dvzdy) * nz * vy + (dvxdz + dvzdx) * nx * vz + (dvydz + dvzdy) * ny * vz +
                     2.0 * dvzdz * nz * vz));
}
#endif

/*!
 * Computes geometry of interfaces, taking into account boundary conditions.
 * Uses this information to calculate the matrix M (equation 8 in Pakmor et al 2016) in each cell.
 * The matrix M is used to compute vx, vy, vz and their gradients at corners.
 * The magnetic field at interfaces is also calculated.
 */
int prepare_stuff(void)
{
  TIMER_START(CPU_BRAGINSKII_PREPARE);
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  refCenter     = (double *)mymalloc("refCenter", Mesh.Ndp * 3 * sizeof(double));
  pIsReflective = (int *)mymalloc("pIsReflective", Mesh.Ndp * sizeof(int));
  memset(pIsReflective, 0, Mesh.Ndp * sizeof(int));
  // 0: no interesting boundary, 1: reflective boundary, 2: outflow boundary

  {
    // we go through all the interfaces, look for boundary points and calculate their mirrored center of mass
    int iface;
    for(iface = 0; iface < Mesh.Nvf; iface++)
      {
        int p1 = Mesh.VF[iface].p1;
        int p2 = Mesh.VF[iface].p2;

        if(p1 < 0 || p2 < 0)
          continue;

        if(Mesh.DP[p1].task == ThisTask && Mesh.DP[p2].task == ThisTask)
          {
            int part1 = Mesh.DP[p1].index;
            int part2 = Mesh.DP[p2].index;

            if((part1 >= NumGas || part2 >= NumGas) && (part1 < NumGas || part2 < NumGas))
              {
                // we found a boundary point
                if(part1 >= NumGas)
                  part1 -= NumGas;
                if(part2 >= NumGas)
                  part2 -= NumGas;

                if(Mesh.VF[iface].area < 1e-5 * fmin(SphP[part1].ActiveArea, SphP[part2].ActiveArea))
                  continue;

                if(P[part1].ID == P[part2].ID)
                  {
                    // it is a reflective of outflow boundary, i.e. not a periodic boundary
                    int pCell, pGhost;

                    if(Mesh.DP[p1].index >= NumGas)
                      {
                        pCell  = p2;
                        pGhost = p1;
                      }
                    else if(Mesh.DP[p2].index >= NumGas)
                      {
                        pCell  = p1;
                        pGhost = p2;
                      }
                    else
                      {
                        // just to get rid of the warning
                        pCell = pGhost = -1;
                      }

                    int particle = Mesh.DP[pCell].index;

                    // calculate normal vector of the interface
                    double nx = Mesh.DP[pGhost].x - P[particle].Pos[0];
                    double ny = Mesh.DP[pGhost].y - P[particle].Pos[1];
                    double nz = Mesh.DP[pGhost].z - P[particle].Pos[2];

                    // perpendicular on the surface
                    double nn = sqrt(nx * nx + ny * ny + nz * nz);
                    nx /= nn;
                    ny /= nn;
                    nz /= nn;
                    double fx = (SphP[particle].Center[0] - Mesh.VF[iface].cx);
                    double fy = (SphP[particle].Center[1] - Mesh.VF[iface].cy);
                    double fz = (SphP[particle].Center[2] - Mesh.VF[iface].cz);
                    double ff = (fx * nx + fy * ny + fz * nz);

                    double px = SphP[particle].Center[0] - ff * nx;
                    double py = SphP[particle].Center[1] - ff * ny;
                    double pz = SphP[particle].Center[2] - ff * nz;

                    refCenter[pGhost * 3 + 0] = 2. * px - SphP[particle].Center[0];
                    refCenter[pGhost * 3 + 1] = 2. * py - SphP[particle].Center[1];
                    refCenter[pGhost * 3 + 2] = 2. * pz - SphP[particle].Center[2];

                    int image_flags = Mesh.DP[pGhost].image_flags;
#if defined(REFLECTIVE_X)
                    if((image_flags & REFL_X_FLAGS) && (image_flags & OUTFLOW_X))
                      pIsReflective[pGhost] = 2;
                    else if(image_flags & REFL_X_FLAGS)
                      pIsReflective[pGhost] = 1;
#endif
#if defined(REFLECTIVE_Y)
                    if((image_flags & REFL_Y_FLAGS) && (image_flags & OUTFLOW_Y))
                      pIsReflective[pGhost] = 2;
                    else if(image_flags & REFL_Y_FLAGS)
                      pIsReflective[pGhost] = 1;
#endif
#if defined(REFLECTIVE_Z)
                    if((image_flags & REFL_Z_FLAGS) && (image_flags & OUTFLOW_Z))
                      pIsReflective[pGhost] = 2;
                    else if(image_flags & REFL_Z_FLAGS)
                      pIsReflective[pGhost] = 1;
#endif
                  }
              }
          }
      }
  }
#endif

  /* for all interfaces, make a list of all cells sharing a corner with the interface
    and save their center positions */

  diff_face_data = (struct diff_face_data *)mymalloc("diff_face_data", Mesh.Nvf * sizeof(struct diff_face_data));
  corner_list    = (struct corner_list *)mymalloc("corner_list", Mesh.Nvf * 20 * sizeof(struct corner_list));
  corner_data    = (struct corner_data *)mymalloc("corner_data", Mesh.Ndt * sizeof(struct corner_data));
  corner_vel     = (struct corner_vel *)mymalloc("corner_vel", Mesh.Ndt * sizeof(struct corner_vel));

  /* compute matrices and magnetic fields at corners */
  int icorner;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    corner_data[icorner].active = 0;

  /* flag interfaces that are needed and compute their timestep */
  int iface;
  int substeps = 1;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      fd->active                = 0;
      fd->sng                   = 0;

      int p1 = Mesh.VF[iface].p1;
      int p2 = Mesh.VF[iface].p2;

      if(p1 < 0 || p2 < 0)
        continue;

      int timebin1 = get_timebin_of_point(p1);
      int timebin2 = get_timebin_of_point(p2);

      if(!TimeBinSynchronized[timebin1] && !TimeBinSynchronized[timebin2])
        continue;

      if(Mesh.DP[p1].ID < Mesh.DP[p2].ID)
        {
          if(TimeBinSynchronized[timebin1])
            {
              /* lower ID is active, its task is responsible */
              if(Mesh.DP[p1].task != ThisTask || Mesh.DP[p1].index >= NumGas)
                continue;
            }
          else
            {
              /* only higher ID is active, its task is responsible */
              if(Mesh.DP[p2].task != ThisTask || Mesh.DP[p2].index >= NumGas)
                continue;
            }
        }
      else if(Mesh.DP[p1].ID > Mesh.DP[p2].ID)
        {
          if(TimeBinSynchronized[timebin2])
            {
              /* lower ID is active, its task is responsible */
              if(Mesh.DP[p2].task != ThisTask || Mesh.DP[p2].index >= NumGas)
                continue;
            }
          else
            {
              /* only higher ID is active, its task is responsible */
              if(Mesh.DP[p1].task != ThisTask || Mesh.DP[p1].index >= NumGas)
                continue;
            }
        }
      else
        {
          /* interface with at least one ghost point */

          if(Mesh.DP[p1].task != ThisTask && Mesh.DP[p2].task != ThisTask)
            continue;

          if(Mesh.DP[p1].task != Mesh.DP[p2].task)
            terminate("This should not happen, I think...");

          if(Mesh.DP[p1].index >= NumGas && Mesh.DP[p2].index >= NumGas)
            continue;

          int p;
          if(Mesh.DP[p1].index < NumGas)
            p = p1;
          else
            p = p2;

          if(!TimeBinSynchronized[P[Mesh.DP[p].index].TimeBinHydro])
            continue;
        }

      double surfacearea = fmax(get_surface_area_of_cell(p1), get_surface_area_of_cell(p2));
      if(Mesh.VF[iface].area <= 1e-5 * surfacearea)
        continue;

      /* if we make it here, the interface needs to be included and we are responsible for it */
      fd->active = 1;

      /* now lets get its timestep */
      int timeBin = imin(timebin1, timebin2);
      if(timeBin == 0)
        {
          fd->dt = 0.;
        }
      else
        {
          fd->dt = (((integertime)1) << timeBin) * All.Timebase_interval / All.cf_hubble_a;
#if defined(BRAGINSKII_VISCOSITY_SUBCYCLE) || defined(BRAGINSKII_RKL2_SUPER_TIME_STEPPING)
          if(timeBin == timebin1)
            {
              int substeps1 = get_substeps_of_point(p1);
              substeps      = imax(substeps1, substeps);
            }
          else
            {
              int substeps2 = get_substeps_of_point(p2);
              substeps      = imax(substeps2, substeps);
            }
#endif
        }
    }

  cornerCount = 0;

#ifndef ONEDIMS
#ifdef TWODIMS
  const int edge_start[3] = {1, 2, 0};
  const int edge_end[3]   = {2, 0, 1};
#else
  const int edge_start[6]     = {0, 0, 0, 1, 1, 2};
  const int edge_end[6]       = {1, 2, 3, 2, 3, 3};
  const int edge_opposite[6]  = {3, 1, 2, 3, 0, 1};
  const int edge_nexttetra[6] = {2, 3, 1, 0, 2, 0};
#endif
#endif

  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];
      int p1                    = Mesh.VF[iface].p1;
      int p2                    = Mesh.VF[iface].p2;

      if(p1 < 0 || p2 < 0 || fd->active == 0)
        {
          fd->area   = 0;
          fd->active = 0;
          continue;
        }

      fd->area        = Mesh.VF[iface].area;
      fd->cornerFirst = -1;
      fd->cornerCount = 0;

#ifdef ONEDIMS
      terminate("not implemented at the moment.");
#else
#ifdef TWODIMS
      int tt   = Mesh.VF[iface].dt_index;
      tetra *t = &Mesh.DT[tt];

      int nr;
      for(nr = 0; nr < 3; nr++)
        {
          int start_index = t->p[edge_start[nr]];
          int end_index   = t->p[edge_end[nr]];

          if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
            break;
        }

      int qq = t->t[nr];

      fd->cornerFirst = cornerCount;

      add_corner(tt, iface);
      add_corner(qq, iface);

      for(icorner = fd->cornerFirst; icorner < fd->cornerFirst + fd->cornerCount; icorner++)
        corner_list[icorner].weight = 0.5;
#else /* 3D */
      int tt   = Mesh.VF[iface].dt_index;
      tetra *t = &Mesh.DT[tt];

      fd->cornerFirst = cornerCount;

      int nr;
      for(nr = 0; nr < 6; nr++)
        {
          int start_index = t->p[edge_start[nr]];
          int end_index   = t->p[edge_end[nr]];

          if((start_index == p1 && end_index == p2) || (start_index == p2 && end_index == p1))
            break;
        }

      tetra *prev, *next;
      int i, j, k, l, m, ii, jj, kk, ll, nn;

      i = edge_start[nr];
      j = edge_end[nr];
      k = edge_opposite[nr];
      l = edge_nexttetra[nr];

      prev = t;

      do
        {
          nn   = prev->t[l];
          next = &Mesh.DT[nn];

          add_corner(nn, iface);

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(next->p[m] == prev->p[k])
                ll = m;
              if(next->p[m] == prev->p[i])
                ii = m;
              if(next->p[m] == prev->p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("inconsistency");

          kk = 6 - (ll + ii + jj);

          prev = next;

          i = ii;
          l = ll;
          j = jj;
          k = kk;
        }
      while(next != t);

      /* fix weights */
#endif
#endif
    }

  // Compute the matrix M (eq. 8 in Pakmor et al 2016)
  compute_least_squares_matrix_at_corners();

  int failedInterfaces = 0;
  int goodInterfaces   = 0;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];

      compute_geometry_of_interface(iface);
      compute_magnetic_field_at_interface(iface);

      if(fd->active)
        {
          // if (fd->failCount == fd->cornerCount){
          if(fd->cornerCount - fd->failCount < 2)
            {
              failedInterfaces++;
            }
          else
            {
              goodInterfaces++;
            }
        }
    }

  MPI_Barrier(MPI_COMM_WORLD);

  int failedInterfacesSum;
  int goodInterfacesSum;
  MPI_Reduce(&failedInterfaces, &failedInterfacesSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&goodInterfaces, &goodInterfacesSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("BRAGINSKII_VISCOSITY: %d/%d (%1.2f %%) interfaces had less than 2 good corners.\n", failedInterfacesSum,
             failedInterfacesSum + goodInterfacesSum, 100.0 * failedInterfacesSum / (failedInterfacesSum + goodInterfacesSum));
  mpi_printf("BRAGINSKII_VISCOSITY: Preparations done.\n");
  // terminate("");
  TIMER_STOP(CPU_BRAGINSKII_PREPARE);

  return substeps;
}

int get_timebin_of_point(int point)
{
  int particle = Mesh.DP[point].index;

  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return P[particle].TimeBinHydro;
    }
  else
    return PrimExch[particle].TimeBinHydro;
}

double get_surface_area_of_cell(int point)
{
  int particle = Mesh.DP[point].index;

  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return SphP[particle].ActiveArea;
    }
  else
    return PrimExch[particle].ActiveArea;
}

#if defined(BRAGINSKII_VISCOSITY_SUBCYCLE) || defined(BRAGINSKII_RKL2_SUPER_TIME_STEPPING)
int get_substeps_of_point(int point)
{
  int particle = Mesh.DP[point].index;

  if(Mesh.DP[point].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;
      return SphP[particle].BragViscositySubsteps;
    }
  else
    return 0;  // PrimExch[particle].BragViscositySubcycles;
}
#endif

void add_corner(int tt, int iface)
{
  if(cornerCount >= Mesh.Nvf * 20)
    terminate("urg");

  corner_list[cornerCount].index = tt;
  cornerCount++;

  corner_data[tt].active = 1;

  struct diff_face_data *fd = &diff_face_data[iface];
  fd->cornerCount++;
}

/*!
 * Calculates the matrix M (equation 8 in Pakmor et al 2016).
 */
void compute_least_squares_matrix_at_corners(void)
{
  TIMER_START(CPU_BRAGINSKII_LEAST_SQUARE_MATRIX);
  int failCount   = 0;
  int activeCount = 0;

  int totMoves = 0;

  int icorner;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    {
      struct corner_data *cd = &corner_data[icorner];

      cd->fail = 0;

      if(!cd->active)
        continue;

      activeCount++;

      cd->tetra = icorner;

      double weights[NUMDIMS + 1];
      int boundary = 0;

      int row, col, k;
      for(row = 0; row < NUMDIMS + 1; row++)
        {
          int point = Mesh.DT[cd->tetra].p[row];

          if(point < 0)
            {
              for(k = 0; k < NUMDIMS + 1; k++)
                for(col = 0; col < NUMDIMS + 1; col++)
                  cd->matrix[k][col] = 0;

              boundary = 1;
              break;
            }

          double Center[3];
          if(!TimeBinSynchronized[Mesh.DP[point].timebin])
            {
              Center[0] = Mesh.DP[point].x;
              Center[1] = Mesh.DP[point].y;
              Center[2] = Mesh.DP[point].z;
            }
          else
            {
              int particle = Mesh.DP[point].index;

              if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
                particle -= NumGas;

              if(Mesh.DP[point].task == ThisTask)
                {
                  for(k = 0; k < 3; k++)
                    Center[k] = SphP[particle].Center[k];
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                  if(pIsReflective[point])
                    for(k = 0; k < 3; k++)
                      Center[k] = refCenter[point * 3 + k];
#endif
                }
              else
                {
                  for(k = 0; k < 3; k++)
                    Center[k] = PrimExch[particle].Center[k];
                }
            }

          double weight = 1.;
          weights[row]  = weight;

          double dist[NUMDIMS];
          dist[0] = NEAREST_X(Center[0] - Mesh.DTC[icorner].cx);
#if NUMDIMS > 1
          dist[1] = NEAREST_Y(Center[1] - Mesh.DTC[icorner].cy);
#endif
#if NUMDIMS > 2
          dist[2] = NEAREST_Z(Center[2] - Mesh.DTC[icorner].cz);
#endif

          cd->matrix[row][0] = weight;
          for(col = 0; col < NUMDIMS; col++)
            cd->matrix[row][col + 1] = dist[col] * weight;
        }

      if(boundary)
        continue;

      double matrixT[(NUMDIMS + 1)][(NUMDIMS + 1)];
      double matrix[(NUMDIMS + 1) * (NUMDIMS + 1)];
      for(row = 0; row < NUMDIMS + 1; row++)
        for(col = 0; col < NUMDIMS + 1; col++)
          {
            matrixT[row][col] = cd->matrix[col][row];
            int idx           = row * (NUMDIMS + 1) + col;
            matrix[idx]       = 0;
            for(k = 0; k < NUMDIMS + 1; k++)
              matrix[idx] += cd->matrix[k][row] * cd->matrix[k][col];
          }

      int s;
      gsl_matrix_view m     = gsl_matrix_view_array(matrix, NUMDIMS + 1, NUMDIMS + 1);
      gsl_permutation *perm = gsl_permutation_alloc(NUMDIMS + 1);
      gsl_linalg_LU_decomp(&m.matrix, perm, &s);

      if(gsl_linalg_LU_det(&m.matrix, s) != 0.)
        {
          double matrix_inverse[(NUMDIMS + 1) * (NUMDIMS + 1)];
          gsl_matrix_view minv = gsl_matrix_view_array(matrix_inverse, NUMDIMS + 1, NUMDIMS + 1);
          gsl_linalg_LU_invert(&m.matrix, perm, &minv.matrix);
          gsl_permutation_free(perm);

          for(row = 0; row < NUMDIMS + 1; row++)
            for(col = 0; col < NUMDIMS + 1; col++)
              {
                cd->matrix[col][row] = 0;
                for(k = 0; k < NUMDIMS + 1; k++)
                  cd->matrix[col][row] += matrix_inverse[row * (NUMDIMS + 1) + k] * matrixT[k][col];

                cd->matrix[col][row] *= weights[row];
              }
        }
      else
        {
          for(row = 0; row < NUMDIMS + 1; row++)
            {
              cd->matrix[row][0] = 1. / (NUMDIMS + 1);
              for(col = 1; col < NUMDIMS + 1; col++)
                cd->matrix[row][col] = 0.;
            }
        }

      for(row = 0; row < NUMDIMS + 1; row++)
        {
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
          if(cd->matrix[row][0] < -0.01 || pIsReflective[Mesh.DT[cd->tetra].p[row]])
#else
          if(cd->matrix[row][0] < -0.01)
#endif
            {
              cd->fail = 1;
              break;
            }
        }

      TIMER_START(CPU_BRAGINSKII_LEAST_SQUARE_MATRIX_GET_VALUES);

      if(cd->fail)
        failCount++;

      for(k = 0; k < 3; k++)
        cd->bfld[k] = 0;

      // cd->density  = 0.0;
      // cd->pressure = 0.0;
      cd->eta_aniso = 0.0;

      double Bmin[3], Bmax[3];

      for(k = 0; k < 3; k++)
        {
          Bmin[k] = +MAX_DOUBLE_NUMBER;
          Bmax[k] = -MAX_DOUBLE_NUMBER;
        }

      for(int p = 0; p < NUMDIMS + 1; p++)
        {
          int point    = Mesh.DT[cd->tetra].p[p];
          int particle = Mesh.DP[point].index;

          if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
            particle -= NumGas;

          double *B, density, pressure, eta_aniso;
          if(Mesh.DP[point].task == ThisTask)
            {
              B        = SphP[particle].B;
              density  = SphP[particle].Density;
              pressure = SphP[particle].Pressure;
              // eta_aniso = SphP[particle].Density * SphP[particle].nu_aniso;
            }
          else
            {
              B        = PrimExch[particle].B;
              density  = PrimExch[particle].Density;
              pressure = PrimExch[particle].Pressure;
              // eta_aniso= PrimExch[particle].Density * PrimExch[particle].nu_aniso;
            }

// TODO: Here we compute nu_aniso (ν_∥). I would instead like to precompute this
// inside timestep.c and store them in e.g. SphP[particle].nu_aniso.
// This would ensure that the ν_∥ values used for the update
// are identical to the ones used for estimating the allowed time step.
#ifndef BRAGINSKII_SPITZER
          double brag_coeff = All.BragViscosityCoefficient;
#else
          // Spitzer coefficient
          double brag_coeff = All.BragViscosityCoefficient * pow(pressure / density, 2.5) / density;
          // Cap the value of the diffusion coefficient
          if(brag_coeff > All.BragViscosityMaximumCoefficient)
            brag_coeff = All.BragViscosityMaximumCoefficient;
#endif

          // eta at center of cell
          eta_aniso = brag_coeff * density;

          struct corner_data *cd = &corner_data[icorner];

          // cd->density += density * cd->matrix[p][0];
          // cd->pressure += pressure * cd->matrix[p][0];

          // eta_aniso at corner
          cd->eta_aniso += eta_aniso * cd->matrix[p][0];

          // B-field at corner
          for(k = 0; k < 3; k++)
            {
              cd->bfld[k] += B[k] * cd->matrix[p][0];

              if(B[k] < Bmin[k])
                Bmin[k] = B[k];
              if(B[k] > Bmax[k])
                Bmax[k] = B[k];
            }
        }

      for(k = 0; k < 3; k++)
        {
          if(cd->bfld[k] > Bmax[k])
            cd->bfld[k] = Bmax[k];
          if(cd->bfld[k] < Bmin[k])
            cd->bfld[k] = Bmin[k];
        }
      TIMER_STOP(CPU_BRAGINSKII_LEAST_SQUARE_MATRIX_GET_VALUES);
    }

  int failCountSum, activeCountSum, totMovesSum;
  MPI_Reduce(&failCount, &failCountSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&activeCount, &activeCountSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&totMoves, &totMovesSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("BRAGINSKII_VISCOSITY: %d/%d (%1.2f %%) active corners were bad, totMoves=%d.\n", failCountSum, activeCountSum,
             100.0 * failCountSum / activeCountSum, totMovesSum);
  TIMER_STOP(CPU_BRAGINSKII_LEAST_SQUARE_MATRIX);
}

static double get_vector_distance(double a[3], double b[3])
{
  double dx = a[0] - b[0];
  double dy = a[1] - b[1];
  double dz = a[2] - b[2];

  return sqrt(dx * dx + dy * dy + dz * dz);
}

static double get_area_triangle(double p1[3], double p2[3], double p3[3])
{
  double a = get_vector_distance(p1, p2);
  double b = get_vector_distance(p2, p3);
  double c = get_vector_distance(p3, p1);

  double s    = 0.5 * (a + b + c);
  double prod = s * (s - a) * (s - b) * (s - c);

  if(prod < 0.)
    return 0.;  // can be a tiny bit below zero for degenerate triangles because of roundoff errors
  else
    return sqrt(prod);
}

void compute_geometry_of_interface(int iface)
{
  struct diff_face_data *fd = &diff_face_data[iface];

  if(!fd->active)
    return;

  fd->failWeight = 0.0;
  fd->failCount  = 0;

  int p1 = Mesh.VF[iface].p1;
  int p2 = Mesh.VF[iface].p2;

  double nx = Mesh.DP[p2].x - Mesh.DP[p1].x;
  double ny = Mesh.DP[p2].y - Mesh.DP[p1].y;
  double nz = Mesh.DP[p2].z - Mesh.DP[p1].z;
  double nn = sqrt(nx * nx + ny * ny + nz * nz);

  nx /= nn;
  ny /= nn;
  nz /= nn;

  fd->nx = nx;
  fd->ny = ny;
  fd->nz = nz;

  // need an ortonormal basis
  if(fd->nx != 0.0 || fd->ny != 0.0)
    {
      fd->mx = -fd->ny;
      fd->my = fd->nx;
      fd->mz = 0.0;
    }
  else
    {
      fd->mx = 1.0;
      fd->my = 0.0;
      fd->mz = 0.0;
    }

  double mm = sqrt(fd->mx * fd->mx + fd->my * fd->my + fd->mz * fd->mz);
  fd->mx /= mm;
  fd->my /= mm;
  fd->mz /= mm;

  fd->px = fd->ny * fd->mz - fd->nz * fd->my;
  fd->py = fd->nz * fd->mx - fd->nx * fd->mz;
  fd->pz = fd->nx * fd->my - fd->ny * fd->mx;

#ifdef TWODIMS
  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      cl->weight             = 0.5;
    }
#else
  // compute weights to get value at center from values at corners
  double cold[3];
  cold[0] =
      0.5 * (Mesh.DTC[corner_list[fd->cornerFirst + fd->cornerCount - 1].index].cx + Mesh.DTC[corner_list[fd->cornerFirst].index].cx);
  cold[1] =
      0.5 * (Mesh.DTC[corner_list[fd->cornerFirst + fd->cornerCount - 1].index].cy + Mesh.DTC[corner_list[fd->cornerFirst].index].cy);
  cold[2] =
      0.5 * (Mesh.DTC[corner_list[fd->cornerFirst + fd->cornerCount - 1].index].cz + Mesh.DTC[corner_list[fd->cornerFirst].index].cz);

  double fc[3];
  fc[0] = Mesh.VF[iface].cx;
  fc[1] = Mesh.VF[iface].cy;
  fc[2] = Mesh.VF[iface].cz;

  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];

      double p[3];
      p[0] = Mesh.DTC[cl->index].cx;
      p[1] = Mesh.DTC[cl->index].cy;
      p[2] = Mesh.DTC[cl->index].cz;

      double area = get_area_triangle(cold, fc, p);

      int icurr = fd->cornerFirst + icorner;
      int inext = fd->cornerFirst + ((icorner + 1) % fd->cornerCount);

      cold[0] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cx + Mesh.DTC[corner_list[inext].index].cx);
      cold[1] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cy + Mesh.DTC[corner_list[inext].index].cy);
      cold[2] = 0.5 * (Mesh.DTC[corner_list[icurr].index].cz + Mesh.DTC[corner_list[inext].index].cz);

      area += get_area_triangle(cold, fc, p);
      cl->weight = area / Mesh.VF[iface].area;
    }

  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_data *cd = &corner_data[cl->index];

      if(cd->fail)
        {
          fd->failWeight += cl->weight;
          fd->failCount++;
        }
    }

  // mpi_printf("fails: %d, corners: %d", fd->failCount, fd->cornerCount);
  // if (fd->failCount == fd->cornerCount){
  // mpi_printf("fails: %d, corners: %d", fd->failCount, fd->cornerCount);
  // terminate( "all corners failed" );
  // }

  // Adjust weights of working corners to have a sum of 1
  // and set bad corners to have a weight of 0.
  if(fd->failCount < fd->cornerCount)
    {
      for(icorner = 0; icorner < fd->cornerCount; icorner++)
        {
          struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
          struct corner_data *cd = &corner_data[cl->index];

          if(cd->fail)
            {
              cl->weight = 0.0;
            }
          else
            {
              cl->weight *= 1.0 / (1.0 - fd->failWeight);
            }
        }
    }

  // Calculate magnetic field components, bm & bp, in (temporary) local coord
  // system
  double bm = 0.;
  double bp = 0.;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_data *cd = &corner_data[cl->index];

      bm += (cd->bfld[0] * fd->mx + cd->bfld[1] * fd->my + cd->bfld[2] * fd->mz) * cl->weight;
      bp += (cd->bfld[0] * fd->px + cd->bfld[1] * fd->py + cd->bfld[2] * fd->pz) * cl->weight;
    }

  // Calculate new unit vectors e_m and e_p such that we
  // orient e_m along the b-component in the interface
  double mx = bm * fd->mx + bp * fd->px;
  double my = bm * fd->my + bp * fd->py;
  double mz = bm * fd->mz + bp * fd->pz;

  // But only store this if e_m can be normalized
  mm = sqrt(mx * mx + my * my + mz * mz);
  if(mm > 0)
    {
      fd->mx = mx / mm;
      fd->my = my / mm;
      fd->mz = mz / mm;

      fd->px = fd->ny * fd->mz - fd->nz * fd->my;
      fd->py = fd->nz * fd->mx - fd->nx * fd->mz;
      fd->pz = fd->nx * fd->my - fd->ny * fd->mx;
    }
#endif
}

// Compute magnetic field and density at interface
void compute_magnetic_field_at_interface(int iface)
{
  struct diff_face_data *fd = &diff_face_data[iface];

  if(!fd->active)
    return;

  int k;
  for(k = 0; k < 3; k++)
    fd->bfld[k] = 0.0;

  // fd->density = 0.0;

  int icorner;
  for(icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_data *cd = &corner_data[cl->index];
      // b \cdot e_1
      fd->bfld[0] += (cd->bfld[0] * fd->nx + cd->bfld[1] * fd->ny + cd->bfld[2] * fd->nz) * cl->weight;
      // b \cdot e_2
      fd->bfld[1] += (cd->bfld[0] * fd->mx + cd->bfld[1] * fd->my + cd->bfld[2] * fd->mz) * cl->weight;
      // b \cdot e_3 = 0 by construction. No need to calculate it here.
      // fd->bfld[2] += 0.0; //(cd->bfld[0] * fd->px + cd->bfld[1] * fd->py + cd->bfld[2] * fd->pz) * cl->weight;
      // fd->density += cd->density * cl->weight;
    }

  double b = sqrt(fd->bfld[0] * fd->bfld[0] + fd->bfld[1] * fd->bfld[1] + fd->bfld[2] * fd->bfld[2]);

#ifdef BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY
  fd->B2 = b * b;
#endif

  // Normalize magnetic field vector
  if(b != 0.)
    for(k = 0; k < 3; k++)
      fd->bfld[k] /= b;
}

#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
void compute_velocity_at_centers(struct sphp_copy *state, struct sphp_copy *stateExch)
{
  int p;
  for(p = 0; p < Mesh.Ndp; p++)
    {
      int part = Mesh.DP[p].index;
      if(Mesh.DP[p].task == ThisTask)
        {
          if(part >= NumGas)
            part -= NumGas;

          center_vel[p].vx = state[part].Momentum[0] / P[part].Mass;
          center_vel[p].vy = state[part].Momentum[1] / P[part].Mass;
          center_vel[p].vz = state[part].Momentum[2] / P[part].Mass;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
          if(pIsReflective[p])
            terminate("Reflective not implemented yet...");
#endif
        }
      else
        {
          double mass      = PrimExch[part].Volume * PrimExch[part].Density;
          center_vel[p].vx = stateExch[part].Momentum[0] / mass;
          center_vel[p].vy = stateExch[part].Momentum[1] / mass;
          center_vel[p].vz = stateExch[part].Momentum[2] / mass;
        }
    }
}
#endif

void compute_quantities_at_corners(struct sphp_copy *sphp_copy, struct sphp_copy *sphp_copyExch)
{
  int icorner, k;
  for(icorner = 0; icorner < Mesh.Ndt; icorner++)
    {
      struct corner_data *cd = &corner_data[icorner];
      struct corner_vel *cv  = &corner_vel[icorner];

      if(!cd->active)
        continue;

      cv->vx = 0.0;
      cv->vy = 0.0;
      cv->vz = 0.0;
      for(k = 0; k < NUMDIMS; k++)
        {
          cv->vxGrad[k] = 0.0;
          cv->vyGrad[k] = 0.0;
          cv->vzGrad[k] = 0.0;
        }

      tetra *tt = &Mesh.DT[cd->tetra];
      int p;
      for(p = 0; p < NUMDIMS + 1; p++)
        {
          int point = tt->p[p];

          if(point < 0)
            break;

          int particle = Mesh.DP[point].index;

          // Map Delauney ghost cell to physical cell
          if(particle >= NumGas && Mesh.DP[point].task == ThisTask)
            particle -= NumGas;

          double vx = 0.0, vy = 0.0, vz = 0.0;
          if(Mesh.DP[point].task == ThisTask)
            {
              // for reflective boundaries the value is constant across the interface, for outflow boundaries we
              // set it to zero. alternatively one could set it to a constant background value here
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
              terminate("Braginskii: Reflective boundaries not implemented!");
#endif
              vx = sphp_copy[particle].Momentum[0] / P[particle].Mass;
              vy = sphp_copy[particle].Momentum[1] / P[particle].Mass;
              vz = sphp_copy[particle].Momentum[2] / P[particle].Mass;
            }
          else
            {
              double mass = PrimExch[particle].Volume * PrimExch[particle].Density;
              vx          = sphp_copyExch[particle].Momentum[0] / mass;
              vy          = sphp_copyExch[particle].Momentum[1] / mass;
              vz          = sphp_copyExch[particle].Momentum[2] / mass;
            }

          cv->vx += vx * cd->matrix[p][0];
          cv->vy += vy * cd->matrix[p][0];
          cv->vz += vz * cd->matrix[p][0];

          for(k = 0; k < NUMDIMS; k++)
            {
              cv->vxGrad[k] += vx * cd->matrix[p][k + 1];
              cv->vyGrad[k] += vy * cd->matrix[p][k + 1];
              cv->vzGrad[k] += vz * cd->matrix[p][k + 1];
            }
        }

      for(k = NUMDIMS; k < 3; k++)
        {
          cv->vxGrad[k] = 0.0;
          cv->vyGrad[k] = 0.0;
          cv->vzGrad[k] = 0.0;
        }
    }
}

void set_diffusion_coefficients(struct diff_face_data *diff_face_data)
{
  // if(All.ComovingIntegrationOn)
  //   coeff *= All.cf_atime * All.cf_atime / All.HubbleParam;

  for(int iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];

      if(fd->active)
        {
#ifndef BRAGINSKII_SPITZER
          fd->eta_aniso = 0.0;
          for(int icorner = 0; icorner < fd->cornerCount; icorner++)
            {
              struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
              struct corner_data *cd = &corner_data[cl->index];
              fd->eta_aniso += cl->weight * cd->eta_aniso;
            }
#else
          // Harmonic mean of viscosity values at corners
          double sum = 0.0;
          for(int icorner = 0; icorner < fd->cornerCount; icorner++)
            {
              struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
              struct corner_data *cd = &corner_data[cl->index];
              sum += cl->weight / cd->eta_aniso;
            }
          if(sum > 0.0)
            fd->eta_aniso = 1.0 / sum;
          else
            fd->eta_aniso = 0.0;
#endif
        }
    }
}

/*
Calculates the velocity and the gradient of the velocity in the local
coordinate system at Voronoi cell interfaces by using the values in Cartesian
coordinate system at corners.
*/
void calculate_grad_v_in_local_coordinates_without_limiters(struct diff_face_data *fd, struct grad_v *grad_v)
{
  double nx = fd->nx;
  double mx = fd->mx;
  double px = fd->px;
  double ny = fd->ny;
  double my = fd->my;
  double py = fd->py;
  double nz = fd->nz;
  double mz = fd->mz;
  double pz = fd->pz;

  // Reset to zero
  grad_v->vn    = 0.0;
  grad_v->vm    = 0.0;
  grad_v->dvndn = 0.0;
  grad_v->dvmdn = 0.0;
  grad_v->dvndm = 0.0;
  grad_v->dvmdm = 0.0;
  grad_v->dvpdp = 0.0;

  // Loop over all corners
  for(int icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_vel *cv  = &corner_vel[cl->index];

      // Velocity in Cartesian system
      double vx = cv->vx;
      double vy = cv->vy;
      double vz = cv->vz;

      // Elements of the gradient tensor nabla v in Cartesian system
      double dvxdx = cv->vxGrad[0], dvxdy = cv->vxGrad[1], dvxdz = cv->vxGrad[2];
      double dvydx = cv->vyGrad[0], dvydy = cv->vyGrad[1], dvydz = cv->vyGrad[2];
      double dvzdx = cv->vzGrad[0], dvzdy = cv->vzGrad[1], dvzdz = cv->vzGrad[2];

      // Components of v in the local coordinate system
      grad_v->vn += cl->weight * (nx * vx + ny * vy + nz * vz);
      grad_v->vm += cl->weight * (mx * vx + my * vy + mz * vz);

      // Elements of the gradient tensor nabla v in local system
      grad_v->dvndn += cl->weight * (nx * (nx * dvxdx + ny * dvydx + nz * dvzdx) + ny * (nx * dvxdy + ny * dvydy + nz * dvzdy) +
                                     nz * (nx * dvxdz + ny * dvydz + nz * dvzdz));
      grad_v->dvmdn += cl->weight * (dvxdx * mx * nx + dvydx * my * nx + dvzdx * mz * nx + dvxdy * mx * ny + dvydy * my * ny +
                                     dvzdy * mz * ny + dvxdz * mx * nz + dvydz * my * nz + dvzdz * mz * nz);
      grad_v->dvndm += cl->weight * (dvxdx * mx * nx + dvxdy * my * nx + dvxdz * mz * nx + dvydx * mx * ny + dvydy * my * ny +
                                     dvydz * mz * ny + dvzdx * mx * nz + dvzdy * my * nz + dvzdz * mz * nz);
      grad_v->dvmdm += cl->weight * (dvxdx * mx * mx + dvxdy * mx * my + dvydx * mx * my + dvydy * my * my + dvxdz * mx * mz +
                                     dvzdx * mx * mz + dvydz * my * mz + dvzdy * my * mz + dvzdz * mz * mz);
      grad_v->dvpdp += cl->weight * (dvxdx * px * px + dvxdy * px * py + dvydx * px * py + dvydy * py * py + dvxdz * px * pz +
                                     dvzdx * px * pz + dvydz * py * pz + dvzdy * py * pz + dvzdz * pz * pz);
    }
}

/*
Calculates the velocity and the gradient of the velocity in the local
coordinate system at Voronoi cell interfaces by using the values in Cartesian
coordinate system at corners.
*/
void calculate_grad_v_in_local_coordinates(struct diff_face_data *fd, struct grad_v *grad_v)
{
  double nx = fd->nx;
  double mx = fd->mx;
  double ny = fd->ny;
  double my = fd->my;
  double nz = fd->nz;
  double mz = fd->mz;

  // Reset to zero
  grad_v->vn = 0.0;
  grad_v->vm = 0.0;

  // Loop over all corners
  for(int icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_vel *cv  = &corner_vel[cl->index];

      // Velocity in Cartesian system
      double vx = cv->vx;
      double vy = cv->vy;
      double vz = cv->vz;

      // Components of v in the local coordinate system
      grad_v->vn += cl->weight * (nx * vx + ny * vy + nz * vz);
      grad_v->vm += cl->weight * (mx * vx + my * vy + mz * vz);
    }

// Calculated limited derivatives
#ifndef BRAGINSKII_SIMPLE_DERIVATIVES
  grad_v->dvndn = limited_derivative(fd, 1);
  grad_v->dvmdn = limited_derivative(fd, 2);
#endif
  grad_v->dvndm = limited_derivative(fd, 3);
  grad_v->dvmdm = limited_derivative(fd, 4);
  grad_v->dvpdp = limited_derivative(fd, 5);
}

#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
void calculate_simple_normal_gradients(struct diff_face_data *fd, int iface, struct grad_v *grad_v)
{
  double Center1[3], Center2[3];
  point_get_center(Mesh.VF[iface].p1, Center1);
  point_get_center(Mesh.VF[iface].p2, Center2);

  double dx = NEAREST_X(Center2[0] - Center1[0]);
  double dy = NEAREST_Y(Center2[1] - Center1[1]);
  double dz = NEAREST_Z(Center2[2] - Center1[2]);

  double dist2 = dx * dx + dy * dy + dz * dz;

  if(dist2 > 0.)
    {
      double nx = fd->nx;
      double mx = fd->mx;
      double ny = fd->ny;
      double my = fd->my;
      double nz = fd->nz;
      double mz = fd->mz;

      double vx1 = center_vel[Mesh.VF[iface].p1].vx;
      double vx2 = center_vel[Mesh.VF[iface].p2].vx;

      double vy1 = center_vel[Mesh.VF[iface].p1].vy;
      double vy2 = center_vel[Mesh.VF[iface].p2].vy;

      double vz1 = center_vel[Mesh.VF[iface].p1].vz;
      double vz2 = center_vel[Mesh.VF[iface].p2].vz;

      double vn1 = nx * vx1 + ny * vy1 + nz * vz1;
      double vn2 = nx * vx2 + ny * vy2 + nz * vz2;

      double vm1 = mx * vx1 + my * vy1 + mz * vz1;
      double vm2 = mx * vx2 + my * vy2 + mz * vz2;

      // grad_v->dvndn = (vn2 - vn1)/sqrt(dist2);
      // grad_v->dvmdn = (vm2 - vm1)/sqrt(dist2);
      grad_v->dvndn = (vn2 - vn1) * (dx * fd->nx + dy * fd->ny + dz * fd->nz) / dist2;
      grad_v->dvmdn = (vm2 - vm1) * (dx * fd->nx + dy * fd->ny + dz * fd->nz) / dist2;
    }
  else
    {
      grad_v->dvndn = 0.0;
      grad_v->dvmdn = 0.0;
    }
}
#endif

double limited_derivative(struct diff_face_data *fd, const int mycase)
{
  double first = 0.0, der, der_limited, sum = 0.0;
  int signs_ok;
  int failCount = 0;

  double nx = fd->nx;
  double mx = fd->mx;
  double px = fd->px;
  double ny = fd->ny;
  double my = fd->my;
  double py = fd->py;
  double nz = fd->nz;
  double mz = fd->mz;
  double pz = fd->pz;

  for(int icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_vel *cv  = &corner_vel[cl->index];
      struct corner_data *cd = &corner_data[cl->index];

      if(cd->fail > 0)
        failCount++;

      // Elements of the gradient tensor nabla v in Cartesian system
      double dvxdx = cv->vxGrad[0], dvxdy = cv->vxGrad[1], dvxdz = cv->vxGrad[2];
      double dvydx = cv->vyGrad[0], dvydy = cv->vyGrad[1], dvydz = cv->vyGrad[2];
      double dvzdx = cv->vzGrad[0], dvzdy = cv->vzGrad[1], dvzdz = cv->vzGrad[2];

      switch(mycase)
        {
          case 1:
            // dvndn
            der = nx * (nx * dvxdx + ny * dvydx + nz * dvzdx) + ny * (nx * dvxdy + ny * dvydy + nz * dvzdy) +
                  nz * (nx * dvxdz + ny * dvydz + nz * dvzdz);
            break;
          case 2:
            // dvmdn
            der = dvxdx * mx * nx + dvydx * my * nx + dvzdx * mz * nx + dvxdy * mx * ny + dvydy * my * ny + dvzdy * mz * ny +
                  dvxdz * mx * nz + dvydz * my * nz + dvzdz * mz * nz;
            break;
          case 3:
            // dvndm
            der = dvxdx * mx * nx + dvxdy * my * nx + dvxdz * mz * nx + dvydx * mx * ny + dvydy * my * ny + dvydz * mz * ny +
                  dvzdx * mx * nz + dvzdy * my * nz + dvzdz * mz * nz;
            break;
          case 4:
            // dvmdm
            der = dvxdx * mx * mx + dvxdy * mx * my + dvydx * mx * my + dvydy * my * my + dvxdz * mx * mz + dvzdx * mx * mz +
                  dvydz * my * mz + dvzdy * my * mz + dvzdz * mz * mz;
            break;
          case 5:
            // dvpdp
            der = dvxdx * px * px + dvxdy * px * py + dvydx * px * py + dvydy * py * py + dvxdz * px * pz + dvzdx * px * pz +
                  dvydz * py * pz + dvzdy * py * pz + dvzdz * pz * pz;
            break;
          default:
            terminate("limited_derivative: mycase should be 1-5");
        }

      // Only execute this code if corner has not failed
      if(!cd->fail)
        {
          // Store the first non-zero estimate
          if(first == 0.0)
            {
              first    = der;
              signs_ok = 1;
            }
          else
            {
              // Check if derivative estimates differ in sign
              if(der * first < 0.0)
                {
                  signs_ok = 0;
                }
            }
          if(der != 0.0)
            {
              sum += cl->weight / der;
            }
        }
    }

  // if (signs_ok == 1 && sum != 0.0 && failCount < fd->cornerCount){
  if(signs_ok == 1 && sum != 0.0 && (fd->cornerCount - failCount >= 2))
    {
      der_limited = 1.0 / sum;
    }
  else
    {
      der_limited = 0.0;
      // TODO: Use a different estimate when this happens.
    }

  return der_limited;
}

double limited_derivative_iso(struct diff_face_data *fd, const int mycase)
{
  double first = 0.0, der, der_limited, sum = 0.0;
  int signs_ok;
  int failCount = 0;

  for(int icorner = 0; icorner < fd->cornerCount; icorner++)
    {
      struct corner_list *cl = &corner_list[fd->cornerFirst + icorner];
      struct corner_vel *cv  = &corner_vel[cl->index];
      struct corner_data *cd = &corner_data[cl->index];

      if(cd->fail > 0)
        failCount++;

      // Elements of the gradient tensor nabla v in Cartesian system
      switch(mycase)
        {
          case 1:
            // dvxdx
            der = cv->vxGrad[0];
            break;
          case 2:
            // dvxdy
            der = cv->vxGrad[1];
            break;
          case 3:
            // dvxdz
            der = cv->vxGrad[2];
            break;
          case 4:
            // dvydx
            der = cv->vyGrad[0];
            break;
          case 5:
            // dvydy
            der = cv->vyGrad[1];
            break;
          case 6:
            // dvydy
            der = cv->vyGrad[2];
            break;
          case 7:
            // dvydz
            der = cv->vzGrad[0];
            break;
          case 8:
            // dvzdy
            der = cv->vzGrad[1];
            break;
          case 9:
            // dvzdz
            der = cv->vzGrad[2];
            break;
          default:
            terminate("limited_derivative: mycase should be 1-9");
        }

      // Only execute this code if corner has not failed
      if(!cd->fail)
        {
          // Store the first non-zero estimate
          if(first == 0.0)
            {
              first    = der;
              signs_ok = 1;
            }
          else
            {
              // Check if derivative estimates differ in sign
              if(der * first < 0.0)
                {
                  signs_ok = 0;
                }
            }
          if(der != 0.0)
            {
              sum += cl->weight / der;
            }
        }
    }

  // if (signs_ok == 1 && sum != 0.0 && failCount < fd->cornerCount){
  if(signs_ok == 1 && sum != 0.0 && (fd->cornerCount - failCount >= 2))
    {
      der_limited = 1.0 / sum;
    }
  else
    {
      der_limited = 0.0;
      // TODO: Use a different estimate when this happens.
    }

  return der_limited;
}

/* Copy state from first input to second input */
void copy_state(struct sphp_copy *state1, struct sphp_copy *state2)
{
  for(int part = 0; part < NumGas; part++)
    {
      state2[part].Energy      = state1[part].Energy;
      state2[part].Momentum[0] = state1[part].Momentum[0];
      state2[part].Momentum[1] = state1[part].Momentum[1];
      state2[part].Momentum[2] = state1[part].Momentum[2];
    }
}

// Copy some of the values in SphP to sphp_copy
void copy_SphP(struct sphp_copy *sphp_copy)
{
  for(int part = 0; part < NumGas; part++)
    {
      sphp_copy[part].Energy      = SphP[part].Energy;
      sphp_copy[part].Momentum[0] = SphP[part].Momentum[0];
      sphp_copy[part].Momentum[1] = SphP[part].Momentum[1];
      sphp_copy[part].Momentum[2] = SphP[part].Momentum[2];
    }
}

/* Copy some of the values from sphp_copy to SphP. */
void set_SphP(struct sphp_copy *sphp_copy)
{
  for(int part = 0; part < NumGas; part++)
    {
      SphP[part].Energy      = sphp_copy[part].Energy;
      SphP[part].Momentum[0] = sphp_copy[part].Momentum[0];
      SphP[part].Momentum[1] = sphp_copy[part].Momentum[1];
      SphP[part].Momentum[2] = sphp_copy[part].Momentum[2];
    }
}

/*
Apply the calculated velocity flux directly on this task and calculate a flux
list for later communication to other tasks. The latter is used by
exchange_and_apply_fluxes.
*/
int apply_flux_or_calculate_fluxlist(struct viscflux_list_data *ViscFluxList, struct sphp_copy *sphp_copy, struct diff_face_data *fd,
                                     double flux_E, double flux_Mx, double flux_My, double flux_Mz, int iface, int *MaxNflux,
                                     int Nflux)
{
  int k;
  for(k = 0; k < 2; k++)
    {
      double dir = +1. - k * 2.; /* switches between -1 and +1 */
      int q      = (k == 0 ? Mesh.VF[iface].p1 : Mesh.VF[iface].p2);
      int p      = Mesh.DP[q].index;

      if(Mesh.DP[q].task == ThisTask)
        {
          if(Mesh.DP[q].index >= NumGas) /* this is a local ghost point */
            {
              if(Mesh.DP[Mesh.VF[iface].p1].ID == Mesh.DP[Mesh.VF[iface].p2].ID) /* this may happen for reflective points */
                continue;
              p -= NumGas;
            }

          sphp_copy[p].Energy += dir * fd->area * flux_E;
          sphp_copy[p].Momentum[0] += dir * fd->area * flux_Mx;
          sphp_copy[p].Momentum[1] += dir * fd->area * flux_My;
          sphp_copy[p].Momentum[2] += dir * fd->area * flux_Mz;
        }
      else
        {
          if(Nflux >= *MaxNflux)
            {
              Mesh.Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
              *MaxNflux = Mesh.Indi.AllocFacNflux;
              ViscFluxList =
                  (struct viscflux_list_data *)myrealloc_movable(ViscFluxList, *MaxNflux * sizeof(struct viscflux_list_data));

              if(Nflux >= *MaxNflux)
                terminate("Nflux >= MaxNflux");
            }

          ViscFluxList[Nflux].task    = Mesh.DP[q].task;
          ViscFluxList[Nflux].index   = Mesh.DP[q].originalindex;
          ViscFluxList[Nflux].dEnergy = dir * fd->area * flux_E;
          ViscFluxList[Nflux].dMx     = dir * fd->area * flux_Mx;
          ViscFluxList[Nflux].dMy     = dir * fd->area * flux_My;
          ViscFluxList[Nflux].dMz     = dir * fd->area * flux_Mz;

          Nflux++;
        }
    }
  return Nflux;
}

/*
Communicate velocity fluxes calculated in apply_flux_or_calculate_fluxlist
to other tasks.
*/
void exchange_and_apply_fluxes(struct viscflux_list_data *ViscFluxList, int Nflux, struct sphp_copy *sphp_copy)
{
  /* exchange and apply fluxes */
  int i, j, nimport, ngrp, recvTask;
  mysort(ViscFluxList, Nflux, sizeof(struct viscflux_list_data), viscflux_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[ViscFluxList[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate("Send_count[ThisTask]=%d", Send_count[ThisTask]);

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct viscflux_list_data *ViscFluxListGet =
      (struct viscflux_list_data *)mymalloc("FluxListGet", nimport * sizeof(struct viscflux_list_data));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&ViscFluxList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct viscflux_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &ViscFluxListGet[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct viscflux_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < nimport; i++)
    {
      int p = ViscFluxListGet[i].index;
      sphp_copy[p].Energy += ViscFluxListGet[i].dEnergy;
      sphp_copy[p].Momentum[0] += ViscFluxListGet[i].dMx;
      sphp_copy[p].Momentum[1] += ViscFluxListGet[i].dMy;
      sphp_copy[p].Momentum[2] += ViscFluxListGet[i].dMz;
    }

  myfree(ViscFluxListGet);
}

#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
/*
Calculate RKL2 weights.
*/
void find_rkl_weights(double *mu, double *mu_til, double *nu, double *gam_til, int s)
{
  int j;
  double b[s + 1];
  double a[s + 1];

  /*Initialize to zero*/
  for(j = 0; j < s + 1; j++)
    {
      mu[j]      = 0.0;
      mu_til[j]  = 0.0;
      nu[j]      = 0.0;
      gam_til[j] = 0.0;
      b[j]       = 0.0;
    }

  double w1 = 4.0 / (s * s + s - 2.0);

  for(j = 0; j < 3; j++)
    {
      b[j] = 1.0 / 3.0;
    }

  for(j = 2; j < s + 1; j++)
    {
      b[j] = (j * j + j - 2.0) / (2.0 * j * (j + 1.0));
    }

  for(j = 0; j < s + 1; j++)
    {
      a[j] = 1.0 - b[j];
    }

  mu_til[1] = b[1] * w1;

  for(j = 2; j < s + 1; j++)
    {
      mu[j]      = (2.0 * j - 1.0) / j * b[j] / b[j - 1];
      mu_til[j]  = mu[j] * w1;
      nu[j]      = -(j - 1.0) / j * b[j] / b[j - 2];
      gam_til[j] = -a[j - 1] * mu_til[j];
    }
}
#endif  // STS

void point_get_center(int p, double *Center)
{
  if(!TimeBinSynchronized[Mesh.DP[p].timebin])
    {
      Center[0] = Mesh.DP[p].x;
      Center[1] = Mesh.DP[p].y;
      Center[2] = Mesh.DP[p].z;
      return;
    }

  int particle = Mesh.DP[p].index;
  int j;
  if(Mesh.DP[p].task == ThisTask)
    {
      if(particle >= NumGas)
        particle -= NumGas;

      for(j = 0; j < 3; j++)
        Center[j] = SphP[particle].Center[j];

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      if(pIsReflective[p])
        for(j = 0; j < 3; j++)
          Center[j] = refCenter[p * 3 + j];
#endif
    }
  else
    {
      for(j = 0; j < 3; j++)
        Center[j] = PrimExch[particle].Center[j];
    }
}

int viscflux_list_data_compare(const void *a, const void *b)
{
  if(((struct viscflux_list_data *)a)->task < (((struct viscflux_list_data *)b)->task))
    return -1;

  if(((struct viscflux_list_data *)a)->task > (((struct viscflux_list_data *)b)->task))
    return +1;

  return 0;
}

/*
This function is used to communicate the values in a state between tasks.
Values not in the present task can be accessed in stateExch after this
function has been called. This function is simply a slightly modified version
of exchange_primitive_variables which can be found in voronoi_exchange.c
*/
void exchange_state_variables(struct sphp_copy *state, struct sphp_copy *stateExch)
{
  int listp;
  struct sphp_copy *tmpStateExch;
  int i, j, p, task, off;
  int ngrp, recvTask, place;

  tmpStateExch = (struct sphp_copy *)mymalloc("tmpStateExch", Mesh_nexport * sizeof(struct sphp_copy));

  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpStateExch[off].Energy      = state[place].Energy;
              tmpStateExch[off].Momentum[0] = state[place].Momentum[0];
              tmpStateExch[off].Momentum[1] = state[place].Momentum[1];
              tmpStateExch[off].Momentum[2] = state[place].Momentum[2];
            }
          listp = ListExports[listp].nextexport;
        }
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&tmpStateExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct sphp_copy), MPI_BYTE,
                           recvTask, TAG_DENS_A, &stateExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct sphp_copy), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpStateExch);
}

void free_stuff(void)
{
  /* free temporary arrays */
  myfree(corner_vel);
  myfree(corner_data);
  myfree(corner_list);
  myfree(diff_face_data);
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  myfree(pIsReflective);
  myfree(refCenter);
#endif
}

#endif
