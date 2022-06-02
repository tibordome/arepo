/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_thermochem.c
 * \date        06/2018
 * \author      Rahul Kannan
 * \brief       Coupled thermochemistry solver - solves for chemistry and cooling of H,HE
 * \details     First uses the PS2009 semi-implicit scheme to solve for ionization rates
 *              and cooling. If the change is larger than 10%, switches to solving the
 *              equations using the CVODE solver from the SUNDIALS library.
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

/* -----------------------------------------------------------------
 * Based on CVODE solver from SUNDIALS library
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------*/

#include <stdio.h>

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"
#include "RT.h"

#ifdef MRT_COUPLED_THERMOCHEMISTRY

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   SM_ELEMENT_D macro in dense.h. SM_ELEMENT_D numbers rows and columns of
   a dense matrix starting from 0. */

#define Ith(v, i) NV_Ith_S(v, i - 1)                /* Ith numbers components 1..NEQ */
#define IJth(A, i, j) SM_ELEMENT_D(A, i - 1, j - 1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */

#define NEQ 4 /* number of equations  */

#define RTOL RCONST(1.0e-6)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-10) /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-10)
#define ATOL3 RCONST(1.0e-10)
#define ATOL4 RCONST(1.0e-10)
#define ATOL5 RCONST(1.0e-10)

#define ZERO RCONST(0.0)
#define NCOEFFS 34

static const int ALPHA_HII         = 0;
static const int SIGMA_HI          = 1;
static const int ALPHA_HEIII       = 2;
static const int SIGMA_HEI         = 3;
static const int ALPHA_HEII        = 4;
static const int SIGMA_HEII        = 5;
static const int COOL_I            = 6;
static const int COOL_II           = 7;
static const int COOL_III          = 8;
static const int COOL_IV           = 9;
static const int COOL_V            = 10;
static const int COOL_VI           = 11;
static const int COOL_VII          = 12;
static const int COOL_VIII         = 13;
static const int COOL_IX           = 14;
static const int COOL_X            = 15;
static const int COOL_XI           = 16;
static const int PHOT_ION_HI       = 17;
static const int PHOT_ION_HEI      = 18;
static const int PHOT_ION_HEII     = 19;
static const int PHOT_COOL_HI      = 20;
static const int PHOT_COOL_HEI     = 21;
static const int PHOT_COOL_HEII    = 22;
static const int IR_PHOT_COEFF_I   = 23;
static const int IR_PHOT_COEFF_II  = 24;
static const int IR_TEMP_COEFF_I   = 25;
static const int IR_TEMP_COEFF_II  = 26;
static const int MOLECULAR_COOLING = 27;
static const int PE_HEATING_I      = 28;
static const int PE_HEATING_II     = 29;
static const int METAL_COOLING     = 30;
/*Inverse Compton cooling of the CMB*/
static const int COOL_INV_COMPTON = 31;
/*Dielectric recombination and cooling for HeII*/
static const int ALPHA_HEII_DIELEC = 32;
static const int COOL_HEII_DIELEC  = 33;

#ifdef MRT_PHOTOELECTRIC_HEATING
static const int FOPT        = 0;    // frequency bin number for photoelctric heating
static const double FRAC_FUV = 0.6;  // photon in the FUV range in optical bin (13.6-5.6)/13.6 ~ 0.6
static double fac_pe;
static double e_pe;
#endif

static double cspeed;
static double fac_chem;
static double fac_cool;
static double coeff[NCOEFFS];

static double clight;

#ifdef MRT_UVB
static double UVBionHI;
static double UVBionHeI;
static double UVBionHeII;
static double UVBcoolHI;
static double UVBcoolHeI;
static double UVBcoolHeII;

static const double ss_alpha1 = -3.98;
static const double ss_alpha2 = -1.09;
static const double ss_beta   = 1.29;
static const double ss_xi     = 1;
static const double ss_n0     = -3.11;
static const double ss_f      = 0.01;
#endif

// static double utherm_particle ;

#ifdef MRT_MOLECULAR_COOLING
static double fratio;
#endif

#ifdef MRT_METAL_COOLING
static double lognH;
#endif

static long int MAXSTEPS = 100000;

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);

/*routine that calls the main function*/

static int calculate_thermochem(int i, double dt);
static void set_coeffs(int i, double nH, double rho, double T, double inv_cv);
static void set_coeffs_photons(int i, double nH, double rho, double T, double inv_cv);
static void set_out(int i, N_Vector y, double nH, double rho, double T, double dt);
static double functional_form(double T, double k1, double k2, double k3, double k4, double k5, double k6, double l, double m,
                              double n);
static double differential_form(double T, double k1, double k2, double k3, double k4, double k5, double k6, double l, double m,
                                double n);
static double functional(int num, double T);
static double differential(int num, double T);
static int mrt_thermochem_semi_implicit(int i, double dt, double T);
static void absorb_photons(int i, double dt, double nH, double T);

static double compton_molecular(double T);
static double compton_differential(double T);

static void set_therm_photons(int i);
static void set_updated_photons(int i);

static double dielectric_functional(double T);
static double dielectric_differential(double T);

#ifdef MRT_MOLECULAR_COOLING
static double molecular_functional(double T);
static double molecular_differential(double T);
#endif

#ifdef MRT_PHOTOELECTRIC_HEATING
static double pe_functional(double T, int mode);
static double pe_differential(double T, int mode);
#endif

#ifdef MRT_METAL_COOLING
static double metal_cooling(double T, int mode);
#endif

static struct thermstruct
{
  double DensPhot[MRT_BINS];
  double Cons_DensPhot[MRT_BINS];

  double rho, nH, inv_cv, temp, initU, mu;

} therm;

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

void mrt_thermochemistry(void)
{
  TIMER_START(CPU_RT_CHEM);
  int idx, i;
  double dt;

  mpi_printf("RT: Updating coupled thermochemistry and cooling with CVODE.....\n");

  fac_chem = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) * All.HubbleParam * All.HubbleParam;
  fac_cool = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs * pow(All.HubbleParam, 3);

  cspeed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s;

#ifdef MRT_DUAL_LIGHTSPEED
  clight = cspeed;
#else
  clight = c_internal_units;
#endif

#ifdef MRT_UVB
  UVBionHI    = 3.76244e-14;
  UVBionHeI   = 2.0821e-14;
  UVBionHeII  = 1.12165e-16;
  UVBcoolHI   = 2.4774e-25;
  UVBcoolHeI  = 2.21352e-25;
  UVBcoolHeII = 5.0084e-27;

  UVBionHI *= All.UnitTime_in_s;
  UVBionHeI *= All.UnitTime_in_s;
  UVBionHeII *= All.UnitTime_in_s;

  UVBcoolHI *= All.UnitTime_in_s / All.UnitEnergy_in_cgs;
  UVBcoolHeI *= All.UnitTime_in_s / All.UnitEnergy_in_cgs;
  UVBcoolHeII *= All.UnitTime_in_s / All.UnitEnergy_in_cgs;
#endif

#ifdef MRT_PHOTOELECTRIC_HEATING
  fac_pe = 3.9e-14 * pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs * All.HubbleParam * All.HubbleParam;
#endif

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn)
        dt /= All.cf_hubble_a;

      if(All.Time == 0.0)
        {
          TIMER_STOP(CPU_RT_CHEM);
          return;
        }

#ifdef MRT_SUBCYCLE
      dt /= ((double)(All.RTNumSubCycles));
#else
      dt *= 0.5;
#endif

#if defined(COOLING) && !defined(GRACKLE)
      if(do_cooling_mrt(i))
#endif

        set_therm_photons(i);
#ifndef MRT_DUAL_LIGHTSPEED
      calculate_thermochem(i, dt);
#else
      int num     = ((int)(clight / c_internal_units));
      double frac = (1.0 / ((double)(num)));
      for(int kk = 0; kk < num; kk++)
        calculate_thermochem(i, dt * frac);
#endif
      set_updated_photons(i);
    }
  TIMER_STOP(CPU_RT_CHEM);
  return;
}

static void set_therm_photons(int i)
{
  therm.rho               = (P[i].Mass / SphP[i].Volume) * All.cf_a3inv;
  therm.nH                = HYDROGEN_MASSFRAC * therm.rho / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
  double molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
  //  molecular_weight = 1.0 ;
#ifdef MRT_LEVITATION_TEST
  molecular_weight = 2.33;
#endif
  therm.mu = molecular_weight;

  therm.inv_cv = GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) /
                 (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

  double vel2 = SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] +
                SphP[i].Momentum[2] * SphP[i].Momentum[2];
  double utherm_particle = (SphP[i].Energy - 0.5 * vel2 / P[i].Mass) / P[i].Mass / (All.cf_atime * All.cf_atime);

  therm.initU = utherm_particle;
  therm.temp  = utherm_particle * therm.inv_cv;

  for(int kk = 0; kk < MRT_BINS; kk++)
    therm.Cons_DensPhot[kk] = SphP[i].DensPhot[kk] * SphP[i].Volume;

  set_coeffs(i, therm.nH, therm.rho, therm.temp, therm.inv_cv);
}

static void set_updated_photons(int i)
{
  double du = therm.temp / therm.inv_cv - therm.initU;

  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

  for(int kk = 0; kk < MRT_BINS; kk++)
    {
      double ratio = therm.Cons_DensPhot[kk] / SphP[i].Volume / SphP[i].DensPhot[kk];
      SphP[i].Cons_DensPhot[kk] *= ratio;
      for(int ll = 0; ll < 3; ll++)
        SphP[i].Cons_RT_F[kk][ll] *= ratio;
    }
}

int calculate_thermochem(int i, double dt)
{
  for(int kk = 0; kk < MRT_BINS; kk++)
    therm.DensPhot[kk] = therm.Cons_DensPhot[kk] / SphP[i].Volume;

  double rho              = therm.rho;
  double nH               = therm.nH;
  double temp             = therm.temp;
  double molecular_weight = therm.mu;
  double inv_cv           = therm.inv_cv;

  if(temp == All.MinGasTemp)
    temp = All.MinGasTemp + 5.0;

  double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

#ifdef MRT_MOLECULAR_COOLING
  fratio = 1.0;
#endif

#ifndef MRT_NO_UV
  absorb_photons(i, dt, nH, temp);
#endif

  set_coeffs_photons(i, nH, rho, temp, inv_cv);

  double nh2 = SphP[i].HII;
#ifdef MRT_INCLUDE_HE
  double nhe2 = SphP[i].HeII;
  double nhe3 = SphP[i].HeIII;
#endif

#ifdef MRT_METAL_COOLING
  lognH = log10(nH * pow(All.UnitLength_in_cm, -3));
#endif

  int cont_flag = mrt_thermochem_semi_implicit(i, dt, temp);
  if(!cont_flag)
    return 1;

  realtype reltol, t, tout;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int flag;
  booleantype check_negative;

  y = abstol = NULL;
  A          = NULL;
  LS         = NULL;
  cvode_mem  = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  if(check_flag((void *)y, "N_VNew_Serial", 0))
    return (1);
  abstol = N_VNew_Serial(NEQ);
  if(check_flag((void *)abstol, "N_VNew_Serial", 0))
    return (1);

  /* Initialize y */
  Ith(y, 1) = SphP[i].HII;
  Ith(y, 2) = SphP[i].HeII;
  Ith(y, 3) = SphP[i].HeIII;
  Ith(y, 4) = temp;

  if(SphP[i].HII < 0.0 || SphP[i].HeII < 0.0 || SphP[i].HeIII < 0.0 || SphP[i].HII > 1.0 || SphP[i].HeII > y_fac ||
     SphP[i].HeIII > y_fac || isnan(temp))
    terminate("%g %g %g %g %g %g %g\n", SphP[i].HII, SphP[i].HeII, SphP[i].HeIII, temp, SphP[i].Ne, temp, SphP[i].Energy);

  // Ith(y,4) = temp ;

  /* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
  Ith(abstol, 1) = ATOL1;
  Ith(abstol, 2) = ATOL2;
  Ith(abstol, 3) = ATOL3;
  Ith(abstol, 4) = ATOL4;

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF);  //, CV_NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0))
    return (1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, 0.0, y);
  if(check_flag(&flag, "CVodeInit", 1))
    return (1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if(check_flag(&flag, "CVodeSVtolerances", 1))
    return (1);

  flag = CVodeSetUserData(cvode_mem, &check_negative);
  flag = CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);

  /* Call CVodeRootInit to specify the root function g with 2 components */
  flag = CVodeRootInit(cvode_mem, 1, g);
  if(check_flag(&flag, "CVodeRootInit", 1))
    return (1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_flag((void *)A, "SUNDenseMatrix", 0))
    return (1);

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNDenseLinearSolver(y, A);
  if(check_flag((void *)LS, "SUNDenseLinearSolver", 0))
    return (1);

  /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if(check_flag(&flag, "CVDlsSetLinearSolver", 1))
    return (1);

  /* Set the user-supplied Jacobian routine Jac */
  flag = CVDlsSetJacFn(cvode_mem, Jac);
  if(check_flag(&flag, "CVDlsSetJacFn", 1))
    return (1);

  check_negative = 0;

  t    = 0.0;
  tout = dt;

  flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
  if(flag != CV_SUCCESS && flag != CV_ROOT_RETURN)
    {
      print_particle_info(i);
      terminate("CVODE ERROR FLAG %d ID %d|| \n Density = %g || \n Old = %g %g %g %g %g|| \n New = %g %g %g %g %g\n", flag, P[i].ID,
                nH * pow(All.UnitLength_in_cm, -3), SphP[i].HII, SphP[i].HeII, SphP[i].HeIII, (y_fac - SphP[i].HeII - SphP[i].HeIII),
                temp, Ith(y, 1), Ith(y, 2), Ith(y, 3), (y_fac - Ith(y, 2)), Ith(y, 3), Ith(y, 4));
    }

  set_out(i, y, nH, rho, temp, dt);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(A);

  return (flag);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * g routine. Compute function f(t,y).
 */

static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  gout[0] = Ith(y, 4) - All.MinGasTemp;
  return (0);
}

/*
 * f routine. Compute function f(t,y).
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3, y4, y5;

  booleantype *check_negative;

  check_negative = (booleantype *)user_data;
  double y_fac   = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

  y1 = Ith(y, 1);
  y2 = Ith(y, 2);
  y3 = Ith(y, 3);
  y4 = Ith(y, 4);

  if(*check_negative && (y4 < 0))
    return (1);

  double val[NCOEFFS];
  for(int i = 0; i < NCOEFFS; i++)
    val[i] = coeff[i] * functional(i, y4);

  val[ALPHA_HEII] += val[ALPHA_HEII_DIELEC];
  val[COOL_V] += val[COOL_HEII_DIELEC];

  double nh2 = y1;
  double nh1 = 1.0 - nh2;

  double nhe2 = y2;
  double nhe3 = y3;
  double nhe1 = y_fac - nhe2 - nhe3;

  double ne = nh2 + nhe2 + 2.0 * nhe3;

  Ith(ydot, 1) = -val[ALPHA_HII] * nh2 * ne + val[SIGMA_HI] * nh1 * ne + val[PHOT_ION_HI] * nh1;
  Ith(ydot, 2) = val[ALPHA_HEIII] * nhe3 * ne + val[SIGMA_HEI] * nhe1 * ne + val[PHOT_ION_HEI] * nhe1 - val[ALPHA_HEII] * nhe2 * ne -
                 val[SIGMA_HEII] * nhe2 * ne - val[PHOT_ION_HEII] * nhe2;
  Ith(ydot, 3) = -val[ALPHA_HEIII] * nhe3 * ne + val[SIGMA_HEII] * nhe2 * ne + val[PHOT_ION_HEII] * nhe2;

  Ith(ydot, 4) = val[PHOT_COOL_HI] * nh1 + val[PHOT_COOL_HEI] * nhe1 + val[PHOT_COOL_HEII] * nhe2 - val[COOL_I] * nh2 * ne -
                 (val[COOL_II] + val[COOL_III]) * nh1 * ne - (val[COOL_V] + val[COOL_VIII] + val[COOL_X]) * nhe2 * ne -
                 val[COOL_VI] * nhe3 * ne - val[COOL_VII] * nhe1 * ne - (val[COOL_IX] + val[COOL_XI]) * nhe2 * ne * ne -
                 val[COOL_IV] * (nh2 + nhe2 + 4.0 * nhe3) * ne - val[COOL_INV_COMPTON] * ne
#ifdef MRT_MOLECULAR_COOLING
                 - val[MOLECULAR_COOLING]
#endif
#ifdef MRT_PHOTOELECTRIC_HEATING
                 + val[PE_HEATING_I] + val[PE_HEATING_II]
#endif
#ifdef MRT_METAL_COOLING
                 - val[METAL_COOLING]
#endif
      ;

  return (0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

  realtype y1, y2, y3, y4, y5;

  y1 = Ith(y, 1);
  y2 = Ith(y, 2);
  y3 = Ith(y, 3);
  y4 = Ith(y, 4), y5 = Ith(y, 5);

  double val[NCOEFFS];
  double diff[NCOEFFS];

  for(int i = 0; i < NCOEFFS; i++)
    {
      val[i]  = coeff[i] * functional(i, y4);
      diff[i] = val[i] * differential(i, y4);
    }

  val[ALPHA_HEII] += val[ALPHA_HEII_DIELEC];
  val[COOL_V] += val[COOL_HEII_DIELEC];

  diff[ALPHA_HEII] += diff[ALPHA_HEII_DIELEC];
  diff[COOL_V] += diff[COOL_HEII_DIELEC];

  double nh2 = y1;
  double nh1 = 1.0 - nh2;

  double nhe2 = y2;
  double nhe3 = y3;
  double nhe1 = y_fac - nhe2 - nhe3;

  double ne = nh2 + nhe2 + 2.0 * nhe3;

  double dnh1_dnh2   = -1.0;
  double dne_dnh2    = 1.0;
  double dnhe1_dnhe2 = -1.0;
  double dnhe1_dnhe3 = -1.0;
  double dne_dnhe2   = 1.0;
  double dne_dnhe3   = 2.0;

  IJth(J, 1, 1) = -val[ALPHA_HII] * ne - val[ALPHA_HII] * nh2 * dne_dnh2 + val[SIGMA_HI] * dnh1_dnh2 * ne +
                  val[SIGMA_HI] * nh1 * dne_dnh2 + val[PHOT_ION_HI] * dnh1_dnh2;
  IJth(J, 1, 2) = -val[ALPHA_HII] * nh2 * dne_dnhe2 + val[SIGMA_HI] * nh1 * dne_dnhe2;
  IJth(J, 1, 3) = -val[ALPHA_HII] * nh2 * dne_dnhe3 + val[SIGMA_HI] * nh1 * dne_dnhe3;
  IJth(J, 1, 4) = -diff[ALPHA_HII] * nh2 * ne + diff[SIGMA_HI] * nh1 * ne;

  IJth(J, 2, 1) = val[ALPHA_HEIII] * nhe3 * dne_dnh2 + val[SIGMA_HEI] * nhe1 * dne_dnh2 - val[ALPHA_HEII] * nhe2 * dne_dnh2 -
                  val[SIGMA_HEII] * nhe2 * dne_dnh2;
  IJth(J, 2, 2) = val[ALPHA_HEIII] * nhe3 * dne_dnhe2 + val[SIGMA_HEI] * dnhe1_dnhe2 * ne + val[SIGMA_HEI] * nhe1 * dne_dnhe2 +
                  val[PHOT_ION_HEI] * dnhe1_dnhe2 - val[ALPHA_HEII] * ne - val[ALPHA_HEII] * nhe2 * dne_dnhe2 - val[SIGMA_HEII] * ne -
                  val[SIGMA_HEII] * nhe2 * dne_dnhe2 - val[PHOT_ION_HEII];
  IJth(J, 2, 3) = val[ALPHA_HEIII] * ne + val[ALPHA_HEIII] * nhe3 * dne_dnhe3 + val[SIGMA_HEI] * dnhe1_dnhe3 * ne +
                  val[SIGMA_HEI] * nhe1 * dne_dnhe3 + val[PHOT_ION_HEI] * dnhe1_dnhe3 - val[ALPHA_HEII] * nhe2 * dne_dnhe3 -
                  val[SIGMA_HEII] * nhe2 * dne_dnhe3;
  IJth(J, 2, 4) =
      diff[ALPHA_HEIII] * nhe3 * ne + diff[SIGMA_HEI] * nhe1 * ne - diff[ALPHA_HEII] * nhe2 * ne - diff[SIGMA_HEII] * nhe2 * ne;

  IJth(J, 3, 1) = -val[ALPHA_HEIII] * nhe3 * dne_dnh2 + val[SIGMA_HEII] * nhe2 * dne_dnh2;
  IJth(J, 3, 2) =
      -val[ALPHA_HEIII] * nhe3 * dne_dnhe2 + val[SIGMA_HEII] * ne + val[SIGMA_HEII] * nhe2 * dne_dnhe2 + val[PHOT_ION_HEII];
  IJth(J, 3, 3) = -val[ALPHA_HEIII] * ne - val[ALPHA_HEIII] * nhe3 * dne_dnhe3 + val[SIGMA_HEII] * nhe2 * dne_dnhe3;
  IJth(J, 3, 4) = -diff[ALPHA_HEIII] * nhe3 * ne + diff[SIGMA_HEII] * nhe2 * ne;

  IJth(J, 4, 1) = val[PHOT_COOL_HI] * dnh1_dnh2 - val[COOL_I] * ne - val[COOL_I] * nh2 * dne_dnh2 -
                  (val[COOL_II] + val[COOL_III]) * dnh1_dnh2 * ne - (val[COOL_II] + val[COOL_III]) * nh1 * dne_dnh2 -
                  (val[COOL_V] + val[COOL_VIII] + val[COOL_X]) * nhe2 * dne_dnh2 - val[COOL_VI] * nhe3 * dne_dnh2 -
                  val[COOL_VII] * nhe1 * dne_dnh2 - (val[COOL_IX] + val[COOL_XI]) * nhe2 * 2.0 * ne * dne_dnh2 - val[COOL_IV] * ne -
                  val[COOL_IV] * (nh2 + nhe2 + 4.0 * nhe3) * dne_dnh2 - val[COOL_INV_COMPTON] * dne_dnh2;

  IJth(J, 4, 2) = val[PHOT_COOL_HEI] * dnhe1_dnhe2 + val[PHOT_COOL_HEII] - val[COOL_I] * nh2 * dne_dnhe2 -
                  (val[COOL_II] + val[COOL_III]) * nh1 * dne_dnhe2 - (val[COOL_V] + val[COOL_VIII] + val[COOL_X]) * ne -
                  (val[COOL_V] + val[COOL_VIII] + val[COOL_X]) * nhe2 * dne_dnhe2 - val[COOL_VI] * nhe3 * dne_dnhe2 -
                  val[COOL_VII] * dnhe1_dnhe2 * ne - val[COOL_VII] * nhe1 * dne_dnhe2 - (val[COOL_IX] + val[COOL_XI]) * ne * ne -
                  (val[COOL_IX] + val[COOL_XI]) * nhe2 * 2.0 * ne * dne_dnhe2 - val[COOL_IV] * ne -
                  val[COOL_IV] * (nh2 + nhe2 + 4.0 * nhe3) * dne_dnhe2 - val[COOL_INV_COMPTON] * dne_dnhe2;

  IJth(J, 4, 3) = val[PHOT_COOL_HEI] * dnhe1_dnhe3 - val[COOL_I] * nh2 * dne_dnhe3 - (val[COOL_II] + val[COOL_III]) * nh1 * dne_dnhe3 -
                  (val[COOL_V] + val[COOL_VIII] + val[COOL_X]) * nhe2 * dne_dnhe3 - val[COOL_VI] * ne -
                  val[COOL_VI] * nhe3 * dne_dnhe3 - val[COOL_VII] * dnhe1_dnhe3 * ne - val[COOL_VII] * nhe1 * dne_dnhe3 -
                  (val[COOL_IX] + val[COOL_XI]) * nhe2 * 2.0 * ne * dne_dnhe3 - val[COOL_IV] * 4.0 * ne -
                  val[COOL_IV] * (nh2 + nhe2 + 4.0 * nhe3) * dne_dnhe3 - val[COOL_INV_COMPTON] * dne_dnhe3;

  IJth(J, 4, 4) = -diff[COOL_I] * nh2 * ne - (diff[COOL_II] + diff[COOL_III]) * nh1 * ne -
                  (diff[COOL_V] + diff[COOL_VIII] + diff[COOL_X]) * nhe2 * ne - diff[COOL_VI] * nhe3 * ne -
                  diff[COOL_VII] * nhe1 * ne - (diff[COOL_IX] + diff[COOL_XI]) * nhe2 * ne * ne -
                  diff[COOL_IV] * (nh2 + nhe2 + 4.0 * nhe3) * ne - diff[COOL_INV_COMPTON] * ne
#ifdef MRT_MOLECULAR_COOLING
                  - diff[MOLECULAR_COOLING]
#endif
#ifdef MRT_PHOTOELECTRIC_HEATING
                  + diff[PE_HEATING_I] + diff[PE_HEATING_II]
#endif
#ifdef MRT_METAL_COOLING
                  - diff[METAL_COOLING]
#endif
      ;

  return (0);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if(opt == 0 && flagvalue == NULL)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
      return (1);
    }

  /* Check if flag < 0 */
  else if(opt == 1)
    {
      errflag = (int *)flagvalue;
      if(*errflag < 0)
        {
          fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
          return (1);
        }
    }

  /* Check if function returned NULL pointer - no memory allocated */
  else if(opt == 2 && flagvalue == NULL)
    {
      fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
      return (1);
    }

  return (0);
}

static void set_coeffs(int i, double nH, double rho, double T, double inv_cv)
{
  for(int kk = 0; kk < NCOEFFS; kk++)
    coeff[kk] = 0.0;

#ifdef MRT_PHOTOELECTRIC_HEATING
  e_pe = 0.0;
#endif

#ifndef MRT_NO_UV
  coeff[ALPHA_HII] = nH * (2.59e-13 * fac_chem);  // alpha_HII
  coeff[SIGMA_HI]  = nH * (5.85e-11 * fac_chem);  // sigma_HI

#ifdef MRT_INCLUDE_HE
  coeff[ALPHA_HEIII]       = nH * (3.36e-10 * fac_chem);  // alpha_HeIII
  coeff[SIGMA_HEI]         = nH * (2.38e-11 * fac_chem);  // sigma_HeI
  coeff[ALPHA_HEII]        = nH * (1.5e-10 * fac_chem);   // alpha_HeII
  coeff[SIGMA_HEII]        = nH * (5.68e-12 * fac_chem);  // sigma_HeII
  coeff[ALPHA_HEII_DIELEC] = nH * (1.9e-3 * fac_chem);    // Dielec alpha_HeII
#endif

#ifdef MRT_COOLING_HEATING

  double val = nH * nH * inv_cv / rho;

  coeff[COOL_I]   = val * 8.7e-27 * fac_cool;
  coeff[COOL_II]  = val * 1.27e-21 * fac_cool;
  coeff[COOL_III] = val * 7.5e-19 * fac_cool;
  coeff[COOL_IV]  = val * 1.42e-27 * fac_cool;

#ifdef MRT_INCLUDE_HE
  coeff[COOL_V]           = val * 1.55e-26 * fac_cool;
  coeff[COOL_VI]          = val * 3.48e-26 * fac_cool;
  coeff[COOL_VII]         = val * 9.38e-22 * fac_cool;
  coeff[COOL_VIII]        = val * 4.95e-22 * fac_cool;
  coeff[COOL_IX]          = val * 5.01e-27 * pow(All.HubbleParam / All.UnitLength_in_cm, 3) * nH * fac_cool;
  coeff[COOL_X]           = val * 5.54e-17 * fac_cool;
  coeff[COOL_XI]          = val * 9.10e-27 * pow(All.HubbleParam / All.UnitLength_in_cm, 3) * nH * fac_cool;
  coeff[COOL_HEII_DIELEC] = val * 1.24e-13 * fac_cool;
#endif

  if(All.ComovingIntegrationOn)
    coeff[COOL_INV_COMPTON] = val * 5.62421e-36 * (All.UnitTime_in_s / All.UnitEnergy_in_cgs) / nH / pow(All.Time, 4);

#ifdef MRT_MOLECULAR_COOLING
  double nh_cgs  = nH / pow(All.UnitLength_in_cm, 3);
  double zfrac   = 1.0;  // correct this if metallicity included
  double fshield = fratio;
  //      printf("fshield =  %g \n", fshield) ;
  double f1 = (0.001 + 0.1 * nh_cgs / (1.0 + nh_cgs) + 0.09 * nh_cgs / (1.0 + 0.1 * nh_cgs) + zfrac * zfrac / (1.0 + nh_cgs));
  double f2 = (1.0 + zfrac) / (1.0 + 0.00143 * nh_cgs);
  double f  = f1 * f2 * fshield;
  coeff[MOLECULAR_COOLING] = 2.896e-26 * fac_cool * f * val;
#endif

#ifdef MRT_PHOTOELECTRIC_HEATING
  e_pe = therm.DensPhot[FOPT] * 1e63 * MeanPhotonEnergy[FOPT] * FRAC_FUV / fac_pe;
  //  e_pe = 2.0 ;
  double z_pe     = 1.0;  // correct this if metallicity included
  double pe_const = 1.3e-24 * fac_cool * val * z_pe * e_pe / nH;
  e_pe /= 0.5 * SphP[i].Ne * nH;
  coeff[PE_HEATING_I] = pe_const * 0.049;
  coeff[PE_HEATING_I] = pe_const * 0.037;
#endif

#ifdef MRT_METAL_COOLING
  double zmetal;
#ifdef GFM_COOLING_METAL
  zmetal = SphP[i].Metallicity / GFM_SOLAR_METALLICITY;  // correct this if metallicity included
#else
  zmetal = 1.0;
#endif
  coeff[METAL_COOLING] = fac_cool * val * zmetal;
#endif

#endif

#endif

  return;
}

static void set_coeffs_photons(int i, double nH, double rho, double T, double inv_cv)
{
  double sum_HI, sum_HeI, sum_HeII, cool_HI, cool_HeI, cool_HeII;
  sum_HI = sum_HeI = sum_HeII = cool_HI = cool_HeI = cool_HeII = 0.0;

  double val = nH * nH * inv_cv / rho;
  for(int kk = 0; kk < UV_BINS; kk++)
    {
      double gamma = therm.DensPhot[kk] * 1e63 / nH / pow(All.cf_atime, 3);

      sum_HI += nH * clight * mrt_sigma_HI[kk] * gamma;
#ifdef MRT_INCLUDE_HE
      sum_HeI += nH * clight * mrt_sigma_HeI[kk] * gamma;
      sum_HeII += nH * clight * mrt_sigma_HeII[kk] * gamma;
#endif

#ifdef MRT_COOLING_HEATING
      cool_HI += val * clight * mrt_sigma_HI[kk] * G_HI[kk] * gamma;
#ifdef MRT_INCLUDE_HE
      cool_HeI += val * clight * mrt_sigma_HeI[kk] * G_HeI[kk] * gamma;
      cool_HeII += val * clight * mrt_sigma_HeII[kk] * G_HeII[kk] * gamma;
#endif
#endif
    }

  coeff[PHOT_ION_HI]   = sum_HI;
  coeff[PHOT_ION_HEI]  = sum_HeI;
  coeff[PHOT_ION_HEII] = sum_HeII;

  coeff[PHOT_COOL_HI]   = cool_HI;
  coeff[PHOT_COOL_HEI]  = cool_HeI;
  coeff[PHOT_COOL_HEII] = cool_HeII;

#ifdef MRT_UVB
  double nh_cgs1 = nH / pow(All.UnitLength_in_cm, 3);
  double ss =
      (1 - ss_f) * pow(1 + pow(nh_cgs1 / pow(10, ss_n0), ss_beta), ss_alpha1) + ss_f * pow(1 + nh_cgs1 / pow(10, ss_n0), ss_alpha2);
  ss = 1 - ss_xi * (1 - ss);

  coeff[PHOT_ION_HI] += UVBionHI * ss;
  coeff[PHOT_ION_HEI] += UVBionHeI * ss;
  coeff[PHOT_ION_HEII] += UVBionHeII * ss;

  coeff[PHOT_COOL_HI] += val / nH * UVBcoolHI * ss;
  coeff[PHOT_COOL_HEI] += val / nH * UVBcoolHeI * ss;
  coeff[PHOT_COOL_HEII] += val / nH * UVBcoolHeII * ss;

#endif

  return;
}

static void set_out(int i, N_Vector y, double nH, double rho, double T, double dt)
{
  double nh2, nhe2, nhe3, temp, EIR;
  nh2  = Ith(y, 1);
  nhe2 = Ith(y, 2);
  nhe3 = Ith(y, 3), temp = Ith(y, 4);

  if(isnan(temp))
    terminate("Temp is nan %g %g %g %g \n", nh2, nhe2, nhe3, temp);

  double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

  if(nh2 < 0.0)
    nh2 = 1e-10;

  if(nh2 > 1.0)
    nh2 = 1.0;

  if(nhe2 < 0.0)
    nhe2 = 1e-10;

  if(nhe2 > y_fac)
    {
      nhe2 = 1.0;
      nhe3 = 1e-10;
    }

  if(nhe3 < 0.0)
    nhe3 = 1e-10;

  if(nhe3 > y_fac)
    {
      nhe3 = 1.0;
      nhe2 = 1e-10;
    }

  if(nhe2 + nhe3 > y_fac)
    {
      nhe2 *= y_fac / (nhe2 + nhe3);
      nhe3 *= y_fac / (nhe2 + nhe3);
    }

  if(temp < All.MinGasTemp)
    temp = All.MinGasTemp;

  SphP[i].nHII   = nh2 * P[i].Mass;
  SphP[i].nHeII  = nhe2 * P[i].Mass;
  SphP[i].nHeIII = nhe3 * P[i].Mass;

  SphP[i].ne   = SphP[i].nHII + SphP[i].nHeII + 2.0 * SphP[i].nHeIII;
  SphP[i].nHI  = (1.0 - nh2) * P[i].Mass;
  SphP[i].nHeI = (y_fac - nhe2 - nhe3) * P[i].Mass;

  SphP[i].HII   = nh2;
  SphP[i].HeII  = nhe2;
  SphP[i].HeIII = nhe3;
  SphP[i].Ne    = nh2 + nhe2 + 2.0 * nhe3;
  SphP[i].HI    = 1.0 - nh2;
  SphP[i].HeI   = y_fac - nhe2 - nhe3;

  therm.temp = temp;
  return;
}

#ifdef MRT_METAL_COOLING
static double metal_cooling(double T, int mode)
{
  double c1, c2, t1, t2, cval;
  c1 = c2 = t1 = t2 = cval = 0.0;

  double metal_cooling_val, metal_cooling_diff;

  metal_cooling_val = metal_cooling_diff = 0.0;

  double temp;

  temp = get_CoolingMetalRate(0.0, lognH, log10(T), &cval, &c1, &c2, &t1, &t2);
  if(temp != 0.0)
    metal_cooling_val = cval;
  if(t1 != t2)
    metal_cooling_diff = 1.0 / cval * (c2 - c1) / (pow(10.0, t2) - pow(10.0, t1));

  if(mode == 1)
    return metal_cooling_val;
  else
    return metal_cooling_diff;
}
#endif

#ifdef MRT_PHOTOELECTRIC_HEATING
static double pe_functional(double T, int mode)
{
  double val;
  double x_pe = e_pe * sqrt(T);
  if(mode == 1)
    val = 1.0 / (1.0 + pow(x_pe / 1925.0, 0.73));
  else
    val = pow(T / 1e4, 0.7) / (1.0 + x_pe / 5000.0);

  return val;
}

static double pe_differential(double T, int mode)
{
  double val;
  double x_pe = e_pe * sqrt(T);
  if(mode == 1)
    val = -0.5 * (0.73 / T) * pow(x_pe / 1925.0, 0.73) / (1.0 + pow(x_pe / 1925.0, 0.73));
  else
    val = (0.7 / T) - (0.5 * x_pe / 5000.0 / T) / (1.0 + x_pe / 5000.0);

  return val;
}
#endif

static double compton_functional(double T) { return (T - 2.727 / All.Time); }

static double compton_differential(double T) { return 1.0 / (T - 2.727 / All.Time); }

static double dielectric_functional(double T)
{
  double a = 1.5;
  double b = 470000.0;
  double c = 0.3;
  double d = 94000.0;

  return pow(T, -a) * exp(-b / T) * (1.0 + c * exp(-d / T));
}

static double dielectric_differential(double T)
{
  double a = 1.5;
  double b = 470000.0;
  double c = 0.3;
  double d = 94000.0;
  return c * d * exp(-d / T) / (T * T * (1.0 + c * exp(-d / T))) + b / T / T - a / T;
}

#ifdef MRT_MOLECULAR_COOLING
static double molecular_functional(double T)
{
  double l1 = -4.9202;
  double k1 = 125.215;
  double l2 = -1.7288;
  double k2 = 1349.86;
  double l3 = -0.3075;
  double k3 = 6450.06;
  double k4 = 158000.0;

  double zi = 1.0 / (pow(T / k1, l1) + pow(T / k2, l2) + pow(T / k3, l3));
  double x  = -(T / k4) * (T / k4);

  return zi * exp(x);
}

static double molecular_differential(double T)
{
  double l1 = -4.9202;
  double k1 = 125.215;
  double l2 = -1.7288;
  double k2 = 1349.86;
  double l3 = -0.3075;
  double k3 = 6450.06;
  double k4 = 158000.0;

  double zi = 1.0 / (pow(T / k1, l1) + pow(T / k2, l2) + pow(T / k3, l3));

  double dziT = (l1 / k1) * pow(T / k1, l1 - 1) + (l2 / k2) * pow(T / k2, l2 - 1) + (l3 / k3) * pow(T / k3, l3 - 1);
  double dx   = 2.0 * T / k4 / k4;

  return -dziT * zi - dx;
}
#endif

static double functional_form(double T, double k1, double k2, double k3, double k4, double k5, double k6, double l, double m, double n)
{
  double l1 = k1;
  double l2, l3, l4, l5;

  if(l == 0.0)
    l2 = 1.0;
  else
    l2 = pow(T / k2, l);

  if(m == 0.0)
    l3 = 1.0;
  else
    l3 = pow(T / k3, m);

  if(k4 == 0.0)
    l4 = 1.0;
  else
    l4 = exp(-k4 / T);

  l5 = 1.0 + k5 * pow(T / k6, n);

  return l1 * l2 * l3 * l4 / l5;
}

static double differential_form(double T, double k1, double k2, double k3, double k4, double k5, double k6, double l, double m,
                                double n)
{
  double l1 = (l + m) / T;
  double l2 = k4 / (T * T);
  double l3 = -k5 * n * pow(T, n - 1) / (pow(k6, n) + k5 * pow(T, n));

  return l1 + l2 + l3;
}

static double functional(int num, double T)
{
  double k1, k2, k3, k4, k5, k6, l, m, n;
  double val;
  k1 = 1.0;
  if(num == ALPHA_HII)
    {
      k2 = 1e4;
      k3 = 1.0;
      k4 = 0.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = -0.7;
      m  = 0.0;
      n  = 0.0;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == SIGMA_HI)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 157809.1;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == ALPHA_HEIII)
    {
      k2 = 1.0;
      k3 = 1e3;
      k4 = 0.0;
      k5 = 1.0;
      k6 = 1e6;
      l  = -0.5;
      m  = -0.2;
      n  = 0.7;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == SIGMA_HEI)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 285335.4;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == ALPHA_HEII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 0.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = -0.6353;
      m  = 0.0;
      n  = 0.0;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == SIGMA_HEII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 631515.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_I)
    {
      k2 = 1.0;
      k3 = 1e3;
      k4 = 0.0;
      k5 = 1.0;
      k6 = 1e6;
      l  = 0.5;
      m  = -0.2;
      n  = 0.7;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_II)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 157809.1;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_III)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 118348.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.0;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_IV)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 0.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = 0.5;
      m  = 0.0;
      n  = 0.0;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_V)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 1.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = 0.3647;
      m  = 0.0;
      n  = 0.0;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_VI)
    {
      k2 = 1.0;
      k3 = 1e3;
      k4 = 0.0;
      k5 = 1.0;
      k6 = 1e6;
      l  = 0.5;
      m  = -0.2;
      n  = 0.7;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_VII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 285335.4;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_VIII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 631515.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_IX)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 55338.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = -0.1687;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_X)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 473638.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = -0.397;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_XI)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 13179.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = -0.1681;
      m  = 0.0;
      n  = 0.5;

      val = functional_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
#ifdef MRT_MOLECULAR_COOLING
  else if(num == MOLECULAR_COOLING)
    {
      val = molecular_functional(T);
    }
#endif
#ifdef MRT_PHOTOELECTRIC_HEATING
  else if(num == PE_HEATING_I)
    {
      val = pe_functional(T, 1);
    }
  else if(num == PE_HEATING_II)
    {
      val = pe_functional(T, 2);
    }
#endif
#ifdef MRT_METAL_COOLING
  else if(num == METAL_COOLING)
    val = metal_cooling(T, 1);
#endif
  else if(num == COOL_INV_COMPTON)
    {
      if(All.ComovingIntegrationOn)
        val = compton_functional(T);
      else
        val = 1.0;
    }
  else if(num == ALPHA_HEII_DIELEC || num == COOL_HEII_DIELEC)
    {
      val = dielectric_functional(T);
    }
  else
    val = 1.0;

  return val;
}

static double differential(int num, double T)
{
  double k1, k2, k3, k4, k5, k6, l, m, n;
  double val;
  k1 = 1.0;
  if(num == ALPHA_HII)
    {
      k2 = 1e4;
      k3 = 1.0;
      k4 = 0.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = -0.7;
      m  = 0.0;
      n  = 0.0;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == SIGMA_HI)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 157809.1;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == ALPHA_HEIII)
    {
      k2 = 1.0;
      k3 = 1e3;
      k4 = 0.0;
      k5 = 1.0;
      k6 = 1e6;
      l  = -0.5;
      m  = -0.2;
      n  = 0.7;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == SIGMA_HEI)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 285335.4;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == ALPHA_HEII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 0.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = -0.6353;
      m  = 0.0;
      n  = 0.0;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == SIGMA_HEII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 631515.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_I)
    {
      k2 = 1.0;
      k3 = 1e3;
      k4 = 0.0;
      k5 = 1.0;
      k6 = 1e6;
      l  = 0.5;
      m  = -0.2;
      n  = 0.7;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_II)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 157809.1;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_III)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 118348.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.0;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_IV)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 0.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = 0.5;
      m  = 0.0;
      n  = 0.0;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_V)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 1.0;
      k5 = 0.0;
      k6 = 1.0;
      l  = 0.3647;
      m  = 0.0;
      n  = 0.0;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_VI)
    {
      k2 = 1.0;
      k3 = 1e3;
      k4 = 0.0;
      k5 = 1.0;
      k6 = 1e6;
      l  = 0.5;
      m  = -0.2;
      n  = 0.7;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_VII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 285335.4;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_VIII)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 631515.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = 0.5;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_IX)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 55338.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = -0.1687;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_X)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 473638.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = -0.397;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
  else if(num == COOL_XI)
    {
      k2 = 1.0;
      k3 = 1.0;
      k4 = 13179.0;
      k5 = 1.0;
      k6 = 1e5;
      l  = -0.1681;
      m  = 0.0;
      n  = 0.5;

      val = differential_form(T, k1, k2, k3, k4, k5, k6, l, m, n);
    }
#ifdef MRT_MOLECULAR_COOLING
  else if(num == MOLECULAR_COOLING)
    {
      val = molecular_differential(T);
    }
#endif
#ifdef MRT_PHOTOELECTRIC_HEATING
  else if(num == PE_HEATING_I)
    {
      val = pe_differential(T, 1);
    }
  else if(num == PE_HEATING_II)
    {
      val = pe_differential(T, 2);
    }
#endif
#ifdef MRT_METAL_COOLING
  else if(num == METAL_COOLING)
    val = metal_cooling(T, 2);
#endif
  else if(num == COOL_INV_COMPTON)
    {
      if(All.ComovingIntegrationOn)
        val = compton_differential(T);
      else
        val = 0.0;
    }
  else if(num == ALPHA_HEII_DIELEC || num == COOL_HEII_DIELEC)
    {
      val = dielectric_differential(T);
    }
  else
    val = 0.0;
  return val;
}

static int mrt_thermochem_semi_implicit(int i, double dt, double T)
{
  double nh1, nh2, nhe1, nhe2, nhe3, ne;

  double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

  double val[NCOEFFS];
  for(int kk = 0; kk < NCOEFFS; kk++)
    val[kk] = coeff[kk] * functional(kk, T);

  val[ALPHA_HEII] += val[ALPHA_HEII_DIELEC];
  val[COOL_V] += val[COOL_HEII_DIELEC];

  double A = dt * val[SIGMA_HI] * SphP[i].Ne;
  double B = dt * val[PHOT_ION_HI];
  double C = dt * val[ALPHA_HII] * SphP[i].Ne;

  nh2 = SphP[i].HII + B + A;
  nh2 /= 1.0 + B + C + A;

  if(nh2 < 0.0 || nh2 > 1.0 || isnan(nh2))
    return 1;

  double D = dt * val[SIGMA_HEII] * SphP[i].Ne;
  double E = dt * val[ALPHA_HEIII] * SphP[i].Ne;
  double F = dt * val[SIGMA_HEI] * SphP[i].Ne;
  double J = dt * val[ALPHA_HEII] * SphP[i].Ne;
  double G = dt * val[PHOT_ION_HEI];
  double L = dt * val[PHOT_ION_HEII];

#ifdef MRT_INCLUDE_HE
  nhe2 = SphP[i].HeII / y_fac;
  nhe3 = SphP[i].HeIII / y_fac;
#else
  nhe2 = 0.0;
  nhe3 = 0.0;
#endif

  nhe2 = nhe2 + G + F - ((G + F - E) / (1.0 + E)) * nhe3;
  nhe2 /= 1.0 + G + F + D + J + L + ((G + F - E) / (1.0 + E)) * (D + L);

  if(nhe2 < 0.0 || nhe2 > 1.0 || isnan(nhe2))
    return 1;

  nhe3 = nhe3 + (D + L) * nhe2;
  nhe3 /= 1.0 + E;

  if(nhe3 < 0.0 || nhe3 > 1.0 || isnan(nhe3))
    return 1;

  if(nhe2 + nhe3 > 1.0)
    return 1;

  nhe2 *= y_fac;
  nhe3 *= y_fac;

  nh1  = 1.0 - nh2;
  nhe1 = y_fac - nhe2 - nhe3;
  ne   = nh2 + nhe2 + 2.0 * nhe3;

  double inv_cv = therm.inv_cv;

  double dT = dt * (val[PHOT_COOL_HI] * nh1 + val[PHOT_COOL_HEI] * nhe1 + val[PHOT_COOL_HEII] * nhe2 - val[COOL_I] * nh2 * ne -
                    (val[COOL_II] + val[COOL_III]) * nh1 * ne - (val[COOL_V] + val[COOL_VIII] + val[COOL_X]) * nhe2 * ne -
                    val[COOL_VI] * nhe3 * ne - val[COOL_VII] * nhe1 * ne - (val[COOL_IX] + val[COOL_XI]) * nhe2 * ne * ne -
                    val[COOL_IV] * (nh2 + nhe2 + 4.0 * nhe3) * ne
#ifdef MRT_MOLECULAR_COOLING
                    - val[MOLECULAR_COOLING]
#endif
#ifdef MRT_PHOTOELECTRIC_HEATING
                    + val[PE_HEATING_I] + val[PE_HEATING_II]
#endif
#ifdef MRT_METAL_COOLING
                    - val[METAL_COOLING]
#endif
                   );

  if(fabs(dT / T) > 0.1 || isnan(dT))
    return 1;

  SphP[i].nHII   = nh2 * P[i].Mass;
  SphP[i].nHeII  = nhe2 * P[i].Mass;
  SphP[i].nHeIII = nhe3 * P[i].Mass;

  SphP[i].ne   = SphP[i].nHII + SphP[i].nHeII + 2.0 * SphP[i].nHeIII;
  SphP[i].nHI  = (1.0 - nh2) * P[i].Mass;
  SphP[i].nHeI = (y_fac - nhe2 - nhe3) * P[i].Mass;

  SphP[i].Ne  = SphP[i].ne / P[i].Mass;
  SphP[i].HI  = SphP[i].nHI / P[i].Mass;
  SphP[i].HeI = SphP[i].nHeI / P[i].Mass;

  SphP[i].HII   = nh2;
  SphP[i].HeII  = nhe2;
  SphP[i].HeIII = nhe3;

  therm.temp = T + dT;

  return 0;
}

#ifndef MRT_NO_UV
static void absorb_photons(int i, double dt, double nH, double T)
{
  double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;

  double nh2 = SphP[i].HII;
#ifdef MRT_INCLUDE_HE
  double nhe2 = SphP[i].HeII;
  double nhe3 = SphP[i].HeIII;
#endif
  double KK = dt * clight * nH;
  for(int j = 0; j < UV_BINS; j++)
    {
      double sum_KK = mrt_sigma_HI[j] * (1.0 - nh2)
#ifdef MRT_INCLUDE_HE
                      + mrt_sigma_HeI[j] * (y_fac - nhe2 - nhe3) + mrt_sigma_HeII[j] * nhe2
#endif
          ;

#ifdef MRT_IR
      double AA = SphP[i].KappaIR_R[j] / nH;
      sum_KK += AA;

      if(KK * AA > 1.0)
        therm.Cons_DensPhot[UV_BINS] += MeanPhotonEnergy[j] * therm.Cons_DensPhot[j] * 1e63;
      else
        therm.Cons_DensPhot[UV_BINS] += KK * AA * MeanPhotonEnergy[j] * therm.Cons_DensPhot[j] * 1e63;
#endif

      double ratio = exp(-KK * sum_KK);

#ifdef MRT_MOLECULAR_COOLING
      if(nu[j] == 13.6)
        {
          double abs1 = clight * mrt_sigma_HI[j] * therm.DensPhot[j] * 1e63;
#ifdef MRT_UVB
          abs1 += UVBionHI;
#endif
          abs1 *= 1e12 / All.UnitTime_in_s;

          double abs2 = 4.4 * 3.0856776e+18 / All.UnitLength_in_cm;
          double tau  = mrt_sigma_HI[j] * nH * abs2 * pow(T / 1e4, -0.173) * pow(abs1, -0.66667);
          fratio      = 1.0 - exp(-tau);
        }
#endif

      therm.Cons_DensPhot[j] *= ratio;
    }
#ifdef MRT_IR_ONLY_CHEMISTRY

  double A      = dt * SphP[i].KappaIR_R[UV_BINS] * clight;
  double ratio1 = 1.0 / (1 + A);

  for(int j = 0; j < 3; j++)
    {
      SphP[i].Cons_RT_F[UV_BINS][j] *= ratio1;
      if(isnan(SphP[i].Cons_RT_F[UV_BINS][j]))
        terminate("IR Cooling - What!!\n");
    }
#endif
  return;
}
#endif

#endif
