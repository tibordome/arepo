/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/braginskii_viscosity/braginskii_viscosity.h
 * \date        12/2017
 * \author      Thomas Berlok
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef BRAGINSKII_VISCOSITY_H
#define BRAGINSKII_VISCOSITY_H

#include "../allvars.h"

/* Braginskii viscosity */
#ifdef BRAGINSKII_VISCOSITY

// Configuration consistency checks
#if !defined(MHD)
#error Braginskii viscosity: Config option MHD is required!
#endif

#if !defined(VORONOI_MESH_KEEP_DT_AND_DTC)
#error Braginskii viscosity: Config option VORONOI_MESH_KEEP_DT_AND_DTC is required.
#endif

#if !defined(TETRA_INDEX_IN_FACE)
#error Braginskii viscosity: Config option TETRA_INDEX_IN_FACE is required.
#endif

#if defined(BRAGINSKII_RKL2_SUPER_TIME_STEPPING) && defined(BRAGINSKII_VISCOSITY_SUBCYCLE)
#error Braginskii viscosity: Subcycling and super timestepping cannot be used at the same time.
#endif

#if !defined(BRAGINSKII_RKL2_SUPER_TIME_STEPPING) && !defined(BRAGINSKII_VISCOSITY_SUBCYCLE)
#ifdef BRAGINSKII_OUTPUT_NUM_SUBSTEPS
#error Braginskii viscosity: BRAGINSKII_OUTPUT_NUM_SUBSTEPS requires either BRAGINSKII_RKL2_SUPER_TIME_STEPPING or BRAGINSKII_VISCOSITY_SUBCYCLE to be turned on.
#endif
#endif

// #if defined(BRAGINSKII_RKL2_SUPER_TIME_STEPPING) && !defined(FORCE_EQUAL_TIMESTEPS)
// #error Braginskii viscosity: Super timestepping is unfortunately incompatible with local timestepping.
// #endif

// #if defined(BRAGINSKII_RKL2_SUPER_TIME_STEPPING) && !defined(BRAGINSKII_VISCOSITY_GLOBAL_TIMESTEP)
// #error Braginskii viscosity: Super time stepping (STS) only implemented for global time steps
// #endif

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
#error Braginskii viscosity: Reflective boundaries are not implemented
#endif

struct diff_face_data
{
  int cornerFirst;
  // TODO: Reduce cornerCount, failCount and failWeight to fewer variables?
  int cornerCount;
  int failCount;
  double bfld[3];  // TODO: We could save a little of memory: double bfld[2];
  int sng;

  int active;
  double dt;
  double area;
  double eta_aniso;
  // double density;

#ifdef BRAGINSKII_LIMIT_PRESSURE_ANISOTROPY
  // Magnetic field strength
  double B2;
#endif

  double nx, ny, nz;
  double mx, my, mz;
  double px, py, pz;

  double failWeight;
};

struct corner_list
{
  int index;
  double weight;
};

// Struct for storing the velocity and its tensor at the corners of cells
// These are not included in corner_data because they change in time
struct corner_vel
{
  // Velocity components
  double vx;
  double vy;
  double vz;
  // Gradients of velocity components
  double vxGrad[3];
  double vyGrad[3];
  double vzGrad[3];
};

// Struct for storing the velocity at the center of cells
struct center_vel
{
  // Velocity components
  double vx;
  double vy;
  double vz;
};

// Struct for communicating viscous fluxes between tasks
struct viscflux_list_data
{
  int task;
  int index;
  double dEnergy;
  double dMx;
  double dMy;
  double dMz;
};

// Struct for storing several copies of the velocities at different time steps.
// For super time stepping, we need to have copies of the velocity at several
// substeps between each super step.
struct sphp_copy
{
  double Energy;
  double Momentum[3];
};

struct grad_v
{
  double vn;
  double vm;
  double dvndn;
  double dvmdn;
  double dvndm;
  double dvmdm;
  double dvpdp;
};

void do_braginskii_viscosity(void);
void set_diffusion_coefficients(struct diff_face_data *diff_face_data);

static void add_corner(int tetra, int iface);
static int get_timebin_of_point(int point);
static double get_surface_area_of_cell(int point);

static int prepare_stuff(void);
static void compute_least_squares_matrix_at_corners(void);
static void compute_geometry_of_interface(int iface);
static void compute_magnetic_field_at_interface(int iface);
static void compute_quantities_at_corners(struct sphp_copy *sphp_copy, struct sphp_copy *sphp_copyExch);
#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
static void compute_velocity_at_centers(struct sphp_copy *state, struct sphp_copy *stateExch);
#endif
static void diffusion_explicit(int substeps);
#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
static void diffusion_explicit_sts(int substeps);
#endif
static void compute_diffusion_coefficients(void);
static int viscflux_list_data_compare(const void *a, const void *b);
static void limit_flux(void);
static void free_stuff(void);
static void point_get_center(int p, double *Center);
static void exchange_and_apply_fluxes(struct viscflux_list_data *ViscFluxList, int Nflux, struct sphp_copy *sphp_copy);
static void calculate_grad_v_in_local_coordinates(struct diff_face_data *fd, struct grad_v *grad_v);
static void calculate_grad_v_in_local_coordinates_without_limiters(struct diff_face_data *fd, struct grad_v *grad_v);
#ifdef BRAGINSKII_SIMPLE_DERIVATIVES
static void calculate_simple_normal_gradients(struct diff_face_data *fd, int iface, struct grad_v *grad_v);
#endif
static int apply_flux_or_calculate_fluxlist(struct viscflux_list_data *ViscFluxList, struct sphp_copy *sphp_copy,
                                            struct diff_face_data *fd, double flux_E, double flux_Mx, double flux_My, double flux_Mz,
                                            int iface, int *MaxNflux, int Nflux);
static void calculate_and_add_fluxes(struct sphp_copy *state1, struct sphp_copy *state2, double factor, double factorE);
static void calculate_flux_across_face(struct diff_face_data *fd, struct grad_v *grad_v, double *flux_E, double *flux_Mx,
                                       double *flux_My, double *flux_Mz, double factor, double factorE);
static void calculate_isotropic_flux_across_face(struct diff_face_data *fd, double *flux_E, double *flux_Mx, double *flux_My,
                                                 double *flux_Mz, double factor, double factorE);
static double limited_derivative(struct diff_face_data *fd, const int mycase);
static double limited_derivative_iso(struct diff_face_data *fd, const int mycase);
static void copy_SphP(struct sphp_copy *sphp_copy);
static void set_SphP(struct sphp_copy *sphp_copy);
static void copy_state(struct sphp_copy *state1, struct sphp_copy *state2);
#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
static void find_rkl_weights(double *mu, double *mu_til, double *nu, double *gam_til, int s);
#endif  // BRAGINSKII_RKL2_SUPER_TIME_STEPPING
static void exchange_state_variables(struct sphp_copy *state, struct sphp_copy *stateExch);
#endif

#endif /*BRAGINSKII_VISCOSITY_H*/
