/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/MRT/RT_proto.h
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef MRT_RT_PROTO_H
#define MRT_RT_PROTO_H

#include "../allvars.h"

#ifdef MRT

void mrt_run(void);
void mrt_run_sources(void);
void exchange_primitive_variables_RT(void);
void exchange_primitive_variables_and_gradients_RT(void);
void update_primitive_variables_RT(void);
void calculate_gradients_RT(void);
void compute_interface_fluxes_RT(tessellation *T);
int face_get_state_RT(tessellation *T, int p, int i, struct state *st);
void state_convert_to_local_frame_RT(struct state *st, double *vel_face, double hubble_a, double atime);
void face_do_time_extrapolation_RT(struct state *delta, struct state *st, double atime);
void face_do_spatial_extrapolation_RT(struct state *delta, struct state *st, struct state *st_other);
void face_add_extrapolations_RT(struct state *st_face, struct state *delta_time, struct state *delta_space);
void face_add_extrapolation_RT(struct state *st_face, struct state *delta);
void face_turn_velocities_RT(struct state *st, struct geometry *geom);
void face_turnback_velocities_RT(struct state_face *st_face, struct geometry *geom);
void face_limit_fluxes_RT(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R,
                          struct fluxes *flux, double dt, double *count, double *count_reduced);
void apply_flux_list_RT(void);
double face_timestep_RT(struct state *state_L, struct state *state_R, double *hubble_a, double *atime);

double get_spectrum(double e, int i_age, int i_metallicity);

void mpi_printf_rt(const int flag, const char *fmt, ...);

void init_RT(void);
void RT_initialize_cell(int i);
void add_source_fluxes(void);
void set_VET_single(int, struct sph_particle_data *);
void mrt_update_chemistry(void);

void init_gradients_RT(void);
void gradient_init_RT(MyFloat *addr, MyFloat *addr_exch, MySingle *addr_grad, int type);
void face_turn_momentum_flux_RT(struct fluxes *flux, struct geometry *geom);

void set_full_ionization_mrt(int);

#ifdef MRT_RIEMANN_ROSUNOV
double godunov_flux_3d_rosunov_RT(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux,
                                  double vel_face_x);
#endif

#ifdef MRT_RIEMANN_ROSUNOV_NEW
double godunov_flux_3d_rosunov_RT_new(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux,
                                      double vel_face);
#endif

#ifdef MRT_RIEMANN_HLLE
double godunov_flux_3d_HLLE_RT(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux,
                               double vel_face_x);
#endif

#ifdef MRT_RIEMANN_HLLE_NEW
double godunov_flux_3d_HLLE_RT_new(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux,
                                   double vel_face);
#endif

#if defined(MRT_RIEMANN_HLLE) || defined(MRT_RIEMANN_HLLE_NEW)
void readin_hlle_eingenvalues(void);
#endif

#if defined(MRT_RIEMANN_ROSUNOV_NEW) || defined(MRT_RIEMANN_HLLE_NEW)
void flux_convert_to_lab_frame_RT_new(const struct state_face *st, double *vel_face, struct fluxes *flux);
#endif

#ifdef MRT_SINGLE_STAR
void read_stellar_table(void);
#endif

#ifdef MRT_COMOVING
void cell_do_lorentz_boost(int, struct sph_particle_data *, double, double, double);
void do_comoving_frame_source_terms(void);
#endif

#ifdef MRT_IR
#ifdef MRT_IR_GRAIN_KAPPA
void read_grain_kappa_data(void);
#endif
double mrt_update_IR_cooling(int, double);
#ifdef MRT_IR_LTE_GSL
int mrt_IR_rate_ODEs(double, const double *, double *, void *);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
#endif
#endif

#ifdef MRT_LOCAL_FEEDBACK
double inject_photons_from_star(int, double, double);
#endif

#if defined(MRT_IR) || defined(MRT_UV_ONLY_DUST) || defined(MRT_CHEM_SG)
void set_kappa_times_rho_IR(int, struct sph_particle_data *);
#endif

#ifdef MRT_UV_ONLY_DUST
void mrt_update_chemistry_only_dust(void);
#endif

#ifndef MRT_NO_UV
void initialize_ionic_species(void);
#endif

#ifdef MRT_RADIATION_PRESSURE
void do_radiation_pressure_source_terms(void);
#endif

#ifdef COOLING
int do_cooling_mrt(int i);
#ifdef MRT_EQUIL_CHEM_COOL
void update_radiation_state_mrt(int i, PhotoCurrent *ppc);
void update_chem_mrt(int i, double dtcool, GasState *pgs);
#endif
#endif

#if defined(MRT_COOLING_HEATING)
double mrt_DoHeating(int, double);
double mrt_DoCooling(int, double);
double mrt_get_cooling_rate(int, double);
double mrt_get_heating_rate(int);
double mrt_GetCoolingTime(int, double, double, double *);
#endif

#ifdef MRT_COUPLED_THERMOCHEMISTRY
void mrt_thermochemistry(void);
#endif

#ifdef MRT_CHEMISTRY_PS2009
void mrt_update_chemistry_ps2009(void);
void mrt_write_stats(void);
void mrt_IR_chemistry(void);
#endif

#ifdef MRT_IR_ONLY_CHEMISTRY
void mrt_IR_chemistry(void);
#endif

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
void mrt_setup(int);
void exchange_vector(void);
#endif

#ifdef MRT_CHEMISTRY_PS2011
int mrt_rate_ODEs(double, const double *, double *, void *);
void mrt_update_chemistry_ps2011(void);
#endif

#ifdef MRT_CHEM_SG
void set_rates_sgchem(void);
#endif

#ifndef MRT_NO_UV
void mrt_get_sigma(void);
void mrt_get_luminosities(void);
#endif

#ifdef MRT_TIME_EXTRAPOLATION
void calculate_div_F(struct state *dl, struct state *st);
#endif

#if defined(MRT_SOURCES) || defined(MRT_LOCAL_FEEDBACK)
void mrt_load_spectrum_table(void);
void mrt_free_spectrum_table(void);

#ifdef MRT_STARS
void start_stellar_sources(void);
void end_stellar_sources(void);
void do_ionizing_stellar_sources(void);
void add_ionizing_stellar_radiation(void);
#endif

#ifdef MRT_BH
void start_blackhole_sources(void);
void end_blackhole_sources(void);
void do_ionizing_blackhole_sources(void);
void add_ionizing_blackhole_radiation(void);
#ifdef MRT_BH_BIPOLAR
int is_cell_within_cone(MyIDType, MyDouble *);
#endif
#endif

#endif /* defined(MRT_SOURCES) || defined(MRT_LOCAL_FEEDBACK) */

static inline double interpolate_2d(double (*table)[101], int i, int j, double dx, double dy)
{
  return (1 - dx) * (1 - dy) * table[i][j] + (1 - dx) * dy * table[i][j + 1] + dx * (1 - dy) * table[i + 1][j] +
         dx * dy * table[i + 1][j + 1];
}

static inline MyFloat interpolate_4d(MyFloat ****table, int i, int j, int k, int l, double dx, double dy, double dz, double dw)
{
  int il = i, jl = j, kl = k, ll = l;
  int ir = i + 1, jr = j + 1, kr = k + 1, lr = l + 1;

  double dxl = 1 - dx, dyl = 1 - dy, dzl = 1 - dz, dwl = 1 - dw;
  double dxr = dx, dyr = dy, dzr = dz, dwr = dw;
  if(dxr == 0)
    ir = i;

  if(dyr == 0)
    jr = j;

  if(dzr == 0)
    kr = k;

  if(dwr == 0)
    lr = l;

  return dxl * dyl * dzl * dwl * table[il][jl][kl][ll] + dxl * dyl * dzl * dwr * table[il][jl][kl][lr] +
         dxl * dyl * dzr * dwl * table[il][jl][kr][ll] + dxl * dyl * dzr * dwr * table[il][jl][kr][lr] +
         dxl * dyr * dzl * dwl * table[il][jr][kl][ll] + dxl * dyr * dzl * dwr * table[il][jr][kl][lr] +
         dxl * dyr * dzr * dwl * table[il][jr][kr][ll] + dxl * dyr * dzr * dwr * table[il][jr][kr][lr] +
         dxr * dyl * dzl * dwl * table[ir][jl][kl][ll] + dxr * dyl * dzl * dwr * table[ir][jl][kl][lr] +
         dxr * dyl * dzr * dwl * table[ir][jl][kr][ll] + dxr * dyl * dzr * dwr * table[ir][jl][kr][lr] +
         dxr * dyr * dzl * dwl * table[ir][jr][kl][ll] + dxr * dyr * dzl * dwr * table[ir][jr][kl][lr] +
         dxr * dyr * dzr * dwl * table[ir][jr][kr][ll] + dxr * dyr * dzr * dwr * table[ir][jr][kr][lr];
}

#endif

#endif
