/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sfr_mcs_proto.h
 * \date        10/2015
 * \author     	Matthew C Smith
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef SFR_MCS_PROTO_H
#define SFR_MCS_PROTO_H

void init_star_formation(void);
void sfr_criteria_announce(void);
#ifndef SN_MCS_INITIAL_DRIVING
void check_AuxDataID_references_mcs(void);
void sfr_mcs_add_star(int i, int j, MyDouble mass_of_star, MyFloat birthtime);
#endif

#ifdef SFR_MCS_LOG
void setup_sf_log(void);
void sf_add_to_log(double rho);
void write_sf_dens_log(void);
#endif

#ifdef IMF_SAMPLING_MCS
void init_imf(void);
#ifdef IMF_SAMPLING_MCS_VERBOSE
MyFloat sample_imf(MyFloat *results, MyFloat target);
#else
void sample_imf(MyFloat *results, MyFloat target);
#endif
void do_imf_sampling(int old_N_star);
void set_star_properties(int old_N_star);
void set_individual_star_properties(int auxid, int nslot);
void get_star_mass_interp_index(MyFloat m, MyFloat *masses, int N_m, int *i_low, MyFloat *delta);
void check_for_dead_stars(void);
void init_star_properties_tables(void);
#endif

#ifdef GFM_COOLING_METAL

#define GFM_MIN_METAL -20 /* This is otherwise defined in stellar_evolution_vars.h */

void metal_cool_mcs_init(double localMetallicityFloor, double localLogMetallicityFloorInSolar);
#ifdef UVB_SELF_SHIELDING
void update_radiation_state(MyFloat rho, MyFloat metallicity, PhotoCurrent *localpc);
#endif
void update_gas_state(MyFloat rho, MyFloat metallicity, GasState *localgs);

/* This function declared in stellar_evolution_proto.h but needed in cooling_metal.c, so reproduced here. */
static inline MyFloat interpol_4d(MyFloat ****table, int i, int j, int k, int l, double dx, double dy, double dz, double dw)
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

#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)
void init_stellar_feedback(void);
void do_stellar_feedback(void);
#endif

#if(defined(SN_MCS) && !defined(SN_MCS_INITIAL_DRIVING)) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)
void update_stellar_ages(void);
#endif

#if((defined(SN_MCS) && !defined(SN_MCS_SINGLE_INJECTION)) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)) && \
    !(defined(SN_MCS_INITIAL_DRIVING) || defined(IMF_SAMPLING_MCS))
void read_sb99_tables(void);
#ifndef SB99_FIXED_Z
int get_sb99_z_index(MyFloat metallicity);
#endif
#endif

#ifdef SN_MCS
void init_sne(void);

#ifndef SN_MCS_INITIAL_DRIVING
#ifndef IMF_SAMPLING_MCS
void check_for_supernovae(int *local_n_sn_event, int *global_n_sn_event);
void get_sb99_t_indicies(MyFloat t, MyFloat *timesteps, int N_t, int *it_low, int *it_high, MyFloat *delta);
#ifndef SB99_FIXED_Z
double get_sn_rate(MyFloat t, int iz);
#else
double get_sn_rate(MyFloat t);
#ifdef SN_MCS_VARIABLE_EJECTA
void get_ejecta_properties(MyFloat t, double *m_ej, double *z_ej);
#endif
#endif
#endif
void find_sn_host_cells_and_distribute(int n_sn_event);
#else
int check_for_supernovae();
#endif

void stellar_feedback_distribute(void);
void sn_finish(void);
void do_supernovae(void);

#ifdef SN_MCS_LOG
void setup_sn_log(void);
void sn_add_to_log(int i);
void write_sn_dens_log(void);
#endif  // SN_MCS_LOG

#ifdef IMF_SAMPLING_MCS
void snII_yields(MyFloat mstar, MyFloat *m_ej, MyFloat *m_Z);
#endif

#ifdef TRACER_MC
int consider_moving_tracers_ejecta_local(int p_from, int p_to, double prob);
int consider_moving_tracers_ejecta(int p, int pother_task, int pother_index, double prob);
void finalise_host_tracers(int p);
void add_tracer_to_TFL(int p_orig, int tr_ind, int pother_task, int pother_index, MyIDType pother_ID);
#endif
#endif  // SN_MCS

#ifdef HII_MCS
void init_hii(void);
#ifndef HII_MCS_TEST
int update_photon_rates(void);
#endif
int set_star_particle_ionization_properties(void);
void place_hii_regions(void);
void place_hii_regions_internal(void);
#ifdef HII_MCS_ANISO
void vec2pix_ring(const double *vec, double *vlen, int *pix);
#endif
#ifdef HII_MCS_LR
double calculate_uv_background_boost_factor(int i);
#endif
#endif

#ifdef PE_MCS
void init_pe(void);
void update_FUV_luminosities(void);
double calculate_pe_heating_rate(int i);
double get_DGR_relative_solar(double met);
#if(PE_MCS_PRESHIELD == 1)
void estimate_dust_column(dust_column_search_data *searchdata_input, int nn, int verbose);
#endif
#endif

#ifdef TURB_APPROX_MCS
void turb_approx_init(void);
void update_turbulent_energy(void);
#ifdef TURB_APPROX_MCS_OUTPUT_RATES
double get_turbulent_production_rate_single(int i);
#endif
#endif 

#endif  // SFR_MCS_PROTO_H
