/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sne/sne_proto.h
 * \date        MM/YYYY
 * \author      Robin Tress, Rowan Smith, Andre Bubel
 * \brief
 * \details     Please contact the authors at robin.tress@uni-heidelberg.de
 *              and rowan.smith@manchester.ac.uk before using this to avoid overlapping
 *              of projects. And please report any issue you may encounter by using this
 *              routine.
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

void sne_feedback(void);

int do_sn_energy_injection(double sne_pos[3], int local_n, int indices[], double radius, double mean_density, double total_mass,
                           double sn_E51);

int decide_injection_scheme(double mean_density, double radius);
void energy_injection(int local_n, int indices[], double total_mass, double sn_E51);
void momentum_injection(double sne_pos[3], int local_n, int indices[], double mean_density, double total_mass);

void update_timesteps_around_injection_region(double sne_pos[3], double injectionradius);

void find_gas_properties_of_injection_region(int local_n, int indices[], double *total_mass, double *total_volume,
                                             double *mean_density, double *mean_utherm, double *mean_temp, double vel_CM[3]);
void distribute_mass(int local_n, int indices[], double mean_density, double mean_utherm, double vel_CM[3]);

void find_particles_within_a_sphere(double center[3], double radius, int *local_n, int *total_n, int **indices);
void find_injection_cells(double center[3], double *radius, int *local_n, int *total_n, int **indices, int sn_size);

enum SNE_type is_it_time_for_a_new_sn(void);
void determine_position_and_injection_region_of_sn(double sne_pos[3], enum SNE_type sne_type, double *radius, int *local_n,
                                                   int *total_n, int **indices, int sn_size);
int is_this_position_viable(double sne_pos[3], double radius);

void remove_sink_particle(int sink_index);

double find_temperature_of_cell(int idx);

void inject_mass(enum SNE_type sne_type, int local_n, int indices[], double total_volume, double *total_mass, double *mean_density,
                 double sn_mass);

enum SNE_type sne_only_once_at_center_TIMING(void);
void sne_only_once_at_center_POSITIONING(double sne_pos[3]);

enum SNE_type sne_random_TIMING(void);
void sne_random_POSITIONING(double sne_pos[3]);

enum SNE_type sne_at_cluster_position_TIMING(void);
void sne_at_cluster_position_POSITIONING(double sne_pos[3]);

enum SNE_type sne_random_in_thin_disc_TIMING(void);
void sne_random_in_thin_disc_POSITIONING(double sne_pos[3]);
int sne_random_in_thin_disc_VALIDATE_POSITION(double sne_pos[3], double radius);

enum SNE_type sne_following_a_given_general_distribution_TIMING(void);
void sne_following_a_given_general_distribution_POSITIONING(double sne_pos[3]);

enum SNE_type sne_at_cluster_position_and_random_TIMING(void);

enum SNE_type sne_sink_feedback_TIMING(void);
void sne_sink_feedback_POSITIONING(double sne_pos[3]);
void sne_sink_feedback_RETURN_MASS(int local_n, int indices[], double total_volume, double *total_mass, double *mean_density,
                                   double sn_mass);

enum SNE_type sne_sinks_and_SNIa_TIMING(void);
enum SNE_type sne_sinks_and_SNIa_McMillan_TIMING(void);
void sne_at_disc_particle_position_POSITIONING(double sne_pos[3]);

#ifdef POPIII_SNE
void sne_sink_feedback_PROPERTIES(double *sn_mass, double *sn_E51);
#endif

enum SNE_type sne_poisson_sampling_SNII_AND_SNIa_TIMING(void);
void sne_poisson_sampling_SNII_AND_SNIa_POSITIONING(double sne_pos[3]);

double random_exponential_distributed(double rate);

#ifdef INJECT_TRACER_INTO_SN
int inject_tracer_particles(int local_n, int indices[], int total_n, enum SNE_type sne_type);
void add_new_tracer(int p, MyIDType newid);
#endif

#ifdef TRACER_MC
void remove_remaining_tracers_from_cell(int p);
#endif

void sne_log(double sne_pos[3], int n, int injection_scheme, double radius, double mean_density, double mean_temperature, int n_tracer,
             double sn_mass, double sn_E51);

void sne_init(void);
void sne_destroy(void);

void sne_initialize_radial_CDF(void);
double inverse_radial_CDF(double x);
double inverse_z_CDF(double x);

#ifdef GALPOT
void find_total_stellar_and_gas_mass(void);
#endif
