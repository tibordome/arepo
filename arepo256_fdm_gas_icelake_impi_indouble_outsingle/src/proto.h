/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/proto.h
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

#ifndef PROTO_H
#define PROTO_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "amr/amr_proto.h"
#include "coord_util.h"
#include "forcetree.h"
#include "helm_eos.h"
#include "opal_eos.h"
#include "timer.h"

#if NUM_THREADS > 1
#include <omp.h>
#endif

#ifdef BLACK_HOLES
#include "blackhole/blackhole_proto.h"
#endif

#ifdef MONOTONE_CONDUCTION
#include "conduction/conduction.h"
#endif

#ifdef VORONOI_PROJ
#include "voronoi_proj.h"
#endif

#ifdef SIDM
#include "sidm/sidm_proto.h"
#endif

#ifdef DVR_RENDER
#include "dvr_render/dvr_render.h"
#endif

#if defined(COOLING) && !defined(GRACKLE) && !defined(CHIMES)
#include "cooling/cooling_proto.h"
#endif

#ifdef ATOMIC_DM /*_COOLING */
#include "atomic_dm/cooling_atomic_dm_proto.h"
#endif

#ifdef GFM_STELLAR_EVOLUTION
#include "GFM/stellar_evolution_proto.h"
#endif

#ifdef GFM_COOLING_METAL
#include "GFM/cooling_metal_proto.h"
#endif

#ifdef GFM_DUST
#include "GFM/stellar_evolution_dust_proto.h"
#endif

#ifdef MRT
#include "MRT/RT.h"
#include "MRT/RT_proto.h"
#endif

#ifdef LOCAL_FEEDBACK
#include "local_feedback/local_feedback.h"
#endif

#ifdef GFM_DUST_COOLING
#include "GFM/cooling_dust_proto.h"
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
#include "GFM/stellar_photometrics_proto.h"
#endif

#ifdef GFM_AGN_RADIATION
#include "GFM/agn_radiation_proto.h"
#endif

#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL) || defined(GFM_WINDS_VARIABLE)
#include "GFM/winds_proto.h"
#endif

#ifdef GFM
#include "GFM/helper_proto.h"
#endif

#ifdef DUST_LIVE
#include "dust_live/dust_proto.h"
#endif

#ifdef ADJ_BOX_POWERSPEC
#include "power_spec/adj_box_powerspec_proto.h"
#endif

#ifdef FLD
#include "fld/fld_proto.h"
#endif

#ifdef SMUGGLE_SFR
#include "SMUGGLE/sfr_proto.h"
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
#include "SMUGGLE/stellar_feedback_proto.h"
#endif

#ifdef SMUGGLE_MOLEC_COOLING
#include "SMUGGLE/cooling_molecules_proto.h"
#endif

#ifdef SMUGGLE_DUST_HEATING_COOLING
#include "SMUGGLE/heating_cooling_dust_proto.h"
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK
#include "SMUGGLE/radiation_stellar_feedback_proto.h"
#endif

#ifdef TEST_COOLING_METAL
#include "Test_Cooling_Metal/test_cooling_metal_proto.h"
#endif

#ifdef SMUGGLE_TEST_SFR
#include "SMUGGLE/test_sfr_proto.h"
#endif

#ifdef OTVET
#include "OTVET/otvet_proto.h"
#ifdef SMUGGLE_RADPRESS_OPT_THICK
#include "SMUGGLE/radpressthick_proto.h"
#endif
#endif

#ifdef GALPOT
#include "galpot/galpot.h"
#endif

#ifdef DIFFUSION
void diffuse(double *data, double *diffusion_coeff, double dt);
#endif

#ifdef TURBULENT_METALDIFFUSION
void turbulent_metal_mixing(void);
void mm_exchange_in_vector(double *in);
double mm_vector_multiply(double *a, double *b);
void calculate_shear_tensor(void);
void mm_gradient(double *in);
void mm_matrix_multiply(double *in, double *out, double *diag, double theta);
void mm_matrix_multiply_isotropic(double *in, double *out, double *diag, double theta);
#endif

/* VITALI  --- CONDUCTION */
#ifdef CONDUCTION
void conduction(void);
void init_spitzer_conductivity(void);
void conduction_matrix_multiply(double *in, double *out, double *diag, double theta);
double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport, int *nsend_local);
#endif

#ifdef BRAGINSKII_VISCOSITY
void init_braginskii_viscosity(void);
#endif

#ifdef SMUGGLE_RADPRESS_OPT_THIN
#include "SMUGGLE/radpressthin_proto.h"
#endif

#ifdef CIRCUMSTELLAR
#include "circumstellar/circumstellar_proto.h"
#endif

#ifdef TGSET
#include "tgset/tgset_proto.h"
#endif

#ifdef TGCHEM
#include "tgchem/tgchem_proto.h"
#endif

#ifdef HEALRAY
#include "healray/healray_proto.h"
#endif

#ifdef SGCHEM
#include "SGChem/sgchem_proto.h"
#endif

#ifdef SNE_FEEDBACK
#include "sne/sne_proto.h"
#endif

#ifdef SINKS
#include "sinks/sinks_proto.h"
#endif

#ifdef SINK_PARTICLES
#include "sink_particles/proto_sink_particles.h"
#endif

#ifdef TREECOLV2
#include "TreeColV2/treecolv2_proto.h"
#endif

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
#include "shock_finder/shock_finder.h"
#endif

#ifdef DG
#include "dg/dg_proto.h"
#endif

#ifdef MHD_CT
#include "constrained_transport/constrained_transport.h"
#endif

#ifdef DM_WINDTUNNEL
#include "dmwindtunnel/dmwindtunnel.h"
#endif

#ifdef BECDM
#include "becdm/becdm.h"
#endif

#ifdef NON_IDEAL_MHD
#include "mhd_nonideal/proto_non_ideal_mhd.h"
#endif

#ifdef MODGRAV
#include "modgrav/modgrav_forcetree.h"
#include "modgrav/modgrav_pm.h"
#endif

#ifdef SIMPLEX
#include "simplex/sx_proto.h"
#endif

#ifdef SFR_MCS
#include "sfr_mcs/sfr_mcs_proto.h"
#endif

#ifdef SOLAR
#include "solar/solar.h"
#endif

void sfr_init(void);
void sfr_create_star_particles(void);

#ifdef SUBBOX_SNAPSHOTS
void read_subbox_coordinates(const char *fname);
#endif

void interpolate_from_wind_table(double t, double *rho, double *vel);
void read_windtunnel_file(void);
void WindtunnelReadIn_InitialiseGlobals(void);
void WindtunnelReadIn_CalculateIndex(float, float, float, int *, int *, int *, float *, float *, float *);
int WindtunnelReadIn_GetFlatIndex(int, int, int, int, int, int);
float WindtunnelReadIn_TrilinearInterpolation(float *, int, int, int, float, float, float);

void ngb_finish_rangebounds_update(int nchanged, int *nodelist);
void ngb_update_rangebounds(int i, int *nchanged, int *nodelist);
double fof_find_nearest_dmparticle(MyIDType *vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask);

int ngb_treefind_fof_nearest(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport,
                             int *nsend_local);
int ngb_treebuild(int npart);
void ngb_insert_pseudo_particles(void);
void ngb_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z);
void ngb_recompute_nodes(void);
void ngb_recompute_node_recursive(int no, int mode);
void ngb_treeallocate(void);
void ngb_treefree(void);
int ngb_treefind_export_node_threads(int no, int target, int thread_id, int image_flag);
int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int mode, int thread_id, int numnodes,
                                  int *firstnode);
void ngb_recompute_nodes_test(void);
void ngb_recompute_nodes_normal(void);
void ngb_treerealloc(int delta_Nodes);

void drift_node(struct NgbNODE *current, integertime time1);
void drift_all_particles(void);

void blackhole_blow_wind(void);
void blackhole_update_wind_affected_cells(void);

#if defined(PMGRID) || defined(POWERSPEC_GRID)
void my_slab_based_fft(fft_plan *plan, void *data, void *workspace, int forward);
void my_slab_based_fft_c2c(fft_plan *plan, void *data, void *workspace, int forward);

void my_slab_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
void my_slab_transposeA(fft_plan *plan, fft_real *field, fft_real *scratch);
void my_slab_transposeB(fft_plan *plan, fft_real *field, fft_real *scratch);

void my_column_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
void my_column_based_fft_init_c2c(fft_plan *plan, int NgridX, int NgridY, int NgridZ);
void my_column_based_fft(fft_plan *plan, void *data, void *workspace, int forward);
void my_column_based_fft_c2c(fft_plan *plan, void *data, void *workspace, int forward);

void my_fft_swap23(fft_plan *plan, fft_real *data, fft_real *out);
void my_fft_swap13(fft_plan *plan, fft_real *data, fft_real *out);
void my_fft_swap23back(fft_plan *plan, fft_real *data, fft_real *out);
void my_fft_swap13back(fft_plan *plan, fft_real *data, fft_real *out);
#endif

#ifdef MHD
void do_mhd_source_terms_first_half(void);
void do_mhd_source_terms_second_half(void);
void do_mhd_source_terms_divb(void);
#endif

#ifdef MHD_DEDNER
double get_dedner_speed(int i);
void compute_dedner_speed(void);
#endif

#ifdef INSPIRAL
void do_inspiral_source_terms_first_half(void);
void do_inspiral_source_terms_second_half(void);
#endif

double get_desired_softening_from_mass(double mass);
void log_restart_debug(void);
int get_thread_num(void);
void report_pinning(void);
void detect_topology(void);
void pin_to_core_set(void);
void get_core_set(void);
void drift_particle_core(int i, integertime time1);
int derefine_should_this_cell_be_merged(int i, int flag);
int can_this_cell_be_split(int i);
void gravity_external(void);
void gravity(int timebin, int fullflag);
int my_ffsll(peanokey i);
void powerspec_vel(int RestartSnapNum);
void set_cosmo_factors_for_current_time(void);
void sub_turb_move_perturbers(double t0, double t1);
void sub_turb_add_forces(void);
void sub_turb_read_table(void);
void sub_turb_parent_halo_accel(double dx, double dy, double dz, double *acc);
double sub_turb_enclosed_mass(double r, double msub, double vmax, double radvmax, double c);
void calc_exact_gravity_for_particle_type(void);
void calculate_non_standard_physics_with_valid_gravity_tree(void);
void calculate_non_standard_physics_with_valid_gravity_tree_always(void);
void calculate_non_standard_physics_prior_gravity_calculation(void);
int get_softeningtype_for_hydro_cell(int i);

#ifdef ACCRETE_ONTO_CENTRAL_POTENTIAL
void accrete_onto_central_potential(void);
#endif

#ifdef PERTURB_VELOCITIES
void perturb_velocities(void);
#endif

double blackhole_friction_dm_density(void);
void blackhole_friction_apply(void);
void gravity_forcetest_testforcelaw(void);

void powersepc_turb_init(void);
int powerspec_turb_find_nearest_evaluate(int target, int mode, int thread_id);
void powerspec_turb_calc_dispersion(void);
double powerspec_turb_obtain_fields(int type);
void powerspec_turb_save(char *fname, double *disp);
void powerspec_turb_collect(void);
void powerspec_turb(int filenr, int type);
void subdivide_evenly(int N, int pieces, int index, int *first, int *count);
void force_evaluate_direct(int target, int result_idx, int nimport);
void gravity_direct(int timebin);
void gravity_force_finalize(int timebin);
void permutate_chunks_in_list(int ncount, int *list);
double get_default_softening_of_particletype(int type);

void get_disk_forces(double RR, double zz, double *f_R, double *f_z);
void growing_disk_init(void);

int imax(int a, int b);
int imin(int a, int b);
size_t smax(size_t a, size_t b);
double max_array(const double *a, int num_elements);
void minimum_large_ints(int n, const long long *src, long long *res);
void sumup_large_ints_comm(int n, const int *src, long long *res, MPI_Comm comm);
void sumup_large_ints(int n, const int *src, long long *res);
void sumup_longs(int n, const long long *src, long long *res);
double second(void);
double measure_time(void);
double timediff(double t0, double t1);
double mysort(void *base, size_t nel, size_t width, int (*compar)(const void *, const void *));
int myflush(FILE *fstream);
int flush_everything(void);
double get_random_number_aux(void);

double get_softening_of_particle(int i);
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
int get_softening_type_from_mass(double mass);
#endif

void ngb_update_velocities(void);
void hello(void);
void sticky_boundary_vertex_velocities(void);

void sticky_vertex_velocities(void);

double blackhole_get_radiomode_mdot(double vvir);
void blackhole_place_nf_bubbles(int num_bubbles, int actioncode);
void blackhole_do_bubbles_nf(int num_activebh);
double blackhole_get_radio_efficiency(double vvir);

void init_turb(void);
void set_turb_ampl(void);
void add_turb_accel(void);
void reset_turb_temp(void);
void log_turb_temp(void);
void init_static_nfw(void);
double get_turb_pot(double x, double y, double z);
void do_turb_driving_step_first_half(void);
void do_turb_driving_step_second_half(void);

void circumstellar_swallow_gas(void);
void find_long_range_step_constraint(void);

void ngb_treemodifylength(int delta_NgbMaxPart);
void domain_resize_storage(int count_get, int count_get_sph, int option_flag);
void domain_resize_storage_stars(int count_get_star);
void domain_resize_storage_blackholes(int count_get_BHs);
void domain_resize_storage_dust(int count_get_dust);
void domain_resize_storage_tracer(int count_get_tracer);

void init_individual_softenings(void);
double get_softening_from_mass(double mass, int type);

void do_derefinements_and_refinements(void);
int refine_criterion_volume(int i);
int refine_criterion_special_boundary(int i);
int refine_criterion_windtunnel(int i);
int refine_criterion_jeans_ref(int i);
int refine_criterion_default(int i);

#ifdef BH_BASED_CGM_ZOOM
double distance_dependent_target_mass(double r);
void bh_based_cgm_zoom_update(void);
#endif

void fill_slice(ray_data *Ray, int Nray, int gradients_flag, int pixels_x, int pixels_y, float *density, float *temperature,
                float *metallicity, float *velocity, float *Bfield, float *vorticity, float *photon_density, float *chem_elements,
                float *density_trmc
#ifdef CHEM_IMAGE
                ,
                float *dust, float *xH2, float *xHP, float *xCO
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 16
                ,
                float *xCHX, float *xOHX, float *xHCOP, float *xCP, float *xMP, float *xHEP
#endif
#if CHEMISTRYNETWORK == 1 || CHEMISTRYNETWORK == 15
                ,
                float *xHEP
#endif
#endif
#ifdef SX_OUTPUT_IMAGE
                ,
                float *rih, float *hrih
#endif
#ifdef SX_OUTPUT_IMAGE_ALL
                ,
                float *sxrates
#endif
);

double get_disk_mass(double time);
void mark_active_timebins(void);
void set_pinning_openmp_threads(void);
void report_pinning_openmp_threads(void);

void voronoi_test(void);

void execute_resubmit_command(void);
void output_compile_time_options(void);
void init_io_fields(void);
void produce_dump(void);
void create_snapshot_if_desired(void);
void output_log_messages(void);
void mpi_report_committable_memory(void);
long long report_comittable_memory(long long *MemTotal, long long *Committed_AS, long long *SwapTotal, long long *SwapFree);

int check_for_interruption_of_run(void);
void set_non_standard_physics_for_current_time(void);
void calculate_non_standard_physics_prior_mesh_construction(void);
void calculate_non_standard_physics_end_of_step(void);
void agn_radiation_info(void);
void compute_statistics(void);
void dvr_render_main(void);

void rt_init_sourceid(void);

double godunov_flux_3d_hllc_states(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
double godunov_flux_3d_hllc_gamma(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face,
                                  struct fluxes *flux);

void face_limit_fluxes(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R,
                       struct fluxes *flux, double dt, long long *count, long long *count_reduced);
double get_sound_speed(int p);
void set_pressure_of_cell(int i);
void gradient_init(MyFloat *addr, MyFloat *addr_exch, MySingle *addr_grad, int type);
void limit_vel_gradient(double *d, MySingle *grad_vx, MySingle *grad_vy, MySingle *grad_vz, double csnd);

/* rt cg method */
void rt_cgmethod(tessellation *T, double dt);
double radtransfer_vector_multiply(double *a, double *b);
double radtransfer_vector_sum(double *a);
void radtransfer_matrix_multiply(tessellation *T, double dt, double *in, double *out, double *sum);

/* radiative transfer */
double rt_DoCooling(int i, double dtime);
double rt_DoHeating(int i, double dtime);
void rt_calc_column_density(double *source_pos);
void rt_find_short_characteristics_column_densities(double *source_pos);
void rt_try_to_compute_column(int i);
void rt_find_upstream_ray_intersections(double *source_pos);
void rt_exchange_column_densities(void);
void rt_inject_photons_single(tessellation *T, double dt);
void rt_advect_radiation(tessellation *T, double dt);
void rt_get_vectors(void);
void rt_select_new_brightest_sources(void);
int radinflux_compare_index_sourceid(const void *a, const void *b);
int radinflux_compare_task(const void *a, const void *b);
int rt_compare_difflist_dphotons(const void *a, const void *b);
int rt_compare_difflist_sourceid(const void *a, const void *b);
int rt_get_cell(tessellation *T, double x, double y, double z);
double rt_get_advect_tistep(void);
void rt_inject_photons_spread(tessellation *T, double dt);
int rt_source_evaluate(int target, int mode, int *nexport, int *nsend_local, int source_id, double dphotons);
int rt_source_evaluate_sfr(int target, int mode, int *nexport, int *nsend_local, double dphotons);
void rt_create_source_list(void);
int source_list_lum_compare(const void *a, const void *b);
double rt_get_cooling_rate(int i, double utherm);
void rt_update_chemistry(double);
void rt_photoheating(int, double);
void rt_set_simple_inits(void);
int rt_rate_ODEs(double t, const double y[], double f[], void *params);

void peano_hilbert_key_inverse(peanokey key, int bits, peano1D *x, peano1D *y, peano1D *z);

void find_nearest_meshpoint_global(mesh_search_data *searchdata, int n, int hsmlguess, int verbose);

int voronoi_ghost_search_quick(void);
int ngb_treefind_ghost_search_quick(MyDouble searchcenter[3], MyDouble refpos[3], MyFloat hsml, int target, int origin, int *startnode,
                                    int bitflags, int mode, int *nexport, int *nsend_local);
int voronoi_ghost_search_evaluate_quick(int target, int mode, int q, int *nexport, int *nsend_local);

void reorder_DP(void);
void peano_hilbert_order_DP(void);

void do_expansion_decay_second_half(void);
void validate_vertex_velocities(void);
void add_mesh_correction_vector(void);
double get_cell_radius(int i);
#ifdef ACTIVE_CELL_SPIN
double get_cell_radius_from_volume(double Volume);
#endif

int voronoi_get_connected_particles(tessellation *T);
void voronoi_init_connectivity(tessellation *T);
void voronoi_update_connectivity(tessellation *T);
int compare_foreign_connection(const void *a, const void *b);
void voronoi_remove_connection(int i);

int pmforce_is_particle_high_res(int type, MyDouble *pos);
void cooling_only(void);

void report_VmRSS(void);
void tree_based_timesteps_setsoundspeeds(void);
void coffee_overide_velocities(void);
void init_aux_fields(void);
void voronoi_update_ghost_velvertex(void);

void healthtest(void);

int should_this_cell_be_split(int i);
int do_refinements(void);
int do_derefinements(void);
void move_collisionless_particle(int new_i, int old_i);

void dump_memory_table(void);

void report_detailed_memory_usage_of_largest_task(void);

void calc_picture_contribution(tessellation *T, int tt, point *p0, point *p1, double *sigma, double *sigmatemp, double *sigmaweight,
                               int weight_flag, int gradients_flag
#if defined(COOLING) && !defined(GRACKLE)
                               ,
                               double *sigmaszy
#endif
#ifdef GFM_STELLAR_EVOLUTION
                               ,
                               double *sigmametal
#endif
#ifdef GFM_AGN_RADIATION
                               ,
                               double *sigmaagnbol
#endif
#ifdef TRACER_MC
                               ,
                               double *sigmatracernum, double *sigmatrweight
#endif
#ifdef CHEM_IMAGE
                               ,
                               double *sigmadust, double *sigmah2, double *sigmahp, double *sigmaco
#endif
#ifdef MRT
                               ,
                               double *sigmaphoton, double *sigmah1
#endif
#ifdef SX_OUTPUT_IMAGE
                               ,
                               double *sigmarih, double *sigmahrih
#endif
);

void setup_rays(int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax,
                double zmin, double zmax);
void exchange_rays(void);
int advance_rays_for_one_cell(int integrate, int weight_flag, int gradients_flag);
void make_3d_voronoi_projected_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis,
                                     double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int weight_flag);

int image_get_next_tetra(tessellation *T, int tt, point *ppstart, point *ppend, int *nexttetra, point *ppexit, int *previous_tetra);
void calculate_vertex_velocity_divergence(void);

void make_voronoi_image(int num);
void make_voronoi_image_slice(int num);
void make_list_of_active_particles(void);

void find_gravity_timesteps_and_do_gravity_step_first_half(void);
void do_gravity_step_second_half(void);

void voronoi_1D_reorder_gas(void);
int voronoi_1D_compare_key(const void *a, const void *b);
void voronoi_1D_order(void);

void voronoi_probe_intersections(tessellation *T, int tt, point *ppstart, point *ppend, int *count, int *f, int *edge, int *corner,
                                 int *orientations);

void pm2d_init_periodic(void);
void pm2d_init_periodic_allocate(void);
void pm2d_init_periodic_free(void);
void pm2d_force_periodic(int mode);
int pm2d_periodic_compare_sortindex(const void *a, const void *b);
void pm2d_mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));

int timestep_evaluate(int target, int mode, int threadid);

void tree_based_timesteps(void);

void coffee_get_velocity(double x, double y, double *vx, double *vy);

double externaldisk_potential(double r);
double externaldisk_dphidR(double r, void *param);

#ifdef SPIRAL
void galaxy_potential(double sp_pos, double cellrad, double sp_t, double sp_acc);
double spiral_potential_calc(double x, double y, double z, double ti);
#endif

#ifdef GRAVITY_TABLE
void grav_table_init(void);
void grav_table_find_grav_acceleration(double xi, double yi, double zi, double *sp_acc);
double grav_table_interpolation(double x, double x1, double x2, double f_x1, double f_x2);
#endif

#ifdef HOST_MEMORY_REPORTING
void check_maxmemsize_setting(void);
#endif

int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);

int MPI_Sizelimited_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount,
                             MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);

int MPI_hypercube_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcount, int *displs,
                             MPI_Datatype recvtype, MPI_Comm comm);

double get_cooling_luminosity(int i);

int get_image_limits(int argc, char **argv, int RestartFlag, int *pixels_x, int *pixels_y, int *xaxis, int *yaxis, int *zaxis,
                     double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax, int *weight_flag);
void extract_position_from_axes(int i, int j, int xaxis, int yaxis, int zaxis, int pixels_x, int pixels_y, double xmin, double xmax,
                                double ymin, double ymax, double zval, double *p);
void make_3d_voronoi_slice_image(int num, int gradients_flag, int pixels_x, int pixels_y, int xaxis, int yaxis, int zaxis, double xmin,
                                 double xmax, double ymin, double ymax, double zval);
void make_3d_voronoi_grid(int num, int pixels_x, int pixels_y, int pixels_z, double xmin, double xmax, double ymin, double ymax,
                          double zmin, double zmax);

void set_special_noh_boundary_conditions(void);

void conduction(void);
void conduction_matrix_multiply(double *in, double *out, double *diag, double theta);

void conduction_PCG(double *in, double *out);
void conduction_PCG_isotropic(double *in, double *out);
void conduction_PCG_anisotropic(double *in, double *out);

double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport, int *nsend_local);

void limit_gradient_special(double *d, double press, double rho, double min_csnd, double max_csnd, MySingle *dpress, MySingle *drho);

double blackhole_get_mdot_radio_from_radiolum(int n);

double parallel_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
double parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *), MPI_Comm comm);

int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);

void test_mpi_alltoall_performance(void);

void drift_particle(int i, integertime time1);

void put_symbol(char *string, double t0, double t1, char c);
void write_cpu_log(void);

void conduction_matrix_multiply_isotropic(double *in, double *out, double *diag, double theta);
void conduction_matrix_multiply_anisotropic(double *in, double *out, double *diag, double theta);

size_t FreeBytes(void);

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line, int clear_flag,
                        const char *callorigin);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line, int clear_flag,
                                const char *callorigin);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);
void *myfree_query_last_block(void);

void mymalloc_init(void);

void calculate_maxid(void);
#if defined(INJECT_TRACER_INTO_SN) || defined(TRACK_ROTATING_HIGHRES_REGION)
void calculate_max_tracer_mc_id(void);
#endif

void determine_compute_nodes(void);

double INLINE_FUNC hubble_function(double a);
#ifdef DARKENERGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef TIMEDEPGRAV
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif
#ifdef EXTERNALHUBBLE
double INLINE_FUNC hubble_function_external(double a);
#endif
#endif

#endif

void set_DMmass(void);

void blackhole_centering(void);

int ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport,
                          int *nsend_local);

void fof_fof(int num);
void fof_subfind_exchange(MPI_Comm Communicator);
void fof_prepare_output_order(void);

void subfind(int num);
double subfind_density(int mode);
double subfind_overdensity(void);
void subfind_density_hsml_guess(void);

void write_file(const char *fname, int readTask, int lastTask, int subbox_flag);

void distribute_file(int nfiles, int *filenr, int *master, int *last);

int get_values_per_blockelement(enum iofields blocknr);

int get_datatype_in_block(enum iofields blocknr, int mode);
void get_dataset_name(enum iofields blocknr, char *buf);

int blockpresent(enum iofields blocknr, int write);
void fill_write_buffer(void *buffer, enum iofields blocknr, int *pindex, int pc, int type, int subbox_flag);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);

int get_particles_in_block(enum iofields blocknr, int *typelist);

int get_bytes_per_blockelement(enum iofields blocknr, int mode);

void read_file(const char *fname, int filenr, int readTask, int lastTask, int readTypes);

void get_Tab_IO_Label(enum iofields blocknr, char *label);

void long_range_init_regionsize(void);

int find_files(const char *fname);

double get_random_number(void);

int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);

void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_idlist(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));

int density_isactive(int n);

void GetMachNumberCR(struct sph_particle_data *Particle);
void GetMachNumber(struct sph_particle_data *Particle);
void GetShock_DtEnergy(struct sph_particle_data *Particle);

void my_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno);

void reconstruct_timebins(void);

peanokey peano_hilbert_key(peano1D x, peano1D y, peano1D z, int bits);
peanokey peano_and_morton_key(peano1D x, peano1D y, peano1D z, int bits, peanokey *morton);
peanokey morton_key(peano1D x, peano1D y, peano1D z, int bits);

peanokey position_to_peanokey(MyDouble pos[3]);
int peanokey_to_topnode(peanokey key);

void enable_core_dumps_and_fpu_exceptions(void);

void find_next_sync_point(void);

void set_units_sfr(void);

void gravity_forcetest(void);

void allocate_memory(void);
void begrun0(void);
void begrun1(void);
void begrun2(void);
int init(void);
void loadrestart(void);
void reread_params_after_loading_restart(void);
void check_omega(void);
void close_logfiles(void);
void compute_grav_accelerations(int timebin, int fullflag);
void compute_global_quantities_of_system(void);
void cooling_and_starformation(void);
void density(void);
void do_box_wrapping(void);
void domain_Decomposition(void);
void domain_shiftPosition(double *pos, enum domain_displace_mode mode);
void domain_displacePosition(double *pos, enum domain_displace_mode mode);
void domain_displacePositions(enum domain_displace_mode mode);
double enclosed_mass(double R);
void endrun(void);
void energy_statistics(void);
#ifdef BINARYLOG
void binary_statistics(void);
#endif
void ensure_neighbours(void);

void every_timestep_stuff(void);
void ewald_corr(double dx, double dy, double dz, double *fper);

void ewald_force(double x, double y, double z, double force[3]);

int my_fls(int x);

void ewald_init(void);
double ewald_psi(double x, double y, double z);
double ewald_pot_corr(double dx, double dy, double dz);
integertime find_next_outputtime(integertime time);

double get_starformation_rate(int i);
double calc_egyeff(int i, double gasdens, double *ne, double *x, double *tsfr, double *factorEVP);

#ifdef MODIFIED_EOS
double effective_eos(double rho);
double invert_effective_eos(double egyeff_old);
void check_modified_eos_parameters(void);
#endif

void gravity_tree(int timebin);
int init(void);
#ifndef LT_STELLAREVOLUTION
void init_clouds(void);
void integrate_sfr(void);
#else
void init_clouds(int, double, double, double, double *, double *);
void integrate_sfr(double, double, double, double, double);
#endif
#define my_fwrite(data, size, nmemb, stream) my_fwrite_fullinfo(data, size, nmemb, stream, __func__, __FILE__, __LINE__)
#define my_fread(data, size, nmemb, stream) my_fread_fullinfo(data, size, nmemb, stream, __func__, __FILE__, __LINE__)
size_t my_fwrite_fullinfo(const void *data, size_t size, size_t nmemb, FILE *stream, const char *func, const char *file, int line);
size_t my_fread_fullinfo(void *data, size_t size, size_t nmemb, FILE *stream, const char *func, const char *file, int line);
int my_system(const char *command);
void open_logfiles(void);
void write_outputfiles_header(void);
void peano_hilbert_order(void);
void read_ic(const char *fname, int readTypes);
void read_header_attributes(FILE *fd);
void read_header_attributes_in_hdf5(const char *fname);
MyIDType determine_ids_offset(void);
void read_ic_cluster(const char *fname);
void read_ic_cluster_gas(const char *fname);
void read_ic_cluster_wimp(const char *fname);
void read_outputlist(const char *fname);
void read_parameter_file(const char *fname);
void check_parameters(void);
void reorder_gas(int *Id);
void reorder_particles(int *Id);
void restart(int mod);
void run(void);
void savepositions(int num, int subbox_flag);
void mpi_printf(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
void mpi_fprintf(FILE *stream, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
void mpi_printf_each(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
void file_path_sprintf(char *buf, const char *fmt, ...) __attribute__((format(printf, 2, 3)));
FILE *open_file(const char *fname);

void open_image_files(char *id1, int num, FILE **dens, FILE **temp, FILE **met, FILE **vel, FILE **mag, FILE **vort, FILE **phot,
                      FILE **chem, FILE **denstr, FILE **dust, FILE **h2, FILE **hp, FILE **co, FILE **chx, FILE **ohx, FILE **hcop,
                      FILE **cp, FILE **mp, FILE **hep, FILE **rih, FILE **hrih, FILE **sxrates);

void write_image_header(FILE *fd, int nx, int ny, int nz);
void write_image_footer(FILE *fd, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
void set_softenings(void);
void set_units(void);
void setup_smoothinglengths(void);

double get_hydrokick_factor(integertime time0, integertime time1);
double get_gravkick_factor(integertime time0, integertime time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1);

#ifdef USE_SFR
void convert_cell_into_star(int i, double birthtime);
void spawn_star_from_cell(int igas, double birthtime, int istar, MyDouble mass_of_star);
int make_star(int idx, int i, double prob, MyDouble mass_of_star, double *sum_mass_stars);
#endif

#ifdef GFM_WINDS
void make_wind(int idx, int i, double prob, MyDouble mass_of_wind, double v_wind, double u_wind);
#endif

#if defined(COOLING) && !defined(GRACKLE)
void cool_cell(int i);
#endif

#ifdef ATOMIC_DM
void ADM_cooling(void);
void ADM_cool_cell(int i);
#endif

/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included! */
#ifdef DECALPHA_NOPOW
double pow(double, double);
#endif

void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(int mode, const int *typelist);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);
int pmforce_nonperiodic(int grnr);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

void pm_setup_nonperiodic_kernel(void);

#ifdef AUTO_SWAP_ENDIAN_READIC
void swap_Nbyte(char *data, int n, int m);
void swap_header(void);
#endif

#ifdef VARIABLE_GAMMA
double godunov_flux_3d_gamma(struct state *st_L, struct state *st_R, struct state_face *st_face);
#endif

#ifdef SECOND_DERIVATIVES
void init_hessians(void);
#endif

#ifdef RT_ADVECT
void rt_advect_main(void);
void rt_init_gradients(void);
int rt_face_get_normals(tessellation *T, int i);
int rt_check_responsibility_of_this_task(const point *DP, int p1, int p2);
#endif

void init_gradients(void);
void init_scalars(void);

#ifdef LOCALLY_ISOTHERM_DISK
double get_isotherm_disk_sound_speed(int p);
int get_isotherm_disk_flag(int p);
#endif

#ifdef SPECIAL_BOUNDARY
void boundary_overide_velocities(struct particle_data *localP, struct sph_particle_data *localSphP, int p);
double usr_defined_get_damping_weight(int i);
void boundary_get_velocity(double x, double y, double z, double *vx, double *vy, double *vz, double dt);
void boundary_x_motion(double x, double y, double z, double *vx, double *vy, double *vz);
void boundary_y_motion(double x, double y, double z, double *vx, double *vy, double *vz);
void boundary_z_motion(double x, double y, double z, double *vx, double *vy, double *vz);
void boundary_circular_motion(double x, double y, double z, double *vx, double *vy, double *vz);
void boundary_double_concentric_motion(double x, double y, double z, double *vx, double *vy, double *vz);
void boundary_disk_motion(double x, double y, double z, double *vx, double *vy, double *vz, double dt);
void get_boundary_cell_state(struct state *state_inside, struct state *state_outside, int orient);
void outside_state_reflective(struct state *state_inside, struct state *state_outside);
void outside_state_noslip(struct state *state_inside, struct state *state_outside);
void outside_state_nrbc(struct state *state_inside, struct state *state_outside);
void outside_state_pressure(struct state *state_inside, struct state *state_outside);
void outside_state_diode(struct state *state_inside, struct state *state_outside, int orient);
#endif

void print_particle_info(int i);
void print_particle_info_from_ID(const MyIDType ID);
void print_state_info(const struct state *st);
void print_state_face_info(const struct state_face *st);

#ifdef VARIABLE_GAMMA
double godunov_flux_3d_gamma(struct state *st_L, struct state *st_R, struct state_face *st_face);
void get_mach_numbers_gamma(struct state *st_L, struct state *st_R, double Press, double *Mach_L, double *Mach_R);
void sample_solution_3d_gamma(double S, struct state *st_L, struct state *st_R, double Press_star, double Vel_star, double W_L,
                              double W_R, struct state_face *st_face);
int riemann_gamma(struct state *st_L, struct state *st_R, double *Press, double *Vel, double *W_L, double *W_R);
void calcW(double Pstar, double GammaFac, struct state *st, double GammaEmin, double GammaEmax, double *W);
#endif

void face_set_scalar_states_and_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);

void face_turn_momentum_flux(struct fluxes *flux, struct geometry *geom);

#if defined(RIEMANN_HLLC) || defined(RIEMANN_ROSUNOV) || defined(RIEMANN_HLLD) || defined(RIEMANN_HLL) || \
    defined(RIEMANN_HLLC_GAMMA) || defined(SPECIAL_RELATIVITY) || defined(GENERAL_RELATIVITY)
void flux_convert_to_lab_frame(struct state *st_L, struct state *st_R, double *vel_face, struct fluxes *flux,
                               struct state_face *st_face);
#endif

void face_clear_fluxes(struct fluxes *flux);
int face_check_responsibility_of_this_task(tessellation *T, int p1, int p2, struct state *st_L, struct state *st_R);
int face_get_normals(tessellation *T, int i, struct geometry *geom);
int face_get_state(tessellation *T, int p, int i, struct state *st);
void face_boundary_check(point *p, double *velx, double *vely, double *velz);
void face_boundary_check_vertex(tessellation *T, int p, double *velx, double *vely, double *velz);
double face_timestep(struct state *state_L, struct state *state_R, double *hubble_a, double *atime);
void state_convert_to_moving_frame(struct state *st, double *vel_face, double hubble_a, double atime);
void face_do_time_extrapolation(struct state *delta, struct state *st, double atime);
void face_do_spatial_extrapolation(struct state *delta, struct state *st, struct state *st_other);
void face_do_spatial_extrapolation_single_quantity(double *delta, double st, double st_other, MySingle *grad, double *dx, double *r);

void face_add_extrapolations(struct state *st_face, struct state *delta_time, struct state *delta_space, struct fvs_stat *stat);
void face_add_extrapolation(struct state *st_face, struct state *delta, struct fvs_stat *stat);

void face_damp_velocities(struct state *st);
void face_turn_velocities(struct state *st, struct geometry *geom);
void solve_advection(struct state *st_L, struct state *st_R, struct state_face *st_face, struct geometry *geom, double *vel_face);
void face_turnback_velocities(struct state_face *st_face, struct geometry *geom);
void face_get_fluxes(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux, struct geometry *geom,
                     double *vel_face);
void face_add_fluxes_advection(struct state_face *st_face, struct fluxes *flux, struct geometry *geom, double *vel_face);

#ifdef TRACER_PARTICLE
void tracer_particle_assign_cell_properties_and_timestep(void);
#endif

#ifdef TRACER_TRAJECTORY
void tracer_init(void);
void tracer_init_output(double tmass);
void tracer_init_output_configuration(void);
void tracer_write_output_if_needed(void);
void tracer_write_output(void);
#endif

#if defined(TRACER_PART_NUM_FLUID_QUANTITIES) || defined(TRACER_MC_NUM_FLUID_QUANTITIES)
void set_tracer_part_indices(void);
#endif

#ifdef TRACER_MC
long long get_total_number_of_tracers(int flag);
int get_max_number_of_tracers(void);
int get_number_of_tracers(int i);
void move_one_tracer(int p, int pother_task, int pother_index, MyIDType pother_ID);
void release_tracer_slot(int itr);

void test_tracer_mc_id_uniqueness(void);
void record_tracer_parent_fluid_properties(void);
void reorder_tracers(void);
void restore_tracer_connectivity(void);
void check_tracer_lists(void);

int consider_moving_tracers(int p, int pother_task, int pother_index, MyIDType pother_ID, double prob);
int consider_moving_tracers_local(int p_from, int p_to, double prob);

int get_free_tracer_slot(void);

void move_tracer_between_parents(int p_from, int p_to, int itracer);
void remove_tracer_from_parent(int p, int itracer);
void add_tracer_to_parent(int p, int itracer);

int exchange_tracer_flux_list(void);
void start_MC_tracer(int nexport_tracer);
void finish_MC_tracer(void);
void load_MC_tracer_start(void);
void load_MC_tracer_reattach(void);
void domain_attach_new_tracers(int offset, int *head_counter);
#endif

#if defined(TRACER_PARTICLE) || defined(TRACER_MC)
void reset_tracer_parent_fluid_properties(void);
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
void special_particle_create_list(void);
void special_particle_update_list(void);
#endif

#ifdef REFINEMENT_AROUND_DM
void dm_particle_create_list(void);
void dm_particle_update_list(void);
#endif

#ifdef GMC_REFINEMENT
int gmc_refinement_criteria(int i);
int gmc_derefinement_criteria(int i);
#endif

#ifdef ONEDIMS_SPHERICAL
void gravity_monopole_1d_spherical(void);
#endif

#ifdef RELAXOBJECT
void relaxobject(void);
#endif

#ifdef RELAXOBJECT_COOLING2
void load_temperature_profil(void);
#endif

#ifdef RELAXOBJECT_BINARY
void do_binary_source_terms_first_half(void);
void do_binary_source_terms_second_half(void);
#endif

int analytic_riemann_solution(int argc, char *argv[]);
double godunov_flux_3d(struct state *st_L, struct state *st_R, struct state_face *st_face);
#if defined(RIEMANN_ROSUNOV) || (defined(RIEMANN_HLLC) && defined(DG))
double godunov_flux_3d_rosunov(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
#endif
#if defined(RIEMANN_HLL) || (defined(RIEMANN_HLLC) && defined(DG))
double godunov_flux_3d_hll(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
#endif
#ifdef RIEMANN_HLLC
double godunov_flux_3d_hllc(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux);
#endif
#ifdef RIEMANN_HLLD
double godunov_flux_3d_hlld(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face, struct fluxes *flux);
#endif
void sample_solution_vaccum_left_3d(double S, struct state *st_R, struct state_face *st_face);
void sample_solution_vaccum_right_3d(double S, struct state *st_L, struct state_face *st_face);
void sample_solution_vaccum_generate_3d(double S, struct state *st_L, struct state *st_R, struct state_face *st_face);
void get_mach_numbers(struct state *st_L, struct state *st_R, double Press);
void sample_solution_3d(double S, struct state *st_L, struct state *st_R, double Press, double Vel, struct state_face *st_face);
int riemann(struct state *st_L, struct state *st_R, double *Press, double *Vel);
void pressure_function(double P, struct state *st, double *F, double *FD);
double guess_for_pressure(struct state *st_L, struct state *st_R);
void riemann_isotherm(struct state *st_L, struct state *st_R, double *Rho, double *Vel, double csnd);
void isothermal_function(double rhostar, double rho, double *F, double *FD);
void sample_solution_isothermal3d(double S, struct state *st_L, struct state *st_R, double Rho, double Vel, struct state_face *st_face,
                                  double csnd);
void apply_flux_list(void);
int flux_list_data_compare(const void *a, const void *b);

void set_vertex_velocities(void);
int scalar_init(MyFloat *addr, MyFloat *addr_mass, int type);
void compute_interface_fluxes(tessellation *T);
void update_primitive_variables(void);
void set_pressure_of_cell_internal(struct particle_data *P, struct sph_particle_data *SphP, int i);
void do_validity_checks(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);
void do_special_boundaries(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);
void update_primitive_variables_single(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);
void update_internal_energy(struct particle_data *P, struct sph_particle_data *SphP, int i, struct pv_update_data *pvd);

#ifdef SPECIAL_RELATIVITY
void update_primitive_variables_special_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i,
                                                   struct pv_update_data *pvd);
double get_cfl_sound_speed_special_relativity(int p);
double godunov_flux_3d_hlle_special_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face,
                                               struct fluxes *flux);
void compute_conserved_quantities_from_ICs_special_relativity(int i);
#endif

#ifdef GENERAL_RELATIVITY
void update_primitive_variables_general_relativity(struct particle_data *P, struct sph_particle_data *SphP, int i,
                                                   struct pv_update_data *pvd);
double get_cfl_sound_speed_general_relativity(int p);
double godunov_flux_3d_hlle_general_relativity(struct state *st_L, struct state *st_R, double *vel_face, struct state_face *st_face,
                                               struct fluxes *flux, struct geometry *geom, double facex, double facey, double facez);
void compute_conserved_quantities_from_ICs_general_relativity(int i);
void get_metric_general_relativity(double x, double y, double z, double *alp, double *psi, double *bx, double *by, double *bz);
void get_metric_and_derivs_general_relativity(double x, double y, double z, double *alp, double *psi, double *bx, double *by,
                                              double *bz, double dalp[3], double dpsi[3], double dbx[3], double dby[3], double dbz[3]);
void do_general_relativity_source_terms(void);
void apply_general_relativity_source_terms(void);
void general_relativity_statistics(void);
#ifdef ATMOSPHERE_GENERAL_RELATIVITY
void general_relativity_determine_maximum_density(void);
#endif
#if METRIC_TYPE == 3 || METRIC_TYPE == 4
int read_fixed_numerical_1d_metric(void);
#endif
#endif

/* MPI utility functions */
void mpi_exchange_buffers(void *send_buf, int *send_count, int *send_offset, void *recv_buf, int *recv_count, int *recv_offset,
                          int item_size, int commtag, int include_self);

int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset, int send_identical);

void *sort_based_on_mesh_search(mesh_search_data *search, void *data, int n_items, int item_size);

void *sort_based_on_field(void *data, int field_offset, int n_items, int item_size);

void mpi_distribute_items_from_search(mesh_search_data *search, void *data, int *n_items, int *max_n, int item_size, int commtag,
                                      int task_offset, int cell_offset);

void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size, int commtag);

void tile_ics(void);

void blackhole_harmonic_force(void);
void reallocate_memory_maxpart(void);
void reallocate_memory_maxpartsph(void);
void reallocate_memory_maxpartBHs(void);
void reallocate_memory_maxpartsinks(void);
void reallocate_memory_maxpartstar(void);
void reallocate_memory_maxpartdust(void);
void reallocate_memory_maxpart_ignore_timebins(void);
void reallocate_memory_maxpartsph_ignore_timebins(void);
void share_particle_number_in_file(const char *fname, int filenr, int readTask, int lastTask, int readTypes);
void reallocate_memory_maxparttracer(int delta);

void calc_memory_checksum(void *base, size_t bytes);

int compare_nodelist_task_index_node(const void *a, const void *b);
int compare_partlist_task_index(const void *a, const void *b);
void report_pinning_threads(void);

void allreduce_sparse_double_sum(const double *loc, double *glob, int N);
void allreduce_sparse_imin(const int *loc, int *glob, int N);

void myMPI_Alltoallv(void *sendb, size_t *sendcounts, size_t *sdispls, void *recvb, size_t *recvcounts, size_t *rdispls, int len,
                     int big_flag, MPI_Comm comm);

int myMPI_Sendrecv(void *sendb, size_t sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvb, size_t recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);
size_t roundup_to_multiple_of_cacheline_size(size_t n);

#ifdef HAVE_HDF5
#include <hdf5.h>

#include "hdf5_util.h"

herr_t my_hdf5_error_handler(void *unused);
void write_header_attributes_in_hdf5(hid_t handle);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_compile_time_options_in_hdf5(hid_t handle);
#endif

void init_cpu_log(void);

void cuda_init(void);
void cuda_finish(void);

void twopoint(void);
void twopoint_save(void);

void calculate_power_spectra(int num, const long long *ntot_type_all);
void calculate_power_spectra_and_ntot(int num, const int *typeflag);

void write_error(int check, size_t nwritten, size_t nmemb);

void init_field(enum iofields field, const char *label, const char *datasetname, enum types_in_memory type_in_memory,
                enum types_in_file type_in_file_output, enum types_in_file type_in_file_input, int values_per_block, enum arrays array,
                void *pointer_to_field, void (*io_func)(int, int, void *, int), int typelist_bitmask);

void init_units(enum iofields field, double a, double h, double L, double M, double V, double c);
void init_snapshot_type(enum iofields field, enum sn_type type);

#ifdef SIDM
void reset_timestep_counter(void);
#endif

#ifdef BAROTROPIC
void apply_barotropic_eos(void);
#endif

#ifdef CHIMES
void ChimesInitCool(void);
void chimes_cooling_only(void);
void chimes_update_gas_vars(int i);
void chimes_update_particle_energies(int i);
double chimes_convert_ucgs_to_temp(double ucgs, int i);
double chimes_convert_temp_to_ucgs(double temp, int i);
double ChimesGetCoolingTime(int i, struct All_rate_variables_structure *this_all_rates);
void chimes_cooling_for_starformation(void);
double chimes_calc_egyeff(int i, double gasdens, double *x, double *tsfr, double *factorEVP);
void init_chimes_parallel(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates,
                          struct Reactions_Structure **this_all_reactions_root,
                          struct Reactions_Structure **this_nonmolecular_reactions_root, double *dustG_arr, double *H2_dissocJ_arr);
void init_chimes_pthreads(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates,
                          struct Reactions_Structure **this_all_reactions_root,
                          struct Reactions_Structure **this_nonmolecular_reactions_root);
void chimes_set_pointers(int i);
void chimes_move_abundances(int i, int j);
void chimes_set_pointers_to_null(int i);
void chimes_update_all_pointers(void);

#ifdef CHIMES_PTHREADS
void chimes_cooling_only_pthreads(void);
void chimes_cooling_for_starformation_pthreads(void);
void chimes_create_pthreads_buffers(int N_active);
void chimes_send_gasvars_to_node_root(int N_active, int *N_active_node, int *buf_offset);
void chimes_run_threads(void);
void *chimes_worker_thread_routine(void *arg);
void chimes_send_abundances_to_original_task(int N_active, int *N_active_node, int *buf_offset);
void chimes_free_pthreads_buffers(int N_active);
#endif

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
void ChimesLoadPhotoIonTables(void);
void chimes_allocate_memory_to_photoion_tables(struct globalVariables *myGlobalVars, struct PhotoIonTables_UVB **myPhotoIonTable,
                                               int N_Elements_in_Bens_tables);
void ChimesReadPhotoIonTables_UVB(struct globalVariables *myGlobalVars, char *photoIonTablePath,
                                  struct PhotoIonTables_UVB *myPhotoIonTable, int N_Elements_in_Bens_tables, int current_spectrum);
void ChimesReadPhotoIonTables_UVB_parallel(struct globalVariables *myGlobalVars, char *photoIonTablePath,
                                           struct PhotoIonTables_UVB *myPhotoIonTable, int N_Elements_in_Bens_tables,
                                           int current_spectrum);
void ChimesCopyPhotoIonTables(struct PhotoIonTables_UVB *myPhotoIonTable, int i, int j);
void ChimesInterpolatePhotoIonTables(struct globalVariables *myGlobalVars, struct PhotoIonTables_UVB *myPhotoIonTable, int z_index,
                                     double dz);
void chimes_redshift_index(double *table, int ntable, double x, int *i, double *dx);
#endif

#ifdef GFM_STELLAR_EVOLUTION
void chimes_update_element_abundances(int i);
#endif
#endif /* CHIMES */

#if defined(BIERMANN_BATTERY) || defined(DURRIVE_BATTERY)
void set_electron_variables_of_cells_for_magnetic_batteries(void);
void do_magnetic_batteries_source_term(void);
#endif
#ifdef DURRIVE_BATTERY
MyDouble f_mt_HI(double e);
MyDouble f_mt_HeI(double e);
MyDouble f_mt_HeII(double e);
#endif

#endif
