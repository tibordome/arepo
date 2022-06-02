/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind.h
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

#ifndef SUBFIND_H
#define SUBFIND_H

#include "../allvars.h"
#include "../domain.h"

#define FIND_SMOOTHING_LENGTHS 0
#define FIND_TOTAL_DENSITIES 1

#define SUBFIND_SO_POT_CALCULATION_PARTICLE_NUMBER 10000

#define SUBFIND_GAL_RADIUS_FAC 2.0 /* for subfind metal calculation */

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
extern int *NodeGrNr;
#endif

extern int GrNr;
extern int NumPartGroup;

extern struct topnode_data *SubTopNodes;
extern struct local_topnode_data *Sub_LocTopNodes;

extern int *SubDomainTask;
extern int *SubDomainNodeIndex;
extern int *SubNextnode;
extern int SubNTopleaves;
extern int SubNTopnodes;

extern int SubTree_MaxPart;
extern int SubTree_NumNodes;
extern int SubTree_MaxNodes;
extern int SubTree_FirstNonTopLevelNode;
extern int SubTree_NumPartImported;
extern int SubTree_NumPartExported;
extern int SubTree_ImportedNodeOffset;
extern int SubTree_NextFreeNode;
extern MyDouble *SubTree_Pos_list;
extern struct NODE *SubNodes;
extern struct ExtNODE *SubExtNodes;

extern double SubTreeAllocFactor;

extern int *SubTree_ResultIndexList;
extern int *SubTree_Task_list;
extern unsigned long long *SubTree_IntPos_list;

extern double SubDomainCorner[3], SubDomainCenter[3], SubDomainLen, SubDomainFac;
extern double SubDomainInverseLen, SubDomainBigFac;

extern MyDouble GrCM[3];

extern int Ncollective;
extern int NprocsCollective;
extern int MaxNsubgroups;
extern int MaxNgbs;
extern int MaxSerialGroupLen;
extern r2type *R2list;

extern int CommSplitColor;
extern MPI_Comm SubComm;

extern int SubNTask, SubThisTask;
extern int SubTagOffset;

#ifdef ADD_GROUP_PROPERTIES
extern MyFloat *SubGroupPos, *SubGroupPosLocal, *SubGroupPosSend;
#endif

extern struct proc_assign_data
{
  int GrNr;
  int Len;
  int FirstTask;
  int NTask;
} * ProcAssign;

extern struct subgroup_properties
{
  int Len;
  int LenType[NTYPES];
  int GrNr;
  int SubNr;
  int SubParent;
  MyIDType SubMostBoundID;
  MyFloat Mass;
  MyFloat MassType[NTYPES];
  MyFloat SubVelDisp;
  MyFloat SubVmax;
  MyFloat SubVmaxRad;
  MyFloat SubHalfMassRad;
  MyFloat SubHalfMassRadType[NTYPES];
  MyFloat SubMassInRad;
  MyFloat SubMassInRadType[NTYPES];
  MyFloat SubMassInHalfRad;
  MyFloat SubMassInHalfRadType[NTYPES];
  MyFloat SubMassInMaxRad;
  MyFloat SubMassInMaxRadType[NTYPES];
  MyFloat Pos[3];
  MyFloat CM[3];
  MyFloat Vel[3];
  MyFloat Spin[3];

#ifdef MHD
  MyFloat Bfld_Halo, Bfld_Disk;
#endif

#ifdef SUBFIND_EXTENDED_PROPERTIES
  MyFloat Ekin, Epot, Ethr;
  MyFloat J[3], Jdm[3], Jgas[3], Jstars[3], CMFrac, CMFracType[NTYPES];
  MyFloat J_inRad[3], Jdm_inRad[3], Jgas_inRad[3], Jstars_inRad[3], CMFrac_inRad, CMFracType_inRad[NTYPES];
  MyFloat J_inHalfRad[3], Jdm_inHalfRad[3], Jgas_inHalfRad[3], Jstars_inHalfRad[3], CMFrac_inHalfRad, CMFracType_inHalfRad[NTYPES];
#endif

#ifdef USE_SFR
  MyFloat Sfr, SfrInRad, SfrInHalfRad, SfrInMaxRad, GasMassSfr;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MyFloat GasMassMetallicity;
  MyFloat GasMassMetallicityHalfRad;
  MyFloat GasMassMetallicityMaxRad;
  MyFloat GasMassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetalsHalfRad[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetalsMaxRad[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetallicity;
  MyFloat StellarMassMetallicityHalfRad;
  MyFloat StellarMassMetallicityMaxRad;
  MyFloat StellarMassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetalsHalfRad[GFM_N_CHEM_ELEMENTS];
  MyFloat StellarMassMetalsMaxRad[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetallicitySfr;
  MyFloat GasMassMetalsSfr[GFM_N_CHEM_ELEMENTS];
  MyFloat GasMassMetallicitySfrWeighted;
  MyFloat GasMassMetalsSfrWeighted[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
  MyFloat GasMassDustMetallicity;
  MyFloat GasMassDustMetallicityHalfRad;
  MyFloat GasMassDustMetallicityMaxRad;
  MyFloat GasMassDustMetallicitySfr;
  MyFloat GasMassDustMetallicitySfrWeighted;
#endif
#endif
#ifdef BLACK_HOLES
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
#endif
#ifdef GFM_WINDS
  MyFloat WindMass;
#endif
#ifdef SUBFIND_MEASURE_H2MASS
  MyFloat H2_Mass;
#endif
#ifdef GFM_STELLAR_PHOTOMETRICS
  MyFloat Magnitude_U, Magnitude_B, Magnitude_V, Magnitude_K;
  MyFloat Magnitude_g, Magnitude_r, Magnitude_i, Magnitude_z;
  MyFloat SurfaceBrightnessLimitRad;
  MyFloat SubMassInPhotRad;
#endif
} * SubGroup;

#ifdef ADD_GROUP_PROPERTIES
extern struct subgroup_properties *SubGroupAll;
#endif

#if defined(MHD) && defined(ADD_MAGNETIC_GROUP_PROPERTIES)
struct rlist_mhd
{
  double r;
  double b_egy;
};
#endif

extern struct nearest_r2_data
{
  double dist[2];
} * R2Loc;

extern struct nearest_ngb_data
{
  long long index[2];
  int count;
} * NgbLoc;

extern int NumPaux;

extern struct paux_data
{
  int TaskOfGr;
  int LocGrIndex;
  unsigned char Type;
  unsigned char SofteningType;
  MyDouble Pos[3];
  MyDouble Mass;
} * Paux;

extern struct submp_data
{
  int index;
  int GrNr;
  int OldIndex;
  MyFloat DM_Density;
#ifdef ADD_GROUP_PROPERTIES
  int OriginalSubNr;
#endif
} * submp;

extern struct cand_dat
{
  int head;
  int len;
  int nsub;
  int rank, subnr, parent;
  int bound_length;
} * candidates;

extern struct coll_cand_dat
{
  long long head;
  long long rank;
  int len;
  int nsub;
  int subnr, parent;
  int bound_length;
} * coll_candidates;

typedef struct
{
  double rho;
#ifdef SUBFIND_CALC_MORE
  double vx, vy, vz;
  double v2;
#endif
} SubDMData;

void subfind_determine_sub_halo_properties(struct unbind_data *d, int num, struct subgroup_properties *subgroup, int grnr, int subnr,
                                           int parallel_flag, int nsubgroups_cat, int grnr_local);
int subfind_treefind_collective_export_node_threads(int no, int i, int thread_id);
void subfind_domain_do_local_refine(int n, int *list);
void assign_group_numbers_based_on_catalogue(int ngroups_cat, int nsubgroups_cat);
void subfind_mark_cgm(void);
double subfind_so_potegy(double *egypot);
double subfind_exchange(void);
void subfind_save_final(int num);

void subfind_distlinklist_get_two_heads(long long ngb_index1, long long ngb_index2, long long *head, long long *head_attach);
int subfind_distlinklist_get_tail_set_tail_increaselen(long long index, long long *tail, long long newtail);
long long subfind_distlinklist_setrank_and_get_next(long long index, long long *rank);
long long subfind_distlinklist_get_rank(long long index);
void subfind_distlinklist_set_next(long long index, long long next);
void subfind_distlinklist_add_particle(long long index);
void subfind_distlinklist_add_bound_particles(long long index, int nsub);
void subfind_distlinklist_mark_particle(long long index, int target, int submark);
long long subfind_distlinklist_get_next(long long index);
long long subfind_distlinklist_get_head(long long index);
void subfind_distlinklist_set_headandnext(long long index, long long head, long long next);
void subfind_distlinklist_set_tailandlen(long long index, long long tail, int len);
void subfind_distlinklist_get_tailandlen(long long index, long long *tail, int *len);
void subfind_distlinklist_set_all(long long index, long long head, long long tail, int len, long long next);
long long subfind_distlinklist_set_head_get_next(long long index, long long head);

void subfind_reorder_according_to_submp(void);
void subfind_reorder_PS(int *Id, int Nstart, int N);
void subfind_reorder_P(int *Id, int Nstart, int N);
void subfind_distribute_particles(MPI_Comm Communicator);

void subfind_process_group_collectively(int nsubgroups_cat);
void subfind_coll_findExtent(void);
void subfind_coll_domain_decomposition(void);
void subfind_coll_domain_combine_topleaves_to_domains(int ncpu, int ndomain);
void subfind_coll_domain_allocate(void);
void subfind_coll_domain_free(void);
int subfind_coll_domain_determineTopTree(void);
void subfind_coll_domain_walktoptree(int no);
void subfind_col_find_coll_candidates(int totgrouplen);
void subfind_coll_treeallocate(int maxpart, int maxindex);
void subfind_coll_treefree(void);
void subfind_coll_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z);
void subfind_coll_exchange_topleafdata(void);
void subfind_coll_update_node_recursive(int no, int sib, int father, int *last);
void subfind_coll_insert_pseudo_particles(void);
int subfind_coll_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, unsigned long long xc, unsigned long long yc,
                                    unsigned long long zc, unsigned long long ilen);
int subfind_coll_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char levels);
int subfind_coll_treebuild_construct(int npart, struct unbind_data *ud);
int subfind_coll_treebuild(int npart, struct unbind_data *ud);

int subfind_unbind(struct unbind_data *ud, int len, int *len_non_gas);
void subfind_unbind_independent_ones(int count);
void subfind_distribute_groups(void);
void subfind_potential_compute(int num, struct unbind_data *d, int phase, double weakly_bound_limit);
int subfind_col_unbind(struct unbind_data *d, int num, int *num_non_gas);
void subfind_find_linkngb(void);
int subfind_loctree_treebuild(int npart, struct unbind_data **mp);
void subfind_loctree_update_node_recursive(int no, int sib, int father);
double subfind_loctree_treeevaluate_potential(int target);
void subfind_loctree_copyExtent(void);
double subfind_locngb_treefind(MyDouble xyz[3], int desngb, double hguess);
void subfind_loctree_findExtent(int npart, struct unbind_data *mp);
int subfind_locngb_treefind_variable(MyDouble searchcenter[3], double hguess);
size_t subfind_loctree_treeallocate(int maxnodes, int maxpart);
void subfind_loctree_treefree(void);
void subfind_find_nearesttwo(void);
int subfind_process_group_serial(int gr, int offset, int nsubgroups_cat);
void subfind_poll_for_requests(void);
double subfind_get_particle_balance(void);

int subfind_compare_rlist_mhd(const void *a, const void *b);
int subfind_compare_procassign_GrNr(const void *a, const void *b);
int subfind_compare_FileOrder(const void *a, const void *b);
int subfind_compare_submp_OldIndex(const void *a, const void *b);
int subfind_compare_submp_GrNr_DM_Density(const void *a, const void *b);
int subfind_compare_densities(const void *a, const void *b);
int subfind_compare_binding_energy(const void *a, const void *b);
int subfind_compare_dist_rotcurve(const void *a, const void *b);
int subfind_compare_coll_candidates_rank(const void *a, const void *b);
int subfind_compare_coll_candidates_boundlength(const void *a, const void *b);
int subfind_compare_coll_candidates_nsubs(const void *a, const void *b);
int subfind_compare_coll_candidates_subnr(const void *a, const void *b);
int subfind_locngb_compare_key(const void *a, const void *b);
int subfind_compare_serial_candidates_subnr(const void *a, const void *b);
int subfind_compare_serial_candidates_rank(const void *a, const void *b);
int subfind_compare_serial_candidates_boundlength(const void *a, const void *b);
int subfind_compare_dist_rotcurve(const void *a, const void *b);
int subfind_compare_binding_energy(const void *a, const void *b);
int subfind_compare_ID_list(const void *a, const void *b);
int subfind_compare_SubGroup_GrNr_SubNr(const void *a, const void *b);
int subfind_compare_dist_rotcurve(const void *a, const void *b);

#ifdef SUBFIND_EXTENDED_PROPERTIES
void subfind_fof_calc_am_collective(int snapnr, int ngroups_cat);
int subfind_fof_calc_am_serial(int gr, int Offs, int snapnr, int ngroups_cat);
void subfind_add_grp_props_calc_fof_angular_momentum(int num, int ngroups_cat);
#endif

#ifdef ADD_GROUP_PROPERTIES
int get_number_of_groups_in_catalogue(int snapnr);
int get_number_of_subgroups_in_catalogue(int snapnr);
int get_number_of_group_catalogue_files(int snapnr);
int get_number_of_groups_in_file(int snapnr, int fnr);
int get_number_of_subgroups_in_file(int snapnr, int fnr);
void get_catalogue_prop(int snapnr, int ngroups_cat);
void fof_read_group_prop(int snapnr, int number_of_files, int *group_range);
void fof_read_subgroup_prop(int snapnr, int number_of_files, int *subgroup_range);
void fof_append_group_properties(int snapnr, int ngroups_cat);
void fof_append_subgroup_properties(int snapnr, int ngroups_cat);
void fof_collect_groups(void);
void fof_collect_subgroups(void);
void subfind_add_grp_props_read_catalogue(int num, int ngroups_cat, int nsubgroups_cat);
void subfind_add_grp_props_distribute_catalogue_subfind(void);
void subfind_add_grp_props_finalize(int num, int ngroups_cat, int nsubgroups_cat);
#endif

#endif
