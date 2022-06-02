/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/allvars.c
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

#include "allvars.h"

struct data_nodelist *DataNodeList; /* to be deleted */

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
MPI_Status mpistat;
#endif

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

#if NUM_THREADS > 1
omp_lock_t *ParticleLocks;
omp_lock_t *Ngb_NodeLocks;
#endif

int ThisTask; /*!< the number of the local processor  */
int NTask;    /*!< number of processors */
int PTask;    /*!< note: "NTask = 2^PTask" (PTask is the smallest number such that
                         NTask <= 2^PTask) */

int ThisNode;        /**< the rank of the current compute node  */
int NumNodes;        /**< the number of compute nodes used  */
int MinTasksPerNode; /**< the minimum number of MPI tasks that is found on any of the nodes  */
int MaxTasksPerNode; /**< the maximum number of MPI tasks that is found on any of the nodes  */
int TasksInThisNode; /**< number of MPI tasks on  current compute node */
int RankInThisNode;  /**< rank of the MPI task on the current compute node */
long long MemoryOnNode;

double CPUThisRun; /*!< Sums CPU time of current process */

int MaxTopNodes; /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag; /*!< taken from command line used to start code. 0 is normal start-up from
                    initial conditions, 1 is resuming a run from a set of restart files, while 2
                    marks a restart from a snapshot file. Values larger than 2 control
                    special postprocessing options. */
int RestartSnapNum;

int Argc;
char **Argv;

size_t AllocatedBytes;

int Nforces;
int *TargetList;
struct thread_data Thread[NUM_THREADS];

#ifdef IMPOSE_PINNING
hwloc_cpuset_t cpuset_thread[NUM_THREADS];
#endif

int *Exportflag,
    *ThreadsExportflag[NUM_THREADS]; /*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset;
int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;

int *TasksThatSend, *TasksThatRecv, NSendTasks, NRecvTasks;

struct send_recv_counts *Send, *Recv;

int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;
int Force_nimport, Force_nexport, *Force_Send_offset, *Force_Send_count, *Force_Recv_count, *Force_Recv_offset;

int TakeLevel;

int TagOffset;

int TimeBinSynchronized[TIMEBINS];
struct TimeBinData TimeBinsHydro, TimeBinsGravity;

#ifdef TRACER_PARTICLE
struct TimeBinData TimeBinsTracer;
#endif

#ifdef USE_SFR
double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
struct TimeBinData TimeBinsBHAccretion;
double TimeBin_BH_mass[TIMEBINS];
double TimeBin_BH_dynamicalmass[TIMEBINS];
double TimeBin_BH_Mdot[TIMEBINS];
double TimeBin_BH_Medd[TIMEBINS];
#endif

#ifdef SINKS
struct TimeBinData TimeBinsSinksAccretion;
#endif

#ifdef DUST_LIVE
struct TimeBinData TimeBinsDust;
#endif

char DumpFlag              = DUMP_BOTH;
char DumpFlagNextSnap      = DUMP_BOTH;
char WroteSnapThisTimestep = 0;

int FlagNyt = 0;

double CPU_Step[CPU_LAST];
double CPU_Step_Stored[CPU_LAST];

double WallclockTime; /*!< This holds the last wallclock time measurement for timings measurements */
double StartOfRun;    /*!< This stores the time of the start of the run for evaluating the elapsed time */

double EgyInjection;

int NumPart; /*!< number of particles on the LOCAL processor */
int NumGas;  /*!< number of gas particles on the LOCAL processor  */
#ifdef TRACER_MC
int N_tracer; /*!< number of tracer particles on the LOCAL processor  */
#endif
#if defined(GFM) || defined(SFR_MCS)
int N_star; /*!< number of star particles on the LOCAL processor  */
#endif
#ifdef BLACK_HOLES
int NumBHs; /*!< number of BH particles on the LOCAL processor  */
#endif
#ifdef SINKS
int NumSinks; /*!< number of sink particles on the LOCAL processor  */
#endif
#ifdef DUST_LIVE
int N_dust; /*!< number of dust particles on the LOCAL processor */
#endif

#if defined(INJECT_TRACER_INTO_SN) || defined(TRACK_ROTATING_HIGHRES_REGION)
MyIDType MaxTracerID; /*!< stores the maximum ID the tracer reached on any processor */
#endif

gsl_rng *random_generator;     /**< a random number generator  */
gsl_rng *random_generator_aux; /**< an auxialiary random number generator for use if one doesn't want to influence the main code's
                                  random numbers  */

#ifdef SNE_FEEDBACK

gsl_rng *sne_rng;
double SNERandomNextTime = 0.;
#ifdef CLUSTERED_SNE
double SNEClusterNextTime = 0.;
#endif

#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
int Sink_about_to_explode = 0;
int Star_about_to_explode = 0;
int Wherewasi             = -1;
#endif

#endif /*SNE_FEEDBACK*/

#ifdef USE_SFR
int Stars_converted; /*!< current number of star particles in gas particle block */
#endif

#ifdef TOLERATE_WRITE_ERROR
int WriteErrorFlag;
char AlternativeOutputDir[MAXLEN_PATH];
#endif

#ifdef GFM_AGN_RADIATION
int CellsWithAGNRadiation;
#endif

#if defined(GFM_WINDS_LOCAL) || defined(GFM_WINDS)
double WindEnergy_Should, WindEnergy_Is;
#endif
#ifdef BLACK_HOLES
double AGNEnergyEM_Is, AGNEnergyEMobs_Is, AGNEnergyM_Should, AGNEnergyM_Is, AGNEnergyT_Should, AGNEnergyT_Is;
#endif
#ifdef SUBBOX_SNAPSHOTS
int Nsubboxes; /**< the number of subboxes for frequent snapshots */
double *SubboxXmin, *SubboxXmax, *SubboxYmin, *SubboxYmax, *SubboxZmin, *SubboxZmax; /**< the coordinates of subboxes */
#endif

#ifdef EOS_OPAL
double opal_rhomax, opal_rhomin; /**< boundaries of eos table */
#endif

#if defined(SN_MCS) && defined(IMF_SAMPLING_MCS)
int NumSNLocal;  /** Number of SN events this timestep on this processor **/
int NumSNGlobal; /** Number of SN events this timestep everywhere **/
#endif

double TimeOfLastDomainConstruction; /*!< holds what it says */

int *Ngblist; /*!< Buffer to hold indices of neighbours retrieved by the neighbour search
                 routines */
#ifdef GRACKLE
code_units my_grackle_units;
#endif

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
double DomainInverseLen, DomainBigFac;
int *DomainStartList, *DomainEndList;
double *DomainCost, *TaskCost;
int *DomainCount, *TaskCount;
struct no_list_data *ListNoData;

int domain_bintolevel[TIMEBINS];
int domain_refbin[TIMEBINS];
int domain_grav_weight[TIMEBINS];
int domain_hydro_weight[TIMEBINS];
int domain_to_be_balanced[TIMEBINS];

int *DomainTask;
int *DomainNewTask;
int *DomainNodeIndex;

peanokey *Key, *KeySorted;

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;

/* variables for input/output , usually only used on process 0 */

char ParameterFile[MAXLEN_PATH]; /*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,                      /*!< file handle for info.txt log-file. */
    *FdEnergy,                     /*!< file handle for energy.txt log-file. */
    *FdTimings,                    /*!< file handle for timings.txt log-file. */
    *FdDomain,                     /*!< file handle for domain.txt log-file. */
    *FdBalance,                    /*!< file handle for balance.txt log-file. */
    *FdMemory, *FdTimebin, *FdCPU; /*!< file handle for cpu.txt log-file. */

#ifdef SIMPLEX
FILE *FdSimplex; /**< file handle for simplex.txt log-file. */
#endif

#ifdef LOCAL_FEEDBACK
FILE *FdLocalFeedback;
#endif

#ifdef DETAILEDTIMINGS
FILE *FdDetailed;
#endif

#ifdef OUTPUT_CPU_CSV
FILE *FdCPUCSV;
#endif

#ifdef RESTART_DEBUG
FILE *FdRestartTest;
#endif

#ifdef USE_SFR
FILE *FdSfr; /*!< file handle for sfr.txt log-file. */
#endif

#ifdef BLACKHOLE_POTMIN_DIAGNOSTIC
FILE *FdBHDiag;
#endif

#if defined(VS_TURB) || defined(AB_TURB)
FILE *FdTurb;
#endif

#ifdef RT_ADVECT
FILE *FdRad; /*!< file handle for radtransfer.txt log-file. */
#endif

struct pair_data *Pairlist;

#ifdef BLACK_HOLES
FILE *FdBlackHoles; /*!< file handle for blackholes.txt log-file. */
FILE *FdBlackHolesDetails;
FILE *FdBlackHolesMergers;
#ifdef BH_SPIN_EVOLUTION
FILE *FdBlackHolesSpin;
#endif
#ifdef BH_NEW_CENTERING
FILE *FdBlackHolesRepos;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
FILE *FdBlackHolesBipolar;
#endif
#endif

#ifdef SINKS
#ifdef SINKS_MERGERS
FILE *FdSinksMergers;
#endif
#endif

#ifdef GFM_STELLAR_EVOLUTION
FILE *FdMetalsGas;   /*!< file handle for metals_gas.txt log-file. */
FILE *FdMetalsStars; /*!< file handle for metals_star.txt log-file. */
FILE *FdMetalsTot;   /*!< file handle for metals_tot.txt log-file. */
FILE *FdSN;          /*!< file handle for SN.txt log-file. */
#endif

#if defined(SMUGGLE_STAR_FEEDBACK) && defined(SMUGGLE_OUTPUT_STELLAR_FEEDBACK)
FILE *FdFeedback; /*!< file handle for stellar_feedback.txt log-file. */
#endif

#ifdef FORCETEST
FILE *FdForceTest; /*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
FILE *FdDE; /*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef OTVET
FILE *FdOTVET;
FILE *FdOTVETStar;
#endif

#ifdef MRT
FILE *FdMRT;
FILE *FdMRTStar;
#endif

#ifdef BINARYLOG
FILE *FdBinary; /**< file handle for binary.txt log-file. */
#endif

#ifdef NUCLEAR_NETWORK
FILE *FdNetwork;
#endif

#ifdef SINK_PARTICLES
FILE *FdSinkPart;
#endif

#ifdef GENERAL_RELATIVITY
FILE *FdGR;
#endif

#ifdef SNE_FEEDBACK
FILE *FdSNe;
#endif

#ifdef COSMIC_RAYS
FILE *FdCREnergy;
#endif

#ifdef GW_SIGNAL
FILE *FdGW;
#endif

#ifdef SGS_TURBULENCE
FILE *FdSgsTurbulence;
#ifdef SGS_TURBULENCE_LOG_PRODUCTION_DISSIPATION
FILE *FdSgsTurbulenceProductionDissipation;
#endif
#endif

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
FILE *FdDust;
#endif

#ifdef SN_MCS
FILE *FdSnr; /*!< file handle for snr.txt log-file. */
#endif

#ifdef HII_MCS_LOG
FILE *FdHii; /*!< file handle for hii.txt log-file. */
#endif

#ifdef SFR_MCS_LOG
FILE *FdSFdens;
#endif
#ifdef SN_MCS_LOG
FILE *FdSNdens;
#endif

#ifdef SFR_MCS_LOG_DETAILS
FILE *FdSFDetails;
#endif
#ifdef SN_MCS_LOG_DETAILS
FILE *FdSNDetails;
#endif

void *CommBuffer; /*!< points to communication buffer, which is used at a few places */

#ifdef TRACER_PARTICLE
int TracerPartTmaxIndex;       /* store maximum temperature for velocity tracers */
int TracerPartTmaxTimeIndex;   /* store time of maximum temperature for velocity tracers */
int TracerPartTmaxRhoIndex;    /* store density at time of maximum temperature for velocity tracers */
int TracerPartRhomaxIndex;     /* store maximum density for velocity tracers */
int TracerPartRhomaxTimeIndex; /* store time of maximum density for velocity tracers */
int TracerPartMachmaxIndex;    /* store maximum Mach from the riemann solution in SphP for the velocity tracers */
int TracerPartEntmaxIndex;     /* store maximum entropy (= P/rho^gamma) for velocity tracers */
int TracerPartEntmaxTimeIndex; /* store maximum entropy (= P/rho^gamma) for velocity tracers */
#endif

#ifdef TRACER_MC
int TracerMCTmaxIndex;                  /* store maximum temperature for MC tracers */
int TracerMCTmaxTimeIndex;              /* store time of maximum temperature for MC tracers */
int TracerMCTmaxRhoIndex;               /* store density at time of maximum temperature for MC tracers */
int TracerMCRhomaxIndex;                /* store maximum density for MC tracers */
int TracerMCRhomaxTimeIndex;            /* store time of maximum density for MC tracers */
int TracerMCMachmaxIndex;               /* store maximum Mach from the riemann solution in SphP for the MC tracers */
int TracerMCEntmaxIndex;                /* store maximum entropy (= P/rho^gamma) for MC tracers */
int TracerMCEntmaxTimeIndex;            /* store time of maximum entropy (= P/rho^gamma) for MC tracers */
int TracerMCLastStarTimeIndex;          /* store the last time a MC tracer belonged to a star/wind particle */
int TracerMCWindCounterIndex;           /* store the number of times a tracer has been kicked to the wind (GFM_WINDS) */
int TracerMCExchangeCounterIndex;       /* store the number of times a tracer has been moved between particles/cells */
int TracerMCExchangeDistanceIndex;      /* store the sum of the radii of the cells between which the tracer was exchanged */
int TracerMCExchangeDistanceErrorIndex; /* store an estimate of the error of the trajectory of the tracer */
int TracerMCShockMachMaxIndex;          /* store maximum Mach number encountered from the on-the-fly shock finder */
#endif

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
struct global_data_all_processes All;

/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P, /*!< holds particle data on local processor */
    *DomainPartBuf;      /*!< buffer for particle data used in domain decomposition */

struct subfind_data *PS;

/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data *SphP, /*!< holds SPH particle data on local processor */
    *DomainSphBuf;              /*!< buffer for SPH particle data in domain decomposition */

#ifdef TRACER_MC
struct tracer_linked_list *TracerLinkedList;
int *TracerLinkedListHeap;
MyIDType *tracer_cellids;
#endif

#ifdef GRAVITY_TABLE
struct grav_table_data *GravT;
#endif

#if defined(GFM) || defined(SFR_MCS)
struct star_particle_data *StarP, *DomainStarBuf;
#endif
#ifdef BLACK_HOLES
struct bh_particle_data *BHP, *DomainBHBuf;
#endif
#ifdef DUST_LIVE
struct dust_particle_data *DustP, *DomainDustBuf;
#endif
#ifdef SINKS
struct sink_particle_data *SinkP, *DomainSinksBuf;
#endif

#if defined(OTVET) && !defined(GFM)
struct otvet_source_data *OtvetSourceP;
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
struct special_particle_data *PartSpecialListGlobal;
#endif

#ifdef REFINEMENT_AROUND_DM
struct refine_dm_data *DMPartListGlobal;
#endif

#if defined(CIRCUMSTELLAR) &&                                                                                 \
    (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && \
    !defined(EXTERNALGRAVITY)
struct source_particle_data *SourcePartListGlobal;
#endif

#ifdef TRACER_MC
/* tracer communication structure for export/import */
struct tracer_flux_list_data *TracerFluxListIn, *TracerFluxListGet;
int TracerFluxListGet_index;
int Ntracerflux;
#endif

#ifdef TGSET
struct TGD_struct TGD;
#endif

#ifdef TGCHEM
struct TGCD_struct TGCD;
#endif

#ifdef HEALRAY
struct HRD_struct HRD;
struct HRSL_struct *HRSL;
#endif

#ifdef TREECOLV2
int TCV2_nang_theta, TCV2_nang_phi;
int TCV2_npos, TCV2_nang, TCV2_nentries;
double TCV2_min_nodesize, TCV2_max_nodesize;
double TCV2_dnodesize, TCV2_dtheta, TCV2_dphi;
int *TCV2_map_theta_phi_to_ipix;
int *TCV2_start_ind;
int *TCV2_npix;
int *TCV2_pixel_list;
double *TCV2_pixel_values;
#endif

#ifdef SINKS
struct SKD_struct SKD;
#endif

#ifdef SINK_PARTICLES
int NSinksThisTask;
int NSinksAllTasks;
int NSinksAllTasksOld;
int NSinkBufferSize;
int LastTaskNewSink;
int SinksFormedSinceLastDomain;
double SinkCreationDensityCodeUnits;
double SinkAccretionRadius;
double SinkFormationRadiusSquared;
double SinkFormationRadius;
double SinkSofteningLength;
double SinkTargetAccretionMass;
double SinkEvolutionDumpRateYears;
double SinksLastEvolutionDumpTime;

double SinkCreationDensityCurrent;
double SinkFormationRadiusCurrent;
double SinkTestsGiveupDensityCurrent;
#ifdef STORE_SINK_PARTICLE_SPIN
double AngularMomentum[3];
#endif
struct global_sink_particle_data *SinkP, /* The sink particle data (i.e. all sinks)*/
    *export_SinkP;                       /* The communication buffer */
#endif

#ifdef RT_ADVECT
double Source_Pos[N_SOURCES][3];
double Source_Lum[N_SOURCES];
int Source_ID[N_SOURCES];
double rt_vec[RT_N_DIR][3];
double rt_vec_new[RT_N_DIR][3];
#endif

#ifdef MRT
#ifndef MRT_NO_UV

#ifdef MRT_CHEM_SG
double mrt_sigma_PE[UV_BINS]; /*Photo electric heating cross section*/
double G_PE[UV_BINS];         /*Photo electric heating rate*/
double P_PE[UV_BINS];         /*Photo electric heating - momentum injection*/

double mrt_sigma_H2[UV_BINS]; /*H2 cross section*/
double G_H2[UV_BINS];         /*H2 heating*/
double P_H2[UV_BINS];         /*H2 momentum injection*/
#endif

double mrt_sigma_HI[UV_BINS];
double mrt_sigma_HeI[UV_BINS];
double mrt_sigma_HeII[UV_BINS];
double G_HI[UV_BINS];
double G_HeI[UV_BINS];
double G_HeII[UV_BINS];
double P_HI[UV_BINS];
double P_HeI[UV_BINS];
double P_HeII[UV_BINS];
double nu[UV_BINS + 1];
#if defined(BIERMANN_BATTERY) || defined(DURRIVE_BATTERY)
double Pelec_HI[UV_BINS];   /* Momentum injection from photo-ionisation cross-section for HI   */
double Pelec_HeI[UV_BINS];  /* Momentum injection from photo-ionisation cross-section for HeI  */
double Pelec_HeII[UV_BINS]; /* Momentum injection from photo-ionisation cross-section for HeII */
#endif
#endif
double lum[MRT_BINS];
double MeanPhotonEnergy[UV_BINS];
double c_internal_units;
#ifdef MRT_IR
double radiation_constant;
#endif
#endif

#ifdef SFR_MCS_LOG
gsl_histogram *sf_dens_hist;
#endif
#ifdef SN_MCS_LOG
gsl_histogram *sn_dens_hist;
#endif

#if defined(SN_MCS) && !(defined(SN_MCS_SINGLE_INJECTION) || defined(IMF_SAMPLING_MCS))
struct sb99_data sb99;
#endif

#if defined(IMF_SAMPLING_MCS)
struct star_properties_table StarProperties;
#endif

peanokey *DomainKeyBuf;

/* global state of system
 */
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;

/* Various structures for communication during the gravity computation.
 */

struct directdata *DirectDataIn, *DirectDataAll;

struct accdata *DirectAccOut, *DirectAccIn;

int ThreadsNexport[NUM_THREADS], ThreadsNexportNodes[NUM_THREADS];

struct data_partlist *PartList, *ThreadsPartList[NUM_THREADS];

struct datanodelist *NodeList, *ThreadsNodeList[NUM_THREADS];

struct potdata_out *PotDataResult, /*!< holds the partial results computed for imported particles. Note: We use GravDataResult =
                                      GravDataGet, such that the result replaces the imported data */
    *PotDataOut; /*!< holds partial results received from other processors. This will overwrite the GravDataIn array */

/*! Header for the standard file format.
 */
struct io_header header; /*!< holds header for snapshot files */
#ifdef NTYPES_ICS
struct io_header_ICs header_ICs; /*!< holds header for IC files */
#endif

char (*Parameters)[MAXLEN_PARAM_TAG];
char (*ParametersValue)[MAXLEN_PARAM_VALUE];
char *ParametersType;

/*
 * Variables for Tree
 * ------------------
 */

/** Variables for gravitational tree */
int Tree_MaxPart;
int Tree_NumNodes;
int Tree_MaxNodes;
int Tree_FirstNonTopLevelNode;
int Tree_NumPartImported;
int Tree_NumPartExported;
int Tree_ImportedNodeOffset;
int Tree_NextFreeNode;
MyDouble *Tree_Pos_list;
unsigned long long *Tree_IntPos_list;
int *Tree_Task_list;
int *Tree_ResultIndexList;

struct treepoint_data *Tree_Points;
struct resultsactiveimported_data *Tree_ResultsActiveImported;

#ifdef BLACK_HOLES
struct treepoint_aux_bh_data *Tree_AuxBH_Points;
int Tree_NumBHExported, Tree_NumBHImported;
#endif

#ifdef SINKS
struct treepoint_aux_sinks_data *Tree_AuxSinks_Points;
int Tree_NumSinksExported, Tree_NumSinksImported;
#endif

int *Nextnode; /*!< gives next node in tree walk  (nodes array) */
int *Father;   /*!< gives parent node in tree (Prenodes array) */

struct NODE *Nodes; /*!< points to the actual memory allocted for the nodes */
                    /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
                       gives the first allocated node */

#ifdef MULTIPLE_NODE_SOFTENING
struct ExtNODE *ExtNodes;
#endif

float *Nodes_GravCost;

/** Variables for neighbor tree */
int Ngb_MaxPart;
int Ngb_NumNodes;
int Ngb_MaxNodes;
int Ngb_FirstNonTopLevelNode;
int Ngb_NextFreeNode;
int *Ngb_Father;
int *Ngb_Marker;
int Ngb_MarkerValue;

int *Ngb_DomainNodeIndex;
int *DomainListOfLocalTopleaves;
int *DomainNLocalTopleave;
int *DomainFirstLocTopleave;
#ifdef AMR
int *Ngb_DomainTask;
#endif
int *Ngb_Nextnode;

/** The ngb-tree data structure
 */
struct NgbNODE *Ngb_Nodes;
struct ExtNgbNODE *ExtNgb_Nodes;

#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif

ray_data *Ray;
int Nray;
int MaxNray;

#ifdef NUM_THREADS
int MaxThreads = NUM_THREADS;
#else
int MaxThreads = 1;
#endif

IO_Field *IO_Fields;
int N_IO_Fields   = 0;
int Max_IO_Fields = 0;
