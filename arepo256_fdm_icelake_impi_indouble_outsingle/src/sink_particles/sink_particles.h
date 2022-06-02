/* PCC - 19.12.2013
        This header holds the external structures for the sink
        particle functions.
*/

#ifdef SINK_PARTICLES_FEEDBACK
#ifdef POPIII_SNE
#define MAXSNE 35 /*maximum number of SNe that every sink can hold*/
#define MAXACCRETIONEVENTS 2
#else
#define MAXSNE 800 /*maximum number of SNe that every sink can hold*/
#define MAXACCRETIONEVENTS 50
#endif
#endif

/* Length of the list of sink_particles */
#ifdef POPIII_SNE
#define SINKLISTLENGTH 30
#else
#define SINKLISTLENGTH 30000
#endif

extern int NSinksThisTask;
extern int NSinksAllTasks;
extern int NSinksAllTasksOld;
extern int NSinkBufferSize;
extern int SinksFormedSinceLastDomain;
extern double SinkCreationDensityCodeUnits;
extern double SinkCreationDensityCurrent;
extern double SinkAccretionRadius;
extern double SinkAccretionCurrent;
extern double SinkFormationRadiusSquared;
extern double SinkFormationRadius;
extern double SinkFormationRadiusCurrent;
extern double SinkSofteningLength;
extern double SinkTargetAccretionMass;
extern double SinkTestsGiveupDensityCurrent;
extern double SinkEvolutionDumpRateYears;
extern double SinksLastEvolutionDumpTime;
extern double SinkEvolutionDumpRateCurrent;

#ifdef SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK
int NSinksLastSnapshot;
#endif

/* This structure is the one that we'll use throughout Arepo */
extern struct global_sink_particle_data
{
  double Pos[3];
  double Vel[3];
  double Accel[3];
  double Mass;
#ifdef SINK_SIMPLEX
  double MassOld;  // used to calculate average accretion over larger time periods
#endif
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
  double AccretionRadius;
#endif
  double FormationMass;
  double FormationTime;
#ifdef SGCHEM_ACCRETION_LUMINOSITY
  double AccretionRate;
  double TimeOld;
#endif
  long long ID;
  int HomeTask;
  int Index;
  int FormationOrder;
#ifdef SINK_PARTICLES_FEEDBACK
  int N_sne; /*number of stars that will go supernova associated to this sink*/
  double StellarMass;
  double explosion_time[MAXSNE]; /*vector storing the time the massive stars will go SN*/
#ifdef POPIII_SNE
  double stellar_mass[MAXSNE]; /*vector storing the masses of massive stars*/
#endif
  /*buffer for when we have too many SNe, then for the ones that do not fit within MAXSNE we store the mass
   *and formation time of the accretion event and draw from the IMF every time-step*/
  double MassStillToConvert[MAXACCRETIONEVENTS];
  double AccretionTime[MAXACCRETIONEVENTS];
#endif
#ifdef SINK_PARTICLES_VARIABLE_CREATION
  double CreationDensity;
#ifdef SINK_PARTICLES_FEEDBACK
  double SFeff;
#endif
#endif
#ifdef STORE_SINK_PARTICLE_SPIN
  double AngularMomentum[3];
#endif
} * SinkP, *export_SinkP;

/*The following keeps track of which Task just formed a sink
 */
extern int LastTaskNewSink;
