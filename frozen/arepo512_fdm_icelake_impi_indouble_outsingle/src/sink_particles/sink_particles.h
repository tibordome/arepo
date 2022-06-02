/* PCC - 19.12.2013
	This header holds the external structures for the sink 
        particle functions.
*/

#ifdef SINK_PARTICLES_FEEDBACK
#define MAXSNE 500 /*maximum number of SNe that every sink can hold*/
#endif


extern int NSinksThisTask;
extern int NSinksAllTasks;
extern int NSinksAllTasksOld;
extern int NSinkBufferSize;
extern int SinksFormedSinceLastDomain;
extern double SinkCreationDensityCodeUnits;
extern double SinkAccretionRadius;
extern double SinkFormationRadiusSquared;
extern double SinkFormationRadius;
extern double SinkSofteningLength;
extern double SinkTargetAccretionMass;
extern double SinkTestsGiveupDensityCodeUnits;
extern double SinkEvolutionDumpRateYears;
extern double SinksLastEvolutionDumpTime;


/* This structure is the one that we'll use throughout Arepo */
extern struct global_sink_particle_data
{
  double Pos[3];
  double Vel[3];
  double Accel[3];
  double Mass;
#ifdef SINK_PARTICLES_VARIABLE_ACC_RADIUS
  double AccretionRadius;
#endif
  double FormationMass;
  double FormationTime;
#ifdef SGCHEM_ACCRETION_LUMINOSITY
  double AccretionRate;
  double TimeOld;
#endif
  int ID;
  int HomeTask;
  int Index;
  int FormationOrder;
#ifdef SINK_PARTICLES_FEEDBACK
  int N_sne;  /*number of stars that will go supernova associated to this sink*/
  double explosion_time[MAXSNE]; /*vector storing the time the massive stars will go SN*/
#endif 
} *SinkP,
  *export_SinkP;

/*The following keeps track of which Task just formed a sink
*/
extern int LastTaskNewSink;


