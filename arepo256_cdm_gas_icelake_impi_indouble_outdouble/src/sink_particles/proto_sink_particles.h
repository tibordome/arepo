/* PCC - 02.01.2014
        Contains the prototypes for the ITA-style sink particles
*/

void sink_particles(void);
void init_sink_particles(int mode);
int get_all_sink_particle_info(int mode);
double accrete_onto_sink_particles(void);
int create_sink_particles(void);
void add_sink_to_particle_structure(int idx, int candidate_index, int candidate_task);
void dump_sink_particle_info(int success);
void dump_sink_particle_snapshot(int snapNum);
int load_sink_particle_snapshot(int snapNum);
int open_sink_particle_evolution_file(void);
int open_old_sink_file_and_read_sink_data(void);

void set_sink_particle_parameters(void); /* Test outside of SINK_PARTICLES_FEEDBACK for now but it was inside before */

#ifdef SINK_MERGERS
void perform_mergers(void);
#endif

#ifdef SINK_PARTICLES_FEEDBACK
void assign_sne(int sink_ID, double mass_to_convert, int candidate_task);
void find_number_of_sne_drawing_from_IMF(double mass_to_convert, int *n_sne, double *massarr);
double stellar_lifetime(double mass_of_star);
double get_atime_of_time(double t_Gyr);
void reinject_gas_mass(int isink);
void eliminate_old_sinks(void);

#ifdef SINK_PHOTOION_FEEDBACK
void photoionisation_feedback_from_sinks(void);
double find_stromgren_radius_around_sink(int isink, int *local_n, int *total_n, int **indices);
void heat_and_ionise_stromgren_cells(int local_n, int indices[]);
#endif

#endif /*SINK_PARTICLES_FEEDBACK*/

#ifdef SINK_PARTICLES_VARIABLE_CREATION
double variable_getcreationdensity(double TargetMass);
double variable_getformationr2(double TargetMass, double CreationDensity);
#ifdef SINK_PARTICLES_FEEDBACK
double variable_geteff(double TargetMass, double Density, double CreationDensity);
#endif
#endif
