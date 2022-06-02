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
#ifdef SINK_PARTICLES_FEEDBACK
void assign_sne(int sink_ID, double mass_to_convert, int candidate_task);
void find_number_of_sne_drawing_from_IMF(double mass_to_convert, int *n_sne, double *massarr);
double stellar_lifetime(double mass_of_star);
#endif
