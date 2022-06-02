/* PCC - 12.06.2018
        Contains the prototypes for the new look-up table TreeCol
*/

int treecolv2_read_and_setup_lookup_table(void);

void treecolv2_add_node_contribution(double posx, double posy, double posz
#ifdef TREECOLV2_VEL
                                     ,
                                     double vx_p, double vy_p, double vz_p, double vx_n, double vy_n, double vz_n, double v_th2
#endif
                                     ,
                                     double size, double gasmass, double total_column_map[NPIX]
#ifdef TREECOLV2_H2
                                     ,
                                     double h2_mass, double h2_column_map[NPIX]
#endif
#ifdef TREECOLV2_CO
                                     ,
                                     double co_mass, double co_column_map[NPIX]
#endif
#ifdef TREECOLV2_C
                                     ,
                                     double c_mass, double c_column_map[NPIX]
#endif
                                     ,
                                     int idebug);
