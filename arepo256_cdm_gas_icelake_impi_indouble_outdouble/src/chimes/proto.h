#include <hdf5.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include "allvars.h"

#define ELECTRON_MASS 9.1093829e-28
#define PROTON_MASS 1.6726218e-24
#define PI 3.1415927
#define LIGHTSPEED 3.0e10          /* In cm/s */
#define BOLTZMANNCGS 1.3806e-16    /* In ergs/K */
#define BOLTZMANN_EVK 8.6173324e-5 /*In eV/K*/
#define G0_GAMMA                                                                    \
  2.77 /* For dust processes involving G0, e.g. photoelectric heating, we attenuate \
        * G0 by exp(- G0_gamma * Av) */

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/* The following structure contains
 * information about each non-eq species */

struct Species_Structure
{
  int include_species;
  double element_abundance;
  double creation_rate;
  double destruction_rate;
};

double k1(int T_index, double dT);
double k2(int T_index, double dT);
double k3(int T_index, double dT);
double k5(int T_index, double dT);
double k6(int T_index, double dT);
double k7(int T_index, double dT);
double k10(int T_index, double dT);
double k11(int T_index, double dT);
double k13(int T_index, double dT, double HI_column_density);
double k15(int T_index, double dT);
double k16(int T_index, double dT);
double k17(int T_index, double dT);
double k24(int T_index, double dT);
double k25rr(int T_index, double dT, double HeI_column_density);
double k25di(int T_index, double dT);
double k26(int T_index, double dT);
double k27(int T_index, double dT);
double k83(int T_index, double dT);
double k87(int T_index, double dT);
double k97(int T_index, double dT);
double k101(int T_index, double dT);
double k106(int T_index, double dT);
double k107(int T_index, double dT);
double k112(int T_index, double dT);
double k113(int T_index, double dT);
double k114(int T_index, double dT);
double k116(int T_index, double dT);
double k117(int T_index, double dT);
double k120(int T_index, double dT);
double k121(int T_index, double dT);
double k122(int T_index, double dT);
double k124(int T_index, double dT);
double k128(int T_index, double dT);
double k129(int T_index, double dT);
double k130(int T_index, double dT);
double k131(int T_index, double dT);
double k133(int T_index, double dT);
double k134(int T_index, double dT);
double k135(int T_index, double dT);
double k136(int T_index, double dT);
double k137(int T_index, double dT);
double k138(int T_index, double dT);
double k139(int T_index, double dT);
double k142(int T_index, double dT);
double k147(int T_index, double dT);
double k150(int T_index, double dT);
double k157(int T_index, double dT);
double k186(int T_index, double dT);
double k187(int T_index, double dT);
double k188(int T_index, double dT);
double k189(int T_index, double dT);
double k190(int T_index, double dT);
double k191(int T_index, double dT);
double k192(int T_index, double dT);
double k193(int T_index, double dT);
double k194(int T_index, double dT);
double k195(int T_index, double dT);
double k196(int T_index, double dT);
double k197(int T_index, double dT);
double k198(int T_index, double dT);
double k199(int T_index, double dT);
double k200(int T_index, double dT);
double k201(int T_index, double dT);
double k202(int T_index, double dT);
double k203(int T_index, double dT);
double k204(int T_index, double dT);
double k205(int T_index, double dT);
double k206(int T_index, double dT);
double k207(int T_index, double dT);
double k220(int T_index, double dT);
double k221(int T_index, double dT);
double k222(int T_index, double dT);
double k223(int T_index, double dT);
double k225(int T_index, double dT);
double k226(int T_index, double dT);
double k227(int T_index, double dT);
double k228(int T_index, double dT);
double k229(int T_index, double dT);
double k230(int T_index, double dT);
double k231(int T_index, double dT);
double k232(int T_index, double dT);
double k233(int T_index, double dT);
double k234(int T_index, double dT);
double k235(int T_index, double dT);
double k236(int T_index, double dT);
double k237(int T_index, double dT);
double k238(int T_index, double dT);
double k283(int T_index, double dT);
double k44(void);
double grain_recomb_jac(double T, double xe, double nH_tot, double dust_G, double extinction, double dust_ratio, double a, double b,
                        double c, double d, double e);
double k61(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k63(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k64(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k65(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k66(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k77(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k79(void);
double k80(void);
double k4(void);
double k8(double T, double nHtotal, double xH, double xH2, double xHe);
double k9(double T, double nHtotal, double xH, double xH2, double xHe);
double R52(double Go, double extinction);
double R53(double NH2, double temperature, double *radiation_intensity, double extinction, double b_5,
           struct globalVariables *myGlobalVars, struct gasVariables *myGasVars);
double k60(double T, double Tgr, double xHI, double dust_ratio);
double R70(double zeta_H);
double k84(int T_index, double dT);
double k85(double T);
double k86(double T);
double k88(double T, double nHtotal, double xH, double xH2, double xHe);
double R89(double zeta_H);
double k100(void);
double k108(void);
double R111(double dust_G, double extinction);
double k115(void);
double k118(void);
double k119(void);
double k123(void);
double k125(void);
double k126(void);
double k127(void);
double k132(void);
double k140(void);
double k141(void);
double k143(void);
double k144(void);
double k145(void);
double k146(void);
double k148(void);
double k149(void);
double k151(void);
double k152(void);
double k153(void);
double k154(void);
double k155(void);
double k156(void);
double k158(void);
double k159(void);
double k160(void);
double k161(void);
double k162(void);
double k163(void);
double k164(void);
double k165(void);
double k166(void);
double k167(void);
double k168(void);
double k169(void);
double k170(void);
double k171(void);
double k172(void);
double k173(void);
double k174(void);
double k175(void);
double k176(void);
double k177(void);
double k178(void);
double k179(void);
double k180(void);
double k181(void);
double k182(void);
double k183(void);
double k184(void);
double k185(void);
double k208(void);
double k209(void);
double k210(void);
double k211(void);
double k212(void);
double k213(void);
double k214(void);
double k215(void);
double k216(void);
double k217(void);
double k218(void);
double k219(void);
double k224(void);
double R239(double Go, double extinction);
double R240(double Go, double extinction);
double R241(double Go, double extinction);
double R242(double Go, double extinction);
double R243(double Go, double extinction);
double R244(double Go, double extinction);
double R245(double Go, double extinction);
double R246(double Go, double extinction);
double R247(double Go, double extinction);
double R248(double Go, double extinction);
double R249(double Go, double extinction);
double R250(double Go, double extinction);
double R251(double Go, double extinction);
double R252(double Go, double extinction);
double R253(double Go, double extinction);
double R254(double Go, double extinction);
double R255(double Go, double extinction);
double R256(double Go, double extinction);
double R257(double Go, double extinction);
double R258(double Go, double extinction);
double R259(double Go, double extinction);
double R260(double Go, double extinction);
double R261(double Go, double extinction);
double R262(double Go, double extinction);
double R263(double Go, double extinction);
double R264(double Go, double extinction);
double R265(double Go, double extinction);
double R266(double Go, double extinction);
double R267(double Go, double NH2, double NCO, double extinction, struct globalVariables *myGlobalVars);
double R268(double zeta_H);
double R269(double zeta_H);
double R270(double zeta_H);
double R272(double zeta_H);
double R273(double zeta_H);
double R274(double zeta_H);
double R275(double zeta_H);
double R276(double zeta_H);
double R277(double zeta_H);
double R278(double zeta_H);
double R279(double zeta_H);
double R280(double zeta_H);
double R281(double zeta_H, double T, double xH2, double xCO);
double k290(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k291(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k292(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k293(double T, double ne, double dust_G, double extinction, double dust_ratio);
double k294(void);
double k295(void);
double k296(void);
double k297(void);
double k303(void);
double k304(void);
double cr_secondary_ionisation(double xHII, int ion_index);

void set_species_arrays(struct Species_Structure *mySpecies, struct gasVariables *myGasVars, int *total_network,
                        int *nonmolecular_network, struct globalVariables *myGlobalVars);
struct Reactions_Structure *add_new_reaction(struct Reactions_Structure *conductor, int reactant1, int reactant2, int reactant3,
                                             int N_reactants, int product1, int product2, int product3, int N_products, double *rate,
                                             int incl_mol, int *speciesIndices);
void set_constraint_equations(void *user_data, double *output_array);

void set_kin_scale_vector(N_Vector scale_vector, int *enum_indices, struct gasVariables *myGasVars,
                          struct globalVariables *myGlobalVars);
int kin_f(N_Vector y, N_Vector ydot, void *user_data);
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

double photoheating(double xi, double Gamma, double epsilon, double nH);
double cosmic_ray_heating(double xi, double zeta_i, double nH);
double compton_cooling(double T, double Tcmb, double xe, double nH);
double photoelectric_grain_heating(double T, double nH, double ne, double n_total, double Z, double dust_G, double extinction);
double gas_grain_transfer(double T, double Tgr, double Z);
double grain_surface_recombination(double T, double Z, double ne, double nH, double dust_G, double extinction);
double H2_rovibrational_cooling(double T, double xH2, double xH, double xHII, double xHe, double xe, double nH_tot);
double H2_collis_deexc_rates(double T, double J, int mode);
double H2_crit_density(double T, double xH, double xH2);
double CO_rotational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                             struct globalVariables *myGlobalVars);
double CO_vibrational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                              struct globalVariables *myGlobalVars);
double H2O_rotational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                              struct globalVariables *myGlobalVars);
double H2O_vibrational_cooling(double T, double N, double xm, double xH2, double xHI, double xe, double nH,
                               struct globalVariables *myGlobalVars);
double OH_rotational_cooling(double T, double N, double dv, double xm, double nH, double n, double tau_d);
double calculate_total_number_density(double *my_abundances, double nH, struct globalVariables *myGlobalVars);

void allocate_molecular_header_arrays(struct globalVariables *myGlobalVars);
void ReadMolecularCoolingHeader(struct globalVariables *myGlobalVars);
void GetMolecularCoolingTable(struct globalVariables *myGlobalVars);
void GetPhotoIonTables(struct globalVariables *myGlobalVars, int N_Elements_in_Bens_tables,
                       struct All_rate_variables_structure **this_all_rates, double *dustG_arr, double *H2_dissocJ_arr);
void initialise_bens_tables(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates,
                            double *dustG_arr, double *H2_dissocJ_arr);
void set_equilibrium_abundances(void *user_data);
void set_equilibrium_abundances_from_tables(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void do_equilibrium_cooling(void *user_data);

void get_index_1d_mydbl(double *table, int ntable, double x, int *i, double *dx);
void get_index_1d_irregular(double *table, int ntable, double x, int *i, double *dx);
double interpol_1d_mydbl(double *table, int i, double dx);
double interpol_1d_fltdbl(float *table, int i, double dx);
double interpol_2d_mydbl(double **table, int i, int j, double dx, double dy);
double interpol_2d_fltdbl(float **table, int i, int j, double dx, double dy);
double interpol_3d_mydbl(double ***table, int i, int j, int k, double dx, double dy, double dz);
double interpol_3d_special(double ****table, int i, int j, int k, int l, double dx, double dy, double dz);
double interpol_4d_mydbl(double ****table, int i, int j, int k, int l, double dx, double dy, double dz, double dw);
double interpol_5d_mydbl(double *****table, int i, int j, int k, int l, int m, double dx, double dy, double dz, double dw, double dv);

double determine_subcyclestep(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, int *abundant_indices,
                              struct Species_Structure *species, struct Reactions_Structure *reactions_root,
                              struct All_rate_variables_structure *this_all_rates);
int is_x_in_array(int x, int *my_array, int len_array);
void update_species_fluxes(double *abundances, double temperature, struct Reactions_Structure *reactions_root,
                           struct Species_Structure *species, double Tmol, int mode, int *speciesIndices,
                           struct globalVariables *myGlobalVars);
double evaluate_reduced_network_size(struct Species_Structure *species, struct gasVariables *myGasVars,
                                     struct globalVariables *myGlobalVars, struct Reactions_Structure *reaction_list_node,
                                     int *myNetworkSize, struct All_rate_variables_structure *this_all_rates,
                                     struct Reactions_Structure *this_all_reactions_root,
                                     struct Reactions_Structure *this_nonmolecular_reactions_root);

typedef struct
{
  struct gasVariables *myGasVars;
  struct globalVariables *myGlobalVars;
  struct Species_Structure *species;
  struct Reactions_Structure *root_reactions;
  struct All_rate_variables_structure *this_all_rates;
  double *HI_column;
  double *H2_column;
  double *HeI_column;
  double *HeII_column;
  double *CO_column;
  double *H2O_column;
  double *OH_column;
  double *extinction;
  void *cvode_mem;
  int network_size;
} UserData;
