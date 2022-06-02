#include <stdio.h>
#include <stdlib.h>
#include "../dtypes.h"
#include "arepoconfig.h"

#ifndef CHIMES_ALLVARS_H
#define CHIMES_ALLVARS_H

#define TOTSIZE 157 /* The total number of species in the full network */
#define MOLSIZE 20  /* Total number of molecules in the full network */
#define EL_NAME_LENGTH 20
#define NONEQION_NAME_LENGTH 29
#define DUSTEFFSIZE 4.0e-22 /* in cm^2 / H atom. Used to convert between Av & N_Htot */
#define DUST_CROSS_SECTION                                    \
  1.0e-21 /* in cm^2 / H atom. Used in H2 formation on dust.  \
           * NOTE: NOT the same as DUSTEFFSIZE. The latter is \
           * Av / NH and thus includes absorption AND scattering. */
#define MAX_TEMPERATURE_FOR_RATES                              \
  2.0e9 /* If the temperature in gasVars exceeds this value,   \
         * we limit the temperature used to calculate reaction \
         * and cooling rates (in set_rates.c and cooling.c) to \
         * this value. */

/* This structure contains variables that
 * are specific to each gas particle/cell. */

struct gasVariables
{
  /* NOTE: all abundances are in
   * the form ni / nHtot, i.e. the
   * ratio of number densities of
   * species i to total hydrogen. */
  int index;
  double element_abundances[10]; /* In the order He, C, N, O, Ne, Mg, Si, S, Ca, Fe */
  double nH_tot;                 /* cm^-3 */
  double temperature;            /* K */
  double TempFloor;
  double divVel;                    /* s^-1 */
  double doppler_broad;             /* km s^-1. NOTE: this is ONLY from turbulence; thermal broadening is added later. */
  double *isotropic_photon_density; /* cm^-3 */
  double *dust_G_parameter;
  double *H2_dissocJ; /* This is n / (isotropic_photon_density * c),
                       * where n is photon number density in the
                       * band 12.24 eV to 13.51 eV (all in cgs units). */
  double cr_rate;
  double metallicity;    /* Z / Z_sol */
  double cell_size;      /* cm; use kernel smoothing length in SPH */
  double hydro_timestep; /* s */
  int ForceEqOn;
  int ThermEvolOn;
  double constant_heating_rate; /* erg s^-1 cm^-3 */
  double *abundances;           /* The size of this array will be set by init_chimes() */
#ifdef CHIMES_ADVECT_ABUNDANCES
  double *IonAdvect;
  double HtotAdvect;
#endif
};

/* This structure contains the global
 * variables */

struct globalVariables
{
  /* The following are parameters that will
   * need to be read in from the parameter file. */
  char BenTablesPath[500];
  char PhotoIonTablePath[500];
  char EqAbundanceTablePath[500];
  char MolecularTablePath[500];
  char AdditionalRatesTablesPath[500];
  int reductionOn;
  int updatePhotonFluxOn;
  int cellSelfShieldingOn;
  int N_spectra; /* The number of UV spectra. */
  int StaticMolCooling;
  double T_EqThresh; /* This is used when ForceEqOn == 1 */
  double time_tolerance;
  double min_subcyclestep;
  double T_mol;
  int n_ions_low;
  int n_ions_med;
  int n_ions_high;
  int InitIonState;
  double grain_temperature;
  double cmb_temperature;
  double relativeTolerance;
  double absoluteTolerance;
  double thermalAbsoluteTolerance;
  int element_included[9]; /* Metals only */
  int speciesIndices[TOTSIZE];
  int totalNumberOfSpecies;
  int print_debug_statements;
  int scale_metal_tolerances;
};

/* The following contain the rates tables */

/* Structures for Ben's tables */

struct NonEq_Ionization
{
  int N_Elements;
  int N_Temperatures;
  int *N_Auger;
  char **ElementName;
  char **FileName;
  int *AtomicNumber;
  int *N_Ions;
  int *IonIndexBegin; /* Use this to denote the enums of each ion */
  double *Temperatures;
  double *AtomicWeights;
  int shielding_dimensions[2]; /* For the molecular shielding */
  double *H2CO_shielding_N;
  double *COself_shielding_N;
  double **CO_shielding_S;
  int secondary_ionisation_dims[1]; /* For secondary ionisation tables of HI and HeI */
  double *x_ion_fraction;
  double *n_ion_HI;
  double *n_ion_HeI;

  int shieldingColumnDimensions[1]; /* This is for the shielding by HI and HeI. */
  double *shieldingColumnDensities;
  struct NonEq_Rates *NonEqRates;
};

struct NonEq_Rates
{
  int N_Ions;
  double **alpharad;
  double **alphadi;
  double **betacoll;
  double ***sigmaphot;
  double *E_thresh;
  double *cosmicRays;
  float ****shieldFactor1D;
  float *****shieldFactor2D;

  double **cool;
  double **epsilon; /* This will store the optically thin epsilons. For optically thick, use the shielding factors. */
  double **CTHrecof;
  double **CTHionof;
  double **CTHerecof;
  double **CTHeionof;

  int *CTHrec_mask;
  int *CTHion_mask;
  int *CTHerec_mask;
  int *CTHeion_mask;
};

/* This holds the rates tables
 * for the additional reactions. */

struct Additional_Rates
{
  int N_Temperatures;
  double *Temperatures;
  double *reaction1;
  double *reaction2;
  double *reaction3;
  double *reaction5;
  double *reaction6;
  double *reaction7;
  double *reaction10;
  double *reaction11;
  double *reaction13A;
  double *reaction13B;
  double *reaction15;
  double *reaction16;
  double *reaction17;
  double *reaction24;
  double *reaction25rrA;
  double *reaction25rrB;
  double *reaction25di;
  double *reaction26;
  double *reaction27;
  double *reaction83;
  double *reaction84;
  double *reaction87;
  double *reaction97;
  double *reaction101;
  double *reaction106;
  double *reaction107;
  double *reaction112;
  double *reaction113;
  double *reaction114;
  double *reaction116;
  double *reaction117;
  double *reaction120;
  double *reaction121;
  double *reaction122;
  double *reaction124;
  double *reaction128;
  double *reaction129;
  double *reaction130;
  double *reaction131;
  double *reaction133;
  double *reaction134;
  double *reaction135;
  double *reaction136;
  double *reaction137;
  double *reaction138;
  double *reaction139;
  double *reaction142;
  double *reaction147;
  double *reaction150;
  double *reaction157;
  double *reaction186;
  double *reaction187;
  double *reaction188;
  double *reaction189;
  double *reaction190;
  double *reaction191;
  double *reaction192;
  double *reaction193;
  double *reaction194;
  double *reaction195;
  double *reaction196;
  double *reaction197;
  double *reaction198;
  double *reaction199;
  double *reaction200;
  double *reaction201;
  double *reaction202;
  double *reaction203;
  double *reaction204;
  double *reaction205;
  double *reaction206;
  double *reaction207;
  double *reaction220;
  double *reaction221;
  double *reaction222;
  double *reaction223;
  double *reaction225;
  double *reaction226;
  double *reaction227;
  double *reaction228;
  double *reaction229;
  double *reaction230;
  double *reaction231;
  double *reaction232;
  double *reaction233;
  double *reaction234;
  double *reaction235;
  double *reaction236;
  double *reaction237;
  double *reaction238;
  double *reaction283;
};

extern struct rate_tables_structure
{
  /* This structure contains Ben's tables */
  struct NonEq_Ionization *NonEqIon;

  /* This structure contains tables for the
   * additional reactions. */
  struct Additional_Rates *RatesTables;

  /* These contain the cooling tables
   * for CI and OI (they be read in using
   * our init_chimes() routines). */
  int nei_cooling_table_dimensions[4];
  double *nei_cooling_temperature;
  double *nei_cooling_HIAbundance;
  double *nei_cooling_ElectronAbundance;
  double *nei_cooling_HIIAbundance;
  double ****nei_cooling_CI;
  double ****nei_cooling_OI;

  /* These contain the cooling tables
   * for CII, NII & SiII calculated
   * using Chianti */
  int chianti_cooling_table_dimensions[2];
  double *chianti_cooling_temperature;
  double *chianti_cooling_ElectronDensity;
  double **chianti_cooling_CII;
  double **chianti_cooling_NII;
  double **chianti_cooling_SiII;
  double **chianti_cooling_FeII;

  /* These are the molecular cooling
   * tables (they will also be read in
   * using our init_chimes() routines). */
  int mol_cooling_table_dimensions[11];
  double *mol_cooling_table_CO_rot_T;
  double *mol_cooling_table_CO_rot_N;
  double *mol_cooling_table_CO_vib_T;
  double *mol_cooling_table_CO_vib_N;
  double *mol_cooling_table_H2O_rot_hiT_T;
  double *mol_cooling_table_H2O_rot_hiT_N;
  double *mol_cooling_table_H2O_rot_lowT_T;
  double *mol_cooling_table_H2O_rot_lowT_N;
  double *mol_cooling_table_H2O_vib_T;
  double *mol_cooling_table_H2O_vib_N;
  double **mol_cooling_table_CO_rot_Llte;
  double **mol_cooling_table_CO_rot_nhalf;
  double **mol_cooling_table_CO_rot_a;
  double *mol_cooling_table_CO_rot_L0;
  double **mol_cooling_table_CO_vib_Llte;
  double **mol_cooling_table_H2O_rot_hiT_Llte;
  double **mol_cooling_table_H2O_rot_hiT_nhalf;
  double **mol_cooling_table_H2O_rot_hiT_a;
  double *mol_cooling_table_H2O_rot_hiT_L0;
  double **mol_cooling_table_H2Oortho_rot_lowT_Llte;
  double **mol_cooling_table_H2Oortho_rot_lowT_nhalf;
  double **mol_cooling_table_H2Oortho_rot_lowT_a;
  double *mol_cooling_table_H2Oortho_rot_lowT_L0;
  double **mol_cooling_table_H2Opara_rot_lowT_Llte;
  double **mol_cooling_table_H2Opara_rot_lowT_nhalf;
  double **mol_cooling_table_H2Opara_rot_lowT_a;
  double *mol_cooling_table_H2Opara_rot_lowT_L0;
  double **mol_cooling_table_H2O_vib_Llte;
  double *mol_cooling_table_H2_lte_temperatures;
  double *mol_cooling_table_H2_lte;
} chimesRateTables;

double calculate_mean_molecular_weight(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void allocate_gas_abundances_memory(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void free_gas_abundances_memory(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void initialise_gas_abundances(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);
void GetNonEqTables(int ns, struct globalVariables *myGlobalVars);
void GetEqAbundancesTables(struct globalVariables *myGlobalVars);
void initialise_cooling(struct globalVariables *myGlobalVars);
void initialise_additional_rates_tables(struct globalVariables *myGlobalVars);
void ReadPhotoIonTables(struct globalVariables *myGlobalVars, char *photoIonTablePath, struct NonEq_Ionization *myNonEqIon,
                        int N_Elements_in_Bens_tables, double *dustG_arr, double *H2_dissocJ_arr, int current_spectrum);
int set_species_index_array(struct globalVariables *myGlobalVariables);

/*  The following are the linked lists that
 * will contain information on each reaction
 * in the network. The routines in
 * init_chimes() will allocate memory
 * to these and set them up */
struct Reactions_Structure
{
  int no_of_reactants;
  int reactants[3];
  int no_of_products;
  int products[3];
  int flag_included;
  int extra_Auger_electrons;
  double *rate;
  struct Reactions_Structure *next_reaction;
};

struct All_rate_variables_structure
{
  /* This will be an array,
   * one for each element. */
  struct Bens_rate_structure *BensRates;

  /* Extra reactions */
  double rate_1;
  double rate_2;
  double rate_3;
  double rate_4;
  double rate_5;
  double rate_6;
  double rate_7;
  double rate_8;
  double rate_9;
  double rate_10;
  double rate_11;
  double rate_13;
  double rate_15;
  double rate_16;
  double rate_17;
  double rate_24;
  double rate_25;
  double rate_26;
  double rate_27;
  double rate_44;
  double rate_52;
  double rate_53;
  double rate_60;
  double rate_61;
  double rate_63;
  double rate_64;
  double rate_65;
  double rate_66;
  double rate_70;
  double rate_77;
  double rate_79;
  double rate_80;
  double rate_83;
  double rate_84;
  double rate_85;
  double rate_86;
  double rate_87;
  double rate_88;
  double rate_89;
  double rate_97;
  double rate_100;
  double rate_101;
  double rate_106;
  double rate_107;
  double rate_108;
  double rate_111;
  double rate_112;
  double rate_113;
  double rate_114;
  double rate_115;
  double rate_116;
  double rate_117;
  double rate_118;
  double rate_119;
  double rate_120;
  double rate_121;
  double rate_122;
  double rate_123;
  double rate_124;
  double rate_125;
  double rate_126;
  double rate_127;
  double rate_128;
  double rate_129;
  double rate_130;
  double rate_131;
  double rate_132;
  double rate_133;
  double rate_134;
  double rate_135;
  double rate_136;
  double rate_137;
  double rate_138;
  double rate_139;
  double rate_140;
  double rate_141;
  double rate_142;
  double rate_143;
  double rate_144;
  double rate_145;
  double rate_146;
  double rate_147;
  double rate_148;
  double rate_149;
  double rate_150;
  double rate_151;
  double rate_152;
  double rate_153;
  double rate_154;
  double rate_155;
  double rate_156;
  double rate_157;
  double rate_158;
  double rate_159;
  double rate_160;
  double rate_161;
  double rate_162;
  double rate_163;
  double rate_164;
  double rate_165;
  double rate_166;
  double rate_167;
  double rate_168;
  double rate_169;
  double rate_170;
  double rate_171;
  double rate_172;
  double rate_173;
  double rate_174;
  double rate_175;
  double rate_176;
  double rate_177;
  double rate_178;
  double rate_179;
  double rate_180;
  double rate_181;
  double rate_182;
  double rate_183;
  double rate_184;
  double rate_185;
  double rate_186;
  double rate_187;
  double rate_188;
  double rate_189;
  double rate_190;
  double rate_191;
  double rate_192;
  double rate_193;
  double rate_194;
  double rate_195;
  double rate_196;
  double rate_197;
  double rate_198;
  double rate_199;
  double rate_200;
  double rate_201;
  double rate_202;
  double rate_203;
  double rate_204;
  double rate_205;
  double rate_206;
  double rate_207;
  double rate_208;
  double rate_209;
  double rate_210;
  double rate_211;
  double rate_212;
  double rate_213;
  double rate_214;
  double rate_215;
  double rate_216;
  double rate_217;
  double rate_218;
  double rate_219;
  double rate_220;
  double rate_221;
  double rate_222;
  double rate_223;
  double rate_224;
  double rate_225;
  double rate_226;
  double rate_227;
  double rate_228;
  double rate_229;
  double rate_230;
  double rate_231;
  double rate_232;
  double rate_233;
  double rate_234;
  double rate_235;
  double rate_236;
  double rate_237;
  double rate_238;
  double rate_239;
  double rate_240;
  double rate_241;
  double rate_242;
  double rate_243;
  double rate_244;
  double rate_245;
  double rate_246;
  double rate_247;
  double rate_248;
  double rate_249;
  double rate_250;
  double rate_251;
  double rate_252;
  double rate_253;
  double rate_254;
  double rate_255;
  double rate_256;
  double rate_257;
  double rate_258;
  double rate_259;
  double rate_260;
  double rate_261;
  double rate_262;
  double rate_263;
  double rate_264;
  double rate_265;
  double rate_266;
  double rate_267;
  double rate_268;
  double rate_269;
  double rate_270;
  double rate_272;
  double rate_273;
  double rate_274;
  double rate_275;
  double rate_276;
  double rate_277;
  double rate_278;
  double rate_279;
  double rate_280;
  double rate_281;
  double rate_283;
  double rate_290;
  double rate_291;
  double rate_292;
  double rate_293;
  double rate_294;
  double rate_295;
  double rate_296;
  double rate_297;
  double rate_303;
  double rate_304;
};

double calculate_total_cooling_rate(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, double HI_column_density,
                                    double HeI_column_density, double HeII_column_density, double H2_column_density,
                                    double CO_column_density, double H2O_column_density, double OH_column_density, double extinction,
                                    struct All_rate_variables_structure *this_all_rates);
void init_chimes(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates,
                 struct Reactions_Structure **this_all_reactions_root, struct Reactions_Structure **this_nonmolecular_reactions_root,
                 double *dustG_arr, double *H2_dissocJ_arr);
void initialise_reactions(struct Reactions_Structure *root_node, int incl_mol, struct globalVariables *myGlobalVars,
                          struct All_rate_variables_structure *this_all_rates);
void set_constant_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars,
                        struct All_rate_variables_structure *this_all_rates);
void update_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, double HI_column_density,
                  double H2_column_density, double HeI_column_density, double HeII_column_density, double CO_column_density,
                  double extinction, struct All_rate_variables_structure *this_all_rates);
void update_T_dependent_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars,
                              struct All_rate_variables_structure *this_all_rates);
void chimes_network(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars,
                    struct All_rate_variables_structure *this_all_rates, struct Reactions_Structure *this_all_reactions_root,
                    struct Reactions_Structure *this_nonmolecular_reactions_root);
void check_constraint_equations(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars);

/* This next structure will hold the current
 * values of the rate coefficients at the
 * current temperature. The reaction list will
 * then point to the variables in this structure. */

struct Bens_rate_structure
{
  /* From Bens model */
  double *CollisIon;
  double **PhotoIon;
  double *cosmicRays;
  double *Recomb;
  double *CTHrec;
  double *CTHion;
  double *CTHerec;
  double *CTHeion;
};

extern struct Equilibrium_abundance_structure
{
  int N_Temperatures;
  int N_Densities;
  int N_Metallicities;
  double *Temperatures;
  double *Densities;
  double *Metallicities;
  double ****EqAbundances;
} EquilibriumAbundances;

enum
{
  elec,    /* 0 */
  HI,      /* 1 */
  HII,     /* 2 */
  Hm,      /* 3 */
  HeI,     /* 4 */
  HeII,    /* 5 */
  HeIII,   /* 6 */
  CI,      /* 7 */
  CII,     /* 8 */
  CIII,    /* 9 */
  CIV,     /* 10 */
  CV,      /* 11 */
  CVI,     /* 12 */
  CVII,    /* 13 */
  Cm,      /* 14 */
  NI,      /* 15 */
  NII,     /* 16 */
  NIII,    /* 17 */
  NIV,     /* 18 */
  NV,      /* 19 */
  NVI,     /* 20 */
  NVII,    /* 21 */
  NVIII,   /* 22 */
  OI,      /* 23 */
  OII,     /* 24 */
  OIII,    /* 25 */
  OIV,     /* 26 */
  OV,      /* 27 */
  OVI,     /* 28 */
  OVII,    /* 29 */
  OVIII,   /* 30 */
  OIX,     /* 31 */
  Om,      /* 32 */
  NeI,     /* 33 */
  NeII,    /* 34 */
  NeIII,   /* 35 */
  NeIV,    /* 36 */
  NeV,     /* 37 */
  NeVI,    /* 38 */
  NeVII,   /* 39 */
  NeVIII,  /* 40 */
  NeIX,    /* 41 */
  NeX,     /* 42 */
  NeXI,    /* 43 */
  MgI,     /* 44 */
  MgII,    /* 45 */
  MgIII,   /* 46 */
  MgIV,    /* 47 */
  MgV,     /* 48 */
  MgVI,    /* 49 */
  MgVII,   /* 50 */
  MgVIII,  /* 51 */
  MgIX,    /* 52 */
  MgX,     /* 53 */
  MgXI,    /* 54 */
  MgXII,   /* 55 */
  MgXIII,  /* 56 */
  SiI,     /* 57 */
  SiII,    /* 58 */
  SiIII,   /* 59 */
  SiIV,    /* 60 */
  SiV,     /* 61 */
  SiVI,    /* 62 */
  SiVII,   /* 63 */
  SiVIII,  /* 64 */
  SiIX,    /* 65 */
  SiX,     /* 66 */
  SiXI,    /* 67 */
  SiXII,   /* 68 */
  SiXIII,  /* 69 */
  SiXIV,   /* 70 */
  SiXV,    /* 71 */
  SI,      /* 72 */
  SII,     /* 73 */
  SIII,    /* 74 */
  SIV,     /* 75 */
  SV,      /* 76 */
  SVI,     /* 77 */
  SVII,    /* 78 */
  SVIII,   /* 79 */
  SIX,     /* 80 */
  SX,      /* 81 */
  SXI,     /* 82 */
  SXII,    /* 83 */
  SXIII,   /* 84 */
  SXIV,    /* 85 */
  SXV,     /* 86 */
  SXVI,    /* 87 */
  SXVII,   /* 88 */
  CaI,     /* 89 */
  CaII,    /* 90 */
  CaIII,   /* 91 */
  CaIV,    /* 92 */
  CaV,     /* 93 */
  CaVI,    /* 94 */
  CaVII,   /* 95 */
  CaVIII,  /* 96 */
  CaIX,    /* 97 */
  CaX,     /* 98 */
  CaXI,    /* 99 */
  CaXII,   /* 100 */
  CaXIII,  /* 101 */
  CaXIV,   /* 102 */
  CaXV,    /* 103 */
  CaXVI,   /* 104 */
  CaXVII,  /* 105 */
  CaXVIII, /* 106 */
  CaXIX,   /* 107 */
  CaXX,    /* 108 */
  CaXXI,   /* 109 */
  FeI,     /* 110 */
  FeII,    /* 111 */
  FeIII,   /* 112 */
  FeIV,    /* 113 */
  FeV,     /* 114 */
  FeVI,    /* 115 */
  FeVII,   /* 116 */
  FeVIII,  /* 117 */
  FeIX,    /* 118 */
  FeX,     /* 119 */
  FeXI,    /* 120 */
  FeXII,   /* 121 */
  FeXIII,  /* 122 */
  FeXIV,   /* 123 */
  FeXV,    /* 124 */
  FeXVI,   /* 125 */
  FeXVII,  /* 126 */
  FeXVIII, /* 127 */
  FeXIX,   /* 128 */
  FeXX,    /* 129 */
  FeXXI,   /* 130 */
  FeXXII,  /* 131 */
  FeXXIII, /* 132 */
  FeXXIV,  /* 133 */
  FeXXV,   /* 134 */
  FeXXVI,  /* 135 */
  FeXXVII, /* 136 */
  H2,      /* 137 */
  H2p,     /* 138 */
  H3p,     /* 139 */
  OH,      /* 140 */
  H2O,     /* 141 */
  C2,      /* 142 */
  O2,      /* 143 */
  HCOp,    /* 144 */
  CH,      /* 145 */
  CH2,     /* 146 */
  CH3p,    /* 147 */
  CO,      /* 148 */
  CHp,     /* 149 */
  CH2p,    /* 150 */
  OHp,     /* 151 */
  H2Op,    /* 152 */
  H3Op,    /* 153 */
  COp,     /* 154 */
  HOCp,    /* 155 */
  O2p      /* 156 */
};

#endif  // !(CHIMES_ALLVARS_H)
