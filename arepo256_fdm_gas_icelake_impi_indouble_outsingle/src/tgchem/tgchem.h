/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#define ELECTRON_VOLT 1.60219e-12
#define TEMP_CMB 2.725
#define TGCHEM_CHI_HI (13.6 * ELECTRON_VOLT)
#define TGCHEM_CHI_H2 (4.48 * ELECTRON_VOLT)
#define TGCHEM_CHI_HM (3.53 * ELECTRON_VOLT)

#define TGCHEM_NUM_ABUNDANCES 3
#define TGCHEM_NUM_SPECIES (TGCHEM_NUM_ABUNDANCES + 1)
#define TGCHEM_NUM_CHEM_RATES 8
#define TGCHEM_NUM_EQ 2
#define TGCHEM_NUM_CHEM 7
#define TGCHEM_NUM_COOL 6
#define TGCHEM_NUM_OPAC 5
#define TGCHEM_NUM_COOLING 5
#define TGCHEM_NUM_TEMP 5000

#define TGCHEM_MAX_ABUNDANCE_HM 1.
#define TGCHEM_MAX_ABUNDANCE_H2 0.5
#define TGCHEM_MAX_ABUNDANCE_HII 1.

#define TGCHEM_TEMP_MIN 1e1
#define TGCHEM_TEMP_MAX 1e8
#define TGCHEM_LOG_DTEMP (log10(TGCHEM_TEMP_MAX / TGCHEM_TEMP_MIN) / TGCHEM_NUM_TEMP)

#define TGCHEM_TOL 1e-5
#define TGCHEM_TOL_EQ_STEPSIZE 1e-3
#define TGCHEM_TOL_TRANS 1e-2
#define TGCHEM_TOL_BISECTION_H2 1e-1
#define TGCHEM_TOL_BISECTION_HII 1e-2
#define TGCHEM_TOL_ABHM 1e-20
#define TGCHEM_TOL_ABH2 1e-20
#define TGCHEM_TOL_ABHII 1e-20
#define TGCHEM_TOL_ABENERGY 0.

#define TGCHEM_MAX_RATE 1e100
#define TGCHEM_MAX_NUM_OUT 10000

#define TGCHEM_H2_ENERGY_SA 8
#define TGCHEM_H2_ENERGY_SB 5
#define TGCHEM_H2_NUM_V 3
#define TGCHEM_H2_NUM_J 20
#define TGCHEM_H2_NUM_T 3
#define TGCHEM_H2_TOT_NUM_LINES (TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_J * TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_T)
#define TGCHEM_H2_NUM_COLUMN 200
#define TGCHEM_H2_COLUMN_MIN 1e21
#define TGCHEM_H2_COLUMN_MAX 1e30
#define TGCHEM_H2_LOG_DCOLUMN (log10(TGCHEM_H2_COLUMN_MAX / TGCHEM_H2_COLUMN_MIN) / TGCHEM_H2_NUM_COLUMN)

#define TGCHEM_PLANCK_OPAC_NUM_RHO 15
#define TGCHEM_PLANCK_OPAC_NUM_TEMP 29

typedef struct
{
  int count;
  long num_neq;
  long num_eq_h2;
  long num_eq_hii;
  long num_rate_calls;
  long num_cells;
  double nh;
  double nh_final;
  double hydro_rate;
  double dt_hydro;
  double dt_neq;
  double dt_eq_h2;
  double dt_eq_hii;
  double temp_cmb;
} step_data;

typedef struct
{
  int cell_idx;
  int cell_id;
  int cell_task;
  double rho_int;
  double nh;
  double nh2;
  double nh3;
  double hydro_rate;
  double dt_hydro;
  double divv;
  double t_ff;
  double vth;
  double ljeans;
  double h2_column;
  double abhm;
  int eq_flag[TGCHEM_NUM_SPECIES];
  double species[TGCHEM_NUM_SPECIES];
  double prate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];
  double nrate[TGCHEM_NUM_SPECIES][TGCHEM_NUM_CHEM_RATES];
  double tot_rate[TGCHEM_NUM_SPECIES];
  double tot_chem_rate;
  double tot_rad_rate;
  double dt_rate[TGCHEM_NUM_SPECIES];
  double dt_chem;
  double dt_rad;
} cell_data;

typedef struct
{
  double abhm;
  double abh2;
  double abhii;
  double energy;
  double abhi;
  double abe;
  double mu;
  double ntot;
  double gamma;
  double temp;
} var_data;

typedef struct
{
  step_data *pstep;
  cell_data *pcell;
  var_data *pvar;
} pointer_data;

typedef struct
{
  double nh[TGCHEM_MAX_NUM_OUT];
  double gamma[TGCHEM_MAX_NUM_OUT];
  double temp[TGCHEM_MAX_NUM_OUT];
  double species[TGCHEM_NUM_SPECIES * TGCHEM_MAX_NUM_OUT];
  double rates[2 * TGCHEM_NUM_SPECIES * TGCHEM_NUM_CHEM_RATES * TGCHEM_MAX_NUM_OUT];
} out_data;

extern struct TGCD_struct
{
  // Input parameters
  int ChemMode;
  int ChemIOMode;
  int ChemH2Mode;
  double ChemInitAbH2;
  double ChemInitAbHII;
  double ChemJ21;

  // Additional Parameters
  int DebugFlag;

  // 0d parameters
  double CollapseFac;
  double TestStepSize;
  double NHInit;
  double NHFinal;
  double TempInit;
  double RedShiftInit;

  // Constants
  double AbMax[TGCHEM_NUM_ABUNDANCES];
  double EnergyConvFac;
  double H2OptThickConst;
  double H2OptThickNHThresh;
  double CIEOptThickNHThresh;
  double LWDissRate;
  double HIPhotonEnergy;

  // Logging
  int NumEq;
  int NumNeq;
  int NumRateCalls;
  int MaxNumRateCalls;
  long TotNumRateCalls;
  int NumSubSteps;
  int MaxNumSubSteps;
  long TotNumSubSteps;
  long NumCellsDone;
  double DtEq;
  double DtNEq;

  // CVODE
  void *CVODEMem;
  N_Vector SpeciesTol;

  // Tables
  double *TempTable;
  double *EqTable;
  double *DEqTable;
  double *ChemTable;
  double *DChemTable;
  double *CoolTable;
  double *DCoolTable;
  double *OpacTable;
  double *DOpacTable;

  // H2 data
  double H2Energy[TGCHEM_H2_NUM_V * TGCHEM_H2_NUM_J];
  double *H2SpontCoeff;
  double *H2ColumnTable;
  double *H2SobEmissTable;
  double *H2DSobEmissTable;
  double *H2TestSobXVal;
  double *H2TestFitEscFrac;
  double *H2TestSobEscFrac;

  double *H2EmissTable;
  double *H2DEmissTable;
  double *H2TotEmissTable;
  double *H2DTotEmissTable;
  double *H2CrossTable;
  double *H2DCrossTable;
  double *H2TotCrossTable;
  double *H2DTotCrossTable;
#ifdef HEALRAY
  int *H2LineIdxTable;
#endif

  // Opacity
  double Opac;
  double *PlanckOpacTable;

} TGCD;
