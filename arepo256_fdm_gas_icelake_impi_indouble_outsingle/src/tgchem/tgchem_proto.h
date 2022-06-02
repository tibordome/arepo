/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgchem/tgchem_proto.h
 * \date        01/2013
 * \author      Primordial chemistry and cooling network
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

// In tgchem.c
void tgchem(void);
void tgchem_compute_step_constants(void);

// In tgchem_init.c
void tgchem_begrun(void);
void tgchem_init_cvode(void);
void tgchem_finish_cvode(void);
void tgchem_init_rates(void);

// In tgchem_rates.c
void tgchem_rates(int mode, step_data* pstep, cell_data* pcell);
void tgchem_rates_aux(cell_data* pcell, double dt_step);

// In tgchem_step.c
void tgchem_step(step_data* pstep, cell_data* pcell);
int tgchem_dspecies(double time, N_Vector species, N_Vector dspecies, void* ppointer);

// In tgchem_utils.c
void tgchem_step_data(step_data* pstep);
void tgchem_cell_data(step_data* pstep, cell_data* pcell, int idx);
void tgchem_var_data(cell_data* pcell, var_data* pvar);
void tgchem_dep_vars(double abh2, double abhii, double nh, double* abhi, double* abe, double* mu, double* ntot);
void tgchem_gamma_temp(double abh2, double abe, double energy, double ntot, double* gamma, double* temp);
double tgchem_eq_coeff(int i, double temp);
double tgchem_chem_coeff(int i, double temp);
double tgchem_cool_coeff(int i, double temp);
double tgchem_opac_coeff(int i, double temp);
int tgchem_compute_temp_idx(double temp);
double tgchem_compute_dtemp(double temp, int temp_idx);
double tgchem_limit_abhm(double abhm);
double tgchem_limit_abh2(double abh2);
double tgchem_limit_abhii(double abhii);
double tgchem_limit_energy(double energy);
double tgchem_limit_abhi(double abhi);
double tgchem_limit_temp(double temp);
double tgchem_limit_h2_column(double h2_column);
double tgchem_t_ff(double nh);
void tgchem_lu_solve(void);
void tgchem_debug(int mode, double dt_step, step_data* pstep, cell_data* pcell);
