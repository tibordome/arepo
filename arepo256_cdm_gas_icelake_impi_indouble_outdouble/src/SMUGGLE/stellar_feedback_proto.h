/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/SMUGGLE/stellar_feedback_proto.h
 * \date        03/2020
 * \author     
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifndef FM_STELLAR_FEED_PROTO_H
#define FM_STELLAR_FEED_PROTO_H

MyDouble compute_SN_energy(MyDouble Number_of_SN);
MyDouble quadratic_equation_solver(MyDouble b, MyDouble c);

double smuggle_sample_poisson_distribution(double lambda);
void smuggle_inject_snfeed_into_cell(int j, double mass, double inj_energy, double *inj_mom, double *boost, double *vel, double coolfac);
void smuggle_inject_windfeed_into_cell(int j, double mass, double inj_energy, double *inj_mom, double *boost, double *vel);
double get_feedback_radius_limiter(double gasdens, double metallicity, double initial_mass);

#ifdef SMUGGLE_OUTPUT_STELLAR_FEEDBACK
void output_stellar_feedback_statistics(void);
#endif

#endif
