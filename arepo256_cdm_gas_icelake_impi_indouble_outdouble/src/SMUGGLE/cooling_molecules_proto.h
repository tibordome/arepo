/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/SMUGGLE/cooling_molecules_proto.h
 * \date        03/2020
 * \author      Federico Marinacci
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */


double calc_self_shielding_factor_molecular_cooling(double nH, char J_UV);
double get_MolecularCoolingRate(double logT, double nH, double logZinSolar, char J_UV);

