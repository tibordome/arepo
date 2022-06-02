/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Interface functions
 * \details     
 */

#ifdef SX_CHEMISTRY==3

// Returns fluxes and other properties of an gass particle
void sx_chem_cell_props( int index, double *flux, double *dens, double *utherm );

// TODO: This will be later an interface between Simplex and codes like Polaris or Warpfield

#endif
