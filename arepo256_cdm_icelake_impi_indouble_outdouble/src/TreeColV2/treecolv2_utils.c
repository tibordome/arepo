/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 *  \file        src/TreeColV2/treecolv2_utils.c
 *  \date        03/07/2018
 *  \author      Paul C Clark
 *  \brief       includes the functions used by the TREECOLV2 module
 *  \details
 *
 *
 *  \par Major modifications and contributions:
 *
 *  03.07.2018  These functions represent the NEW TreeCol (v2), which uses a lookup table
 *  to calculate the column subtended by a treenode/particle during the tree walk. Lookup
 *  tables can be obtained from Paul Clark via email (paul.clark@astro.cf.ac.uk)
 *
 */

#include "../allvars.h"
#include "../proto.h"

int treecolv2_read_and_setup_lookup_table(void)
{
  const char *fname = "TreeCol_lookup.dat";
  mpi_printf("TREECOL-V2: About to try and read the look-up table info from file %s\n", fname);
  FILE *fd;
  if(!(fd = fopen(fname, "r")))
    {
      mpi_printf("TREECOL-V2: Could not find file %s. Aborting!\n", fname);
      return 0;
    }

  /* Read in the data from file */

  my_fread(&TCV2_nang_theta, sizeof(int), 1, fd);
  my_fread(&TCV2_nang_phi, sizeof(int), 1, fd);
  mpi_printf("TREECOL-V2: angle to pix map nang_theta %d nang_phi %d total num angles %d  \n", TCV2_nang_theta, TCV2_nang_phi,
             TCV2_nang_theta * TCV2_nang_phi);

  TCV2_map_theta_phi_to_ipix = (int *)mymalloc("TCV2_map_theta_phi_to_ipix", TCV2_nang_theta * TCV2_nang_phi * sizeof(int));
  my_fread(&TCV2_map_theta_phi_to_ipix[0], sizeof(int), TCV2_nang_theta * TCV2_nang_phi, fd);
  mpi_printf("TREECOL-V2: allocated and read in TCV2_map_theta_phi_to_ipix \n");

  my_fread(&TCV2_npos, sizeof(int), 1, fd);
  my_fread(&TCV2_nang, sizeof(int), 1, fd);
  my_fread(&TCV2_nentries, sizeof(int), 1, fd);
  TCV2_nentries++;
  my_fread(&TCV2_min_nodesize, sizeof(double), 1, fd);
  my_fread(&TCV2_max_nodesize, sizeof(double), 1, fd);

  mpi_printf("TREECOL-V2: npos %d nang %d nentries %d min nodesize %g max nodesize %g \n", TCV2_npos, TCV2_nang, TCV2_nentries,
             TCV2_min_nodesize, TCV2_max_nodesize);

  TCV2_start_ind = (int *)mymalloc("TCV2_start_ind", TCV2_npos * TCV2_nang * sizeof(int));
  my_fread(&TCV2_start_ind[0], sizeof(int), TCV2_npos * TCV2_nang, fd);

  TCV2_npix = (int *)mymalloc("TCV2_npix", TCV2_npos * TCV2_nang * sizeof(int));
  my_fread(&TCV2_npix[0], sizeof(int), TCV2_npos * TCV2_nang, fd);

  TCV2_pixel_list = (int *)mymalloc("TCV2_pixel_list", TCV2_nentries * sizeof(int));
  my_fread(&TCV2_pixel_list[0], sizeof(int), TCV2_nentries, fd);

  TCV2_pixel_values = (double *)mymalloc("TCV2_pixel_values", TCV2_nentries * sizeof(double));
  my_fread(&TCV2_pixel_values[0], sizeof(double), TCV2_nentries, fd);

  fclose(fd);

  mpi_printf("TREECOL-V2: Testing pixel values: TCV2_pixel_list[TCV2_nentries-1] %d\n", TCV2_pixel_list[TCV2_nentries - 1]);
  mpi_printf("TREECOL-V2: Testing pixel values: TCV2_pixel_values[67922] %g TCV2_pixel_values[67923] %g \n", TCV2_pixel_values[67922],
             TCV2_pixel_values[67923]);

  /* Set some other values */

  TCV2_dnodesize = (TCV2_max_nodesize - TCV2_min_nodesize) / ((double)(TCV2_nang - 1));
  TCV2_dtheta    = M_PI / ((double)(TCV2_nang_theta - 1));
  TCV2_dphi      = 2.0 * M_PI / ((double)(TCV2_nang_phi - 1));
  mpi_printf("TREECOL-V2: TCV2_dnodesize %g TCV2_dtheta %g TCV2_dphi %g \n", TCV2_dnodesize, TCV2_dtheta, TCV2_dphi);

  return 1;
}

void treecolv2_add_node_contribution(double posx, double posy, double posz
#ifdef TREECOLV2_VEL
                                     ,
                                     double vx_p, double vy_p, double vz_p, double vx_n, double vy_n, double vz_n, double v_th2
#endif
                                     ,
                                     double size, double gas_mass, double total_column_map[NPIX]
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
                                     int idebug)
{
  /* return if node mass is empty */
  if(gas_mass <= 0)
    return;

  /* get distance to node and return if self! */
  double rad = sqrt(posx * posx + posy * posy + posz * posz);
  if(rad <= 0)
    return;

  /* if node is further than some maxium distance, return*/
  double ascale = 1.0;
  if(All.ComovingIntegrationOn)
    ascale = All.Time;
  if(rad > All.TreeColMaxDistance * All.HubbleParam / ascale)
    return;

  /* Work out the column densities */
  double recip_node_size_squ = 1.0 / size / size;
  double total_column        = gas_mass * recip_node_size_squ;
#ifdef TREECOLV2_H2
  double h2_column = h2_mass * recip_node_size_squ;
#endif
#ifdef TREECOLV2_CO
  double co_column = co_mass * recip_node_size_squ;
#endif
#ifdef TREECOLV2_C
  double c_column = c_mass * recip_node_size_squ;
#endif
  double i_projection_h2 = 1; /* Default is to project everything */
  double i_projection_co = 1;
  double i_projection_c  = 1;

  /* Allow for doppler shifting of the lines */
#ifdef TREECOLV2_VEL
  /* Compute squared velocity difference between particle and node. In comoving runs, we assume that the peculiar velocity
     dominates over the effects of the Hubble flow; this should generally be a good approximation for the gas
     that dominates the local shielding
  */
  double v_diff2 = (vx_p - vx_n) * (vx_p - vx_n) + (vy_p - vy_n) * (vy_p - vy_n) + (vz_p - vz_n) * (vz_p - vz_n);
  v_diff2 *= All.cf_a2inv; /* Convert to physical velocity; note that in non-comoving runs, All.cf_a2inv=1 */
  double f_overl       = All.FracOverlap;
  double v_th2_compare = v_th2 * f_overl * f_overl;

  /* Now need to rescale thermal velocity for H2 or CO */
  double v_th2_compare_H2 = v_th2_compare / 2.0;
  double v_th2_compare_CO = v_th2_compare / 28.0;
  if(v_diff2 > v_th2_compare_H2) /* Add matter only if relative velocities are small enough */
    {
      i_projection_h2 = 0; /* Do not project */
      i_projection_co = 0; /* Thermal velocity of CO < thermal velocity of H2  */
    }
  else if(v_diff2 > v_th2_compare_CO)
    i_projection_co = 0;
#endif

  /* Select one of the predetemined node sizes */
  int isize = round((size / rad - TCV2_min_nodesize) / TCV2_dnodesize);
  isize     = fmax(isize, 0);
  isize     = fmin(isize, TCV2_nang - 1);

  /* Select one of the predetermined node positions from the look-up table */
  posx = posx / rad;
  posy = posy / rad;
  posz = posz / rad;

  double newrad = sqrt(posx * posx + posy * posy + posz * posz);

  double theta = acos(posz);

  double phi = atan2(posy, posx);
  if(phi < 0)
    phi += 2 * M_PI;

  int itheta = round(theta / TCV2_dtheta);
  if(itheta < 0 && itheta >= TCV2_nang_theta)
    terminate("Problem with itheta %d pos %g %g %g size %g", itheta, posx, posy, posz, size);

  int iphi = round(phi / TCV2_dphi);
  if(iphi < 0 && iphi >= TCV2_nang_phi)
    terminate("Problem with iphi %d pos %g %g %g size %g", iphi, posx, posy, posz, size);

  int iangle = TCV2_nang_phi * itheta + iphi;

  int ipos = TCV2_map_theta_phi_to_ipix[iangle];
  if(ipos < 0 && ipos >= TCV2_npos)
    terminate("Problem with ipos %d pos %g %g %g size %g", ipos, posx, posy, posz, size);

  /* retrieve and add this precomuted map from the lookup table to the column map(s) */
  int ilook = ipos * TCV2_nang + isize;

  if(ilook < 0 && ilook >= TCV2_npos * TCV2_nang)
    printf("Problem with ilook %d pos %g %g %g size %g \n", ilook, posx, posy, posz, size);
  int start_index = TCV2_start_ind[ilook];

  int ncount          = 0;
  double column_added = 0;
  int index;
  int ipix;
  double weighting;
  for(int i = 0; i < TCV2_npix[ilook]; i++)
    {
      ncount++;
      /* retrieve */
      index = start_index + i;
      if(index < 0 || index >= TCV2_nentries)
        printf(
            "TCV2: Problem with index %d pos %g %g %g size %g itheta %d iphi %d iangle %d ipos %d isize %d ilook %d npix %d "
            "start_index %d \n",
            index, posx * rad, posy * rad, posz * rad, size, itheta, iphi, iangle, ipos, isize, ilook, TCV2_npix[ilook], start_index);
      ipix = TCV2_pixel_list[index];
      if(ipix < 0 || ipix >= NPIX)
        printf(
            "TCV2: Problem with ipix %d os %g %g %g size %g itheta %d iphi %d iangle %d ipos %d isize %d ilook %d npix %d start_index "
            "%d \n",
            ipix, posx * rad, posy * rad, posz * rad, size, itheta, iphi, iangle, ipos, isize, ilook, TCV2_npix[ilook], start_index);
      weighting = TCV2_pixel_values[index];

      /* add to maps */
      column_added += weighting;
      total_column_map[ipix] += weighting * total_column;
#ifdef TREECOLV2_H2
      h2_column_map[ipix] += weighting * h2_column * i_projection_h2;
#endif
#ifdef TREECOLV2_CO
      co_column_map[ipix] += weighting * co_column * i_projection_co;
#endif
#ifdef TREECOLV2_C
      c_column_map[ipix] += weighting * c_column * i_projection_c;
#endif
    }

  if(ncount == 0 || column_added < 1e-30)
    {
      mpi_printf("TCV2: problem -- didn't actually add any column! column added %g ncount %d \n", column_added, ncount);
      terminate("Stopping...");
    }

  return;
}
