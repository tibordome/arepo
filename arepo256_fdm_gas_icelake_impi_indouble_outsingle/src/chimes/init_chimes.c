#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../allvars.h"
#include "../proto.h"
#include "proto.h"

struct Equilibrium_abundance_structure EquilibriumAbundances;
struct rate_tables_structure chimesRateTables;

void GetEqAbundancesTables(struct globalVariables *myGlobalVars)
{
  hid_t file_id, dataset, dataspace_id, memspace_id;
  herr_t status;
  hsize_t dims[2];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int rank, i, j, k, l;
  char set_name[500];
  double *T, *nH, *Z;
  double *array1D;
  dims[1] = 1;

  file_id = H5Fopen(myGlobalVars->EqAbundanceTablePath, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Read in size of each dimension */
  sprintf(set_name, "/N_Temperatures");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(EquilibriumAbundances.N_Temperatures));
  status  = H5Dclose(dataset);

  sprintf(set_name, "/N_Densities");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(EquilibriumAbundances.N_Densities));
  status  = H5Dclose(dataset);

  sprintf(set_name, "/N_Metallicities");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(EquilibriumAbundances.N_Metallicities));
  status  = H5Dclose(dataset);

  /* Allocate memory to tables. */
  EquilibriumAbundances.Temperatures  = (double *)mymalloc("ChimesEqmT", EquilibriumAbundances.N_Temperatures * sizeof(double));
  EquilibriumAbundances.Densities     = (double *)mymalloc("ChimesEqmRho", EquilibriumAbundances.N_Densities * sizeof(double));
  EquilibriumAbundances.Metallicities = (double *)mymalloc("ChimesEqmZ", EquilibriumAbundances.N_Metallicities * sizeof(double));

  EquilibriumAbundances.EqAbundances =
      (double ****)mymalloc("ChimesEqm_dim0", EquilibriumAbundances.N_Temperatures * sizeof(double ***));
  EquilibriumAbundances.EqAbundances[0] = (double ***)mymalloc(
      "ChimesEqm_dim1", EquilibriumAbundances.N_Temperatures * EquilibriumAbundances.N_Densities * sizeof(double **));
  EquilibriumAbundances.EqAbundances[0][0] =
      (double **)mymalloc("ChimesEqm_dim2", EquilibriumAbundances.N_Temperatures * EquilibriumAbundances.N_Densities *
                                                EquilibriumAbundances.N_Metallicities * sizeof(double *));
  EquilibriumAbundances.EqAbundances[0][0][0] = (double *)mymalloc(
      "ChimesEqm_dim3", EquilibriumAbundances.N_Temperatures * EquilibriumAbundances.N_Densities *
                            EquilibriumAbundances.N_Metallicities * myGlobalVars->totalNumberOfSpecies * sizeof(double));
  for(i = 0; i < EquilibriumAbundances.N_Temperatures; i++)
    {
      EquilibriumAbundances.EqAbundances[i] = &(EquilibriumAbundances.EqAbundances[0][i * EquilibriumAbundances.N_Densities]);
      for(j = 0; j < EquilibriumAbundances.N_Densities; j++)
        {
          EquilibriumAbundances.EqAbundances[i][j] =
              &(EquilibriumAbundances
                    .EqAbundances[0][0][((i * EquilibriumAbundances.N_Densities) + j) * EquilibriumAbundances.N_Metallicities]);
          for(k = 0; k < EquilibriumAbundances.N_Metallicities; k++)
            EquilibriumAbundances.EqAbundances[i][j][k] =
                &(EquilibriumAbundances
                      .EqAbundances[0][0][0]
                                   [((((i * EquilibriumAbundances.N_Densities) + j) * EquilibriumAbundances.N_Metallicities) + k) *
                                    myGlobalVars->totalNumberOfSpecies]);
        }
    }

  T       = (double *)mymalloc("Chimes_T", EquilibriumAbundances.N_Temperatures * sizeof(double));
  nH      = (double *)mymalloc("Chimes_nH", EquilibriumAbundances.N_Densities * sizeof(double));
  Z       = (double *)mymalloc("Chimes_Z", EquilibriumAbundances.N_Metallicities * sizeof(double));
  array1D = (double *)mymalloc("Chimes_array1D", EquilibriumAbundances.N_Temperatures * sizeof(double));

  /* Read in tables. */
  sprintf(set_name, "/Temperatures");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, T);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Densities");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nH);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Metallicities");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Z);
  status  = H5Dclose(dataset);

  for(i = 0; i < EquilibriumAbundances.N_Temperatures; i++)
    EquilibriumAbundances.Temperatures[i] = log10(T[i]);

  for(i = 0; i < EquilibriumAbundances.N_Densities; i++)
    EquilibriumAbundances.Densities[i] = log10(nH[i]);

  for(i = 0; i < EquilibriumAbundances.N_Metallicities; i++)
    EquilibriumAbundances.Metallicities[i] = log10(Z[i]);

  dataset      = H5Dopen(file_id, "Abundances");
  dataspace_id = H5Dget_space(dataset);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);

  dims[0]     = EquilibriumAbundances.N_Temperatures;
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < EquilibriumAbundances.N_Densities; i++)
    {
      for(j = 0; j < EquilibriumAbundances.N_Metallicities; j++)
        {
          for(k = 0; k < myGlobalVars->totalNumberOfSpecies; k++)
            {
              /* Note that if we switch off individual
               * elements, you will need to create a
               * new EqAbundances table that does NOT
               * include these excluded elements. */
              offset4D[0] = 0;
              offset4D[1] = i;
              offset4D[2] = j;
              offset4D[3] = k;
              count4D[0]  = EquilibriumAbundances.N_Temperatures;
              count4D[1]  = 1;
              count4D[2]  = 1;
              count4D[3]  = 1;
              status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);

              status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array1D);

              for(l = 0; l < EquilibriumAbundances.N_Temperatures; l++)
                EquilibriumAbundances.EqAbundances[l][i][j][k] =
                    log10(max(array1D[l], 1.0e-95)); /* If the table accidentally contains negative values, we need to force it to be
                                                        positive before taking the log. */
            }
        }
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset);

  status = H5Fclose(file_id);

  myfree(array1D);
  myfree(Z);
  myfree(nH);
  myfree(T);
}

/* This routine sets up the initial photoionisation
 * table in the chimesRateTables structure. */

void GetPhotoIonTables(struct globalVariables *myGlobalVars, int N_Elements_in_Bens_tables,
                       struct All_rate_variables_structure **this_all_rates, double *dustG_arr, double *H2_dissocJ_arr)
{
  int j, l, ns, N_arrayCells;

  /* Read in list of paths to the PhotoIon tables,
   * and read in each table in turn. */
  FILE *fdin;
  char buffer[500];
  char table_path[500];

  if((fdin = fopen(myGlobalVars->PhotoIonTablePath, "r")))
    {
      int current_spectrum = 0;
      while(fgets(buffer, 500, fdin))
        {
          if(current_spectrum > myGlobalVars->N_spectra - 1)
            {
              printf("ERROR: too many UV spectra specified in PhotoIon table list. Aborting. \n");
              exit(-1);
            }
          sscanf(buffer, "%s", table_path);
          printf("Reading PhotoIon table: %s \n", table_path);
          ReadPhotoIonTables(myGlobalVars, table_path, chimesRateTables.NonEqIon, N_Elements_in_Bens_tables, dustG_arr, H2_dissocJ_arr,
                             current_spectrum);
          current_spectrum += 1;
        }
      if(current_spectrum < myGlobalVars->N_spectra)
        {
          printf("ERROR: too few UV spectra specified in PhotoIon table list. Aborting. \n");
          exit(-1);
        }
    }
  else
    {
      printf("ERROR: PhotoIon table list %s not found. \n", myGlobalVars->PhotoIonTablePath);
      exit(-1);
    }

  /* Allocate memory for the photoionisation rates
   * in the AllRates structure. */
  ns = 0;
  for(l = 0; l < N_Elements_in_Bens_tables; l++)
    {
      if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
          else
            N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          (*this_all_rates)->BensRates[ns].PhotoIon = (double **)mymalloc("Chimes_PhotoIon_dim0", N_arrayCells * sizeof(double *));
          (*this_all_rates)->BensRates[ns].PhotoIon[0] =
              (double *)mymalloc("Chimes_PhotoIon_dim1", N_arrayCells * chimesRateTables.NonEqIon->N_Auger[ns] * sizeof(double));
          for(j = 0; j < N_arrayCells; j++)
            (*this_all_rates)->BensRates[ns].PhotoIon[j] =
                &((*this_all_rates)->BensRates[ns].PhotoIon[0][j * chimesRateTables.NonEqIon->N_Auger[ns]]);

          ns += 1;
        }
    }
}

/* The following routine reads a photoionisation table into
 * the given NonEqIon structure. */
void ReadPhotoIonTables(struct globalVariables *myGlobalVars, char *photoIonTablePath, struct NonEq_Ionization *myNonEqIon,
                        int N_Elements_in_Bens_tables, double *dustG_arr, double *H2_dissocJ_arr, int current_spectrum)
{
  hid_t file_id, dataset, dataspace_id, memspace_id;
  herr_t status;
  hsize_t dims[2];
  hsize_t dims2D[2], count2D[2], offset2D[2];
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int rank, i, j, k, l, m, ns, r, included_index, N_arrayCells;
  char set_name[500];
  double *sigmaphot, *E_thresh, *epsilon, *cosmicrays;
  double *H2CO_N, *COself_N, *CO_S, *x_ion, *n_ion_HI, *n_ion_HeI, *shieldColumn;
  double H2_dissocJ, dust_G_parameter;
  float *shieldFac_1D, *shieldFac_2D;
  dims[1] = 1;

  file_id = H5Fopen(photoIonTablePath, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(current_spectrum == 0)
    {
      /* First time this routine is called, read in tables and
       * parameters that are only needed once, not for each
       * UV spectrum. */
      int N_Auger[N_Elements_in_Bens_tables];
      myNonEqIon->N_Auger = (int *)mymalloc("Chimes_N_Auger", chimesRateTables.NonEqIon->N_Elements * sizeof(int));

      sprintf(set_name, "/N_Auger");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_Auger);
      status  = H5Dclose(dataset);

      included_index = 0;
      for(i = 0; i < N_Elements_in_Bens_tables; i++)
        {
          if(i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
            {
              myNonEqIon->N_Auger[included_index] = N_Auger[i];
              included_index += 1;
            }
        }

      sprintf(set_name, "/Column_density_dimension");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, myNonEqIon->shieldingColumnDimensions);
      status  = H5Dclose(dataset);

      sprintf(set_name, "/shielding_dimensions");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, myNonEqIon->shielding_dimensions);
      status  = H5Dclose(dataset);

      sprintf(set_name, "/secondary_ionisation/bin_dimensions");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, myNonEqIon->secondary_ionisation_dims);
      status  = H5Dclose(dataset);

      myNonEqIon->shieldingColumnDensities =
          (double *)mymalloc("ChimesShieldCol", myNonEqIon->shieldingColumnDimensions[0] * sizeof(double));
      myNonEqIon->COself_shielding_N = (double *)mymalloc("ChimesCOself_N", myNonEqIon->shielding_dimensions[0] * sizeof(double));
      myNonEqIon->H2CO_shielding_N   = (double *)mymalloc("ChimesH2CO_N", myNonEqIon->shielding_dimensions[1] * sizeof(double));
      myNonEqIon->CO_shielding_S     = (double **)mymalloc("ChimesCO_S_dim0", myNonEqIon->shielding_dimensions[0] * sizeof(double *));
      myNonEqIon->CO_shielding_S[0]  = (double *)mymalloc(
           "ChimesCO_S_dim1", myNonEqIon->shielding_dimensions[0] * myNonEqIon->shielding_dimensions[1] * sizeof(double));
      for(i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
        myNonEqIon->CO_shielding_S[i] = &(myNonEqIon->CO_shielding_S[0][i * myNonEqIon->shielding_dimensions[1]]);

      myNonEqIon->x_ion_fraction = (double *)mymalloc("Chimes_x_ion", myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      myNonEqIon->n_ion_HI       = (double *)mymalloc("Chimes_n_HI", myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      myNonEqIon->n_ion_HeI      = (double *)mymalloc("Chimes_n_HeI", myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));

      shieldColumn = (double *)mymalloc("Chimes_col", myNonEqIon->shieldingColumnDimensions[0] * sizeof(double));
      COself_N     = (double *)mymalloc("Chimes_CO", myNonEqIon->shielding_dimensions[0] * sizeof(double));
      H2CO_N       = (double *)mymalloc("Chimes_H2CO", myNonEqIon->shielding_dimensions[1] * sizeof(double));
      CO_S         = (double *)mymalloc("Chimes_CO_S", myNonEqIon->shielding_dimensions[0] * sizeof(double));
      x_ion        = (double *)mymalloc("ChimesXion", myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      n_ion_HI     = (double *)mymalloc("ChimesnHI", myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      n_ion_HeI    = (double *)mymalloc("ChimesnHeI", myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));

      sprintf(set_name, "/Column_density_bins");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, shieldColumn);
      status  = H5Dclose(dataset);

      for(i = 0; i < myNonEqIon->shieldingColumnDimensions[0]; i++)
        myNonEqIon->shieldingColumnDensities[i] = log10(shieldColumn[i]);

      sprintf(set_name, "/COself_shielding_N");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, COself_N);
      status  = H5Dclose(dataset);

      for(i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
        myNonEqIon->COself_shielding_N[i] = COself_N[i];

      sprintf(set_name, "/H2CO_shielding_N");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, H2CO_N);
      status  = H5Dclose(dataset);

      for(i = 0; i < myNonEqIon->shielding_dimensions[1]; i++)
        myNonEqIon->H2CO_shielding_N[i] = H2CO_N[i];

      sprintf(set_name, "/CO_shielding_S");
      dataset      = H5Dopen(file_id, set_name);
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
      dims[0]      = myNonEqIon->shielding_dimensions[0];
      rank         = 1;
      memspace_id  = H5Screate_simple(rank, dims, NULL);

      for(j = 0; j < myNonEqIon->shielding_dimensions[1]; j++)
        {
          offset2D[0] = 0;
          offset2D[1] = j;
          count2D[0]  = myNonEqIon->shielding_dimensions[0];
          count2D[1]  = 1;
          status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

          status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, CO_S);

          for(i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
            myNonEqIon->CO_shielding_S[i][j] = CO_S[i];
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      status = H5Dclose(dataset);

      sprintf(set_name, "/secondary_ionisation/x_ion");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_ion);
      status  = H5Dclose(dataset);

      for(i = 0; i < myNonEqIon->secondary_ionisation_dims[0]; i++)
        myNonEqIon->x_ion_fraction[i] = log10(x_ion[i]);

      sprintf(set_name, "/secondary_ionisation/n_ion_HI");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_ion_HI);
      status  = H5Dclose(dataset);

      for(i = 0; i < myNonEqIon->secondary_ionisation_dims[0]; i++)
        myNonEqIon->n_ion_HI[i] = log10(n_ion_HI[i]);

      sprintf(set_name, "/secondary_ionisation/n_ion_HeI");
      dataset = H5Dopen(file_id, set_name);
      status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_ion_HeI);
      status  = H5Dclose(dataset);

      for(i = 0; i < myNonEqIon->secondary_ionisation_dims[0]; i++)
        myNonEqIon->n_ion_HeI[i] = log10(n_ion_HeI[i]);

      myfree(n_ion_HeI);
      myfree(n_ion_HI);
      myfree(x_ion);
      myfree(CO_S);
      myfree(H2CO_N);
      myfree(COself_N);
      myfree(shieldColumn);

      /* Also allocate memory for tables first time this routine is called. */
      ns = 0;
      for(l = 0; l < N_Elements_in_Bens_tables; l++)
        {
          if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
            {
              if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
                N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
              else
                N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

              myNonEqIon->NonEqRates[ns].sigmaphot =
                  (double ***)mymalloc("Chimes_sigmaphot_dim0", myGlobalVars->N_spectra * sizeof(double **));
              myNonEqIon->NonEqRates[ns].sigmaphot[0] =
                  (double **)mymalloc("Chimes_sigmaphot_dim1", myGlobalVars->N_spectra * N_arrayCells * sizeof(double *));
              myNonEqIon->NonEqRates[ns].sigmaphot[0][0] = (double *)mymalloc(
                  "Chimes_sigmaphot_dim2", myGlobalVars->N_spectra * N_arrayCells * myNonEqIon->N_Auger[ns] * sizeof(double));
              for(j = 0; j < myGlobalVars->N_spectra; j++)
                {
                  myNonEqIon->NonEqRates[ns].sigmaphot[j] = &(myNonEqIon->NonEqRates[ns].sigmaphot[0][j * N_arrayCells]);
                  for(k = 0; k < N_arrayCells; k++)
                    myNonEqIon->NonEqRates[ns].sigmaphot[j][k] =
                        &(myNonEqIon->NonEqRates[ns].sigmaphot[0][0][((j * N_arrayCells) + k) * myNonEqIon->N_Auger[ns]]);
                }

              myNonEqIon->NonEqRates[ns].shieldFactor1D =
                  (float ****)mymalloc("Chimes_shield1D_dim0", myGlobalVars->N_spectra * sizeof(float ***));
              myNonEqIon->NonEqRates[ns].shieldFactor1D[0] =
                  (float ***)mymalloc("Chimes_shield1D_dim1", myGlobalVars->N_spectra * 3 * sizeof(float **));
              myNonEqIon->NonEqRates[ns].shieldFactor1D[0][0] =
                  (float **)mymalloc("Chimes_shield1D_dim2", myGlobalVars->N_spectra * 3 * N_arrayCells * sizeof(float *));
              myNonEqIon->NonEqRates[ns].shieldFactor1D[0][0][0] =
                  (float *)mymalloc("Chimes_shield1D_dim3", myGlobalVars->N_spectra * 3 * N_arrayCells *
                                                                myNonEqIon->shieldingColumnDimensions[0] * sizeof(float));
              for(j = 0; j < myGlobalVars->N_spectra; j++)
                {
                  myNonEqIon->NonEqRates[ns].shieldFactor1D[j] = &(myNonEqIon->NonEqRates[ns].shieldFactor1D[0][j * 3]);
                  for(k = 0; k < 3; k++)
                    {
                      myNonEqIon->NonEqRates[ns].shieldFactor1D[j][k] =
                          &(myNonEqIon->NonEqRates[ns].shieldFactor1D[0][0][((j * 3) + k) * N_arrayCells]);
                      for(m = 0; m < N_arrayCells; m++)
                        myNonEqIon->NonEqRates[ns].shieldFactor1D[j][k][m] =
                            &(myNonEqIon->NonEqRates[ns].shieldFactor1D[0][0][0][((((j * 3) + k) * N_arrayCells) + m) *
                                                                                 myNonEqIon->shieldingColumnDimensions[0]]);
                    }
                }

              myNonEqIon->NonEqRates[ns].shieldFactor2D =
                  (float *****)mymalloc("Chimes_shield2D_dim0", myGlobalVars->N_spectra * sizeof(float ****));
              myNonEqIon->NonEqRates[ns].shieldFactor2D[0] =
                  (float ****)mymalloc("Chimes_shield2D_dim1", myGlobalVars->N_spectra * 6 * sizeof(float ***));
              myNonEqIon->NonEqRates[ns].shieldFactor2D[0][0] =
                  (float ***)mymalloc("Chimes_shield2D_dim2", myGlobalVars->N_spectra * 6 * N_arrayCells * sizeof(float **));
              myNonEqIon->NonEqRates[ns].shieldFactor2D[0][0][0] =
                  (float **)mymalloc("Chimes_shield2D_dim3", myGlobalVars->N_spectra * 6 * N_arrayCells *
                                                                 myNonEqIon->shieldingColumnDimensions[0] * sizeof(float *));
              myNonEqIon->NonEqRates[ns].shieldFactor2D[0][0][0][0] = (float *)mymalloc(
                  "Chimes_shield2D_dim4", myGlobalVars->N_spectra * 6 * N_arrayCells * myNonEqIon->shieldingColumnDimensions[0] *
                                              myNonEqIon->shieldingColumnDimensions[0] * sizeof(float));
              for(m = 0; m < myGlobalVars->N_spectra; m++)
                {
                  myNonEqIon->NonEqRates[ns].shieldFactor2D[m] = &(myNonEqIon->NonEqRates[ns].shieldFactor2D[0][m * 6]);
                  for(j = 0; j < 6; j++)
                    {
                      myNonEqIon->NonEqRates[ns].shieldFactor2D[m][j] =
                          &(myNonEqIon->NonEqRates[ns].shieldFactor2D[0][0][((m * 6) + j) * N_arrayCells]);
                      for(k = 0; k < N_arrayCells; k++)
                        {
                          myNonEqIon->NonEqRates[ns].shieldFactor2D[m][j][k] =
                              &(myNonEqIon->NonEqRates[ns].shieldFactor2D[0][0][0][((((m * 6) + j) * N_arrayCells) + k) *
                                                                                   myNonEqIon->shieldingColumnDimensions[0]]);
                          for(r = 0; r < myNonEqIon->shieldingColumnDimensions[0]; r++)
                            myNonEqIon->NonEqRates[ns].shieldFactor2D[m][j][k][r] =
                                &(myNonEqIon->NonEqRates[ns].shieldFactor2D[0][0][0][0][((((((m * 6) + j) * N_arrayCells) + k) *
                                                                                          myNonEqIon->shieldingColumnDimensions[0]) +
                                                                                         r) *
                                                                                        myNonEqIon->shieldingColumnDimensions[0]]);
                        }
                    }
                }

              myNonEqIon->NonEqRates[ns].epsilon = (double **)mymalloc("Chimes_eps_dim0", myGlobalVars->N_spectra * sizeof(double *));
              myNonEqIon->NonEqRates[ns].epsilon[0] =
                  (double *)mymalloc("Chimes_eps_dim1", myGlobalVars->N_spectra * N_arrayCells * sizeof(double));
              for(j = 0; j < myGlobalVars->N_spectra; j++)
                myNonEqIon->NonEqRates[ns].epsilon[j] = &(myNonEqIon->NonEqRates[ns].epsilon[0][j * N_arrayCells]);

              /* Read in E_thresh and cosmicrays from table. */
              E_thresh = (double *)mymalloc("Chimes_E_thresh", N_arrayCells * sizeof(double));
              cosmicrays =
                  (double *)mymalloc("Chimes_cosmic", (chimesRateTables.NonEqIon->N_Ions[ns] - 1) *
                                                          sizeof(double)); /* For hydrogen, H- & H2 cr rates are stored elsewhere */

              sprintf(set_name, "/%s/E_thresh", chimesRateTables.NonEqIon->ElementName[ns]);
              dataset = H5Dopen(file_id, set_name);
              status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, E_thresh);
              status  = H5Dclose(dataset);

              sprintf(set_name, "/%s/cosmicRays", chimesRateTables.NonEqIon->ElementName[ns]);
              dataset = H5Dopen(file_id, set_name);
              status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cosmicrays);
              status  = H5Dclose(dataset);

              for(i = 0; i < N_arrayCells; i++)
                {
                  myNonEqIon->NonEqRates[ns].E_thresh[i]   = E_thresh[i];
                  myNonEqIon->NonEqRates[ns].cosmicRays[i] = cosmicrays[i];
                }
              myfree(cosmicrays);
              myfree(E_thresh);

              ns += 1;
            }
        }
    }  // if current_spectrum == 0

  // The following need to be read in for each UV spectrum.
  ns = 0;
  for(l = 0; l < N_Elements_in_Bens_tables; l++)
    {
      if(l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
        {
          if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
            N_arrayCells = 3; /* Hydrogen photoionisation also includes H- & H2*/
          else
            N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

          sigmaphot    = (double *)mymalloc("Chimes_sigma", N_arrayCells * sizeof(double));
          epsilon      = (double *)mymalloc("ChimesEps", N_arrayCells * sizeof(double));
          shieldFac_1D = (float *)mymalloc("Chimes_shield1", 3 * sizeof(float));
          shieldFac_2D = (float *)mymalloc("Chimes_shield2", 6 * sizeof(float));

          sprintf(set_name, "/%s/sigmaPhot", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset = H5Dopen(file_id, set_name);

          if(ns == 0 || ns == 1)
            {
              status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigmaphot);
              for(i = 0; i < N_arrayCells; i++)
                myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][0] = sigmaphot[i];
              status = H5Dclose(dataset);
            }
          else
            {
              dataspace_id = H5Dget_space(dataset);
              status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
              dims[0]      = N_arrayCells;
              rank         = 1;
              memspace_id  = H5Screate_simple(rank, dims, NULL);

              for(j = 0; j < myNonEqIon->N_Auger[ns]; j++)
                {
                  offset2D[0] = 0;
                  offset2D[1] = j;
                  count2D[0]  = N_arrayCells;
                  count2D[1]  = 1;
                  status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

                  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, sigmaphot);

                  for(i = 0; i < N_arrayCells; i++)
                    myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][j] = sigmaphot[i];
                }
              H5Sclose(memspace_id);
              H5Sclose(dataspace_id);
              status = H5Dclose(dataset);
            }

          sprintf(set_name, "/%s/epsilonPhot", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset = H5Dopen(file_id, set_name);
          status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, epsilon);
          status  = H5Dclose(dataset);

          sprintf(set_name, "/%s/shieldFactor_1D", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset      = H5Dopen(file_id, set_name);
          dataspace_id = H5Dget_space(dataset);
          status       = H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
          dims[0]      = 3;
          rank         = 1;
          memspace_id  = H5Screate_simple(rank, dims, NULL);

          for(i = 0; i < N_arrayCells; i++)
            {
              for(j = 0; j < myNonEqIon->shieldingColumnDimensions[0]; j++)
                {
                  offset3D[0] = 0;
                  offset3D[1] = i;
                  offset3D[2] = j;
                  count3D[0]  = 3;
                  count3D[1]  = 1;
                  count3D[2]  = 1;
                  status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D, NULL);
                  status      = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_1D);

                  for(k = 0; k < 3; k++)
                    myNonEqIon->NonEqRates[ns].shieldFactor1D[current_spectrum][k][i][j] = (float)max(shieldFac_1D[k], -300.0);
                }
            }
          H5Sclose(memspace_id);
          H5Sclose(dataspace_id);
          H5Dclose(dataset);

          sprintf(set_name, "/%s/shieldFactor_2D", chimesRateTables.NonEqIon->ElementName[ns]);
          dataset      = H5Dopen(file_id, set_name);
          dataspace_id = H5Dget_space(dataset);
          status       = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
          dims[0]      = 6;
          rank         = 1;
          memspace_id  = H5Screate_simple(rank, dims, NULL);

          for(m = 0; m < N_arrayCells; m++)
            {
              for(i = 0; i < myNonEqIon->shieldingColumnDimensions[0]; i++)
                {
                  for(j = 0; j < myNonEqIon->shieldingColumnDimensions[0]; j++)
                    {
                      offset4D[0] = 0;
                      offset4D[1] = m;
                      offset4D[2] = i;
                      offset4D[3] = j;
                      count4D[0]  = 6;
                      count4D[1]  = 1;
                      count4D[2]  = 1;
                      count4D[3]  = 1;
                      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);
                      status      = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_2D);

                      for(k = 0; k < 6; k++)
                        myNonEqIon->NonEqRates[ns].shieldFactor2D[current_spectrum][k][m][i][j] = (float)max(shieldFac_2D[k], -300.0);
                    }
                }
            }
          H5Sclose(memspace_id);
          H5Sclose(dataspace_id);
          H5Dclose(dataset);

          for(i = 0; i < N_arrayCells; i++)
            {
              myNonEqIon->NonEqRates[ns].epsilon[current_spectrum][i] = epsilon[i];
            }

          myfree(shieldFac_2D);
          myfree(shieldFac_1D);
          myfree(epsilon);
          myfree(sigmaphot);

          ns += 1;
        }
    }

  /* Read in H2_dissocJ and dust_G_parameter for each UV spectrum. */
  sprintf(set_name, "/H2_dissocJ");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &H2_dissocJ);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/dust_G_parameter");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dust_G_parameter);
  status  = H5Dclose(dataset);

  H2_dissocJ_arr[current_spectrum] = H2_dissocJ;
  dustG_arr[current_spectrum]      = dust_G_parameter;

  status = H5Fclose(file_id);
}

void GetNonEqTables(int ns, struct globalVariables *myGlobalVars)
{
  hid_t file_id, dataset, dataspace_id, memspace_id;
  herr_t status;
  int rank, i, j, k, l;
  hsize_t dims[2];
  hsize_t dims2D[2], offset2D[2], count2D[2];
  hsize_t dims4D[4], offset4D[4], count4D[4];
  dims[1] = 1;

  double cooling_array_1D[chimesRateTables.nei_cooling_table_dimensions[0]];
  double chianti_cooling_array_1D[chimesRateTables.chianti_cooling_table_dimensions[0]];

  char fname[500], set_name[500];

  double alpharad[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  double alphadi[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  double betacoll[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  double cool[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];

  double CTHrecof[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  double CTHionof[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  double CTHerecof[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  double CTHeionof[chimesRateTables.NonEqIon->N_Ions[ns]][chimesRateTables.NonEqIon->N_Temperatures];
  int CTHion_mask[chimesRateTables.NonEqIon->N_Ions[ns]];
  int CTHrec_mask[chimesRateTables.NonEqIon->N_Ions[ns]];
  int CTHeion_mask[chimesRateTables.NonEqIon->N_Ions[ns]];
  int CTHerec_mask[chimesRateTables.NonEqIon->N_Ions[ns]];

  sprintf(fname, "%s%s", myGlobalVars->BenTablesPath, chimesRateTables.NonEqIon->FileName[ns]);

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  sprintf(set_name, "/Rates/AlphaRad");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, alpharad);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/AlphaDi");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, alphadi);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/BetaColl");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, betacoll);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/Cool");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cool);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/CTHrecof");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHrecof);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/CTHionof");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHionof);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/CTHerecof");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHerecof);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Rates/CTHeionof");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHeionof);
  status  = H5Dclose(dataset);

  /* CT masks */
  sprintf(set_name, "/Header/CTHion_mask");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHion_mask);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Header/CTHrec_mask");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHrec_mask);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Header/CTHeion_mask");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHeion_mask);
  status  = H5Dclose(dataset);

  sprintf(set_name, "/Header/CTHerec_mask");
  dataset = H5Dopen(file_id, set_name);
  status  = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CTHerec_mask);
  status  = H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.NonEqIon->N_Ions[ns]; i++)
    {
      chimesRateTables.NonEqIon->NonEqRates[ns].CTHion_mask[i]  = CTHion_mask[i];
      chimesRateTables.NonEqIon->NonEqRates[ns].CTHrec_mask[i]  = CTHrec_mask[i];
      chimesRateTables.NonEqIon->NonEqRates[ns].CTHeion_mask[i] = CTHeion_mask[i];
      chimesRateTables.NonEqIon->NonEqRates[ns].CTHerec_mask[i] = CTHerec_mask[i];
      for(j = 0; j < chimesRateTables.NonEqIon->N_Temperatures; j++)
        {
          /* Note that we want to interpolate all of the rates
           * linearly in LOG space, hence take logs here */
          chimesRateTables.NonEqIon->NonEqRates[ns].alpharad[i][j] = log10(max(alpharad[i][j], 1.0e-300));
          chimesRateTables.NonEqIon->NonEqRates[ns].alphadi[i][j]  = log10(max(alphadi[i][j], 1.0e-300));
          chimesRateTables.NonEqIon->NonEqRates[ns].betacoll[i][j] = log10(max(betacoll[i][j], 1.0e-300));
          chimesRateTables.NonEqIon->NonEqRates[ns].cool[i][j] = log10(max(cool[i][j], 1.0e-300)); /* COOLING RATES ARE POSITIVE!!!*/

          chimesRateTables.NonEqIon->NonEqRates[ns].CTHrecof[i][j]  = log10(max(CTHrecof[i][j], 1.0e-300));
          chimesRateTables.NonEqIon->NonEqRates[ns].CTHerecof[i][j] = log10(max(CTHerecof[i][j], 1.0e-300));

          chimesRateTables.NonEqIon->NonEqRates[ns].CTHionof[i][j]  = log10(max(CTHionof[i][j], 1.0e-300));
          chimesRateTables.NonEqIon->NonEqRates[ns].CTHeionof[i][j] = log10(max(CTHeionof[i][j], 1.0e-300));
        }
    }

  if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == CI)
    {
      /* CI nei cooling */
      dataset      = H5Dopen(file_id, "/Rates/CI_cool");
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);

      dims[0]     = chimesRateTables.nei_cooling_table_dimensions[0];
      rank        = 1;
      memspace_id = H5Screate_simple(rank, dims, NULL);

      for(i = 0; i < chimesRateTables.nei_cooling_table_dimensions[1]; i++)
        {
          for(j = 0; j < chimesRateTables.nei_cooling_table_dimensions[2]; j++)
            {
              for(k = 0; k < chimesRateTables.nei_cooling_table_dimensions[3]; k++)
                {
                  offset4D[0] = 0;
                  offset4D[1] = i;
                  offset4D[2] = j;
                  offset4D[3] = k;
                  count4D[0]  = chimesRateTables.nei_cooling_table_dimensions[0];
                  count4D[1]  = 1;
                  count4D[2]  = 1;
                  count4D[3]  = 1;
                  status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);

                  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, cooling_array_1D);

                  for(l = 0; l < chimesRateTables.nei_cooling_table_dimensions[0]; l++)
                    chimesRateTables.nei_cooling_CI[l][i][j][k] = cooling_array_1D[l];
                }
            }
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset);

      /* CII chianti cooling */
      dataset      = H5Dopen(file_id, "/Rates/CII_cool");
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

      dims[0]     = chimesRateTables.chianti_cooling_table_dimensions[0];
      rank        = 1;
      memspace_id = H5Screate_simple(rank, dims, NULL);

      for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[1]; i++)
        {
          offset2D[0] = 0;
          offset2D[1] = i;
          count2D[0]  = chimesRateTables.chianti_cooling_table_dimensions[0];
          count2D[1]  = 1;
          status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);
          status      = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, chianti_cooling_array_1D);
          for(l = 0; l < chimesRateTables.chianti_cooling_table_dimensions[0]; l++)
            chimesRateTables.chianti_cooling_CII[l][i] = chianti_cooling_array_1D[l];
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset);
    }

  else if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == NI)
    {
      /* NII chianti cooling */
      dataset      = H5Dopen(file_id, "/Rates/NII_cool");
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

      dims[0]     = chimesRateTables.chianti_cooling_table_dimensions[0];
      rank        = 1;
      memspace_id = H5Screate_simple(rank, dims, NULL);

      for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[1]; i++)
        {
          offset2D[0] = 0;
          offset2D[1] = i;
          count2D[0]  = chimesRateTables.chianti_cooling_table_dimensions[0];
          count2D[1]  = 1;
          status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);
          status      = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, chianti_cooling_array_1D);
          for(l = 0; l < chimesRateTables.chianti_cooling_table_dimensions[0]; l++)
            chimesRateTables.chianti_cooling_NII[l][i] = chianti_cooling_array_1D[l];
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset);
    }

  else if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == OI)
    {
      dataset      = H5Dopen(file_id, "/Rates/OI_cool");
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);

      dims[0]     = chimesRateTables.nei_cooling_table_dimensions[0];
      rank        = 1;
      memspace_id = H5Screate_simple(rank, dims, NULL);

      for(i = 0; i < chimesRateTables.nei_cooling_table_dimensions[1]; i++)
        {
          for(j = 0; j < chimesRateTables.nei_cooling_table_dimensions[2]; j++)
            {
              for(k = 0; k < chimesRateTables.nei_cooling_table_dimensions[3]; k++)
                {
                  offset4D[0] = 0;
                  offset4D[1] = i;
                  offset4D[2] = j;
                  offset4D[3] = k;
                  count4D[0]  = chimesRateTables.nei_cooling_table_dimensions[0];
                  count4D[1]  = 1;
                  count4D[2]  = 1;
                  count4D[3]  = 1;
                  status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);

                  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, cooling_array_1D);

                  for(l = 0; l < chimesRateTables.nei_cooling_table_dimensions[0]; l++)
                    chimesRateTables.nei_cooling_OI[l][i][j][k] = cooling_array_1D[l];
                }
            }
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset);
    }

  else if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == SiI)
    {
      /* SiII chianti cooling */
      dataset      = H5Dopen(file_id, "/Rates/SiII_cool");
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

      dims[0]     = chimesRateTables.chianti_cooling_table_dimensions[0];
      rank        = 1;
      memspace_id = H5Screate_simple(rank, dims, NULL);

      for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[1]; i++)
        {
          offset2D[0] = 0;
          offset2D[1] = i;
          count2D[0]  = chimesRateTables.chianti_cooling_table_dimensions[0];
          count2D[1]  = 1;
          status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);
          status      = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, chianti_cooling_array_1D);
          for(l = 0; l < chimesRateTables.chianti_cooling_table_dimensions[0]; l++)
            chimesRateTables.chianti_cooling_SiII[l][i] = chianti_cooling_array_1D[l];
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset);
    }

  else if(chimesRateTables.NonEqIon->IonIndexBegin[ns] == FeI)
    {
      /* FeII chianti cooling */
      dataset      = H5Dopen(file_id, "/Rates/FeII_cool");
      dataspace_id = H5Dget_space(dataset);
      status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

      dims[0]     = chimesRateTables.chianti_cooling_table_dimensions[0];
      rank        = 1;
      memspace_id = H5Screate_simple(rank, dims, NULL);

      for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[1]; i++)
        {
          offset2D[0] = 0;
          offset2D[1] = i;
          count2D[0]  = chimesRateTables.chianti_cooling_table_dimensions[0];
          count2D[1]  = 1;
          status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);
          status      = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, chianti_cooling_array_1D);
          for(l = 0; l < chimesRateTables.chianti_cooling_table_dimensions[0]; l++)
            chimesRateTables.chianti_cooling_FeII[l][i] = chianti_cooling_array_1D[l];
        }
      H5Sclose(memspace_id);
      H5Sclose(dataspace_id);
      H5Dclose(dataset);
    }

  status = H5Fclose(file_id);
}

void initialise_bens_tables(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates,
                            double *dustG_arr, double *H2_dissocJ_arr)
{
  int i, j, k;

  hid_t file_id, dataset, datatype;

  hsize_t el_name_length   = EL_NAME_LENGTH;
  hsize_t file_name_length = NONEQION_NAME_LENGTH;
  int N_Elements_in_Bens_tables;
  int included_index;
  char fname[256];

  chimesRateTables.NonEqIon = (struct NonEq_Ionization *)mymalloc("Chimes_NonEqIon", sizeof(struct NonEq_Ionization));

  sprintf(
      fname, "%sIonRates_Head.HM01.hdf5",
      myGlobalVars->BenTablesPath); /* Note that we ignore everything that depends on the UV spectrum, so we can pick any of them */

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    printf("[read_noneq_tables()]: unable to open file %s\n", fname);

  dataset = H5Dopen(file_id, "N_Elements");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_Elements_in_Bens_tables);
  H5Dclose(dataset);

#ifdef CHIMES_REDSHIFT_DEPENDENT_UVB
  Chimes_N_Elements_in_Bens_tables = N_Elements_in_Bens_tables;
#endif

  dataset = H5Dopen(file_id, "N_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.NonEqIon->N_Temperatures);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "/nei_headers/BinSizes");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.nei_cooling_table_dimensions[0]);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "/chianti_headers/BinSizes");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.chianti_cooling_table_dimensions[0]);
  H5Dclose(dataset);

  /* First we read in the element names and filenames to
   * temporary arrays. */

  char element_names[N_Elements_in_Bens_tables][EL_NAME_LENGTH];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, el_name_length);
  dataset = H5Dopen(file_id, "Element_Names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  H5Dclose(dataset);
  H5Tclose(datatype);

  char file_names[N_Elements_in_Bens_tables][NONEQION_NAME_LENGTH];

  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, file_name_length);
  dataset = H5Dopen(file_id, "File_Names");
  H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, file_names);
  H5Dclose(dataset);
  H5Tclose(datatype);

  /* Now we need to determine how many elements are actually
   * included in our network. */
  chimesRateTables.NonEqIon->N_Elements = 2; /* H & He */
  for(i = 0; i < 9; i++)
    if(myGlobalVars->element_included[i] == 1)
      chimesRateTables.NonEqIon->N_Elements += 1;

  /* Finally, set the element names and file names
   * in the arrays for bens tables. */
  chimesRateTables.NonEqIon->ElementName =
      (char **)mymalloc("Chimes_ElName_dim0", chimesRateTables.NonEqIon->N_Elements * sizeof(char *));
  chimesRateTables.NonEqIon->ElementName[0] =
      (char *)mymalloc("Chimes_ElName_dim1", chimesRateTables.NonEqIon->N_Elements * EL_NAME_LENGTH * sizeof(char));
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    chimesRateTables.NonEqIon->ElementName[i] = &(chimesRateTables.NonEqIon->ElementName[0][i * EL_NAME_LENGTH]);

  included_index = 0;
  for(i = 0; i < N_Elements_in_Bens_tables; i++)
    {
      if(i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
        {
          strcpy(chimesRateTables.NonEqIon->ElementName[included_index], element_names[i]);
          included_index += 1;
        }
    }

  chimesRateTables.NonEqIon->FileName =
      (char **)mymalloc("Chimes_filename_dim0", chimesRateTables.NonEqIon->N_Elements * sizeof(char *));
  chimesRateTables.NonEqIon->FileName[0] =
      (char *)mymalloc("Chimes_filename_dim1", chimesRateTables.NonEqIon->N_Elements * NONEQION_NAME_LENGTH * sizeof(char));
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    chimesRateTables.NonEqIon->FileName[i] = &(chimesRateTables.NonEqIon->FileName[0][i * NONEQION_NAME_LENGTH]);

  included_index = 0;
  for(i = 0; i < N_Elements_in_Bens_tables; i++)
    {
      if(i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
        {
          strcpy(chimesRateTables.NonEqIon->FileName[included_index], file_names[i]);
          included_index += 1;
        }
    }

  chimesRateTables.NonEqIon->Temperatures = (double *)mymalloc("Chimes_T", chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));

  dataset = H5Dopen(file_id, "Temperatures");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.NonEqIon->Temperatures[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.NonEqIon->N_Temperatures; i++)
    chimesRateTables.NonEqIon->Temperatures[i] = log10(chimesRateTables.NonEqIon->Temperatures[i]);

  /* From here and below, first read in tables to temporary arrays,
   * and THEN step through them selecting the elements that
   * are actually present in the network. */

  int atomic_number[N_Elements_in_Bens_tables];
  int N_Ions[N_Elements_in_Bens_tables];
  double atomic_weights[N_Elements_in_Bens_tables];

  dataset = H5Dopen(file_id, "AtomicNumbers");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_number);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "N_Ions");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_Ions);
  H5Dclose(dataset);

  dataset = H5Dopen(file_id, "AtomicWeights");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_weights);
  H5Dclose(dataset);

  chimesRateTables.NonEqIon->AtomicNumber = (int *)mymalloc("Chimes_atomNo", chimesRateTables.NonEqIon->N_Elements * sizeof(int));
  chimesRateTables.NonEqIon->N_Ions       = (int *)mymalloc("Chimes_N_Ions", chimesRateTables.NonEqIon->N_Elements * sizeof(int));
  chimesRateTables.NonEqIon->AtomicWeights =
      (double *)mymalloc("Chimes_atomWeight", chimesRateTables.NonEqIon->N_Elements * sizeof(double));

  included_index = 0;
  for(i = 0; i < N_Elements_in_Bens_tables; i++)
    {
      if(i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
        {
          chimesRateTables.NonEqIon->AtomicNumber[included_index]  = atomic_number[i];
          chimesRateTables.NonEqIon->N_Ions[included_index]        = N_Ions[i];
          chimesRateTables.NonEqIon->AtomicWeights[included_index] = atomic_weights[i];
          included_index += 1;
        }
    }

  chimesRateTables.nei_cooling_temperature =
      (double *)mymalloc("Chimes_nei_T", chimesRateTables.nei_cooling_table_dimensions[0] * sizeof(double));

  dataset = H5Dopen(file_id, "/nei_headers/Temperatures");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_temperature);
  H5Dclose(dataset);

  chimesRateTables.nei_cooling_HIAbundance =
      (double *)mymalloc("Chimes_nei_HI", chimesRateTables.nei_cooling_table_dimensions[1] * sizeof(double));

  dataset = H5Dopen(file_id, "/nei_headers/HIDensities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_HIAbundance);
  H5Dclose(dataset);

  chimesRateTables.nei_cooling_ElectronAbundance =
      (double *)mymalloc("Chimes_nei_elec", chimesRateTables.nei_cooling_table_dimensions[2] * sizeof(double));

  dataset = H5Dopen(file_id, "/nei_headers/ElectronDensities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_ElectronAbundance);
  H5Dclose(dataset);

  chimesRateTables.nei_cooling_HIIAbundance =
      (double *)mymalloc("Chimes_nei_HII", chimesRateTables.nei_cooling_table_dimensions[3] * sizeof(double));

  dataset = H5Dopen(file_id, "/nei_headers/HIIDensities");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_HIIAbundance);
  H5Dclose(dataset);

  chimesRateTables.chianti_cooling_temperature =
      (double *)mymalloc("Chimes_chianti_T", chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double));

  dataset = H5Dopen(file_id, "/chianti_headers/Temperatures");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.chianti_cooling_temperature);
  H5Dclose(dataset);

  chimesRateTables.chianti_cooling_ElectronDensity =
      (double *)mymalloc("Chimes_chianti_elec", chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));

  dataset = H5Dopen(file_id, "/chianti_headers/ElectronDensity");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.chianti_cooling_ElectronDensity);
  H5Dclose(dataset);

  chimesRateTables.NonEqIon->NonEqRates =
      (struct NonEq_Rates *)mymalloc("Chimes_NonEqRates", chimesRateTables.NonEqIon->N_Elements * sizeof(struct NonEq_Rates));

  H5Fclose(file_id);

  chimesRateTables.nei_cooling_CI =
      (double ****)mymalloc("Chimes_nei_CI_dim0", chimesRateTables.nei_cooling_table_dimensions[0] * sizeof(double ***));
  chimesRateTables.nei_cooling_CI[0] =
      (double ***)mymalloc("Chimes_nei_CI_dim1", chimesRateTables.nei_cooling_table_dimensions[0] *
                                                     chimesRateTables.nei_cooling_table_dimensions[1] * sizeof(double **));
  chimesRateTables.nei_cooling_CI[0][0] = (double **)mymalloc(
      "Chimes_nei_CI_dim2", chimesRateTables.nei_cooling_table_dimensions[0] * chimesRateTables.nei_cooling_table_dimensions[1] *
                                chimesRateTables.nei_cooling_table_dimensions[2] * sizeof(double *));
  chimesRateTables.nei_cooling_CI[0][0][0] = (double *)mymalloc(
      "Chimes_nei_CI_dim3", chimesRateTables.nei_cooling_table_dimensions[0] * chimesRateTables.nei_cooling_table_dimensions[1] *
                                chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[3] *
                                sizeof(double));

  for(i = 0; i < chimesRateTables.nei_cooling_table_dimensions[0]; i++)
    {
      chimesRateTables.nei_cooling_CI[i] = &(chimesRateTables.nei_cooling_CI[0][i * chimesRateTables.nei_cooling_table_dimensions[1]]);
      for(j = 0; j < chimesRateTables.nei_cooling_table_dimensions[1]; j++)
        {
          chimesRateTables.nei_cooling_CI[i][j] =
              &(chimesRateTables.nei_cooling_CI[0][0][((i * chimesRateTables.nei_cooling_table_dimensions[1]) + j) *
                                                      chimesRateTables.nei_cooling_table_dimensions[2]]);
          for(k = 0; k < chimesRateTables.nei_cooling_table_dimensions[2]; k++)
            chimesRateTables.nei_cooling_CI[i][j][k] =
                &(chimesRateTables.nei_cooling_CI[0][0][0][((((i * chimesRateTables.nei_cooling_table_dimensions[1]) + j) *
                                                             chimesRateTables.nei_cooling_table_dimensions[2]) +
                                                            k) *
                                                           chimesRateTables.nei_cooling_table_dimensions[3]]);
        }
    }

  chimesRateTables.nei_cooling_OI =
      (double ****)mymalloc("Chimes_nei_OI_dim0", chimesRateTables.nei_cooling_table_dimensions[0] * sizeof(double ***));
  chimesRateTables.nei_cooling_OI[0] =
      (double ***)mymalloc("Chimes_nei_OI_dim1", chimesRateTables.nei_cooling_table_dimensions[0] *
                                                     chimesRateTables.nei_cooling_table_dimensions[1] * sizeof(double **));
  chimesRateTables.nei_cooling_OI[0][0] = (double **)mymalloc(
      "Chimes_nei_OI_dim2", chimesRateTables.nei_cooling_table_dimensions[0] * chimesRateTables.nei_cooling_table_dimensions[1] *
                                chimesRateTables.nei_cooling_table_dimensions[2] * sizeof(double *));
  chimesRateTables.nei_cooling_OI[0][0][0] = (double *)mymalloc(
      "Chimes_nei_OI_dim3", chimesRateTables.nei_cooling_table_dimensions[0] * chimesRateTables.nei_cooling_table_dimensions[1] *
                                chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[3] *
                                sizeof(double));

  for(i = 0; i < chimesRateTables.nei_cooling_table_dimensions[0]; i++)
    {
      chimesRateTables.nei_cooling_OI[i] = &(chimesRateTables.nei_cooling_OI[0][i * chimesRateTables.nei_cooling_table_dimensions[1]]);
      for(j = 0; j < chimesRateTables.nei_cooling_table_dimensions[1]; j++)
        {
          chimesRateTables.nei_cooling_OI[i][j] =
              &(chimesRateTables.nei_cooling_OI[0][0][((i * chimesRateTables.nei_cooling_table_dimensions[1]) + j) *
                                                      chimesRateTables.nei_cooling_table_dimensions[2]]);
          for(k = 0; k < chimesRateTables.nei_cooling_table_dimensions[2]; k++)
            chimesRateTables.nei_cooling_OI[i][j][k] =
                &(chimesRateTables.nei_cooling_OI[0][0][0][((((i * chimesRateTables.nei_cooling_table_dimensions[1]) + j) *
                                                             chimesRateTables.nei_cooling_table_dimensions[2]) +
                                                            k) *
                                                           chimesRateTables.nei_cooling_table_dimensions[3]]);
        }
    }

  chimesRateTables.chianti_cooling_CII =
      (double **)mymalloc("Chimes_chianti_CII_dim0", chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.chianti_cooling_CII[0] =
      (double *)mymalloc("Chimes_chianti_CII_dim1", chimesRateTables.chianti_cooling_table_dimensions[0] *
                                                        chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_CII[i] =
        &(chimesRateTables.chianti_cooling_CII[0][i * chimesRateTables.chianti_cooling_table_dimensions[1]]);

  chimesRateTables.chianti_cooling_NII =
      (double **)mymalloc("Chimes_chianti_NII_dim0", chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.chianti_cooling_NII[0] =
      (double *)mymalloc("Chimes_chianti_NII_dim1", chimesRateTables.chianti_cooling_table_dimensions[0] *
                                                        chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_NII[i] =
        &(chimesRateTables.chianti_cooling_NII[0][i * chimesRateTables.chianti_cooling_table_dimensions[1]]);

  chimesRateTables.chianti_cooling_SiII =
      (double **)mymalloc("Chimes_chianti_SiII_dim0", chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.chianti_cooling_SiII[0] =
      (double *)mymalloc("Chimes_chianti_SiII_dim1", chimesRateTables.chianti_cooling_table_dimensions[0] *
                                                         chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_SiII[i] =
        &(chimesRateTables.chianti_cooling_SiII[0][i * chimesRateTables.chianti_cooling_table_dimensions[1]]);

  chimesRateTables.chianti_cooling_FeII =
      (double **)mymalloc("Chimes_chianti_FeII_dim0", chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.chianti_cooling_FeII[0] =
      (double *)mymalloc("Chimes_chianti_FeII_dim1", chimesRateTables.chianti_cooling_table_dimensions[0] *
                                                         chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_FeII[i] =
        &(chimesRateTables.chianti_cooling_FeII[0][i * chimesRateTables.chianti_cooling_table_dimensions[1]]);

  *this_all_rates = (struct All_rate_variables_structure *)mymalloc("Chimes_allRates", sizeof(struct All_rate_variables_structure));

  (*this_all_rates)->BensRates = (struct Bens_rate_structure *)mymalloc(
      "Chimes_bensRates", chimesRateTables.NonEqIon->N_Elements * sizeof(struct Bens_rate_structure));

  chimesRateTables.NonEqIon->IonIndexBegin = (int *)mymalloc("Chimes_ionIndex", chimesRateTables.NonEqIon->N_Elements * sizeof(int));

  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Hydrogen") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = HI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Helium") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = HeI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Carbon") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = CI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Nitrogen") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = NI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Oxygen") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = OI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Neon") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = NeI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Magnesium") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = MgI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Silicon") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = SiI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Sulphur") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = SI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Calcium") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = CaI;
      else if(strcmp(chimesRateTables.NonEqIon->ElementName[i], "Iron") == 0)
        chimesRateTables.NonEqIon->IonIndexBegin[i] = FeI;

      chimesRateTables.NonEqIon->NonEqRates[i].N_Ions = chimesRateTables.NonEqIon->N_Ions[i];

      chimesRateTables.NonEqIon->NonEqRates[i].alpharad =
          (double **)mymalloc("Chimes_alpharad_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].alpharad[0] = (double *)mymalloc(
          "Chimes_alpharad_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].alpharad[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].alpharad[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].alphadi =
          (double **)mymalloc("Chimes_alphadi_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].alphadi[0] = (double *)mymalloc(
          "Chimes_alphadi_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].alphadi[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].alphadi[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].betacoll =
          (double **)mymalloc("Chimes_betacoll_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].betacoll[0] = (double *)mymalloc(
          "Chimes_betacoll_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].betacoll[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].betacoll[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].cool =
          (double **)mymalloc("Chimes_cool_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].cool[0] = (double *)mymalloc(
          "Chimes_cool_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].cool[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].cool[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].E_thresh =
          (double *)mymalloc("Chimes_Ethresh", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));

      chimesRateTables.NonEqIon->NonEqRates[i].cosmicRays =
          (double *)mymalloc("Chimes_cosmic", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));

      chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof =
          (double **)mymalloc("Chimes_CTHrecof_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[0] = (double *)mymalloc(
          "Chimes_CTHrecof_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].CTHionof =
          (double **)mymalloc("Chimes_CTHionof_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[0] = (double *)mymalloc(
          "Chimes_CTHionof_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof =
          (double **)mymalloc("Chimes_CTHerecof_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[0] = (double *)mymalloc(
          "Chimes_CTHerecof_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof =
          (double **)mymalloc("Chimes_CTHeionof_dim0", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(double *));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[0] = (double *)mymalloc(
          "Chimes_CTHeionof_dim1", chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
        chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[j] =
            &(chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[0][j * chimesRateTables.NonEqIon->N_Temperatures]);

      chimesRateTables.NonEqIon->NonEqRates[i].CTHion_mask =
          (int *)mymalloc("Chimes_CTHion_mask", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHrec_mask =
          (int *)mymalloc("Chimes_CTHrec_mask", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHeion_mask =
          (int *)mymalloc("Chimes_CTHeion_mask", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHerec_mask =
          (int *)mymalloc("Chimes_CTHerec_mask", chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));

      GetNonEqTables(i, myGlobalVars);

      (*this_all_rates)->BensRates[i].CollisIon =
          (double *)mymalloc("Chimes_collis", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].Recomb =
          (double *)mymalloc("Chimes_rec", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHrec =
          (double *)mymalloc("Chimes_CTH_rec", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHion =
          (double *)mymalloc("Chimes_CTH_ion", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHerec =
          (double *)mymalloc("Chimes_CTHe_rec", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHeion =
          (double *)mymalloc("Chimes_CTHe_ion", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].cosmicRays =
          (double *)mymalloc("Chimes_CR", (chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
    }
  GetPhotoIonTables(myGlobalVars, N_Elements_in_Bens_tables, this_all_rates, dustG_arr, H2_dissocJ_arr);
}

void initialise_additional_rates_tables(struct globalVariables *myGlobalVars)
{
  int i;
  hid_t file_id, dataset;
  char fname[256];

  sprintf(fname, "%s", myGlobalVars->AdditionalRatesTablesPath);

  chimesRateTables.RatesTables = (struct Additional_Rates *)mymalloc("Chimes_RatesTables", sizeof(struct Additional_Rates));

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0)
    printf("[read_noneq_tables()]: unable to open file %s\n", fname);

  dataset = H5Dopen(file_id, "N_Temperatures");
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->N_Temperatures);
  H5Dclose(dataset);

  /* Allocate memory to all of the tables
   * in the additional rates structure. */
  chimesRateTables.RatesTables->Temperatures =
      (double *)mymalloc("Chimes_T", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction1 =
      (double *)mymalloc("Chimes_r1", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction2 =
      (double *)mymalloc("Chimes_r2", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction3 =
      (double *)mymalloc("Chimes_r3", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction5 =
      (double *)mymalloc("Chimes_r5", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction6 =
      (double *)mymalloc("Chimes_r6", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction7 =
      (double *)mymalloc("Chimes_r7", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction10 =
      (double *)mymalloc("Chimes_r10", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction11 =
      (double *)mymalloc("Chimes_r11", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction13A =
      (double *)mymalloc("Chimes_r13A", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction13B =
      (double *)mymalloc("Chimes_r13B", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction15 =
      (double *)mymalloc("Chimes_r15", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction16 =
      (double *)mymalloc("Chimes_r16", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction17 =
      (double *)mymalloc("Chimes_r17", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction24 =
      (double *)mymalloc("Chimes_r24", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction25rrA =
      (double *)mymalloc("Chimes_r25rrA", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction25rrB =
      (double *)mymalloc("Chimes_r25rrB", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction25di =
      (double *)mymalloc("Chimes_r25di", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction26 =
      (double *)mymalloc("Chimes_r26", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction27 =
      (double *)mymalloc("Chimes_r27", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction83 =
      (double *)mymalloc("Chimes_r83", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction84 =
      (double *)mymalloc("Chimes_r84", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction87 =
      (double *)mymalloc("Chimes_r87", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction97 =
      (double *)mymalloc("Chimes_r97", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction101 =
      (double *)mymalloc("Chimes_r101", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction106 =
      (double *)mymalloc("Chimes_r106", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction107 =
      (double *)mymalloc("Chimes_r107", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction112 =
      (double *)mymalloc("Chimes_r112", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction113 =
      (double *)mymalloc("Chimes_r113", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction114 =
      (double *)mymalloc("Chimes_r114", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction116 =
      (double *)mymalloc("Chimes_r116", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction117 =
      (double *)mymalloc("Chimes_r117", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction120 =
      (double *)mymalloc("Chimes_r120", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction121 =
      (double *)mymalloc("Chimes_r121", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction122 =
      (double *)mymalloc("Chimes_r122", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction124 =
      (double *)mymalloc("Chimes_r124", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction128 =
      (double *)mymalloc("Chimes_r128", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction129 =
      (double *)mymalloc("Chimes_r129", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction130 =
      (double *)mymalloc("Chimes_r130", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction131 =
      (double *)mymalloc("Chimes_r131", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction133 =
      (double *)mymalloc("Chimes_r133", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction134 =
      (double *)mymalloc("Chimes_r134", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction135 =
      (double *)mymalloc("Chimes_r135", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction136 =
      (double *)mymalloc("Chimes_r136", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction137 =
      (double *)mymalloc("Chimes_r137", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction138 =
      (double *)mymalloc("Chimes_r138", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction139 =
      (double *)mymalloc("Chimes_r139", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction142 =
      (double *)mymalloc("Chimes_r142", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction147 =
      (double *)mymalloc("Chimes_r147", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction150 =
      (double *)mymalloc("Chimes_r150", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction157 =
      (double *)mymalloc("Chimes_r157", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction186 =
      (double *)mymalloc("Chimes_r186", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction187 =
      (double *)mymalloc("Chimes_r187", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction188 =
      (double *)mymalloc("Chimes_r188", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction189 =
      (double *)mymalloc("Chimes_r189", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction190 =
      (double *)mymalloc("Chimes_r190", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction191 =
      (double *)mymalloc("Chimes_r191", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction192 =
      (double *)mymalloc("Chimes_r192", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction193 =
      (double *)mymalloc("Chimes_r193", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction194 =
      (double *)mymalloc("Chimes_r194", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction195 =
      (double *)mymalloc("Chimes_r195", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction196 =
      (double *)mymalloc("Chimes_r196", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction197 =
      (double *)mymalloc("Chimes_r197", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction198 =
      (double *)mymalloc("Chimes_r198", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction199 =
      (double *)mymalloc("Chimes_r199", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction200 =
      (double *)mymalloc("Chimes_r200", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction201 =
      (double *)mymalloc("Chimes_r201", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction202 =
      (double *)mymalloc("Chimes_r202", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction203 =
      (double *)mymalloc("Chimes_r203", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction204 =
      (double *)mymalloc("Chimes_r204", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction205 =
      (double *)mymalloc("Chimes_r205", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction206 =
      (double *)mymalloc("Chimes_r206", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction207 =
      (double *)mymalloc("Chimes_r207", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction220 =
      (double *)mymalloc("Chimes_r220", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction221 =
      (double *)mymalloc("Chimes_r221", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction222 =
      (double *)mymalloc("Chimes_r222", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction223 =
      (double *)mymalloc("Chimes_r223", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction225 =
      (double *)mymalloc("Chimes_r225", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction226 =
      (double *)mymalloc("Chimes_r226", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction227 =
      (double *)mymalloc("Chimes_r227", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction228 =
      (double *)mymalloc("Chimes_r228", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction229 =
      (double *)mymalloc("Chimes_r229", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction230 =
      (double *)mymalloc("Chimes_r230", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction231 =
      (double *)mymalloc("Chimes_r231", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction232 =
      (double *)mymalloc("Chimes_r232", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction233 =
      (double *)mymalloc("Chimes_r233", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction234 =
      (double *)mymalloc("Chimes_r234", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction235 =
      (double *)mymalloc("Chimes_r235", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction236 =
      (double *)mymalloc("Chimes_r236", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction237 =
      (double *)mymalloc("Chimes_r237", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction238 =
      (double *)mymalloc("Chimes_r238", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));
  chimesRateTables.RatesTables->reaction283 =
      (double *)mymalloc("Chimes_r283", chimesRateTables.RatesTables->N_Temperatures * sizeof(double));

  dataset = H5Dopen(file_id, "Temperatures");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->Temperatures[0]);
  H5Dclose(dataset);
  /* Note that the temperature bins in the tables are
   * already log. */

  dataset = H5Dopen(file_id, "/Rates/k1");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction1[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction1[i] = log10(chimesRateTables.RatesTables->reaction1[i]);

  dataset = H5Dopen(file_id, "/Rates/k2");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction2[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction2[i] = log10(chimesRateTables.RatesTables->reaction2[i]);

  dataset = H5Dopen(file_id, "/Rates/k3");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction3[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction3[i] = log10(chimesRateTables.RatesTables->reaction3[i]);

  dataset = H5Dopen(file_id, "/Rates/k5");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction5[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction5[i] = log10(chimesRateTables.RatesTables->reaction5[i]);

  dataset = H5Dopen(file_id, "/Rates/k6");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction6[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction6[i] = log10(chimesRateTables.RatesTables->reaction6[i]);

  dataset = H5Dopen(file_id, "/Rates/k7");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction7[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction7[i] = log10(chimesRateTables.RatesTables->reaction7[i]);

  dataset = H5Dopen(file_id, "/Rates/k10");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction10[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction10[i] = log10(chimesRateTables.RatesTables->reaction10[i]);

  dataset = H5Dopen(file_id, "/Rates/k11");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction11[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction11[i] = log10(chimesRateTables.RatesTables->reaction11[i]);

  dataset = H5Dopen(file_id, "/Rates/k13A");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction13A[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction13A[i] = log10(chimesRateTables.RatesTables->reaction13A[i]);

  dataset = H5Dopen(file_id, "/Rates/k13B");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction13B[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction13B[i] = log10(chimesRateTables.RatesTables->reaction13B[i]);

  dataset = H5Dopen(file_id, "/Rates/k15");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction15[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction15[i] = log10(chimesRateTables.RatesTables->reaction15[i]);

  dataset = H5Dopen(file_id, "/Rates/k16");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction16[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction16[i] = log10(chimesRateTables.RatesTables->reaction16[i]);

  dataset = H5Dopen(file_id, "/Rates/k17");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction17[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction17[i] = log10(chimesRateTables.RatesTables->reaction17[i]);

  dataset = H5Dopen(file_id, "/Rates/k24");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction24[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction24[i] = log10(chimesRateTables.RatesTables->reaction24[i]);

  dataset = H5Dopen(file_id, "/Rates/k25rrA");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction25rrA[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction25rrA[i] = log10(chimesRateTables.RatesTables->reaction25rrA[i]);

  dataset = H5Dopen(file_id, "/Rates/k25rrB");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction25rrB[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction25rrB[i] = log10(chimesRateTables.RatesTables->reaction25rrB[i]);

  dataset = H5Dopen(file_id, "/Rates/k25di");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction25di[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction25di[i] = log10(chimesRateTables.RatesTables->reaction25di[i]);

  dataset = H5Dopen(file_id, "/Rates/k26");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction26[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction26[i] = log10(chimesRateTables.RatesTables->reaction26[i]);

  dataset = H5Dopen(file_id, "/Rates/k27");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction27[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction27[i] = log10(chimesRateTables.RatesTables->reaction27[i]);

  dataset = H5Dopen(file_id, "/Rates/k83");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction83[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction83[i] = log10(chimesRateTables.RatesTables->reaction83[i]);

  dataset = H5Dopen(file_id, "/Rates/k84");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction84[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction84[i] = log10(max(chimesRateTables.RatesTables->reaction84[i], 1.0e-300));

  dataset = H5Dopen(file_id, "/Rates/k87");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction87[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction87[i] = log10(chimesRateTables.RatesTables->reaction87[i]);

  dataset = H5Dopen(file_id, "/Rates/k97");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction97[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction97[i] = log10(chimesRateTables.RatesTables->reaction97[i]);

  dataset = H5Dopen(file_id, "/Rates/k101");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction101[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction101[i] = log10(chimesRateTables.RatesTables->reaction101[i]);

  dataset = H5Dopen(file_id, "/Rates/k106");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction106[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction106[i] = log10(chimesRateTables.RatesTables->reaction106[i]);

  dataset = H5Dopen(file_id, "/Rates/k107");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction107[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction107[i] = log10(chimesRateTables.RatesTables->reaction107[i]);

  dataset = H5Dopen(file_id, "/Rates/k112");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction112[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction112[i] = log10(chimesRateTables.RatesTables->reaction112[i]);

  dataset = H5Dopen(file_id, "/Rates/k113");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction113[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction113[i] = log10(chimesRateTables.RatesTables->reaction113[i]);

  dataset = H5Dopen(file_id, "/Rates/k114");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction114[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction114[i] = log10(chimesRateTables.RatesTables->reaction114[i]);

  dataset = H5Dopen(file_id, "/Rates/k116");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction116[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction116[i] = log10(chimesRateTables.RatesTables->reaction116[i]);

  dataset = H5Dopen(file_id, "/Rates/k117");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction117[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction117[i] = log10(chimesRateTables.RatesTables->reaction117[i]);

  dataset = H5Dopen(file_id, "/Rates/k120");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction120[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction120[i] = log10(chimesRateTables.RatesTables->reaction120[i]);

  dataset = H5Dopen(file_id, "/Rates/k121");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction121[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction121[i] = log10(chimesRateTables.RatesTables->reaction121[i]);

  dataset = H5Dopen(file_id, "/Rates/k122");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction122[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction122[i] = log10(chimesRateTables.RatesTables->reaction122[i]);

  dataset = H5Dopen(file_id, "/Rates/k124");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction124[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction124[i] = log10(chimesRateTables.RatesTables->reaction124[i]);

  dataset = H5Dopen(file_id, "/Rates/k128");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction128[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction128[i] = log10(chimesRateTables.RatesTables->reaction128[i]);

  dataset = H5Dopen(file_id, "/Rates/k129");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction129[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction129[i] = log10(chimesRateTables.RatesTables->reaction129[i]);

  dataset = H5Dopen(file_id, "/Rates/k130");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction130[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction130[i] = log10(chimesRateTables.RatesTables->reaction130[i]);

  dataset = H5Dopen(file_id, "/Rates/k131");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction131[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction131[i] = log10(chimesRateTables.RatesTables->reaction131[i]);

  dataset = H5Dopen(file_id, "/Rates/k133");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction133[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction133[i] = log10(chimesRateTables.RatesTables->reaction133[i]);

  dataset = H5Dopen(file_id, "/Rates/k134");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction134[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction134[i] = log10(chimesRateTables.RatesTables->reaction134[i]);

  dataset = H5Dopen(file_id, "/Rates/k135");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction135[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction135[i] = log10(chimesRateTables.RatesTables->reaction135[i]);

  dataset = H5Dopen(file_id, "/Rates/k136");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction136[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction136[i] = log10(chimesRateTables.RatesTables->reaction136[i]);

  dataset = H5Dopen(file_id, "/Rates/k137");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction137[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction137[i] = log10(chimesRateTables.RatesTables->reaction137[i]);

  dataset = H5Dopen(file_id, "/Rates/k138");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction138[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction138[i] = log10(chimesRateTables.RatesTables->reaction138[i]);

  dataset = H5Dopen(file_id, "/Rates/k139");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction139[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction139[i] = log10(chimesRateTables.RatesTables->reaction139[i]);

  dataset = H5Dopen(file_id, "/Rates/k142");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction142[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction142[i] = log10(chimesRateTables.RatesTables->reaction142[i]);

  dataset = H5Dopen(file_id, "/Rates/k147");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction147[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction147[i] = log10(chimesRateTables.RatesTables->reaction147[i]);

  dataset = H5Dopen(file_id, "/Rates/k150");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction150[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction150[i] = log10(chimesRateTables.RatesTables->reaction150[i]);

  dataset = H5Dopen(file_id, "/Rates/k157");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction157[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction157[i] = log10(chimesRateTables.RatesTables->reaction157[i]);

  dataset = H5Dopen(file_id, "/Rates/k186");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction186[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction186[i] = log10(chimesRateTables.RatesTables->reaction186[i]);

  dataset = H5Dopen(file_id, "/Rates/k187");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction187[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction187[i] = log10(chimesRateTables.RatesTables->reaction187[i]);

  dataset = H5Dopen(file_id, "/Rates/k188");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction188[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction188[i] = log10(chimesRateTables.RatesTables->reaction188[i]);

  dataset = H5Dopen(file_id, "/Rates/k189");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction189[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction189[i] = log10(chimesRateTables.RatesTables->reaction189[i]);

  dataset = H5Dopen(file_id, "/Rates/k190");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction190[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction190[i] = log10(chimesRateTables.RatesTables->reaction190[i]);

  dataset = H5Dopen(file_id, "/Rates/k191");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction191[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction191[i] = log10(chimesRateTables.RatesTables->reaction191[i]);

  dataset = H5Dopen(file_id, "/Rates/k192");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction192[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction192[i] = log10(chimesRateTables.RatesTables->reaction192[i]);

  dataset = H5Dopen(file_id, "/Rates/k193");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction193[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction193[i] = log10(chimesRateTables.RatesTables->reaction193[i]);

  dataset = H5Dopen(file_id, "/Rates/k194");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction194[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction194[i] = log10(chimesRateTables.RatesTables->reaction194[i]);

  dataset = H5Dopen(file_id, "/Rates/k195");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction195[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction195[i] = log10(chimesRateTables.RatesTables->reaction195[i]);

  dataset = H5Dopen(file_id, "/Rates/k196");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction196[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction196[i] = log10(chimesRateTables.RatesTables->reaction196[i]);

  dataset = H5Dopen(file_id, "/Rates/k197");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction197[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction197[i] = log10(chimesRateTables.RatesTables->reaction197[i]);

  dataset = H5Dopen(file_id, "/Rates/k198");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction198[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction198[i] = log10(chimesRateTables.RatesTables->reaction198[i]);

  dataset = H5Dopen(file_id, "/Rates/k199");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction199[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction199[i] = log10(chimesRateTables.RatesTables->reaction199[i]);

  dataset = H5Dopen(file_id, "/Rates/k200");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction200[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction200[i] = log10(chimesRateTables.RatesTables->reaction200[i]);

  dataset = H5Dopen(file_id, "/Rates/k201");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction201[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction201[i] = log10(chimesRateTables.RatesTables->reaction201[i]);

  dataset = H5Dopen(file_id, "/Rates/k202");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction202[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction202[i] = log10(chimesRateTables.RatesTables->reaction202[i]);

  dataset = H5Dopen(file_id, "/Rates/k203");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction203[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction203[i] = log10(chimesRateTables.RatesTables->reaction203[i]);

  dataset = H5Dopen(file_id, "/Rates/k204");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction204[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction204[i] = log10(chimesRateTables.RatesTables->reaction204[i]);

  dataset = H5Dopen(file_id, "/Rates/k205");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction205[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction205[i] = log10(chimesRateTables.RatesTables->reaction205[i]);

  dataset = H5Dopen(file_id, "/Rates/k206");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction206[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction206[i] = log10(chimesRateTables.RatesTables->reaction206[i]);

  dataset = H5Dopen(file_id, "/Rates/k207");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction207[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction207[i] = log10(chimesRateTables.RatesTables->reaction207[i]);

  dataset = H5Dopen(file_id, "/Rates/k220");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction220[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction220[i] = log10(chimesRateTables.RatesTables->reaction220[i]);

  dataset = H5Dopen(file_id, "/Rates/k221");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction221[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction221[i] = log10(chimesRateTables.RatesTables->reaction221[i]);

  dataset = H5Dopen(file_id, "/Rates/k222");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction222[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction222[i] = log10(chimesRateTables.RatesTables->reaction222[i]);

  dataset = H5Dopen(file_id, "/Rates/k223");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction223[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction223[i] = log10(chimesRateTables.RatesTables->reaction223[i]);

  dataset = H5Dopen(file_id, "/Rates/k225");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction225[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction225[i] = log10(chimesRateTables.RatesTables->reaction225[i]);

  dataset = H5Dopen(file_id, "/Rates/k226");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction226[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction226[i] = log10(chimesRateTables.RatesTables->reaction226[i]);

  dataset = H5Dopen(file_id, "/Rates/k227");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction227[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction227[i] = log10(chimesRateTables.RatesTables->reaction227[i]);

  dataset = H5Dopen(file_id, "/Rates/k228");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction228[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction228[i] = log10(chimesRateTables.RatesTables->reaction228[i]);

  dataset = H5Dopen(file_id, "/Rates/k229");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction229[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction229[i] = log10(chimesRateTables.RatesTables->reaction229[i]);

  dataset = H5Dopen(file_id, "/Rates/k230");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction230[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction230[i] = log10(chimesRateTables.RatesTables->reaction230[i]);

  dataset = H5Dopen(file_id, "/Rates/k231");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction231[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction231[i] = log10(chimesRateTables.RatesTables->reaction231[i]);

  dataset = H5Dopen(file_id, "/Rates/k232");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction232[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction232[i] = log10(chimesRateTables.RatesTables->reaction232[i]);

  dataset = H5Dopen(file_id, "/Rates/k233");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction233[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction233[i] = log10(chimesRateTables.RatesTables->reaction233[i]);

  dataset = H5Dopen(file_id, "/Rates/k234");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction234[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction234[i] = log10(chimesRateTables.RatesTables->reaction234[i]);

  dataset = H5Dopen(file_id, "/Rates/k235");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction235[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction235[i] = log10(chimesRateTables.RatesTables->reaction235[i]);

  dataset = H5Dopen(file_id, "/Rates/k236");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction236[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction236[i] = log10(chimesRateTables.RatesTables->reaction236[i]);

  dataset = H5Dopen(file_id, "/Rates/k237");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction237[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction237[i] = log10(chimesRateTables.RatesTables->reaction237[i]);

  dataset = H5Dopen(file_id, "/Rates/k238");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction238[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction238[i] = log10(chimesRateTables.RatesTables->reaction238[i]);

  dataset = H5Dopen(file_id, "/Rates/k283");
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.RatesTables->reaction283[0]);
  H5Dclose(dataset);

  for(i = 0; i < chimesRateTables.RatesTables->N_Temperatures; i++)
    chimesRateTables.RatesTables->reaction283[i] = log10(chimesRateTables.RatesTables->reaction283[i]);

  H5Fclose(file_id);
}

/*
 * ----------------------------------------------------------------------
 * This routine allocates header arrays for the  molecular cooling tables
 * ----------------------------------------------------------------------
 */

void allocate_molecular_header_arrays(struct globalVariables *myGlobalVars)
{
  chimesRateTables.mol_cooling_table_CO_rot_T =
      (double *)mymalloc("Chimes_CO_T1", chimesRateTables.mol_cooling_table_dimensions[0] * sizeof(double));
  chimesRateTables.mol_cooling_table_CO_rot_N =
      (double *)mymalloc("Chimes_CO_N1", chimesRateTables.mol_cooling_table_dimensions[1] * sizeof(double));
  chimesRateTables.mol_cooling_table_CO_vib_T =
      (double *)mymalloc("Chimes_CO_T2", chimesRateTables.mol_cooling_table_dimensions[2] * sizeof(double));
  chimesRateTables.mol_cooling_table_CO_vib_N =
      (double *)mymalloc("Chimes_CO_N2", chimesRateTables.mol_cooling_table_dimensions[3] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_rot_hiT_T =
      (double *)mymalloc("Chimes_H2O_T1", chimesRateTables.mol_cooling_table_dimensions[4] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_rot_hiT_N =
      (double *)mymalloc("Chimes_H2O_N1", chimesRateTables.mol_cooling_table_dimensions[5] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_rot_lowT_T =
      (double *)mymalloc("Chimes_H2O_T2", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_rot_lowT_N =
      (double *)mymalloc("Chimes_H2O_N2", chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_vib_T =
      (double *)mymalloc("Chimes_H2O_T3", chimesRateTables.mol_cooling_table_dimensions[8] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_vib_N =
      (double *)mymalloc("Chimes_H2O_N3", chimesRateTables.mol_cooling_table_dimensions[9] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2_lte_temperatures =
      (double *)mymalloc("Chimes_H2_LTE", chimesRateTables.mol_cooling_table_dimensions[10] * sizeof(double));
}

/*
 * ----------------------------------------------------------------------
 * This routine reads in the header of the ion cooling table files
 * ----------------------------------------------------------------------
 */

void ReadMolecularCoolingHeader(struct globalVariables *myGlobalVars)
{
  hid_t tempfile_id, dataset_id;
  herr_t status;

  tempfile_id = H5Fopen(myGlobalVars->MolecularTablePath, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* Determine the length of each header array */
  dataset_id = H5Dopen(tempfile_id, "/Bin_Sizes");
  status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.mol_cooling_table_dimensions);
  status     = H5Dclose(dataset_id);

  /* allocate arrays for cooling table header */
  allocate_molecular_header_arrays(myGlobalVars);

  /* Fill the header arrays */
  dataset_id = H5Dopen(tempfile_id, "/CO_rot_T");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_CO_rot_T);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/CO_rot_N");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_CO_rot_N);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/CO_vib_T");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_CO_vib_T);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/CO_vib_N");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_CO_vib_N);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2O_rot_hiT_T");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2O_rot_hiT_T);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2O_rot_hiT_N");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2O_rot_hiT_N);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2O_rot_lowT_T");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2O_rot_lowT_T);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2O_rot_lowT_N");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2O_rot_lowT_N);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2O_vib_T");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2O_vib_T);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2O_vib_N");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2O_vib_N);
  status     = H5Dclose(dataset_id);

  dataset_id = H5Dopen(tempfile_id, "/H2_lte_temperatures");
  status =
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.mol_cooling_table_H2_lte_temperatures);
  status = H5Dclose(dataset_id);

  status = H5Fclose(tempfile_id);
}

/*
 * -------------------------------------------------
 * This routine reads in the molecular cooling table
 * -------------------------------------------------
 */

void GetMolecularCoolingTable(struct globalVariables *myGlobalVars)
{
  hid_t file_id, dataset_id, dataspace_id, memspace_id;
  herr_t status;
  int rank, i, j;
  hsize_t dims[2];
  hsize_t dims2D[2], offset2D[2], count2D[2];
  dims[1] = 1;

  double array0[chimesRateTables.mol_cooling_table_dimensions[0]];
  double array2[chimesRateTables.mol_cooling_table_dimensions[2]];
  double array4[chimesRateTables.mol_cooling_table_dimensions[4]];
  double array6[chimesRateTables.mol_cooling_table_dimensions[6]];
  double array8[chimesRateTables.mol_cooling_table_dimensions[8]];
  double array10[chimesRateTables.mol_cooling_table_dimensions[10]];

  file_id = H5Fopen(myGlobalVars->MolecularTablePath, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "/CO_rot_L0");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array0);
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[0]; i++)
    chimesRateTables.mol_cooling_table_CO_rot_L0[i] = array0[i];
  H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/H2O_rot_hiT_L0");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array4);
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[4]; i++)
    chimesRateTables.mol_cooling_table_H2O_rot_hiT_L0[i] = array4[i];
  H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/H2Oortho_rot_lowT_L0");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array6);
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_L0[i] = array6[i];
  H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/H2Opara_rot_lowT_L0");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array6);
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_L0[i] = array6[i];
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/CO_rot_Llte");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  dims[0]     = chimesRateTables.mol_cooling_table_dimensions[0];
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[1]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[0];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array0);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[0]; j++)
        chimesRateTables.mol_cooling_table_CO_rot_Llte[j][i] = array0[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/CO_rot_nhalf");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[1]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[0];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array0);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[0]; j++)
        chimesRateTables.mol_cooling_table_CO_rot_nhalf[j][i] = array0[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/CO_rot_a");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[1]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[0];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array0);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[0]; j++)
        chimesRateTables.mol_cooling_table_CO_rot_a[j][i] = array0[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/CO_vib_Llte");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  dims[0]     = chimesRateTables.mol_cooling_table_dimensions[2];
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[3]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[2];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array2);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[2]; j++)
        chimesRateTables.mol_cooling_table_CO_vib_Llte[j][i] = array2[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2O_rot_hiT_Llte");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  dims[0]     = chimesRateTables.mol_cooling_table_dimensions[4];
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[5]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[4];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array4);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[4]; j++)
        chimesRateTables.mol_cooling_table_H2O_rot_hiT_Llte[j][i] = array4[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2O_rot_hiT_nhalf");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[5]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[4];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array4);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[4]; j++)
        chimesRateTables.mol_cooling_table_H2O_rot_hiT_nhalf[j][i] = array4[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2O_rot_hiT_a");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[5]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[4];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array4);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[4]; j++)
        chimesRateTables.mol_cooling_table_H2O_rot_hiT_a[j][i] = array4[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2Oortho_rot_lowT_Llte");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  dims[0]     = chimesRateTables.mol_cooling_table_dimensions[6];
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[7]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[6];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array6);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[6]; j++)
        chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_Llte[j][i] = array6[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2Oortho_rot_lowT_nhalf");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[7]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[6];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array6);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[6]; j++)
        chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_nhalf[j][i] = array6[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2Oortho_rot_lowT_a");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[7]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[6];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array6);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[6]; j++)
        chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_a[j][i] = array6[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2Opara_rot_lowT_Llte");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  dims[0]     = chimesRateTables.mol_cooling_table_dimensions[6];
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[7]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[6];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array6);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[6]; j++)
        chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_Llte[j][i] = array6[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2Opara_rot_lowT_nhalf");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[7]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[6];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array6);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[6]; j++)
        chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_nhalf[j][i] = array6[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2Opara_rot_lowT_a");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[7]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[6];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array6);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[6]; j++)
        chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_a[j][i] = array6[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id   = H5Dopen(file_id, "/H2O_vib_Llte");
  dataspace_id = H5Dget_space(dataset_id);
  status       = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);

  dims[0]     = chimesRateTables.mol_cooling_table_dimensions[8];
  rank        = 1;
  memspace_id = H5Screate_simple(rank, dims, NULL);

  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[9]; i++)
    {
      offset2D[0] = 0;
      offset2D[1] = i;
      count2D[0]  = chimesRateTables.mol_cooling_table_dimensions[8];
      count2D[1]  = 1;
      status      = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);

      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, array8);

      for(j = 0; j < chimesRateTables.mol_cooling_table_dimensions[8]; j++)
        chimesRateTables.mol_cooling_table_H2O_vib_Llte[j][i] = array8[j];
    }
  H5Sclose(memspace_id);
  H5Sclose(dataspace_id);
  H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "/H2_cooling_lte");
  status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array10);
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[10]; i++)
    chimesRateTables.mol_cooling_table_H2_lte[i] = array10[i];
  H5Dclose(dataset_id);

  H5Fclose(file_id);
}

/*
 * -----------------------------------------------------------------------
 * This routine initialises the cooling at the beginning of the simulation
 * -----------------------------------------------------------------------
 */
void initialise_cooling(struct globalVariables *myGlobalVars)
{
  int i;

  /* Read in the molecular cooling table header */
  ReadMolecularCoolingHeader(myGlobalVars);

  /* Allocate memory for the molecular cooling tables */
  chimesRateTables.mol_cooling_table_CO_rot_L0 =
      (double *)mymalloc("Chimes_CO_L0", chimesRateTables.mol_cooling_table_dimensions[0] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2O_rot_hiT_L0 =
      (double *)mymalloc("Chimes_H2O_L0", chimesRateTables.mol_cooling_table_dimensions[4] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_L0 =
      (double *)mymalloc("Chimes_H2O_oL0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_L0 =
      (double *)mymalloc("Chimes_H20_pL0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double));
  chimesRateTables.mol_cooling_table_H2_lte =
      (double *)mymalloc("Chimes_H2_LTE", chimesRateTables.mol_cooling_table_dimensions[10] * sizeof(double));

  chimesRateTables.mol_cooling_table_CO_rot_Llte =
      (double **)mymalloc("Chimes_CO_rot_Llte_dim0", chimesRateTables.mol_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.mol_cooling_table_CO_rot_Llte[0] =
      (double *)mymalloc("Chimes_CO_rot_Llte_dim1", chimesRateTables.mol_cooling_table_dimensions[0] *
                                                        chimesRateTables.mol_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[0]; i++)
    chimesRateTables.mol_cooling_table_CO_rot_Llte[i] =
        &(chimesRateTables.mol_cooling_table_CO_rot_Llte[0][i * chimesRateTables.mol_cooling_table_dimensions[1]]);

  chimesRateTables.mol_cooling_table_CO_rot_nhalf =
      (double **)mymalloc("Chimes_CO_rot_nhalf_dim0", chimesRateTables.mol_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.mol_cooling_table_CO_rot_nhalf[0] =
      (double *)mymalloc("Chimes_CO_rot_nhalf_dim1", chimesRateTables.mol_cooling_table_dimensions[0] *
                                                         chimesRateTables.mol_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[0]; i++)
    chimesRateTables.mol_cooling_table_CO_rot_nhalf[i] =
        &(chimesRateTables.mol_cooling_table_CO_rot_nhalf[0][i * chimesRateTables.mol_cooling_table_dimensions[1]]);

  chimesRateTables.mol_cooling_table_CO_rot_a =
      (double **)mymalloc("Chimes_CO_rot_a_dim0", chimesRateTables.mol_cooling_table_dimensions[0] * sizeof(double *));
  chimesRateTables.mol_cooling_table_CO_rot_a[0] =
      (double *)mymalloc("Chimes_CO_rot_a_dim1", chimesRateTables.mol_cooling_table_dimensions[0] *
                                                     chimesRateTables.mol_cooling_table_dimensions[1] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[0]; i++)
    chimesRateTables.mol_cooling_table_CO_rot_a[i] =
        &(chimesRateTables.mol_cooling_table_CO_rot_a[0][i * chimesRateTables.mol_cooling_table_dimensions[1]]);

  chimesRateTables.mol_cooling_table_CO_vib_Llte =
      (double **)mymalloc("Chimes_CO_vib_Llte_dim0", chimesRateTables.mol_cooling_table_dimensions[2] * sizeof(double *));
  chimesRateTables.mol_cooling_table_CO_vib_Llte[0] =
      (double *)mymalloc("Chimes_CO_vib_Llte_dim1", chimesRateTables.mol_cooling_table_dimensions[2] *
                                                        chimesRateTables.mol_cooling_table_dimensions[3] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[2]; i++)
    chimesRateTables.mol_cooling_table_CO_vib_Llte[i] =
        &(chimesRateTables.mol_cooling_table_CO_vib_Llte[0][i * chimesRateTables.mol_cooling_table_dimensions[3]]);

  chimesRateTables.mol_cooling_table_H2O_rot_hiT_Llte =
      (double **)mymalloc("Chimes_H2O_Llte_dim0", chimesRateTables.mol_cooling_table_dimensions[4] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2O_rot_hiT_Llte[0] =
      (double *)mymalloc("Chimes_H2O_Llte_dim1", chimesRateTables.mol_cooling_table_dimensions[4] *
                                                     chimesRateTables.mol_cooling_table_dimensions[5] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[4]; i++)
    chimesRateTables.mol_cooling_table_H2O_rot_hiT_Llte[i] =
        &(chimesRateTables.mol_cooling_table_H2O_rot_hiT_Llte[0][i * chimesRateTables.mol_cooling_table_dimensions[5]]);

  chimesRateTables.mol_cooling_table_H2O_rot_hiT_nhalf =
      (double **)mymalloc("Chimes_H2O_nhalf_dim0", chimesRateTables.mol_cooling_table_dimensions[4] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2O_rot_hiT_nhalf[0] =
      (double *)mymalloc("Chimes_H2O_nhalf_dim1", chimesRateTables.mol_cooling_table_dimensions[4] *
                                                      chimesRateTables.mol_cooling_table_dimensions[5] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[4]; i++)
    chimesRateTables.mol_cooling_table_H2O_rot_hiT_nhalf[i] =
        &(chimesRateTables.mol_cooling_table_H2O_rot_hiT_nhalf[0][i * chimesRateTables.mol_cooling_table_dimensions[5]]);

  chimesRateTables.mol_cooling_table_H2O_rot_hiT_a =
      (double **)mymalloc("Chimes_H2O_a_dim0", chimesRateTables.mol_cooling_table_dimensions[4] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2O_rot_hiT_a[0] =
      (double *)mymalloc("Chimes_H2O_a_dim1", chimesRateTables.mol_cooling_table_dimensions[4] *
                                                  chimesRateTables.mol_cooling_table_dimensions[5] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[4]; i++)
    chimesRateTables.mol_cooling_table_H2O_rot_hiT_a[i] =
        &(chimesRateTables.mol_cooling_table_H2O_rot_hiT_a[0][i * chimesRateTables.mol_cooling_table_dimensions[5]]);

  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_Llte =
      (double **)mymalloc("Chimes_oH2O_Llte_dim0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_Llte[0] =
      (double *)mymalloc("Chimes_oH2O_Llte_dim1", chimesRateTables.mol_cooling_table_dimensions[6] *
                                                      chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_Llte[i] =
        &(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_Llte[0][i * chimesRateTables.mol_cooling_table_dimensions[7]]);

  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_nhalf =
      (double **)mymalloc("Chimes_oH2O_nhalf_dim0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_nhalf[0] =
      (double *)mymalloc("Chimes_oH2O_nhalf_dim1", chimesRateTables.mol_cooling_table_dimensions[6] *
                                                       chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_nhalf[i] =
        &(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_nhalf[0][i * chimesRateTables.mol_cooling_table_dimensions[7]]);

  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_a =
      (double **)mymalloc("Chimes_oH2O_a_dim0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_a[0] =
      (double *)mymalloc("Chimes_oH2O_a_dim1", chimesRateTables.mol_cooling_table_dimensions[6] *
                                                   chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_a[i] =
        &(chimesRateTables.mol_cooling_table_H2Oortho_rot_lowT_a[0][i * chimesRateTables.mol_cooling_table_dimensions[7]]);

  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_Llte =
      (double **)mymalloc("Chimes_pH2O_Llte_dim0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_Llte[0] =
      (double *)mymalloc("Chimes_pH2O_Llte_dim1", chimesRateTables.mol_cooling_table_dimensions[6] *
                                                      chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_Llte[i] =
        &(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_Llte[0][i * chimesRateTables.mol_cooling_table_dimensions[7]]);

  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_nhalf =
      (double **)mymalloc("Chimes_pH2O_nhalf_dim0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_nhalf[0] =
      (double *)mymalloc("Chimes_pH2O_nhalf_dim1", chimesRateTables.mol_cooling_table_dimensions[6] *
                                                       chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_nhalf[i] =
        &(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_nhalf[0][i * chimesRateTables.mol_cooling_table_dimensions[7]]);

  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_a =
      (double **)mymalloc("Chimes_pH2O_a_dim0", chimesRateTables.mol_cooling_table_dimensions[6] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_a[0] =
      (double *)mymalloc("Chimes_pH2O_a_dim1", chimesRateTables.mol_cooling_table_dimensions[6] *
                                                   chimesRateTables.mol_cooling_table_dimensions[7] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[6]; i++)
    chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_a[i] =
        &(chimesRateTables.mol_cooling_table_H2Opara_rot_lowT_a[0][i * chimesRateTables.mol_cooling_table_dimensions[7]]);

  chimesRateTables.mol_cooling_table_H2O_vib_Llte =
      (double **)mymalloc("Chimes_H2O_vib_dim0", chimesRateTables.mol_cooling_table_dimensions[8] * sizeof(double *));
  chimesRateTables.mol_cooling_table_H2O_vib_Llte[0] =
      (double *)mymalloc("Chimes_H2Ovib_dim1", chimesRateTables.mol_cooling_table_dimensions[8] *
                                                   chimesRateTables.mol_cooling_table_dimensions[9] * sizeof(double));
  for(i = 0; i < chimesRateTables.mol_cooling_table_dimensions[8]; i++)
    chimesRateTables.mol_cooling_table_H2O_vib_Llte[i] =
        &(chimesRateTables.mol_cooling_table_H2O_vib_Llte[0][i * chimesRateTables.mol_cooling_table_dimensions[9]]);

  /* Read in the molecular cooling table */
  GetMolecularCoolingTable(myGlobalVars);
}
struct Reactions_Structure *add_new_reaction(struct Reactions_Structure *prev_reaction, int reactant1, int reactant2, int reactant3,
                                             int N_reactants, int product1, int product2, int product3, int N_products, double *rate,
                                             int incl_mol, int *speciesIndices)
{
  int i;
  int incl_reaction = 1;
  int reactants[3]  = {reactant1, reactant2, reactant3};
  int products[3]   = {product1, product2, product3};
  struct Reactions_Structure *conductor;
  /* Note: if there are less than 3 reactants/species,
   * you can put anything for the later ones */

  /* The argument 'incl_mol' tells us whether we are
   * including molecules (species with enum >= 137).
   * Also, the argument 'speciesIndices' is used to
   * determine whether all species are included in
   * the complete network, as their index will be
   * -1 if they are not included. */
  for(i = 0; i < N_reactants; i++)
    {
      if((reactants[i] >= H2 && incl_mol == 0) || speciesIndices[reactants[i]] == -1)
        {
          incl_reaction = 0;
          break;
        }
    }
  for(i = 0; i < N_products; i++)
    {
      if((products[i] >= H2 && incl_mol == 0) || speciesIndices[products[i]] == -1)
        {
          incl_reaction = 0;
          break;
        }
    }

  if(incl_reaction == 1)
    {
      conductor                    = mymalloc("Chimes_reaction", sizeof(struct Reactions_Structure));
      prev_reaction->next_reaction = conductor;
      conductor->no_of_reactants   = N_reactants;
      conductor->no_of_products    = N_products;
      if(product3 < 0) /* Auger ionisation */
        conductor->extra_Auger_electrons = -product3;
      else
        conductor->extra_Auger_electrons = 0;
      for(i = 0; i < N_reactants; i++)
        conductor->reactants[i] = speciesIndices[reactants[i]];
      for(i = 0; i < N_products; i++)
        conductor->products[i] = speciesIndices[products[i]];
      conductor->rate          = rate;
      conductor->flag_included = 1;
      conductor->next_reaction = NULL;
      return conductor;
    }
  else
    return prev_reaction;
}

void initialise_reactions(struct Reactions_Structure *root_node, int incl_mol, struct globalVariables *myGlobalVars,
                          struct All_rate_variables_structure *this_all_rates)
{
  int i, j, k;
  /* Starting at the root node, go through every reaction, adding them to
   * the linked list in turn. Note that the root node is blank - it does
   * not contain a reaction. In the RhsFn, this node gets skipped anyway. */
  struct Reactions_Structure *last_reaction;
  last_reaction = root_node;
  last_reaction =
      add_new_reaction(last_reaction, HI, elec, 0, 2, Hm, 0, 0, 1, &this_all_rates->rate_1, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Hm, HII, 0, 2, HI, HI, 0, 2, &this_all_rates->rate_5, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HI, elec, 0, 2, HII, elec, elec, 3, &this_all_rates->rate_11, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HII, elec, 0, 2, HI, 0, 0, 1, &this_all_rates->rate_13, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, Hm, elec, 0, 2, HI, elec, elec, 3, &this_all_rates->rate_15, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Hm, HI, 0, 2, HI, HI, elec, 3, &this_all_rates->rate_16, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeI, elec, 0, 2, HeII, elec, elec, 3, &this_all_rates->rate_24, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeII, elec, 0, 2, HeI, 0, 0, 1, &this_all_rates->rate_25, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeII, HI, 0, 2, HeI, HII, 0, 2, &this_all_rates->rate_26, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeI, HII, 0, 2, HeII, HI, 0, 2, &this_all_rates->rate_27, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CII, SiI, 0, 2, CI, SiII, 0, 2, &this_all_rates->rate_44, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HII, elec, 0, 2, HI, 0, 0, 1, &this_all_rates->rate_61, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeII, elec, 0, 2, HeI, 0, 0, 1, &this_all_rates->rate_63, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CII, elec, 0, 2, CI, 0, 0, 1, &this_all_rates->rate_64, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OII, elec, 0, 2, OI, 0, 0, 1, &this_all_rates->rate_65, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, SiII, elec, 0, 2, SiI, 0, 0, 1, &this_all_rates->rate_66, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, FeII, elec, 0, 2, FeI, 0, 0, 1, &this_all_rates->rate_77, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, FeI, CII, 0, 2, FeII, CI, 0, 2, &this_all_rates->rate_79, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, FeI, SiII, 0, 2, FeII, SiI, 0, 2, &this_all_rates->rate_80, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Hm, HI, 0, 2, H2, elec, 0, 2, &this_all_rates->rate_2, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HI, HII, 0, 2, H2p, 0, 0, 1, &this_all_rates->rate_3, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HI, H2p, 0, 2, H2, HII, 0, 2, &this_all_rates->rate_4, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2p, elec, 0, 2, HI, HI, 0, 2, &this_all_rates->rate_6, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, HII, 0, 2, H2p, HI, 0, 2, &this_all_rates->rate_7, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2, elec, 0, 2, HI, HI, elec, 3, &this_all_rates->rate_8, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, HI, 0, 2, HI, HI, HI, 3, &this_all_rates->rate_9, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, H2, 0, 2, H2, HI, HI, 3, &this_all_rates->rate_10, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, Hm, HII, 0, 2, H2p, elec, 0, 2, &this_all_rates->rate_17, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2p, Hm, 0, 2, H2, HI, 0, 2, &this_all_rates->rate_87, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, HeI, 0, 2, HI, HI, HeI, 3, &this_all_rates->rate_88, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2p, Hm, 0, 2, HI, HI, HI, 3, &this_all_rates->rate_83, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, elec, 0, 2, Hm, HI, 0, 2, &this_all_rates->rate_84, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2, HeII, 0, 2, HeI, HI, HII, 3, &this_all_rates->rate_85, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2, HeII, 0, 2, H2p, HeI, 0, 2, &this_all_rates->rate_86, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2p, 0, 0, 1, HI, HII, 0, 2, &this_all_rates->rate_52, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, 0, 0, 1, HI, HI, 0, 2, &this_all_rates->rate_53, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, 0, 0, 1, H2p, elec, 0, 2, &this_all_rates->rate_70, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HI, HI, 0, 2, H2, 0, 0, 1, &this_all_rates->rate_60, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, 0, 0, 1, HI, HI, 0, 2, &this_all_rates->rate_89, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2p, H2, 0, 2, H3p, HI, 0, 2, &this_all_rates->rate_97, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CO, H3p, 0, 2, HCOp, H2, 0, 2, &this_all_rates->rate_100, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeII, CO, 0, 2, CII, OI, HeI, 3, &this_all_rates->rate_101, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3p, elec, 0, 2, H2, HI, 0, 2, &this_all_rates->rate_106, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HCOp, elec, 0, 2, CO, HI, 0, 2, &this_all_rates->rate_107, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3p, FeI, 0, 2, FeII, H2, HI, 3, &this_all_rates->rate_108, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HCOp, 0, 0, 1, CO, HII, 0, 2, &this_all_rates->rate_111, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H3p, HI, 0, 2, H2p, H2, 0, 2, &this_all_rates->rate_112, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CO, HeII, 0, 2, CI, OII, HeI, 3, &this_all_rates->rate_113, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3p, elec, 0, 2, HI, HI, HI, 3, &this_all_rates->rate_114, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, HII, 0, 2, H3p, 0, 0, 1, &this_all_rates->rate_115, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, OI, 0, 2, CO, 0, 0, 1, &this_all_rates->rate_116, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, HI, 0, 2, OI, HI, HI, 3, &this_all_rates->rate_117, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HOCp, H2, 0, 2, HCOp, H2, 0, 2, &this_all_rates->rate_118, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HOCp, CO, 0, 2, HCOp, CO, 0, 2, &this_all_rates->rate_119, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, H2, 0, 2, CH, HI, 0, 2, &this_all_rates->rate_120, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, HI, 0, 2, CI, H2, 0, 2, &this_all_rates->rate_121, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, H2, 0, 2, CH2, HI, 0, 2, &this_all_rates->rate_122, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, CI, 0, 2, C2, HI, 0, 2, &this_all_rates->rate_123, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, OI, 0, 2, CO, HI, 0, 2, &this_all_rates->rate_124, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH2, HI, 0, 2, CH, H2, 0, 2, &this_all_rates->rate_125, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH2, OI, 0, 2, CO, HI, HI, 3, &this_all_rates->rate_126, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH2, OI, 0, 2, CO, H2, 0, 2, &this_all_rates->rate_127, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, C2, OI, 0, 2, CO, CI, 0, 2, &this_all_rates->rate_128, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OI, H2, 0, 2, OH, HI, 0, 2, &this_all_rates->rate_129, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, HI, 0, 2, OI, H2, 0, 2, &this_all_rates->rate_130, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, H2, 0, 2, H2O, HI, 0, 2, &this_all_rates->rate_131, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, CI, 0, 2, CO, HI, 0, 2, &this_all_rates->rate_132, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, OI, 0, 2, O2, HI, 0, 2, &this_all_rates->rate_133, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, OH, 0, 2, H2O, OI, 0, 2, &this_all_rates->rate_134, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2O, HI, 0, 2, H2, OH, 0, 2, &this_all_rates->rate_135, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, HI, 0, 2, OH, OI, 0, 2, &this_all_rates->rate_136, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, H2, 0, 2, OH, OH, 0, 2, &this_all_rates->rate_137, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, CI, 0, 2, CO, OI, 0, 2, &this_all_rates->rate_138, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CO, HI, 0, 2, CI, OH, 0, 2, &this_all_rates->rate_139, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, H2p, 0, 2, CHp, HI, 0, 2, &this_all_rates->rate_140, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, H3p, 0, 2, CHp, H2, 0, 2, &this_all_rates->rate_141, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CII, H2, 0, 2, CHp, HI, 0, 2, &this_all_rates->rate_142, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CHp, HI, 0, 2, CII, H2, 0, 2, &this_all_rates->rate_143, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CHp, H2, 0, 2, CH2p, HI, 0, 2, &this_all_rates->rate_144, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CHp, OI, 0, 2, COp, HI, 0, 2, &this_all_rates->rate_145, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2, HII, 0, 2, CHp, H2, 0, 2, &this_all_rates->rate_146, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2p, HI, 0, 2, CHp, H2, 0, 2, &this_all_rates->rate_147, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2p, H2, 0, 2, CH3p, HI, 0, 2, &this_all_rates->rate_148, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2p, OI, 0, 2, HCOp, HI, 0, 2, &this_all_rates->rate_149, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH3p, HI, 0, 2, CH2p, H2, 0, 2, &this_all_rates->rate_150, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH3p, OI, 0, 2, HCOp, H2, 0, 2, &this_all_rates->rate_151, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, C2, OII, 0, 2, COp, CI, 0, 2, &this_all_rates->rate_152, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OII, H2, 0, 2, OHp, HI, 0, 2, &this_all_rates->rate_153, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OI, H2p, 0, 2, OHp, HI, 0, 2, &this_all_rates->rate_154, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OI, H3p, 0, 2, OHp, H2, 0, 2, &this_all_rates->rate_155, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OH, H3p, 0, 2, H2Op, H2, 0, 2, &this_all_rates->rate_156, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, CII, 0, 2, COp, HI, 0, 2, &this_all_rates->rate_157, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OHp, H2, 0, 2, H2Op, HI, 0, 2, &this_all_rates->rate_158, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2Op, H2, 0, 2, H3Op, HI, 0, 2, &this_all_rates->rate_159, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, H3p, 0, 2, H3Op, H2, 0, 2, &this_all_rates->rate_160, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, CII, 0, 2, HCOp, HI, 0, 2, &this_all_rates->rate_161, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, CII, 0, 2, HOCp, HI, 0, 2, &this_all_rates->rate_162, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, CI, 0, 2, HCOp, H2, 0, 2, &this_all_rates->rate_163, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, CII, 0, 2, COp, OI, 0, 2, &this_all_rates->rate_164, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, CII, 0, 2, CO, OII, 0, 2, &this_all_rates->rate_165, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, O2, CH2p, 0, 2, HCOp, OH, 0, 2, &this_all_rates->rate_166, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2p, CI, 0, 2, COp, OI, 0, 2, &this_all_rates->rate_167, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CO, H3p, 0, 2, HOCp, H2, 0, 2, &this_all_rates->rate_168, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HCOp, CI, 0, 2, CO, CHp, 0, 2, &this_all_rates->rate_169, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HCOp, H2O, 0, 2, CO, H3Op, 0, 2, &this_all_rates->rate_170, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, HII, 0, 2, CHp, HI, 0, 2, &this_all_rates->rate_171, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2, HII, 0, 2, CH2p, HI, 0, 2, &this_all_rates->rate_172, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2, HeII, 0, 2, CII, HeI, H2, 3, &this_all_rates->rate_173, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, C2, HeII, 0, 2, CII, CI, HeI, 3, &this_all_rates->rate_174, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, HII, 0, 2, OHp, HI, 0, 2, &this_all_rates->rate_175, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OH, HeII, 0, 2, OII, HeI, HI, 3, &this_all_rates->rate_176, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, HII, 0, 2, H2Op, HI, 0, 2, &this_all_rates->rate_177, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, HeII, 0, 2, OH, HeI, HII, 3, &this_all_rates->rate_178, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, HeII, 0, 2, OHp, HeI, HI, 3, &this_all_rates->rate_179, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, HeII, 0, 2, H2Op, HeI, 0, 2, &this_all_rates->rate_180, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, HII, 0, 2, O2p, HI, 0, 2, &this_all_rates->rate_181, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, O2, HeII, 0, 2, O2p, HeI, 0, 2, &this_all_rates->rate_182, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, O2, HeII, 0, 2, OII, OI, HeI, 3, &this_all_rates->rate_183, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2p, CI, 0, 2, O2, CII, 0, 2, &this_all_rates->rate_184, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, COp, HI, 0, 2, CO, HII, 0, 2, &this_all_rates->rate_185, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Cm, HII, 0, 2, CI, HI, 0, 2, &this_all_rates->rate_186, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Om, HII, 0, 2, OI, HI, 0, 2, &this_all_rates->rate_187, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HeII, Hm, 0, 2, HeI, HI, 0, 2, &this_all_rates->rate_188, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CHp, elec, 0, 2, CI, HI, 0, 2, &this_all_rates->rate_189, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2p, elec, 0, 2, CH, HI, 0, 2, &this_all_rates->rate_190, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2p, elec, 0, 2, CI, HI, HI, 3, &this_all_rates->rate_191, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2p, elec, 0, 2, CI, H2, 0, 2, &this_all_rates->rate_192, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH3p, elec, 0, 2, CH2, HI, 0, 2, &this_all_rates->rate_193, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH3p, elec, 0, 2, CH, H2, 0, 2, &this_all_rates->rate_194, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH3p, elec, 0, 2, CH, HI, HI, 3, &this_all_rates->rate_195, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OHp, elec, 0, 2, OI, HI, 0, 2, &this_all_rates->rate_196, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2Op, elec, 0, 2, OI, HI, HI, 3, &this_all_rates->rate_197, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2Op, elec, 0, 2, OI, H2, 0, 2, &this_all_rates->rate_198, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2Op, elec, 0, 2, OH, HI, 0, 2, &this_all_rates->rate_199, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, elec, 0, 2, HI, H2O, 0, 2, &this_all_rates->rate_200, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, elec, 0, 2, OH, H2, 0, 2, &this_all_rates->rate_201, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, elec, 0, 2, OH, HI, HI, 3, &this_all_rates->rate_202, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, elec, 0, 2, OI, HI, H2, 3, &this_all_rates->rate_203, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, O2p, elec, 0, 2, OI, OI, 0, 2, &this_all_rates->rate_204, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, COp, elec, 0, 2, CI, OI, 0, 2, &this_all_rates->rate_205, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HCOp, elec, 0, 2, OH, CI, 0, 2, &this_all_rates->rate_206, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HOCp, elec, 0, 2, CO, HI, 0, 2, &this_all_rates->rate_207, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Hm, CI, 0, 2, CH, elec, 0, 2, &this_all_rates->rate_208, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Hm, OI, 0, 2, OH, elec, 0, 2, &this_all_rates->rate_209, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, Hm, OH, 0, 2, H2O, elec, 0, 2, &this_all_rates->rate_210, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Cm, HI, 0, 2, CH, elec, 0, 2, &this_all_rates->rate_211, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, Cm, H2, 0, 2, CH2, elec, 0, 2, &this_all_rates->rate_212, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Cm, OI, 0, 2, CO, elec, 0, 2, &this_all_rates->rate_213, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Om, HI, 0, 2, OH, elec, 0, 2, &this_all_rates->rate_214, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, Om, H2, 0, 2, H2O, elec, 0, 2, &this_all_rates->rate_215, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Om, CI, 0, 2, CO, elec, 0, 2, &this_all_rates->rate_216, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, elec, 0, 2, Cm, 0, 0, 1, &this_all_rates->rate_217, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, HI, 0, 2, CH, 0, 0, 1, &this_all_rates->rate_218, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, H2, 0, 2, CH2, 0, 0, 1, &this_all_rates->rate_219, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CI, CI, 0, 2, C2, 0, 0, 1, &this_all_rates->rate_220, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CII, HI, 0, 2, CHp, 0, 0, 1, &this_all_rates->rate_221, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CII, H2, 0, 2, CH2p, 0, 0, 1, &this_all_rates->rate_222, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CII, OI, 0, 2, COp, 0, 0, 1, &this_all_rates->rate_223, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OI, elec, 0, 2, Om, 0, 0, 1, &this_all_rates->rate_224, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OI, HI, 0, 2, OH, 0, 0, 1, &this_all_rates->rate_225, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OI, OI, 0, 2, O2, 0, 0, 1, &this_all_rates->rate_226, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, HI, 0, 2, H2O, 0, 0, 1, &this_all_rates->rate_227, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HI, HI, HI, 3, H2, HI, 0, 2, &this_all_rates->rate_228, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, HI, HI, H2, 3, H2, H2, 0, 2, &this_all_rates->rate_229, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, HI, HI, HeI, 3, H2, HeI, 0, 2, &this_all_rates->rate_230, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CI, CI, HeI, 3, C2, HeI, 0, 2, &this_all_rates->rate_231, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CI, OI, HeI, 3, CO, HeI, 0, 2, &this_all_rates->rate_232, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CII, OI, HeI, 3, COp, HeI, 0, 2, &this_all_rates->rate_233, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CI, OII, HeI, 3, COp, HeI, 0, 2, &this_all_rates->rate_234, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OI, HI, HeI, 3, OH, HeI, 0, 2, &this_all_rates->rate_235, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OH, HI, HeI, 3, H2O, HeI, 0, 2, &this_all_rates->rate_236, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OI, OI, HeI, 3, O2, HeI, 0, 2, &this_all_rates->rate_237, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, OI, CH, 0, 2, HCOp, elec, 0, 2, &this_all_rates->rate_238, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H3p, 0, 0, 1, H2, HII, 0, 2, &this_all_rates->rate_239, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H3p, 0, 0, 1, H2p, HI, 0, 2, &this_all_rates->rate_240, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Cm, 0, 0, 1, CI, elec, 0, 2, &this_all_rates->rate_241, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, 0, 0, 1, CI, HI, 0, 2, &this_all_rates->rate_242, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, 0, 0, 1, CHp, elec, 0, 2, &this_all_rates->rate_243, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CHp, 0, 0, 1, CI, HII, 0, 2, &this_all_rates->rate_244, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH2, 0, 0, 1, CH, HI, 0, 2, &this_all_rates->rate_245, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2, 0, 0, 1, CH2p, elec, 0, 2, &this_all_rates->rate_246, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH2p, 0, 0, 1, CHp, HI, 0, 2, &this_all_rates->rate_247, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH3p, 0, 0, 1, CH2p, HI, 0, 2, &this_all_rates->rate_248, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH3p, 0, 0, 1, CHp, H2, 0, 2, &this_all_rates->rate_249, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, C2, 0, 0, 1, CI, CI, 0, 2, &this_all_rates->rate_250, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, Om, 0, 0, 1, OI, elec, 0, 2, &this_all_rates->rate_251, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, 0, 0, 1, OI, HI, 0, 2, &this_all_rates->rate_252, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, 0, 0, 1, OHp, elec, 0, 2, &this_all_rates->rate_253, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OHp, 0, 0, 1, OI, HII, 0, 2, &this_all_rates->rate_254, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2O, 0, 0, 1, OH, HI, 0, 2, &this_all_rates->rate_255, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2O, 0, 0, 1, H2Op, elec, 0, 2, &this_all_rates->rate_256, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2Op, 0, 0, 1, H2p, OI, 0, 2, &this_all_rates->rate_257, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2Op, 0, 0, 1, HII, OH, 0, 2, &this_all_rates->rate_258, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2Op, 0, 0, 1, OII, H2, 0, 2, &this_all_rates->rate_259, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2Op, 0, 0, 1, OHp, HI, 0, 2, &this_all_rates->rate_260, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, 0, 0, 1, HII, H2O, 0, 2, &this_all_rates->rate_261, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H3Op, 0, 0, 1, H2p, OH, 0, 2, &this_all_rates->rate_262, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H3Op, 0, 0, 1, H2Op, HI, 0, 2, &this_all_rates->rate_263, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H3Op, 0, 0, 1, OHp, H2, 0, 2, &this_all_rates->rate_264, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, 0, 0, 1, O2p, elec, 0, 2, &this_all_rates->rate_265, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, 0, 0, 1, OI, OI, 0, 2, &this_all_rates->rate_266, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CO, 0, 0, 1, CI, OI, 0, 2, &this_all_rates->rate_267, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, H2, 0, 0, 1, HII, HI, elec, 3, &this_all_rates->rate_268, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2, 0, 0, 1, HII, Hm, 0, 2, &this_all_rates->rate_269, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CO, 0, 0, 1, COp, elec, 0, 2, &this_all_rates->rate_270, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH, 0, 0, 1, CI, HI, 0, 2, &this_all_rates->rate_272, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CHp, 0, 0, 1, CII, HI, 0, 2, &this_all_rates->rate_273, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CH2, 0, 0, 1, CH2p, elec, 0, 2, &this_all_rates->rate_274, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CH2, 0, 0, 1, CH, HI, 0, 2, &this_all_rates->rate_275, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, C2, 0, 0, 1, CI, CI, 0, 2, &this_all_rates->rate_276, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, OH, 0, 0, 1, OI, HI, 0, 2, &this_all_rates->rate_277, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, H2O, 0, 0, 1, OH, HI, 0, 2, &this_all_rates->rate_278, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, 0, 0, 1, OI, OI, 0, 2, &this_all_rates->rate_279, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, O2, 0, 0, 1, O2p, elec, 0, 2, &this_all_rates->rate_280, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CO, 0, 0, 1, CI, OI, 0, 2, &this_all_rates->rate_281, incl_mol, myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, CII, OH, 0, 2, CO, HII, 0, 2, &this_all_rates->rate_283, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, MgII, elec, 0, 2, MgI, 0, 0, 1, &this_all_rates->rate_290, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction =
      add_new_reaction(last_reaction, SII, elec, 0, 2, SI, 0, 0, 1, &this_all_rates->rate_291, incl_mol, myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CaII, elec, 0, 2, CaI, 0, 0, 1, &this_all_rates->rate_292, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CaIII, elec, 0, 2, CaII, 0, 0, 1, &this_all_rates->rate_293, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, CII, MgI, 0, 2, CI, MgII, 0, 2, &this_all_rates->rate_294, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, NII, MgI, 0, 2, NI, MgII, 0, 2, &this_all_rates->rate_295, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, SiII, MgI, 0, 2, SiI, MgII, 0, 2, &this_all_rates->rate_296, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, SII, MgI, 0, 2, SI, MgII, 0, 2, &this_all_rates->rate_297, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, COp, H2, 0, 2, HCOp, HI, 0, 2, &this_all_rates->rate_303, incl_mol,
                                   myGlobalVars->speciesIndices);
  last_reaction = add_new_reaction(last_reaction, COp, H2, 0, 2, HOCp, HI, 0, 2, &this_all_rates->rate_304, incl_mol,
                                   myGlobalVars->speciesIndices);

  /* Now we add the reactions from Ben's model */
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      if(chimesRateTables.NonEqIon->IonIndexBegin[i] == HI)
        {
          last_reaction = add_new_reaction(last_reaction, HI, 0, 0, 1, HII, elec, 0, 2, &(this_all_rates->BensRates[i].cosmicRays[0]),
                                           incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, HI, 0, 0, 1, HII, elec, 0, 2, &(this_all_rates->BensRates[i].PhotoIon[0][0]),
                                           incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, Hm, 0, 0, 1, HI, elec, 0, 2, &(this_all_rates->BensRates[i].PhotoIon[1][0]),
                                           incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, H2, 0, 0, 1, H2p, elec, 0, 2, &(this_all_rates->BensRates[i].PhotoIon[2][0]),
                                           incl_mol, myGlobalVars->speciesIndices);
        }
      else if(chimesRateTables.NonEqIon->IonIndexBegin[i] == HeI)
        {
          last_reaction = add_new_reaction(last_reaction, HeII, elec, 0, 2, HeIII, elec, elec, 3,
                                           &(this_all_rates->BensRates[i].CollisIon[1]), incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, HeI, 0, 0, 1, HeII, elec, 0, 2,
                                           &(this_all_rates->BensRates[i].PhotoIon[0][0]), incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, HeII, 0, 0, 1, HeIII, elec, 0, 2,
                                           &(this_all_rates->BensRates[i].PhotoIon[1][0]), incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, HeIII, elec, 0, 2, HeII, 0, 0, 1, &(this_all_rates->BensRates[i].Recomb[1]),
                                           incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, HeI, 0, 0, 1, HeII, elec, 0, 2,
                                           &(this_all_rates->BensRates[i].cosmicRays[0]), incl_mol, myGlobalVars->speciesIndices);
          last_reaction = add_new_reaction(last_reaction, HeII, 0, 0, 1, HeIII, elec, 0, 2,
                                           &(this_all_rates->BensRates[i].cosmicRays[1]), incl_mol, myGlobalVars->speciesIndices);
          if(chimesRateTables.NonEqIon->NonEqRates[i].CTHrec_mask[2] == 1)
            last_reaction = add_new_reaction(last_reaction, HeIII, HI, 0, 2, HeII, HII, 0, 2,
                                             &(this_all_rates->BensRates[i].CTHrec[1]), incl_mol, myGlobalVars->speciesIndices);
          if(chimesRateTables.NonEqIon->NonEqRates[i].CTHion_mask[1] == 1)
            last_reaction = add_new_reaction(last_reaction, HeII, HII, 0, 2, HeIII, HI, 0, 2,
                                             &(this_all_rates->BensRates[i].CTHion[1]), incl_mol, myGlobalVars->speciesIndices);
        }
      else
        {
          for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i] - 1; j++)
            {
              last_reaction = add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), elec, 0, 2,
                                               (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), elec, elec, 3,
                                               &(this_all_rates->BensRates[i].CollisIon[j]), incl_mol, myGlobalVars->speciesIndices);
              last_reaction = add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), elec, 0, 2,
                                               (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), 0, 0, 1,
                                               &(this_all_rates->BensRates[i].Recomb[j]), incl_mol, myGlobalVars->speciesIndices);
              last_reaction = add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), 0, 0, 1,
                                               (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), elec, 0, 2,
                                               &(this_all_rates->BensRates[i].cosmicRays[j]), incl_mol, myGlobalVars->speciesIndices);

              /* Note that we only include CT for ions j <= 4*/
              if(j <= 4)
                {
                  if(chimesRateTables.NonEqIon->NonEqRates[i].CTHion_mask[j] == 1)
                    last_reaction =
                        add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), HII, 0, 2,
                                         (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), HI, 0, 2,
                                         &(this_all_rates->BensRates[i].CTHion[j]), incl_mol, myGlobalVars->speciesIndices);
                  if(chimesRateTables.NonEqIon->NonEqRates[i].CTHeion_mask[j] == 1)
                    last_reaction =
                        add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), HeII, 0, 2,
                                         (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), HeI, 0, 2,
                                         &(this_all_rates->BensRates[i].CTHeion[j]), incl_mol, myGlobalVars->speciesIndices);
                  if(chimesRateTables.NonEqIon->NonEqRates[i].CTHrec_mask[j + 1] == 1)
                    last_reaction =
                        add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), HI, 0, 2,
                                         (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), HII, 0, 2,
                                         &(this_all_rates->BensRates[i].CTHrec[j]), incl_mol, myGlobalVars->speciesIndices);
                  if(chimesRateTables.NonEqIon->NonEqRates[i].CTHerec_mask[j + 1] == 1)
                    last_reaction =
                        add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + 1), HeI, 0, 2,
                                         (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), HeII, 0, 2,
                                         &(this_all_rates->BensRates[i].CTHerec[j]), incl_mol, myGlobalVars->speciesIndices);
                }

              /* Note that for Auger ionisations, the 3rd product
               * index is set to -(no of EXTRA electrons). */
              for(k = 0; k < chimesRateTables.NonEqIon->N_Auger[i]; k++)
                {
                  if(chimesRateTables.NonEqIon->NonEqRates[i].sigmaphot[0][j][k] > 0.0)
                    last_reaction =
                        add_new_reaction(last_reaction, (chimesRateTables.NonEqIon->IonIndexBegin[i] + j), 0, 0, 1,
                                         (chimesRateTables.NonEqIon->IonIndexBegin[i] + j + k + 1), elec, -k, 2,
                                         &(this_all_rates->BensRates[i].PhotoIon[j][k]), incl_mol, myGlobalVars->speciesIndices);
                }
            }
        }
    }
}

int set_species_index_array(struct globalVariables *myGlobalVariables)
{
  int i;
  int current_index = 0;

  /* First, include electrons and all
   * H & He ions. */
  for(i = elec; i <= HeIII; i++)
    {
      myGlobalVariables->speciesIndices[i] = current_index;
      current_index += 1;
    }

  /* For each metal, check whether it
   * is included. If it is, add these
   * the speciesIndices, otherwise
   * make that index -1. */
  for(i = CI; i <= Cm; i++)
    {
      if(myGlobalVariables->element_included[0] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = NI; i <= NVIII; i++)
    {
      if(myGlobalVariables->element_included[1] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = OI; i <= Om; i++)
    {
      if(myGlobalVariables->element_included[2] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = NeI; i <= NeXI; i++)
    {
      if(myGlobalVariables->element_included[3] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = MgI; i <= MgXIII; i++)
    {
      if(myGlobalVariables->element_included[4] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = SiI; i <= SiXV; i++)
    {
      if(myGlobalVariables->element_included[5] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = SI; i <= SXVII; i++)
    {
      if(myGlobalVariables->element_included[6] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = CaI; i <= CaXXI; i++)
    {
      if(myGlobalVariables->element_included[7] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  for(i = FeI; i <= FeXXVII; i++)
    {
      if(myGlobalVariables->element_included[8] == 1)
        {
          myGlobalVariables->speciesIndices[i] = current_index;
          current_index += 1;
        }
      else
        myGlobalVariables->speciesIndices[i] = -1;
    }

  /* Determine which molecules are present. */
  for(i = H2; i <= H3p; i++)
    {
      myGlobalVariables->speciesIndices[i] = current_index;
      current_index += 1;
    }

  if(myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[OH] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[H2O] = current_index;
      current_index += 1;
    }
  else
    {
      myGlobalVariables->speciesIndices[OH]  = -1;
      myGlobalVariables->speciesIndices[H2O] = -1;
    }

  if(myGlobalVariables->element_included[0] == 1)
    {
      myGlobalVariables->speciesIndices[C2] = current_index;
      current_index += 1;
    }
  else
    myGlobalVariables->speciesIndices[C2] = -1;

  if(myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[O2] = current_index;
      current_index += 1;
    }
  else
    myGlobalVariables->speciesIndices[O2] = -1;

  if(myGlobalVariables->element_included[0] == 1 && myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[HCOp] = current_index;
      current_index += 1;
    }
  else
    myGlobalVariables->speciesIndices[HCOp] = -1;

  if(myGlobalVariables->element_included[0] == 1)
    {
      myGlobalVariables->speciesIndices[CH] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[CH2] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[CH3p] = current_index;
      current_index += 1;
    }
  else
    {
      myGlobalVariables->speciesIndices[CH]   = -1;
      myGlobalVariables->speciesIndices[CH2]  = -1;
      myGlobalVariables->speciesIndices[CH3p] = -1;
    }

  if(myGlobalVariables->element_included[0] == 1 && myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[CO] = current_index;
      current_index += 1;
    }
  else
    myGlobalVariables->speciesIndices[CO] = -1;

  if(myGlobalVariables->element_included[0] == 1)
    {
      myGlobalVariables->speciesIndices[CHp] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[CH2p] = current_index;
      current_index += 1;
    }
  else
    {
      myGlobalVariables->speciesIndices[CHp]  = -1;
      myGlobalVariables->speciesIndices[CH2p] = -1;
    }

  if(myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[OHp] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[H2Op] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[H3Op] = current_index;
      current_index += 1;
    }
  else
    {
      myGlobalVariables->speciesIndices[OHp]  = -1;
      myGlobalVariables->speciesIndices[H2Op] = -1;
      myGlobalVariables->speciesIndices[H3Op] = -1;
    }

  if(myGlobalVariables->element_included[0] == 1 && myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[COp] = current_index;
      current_index += 1;
      myGlobalVariables->speciesIndices[HOCp] = current_index;
      current_index += 1;
    }
  else
    {
      myGlobalVariables->speciesIndices[COp]  = -1;
      myGlobalVariables->speciesIndices[HOCp] = -1;
    }

  if(myGlobalVariables->element_included[2] == 1)
    {
      myGlobalVariables->speciesIndices[O2p] = current_index;
      current_index += 1;
    }
  else
    myGlobalVariables->speciesIndices[O2p] = -1;

  return current_index;
}

void init_chimes(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates,
                 struct Reactions_Structure **this_all_reactions_root, struct Reactions_Structure **this_nonmolecular_reactions_root,
                 double *dustG_arr, double *H2_dissocJ_arr)
{
  int i;

  myGlobalVars->totalNumberOfSpecies = set_species_index_array(myGlobalVars);

  /* Read in Bens tables and the photoionisation
   * tables. */
  initialise_bens_tables(myGlobalVars, this_all_rates, dustG_arr, H2_dissocJ_arr);

  /* Read in the tables of the additional rates. */
  initialise_additional_rates_tables(myGlobalVars);

  /* Read in tables of equilibrium abundances */
  GetEqAbundancesTables(myGlobalVars);

  /* Read in tables for the molecular, CI and OI
   * cooling. */
  initialise_cooling(myGlobalVars);

  /* Allocate memory to the reaction lists
   * and initialise them.*/
  *this_all_reactions_root = (struct Reactions_Structure *)mymalloc("Chimes_all_reactions", sizeof(struct Reactions_Structure));
  initialise_reactions(*this_all_reactions_root, 1, myGlobalVars, *this_all_rates);

  *this_nonmolecular_reactions_root =
      (struct Reactions_Structure *)mymalloc("Chimes_nonmolecular_reactions", sizeof(struct Reactions_Structure));
  initialise_reactions(*this_nonmolecular_reactions_root, 0, myGlobalVars, *this_all_rates);

  // BEGIN joki============================================================
  if(myGlobalVars->print_debug_statements == 1)
    {
      printf("PRINTING GLOBALVARS:\n");
      printf("myGlobalVars->reductionOn = %i \n ", myGlobalVars->reductionOn);
      printf("myGlobalVars->updatePhotonFluxOn = %i \n ", myGlobalVars->updatePhotonFluxOn);
      printf("myGlobalVars->cellSelfShieldingOn = %i \n ", myGlobalVars->cellSelfShieldingOn);
      printf("myGlobalVars->N_spectra = %i \n", myGlobalVars->N_spectra);
      printf("myGlobalVars->StaticMolCooling = %i \n ", myGlobalVars->StaticMolCooling);
      printf("myGlobalVars->T_EqThresh = %e \n ", myGlobalVars->T_EqThresh);
      printf("myGlobalVars->time_tolerance = %e \n ", myGlobalVars->time_tolerance);
      printf("myGlobalVars-> min_subcyclestep = %e \n ", myGlobalVars->min_subcyclestep);
      printf("myGlobalVars->T_mol = %e \n ", myGlobalVars->T_mol);
      printf("myGlobalVars->n_ions_low = %i \n ", myGlobalVars->n_ions_low);
      printf("myGlobalVars->n_ions_med = %i \n ", myGlobalVars->n_ions_med);
      printf("myGlobalVars->n_ions_high = %i \n ", myGlobalVars->n_ions_high);
      printf("myGlobalVars->InitIonState = %i \n ", myGlobalVars->InitIonState);
      printf("myGlobalVars->grain_temperature = %e \n ", myGlobalVars->grain_temperature);
      printf("myGlobalVars->cmb_temperature = %e \n ", myGlobalVars->cmb_temperature);
      printf("myGlobalVars->relativeTolerance = %e \n ", myGlobalVars->relativeTolerance);
      printf("myGlobalVars->absoluteTolerance = %e \n ", myGlobalVars->absoluteTolerance);
      printf("myGlobalVars->thermalAbsoluteTolerance = %e \n ", myGlobalVars->thermalAbsoluteTolerance);
      printf("myGlobalVars->totalNumberOfSpecies = %i \n ", myGlobalVars->totalNumberOfSpecies);
      for(i = 0; i < 8; i++)
        printf("myGlobalVars->element_included[%i] = %i \n ", i, myGlobalVars->element_included[i]);
      printf("==================================\n");
    }
  // END joki============================================================
}

void allocate_gas_abundances_memory(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  myGasVars->abundances =
      (double *)mymalloc_movable(&(myGasVars->abundances), "Chimes_abun", myGlobalVars->totalNumberOfSpecies * sizeof(double));
#ifdef CHIMES_ADVECT_ABUNDANCES
  myGasVars->IonAdvect =
      (double *)mymalloc_movable(&(myGasVars->IonAdvect), "Chimes_advect", myGlobalVars->totalNumberOfSpecies * sizeof(double));
#endif

  /* We also allocate memory for radiation fields here. */
  myGasVars->isotropic_photon_density =
      (double *)mymalloc_movable(&(myGasVars->isotropic_photon_density), "Chimes_iso_phot", myGlobalVars->N_spectra * sizeof(double));
  myGasVars->dust_G_parameter =
      (double *)mymalloc_movable(&(myGasVars->dust_G_parameter), "Chimes_dustG", myGlobalVars->N_spectra * sizeof(double));
  myGasVars->H2_dissocJ =
      (double *)mymalloc_movable(&(myGasVars->H2_dissocJ), "Chimes_H2_dissocJ", myGlobalVars->N_spectra * sizeof(double));
}

void free_gas_abundances_memory(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  myfree_movable(myGasVars->H2_dissocJ);
  myfree_movable(myGasVars->dust_G_parameter);
  myfree_movable(myGasVars->isotropic_photon_density);

#ifdef CHIMES_ADVECT_ABUNDANCES
  myfree_movable(myGasVars->IonAdvect);
#endif
  myfree_movable(myGasVars->abundances);
}

/* The following routine sets the initial
 * abundances of each species. */
void initialise_gas_abundances(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  int i, init_ion_state;
  /* mode determines in what state
   * we initially set the gas. */
  if(myGlobalVars->InitIonState >= 0 && myGlobalVars->InitIonState <= 26)
    init_ion_state = myGlobalVars->InitIonState;
  else
    {
      printf("WARNING: initialise_gas_abundances() mode not recognised. Assuming fully neutral\n");
      fflush(stdout);
      init_ion_state = 0;
    }

  /* First, set all abundances to zero */
  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    myGasVars->abundances[i] = 0.0;

  /* Now set the abundances of the initial
   * species.  */
  myGasVars->abundances[myGlobalVars->speciesIndices[min(HI + init_ion_state, HII)]] = 1.0;
  myGasVars->abundances[myGlobalVars->speciesIndices[elec]] = myGasVars->abundances[myGlobalVars->speciesIndices[HII]];

  myGasVars->abundances[myGlobalVars->speciesIndices[min(HeI + init_ion_state, HeIII)]] = myGasVars->element_abundances[0];
  for(i = 1; i <= 2; i++)
    myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[HeI + i]] * i;

  if(myGlobalVars->element_included[0] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(CI + init_ion_state, CVII)]] = myGasVars->element_abundances[1];
      for(i = 1; i <= 6; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[CI + i]] * i;
    }

  if(myGlobalVars->element_included[1] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(NI + init_ion_state, NVIII)]] = myGasVars->element_abundances[2];
      for(i = 1; i <= 7; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[NI + i]] * i;
    }

  if(myGlobalVars->element_included[2] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(OI + init_ion_state, OIX)]] = myGasVars->element_abundances[3];
      for(i = 1; i <= 8; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[OI + i]] * i;
    }

  if(myGlobalVars->element_included[3] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(NeI + init_ion_state, NeXI)]] = myGasVars->element_abundances[4];
      for(i = 1; i <= 10; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[NeI + i]] * i;
    }

  if(myGlobalVars->element_included[4] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(MgI + init_ion_state, MgXIII)]] = myGasVars->element_abundances[5];
      for(i = 1; i <= 12; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[MgI + i]] * i;
    }

  if(myGlobalVars->element_included[5] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(SiI + init_ion_state, SiXV)]] = myGasVars->element_abundances[6];
      for(i = 1; i <= 14; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[SiI + i]] * i;
    }

  if(myGlobalVars->element_included[6] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(SI + init_ion_state, SXVII)]] = myGasVars->element_abundances[7];
      for(i = 1; i <= 16; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[SI + i]] * i;
    }

  if(myGlobalVars->element_included[7] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(CaI + init_ion_state, CaXXI)]] = myGasVars->element_abundances[8];
      for(i = 1; i <= 20; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[CaI + i]] * i;
    }

  if(myGlobalVars->element_included[8] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[min(FeI + init_ion_state, FeXXVII)]] = myGasVars->element_abundances[9];
      for(i = 1; i <= 26; i++)
        myGasVars->abundances[myGlobalVars->speciesIndices[elec]] += myGasVars->abundances[myGlobalVars->speciesIndices[FeI + i]] * i;
    }

  /* Make sure that the abundances of any elements
   * not included are set to zero. */
  for(i = 0; i < 9; i++)
    if(myGlobalVars->element_included[i] != 1)
      myGasVars->element_abundances[i + 1] = 0.0;

  // Set constant heating rate to zero.
  myGasVars->constant_heating_rate = 0.0;

  // BEGIN joki============================================================
  if(myGlobalVars->print_debug_statements == 1)
    {
      printf("PRINTING MYGASVARS:\n");
      for(i = 0; i < 9; i++)
        printf("myGasVars->element_abundances[%i] = %e \n", i, myGasVars->element_abundances[i]);
      printf("myGasVars->nH_tot = %e \n ", myGasVars->nH_tot);
      printf("myGasVars->temperature = %e \n ", myGasVars->temperature);
      printf("myGasVars->TempFloor = %e \n ", myGasVars->TempFloor);
      printf("myGasVars->divVel = %e \n ", myGasVars->divVel);
      printf("myGasVars->doppler_broad = %e \n ", myGasVars->doppler_broad);
      for(i = 0; i < myGlobalVars->N_spectra; i++)
        printf("myGasVars->isotropic_photon_density_%i = %e \n ", i, myGasVars->isotropic_photon_density[i]);
      printf("myGasVars->metallicity = %e \n ", myGasVars->metallicity);
      printf("myGasVars->cell_size = %e \n ", myGasVars->cell_size);
      printf("myGasVars->hydro_timestep = %e \n ", myGasVars->hydro_timestep);
      printf("myGasVars->ForceEqOn = %i \n ", myGasVars->ForceEqOn);
      for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
        printf("myGasVars->abundances[%i]= %f \n ", i, myGasVars->abundances[i]);
      printf("==================================\n");
    }
  // END joki============================================================
}
