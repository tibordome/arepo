#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include "proto.h"

double determine_subcyclestep(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, int *abundant_indices,
                              struct Species_Structure *species, struct Reactions_Structure *reactions_root,
                              struct All_rate_variables_structure *this_all_rates)
{
  double dt_new, dt_temp;
  double HI_column_density, HeI_column_density, HeII_column_density, H2_column_density, CO_column_density, H2O_column_density,
      OH_column_density, extinction;
  int i;
  /* Note that the abundant_indices array contains the indices from
   * the enumerated list - you will need to convert these to the
   * indices that only run up to totalNumberOfSpecies, that are
   * used in e.g. the abundances and species arrays. */

  dt_new = 1.0e30;

  update_species_fluxes(myGasVars->abundances, myGasVars->temperature, reactions_root, species, myGlobalVars->T_mol, 1,
                        myGlobalVars->speciesIndices, myGlobalVars);

  if(myGlobalVars->cellSelfShieldingOn == 1)
    {
      /* WARNING: This routine is NOT consistent with using
       * cellSelfShieldingOn == 2. If you intend to use this
       * again, you will need to give this routine the
       * column densities from the user_data structure, as they
       * cannot subsequently be re-computed during the course of
       * the integration. */
      HI_column_density   = myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * myGasVars->cell_size * myGasVars->nH_tot;
      HeI_column_density  = myGasVars->abundances[myGlobalVars->speciesIndices[HeI]] * myGasVars->cell_size * myGasVars->nH_tot;
      HeII_column_density = myGasVars->abundances[myGlobalVars->speciesIndices[HeII]] * myGasVars->cell_size * myGasVars->nH_tot;
      H2_column_density   = myGasVars->abundances[myGlobalVars->speciesIndices[H2]] * myGasVars->cell_size * myGasVars->nH_tot;
      CO_column_density = max(myGasVars->abundances[myGlobalVars->speciesIndices[CO]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      H2O_column_density =
          max(myGasVars->abundances[myGlobalVars->speciesIndices[H2O]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      OH_column_density = max(myGasVars->abundances[myGlobalVars->speciesIndices[OH]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      extinction        = DUSTEFFSIZE * myGasVars->cell_size * myGasVars->nH_tot * myGasVars->metallicity;
    }
  else
    {
      HI_column_density   = 0.0;
      HeI_column_density  = 0.0;
      HeII_column_density = 0.0;
      H2_column_density   = 0.0;
      CO_column_density   = 0.0;
      H2O_column_density  = 0.0;
      OH_column_density   = 0.0;
      extinction          = 0.0;
    }

  /* From rate of change of maximum abundances */
  for(i = 0; i < 9; i++)
    {
      if(abundant_indices[i] != -1)
        {
          dt_temp = myGlobalVars->time_tolerance * fabs(myGasVars->abundances[myGlobalVars->speciesIndices[abundant_indices[i]]] /
                                                        (species[myGlobalVars->speciesIndices[abundant_indices[i]]].creation_rate -
                                                         species[myGlobalVars->speciesIndices[abundant_indices[i]]].destruction_rate));
          if(dt_temp < dt_new)
            dt_new = dt_temp;
        }
    }

  /* From rate of change of T */
  /* This gives the estimated time
   * until T will cross the T_mol
   * threshold. */
  if(myGasVars->ThermEvolOn == 1)
    {
      dt_temp = (myGlobalVars->T_mol - myGasVars->temperature) * 1.5 *
                calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS /
                (-calculate_total_cooling_rate(myGasVars, myGlobalVars, HI_column_density, HeI_column_density, HeII_column_density,
                                               H2_column_density, CO_column_density, H2O_column_density, OH_column_density, extinction,
                                               this_all_rates));
      if(dt_temp > myGlobalVars->min_subcyclestep && dt_temp < dt_new)
        dt_new = dt_temp;
    }

  return max(dt_new, myGlobalVars->min_subcyclestep);
}

int is_x_in_array(int x, int *my_array, int len_array)
{
  int j;

  for(j = 0; j < len_array; j++)
    {
      if(x == my_array[j])
        return 1;
    }

  return 0;
}

void update_species_fluxes(double *abundances, double temperature, struct Reactions_Structure *reactions_root,
                           struct Species_Structure *species, double Tmol, int mode, int *speciesIndices,
                           struct globalVariables *myGlobalVars)
{
  int i;
  /* Note that this updates the creation and destruction
   * rates in the species array using the FULL network,
   * and NOT just the reduced network. */
  struct Reactions_Structure *current_reaction;
  double current_rate;

  current_reaction = reactions_root;

  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    {
      species[i].creation_rate    = 0.0;
      species[i].destruction_rate = 0.0;
    }

  while(current_reaction->next_reaction != NULL)
    {
      current_reaction = current_reaction->next_reaction;
      if(current_reaction->flag_included == 1 || mode == 1)
        {
          current_rate = *(current_reaction->rate);
          for(i = 0; i < current_reaction->no_of_reactants; i++)
            current_rate *= abundances[current_reaction->reactants[i]];
          for(i = 0; i < current_reaction->no_of_reactants; i++)
            species[current_reaction->reactants[i]].destruction_rate += current_rate;
          for(i = 0; i < current_reaction->no_of_products; i++)
            species[current_reaction->products[i]].creation_rate += current_rate;
          /* Auger Ionisations */
          if(current_reaction->extra_Auger_electrons > 0)
            species[speciesIndices[elec]].creation_rate += current_rate * current_reaction->extra_Auger_electrons;
        }
    }
}

void determine_included_reactions(struct Reactions_Structure *reactions_root, struct Species_Structure *species, int *speciesIndices)
{
  int incl_flag, i;
  struct Reactions_Structure *current_reaction;
  current_reaction = reactions_root;

  while(current_reaction->next_reaction != NULL)
    {
      current_reaction = current_reaction->next_reaction;
      incl_flag        = 1;
      for(i = 0; i < current_reaction->no_of_reactants; i++)
        {
          if(species[current_reaction->reactants[i]].include_species == 0)
            {
              incl_flag = 0;
              break;
            }
        }
      for(i = 0; i < current_reaction->no_of_products; i++)
        {
          if(species[current_reaction->products[i]].include_species == 0)
            {
              incl_flag = 0;
              break;
            }
        }
      current_reaction->flag_included = incl_flag;
    }
}

double evaluate_reduced_network_size(struct Species_Structure *species, struct gasVariables *myGasVars,
                                     struct globalVariables *myGlobalVars, struct Reactions_Structure *reaction_list_node,
                                     int *myNetworkSize, struct All_rate_variables_structure *this_all_rates,
                                     struct Reactions_Structure *this_all_reactions_root,
                                     struct Reactions_Structure *this_nonmolecular_reactions_root)
{
  int i, abundant_index;
  int max_indices[9];
  double dt;
  double O_abund = myGasVars->element_abundances[3];
  double C_abund = myGasVars->element_abundances[1];
  int *speciesIndices;
  speciesIndices = myGlobalVars->speciesIndices;

  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    species[i].include_species = 0;

  /* First, determine whether molecules are included */
  if(myGasVars->temperature < myGlobalVars->T_mol)
    {
      species[speciesIndices[H2]].include_species   = 1;
      species[speciesIndices[H2p]].include_species  = 1;
      species[speciesIndices[H3p]].include_species  = 1;
      species[speciesIndices[OH]].include_species   = (O_abund > 0.0);
      species[speciesIndices[H2O]].include_species  = (O_abund > 0.0);
      species[speciesIndices[C2]].include_species   = (C_abund > 0.0);
      species[speciesIndices[O2]].include_species   = (O_abund > 0.0);
      species[speciesIndices[HCOp]].include_species = (O_abund > 0.0 && C_abund > 0.0);
      species[speciesIndices[CH]].include_species   = (C_abund > 0.0);
      species[speciesIndices[CH2]].include_species  = (C_abund > 0.0);
      species[speciesIndices[CH3p]].include_species = (C_abund > 0.0);
      species[speciesIndices[CO]].include_species   = (O_abund > 0.0 && C_abund > 0.0);
      species[speciesIndices[CHp]].include_species  = (C_abund > 0.0);
      species[speciesIndices[CH2p]].include_species = (C_abund > 0.0);
      species[speciesIndices[OHp]].include_species  = (O_abund > 0.0);
      species[speciesIndices[H2Op]].include_species = (O_abund > 0.0);
      species[speciesIndices[H3Op]].include_species = (O_abund > 0.0);
      species[speciesIndices[COp]].include_species  = (O_abund > 0.0 && C_abund > 0.0);
      species[speciesIndices[HOCp]].include_species = (O_abund > 0.0 && C_abund > 0.0);
      species[speciesIndices[O2p]].include_species  = (O_abund > 0.0);
      species[speciesIndices[Cm]].include_species   = (C_abund > 0.0);
      species[speciesIndices[Om]].include_species   = (O_abund > 0.0);

      reaction_list_node = this_all_reactions_root;
      *myNetworkSize     = 3;

      if(O_abund > 0.0 && C_abund > 0.0)
        *myNetworkSize += 19;
      else if(O_abund > 0.0)
        *myNetworkSize += 8;
      else if(C_abund > 0.0)
        *myNetworkSize += 7;
    }
  else
    {
      reaction_list_node = this_nonmolecular_reactions_root;
      *myNetworkSize     = 0;
    }

  /* Next, add electrons and all H and He ions*/
  for(i = 0; i < 7; i++)
    {
      species[speciesIndices[i]].include_species = 1;
      *myNetworkSize += 1;
    }

  /* Now go through each element, adding in
   * the most abundant ion +/- n_ion. We will
   * start everything in max_indices at -1, so
   * if an element is not included in the total
   * network its max_index will be left at -1. */

  for(i = 0; i < 9; i++)
    max_indices[i] = -1;

  /* Carbon */
  if(myGasVars->element_abundances[1] > 0.0)
    {
      abundant_index = CI;
      for(i = CII; i < Cm; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[0] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_low); i <= (abundant_index + myGlobalVars->n_ions_low); i++)
        {
          if(i >= CI && i <= CVII)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Nitrogen */
  if(myGasVars->element_abundances[2] > 0.0)
    {
      abundant_index = NI;
      for(i = NII; i <= NVIII; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[1] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_low); i <= (abundant_index + myGlobalVars->n_ions_low); i++)
        {
          if(i >= NI && i <= NVIII)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Oxygen */
  if(myGasVars->element_abundances[3] > 0.0)
    {
      abundant_index = OI;
      for(i = OII; i < Om; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[2] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_low); i <= (abundant_index + myGlobalVars->n_ions_low); i++)
        {
          if(i >= OI && i <= OIX)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Neon */
  if(myGasVars->element_abundances[4] > 0.0)
    {
      abundant_index = NeI;
      for(i = NeII; i <= NeXI; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[3] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_med); i <= (abundant_index + myGlobalVars->n_ions_med); i++)
        {
          if(i >= NeI && i <= NeXI)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Magnesium */
  if(myGasVars->element_abundances[5] > 0.0)
    {
      abundant_index = MgI;
      for(i = MgII; i <= MgXIII; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[4] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_med); i <= (abundant_index + myGlobalVars->n_ions_med); i++)
        {
          if(i >= MgI && i <= MgXIII)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Silicon */
  if(myGasVars->element_abundances[6] > 0.0)
    {
      abundant_index = SiI;
      for(i = SiII; i < SiXV; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[5] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_med); i <= (abundant_index + myGlobalVars->n_ions_med); i++)
        {
          if(i >= SiI && i <= SiXV)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Sulphur */
  if(myGasVars->element_abundances[7] > 0.0)
    {
      abundant_index = SI;
      for(i = SII; i <= SXVII; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[6] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_high); i <= (abundant_index + myGlobalVars->n_ions_high); i++)
        {
          if(i >= SI && i <= SXVII)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Calcium */
  if(myGasVars->element_abundances[8] > 0.0)
    {
      abundant_index = CaI;
      for(i = CaII; i <= CaXXI; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[7] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_high); i <= (abundant_index + myGlobalVars->n_ions_high); i++)
        {
          if(i >= CaI && i <= CaXXI)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  /* Iron */
  if(myGasVars->element_abundances[9] > 0.0)
    {
      abundant_index = FeI;
      for(i = FeII; i <= FeXXVII; i++)
        {
          if(myGasVars->abundances[speciesIndices[i]] > myGasVars->abundances[speciesIndices[abundant_index]])
            abundant_index = i;
        }
      max_indices[8] = abundant_index;
      for(i = (abundant_index - myGlobalVars->n_ions_high); i <= (abundant_index + myGlobalVars->n_ions_high); i++)
        {
          if(i >= FeI && i <= FeXXVII)
            {
              species[speciesIndices[i]].include_species = 1;
              *myNetworkSize += 1;
            }
        }
    }

  determine_included_reactions(reaction_list_node, species, myGlobalVars->speciesIndices);

  dt = determine_subcyclestep(myGasVars, myGlobalVars, max_indices, species, reaction_list_node, this_all_rates);

  return dt;
}
