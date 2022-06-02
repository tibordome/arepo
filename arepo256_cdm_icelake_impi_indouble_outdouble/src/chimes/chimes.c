#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_spgmr.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sys/types.h>
#include <time.h>
#include "../allvars.h"
#include "../proto.h"
#include "proto.h"

#define MAXSTEPS 1e5
#define FTOL 1.0e-13
#define STOL 1.0e-25

void set_equilibrium_abundances(void *user_data)
{
  /* This is used when ForceEqOn == 2 */
  realtype reltol, abstol_scalar;
  void *cvode_mem;
  int flag, i, j, thermEvolFlag;
  int *enum_indices;
  realtype t;
  N_Vector constr_vector, scale_vector, y;
  double dt_hydro;
  void *k_mem;
  struct gasVariables *myGasVars;
  struct globalVariables *myGlobalVars;
  struct Species_Structure *species;
  double network_size;
  UserData *data;
  data = (UserData *)user_data;

  myGasVars    = data->myGasVars;
  myGlobalVars = data->myGlobalVars;
  species      = data->species;
  network_size = data->network_size;

  double init_abundances[myGlobalVars->totalNumberOfSpecies];

  /* Set up solution, constraint & scale vectors */
  y             = N_VNew_Serial(network_size);
  constr_vector = N_VNew_Serial(network_size);
  scale_vector  = N_VNew_Serial(network_size);
  N_VConst_Serial(1.0, constr_vector);
  N_VConst_Serial(1.0, scale_vector);

  /* Set up KINSol solver */
  k_mem = KINCreate();
  flag  = KINInit(k_mem, kin_f, y);
  flag  = KINSetUserData(k_mem, data);
  flag  = KINSetConstraints(k_mem, constr_vector);
  flag  = KINSetFuncNormTol(k_mem, FTOL);
  flag  = KINSetScaledStepTol(k_mem, STOL);
  flag  = KINDense(k_mem, network_size);
  flag  = KINSetMaxSetupCalls(k_mem, 1);
  flag  = KINSetNumMaxIters(k_mem, 10000);

  /* Set initial guess */
  enum_indices = (int *)mymalloc_movable(&enum_indices, "Chimes_enum",
                                         myGlobalVars->totalNumberOfSpecies * sizeof(int)); /* Gives position in y for given enum. */
  i            = 0;
  for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
    {
      if(species[j].include_species == 1)
        {
          NV_Ith_S(y, i)  = myGasVars->abundances[j];
          enum_indices[j] = i;
          i++;
        }
      else
        enum_indices[j] = -1;
    }

  /* Scale each species by its corresponding element abundance. */
  set_kin_scale_vector(scale_vector, enum_indices, myGasVars, myGlobalVars);
  myfree_movable(enum_indices);

  for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
    init_abundances[j] = myGasVars->abundances[j];

  flag = KINSol(k_mem, y, KIN_NONE, scale_vector, scale_vector);

  /* If necessary, update the column densities. */
  if(myGlobalVars->cellSelfShieldingOn == 1)
    {
      *(data->HI_column)   = myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * myGasVars->cell_size * myGasVars->nH_tot;
      *(data->H2_column)   = myGasVars->abundances[myGlobalVars->speciesIndices[H2]] * myGasVars->cell_size * myGasVars->nH_tot;
      *(data->HeI_column)  = myGasVars->abundances[myGlobalVars->speciesIndices[HeI]] * myGasVars->cell_size * myGasVars->nH_tot;
      *(data->HeII_column) = myGasVars->abundances[myGlobalVars->speciesIndices[HeII]] * myGasVars->cell_size * myGasVars->nH_tot;
      if(myGlobalVars->speciesIndices[CO] > -1)
        *(data->CO_column) =
            max(myGasVars->abundances[myGlobalVars->speciesIndices[CO]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else
        *(data->CO_column) = 0.0;
      if(myGlobalVars->speciesIndices[H2O] > -1)
        *(data->H2O_column) =
            max(myGasVars->abundances[myGlobalVars->speciesIndices[H2O]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else
        *(data->H2O_column) = 0.0;
      if(myGlobalVars->speciesIndices[OH] > -1)
        *(data->OH_column) =
            max(myGasVars->abundances[myGlobalVars->speciesIndices[OH]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else
        *(data->OH_column) = 0.0;
      *(data->extinction) = DUSTEFFSIZE * myGasVars->cell_size * myGasVars->nH_tot * myGasVars->metallicity;
    }

  /* If flag == 2, KINSol may have had problems with
   * the initial conditions being too far from equilibrium.
   * Try evolving the abundances from InitIonState for
   * 10,000 years and then redo KINSol. */
  /* If flag == 1, the initial guess was OK. Still, integrate
   * for 10 Myr to make sure. */
  if(flag == 1 || flag == 2)
    {
      /* Start by integrating init abundances for
       * 10 kyr, followed by a call to
       * KINsol to finally get it into equilibrium.
       * Then repeat this with a factor 10 greater
       * integration time until the kinsol flag is
       * no longer 2. */

      thermEvolFlag             = myGasVars->ThermEvolOn;
      myGasVars->ThermEvolOn    = 0;
      dt_hydro                  = myGasVars->hydro_timestep;
      myGasVars->hydro_timestep = 3.17e11; /* 10 kyr ; */
      do
        {
          /* Update the y vector with the abundances before KINSol was called. */
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            myGasVars->abundances[j] = init_abundances[j];

          i = 0;
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            {
              if(species[j].include_species == 1)
                {
                  NV_Ith_S(y, i) = myGasVars->abundances[j];
                  i++;
                }
            }

          /* Update all rates. */
          set_constant_rates(myGasVars, myGlobalVars, data->this_all_rates);
          update_rates(myGasVars, myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
                       *(data->CO_column), *(data->extinction), data->this_all_rates);
          update_T_dependent_rates(myGasVars, myGlobalVars, data->this_all_rates);

          /* Set up CVode */
          reltol          = myGlobalVars->relativeTolerance;
          abstol_scalar   = myGlobalVars->absoluteTolerance;
          cvode_mem       = CVodeCreate(CV_BDF, CV_NEWTON);
          data->cvode_mem = cvode_mem;
          flag            = CVodeSetUserData(cvode_mem, data);
          flag            = CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);
          flag            = CVodeInit(cvode_mem, f, 0.0, y);
          flag            = CVodeSStolerances(cvode_mem, reltol, abstol_scalar);
          flag            = CVDense(cvode_mem, network_size);
          flag            = CVodeSetMaxConvFails(cvode_mem, 5000);

          /* Integrate for dt_hydro. */
          flag = CVode(cvode_mem, myGasVars->hydro_timestep, y, &t, CV_NORMAL);

          /* Update abundance array with solution vector. */
          i = 0;
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            {
              if(species[j].include_species == 1)
                {
                  myGasVars->abundances[j] = NV_Ith_S(y, i);
                  i++;
                }
            }

          check_constraint_equations(myGasVars, myGlobalVars);

          /* Some abundances may have been updated by the constraint
           * equations. Update y accordingly. */
          i = 0;
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            {
              if(species[j].include_species == 1)
                {
                  NV_Ith_S(y, i) = myGasVars->abundances[j];
                  i++;
                }
            }

          /* Record abundances in case KINSol crashes again. */
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            init_abundances[j] = myGasVars->abundances[j];

          /* Call KINSol to try and find equilibrium. */
          flag = KINSol(k_mem, y, KIN_NONE, scale_vector, scale_vector);

          /* If necessary, update the column densities. */
          if(myGlobalVars->cellSelfShieldingOn == 1)
            {
              *(data->HI_column) = myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * myGasVars->cell_size * myGasVars->nH_tot;
              *(data->H2_column) = myGasVars->abundances[myGlobalVars->speciesIndices[H2]] * myGasVars->cell_size * myGasVars->nH_tot;
              *(data->HeI_column) =
                  myGasVars->abundances[myGlobalVars->speciesIndices[HeI]] * myGasVars->cell_size * myGasVars->nH_tot;
              *(data->HeII_column) =
                  myGasVars->abundances[myGlobalVars->speciesIndices[HeII]] * myGasVars->cell_size * myGasVars->nH_tot;
              if(myGlobalVars->speciesIndices[CO] > -1)
                *(data->CO_column) =
                    max(myGasVars->abundances[myGlobalVars->speciesIndices[CO]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
              else
                *(data->CO_column) = 0.0;
              if(myGlobalVars->speciesIndices[H2O] > -1)
                *(data->H2O_column) =
                    max(myGasVars->abundances[myGlobalVars->speciesIndices[H2O]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
              else
                *(data->H2O_column) = 0.0;
              if(myGlobalVars->speciesIndices[OH] > -1)
                *(data->OH_column) =
                    max(myGasVars->abundances[myGlobalVars->speciesIndices[OH]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
              else
                *(data->OH_column) = 0.0;
              *(data->extinction) = DUSTEFFSIZE * myGasVars->cell_size * myGasVars->nH_tot * myGasVars->metallicity;
            }

          /* If we need to repeat these steps, increase dt_hydro by a factor of 10. */
          myGasVars->hydro_timestep *= 10.0;

          CVodeFree(&cvode_mem);
          data->cvode_mem = NULL;
          /* Quit once dt_hydro > 1000 Gyr */
          if(myGasVars->hydro_timestep > 3.16e19)
            break;
        }
      while(flag == 2); /* For flag == 1, only integrate once. */

      /* Reset ThermEvol flag & dt_hydro. */
      myGasVars->ThermEvolOn    = thermEvolFlag;
      myGasVars->hydro_timestep = dt_hydro;
    }

  if(flag == 2)
    {
      /* Set abundances to what they were after integrating for 1000 Gyr. */
      for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
        myGasVars->abundances[j] = init_abundances[j];
    }
  else
    {
      i = 0;
      for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
        {
          if(species[j].include_species == 1)
            {
              myGasVars->abundances[j] = NV_Ith_S(y, i);
              i++;
            }
        }
    }

  N_VDestroy_Serial(y);
  N_VDestroy_Serial(scale_vector);
  N_VDestroy_Serial(constr_vector);
  KINFree(&k_mem);
  return;
}

void set_equilibrium_abundances_from_tables(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  /* This is used when ForceEqOn == 1 */
  int i, j, k, l;
  double dT, dnH, dZ;

  if(myGasVars->temperature <= pow(10.0, EquilibriumAbundances.Temperatures[0]))
    {
      i  = 0;
      dT = 0.0;
    }
  else
    get_index_1d_mydbl(EquilibriumAbundances.Temperatures, EquilibriumAbundances.N_Temperatures, log10(myGasVars->temperature), &i,
                       &dT);

  if(myGasVars->nH_tot <= pow(10.0, EquilibriumAbundances.Densities[0]))
    {
      j   = 0;
      dnH = 0.0;
    }
  else
    get_index_1d_mydbl(EquilibriumAbundances.Densities, EquilibriumAbundances.N_Densities, log10(myGasVars->nH_tot), &j, &dnH);

  if(myGasVars->metallicity <= pow(10.0, EquilibriumAbundances.Metallicities[0]))
    {
      k  = 0;
      dZ = 0.0;
    }
  else
    get_index_1d_mydbl(EquilibriumAbundances.Metallicities, EquilibriumAbundances.N_Metallicities, log10(myGasVars->metallicity), &k,
                       &dZ);

  /* Note that the equilibrium tables tabulate
   * ionisation (or molecular) fraction, and
   * NOT the abundance wrt H. Now we need to
   * multiply by the appropriate element abundance. */
  for(l = myGlobalVars->speciesIndices[elec]; l <= myGlobalVars->speciesIndices[Hm]; l++)
    myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ));
  for(l = myGlobalVars->speciesIndices[HeI]; l <= myGlobalVars->speciesIndices[HeIII]; l++)
    myGasVars->abundances[l] =
        pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) * myGasVars->element_abundances[0];
  if(myGlobalVars->element_included[0] == 1)
    for(l = myGlobalVars->speciesIndices[CI]; l <= myGlobalVars->speciesIndices[Cm]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[1];
  if(myGlobalVars->element_included[1] == 1)
    for(l = myGlobalVars->speciesIndices[NI]; l <= myGlobalVars->speciesIndices[NVIII]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[2];
  if(myGlobalVars->element_included[2] == 1)
    for(l = myGlobalVars->speciesIndices[OI]; l <= myGlobalVars->speciesIndices[Om]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[3];
  if(myGlobalVars->element_included[3] == 1)
    for(l = myGlobalVars->speciesIndices[NeI]; l <= myGlobalVars->speciesIndices[NeXI]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[4];
  if(myGlobalVars->element_included[4] == 1)
    for(l = myGlobalVars->speciesIndices[MgI]; l <= myGlobalVars->speciesIndices[MgXIII]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[5];
  if(myGlobalVars->element_included[5] == 1)
    for(l = myGlobalVars->speciesIndices[SiI]; l <= myGlobalVars->speciesIndices[SiXV]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[6];
  if(myGlobalVars->element_included[6] == 1)
    for(l = myGlobalVars->speciesIndices[SI]; l <= myGlobalVars->speciesIndices[SXVII]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[7];
  if(myGlobalVars->element_included[7] == 1)
    for(l = myGlobalVars->speciesIndices[CaI]; l <= myGlobalVars->speciesIndices[CaXXI]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[8];
  if(myGlobalVars->element_included[8] == 1)
    for(l = myGlobalVars->speciesIndices[FeI]; l <= myGlobalVars->speciesIndices[FeXXVII]; l++)
      myGasVars->abundances[l] = pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, l, dT, dnH, dZ)) *
                                 myGasVars->element_abundances[9];

  myGasVars->abundances[myGlobalVars->speciesIndices[H2]] =
      pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[H2], dT, dnH, dZ));
  myGasVars->abundances[myGlobalVars->speciesIndices[H2p]] =
      pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[H2p], dT, dnH, dZ));
  myGasVars->abundances[myGlobalVars->speciesIndices[H3p]] =
      pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[H3p], dT, dnH, dZ));

  if(myGlobalVars->element_included[0] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[C2]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[C2], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[CH]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[CH], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[CH2]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[CH2], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[CH3p]] =
          pow(10.0,
              interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[CH3p], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[CHp]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[CHp], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[CH2p]] =
          pow(10.0,
              interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[CH2p], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
    }
  if(myGlobalVars->element_included[2] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[OH]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[OH], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
      myGasVars->abundances[myGlobalVars->speciesIndices[H2O]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[H2O], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
      myGasVars->abundances[myGlobalVars->speciesIndices[O2]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[O2], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
      myGasVars->abundances[myGlobalVars->speciesIndices[OHp]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[OHp], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
      myGasVars->abundances[myGlobalVars->speciesIndices[H2Op]] =
          pow(10.0,
              interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[H2Op], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
      myGasVars->abundances[myGlobalVars->speciesIndices[H3Op]] =
          pow(10.0,
              interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[H3Op], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
      myGasVars->abundances[myGlobalVars->speciesIndices[O2p]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[O2p], dT, dnH, dZ)) *
          myGasVars->element_abundances[3];
    }
  if(myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] =
          pow(10.0,
              interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[HCOp], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[CO]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[CO], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[COp]] =
          pow(10.0, interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[COp], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
      myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]] =
          pow(10.0,
              interpol_3d_special(EquilibriumAbundances.EqAbundances, i, j, k, myGlobalVars->speciesIndices[HOCp], dT, dnH, dZ)) *
          myGasVars->element_abundances[1];
    }

  return;
}

void chimes_network(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars,
                    struct All_rate_variables_structure *this_all_rates, struct Reactions_Structure *this_all_reactions_root,
                    struct Reactions_Structure *this_nonmolecular_reactions_root)
{
  realtype reltol, abstol_scalar, t;
  N_Vector abstol_vector, y;
  void *cvode_mem;

  int flag, i, j, old_network_size;
  double t_current_substep, dt_subcycle;
  double HI_column_density, H2_column_density, HeI_column_density, HeII_column_density, CO_column_density, H2O_column_density,
      OH_column_density, extinction;
  double internal_energy;
  int network_size, total_network_size, nonmolecular_network_size;
  struct Reactions_Structure *root_node_reaction_list;
  struct Species_Structure species[myGlobalVars->totalNumberOfSpecies];
  UserData data;

  if(myGlobalVars->cellSelfShieldingOn > 0)
    {
      /* if cellSelfShieldingOn == 1, these column densities will
       * subsequently be updated throughout the course of the
       * integration. If cellSelfShieldingOn == 2, we will set
       * these column densities here but they will NOT subsequently
       * be updated. */
      HI_column_density   = myGasVars->abundances[myGlobalVars->speciesIndices[HI]] * myGasVars->cell_size * myGasVars->nH_tot;
      H2_column_density   = myGasVars->abundances[myGlobalVars->speciesIndices[H2]] * myGasVars->cell_size * myGasVars->nH_tot;
      HeI_column_density  = myGasVars->abundances[myGlobalVars->speciesIndices[HeI]] * myGasVars->cell_size * myGasVars->nH_tot;
      HeII_column_density = myGasVars->abundances[myGlobalVars->speciesIndices[HeII]] * myGasVars->cell_size * myGasVars->nH_tot;
      if(myGlobalVars->speciesIndices[CO] > -1)
        CO_column_density =
            max(myGasVars->abundances[myGlobalVars->speciesIndices[CO]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else
        CO_column_density = 0.0;
      if(myGlobalVars->speciesIndices[H2O] > -1)
        H2O_column_density =
            max(myGasVars->abundances[myGlobalVars->speciesIndices[H2O]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else
        H2O_column_density = 0.0;
      if(myGlobalVars->speciesIndices[OH] > -1)
        OH_column_density =
            max(myGasVars->abundances[myGlobalVars->speciesIndices[OH]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else
        OH_column_density = 0.0;
      extinction = DUSTEFFSIZE * myGasVars->cell_size * myGasVars->nH_tot * myGasVars->metallicity;
    }
  else
    {
      HI_column_density   = 0.0;
      H2_column_density   = 0.0;
      HeI_column_density  = 0.0;
      HeII_column_density = 0.0;
      CO_column_density   = 0.0;
      H2O_column_density  = 0.0;
      OH_column_density   = 0.0;
      extinction          = 0.0;
    }

  set_constant_rates(myGasVars, myGlobalVars, this_all_rates);
  update_rates(myGasVars, myGlobalVars, HI_column_density, H2_column_density, HeI_column_density, HeII_column_density,
               CO_column_density, extinction, this_all_rates);
  update_T_dependent_rates(myGasVars, myGlobalVars, this_all_rates);

  set_species_arrays(species, myGasVars, &total_network_size, &nonmolecular_network_size, myGlobalVars);
  if(myGasVars->temperature <= myGlobalVars->T_mol)
    {
      root_node_reaction_list = this_all_reactions_root;
      network_size            = total_network_size;
    }
  else
    {
      /* Exclude all molecular species and set their abundances to zero. */
      for(i = H2; i <= O2p; i++)
        {
          if(myGlobalVars->speciesIndices[i] > -1)
            {
              species[myGlobalVars->speciesIndices[i]].include_species = 0;
              myGasVars->abundances[myGlobalVars->speciesIndices[i]]   = 0.0;
            }
        }
      check_constraint_equations(myGasVars, myGlobalVars);
      root_node_reaction_list = this_nonmolecular_reactions_root;
      network_size            = nonmolecular_network_size;
    }

  /* Set up structure to pass user data to the solver. */
  data.myGasVars      = myGasVars;
  data.myGlobalVars   = myGlobalVars;
  data.species        = species;
  data.root_reactions = root_node_reaction_list;
  data.this_all_rates = this_all_rates;
  data.HI_column      = &HI_column_density;
  data.H2_column      = &H2_column_density;
  data.HeI_column     = &HeI_column_density;
  data.HeII_column    = &HeII_column_density;
  data.CO_column      = &CO_column_density;
  data.H2O_column     = &H2O_column_density;
  data.OH_column      = &OH_column_density;
  data.extinction     = &extinction;
  data.network_size   = network_size;

  if(myGasVars->ForceEqOn > 0)
    {
      if(myGasVars->ThermEvolOn == 0)
        {
          if(myGasVars->ForceEqOn == 1)
            set_equilibrium_abundances_from_tables(myGasVars, myGlobalVars);
          else
            set_equilibrium_abundances(&data);
        }
      else
        do_equilibrium_cooling(&data);

      return;
    }

  if(myGlobalVars->reductionOn == 1)
    dt_subcycle = evaluate_reduced_network_size(species, myGasVars, myGlobalVars, root_node_reaction_list, &network_size,
                                                this_all_rates, this_all_reactions_root, this_nonmolecular_reactions_root);

  /* Create a serial vector of length network_size for the initial conditions. */
  if(myGasVars->ThermEvolOn == 0)
    {
      y = N_VNew_Serial(network_size);
      if(myGlobalVars->scale_metal_tolerances == 1)
        abstol_vector = N_VNew_Serial(network_size);
    }
  else
    {
      y               = N_VNew_Serial(network_size + 1);
      internal_energy = myGasVars->temperature * 1.5 *
                        calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS;
      NV_Ith_S(y, network_size)             = internal_energy;
      abstol_vector                         = N_VNew_Serial(network_size + 1);
      NV_Ith_S(abstol_vector, network_size) = myGlobalVars->thermalAbsoluteTolerance;
    }

  i = 0;
  for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
    {
      if(species[j].include_species == 1)
        {
          NV_Ith_S(y, i) = myGasVars->abundances[j];

          if(myGlobalVars->scale_metal_tolerances == 1)
            NV_Ith_S(abstol_vector, i) = myGlobalVars->absoluteTolerance * species[j].element_abundance;
          else if(myGasVars->ThermEvolOn == 1)
            NV_Ith_S(abstol_vector, i) = myGlobalVars->absoluteTolerance;

          i++;
        }
    }

  /* Set up the solver */
  /* Set the tolerances */
  reltol        = myGlobalVars->relativeTolerance;
  abstol_scalar = myGlobalVars->absoluteTolerance;

  /* Use CVodeCreate to create the solver
   * memory and specify the Backward Differentiation
   * Formula and Newton iteration. */
  cvode_mem      = CVodeCreate(CV_BDF, CV_NEWTON);
  data.cvode_mem = cvode_mem;

  /* Set the user data for CVode */
  flag = CVodeSetUserData(cvode_mem, &data);

  /* Use CVodeSetMaxNumSteps to set the maximum number
   * of steps CVode takes. */
  flag = CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);

  /* Use CVodeInit to initialise the integrator
   * memory and specify the right hand side
   * function in y' = f(t,y) (i.e. the rate
   * equations), the initial time 0.0 and the
   * initial conditions, in y. */
  flag = CVodeInit(cvode_mem, f, 0.0, y);

  /* Use CVodeSVtolerances to specify the scalar relative and absolute tolerances. */
  if(myGasVars->ThermEvolOn == 0 && myGlobalVars->scale_metal_tolerances == 0)
    flag = CVodeSStolerances(cvode_mem, reltol, abstol_scalar);
  else
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vector);

  /* Use CVDense to specify the CVDENSE dense linear solver. */
  if(myGasVars->ThermEvolOn == 0)
    flag = CVDense(cvode_mem, network_size);
  else
    flag = CVDense(cvode_mem, network_size + 1);

  /* Specify the maximum number of convergence test failures. */
  flag = CVodeSetMaxConvFails(cvode_mem, 5000);

  if(myGlobalVars->reductionOn == 1 && dt_subcycle < myGasVars->hydro_timestep)
    {
      /* This means that we will have to re-evaluate
       * the reduced network more than once before
       * the end of the hydro time step */
      t_current_substep = 0.0;
      while(t_current_substep < myGasVars->hydro_timestep)
        {
          dt_subcycle = min(dt_subcycle, myGasVars->hydro_timestep - t_current_substep);
          dt_subcycle = max(dt_subcycle, myGlobalVars->min_subcyclestep);

          flag = CVode(cvode_mem, dt_subcycle, y, &t, CV_NORMAL);

          /* Write the output abundances to the gas cell .
           * Note that species not included in the reduced
           * network are kept constant in the GasVars struct. */
          i = 0;
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            {
              if(species[j].include_species == 1)
                {
                  myGasVars->abundances[j] = NV_Ith_S(y, i);
                  i++;
                }
            }

          check_constraint_equations(myGasVars, myGlobalVars);

          if(myGasVars->ThermEvolOn == 1)
            {
              myGasVars->temperature = max(
                  NV_Ith_S(y, network_size) /
                      (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS),
                  myGasVars->TempFloor);
              /* If T has reached TempFloor, ensure that the thermal
               * energy is put back onto this floor */
              if(myGasVars->temperature == myGasVars->TempFloor)
                NV_Ith_S(y, network_size) = myGasVars->temperature * 1.5 *
                                            calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) *
                                            BOLTZMANNCGS;
            }

          old_network_size = network_size; /* This variable is used to see when the size of the network changes. */
          if(myGlobalVars->reductionOn == 1)
            dt_subcycle = evaluate_reduced_network_size(species, myGasVars, myGlobalVars, root_node_reaction_list, &network_size,
                                                        this_all_rates, this_all_reactions_root, this_nonmolecular_reactions_root);

          data.root_reactions = root_node_reaction_list;

          /* If the network size has changed we will need
           * to recreate y */
          if(old_network_size != network_size)
            {
              N_VDestroy_Serial(y);
              if(myGasVars->ThermEvolOn == 0)
                {
                  y = N_VNew_Serial(network_size);
                  if(myGlobalVars->scale_metal_tolerances == 1)
                    {
                      N_VDestroy_Serial(abstol_vector);
                      abstol_vector = N_VNew_Serial(network_size);
                    }
                }
              else
                {
                  N_VDestroy_Serial(abstol_vector);
                  y                         = N_VNew_Serial(network_size + 1);
                  NV_Ith_S(y, network_size) = myGasVars->temperature * 1.5 *
                                              calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) *
                                              BOLTZMANNCGS;
                  abstol_vector                         = N_VNew_Serial(network_size + 1);
                  NV_Ith_S(abstol_vector, network_size) = myGlobalVars->thermalAbsoluteTolerance;
                }

              i = 0;
              for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
                {
                  if(species[j].include_species == 1)
                    {
                      if(myGlobalVars->scale_metal_tolerances == 1)
                        NV_Ith_S(abstol_vector, i) = myGlobalVars->absoluteTolerance * species[j].element_abundance;
                      else if(myGasVars->ThermEvolOn == 1)
                        NV_Ith_S(abstol_vector, i) = myGlobalVars->absoluteTolerance;

                      i++;
                    }
                }
            }

          /* The check_constrain_equations routine
           * may have changed the abundances -
           * need to update the vector y */
          i = 0;
          for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
            {
              if(species[j].include_species == 1)
                {
                  NV_Ith_S(y, i) = myGasVars->abundances[j];
                  i++;
                }
            }

          /* Now update data.network_size & reset the CVode memory */
          data.network_size = network_size;

          CVodeFree(&cvode_mem);
          cvode_mem      = CVodeCreate(CV_BDF, CV_NEWTON);
          data.cvode_mem = cvode_mem;
          flag           = CVodeSetUserData(cvode_mem, &data);
          flag           = CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);
          flag           = CVodeInit(cvode_mem, f, 0.0, y);
          if(myGasVars->ThermEvolOn == 0)
            {
              if(myGlobalVars->scale_metal_tolerances == 0)
                flag = CVodeSStolerances(cvode_mem, reltol, abstol_scalar);
              else
                flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vector);
              flag = CVDense(cvode_mem, network_size);
            }
          else
            {
              flag = CVodeSVtolerances(cvode_mem, reltol, abstol_vector);
              flag = CVDense(cvode_mem, network_size + 1);
            }
          flag = CVodeSetMaxConvFails(cvode_mem, 5000);
          t_current_substep += t;
        }
    }
  else
    {
      /* Either reduction is off, or the subcyclestep
       * is long enough that we can just integrate
       * to the end of the hydro timestep, without
       * needing to re-evaluate the reduced network. */
      flag = CVode(cvode_mem, myGasVars->hydro_timestep, y, &t, CV_NORMAL);

      /* Write the output abundances to the gas cell
       * Note that species not included in the reduced
       * network are kept constant in the GasVars struct. */
      i = 0;
      for(j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
        {
          if(species[j].include_species == 1)
            {
              myGasVars->abundances[j] = NV_Ith_S(y, i);
              i++;
            }
        }

      check_constraint_equations(myGasVars, myGlobalVars);

      if(myGasVars->ThermEvolOn == 1)
        {
          myGasVars->temperature =
              max(NV_Ith_S(y, network_size) /
                      (1.5 * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS),
                  myGasVars->TempFloor);
        }
    }

  N_VDestroy_Serial(y);
  if(myGasVars->ThermEvolOn == 1 || myGlobalVars->scale_metal_tolerances == 1)
    N_VDestroy_Serial(abstol_vector);
  CVodeFree(&cvode_mem);

  return;
}
