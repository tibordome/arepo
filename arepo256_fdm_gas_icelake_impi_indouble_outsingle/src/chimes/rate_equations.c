#include <cvode/cvode.h>
#include <kinsol/kinsol.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sys/types.h>
#include <time.h>
#include "../allvars.h"
#include "../proto.h"
#include "proto.h"

void check_constraint_equations(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  double x;
  int i;

  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    myGasVars->abundances[i] = max(myGasVars->abundances[i], 0.0);

  /* Helium */
  if(myGasVars->element_abundances[0] > 0.0)
    {
      x = myGasVars->abundances[myGlobalVars->speciesIndices[HeI]] + myGasVars->abundances[myGlobalVars->speciesIndices[HeII]] +
          myGasVars->abundances[myGlobalVars->speciesIndices[HeIII]];
      if(fabs((x - myGasVars->element_abundances[0]) / myGasVars->element_abundances[0]) > 0.01)
        {
          for(i = 0; i < 3; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[HeI + i]] *= myGasVars->element_abundances[0] / x;
        }
    }

  /* Nitrogen */
  if(myGasVars->element_abundances[2] > 0.0 && myGlobalVars->element_included[1] == 1)
    {
      x = 0.0;
      for(i = 0; i < 8; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[NI + i]];
      if(fabs((x - myGasVars->element_abundances[2]) / myGasVars->element_abundances[2]) > 0.01)
        {
          for(i = 0; i < 8; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[NI + i]] *= myGasVars->element_abundances[2] / x;
        }
    }

  /* Neon */
  if(myGasVars->element_abundances[4] > 0.0 && myGlobalVars->element_included[3] == 1)
    {
      x = 0.0;
      for(i = 0; i < 11; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[NeI + i]];
      if(fabs((x - myGasVars->element_abundances[4]) / myGasVars->element_abundances[4]) > 0.01)
        {
          for(i = 0; i < 11; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[NeI + i]] *= myGasVars->element_abundances[4] / x;
        }
    }

  /* Magnesium */
  if(myGasVars->element_abundances[5] > 0.0 && myGlobalVars->element_included[4] == 1)
    {
      x = 0.0;
      for(i = 0; i < 13; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[MgI + i]];
      if(fabs((x - myGasVars->element_abundances[5]) / myGasVars->element_abundances[5]) > 0.01)
        {
          for(i = 0; i < 13; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[MgI + i]] *= myGasVars->element_abundances[5] / x;
        }
    }

  /* Silicon */
  if(myGasVars->element_abundances[6] > 0.0 && myGlobalVars->element_included[5] == 1)
    {
      x = 0.0;
      for(i = 0; i < 15; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[SiI + i]];
      if(fabs((x - myGasVars->element_abundances[6]) / myGasVars->element_abundances[6]) > 0.01)
        {
          for(i = 0; i < 15; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[SiI + i]] *= myGasVars->element_abundances[6] / x;
        }
    }
  /* Sulphur */
  if(myGasVars->element_abundances[7] > 0.0 && myGlobalVars->element_included[6] == 1)
    {
      x = 0.0;
      for(i = 0; i < 17; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[SI + i]];
      if(fabs((x - myGasVars->element_abundances[7]) / myGasVars->element_abundances[7]) > 0.01)
        {
          for(i = 0; i < 17; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[SI + i]] *= myGasVars->element_abundances[7] / x;
        }
    }

  /* Calcium */
  if(myGasVars->element_abundances[8] > 0.0 && myGlobalVars->element_included[7] == 1)
    {
      x = 0.0;
      for(i = 0; i < 21; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[CaI + i]];
      if(fabs((x - myGasVars->element_abundances[8]) / myGasVars->element_abundances[8]) > 0.01)
        {
          for(i = 0; i < 21; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[CaI + i]] *= myGasVars->element_abundances[8] / x;
        }
    }

  /* Iron */
  if(myGasVars->element_abundances[9] > 0.0 && myGlobalVars->element_included[8] == 1)
    {
      x = 0.0;
      for(i = 0; i < 27; i++)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[FeI + i]];
      if(fabs((x - myGasVars->element_abundances[9]) / myGasVars->element_abundances[9]) > 0.01)
        {
          for(i = 0; i < 27; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[FeI + i]] *= myGasVars->element_abundances[9] / x;
        }
    }

  /* Carbon */
  if(myGasVars->element_abundances[1] > 0.0 && myGlobalVars->element_included[0] == 1)
    {
      x = 0.0;
      for(i = 0; i < 8; i++) /* Includes Cm */
        x += myGasVars->abundances[myGlobalVars->speciesIndices[CI + i]];
      x += 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[C2]] + myGasVars->abundances[myGlobalVars->speciesIndices[CH]] +
           myGasVars->abundances[myGlobalVars->speciesIndices[CH2]] + myGasVars->abundances[myGlobalVars->speciesIndices[CH3p]] +
           myGasVars->abundances[myGlobalVars->speciesIndices[CHp]] + myGasVars->abundances[myGlobalVars->speciesIndices[CH2p]];
      if(myGlobalVars->element_included[2] == 1)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] + myGasVars->abundances[myGlobalVars->speciesIndices[CO]] +
             myGasVars->abundances[myGlobalVars->speciesIndices[COp]] + myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]];
      if(fabs((x - myGasVars->element_abundances[1]) / myGasVars->element_abundances[1]) > 0.01)
        {
          for(i = 0; i < 8; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[CI + i]] *= myGasVars->element_abundances[1] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[C2]] *= myGasVars->element_abundances[1] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH]] *= myGasVars->element_abundances[1] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH2]] *= myGasVars->element_abundances[1] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH3p]] *= myGasVars->element_abundances[1] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CHp]] *= myGasVars->element_abundances[1] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH2p]] *= myGasVars->element_abundances[1] / x;
          if(myGlobalVars->element_included[2] == 1)
            {
              myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] *= myGasVars->element_abundances[1] / x;
              myGasVars->abundances[myGlobalVars->speciesIndices[CO]] *= myGasVars->element_abundances[1] / x;
              myGasVars->abundances[myGlobalVars->speciesIndices[COp]] *= myGasVars->element_abundances[1] / x;
              myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]] *= myGasVars->element_abundances[1] / x;
            }
        }
    }

  /* Oxygen */
  if(myGasVars->element_abundances[3] > 0.0 && myGlobalVars->element_included[2] == 1)
    {
      x = 0.0;
      for(i = 0; i < 10; i++) /* Includes Om */
        x += myGasVars->abundances[myGlobalVars->speciesIndices[OI + i]];
      x += myGasVars->abundances[myGlobalVars->speciesIndices[OH]] + myGasVars->abundances[myGlobalVars->speciesIndices[H2O]] +
           2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[O2]] + myGasVars->abundances[myGlobalVars->speciesIndices[OHp]] +
           myGasVars->abundances[myGlobalVars->speciesIndices[H2Op]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[O2p]] +
           myGasVars->abundances[myGlobalVars->speciesIndices[H3Op]];
      if(myGlobalVars->element_included[0] == 1)
        x += myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] + myGasVars->abundances[myGlobalVars->speciesIndices[CO]] +
             myGasVars->abundances[myGlobalVars->speciesIndices[COp]] + myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]];
      if(fabs((x - myGasVars->element_abundances[3]) / myGasVars->element_abundances[3]) > 0.01)
        {
          for(i = 0; i < 10; i++)
            myGasVars->abundances[myGlobalVars->speciesIndices[OI + i]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[OH]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[H2O]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[O2]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[OHp]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[H2Op]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[O2p]] *= myGasVars->element_abundances[3] / x;
          myGasVars->abundances[myGlobalVars->speciesIndices[H3Op]] *= myGasVars->element_abundances[3] / x;
          if(myGlobalVars->element_included[0] == 1)
            {
              myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] *= myGasVars->element_abundances[3] / x;
              myGasVars->abundances[myGlobalVars->speciesIndices[CO]] *= myGasVars->element_abundances[3] / x;
              myGasVars->abundances[myGlobalVars->speciesIndices[COp]] *= myGasVars->element_abundances[3] / x;
              myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]] *= myGasVars->element_abundances[3] / x;
            }
        }
    }

  /* Hydrogen */
  x = myGasVars->abundances[myGlobalVars->speciesIndices[HI]] + myGasVars->abundances[myGlobalVars->speciesIndices[HII]] +
      myGasVars->abundances[myGlobalVars->speciesIndices[Hm]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[H2]] +
      2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[H2p]] + 3.0 * myGasVars->abundances[myGlobalVars->speciesIndices[H3p]];

  if(myGlobalVars->element_included[0] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[CH]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[CH2]] +
         3.0 * myGasVars->abundances[myGlobalVars->speciesIndices[CH3p]] + myGasVars->abundances[myGlobalVars->speciesIndices[CHp]] +
         2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[CH2p]];
  if(myGlobalVars->element_included[2] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[OH]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[H2O]] +
         myGasVars->abundances[myGlobalVars->speciesIndices[OHp]] + 2.0 * myGasVars->abundances[myGlobalVars->speciesIndices[H2Op]] +
         3.0 * myGasVars->abundances[myGlobalVars->speciesIndices[H3Op]];
  if(myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] + myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]];

  if(fabs(x - 1.0) > 0.01)
    {
      myGasVars->abundances[myGlobalVars->speciesIndices[HI]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[HII]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[Hm]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[H2]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[H2p]] /= x;
      myGasVars->abundances[myGlobalVars->speciesIndices[H3p]] /= x;
      if(myGlobalVars->element_included[0] == 1)
        {
          myGasVars->abundances[myGlobalVars->speciesIndices[CH]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH2]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH3p]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CHp]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[CH2p]] /= x;
        }
      if(myGlobalVars->element_included[2] == 1)
        {
          myGasVars->abundances[myGlobalVars->speciesIndices[OH]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[H2O]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[OHp]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[H2Op]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[H3Op]] /= x;
        }
      if(myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
        {
          myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] /= x;
          myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]] /= x;
        }
    }

  /* Electrons */
  x = myGasVars->abundances[myGlobalVars->speciesIndices[HII]];
  x -= myGasVars->abundances[myGlobalVars->speciesIndices[Hm]];
  x += myGasVars->abundances[myGlobalVars->speciesIndices[HeII]];
  x += myGasVars->abundances[myGlobalVars->speciesIndices[HeIII]] * 2.0;

  if(myGlobalVars->element_included[0] == 1)
    {
      for(i = 1; i <= 6; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[CI + i]];
      x -= myGasVars->abundances[myGlobalVars->speciesIndices[Cm]];
      x += myGasVars->abundances[myGlobalVars->speciesIndices[CH3p]] + myGasVars->abundances[myGlobalVars->speciesIndices[CHp]] +
           myGasVars->abundances[myGlobalVars->speciesIndices[CH2p]];
    }

  if(myGlobalVars->element_included[1] == 1)
    {
      for(i = 1; i <= 7; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[NI + i]];
    }

  if(myGlobalVars->element_included[2] == 1)
    {
      for(i = 1; i <= 8; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[OI + i]];
      x -= myGasVars->abundances[myGlobalVars->speciesIndices[Om]];
      x += myGasVars->abundances[myGlobalVars->speciesIndices[OHp]] + myGasVars->abundances[myGlobalVars->speciesIndices[H2Op]] +
           myGasVars->abundances[myGlobalVars->speciesIndices[H3Op]] + myGasVars->abundances[myGlobalVars->speciesIndices[O2p]];
    }

  if(myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    x += myGasVars->abundances[myGlobalVars->speciesIndices[HCOp]] + +myGasVars->abundances[myGlobalVars->speciesIndices[COp]] +
         myGasVars->abundances[myGlobalVars->speciesIndices[HOCp]];

  if(myGlobalVars->element_included[3] == 1)
    {
      for(i = 1; i <= 10; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[NeI + i]];
    }

  if(myGlobalVars->element_included[4] == 1)
    {
      for(i = 1; i <= 12; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[MgI + i]];
    }

  if(myGlobalVars->element_included[5] == 1)
    {
      for(i = 1; i <= 14; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[SiI + i]];
    }

  if(myGlobalVars->element_included[6] == 1)
    {
      for(i = 1; i <= 16; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[SI + i]];
    }

  if(myGlobalVars->element_included[7] == 1)
    {
      for(i = 1; i <= 20; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[CaI + i]];
    }

  if(myGlobalVars->element_included[8] == 1)
    {
      for(i = 1; i <= 26; i++)
        x += i * myGasVars->abundances[myGlobalVars->speciesIndices[FeI + i]];
    }

  x += myGasVars->abundances[myGlobalVars->speciesIndices[H2p]] + myGasVars->abundances[myGlobalVars->speciesIndices[H3p]];

  if(fabs((x - myGasVars->abundances[myGlobalVars->speciesIndices[elec]]) /
          myGasVars->abundances[myGlobalVars->speciesIndices[elec]]) > 0.01)
    myGasVars->abundances[myGlobalVars->speciesIndices[elec]] = x;
}

void set_constraint_equations(void *user_data, double *output_array)
{
  /* This routine calculates the electron and
   * singly ionised atomic abundances from their
   * constraint equations and writes these to
   * the output array in the order elec, HI,
   * HeI etc.*/
  double x_temp;
  int i, output_index;
  UserData *data;
  data = (UserData *)user_data;

  output_index = 0;

  /* Electrons */
  x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HII]];
  x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[Hm]];
  x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeII]];
  x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeIII]] * 2.0;

  if(data->myGlobalVars->element_included[0] == 1 && data->myGasVars->element_abundances[1] > 0.0)
    {
      for(i = 1; i <= 6; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CI + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[Cm]];
      x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH3p]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CHp]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH2p]];
    }

  if(data->myGlobalVars->element_included[1] == 1 && data->myGasVars->element_abundances[2] > 0.0)
    {
      for(i = 1; i <= 7; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NI + i]];
    }

  if(data->myGlobalVars->element_included[2] == 1 && data->myGasVars->element_abundances[3] > 0.0)
    {
      for(i = 1; i <= 8; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OI + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[Om]];
      x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OHp]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2Op]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H3Op]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[O2p]];
    }

  if(data->myGlobalVars->element_included[0] == 1 && data->myGasVars->element_abundances[1] > 0.0 &&
     data->myGlobalVars->element_included[2] == 1 && data->myGasVars->element_abundances[3] > 0.0)
    x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HCOp]] +
              +data->myGasVars->abundances[data->myGlobalVars->speciesIndices[COp]] +
              data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HOCp]];

  if(data->myGlobalVars->element_included[3] == 1 && data->myGasVars->element_abundances[4] > 0.0)
    {
      for(i = 1; i <= 10; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NeI + i]];
    }

  if(data->myGlobalVars->element_included[4] == 1 && data->myGasVars->element_abundances[5] > 0.0)
    {
      for(i = 1; i <= 12; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[MgI + i]];
    }

  if(data->myGlobalVars->element_included[5] == 1 && data->myGasVars->element_abundances[6] > 0.0)
    {
      for(i = 1; i <= 14; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SiI + i]];
    }

  if(data->myGlobalVars->element_included[6] == 1 && data->myGasVars->element_abundances[7] > 0.0)
    {
      for(i = 1; i <= 16; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SI + i]];
    }

  if(data->myGlobalVars->element_included[7] == 1 && data->myGasVars->element_abundances[8] > 0.0)
    {
      for(i = 1; i <= 20; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CaI + i]];
    }

  if(data->myGlobalVars->element_included[8] == 1 && data->myGasVars->element_abundances[9] > 0.0)
    {
      for(i = 1; i <= 26; i++)
        x_temp += i * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[FeI + i]];
    }
  x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2p]] +
            data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H3p]];

  output_array[output_index] = x_temp;
  output_index += 1;

  /* Hydrogen */
  x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HI]] +
           data->myGasVars->abundances[data->myGlobalVars->speciesIndices[Hm]] +
           2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2]] +
           2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2p]] +
           3.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H3p]];

  if(data->myGlobalVars->element_included[0] == 1)
    x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH]] +
              2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH2]] +
              3.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH3p]] +
              data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CHp]] +
              2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH2p]];
  if(data->myGlobalVars->element_included[2] == 1)
    x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OH]] +
              2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2O]] +
              data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OHp]] +
              2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2Op]] +
              3.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H3Op]];
  if(data->myGlobalVars->element_included[0] == 1 && data->myGlobalVars->element_included[2] == 1)
    x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HCOp]] +
              data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HOCp]];

  output_array[output_index] = 1.0 - x_temp;
  output_index += 1;

  /* Helium */
  if(data->myGasVars->element_abundances[0] > 0.0)
    {
      output_array[output_index] =
          data->myGasVars->element_abundances[0] - (data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeI]] +
                                                    data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeIII]]);
      output_index += 1;
    }

  /* Carbon */
  if(data->myGasVars->element_abundances[1] > 0.0 && data->myGlobalVars->element_included[0] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CI]];
      for(i = 0; i < 6; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CII]];
      x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[Cm]];

      x_temp += 2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[C2]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH2]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH3p]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CHp]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CH2p]];
      if(data->myGlobalVars->element_included[2] == 1)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HCOp]] +
                  data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CO]] +
                  data->myGasVars->abundances[data->myGlobalVars->speciesIndices[COp]] +
                  data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HOCp]];

      output_array[output_index] = data->myGasVars->element_abundances[1] - x_temp;
      output_index += 1;
    }

  /* Nitrogen */
  if(data->myGasVars->element_abundances[2] > 0.0 && data->myGlobalVars->element_included[1] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NI]];
      for(i = 0; i < 7; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NII]];
      output_array[output_index] = data->myGasVars->element_abundances[2] - x_temp;
      output_index += 1;
    }

  /* Oxygen */
  if(data->myGasVars->element_abundances[3] > 0.0 && data->myGlobalVars->element_included[2] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OI]];
      for(i = 0; i < 8; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OII + i]];
      x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[Om]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OII]];

      x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OH]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2O]] +
                2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[O2]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OHp]] +
                data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2Op]] +
                2.0 * data->myGasVars->abundances[data->myGlobalVars->speciesIndices[O2p]];
      if(data->myGlobalVars->element_included[0] == 1)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HCOp]] +
                  data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CO]] +
                  data->myGasVars->abundances[data->myGlobalVars->speciesIndices[COp]] +
                  data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HOCp]];

      output_array[output_index] = data->myGasVars->element_abundances[3] - x_temp;
      output_index += 1;
    }

  /* Neon */
  if(data->myGasVars->element_abundances[4] > 0.0 && data->myGlobalVars->element_included[3] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NeI]];
      for(i = 0; i < 10; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NeII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[NeII]];
      output_array[output_index] = data->myGasVars->element_abundances[4] - x_temp;
      output_index += 1;
    }

  /* Magnesium */
  if(data->myGasVars->element_abundances[5] > 0.0 && data->myGlobalVars->element_included[4] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[MgI]];
      for(i = 0; i < 12; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[MgII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[MgII]];
      output_array[output_index] = data->myGasVars->element_abundances[5] - x_temp;
      output_index += 1;
    }

  /* Silicon */
  if(data->myGasVars->element_abundances[6] > 0.0 && data->myGlobalVars->element_included[5] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SiI]];
      for(i = 0; i < 14; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SiII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SiII]];
      output_array[output_index] = data->myGasVars->element_abundances[6] - x_temp;
      output_index += 1;
    }

  /* Sulphur */
  if(data->myGasVars->element_abundances[7] > 0.0 && data->myGlobalVars->element_included[6] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SI]];
      for(i = 0; i < 16; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[SII]];
      output_array[output_index] = data->myGasVars->element_abundances[7] - x_temp;
      output_index += 1;
    }

  /* Calcium */
  if(data->myGasVars->element_abundances[8] > 0.0 && data->myGlobalVars->element_included[7] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CaI]];
      for(i = 0; i < 20; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CaII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CaII]];
      output_array[output_index] = data->myGasVars->element_abundances[8] - x_temp;
      output_index += 1;
    }

  /* Iron */
  if(data->myGasVars->element_abundances[9] > 0.0 && data->myGlobalVars->element_included[8] == 1)
    {
      x_temp = data->myGasVars->abundances[data->myGlobalVars->speciesIndices[FeI]];
      for(i = 0; i < 26; i++)
        x_temp += data->myGasVars->abundances[data->myGlobalVars->speciesIndices[FeII + i]];
      x_temp -= data->myGasVars->abundances[data->myGlobalVars->speciesIndices[FeII]];
      output_array[output_index] = data->myGasVars->element_abundances[9] - x_temp;
      output_index += 1;
    }

  return;
}

void set_kin_scale_vector(N_Vector scale_vector, int *enum_indices, struct gasVariables *myGasVars,
                          struct globalVariables *myGlobalVars)
{
  int i;

  for(i = HeI; i <= HeIII; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[0];
        }
    }
  for(i = CI; i <= Cm; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[1];
        }
    }
  for(i = NI; i <= NVIII; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[2];
        }
    }
  for(i = OI; i <= OIX; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[3];
        }
    }
  for(i = NeI; i <= NeXI; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[4];
        }
    }
  for(i = MgI; i <= MgXIII; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[5];
        }
    }
  for(i = SiI; i <= SiXV; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[6];
        }
    }
  for(i = SI; i <= SXVII; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[7];
        }
    }
  for(i = CaI; i <= CaXXI; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[8];
        }
    }
  for(i = FeI; i <= FeXXVII; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[i]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[i]]) = 1.0 / myGasVars->element_abundances[9];
        }
    }
  if(myGasVars->temperature <= myGlobalVars->T_mol)
    {
      if(myGlobalVars->speciesIndices[OH] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[OH]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[OH]]) = 1.0 / myGasVars->element_abundances[3];
        }
      if(myGlobalVars->speciesIndices[H2O] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[H2O]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[H2O]]) = 1.0 / myGasVars->element_abundances[3];
        }
      if(myGlobalVars->speciesIndices[C2] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[C2]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[C2]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[O2] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[O2]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[O2]]) = 1.0 / myGasVars->element_abundances[3];
        }
      if(myGlobalVars->speciesIndices[HCOp] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[HCOp]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[HCOp]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[CH] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[CH]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[CH]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[CH2] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[CH2]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[CH2]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[CH3p] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[CH3p]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[CH3p]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[CO] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[CO]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[CO]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[CHp] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[CHp]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[CHp]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[CH2p] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[CH2p]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[CH2p]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[OHp] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[OHp]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[OHp]]) = 1.0 / myGasVars->element_abundances[3];
        }
      if(myGlobalVars->speciesIndices[H2Op] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[H2Op]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[H2Op]]) = 1.0 / myGasVars->element_abundances[3];
        }
      if(myGlobalVars->speciesIndices[H3Op] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[H3Op]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[H3Op]]) = 1.0 / myGasVars->element_abundances[3];
        }
      if(myGlobalVars->speciesIndices[COp] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[COp]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[COp]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[HOCp] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[HOCp]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[HOCp]]) = 1.0 / myGasVars->element_abundances[1];
        }
      if(myGlobalVars->speciesIndices[O2p] > -1)
        {
          if(enum_indices[myGlobalVars->speciesIndices[O2p]] > -1)
            NV_Ith_S(scale_vector, enum_indices[myGlobalVars->speciesIndices[O2p]]) = 1.0 / myGasVars->element_abundances[3];
        }
    }

  return;
}

int kin_f(N_Vector y, N_Vector ydot, void *user_data)
{
  double current_rate;
  int i, j, constr_index;
  UserData *data;

  data = (UserData *)user_data;
  int *indices; /* We will use this array to relate the enum types of
                 * each (non-eq) species to their position in y */
  double constraint_abundances[12];
  int *enum_indices;

  indices      = (int *)mymalloc("Chimes_indices", data->network_size * sizeof(int));
  enum_indices = (int *)mymalloc("Chimes_enum", data->myGlobalVars->totalNumberOfSpecies * sizeof(int));

  /* First, loop through the enum types of all
   * non-eq species. If they are included in
   * the network then their abundance is in
   * the vector y. */

  i = 0; /* We use this to keep track of where we are in the vector y */
  for(j = 0; j < data->myGlobalVars->totalNumberOfSpecies; j++)
    {
      if(data->species[j].include_species == 1)
        {
          data->myGasVars->abundances[j] = NV_Ith_S(y, i);
          indices[i]                     = j;
          enum_indices[j]                = i;
          i++;
        }
      else
        enum_indices[j] = -1;
    }

  /* Constraint equations */
  set_constraint_equations(data, constraint_abundances);

  /* If cell self-shielding is on, update
   * column densities */

  if(data->myGlobalVars->cellSelfShieldingOn == 1)
    {
      *(data->HI_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HI]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      *(data->H2_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      *(data->HeI_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeI]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      *(data->HeII_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeII]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      if(data->myGlobalVars->speciesIndices[CO] > -1)
        *(data->CO_column) = max(data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CO]], 0.0) *
                             data->myGasVars->cell_size * data->myGasVars->nH_tot;
      else
        *(data->CO_column) = 0.0;
      if(data->myGlobalVars->speciesIndices[H2O] > -1)
        *(data->H2O_column) = max(data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2O]], 0.0) *
                              data->myGasVars->cell_size * data->myGasVars->nH_tot;
      else
        *(data->H2O_column) = 0.0;
      if(data->myGlobalVars->speciesIndices[OH] > -1)
        *(data->OH_column) = max(data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OH]], 0.0) *
                             data->myGasVars->cell_size * data->myGasVars->nH_tot;
      else
        *(data->OH_column) = 0.0;
      *(data->extinction) = DUSTEFFSIZE * data->myGasVars->cell_size * data->myGasVars->nH_tot * data->myGasVars->metallicity;
    }

  /* Update the rate coefficients that depend
   * on other variables, e.g. abundances */
  update_rates(data->myGasVars, data->myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
               *(data->CO_column), *(data->extinction), data->this_all_rates);

  /* Zero all species rates */
  for(i = 0; i < data->network_size; i++)
    {
      data->species[indices[i]].creation_rate    = 0.0;
      data->species[indices[i]].destruction_rate = 0.0;
    }

  /* Go through the list of reactions, updating the creation
   * and destruction rates of the non-eq species as we go */

  struct Reactions_Structure *current_reaction;
  current_reaction = data->root_reactions;

  while(current_reaction->next_reaction != NULL)
    {
      current_reaction = current_reaction->next_reaction;
      if(current_reaction->flag_included == 1)
        {
          current_rate = *(current_reaction->rate);
          for(i = 0; i < current_reaction->no_of_reactants; i++)
            current_rate *= data->myGasVars->abundances[current_reaction->reactants[i]];

          for(i = 0; i < current_reaction->no_of_reactants; i++)
            data->species[current_reaction->reactants[i]].destruction_rate += current_rate;
          for(i = 0; i < current_reaction->no_of_products; i++)
            data->species[current_reaction->products[i]].creation_rate += current_rate;

          /* For Auger ionisations, we do not have enough
           * room to fit in up to 11 products in the reactions
           * structure. Therefore, for these the third product
           * index is -(N_auger - 1), i.e. the extra number of
           * electrons that are produced. Therefore, we now
           * need to add these to the production rate of e- */
          if(current_reaction->extra_Auger_electrons > 0)
            data->species[data->myGlobalVars->speciesIndices[elec]].creation_rate +=
                current_rate * current_reaction->extra_Auger_electrons;
        }
    }

  /* Now set the output ydot vector for the chemical abundances */
  for(i = 0; i < data->network_size; i++)
    NV_Ith_S(ydot, i) = data->species[indices[i]].creation_rate - data->species[indices[i]].destruction_rate;

  /* Finally, use the constraint equations */
  constr_index = 0;
  NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[elec]]) =
      NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[elec]]) - constraint_abundances[constr_index++];
  NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[HII]]) =
      NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[HII]]) - constraint_abundances[constr_index++];
  NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[HeII]]) =
      NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[HeII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[1] > 0.0 && data->myGlobalVars->element_included[0] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[CII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[CII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[2] > 0.0 && data->myGlobalVars->element_included[1] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[NII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[NII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[3] > 0.0 && data->myGlobalVars->element_included[2] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[OII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[OII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[4] > 0.0 && data->myGlobalVars->element_included[3] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[NeII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[NeII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[5] > 0.0 && data->myGlobalVars->element_included[4] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[MgII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[MgII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[6] > 0.0 && data->myGlobalVars->element_included[5] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[SiII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[SiII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[7] > 0.0 && data->myGlobalVars->element_included[6] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[SII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[SII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[8] > 0.0 && data->myGlobalVars->element_included[7] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[CaII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[CaII]]) - constraint_abundances[constr_index++];
  if(data->myGasVars->element_abundances[9] > 0.0 && data->myGlobalVars->element_included[8] == 1)
    NV_Ith_S(ydot, enum_indices[data->myGlobalVars->speciesIndices[FeII]]) =
        NV_Ith_S(y, enum_indices[data->myGlobalVars->speciesIndices[FeII]]) - constraint_abundances[constr_index++];

  myfree(enum_indices);
  myfree(indices);

  return 0;
}

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  double current_rate;
  int i, j;
  UserData *data;

  data = (UserData *)user_data;
  int indices[data->network_size]; /* We will use this array to relate the enum types of
                                    * each (non-eq) species to their position in y */

  /* First, loop through the enum types of all
   * non-eq species. If they are included in
   * the network then their abundance is in
   * the vector y. */

  i = 0; /* We use this to keep track of where we are in the vector y */
  for(j = 0; j < data->myGlobalVars->totalNumberOfSpecies; j++)
    {
      if(data->species[j].include_species == 1)
        {
          data->myGasVars->abundances[j] = NV_Ith_S(y, i);
          indices[i]                     = j;
          i++;
        }
    }

  /* If Thermal Evolution is switched on, the final element in the
   * vector y is the internal energy (per unit volume). Use this
   * to update the temperature, and also the rates that depend on T */
  if(data->myGasVars->ThermEvolOn == 1)
    {
      data->myGasVars->temperature =
          max(NV_Ith_S(y, data->network_size) /
                  (1.5 * calculate_total_number_density(data->myGasVars->abundances, data->myGasVars->nH_tot, data->myGlobalVars) *
                   BOLTZMANNCGS),
              10.1); /* The rates are not defined below ~10 K */
      update_T_dependent_rates(data->myGasVars, data->myGlobalVars, data->this_all_rates);
    }

  /* If cell self-shielding is on, update
   * column densities */

  if(data->myGlobalVars->cellSelfShieldingOn == 1)
    {
      *(data->HI_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HI]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      *(data->H2_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      *(data->HeI_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeI]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      *(data->HeII_column) =
          data->myGasVars->abundances[data->myGlobalVars->speciesIndices[HeII]] * data->myGasVars->cell_size * data->myGasVars->nH_tot;
      if(data->myGlobalVars->speciesIndices[CO] > -1)
        *(data->CO_column) = max(data->myGasVars->abundances[data->myGlobalVars->speciesIndices[CO]], 0.0) *
                             data->myGasVars->cell_size * data->myGasVars->nH_tot;
      else
        *(data->CO_column) = 0.0;
      if(data->myGlobalVars->speciesIndices[H2O] > -1)
        *(data->H2O_column) = max(data->myGasVars->abundances[data->myGlobalVars->speciesIndices[H2O]], 0.0) *
                              data->myGasVars->cell_size * data->myGasVars->nH_tot;
      else
        *(data->H2O_column) = 0.0;
      if(data->myGlobalVars->speciesIndices[OH] > -1)
        *(data->OH_column) = max(data->myGasVars->abundances[data->myGlobalVars->speciesIndices[OH]], 0.0) *
                             data->myGasVars->cell_size * data->myGasVars->nH_tot;
      else
        *(data->OH_column) = 0.0;
      *(data->extinction) = DUSTEFFSIZE * data->myGasVars->cell_size * data->myGasVars->nH_tot * data->myGasVars->metallicity;
    }

  /* Update the rate coefficients that depend
   * on other variables, e.g. abundances */
  update_rates(data->myGasVars, data->myGlobalVars, *(data->HI_column), *(data->H2_column), *(data->HeI_column), *(data->HeII_column),
               *(data->CO_column), *(data->extinction), data->this_all_rates);

  /* Zero all species rates */
  for(i = 0; i < data->network_size; i++)
    {
      data->species[indices[i]].creation_rate    = 0.0;
      data->species[indices[i]].destruction_rate = 0.0;
    }

  /* Go through the list of reactions, updating the creation
   * and destruction rates of the non-eq species as we go */

  struct Reactions_Structure *current_reaction;
  current_reaction = data->root_reactions;

  while(current_reaction->next_reaction != NULL)
    {
      current_reaction = current_reaction->next_reaction;
      if(current_reaction->flag_included == 1)
        {
          current_rate = *(current_reaction->rate);
          for(i = 0; i < current_reaction->no_of_reactants; i++)
            current_rate *= data->myGasVars->abundances[current_reaction->reactants[i]];

          for(i = 0; i < current_reaction->no_of_reactants; i++)
            data->species[current_reaction->reactants[i]].destruction_rate += current_rate;
          for(i = 0; i < current_reaction->no_of_products; i++)
            data->species[current_reaction->products[i]].creation_rate += current_rate;

          /* For Auger ionisations, we do not have enough
           * room to fit in up to 11 products in the reactions
           * structure. Therefore, for these the third product
           * index is -(N_auger - 1), i.e. the extra number of
           * electrons that are produced. Therefore, we now
           * need to add these to the production rate of e- */
          if(current_reaction->extra_Auger_electrons > 0)
            data->species[data->myGlobalVars->speciesIndices[elec]].creation_rate +=
                current_rate * current_reaction->extra_Auger_electrons;
        }
    }

  /* Now set the output ydot vector for the chemical abundances */
  for(i = 0; i < data->network_size; i++)
    NV_Ith_S(ydot, i) = data->species[indices[i]].creation_rate - data->species[indices[i]].destruction_rate;

  /* Finally, if Thermal Evolution is switched on, calculate the cooling rate */
  if(data->myGasVars->ThermEvolOn == 1)
    {
      if(data->myGasVars->temperature > data->myGasVars->TempFloor)
        NV_Ith_S(ydot, data->network_size) = -calculate_total_cooling_rate(
            data->myGasVars, data->myGlobalVars, *(data->HI_column), *(data->HeI_column), *(data->HeII_column), *(data->H2_column),
            *(data->CO_column), *(data->H2O_column), *(data->OH_column), *(data->extinction),
            data->this_all_rates); /* Note that network_size is the number of chemcial species, hence No. of eqns = network_size + 1
                                      when ThermEvol is on */
      else
        NV_Ith_S(ydot, data->network_size) =
            max(-calculate_total_cooling_rate(data->myGasVars, data->myGlobalVars, *(data->HI_column), *(data->HeI_column),
                                              *(data->HeII_column), *(data->H2_column), *(data->CO_column), *(data->H2O_column),
                                              *(data->OH_column), *(data->extinction), data->this_all_rates),
                0.0); /* Once T falls below T_floor, set T_dot >= 0 */
    }

  return 0;
}
