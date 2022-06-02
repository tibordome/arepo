/* This file contains functions to
 * set up the reaction and
 * photoionisation rates at the
 * beginning of the run, based
 * on the given parameters.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "proto.h"

void set_constant_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars,
                        struct All_rate_variables_structure *this_all_rates)
{
  /* Note that this function also includes all the
   * rates in update_T_dependent_rates, but NOT those
   * in update_rates */
  int i, j;
  double nH_tot  = myGasVars->nH_tot;
  double cr_rate = myGasVars->cr_rate;

  this_all_rates->rate_4   = nH_tot * k4();
  this_all_rates->rate_44  = nH_tot * k44();
  this_all_rates->rate_79  = nH_tot * k79();
  this_all_rates->rate_80  = nH_tot * k80();
  this_all_rates->rate_70  = R70(cr_rate);
  this_all_rates->rate_89  = R89(cr_rate);
  this_all_rates->rate_100 = nH_tot * k100();
  this_all_rates->rate_108 = nH_tot * k108();
  this_all_rates->rate_115 = nH_tot * k115();
  this_all_rates->rate_118 = nH_tot * k118();
  this_all_rates->rate_119 = nH_tot * k119();
  this_all_rates->rate_123 = nH_tot * k123();
  this_all_rates->rate_125 = nH_tot * k125();
  this_all_rates->rate_126 = nH_tot * k126();
  this_all_rates->rate_127 = nH_tot * k127();
  this_all_rates->rate_132 = nH_tot * k132();
  this_all_rates->rate_140 = nH_tot * k140();
  this_all_rates->rate_141 = nH_tot * k141();
  this_all_rates->rate_143 = nH_tot * k143();
  this_all_rates->rate_144 = nH_tot * k144();
  this_all_rates->rate_145 = nH_tot * k145();
  this_all_rates->rate_146 = nH_tot * k146();
  this_all_rates->rate_148 = nH_tot * k148();
  this_all_rates->rate_149 = nH_tot * k149();
  this_all_rates->rate_151 = nH_tot * k151();
  this_all_rates->rate_152 = nH_tot * k152();
  this_all_rates->rate_153 = nH_tot * k153();
  this_all_rates->rate_154 = nH_tot * k154();
  this_all_rates->rate_155 = nH_tot * k155();
  this_all_rates->rate_156 = nH_tot * k156();
  this_all_rates->rate_158 = nH_tot * k158();
  this_all_rates->rate_159 = nH_tot * k159();
  this_all_rates->rate_160 = nH_tot * k160();
  this_all_rates->rate_161 = nH_tot * k161();
  this_all_rates->rate_162 = nH_tot * k162();
  this_all_rates->rate_163 = nH_tot * k163();
  this_all_rates->rate_164 = nH_tot * k164();
  this_all_rates->rate_165 = nH_tot * k165();
  this_all_rates->rate_166 = nH_tot * k166();
  this_all_rates->rate_167 = nH_tot * k167();
  this_all_rates->rate_168 = nH_tot * k168();
  this_all_rates->rate_169 = nH_tot * k169();
  this_all_rates->rate_170 = nH_tot * k170();
  this_all_rates->rate_171 = nH_tot * k171();
  this_all_rates->rate_172 = nH_tot * k172();
  this_all_rates->rate_173 = nH_tot * k173();
  this_all_rates->rate_174 = nH_tot * k174();
  this_all_rates->rate_175 = nH_tot * k175();
  this_all_rates->rate_176 = nH_tot * k176();
  this_all_rates->rate_177 = nH_tot * k177();
  this_all_rates->rate_178 = nH_tot * k178();
  this_all_rates->rate_179 = nH_tot * k179();
  this_all_rates->rate_180 = nH_tot * k180();
  this_all_rates->rate_181 = nH_tot * k181();
  this_all_rates->rate_182 = nH_tot * k182();
  this_all_rates->rate_183 = nH_tot * k183();
  this_all_rates->rate_184 = nH_tot * k184();
  this_all_rates->rate_185 = nH_tot * k185();
  this_all_rates->rate_208 = nH_tot * k208();
  this_all_rates->rate_209 = nH_tot * k209();
  this_all_rates->rate_210 = nH_tot * k210();
  this_all_rates->rate_211 = nH_tot * k211();
  this_all_rates->rate_212 = nH_tot * k212();
  this_all_rates->rate_213 = nH_tot * k213();
  this_all_rates->rate_214 = nH_tot * k214();
  this_all_rates->rate_215 = nH_tot * k215();
  this_all_rates->rate_216 = nH_tot * k216();
  this_all_rates->rate_217 = nH_tot * k217();
  this_all_rates->rate_218 = nH_tot * k218();
  this_all_rates->rate_219 = nH_tot * k219();
  this_all_rates->rate_224 = nH_tot * k224();
  this_all_rates->rate_268 = R268(cr_rate);
  this_all_rates->rate_269 = R269(cr_rate);
  this_all_rates->rate_270 = R270(cr_rate);
  this_all_rates->rate_272 = R272(cr_rate);
  this_all_rates->rate_273 = R273(cr_rate);
  this_all_rates->rate_274 = R274(cr_rate);
  this_all_rates->rate_275 = R275(cr_rate);
  this_all_rates->rate_276 = R276(cr_rate);
  this_all_rates->rate_277 = R277(cr_rate);
  this_all_rates->rate_278 = R278(cr_rate);
  this_all_rates->rate_279 = R279(cr_rate);
  this_all_rates->rate_280 = R280(cr_rate);
  this_all_rates->rate_294 = nH_tot * k294();
  this_all_rates->rate_295 = nH_tot * k295();
  this_all_rates->rate_296 = nH_tot * k296();
  this_all_rates->rate_297 = nH_tot * k297();
  this_all_rates->rate_303 = nH_tot * k303();
  this_all_rates->rate_304 = nH_tot * k304();

  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      for(j = 0; j < (chimesRateTables.NonEqIon->N_Ions[i] - 1); j++)
        this_all_rates->BensRates[i].cosmicRays[j] = chimesRateTables.NonEqIon->NonEqRates[i].cosmicRays[j] * cr_rate;
    }
}

void update_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, double HI_column_density,
                  double H2_column_density, double HeI_column_density, double HeII_column_density, double CO_column_density,
                  double extinction, struct All_rate_variables_structure *this_all_rates)
{
  /* Some of the rate coefficients depend on
   * other variables, e.g. abundances etc.
   * These are updated here. */
  int i, j, k, l, T_index, N_species;
  double dT;
  double nHtotal = myGasVars->nH_tot;
  double xe      = myGasVars->abundances[myGlobalVars->speciesIndices[elec]];
  double xH      = myGasVars->abundances[myGlobalVars->speciesIndices[HI]];
  double xHII    = myGasVars->abundances[myGlobalVars->speciesIndices[HII]];
  double xH2     = myGasVars->abundances[myGlobalVars->speciesIndices[H2]];
  double xHe     = myGasVars->abundances[myGlobalVars->speciesIndices[HeI]];
  double dust_G;
  double cr_rate = myGasVars->cr_rate;
  double T       = myGasVars->temperature;
  int NHI_index, NHeI_index, NH_eff_index, NHe_eff_index;
  double dNHI, dNHeI, dNH_eff, dNHe_eff;
  double shieldFactor, S1, S2, S3;

  if(T > MAX_TEMPERATURE_FOR_RATES)
    T = MAX_TEMPERATURE_FOR_RATES;

  dust_G = 0.0;
  for(l = 0; l < myGlobalVars->N_spectra; l++)
    dust_G += myGasVars->isotropic_photon_density[l] * LIGHTSPEED * myGasVars->dust_G_parameter[l];

  /* T_index for additional rates. */
  get_index_1d_mydbl(chimesRateTables.RatesTables->Temperatures, chimesRateTables.RatesTables->N_Temperatures, log10(T), &T_index,
                     &dT);

  this_all_rates->rate_13  = nHtotal * k13(T_index, dT, HI_column_density);
  this_all_rates->rate_25  = nHtotal * (k25rr(T_index, dT, HeI_column_density) + k25di(T_index, dT));
  this_all_rates->rate_52  = R52(dust_G, extinction);
  this_all_rates->rate_111 = R111(dust_G, extinction);
  this_all_rates->rate_61  = nHtotal * k61(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_63  = nHtotal * k63(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_64  = nHtotal * k64(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_65  = nHtotal * k65(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_66  = nHtotal * k66(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_77  = nHtotal * k77(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_290 = nHtotal * k290(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_291 = nHtotal * k291(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_292 = nHtotal * k292(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_293 = nHtotal * k293(T, xe * nHtotal, dust_G, extinction, myGasVars->metallicity);
  this_all_rates->rate_8   = nHtotal * k8(T, nHtotal, xH, xH2, xHe);
  this_all_rates->rate_9   = nHtotal * k9(T, nHtotal, xH, xH2, xHe);
  this_all_rates->rate_60  = nHtotal * k60(T, myGlobalVars->grain_temperature, xH, myGasVars->metallicity);
  this_all_rates->rate_88  = nHtotal * k88(T, nHtotal, xH, xH2, xHe);
  this_all_rates->rate_281 = R281(cr_rate, T, xH2, myGasVars->abundances[myGlobalVars->speciesIndices[CO]]);
  this_all_rates->rate_53 =
      R53(H2_column_density, T, myGasVars->isotropic_photon_density, extinction, myGasVars->doppler_broad, myGlobalVars, myGasVars);
  this_all_rates->rate_239 = R239(dust_G, extinction);
  this_all_rates->rate_240 = R240(dust_G, extinction);
  this_all_rates->rate_241 = R241(dust_G, extinction);
  this_all_rates->rate_242 = R242(dust_G, extinction);
  this_all_rates->rate_243 = R243(dust_G, extinction);
  this_all_rates->rate_244 = R244(dust_G, extinction);
  this_all_rates->rate_245 = R245(dust_G, extinction);
  this_all_rates->rate_246 = R246(dust_G, extinction);
  this_all_rates->rate_247 = R247(dust_G, extinction);
  this_all_rates->rate_248 = R248(dust_G, extinction);
  this_all_rates->rate_249 = R249(dust_G, extinction);
  this_all_rates->rate_250 = R250(dust_G, extinction);
  this_all_rates->rate_251 = R251(dust_G, extinction);
  this_all_rates->rate_252 = R252(dust_G, extinction);
  this_all_rates->rate_253 = R253(dust_G, extinction);
  this_all_rates->rate_254 = R254(dust_G, extinction);
  this_all_rates->rate_255 = R255(dust_G, extinction);
  this_all_rates->rate_256 = R256(dust_G, extinction);
  this_all_rates->rate_257 = R257(dust_G, extinction);
  this_all_rates->rate_258 = R258(dust_G, extinction);
  this_all_rates->rate_259 = R259(dust_G, extinction);
  this_all_rates->rate_260 = R260(dust_G, extinction);
  this_all_rates->rate_261 = R261(dust_G, extinction);
  this_all_rates->rate_262 = R262(dust_G, extinction);
  this_all_rates->rate_263 = R263(dust_G, extinction);
  this_all_rates->rate_264 = R264(dust_G, extinction);
  this_all_rates->rate_265 = R265(dust_G, extinction);
  this_all_rates->rate_266 = R266(dust_G, extinction);
  this_all_rates->rate_267 = R267(dust_G, H2_column_density, CO_column_density, extinction, myGlobalVars);

  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HI_column_density, 1.0e-50)), &NHI_index, &dNHI);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HI_column_density + 3.0 * H2_column_density, 1.0e-50)), &NH_eff_index, &dNH_eff);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HeI_column_density, 1.0e-50)), &NHeI_index, &dNHeI);
  get_index_1d_mydbl(chimesRateTables.NonEqIon->shieldingColumnDensities, chimesRateTables.NonEqIon->shieldingColumnDimensions[0],
                     log10(max(HeI_column_density + 0.75 * HeII_column_density, 1.0e-50)), &NHe_eff_index, &dNHe_eff);
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      /* Secondary ionisations from Cosmic Rays; HI & HeI only */
      if(chimesRateTables.NonEqIon->IonIndexBegin[i] == HI || chimesRateTables.NonEqIon->IonIndexBegin[i] == HeI)
        this_all_rates->BensRates[i].cosmicRays[0] =
            chimesRateTables.NonEqIon->NonEqRates[i].cosmicRays[0] * cr_rate *
            (1.0 + cr_secondary_ionisation(xHII, chimesRateTables.NonEqIon->IonIndexBegin[i]));

      if(chimesRateTables.NonEqIon->IonIndexBegin[i] == HI)
        N_species = 3;
      else
        N_species = chimesRateTables.NonEqIon->N_Ions[i] - 1;

      for(j = 0; j < N_species; j++)
        {
          /* Zero the PhotoIon rates. */
          this_all_rates->BensRates[i].PhotoIon[j][0] = 0.0;
          if(chimesRateTables.NonEqIon->N_Auger[i] > 1)
            {
              for(k = 1; k < chimesRateTables.NonEqIon->N_Auger[i]; k++)
                this_all_rates->BensRates[i].PhotoIon[j][k] = 0.0;
            }

          /* Loop over UV spectra */
          for(l = 0; l < myGlobalVars->N_spectra; l++)
            {
              /* Standard photoionisation */
              if(chimesRateTables.NonEqIon->NonEqRates[i].E_thresh[j] >= 13.6)
                {
                  if(HI_column_density == 0.0 && H2_column_density == 0.0 && HeI_column_density == 0.0 && HeII_column_density == 0.0)
                    shieldFactor = 1.0;
                  else
                    {
                      if(chimesRateTables.NonEqIon->NonEqRates[i].E_thresh[j] < 15.4)
                        {
                          S1 = pow(10.0, interpol_1d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor1D[l][0][j],
                                                            NHI_index, dNHI));
                          S2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor2D[l][0][j],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor2D[l][1][j],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      else if(chimesRateTables.NonEqIon->NonEqRates[i].E_thresh[j] < 54.42)
                        {
                          S1 = 0.0;
                          S2 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor2D[l][0][j],
                                                            NH_eff_index, NHeI_index, dNH_eff, dNHeI));
                          S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor2D[l][1][j],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      else
                        {
                          S1 = 0.0;
                          S2 = 0.0;
                          S3 = pow(10.0, interpol_2d_fltdbl(chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor2D[l][1][j],
                                                            NH_eff_index, NHe_eff_index, dNH_eff, dNHe_eff));
                        }
                      shieldFactor = S1 + S2 + S3;
                    }
                  this_all_rates->BensRates[i].PhotoIon[j][0] += chimesRateTables.NonEqIon->NonEqRates[i].sigmaphot[l][j][0] *
                                                                 myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                }
              else
                {
                  /* For these species, UV attenuation is
                   * by dust. The gamma_d values are stored
                   * in the shieldFactor tables. */
                  shieldFactor =
                      exp(-pow(10.0, ((double)chimesRateTables.NonEqIon->NonEqRates[i].shieldFactor1D[l][0][j][0])) * extinction);
                  this_all_rates->BensRates[i].PhotoIon[j][0] += chimesRateTables.NonEqIon->NonEqRates[i].sigmaphot[l][j][0] *
                                                                 myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                }

              /* Auger ionisation */
              if(chimesRateTables.NonEqIon->N_Auger[i] > 1)
                {
                  for(k = 1; k < chimesRateTables.NonEqIon->N_Auger[i]; k++)
                    this_all_rates->BensRates[i].PhotoIon[j][k] += chimesRateTables.NonEqIon->NonEqRates[i].sigmaphot[l][j][k] *
                                                                   myGasVars->isotropic_photon_density[l] * LIGHTSPEED * shieldFactor;
                }
            }
        }
    }
}

void update_T_dependent_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars,
                              struct All_rate_variables_structure *this_all_rates)
{
  /* These are rates that depend only on T
   * - they only need to be updated in the
   * RhsFn if ThermEvol is On. */
  int T_index, i, j;
  double dT;
  double temperature = myGasVars->temperature;
  double nH_tot      = myGasVars->nH_tot;

  if(temperature > MAX_TEMPERATURE_FOR_RATES)
    temperature = MAX_TEMPERATURE_FOR_RATES;

  /* T_index for additional rates. */
  get_index_1d_mydbl(chimesRateTables.RatesTables->Temperatures, chimesRateTables.RatesTables->N_Temperatures, log10(temperature),
                     &T_index, &dT);

  this_all_rates->rate_1   = nH_tot * k1(T_index, dT);
  this_all_rates->rate_5   = nH_tot * k5(T_index, dT);
  this_all_rates->rate_11  = nH_tot * k11(T_index, dT);
  this_all_rates->rate_15  = nH_tot * k15(T_index, dT);
  this_all_rates->rate_16  = nH_tot * k16(T_index, dT);
  this_all_rates->rate_24  = nH_tot * k24(T_index, dT);
  this_all_rates->rate_26  = nH_tot * k26(T_index, dT);
  this_all_rates->rate_27  = nH_tot * k27(T_index, dT);
  this_all_rates->rate_2   = nH_tot * k2(T_index, dT);
  this_all_rates->rate_3   = nH_tot * k3(T_index, dT);
  this_all_rates->rate_6   = nH_tot * k6(T_index, dT);
  this_all_rates->rate_7   = nH_tot * k7(T_index, dT);
  this_all_rates->rate_10  = nH_tot * k10(T_index, dT);
  this_all_rates->rate_17  = nH_tot * k17(T_index, dT);
  this_all_rates->rate_83  = nH_tot * k83(T_index, dT);
  this_all_rates->rate_84  = nH_tot * k84(T_index, dT);
  this_all_rates->rate_85  = nH_tot * k85(temperature);
  this_all_rates->rate_86  = nH_tot * k86(temperature);
  this_all_rates->rate_87  = nH_tot * k87(T_index, dT);
  this_all_rates->rate_97  = nH_tot * k97(T_index, dT);
  this_all_rates->rate_101 = nH_tot * k101(T_index, dT);
  this_all_rates->rate_106 = nH_tot * k106(T_index, dT);
  this_all_rates->rate_107 = nH_tot * k107(T_index, dT);
  this_all_rates->rate_112 = nH_tot * k112(T_index, dT);
  this_all_rates->rate_113 = nH_tot * k113(T_index, dT);
  this_all_rates->rate_114 = nH_tot * k114(T_index, dT);
  this_all_rates->rate_116 = nH_tot * k116(T_index, dT);
  this_all_rates->rate_117 = nH_tot * k117(T_index, dT);
  this_all_rates->rate_120 = nH_tot * k120(T_index, dT);
  this_all_rates->rate_121 = nH_tot * k121(T_index, dT);
  this_all_rates->rate_122 = nH_tot * k122(T_index, dT);
  this_all_rates->rate_124 = nH_tot * k124(T_index, dT);
  this_all_rates->rate_128 = nH_tot * k128(T_index, dT);
  this_all_rates->rate_129 = nH_tot * k129(T_index, dT);
  this_all_rates->rate_130 = nH_tot * k130(T_index, dT);
  this_all_rates->rate_131 = nH_tot * k131(T_index, dT);
  this_all_rates->rate_133 = nH_tot * k133(T_index, dT);
  this_all_rates->rate_134 = nH_tot * k134(T_index, dT);
  this_all_rates->rate_135 = nH_tot * k135(T_index, dT);
  this_all_rates->rate_136 = nH_tot * k136(T_index, dT);
  this_all_rates->rate_137 = nH_tot * k137(T_index, dT);
  this_all_rates->rate_138 = nH_tot * k138(T_index, dT);
  this_all_rates->rate_139 = nH_tot * k139(T_index, dT);
  this_all_rates->rate_142 = nH_tot * k142(T_index, dT);
  this_all_rates->rate_147 = nH_tot * k147(T_index, dT);
  this_all_rates->rate_150 = nH_tot * k150(T_index, dT);
  this_all_rates->rate_157 = nH_tot * k157(T_index, dT);
  this_all_rates->rate_186 = nH_tot * k186(T_index, dT);
  this_all_rates->rate_187 = nH_tot * k187(T_index, dT);
  this_all_rates->rate_188 = nH_tot * k188(T_index, dT);
  this_all_rates->rate_189 = nH_tot * k189(T_index, dT);
  this_all_rates->rate_190 = nH_tot * k190(T_index, dT);
  this_all_rates->rate_191 = nH_tot * k191(T_index, dT);
  this_all_rates->rate_192 = nH_tot * k192(T_index, dT);
  this_all_rates->rate_193 = nH_tot * k193(T_index, dT);
  this_all_rates->rate_194 = nH_tot * k194(T_index, dT);
  this_all_rates->rate_195 = nH_tot * k195(T_index, dT);
  this_all_rates->rate_196 = nH_tot * k196(T_index, dT);
  this_all_rates->rate_197 = nH_tot * k197(T_index, dT);
  this_all_rates->rate_198 = nH_tot * k198(T_index, dT);
  this_all_rates->rate_199 = nH_tot * k199(T_index, dT);
  this_all_rates->rate_200 = nH_tot * k200(T_index, dT);
  this_all_rates->rate_201 = nH_tot * k201(T_index, dT);
  this_all_rates->rate_202 = nH_tot * k202(T_index, dT);
  this_all_rates->rate_203 = nH_tot * k203(T_index, dT);
  this_all_rates->rate_204 = nH_tot * k204(T_index, dT);
  this_all_rates->rate_205 = nH_tot * k205(T_index, dT);
  this_all_rates->rate_206 = nH_tot * k206(T_index, dT);
  this_all_rates->rate_207 = nH_tot * k207(T_index, dT);
  this_all_rates->rate_220 = nH_tot * k220(T_index, dT);
  this_all_rates->rate_221 = nH_tot * k221(T_index, dT);
  this_all_rates->rate_222 = nH_tot * k222(T_index, dT);
  this_all_rates->rate_223 = nH_tot * k223(T_index, dT);
  this_all_rates->rate_225 = nH_tot * k225(T_index, dT);
  this_all_rates->rate_226 = nH_tot * k226(T_index, dT);
  this_all_rates->rate_227 = nH_tot * k227(T_index, dT);
  this_all_rates->rate_228 = nH_tot * nH_tot * k228(T_index, dT);
  this_all_rates->rate_229 = nH_tot * nH_tot * k229(T_index, dT);
  this_all_rates->rate_230 = nH_tot * nH_tot * k230(T_index, dT);
  this_all_rates->rate_231 = nH_tot * nH_tot * k231(T_index, dT);
  this_all_rates->rate_232 = nH_tot * nH_tot * k232(T_index, dT);
  this_all_rates->rate_233 = nH_tot * nH_tot * k233(T_index, dT);
  this_all_rates->rate_234 = nH_tot * nH_tot * k234(T_index, dT);
  this_all_rates->rate_235 = nH_tot * nH_tot * k235(T_index, dT);
  this_all_rates->rate_236 = nH_tot * nH_tot * k236(T_index, dT);
  this_all_rates->rate_237 = nH_tot * nH_tot * k237(T_index, dT);
  this_all_rates->rate_238 = nH_tot * k238(T_index, dT);
  this_all_rates->rate_283 = nH_tot * k283(T_index, dT);

  /* T index for Bens tables. */
  get_index_1d_mydbl(chimesRateTables.NonEqIon->Temperatures, chimesRateTables.NonEqIon->N_Temperatures, log10(temperature), &T_index,
                     &dT);
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      for(j = 0; j < (chimesRateTables.NonEqIon->N_Ions[i] - 1); j++)
        {
          this_all_rates->BensRates[i].CollisIon[j] =
              pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].betacoll[j], T_index, dT)) * nH_tot;
          this_all_rates->BensRates[i].Recomb[j] =
              (pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].alpharad[j + 1], T_index, dT)) +
               pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].alphadi[j + 1], T_index, dT))) *
              nH_tot;

          if(chimesRateTables.NonEqIon->NonEqRates[i].CTHrec_mask[j + 1] == 1)
            this_all_rates->BensRates[i].CTHrec[j] =
                pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[j + 1], T_index, dT)) * nH_tot;
          if(chimesRateTables.NonEqIon->NonEqRates[i].CTHerec_mask[j + 1] == 1)
            this_all_rates->BensRates[i].CTHerec[j] =
                pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[j + 1], T_index, dT)) * nH_tot;
          if(chimesRateTables.NonEqIon->NonEqRates[i].CTHion_mask[j] == 1)
            this_all_rates->BensRates[i].CTHion[j] =
                pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[j], T_index, dT)) * nH_tot;
          if(chimesRateTables.NonEqIon->NonEqRates[i].CTHeion_mask[j] == 1)
            this_all_rates->BensRates[i].CTHeion[j] =
                pow(10.0, interpol_1d_mydbl(chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[j], T_index, dT)) * nH_tot;
        }
    }
}

void set_species_arrays(struct Species_Structure *mySpecies, struct gasVariables *myGasVars, int *total_network,
                        int *nonmolecular_network, struct globalVariables *myGlobalVars)
{
  /* This function determines which
   * species are to be included in the
   * network. Here we exclude species
   * that contain elements whose
   * abundances are zero. */
  int i;
  int inclSpeciesFlags[10];
  int inclSpeciesFlag_CO;
  double min_C_O_abundance;

  for(i = 0; i < 10; i++)
    {
      if(myGasVars->element_abundances[i] > 0.0)
        inclSpeciesFlags[i] = 1;
      else
        inclSpeciesFlags[i] = 0;
    }

  if(myGasVars->element_abundances[1] > 0.0 && myGasVars->element_abundances[3] > 0.0)
    inclSpeciesFlag_CO = 1;
  else
    inclSpeciesFlag_CO = 0;

  /* Determine which species are included */
  for(i = myGlobalVars->speciesIndices[elec]; i <= myGlobalVars->speciesIndices[Hm]; i++)
    {
      mySpecies[i].include_species   = 1;
      mySpecies[i].element_abundance = 1.0;
    }

  for(i = myGlobalVars->speciesIndices[HeI]; i <= myGlobalVars->speciesIndices[HeIII]; i++)
    {
      mySpecies[i].include_species   = inclSpeciesFlags[0];
      mySpecies[i].element_abundance = myGasVars->element_abundances[0];
    }

  if(myGlobalVars->element_included[0] == 1)
    {
      for(i = myGlobalVars->speciesIndices[CI]; i <= myGlobalVars->speciesIndices[Cm]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[1];
          mySpecies[i].element_abundance = myGasVars->element_abundances[1];
        }
    }

  if(myGlobalVars->element_included[1] == 1)
    {
      for(i = myGlobalVars->speciesIndices[NI]; i <= myGlobalVars->speciesIndices[NVIII]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[2];
          mySpecies[i].element_abundance = myGasVars->element_abundances[2];
        }
    }

  if(myGlobalVars->element_included[2] == 1)
    {
      for(i = myGlobalVars->speciesIndices[OI]; i <= myGlobalVars->speciesIndices[Om]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[3];
          mySpecies[i].element_abundance = myGasVars->element_abundances[3];
        }
    }

  if(myGlobalVars->element_included[3] == 1)
    {
      for(i = myGlobalVars->speciesIndices[NeI]; i <= myGlobalVars->speciesIndices[NeXI]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[4];
          mySpecies[i].element_abundance = myGasVars->element_abundances[4];
        }
    }

  if(myGlobalVars->element_included[4] == 1)
    {
      for(i = myGlobalVars->speciesIndices[MgI]; i <= myGlobalVars->speciesIndices[MgXIII]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[5];
          mySpecies[i].element_abundance = myGasVars->element_abundances[5];
        }
    }

  if(myGlobalVars->element_included[5] == 1)
    {
      for(i = myGlobalVars->speciesIndices[SiI]; i <= myGlobalVars->speciesIndices[SiXV]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[6];
          mySpecies[i].element_abundance = myGasVars->element_abundances[6];
        }
    }

  if(myGlobalVars->element_included[6] == 1)
    {
      for(i = myGlobalVars->speciesIndices[SI]; i <= myGlobalVars->speciesIndices[SXVII]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[7];
          mySpecies[i].element_abundance = myGasVars->element_abundances[7];
        }
    }

  if(myGlobalVars->element_included[7] == 1)
    {
      for(i = myGlobalVars->speciesIndices[CaI]; i <= myGlobalVars->speciesIndices[CaXXI]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[8];
          mySpecies[i].element_abundance = myGasVars->element_abundances[8];
        }
    }

  if(myGlobalVars->element_included[8] == 1)
    {
      for(i = myGlobalVars->speciesIndices[FeI]; i <= myGlobalVars->speciesIndices[FeXXVII]; i++)
        {
          mySpecies[i].include_species   = inclSpeciesFlags[9];
          mySpecies[i].element_abundance = myGasVars->element_abundances[9];
        }
    }

  mySpecies[myGlobalVars->speciesIndices[H2]].include_species  = 1;
  mySpecies[myGlobalVars->speciesIndices[H2p]].include_species = 1;
  mySpecies[myGlobalVars->speciesIndices[H3p]].include_species = 1;

  mySpecies[myGlobalVars->speciesIndices[H2]].element_abundance  = 1.0;
  mySpecies[myGlobalVars->speciesIndices[H2p]].element_abundance = 1.0;
  mySpecies[myGlobalVars->speciesIndices[H3p]].element_abundance = 1.0;

  if(myGlobalVars->element_included[2] == 1)
    {
      mySpecies[myGlobalVars->speciesIndices[OH]].include_species   = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[H2O]].include_species  = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[O2]].include_species   = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[OHp]].include_species  = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[H2Op]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[H3Op]].include_species = inclSpeciesFlags[3];
      mySpecies[myGlobalVars->speciesIndices[O2p]].include_species  = inclSpeciesFlags[3];

      mySpecies[myGlobalVars->speciesIndices[OH]].element_abundance   = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[H2O]].element_abundance  = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[O2]].element_abundance   = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[OHp]].element_abundance  = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[H2Op]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[H3Op]].element_abundance = myGasVars->element_abundances[3];
      mySpecies[myGlobalVars->speciesIndices[O2p]].element_abundance  = myGasVars->element_abundances[3];
    }

  if(myGlobalVars->element_included[0] == 1)
    {
      mySpecies[myGlobalVars->speciesIndices[C2]].include_species   = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[CH]].include_species   = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[CH2]].include_species  = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[CH3p]].include_species = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[CHp]].include_species  = inclSpeciesFlags[1];
      mySpecies[myGlobalVars->speciesIndices[CH2p]].include_species = inclSpeciesFlags[1];

      mySpecies[myGlobalVars->speciesIndices[C2]].element_abundance   = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[CH]].element_abundance   = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[CH2]].element_abundance  = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[CH3p]].element_abundance = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[CHp]].element_abundance  = myGasVars->element_abundances[1];
      mySpecies[myGlobalVars->speciesIndices[CH2p]].element_abundance = myGasVars->element_abundances[1];
    }

  if(myGlobalVars->element_included[0] == 1 && myGlobalVars->element_included[2] == 1)
    {
      min_C_O_abundance = min(myGasVars->element_abundances[1], myGasVars->element_abundances[3]);

      mySpecies[myGlobalVars->speciesIndices[HCOp]].include_species = inclSpeciesFlag_CO;
      mySpecies[myGlobalVars->speciesIndices[CO]].include_species   = inclSpeciesFlag_CO;
      mySpecies[myGlobalVars->speciesIndices[COp]].include_species  = inclSpeciesFlag_CO;
      mySpecies[myGlobalVars->speciesIndices[HOCp]].include_species = inclSpeciesFlag_CO;

      mySpecies[myGlobalVars->speciesIndices[HCOp]].element_abundance = min_C_O_abundance;
      mySpecies[myGlobalVars->speciesIndices[CO]].element_abundance   = min_C_O_abundance;
      mySpecies[myGlobalVars->speciesIndices[COp]].element_abundance  = min_C_O_abundance;
      mySpecies[myGlobalVars->speciesIndices[HOCp]].element_abundance = min_C_O_abundance;
    }

  /* Now loop through this array and determine the total
   * number of species that are included in the network. */
  *total_network = 0;
  for(i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    {
      if(mySpecies[i].include_species == 1)
        *total_network += 1;
    }

  /* Now subtract from this the number of
   * molecules included in the network. */
  *nonmolecular_network = *total_network;
  for(i = H2; i <= O2p; i++)
    {
      if(myGlobalVars->speciesIndices[i] > -1)
        {
          if(mySpecies[myGlobalVars->speciesIndices[i]].include_species == 1)
            *nonmolecular_network -= 1;
        }
    }
}
