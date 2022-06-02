/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        22/10/2018
 * \author      Ondrej Jaura
 * \brief       Hard-coded ideal photon sources
 * \details     
 */

#if (SX_SOURCES==11)  // simple source in the middle of the box

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5 } 
};       
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 5e48, 0.0} 
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==12) // baczynski ICs (old paper), photo-evaporization of a dense clump

double sxTS_num = 2;
double sxTS_positions[][3] = {
  {0.5, 0.0625, 0.5 }, 
  {0.0625, 0.5, 0.5 } 
};       
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 5e48, 0.0},
  {0.0, 0.0, 5e48, 0.0} 
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==13) // one ray source

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.05, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e48, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==14) // source in the center of the cosmological box

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {100.0, 100.0, 100.0}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e48, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==15) // Baczinsky's R/D-front test source

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e49, 0.0, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 }; // [erg]
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0, 0, 0, 0 };              // [cm^2]

#elif (SX_SOURCES==16) // Baczinsky's Photoevaporation of a dense clump by two sources

double sxTS_num = 2;
double sxTS_positions[][3] = {
  {0.5, 0.0625, 0.5},
  {0.0625, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1.61e48, 0.0, 0.0},
  {0.0, 0.0, 1.61e48, 0.0, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 0.72*ELECTRONVOLT_IN_ERGS, 
				   2.94*ELECTRONVOLT_IN_ERGS, 
				   4.44*ELECTRONVOLT_IN_ERGS, 0, 0, 0 };        // [erg]
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 5.38e-18, 2.43e-18, 6.01e-18, 0, 0, 0 }; // [cm^2]

#elif (SX_SOURCES==17) // Ten sources in the center, distributed on the sphere ~0.05 [cu]

double sxTS_num = 10;
double sxTS_positions[][3] = {
  {  0.521794494718 ,  0.5 ,  0.545  },
  {  0.473670664659 ,  0.524119827953 ,  0.535  },
  {  0.503785644927 ,  0.456864528605 ,  0.525  },
  {  0.529020684058 ,  0.537852343347 ,  0.515  },
  {  0.451011122648 ,  0.49133455738 ,  0.505  },
  {  0.541976295917 ,  0.473298116525 ,  0.495  },
  {  0.487617663835 ,  0.546061673342 ,  0.485  },
  {  0.480042140391 ,  0.461572355786 ,  0.475  },
  {  0.533540479046 ,  0.512248929153 ,  0.465  },
  {  0.479854355659 ,  0.508315829129 ,  0.455  }
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 1e48, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==18) // Hundred sources in the center, distributed on the sphere ~0.06 [cu]

double sxTS_num = 100;
double sxTS_positions[][3] = {
  {  0.508464041588 ,  0.5 ,  0.5594  },
  {  0.489244520067 ,  0.509852900659 ,  0.5582  },
  {  0.501637920428 ,  0.481336741531 ,  0.557  },
  {  0.513418249705 ,  0.517501730625 ,  0.5558  },
  {  0.475503779975 ,  0.495666963594 ,  0.5546  },
  {  0.523083146412 ,  0.485316391733 ,  0.5534  },
  {  0.492320099035 ,  0.52856884879 ,  0.5522  },
  {  0.485432129526 ,  0.471950451878 ,  0.551  },
  {  0.531435141871 ,  0.511480063394 ,  0.5498  },
  {  0.467476168983 ,  0.513425364649 ,  0.5486  },
  {  0.515591772777 ,  0.466681287215 ,  0.5474  },
  {  0.511457374387 ,  0.536527914971 ,  0.5462  },
  {  0.465662994612 ,  0.480101003519 ,  0.545  },
  {  0.540050376144 ,  0.491195037155 ,  0.5438  },
  {  0.475699573537 ,  0.534564856049 ,  0.5426  },
  {  0.494418965537 ,  0.456931542234 ,  0.5414  },
  {  0.534058740077 ,  0.528704742194 ,  0.5402  },
  {  0.454442884106 ,  0.501883929771 ,  0.539  },
  {  0.53302840838 ,  0.467132322262 ,  0.5378  },
  {  0.497803869977 ,  0.547493336511 ,  0.5366  },
  {  0.468961361436 ,  0.462805337532 ,  0.5354  },
  {  0.548858415346 ,  0.506573830682 ,  0.5342  },
  {  0.458866887841 ,  0.528619348072 ,  0.533  },
  {  0.511167181009 ,  0.4503607608 ,  0.5318  },
  {  0.525659735676 ,  0.544779660171 ,  0.5306  },
  {  0.450170980068 ,  0.48410318357 ,  0.5294  },
  {  0.548076604722 ,  0.477787389203 ,  0.5282  },
  {  0.479314185778 ,  0.549427695576 ,  0.527  },
  {  0.481666130304 ,  0.449027171729 ,  0.5258  },
  {  0.548442108221 ,  0.525459814435 ,  0.5246  },
  {  0.446575979203 ,  0.51408240043 ,  0.5234  },
  {  0.530147325963 ,  0.453113981431 ,  0.5222  },
  {  0.509519784724 ,  0.555392902964 ,  0.521  },
  {  0.455220007907 ,  0.465319857149 ,  0.5198  },
  {  0.556849425949 ,  0.495290141266 ,  0.5186  },
  {  0.461006178266 ,  0.54215117871 ,  0.5174  },
  {  0.500281825971 ,  0.442229068087 ,  0.5162  },
  {  0.539033000932 ,  0.543028186556 ,  0.515  },
  {  0.441857736016 ,  0.494611388042 ,  0.5138  },
  {  0.546727836209 ,  0.464535238289 ,  0.5126  },
  {  0.48945660594 ,  0.557955818014 ,  0.5114  },
  {  0.46850881879 ,  0.449957363119 ,  0.5102  },
  {  0.557211530987 ,  0.515679308727 ,  0.509  },
  {  0.447071814896 ,  0.527161870731 ,  0.5078  },
  {  0.520730440182 ,  0.444083197071 ,  0.5066  },
  {  0.522530751621 ,  0.555346230508 ,  0.5054  },
  {  0.44591364911 ,  0.474367468962 ,  0.5042  },
  {  0.557265049582 ,  0.482344573174 ,  0.503  },
  {  0.469676313568 ,  0.551741994948 ,  0.5018  },
  {  0.487403487133 ,  0.441340236417 ,  0.5006  },
  {  0.548912377524 ,  0.534745061878 ,  0.4994  },
  {  0.440487404703 ,  0.507416940137 ,  0.4982  },
  {  0.538841525933 ,  0.454367381588 ,  0.497  },
  {  0.502181229659 ,  0.559813060757 ,  0.4958  },
  {  0.458056090928 ,  0.45743771045 ,  0.4946  },
  {  0.559558093327 ,  0.503045245352 ,  0.4934  },
  {  0.454138505242 ,  0.537893050795 ,  0.4922  },
  {  0.508197004527 ,  0.441247901171 ,  0.491  },
  {  0.533531942415 ,  0.54869875602 ,  0.4898  },
  {  0.442593058324 ,  0.486790796866 ,  0.4886  },
  {  0.551039648336 ,  0.471082975641 ,  0.4874  },
  {  0.48198163652 ,  0.555541863288 ,  0.4862  },
  {  0.475891312975 ,  0.447143863082 ,  0.485  },
  {  0.55318333552 ,  0.522563085405 ,  0.4838  },
  {  0.445873055978 ,  0.519170652853 ,  0.4826  },
  {  0.52678472586 ,  0.449635146574 ,  0.4814  },
  {  0.514169343124 ,  0.554837849296 ,  0.4802  },
  {  0.452873326401 ,  0.469371963244 ,  0.479  },
  {  0.554981920839 ,  0.49082675734 ,  0.4778  },
  {  0.465958250876 ,  0.543515506622 ,  0.4766  },
  {  0.495747721443 ,  0.445440325083 ,  0.4754  },
  {  0.539584175242 ,  0.536979089638 ,  0.4742  },
  {  0.446420836933 ,  0.499477233324 ,  0.473  },
  {  0.539398212393 ,  0.464608745993 ,  0.4718  },
  {  0.494919219241 ,  0.552055985889 ,  0.4706  },
  {  0.468999287813 ,  0.45873747652 ,  0.4694  },
  {  0.550013233169 ,  0.509350749058 ,  0.4682  },
  {  0.457459133301 ,  0.526481590973 ,  0.467  },
  {  0.513262030013 ,  0.452518650399 ,  0.4658  },
  {  0.521907760821 ,  0.543207522676 ,  0.4646  },
  {  0.455502076406 ,  0.483255604047 ,  0.4634  },
  {  0.543241944287 ,  0.482642170231 ,  0.4622  },
  {  0.480272295102 ,  0.541107391786 ,  0.461  },
  {  0.487084681782 ,  0.457371904156 ,  0.4598  },
  {  0.537360639864 ,  0.522140970823 ,  0.4586  },
  {  0.458646847317 ,  0.508669300041 ,  0.4574  },
  {  0.523910402386 ,  0.466685548815 ,  0.4562  },
  {  0.504715881098 ,  0.539405081721 ,  0.455  },
  {  0.470969318601 ,  0.475044448763 ,  0.4538  },
  {  0.536768076713 ,  0.498838735682 ,  0.4526  },
  {  0.474818313829 ,  0.524574838384 ,  0.4514  },
  {  0.501871918809 ,  0.466586590716 ,  0.4502  },
  {  0.520013142775 ,  0.524463730629 ,  0.449  },
  {  0.470720974009 ,  0.49576932192 ,  0.4478  },
  {  0.522608113584 ,  0.484595026771 ,  0.4466  },
  {  0.494303544686 ,  0.524215499104 ,  0.4454  },
  {  0.489222594714 ,  0.480759222071 ,  0.4442  },
  {  0.517792288642 ,  0.50586808869 ,  0.443  },
  {  0.486699632867 ,  0.505988341516 ,  0.4418  },
  {  0.503343658354 ,  0.492224400421 ,  0.4406  }
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0},
  {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}, {0.0, 0.0, 1e47, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==19) // R/D-front  multi-frequency test source                                                                                                                                               
double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e49, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==20) // Star particle put into the center of the FiBY simulation

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {1.42, 1.42, 1.42}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 1e53, 0.0, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0, 0, 0, 0 };

#elif (SX_SOURCES==21) // Ten sources in the center, distributed on the sphere ~0.05 [cu]

double sxTS_num = 6;
double sxTS_positions[][3] = {
  {  0.5 ,  0.5 ,  0.5  },
  {  0.34,  0.67,  0.14 },
  {  0.45,  0.25,  0.67 },
  {  0.67,  0.35,  0.24 },
  {  0.87,  0.43,  0.53 },
  {  0.89,  0.56,  0.75 }
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {0.0, 0.0, 0e48, 0.0},
  {0.0, 0.0, 1e48, 0.0},
  {0.0, 0.0, 2e48, 0.0},
  {0.0, 0.0, 3e48, 0.0},
  {0.0, 0.0, 4e48, 0.0},
  {0.0, 0.0, 5e48, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0, 0 };

#elif (SX_SOURCES==22) // Test source for the FiBY chemistry

double sxTS_num = 4;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5},
  {0.3, 0.3, 0.5},
  {0.301, 0.301, 0.5},
  {0.7, 0.7, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  { 0e49, 1e49, 2e49, 3e49, 4e49, 5e49, 6e49, 7e49, 8e49, 9e49 },
  { 0e50, 1e50, 2e50, 3e50, 4e50, 5e50, 6e50, 7e50, 8e50, 9e50 },
  { 0e50, 1e50, 2e50, 3e50, 4e50, 5e50, 6e50, 7e50, 8e50, 9e50 },
  { 0e51, 1e51, 2e51, 3e51, 4e51, 5e51, 6e51, 7e51, 8e51, 9e51 }
};

#elif (SX_SOURCES==23) // Monochromatic source for FiBY chemistry

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  { 1e49, 1e49, 1e49, 1e49, 1e49, 1e49, 1e49, 1e49, 1e49, 1e49 }
};

#elif (SX_SOURCES==24) // Monochromatic source for FiBY chemistry

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  { 1e49, 0,0,0,0, 0,0,0,0,0 }
};

#elif (SX_SOURCES==25) // Monochromatic source for FiBY chemistry

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  { 0, 0,0,0,0, 1e49,0,0,0,0 }
};

#elif (SX_SOURCES==26) 
// Black body emission of a star with T = 1e5 and R = 1.5552*Rsol
// Total emission rate is 5e48 ph/s
// 10 bins for FiBY simulation

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {1.18109476e+48,   1.05822280e+48,   8.23215377e+47,   6.31144628e+47,
   6.70289477e+47,   3.42722142e+47,   1.63403804e+47,   7.47440109e+46,
   3.19118619e+46,   2.27303493e+46}
};

#elif (SX_SOURCES==27)
// the same as above but for SGChem bins
// bins 5.6-13.6 eV are switched off

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  //  {  1.08943252e+48, 5.29946691e+47, 3.55338124e+47, 1.88124292e+48, 2.76289817e+48 }
  {  0.0, 0.0, 3.55338124e+47, 1.88124292e+48, 2.76289817e+48 }  // ~5e48
};
double sxTS_energy[SX_NENERGY] = { 1.2804e-12, 9.77036e-12, 7.20688e-12, 2.19735e-11, 3.95975e-11, 3.7034e-11 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 5.33922e-18, 2.3221e-18, 6.48013e-18, 4.12921e-18, 4.3103e-19, 1.5725e-18 };

#elif (SX_SOURCES==28)
// testing source for Anna's halo

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {  0.0, 0.0, 5e+49, 0.0, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 5e-18, 0.0, 0.0, 0.0, 0.0, 0 };

#elif (SX_SOURCES==29)
// testing source for Bonnor-Ebert sphere cut

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.15, 0.15, 0.15}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {  0.0, 0.0, 1e+52, 0.0, 0.0}
};
double sxTS_energy[SX_NENERGY] = { 2.0*ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 };
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 5e-18, 0.0, 0.0, 0.0, 0.0, 0 };

#elif (SX_SOURCES==30)
// The Radiation Pressure test (Section 3.3, Rosdahl & Teyssier 2015)

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {  0.0, 0.0, 1.8e+50, 0.0, 0.0}   // phot/s
};
double sxTS_energy[SX_NENERGY] = { (15.0-13.6)*ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 };  // erg
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 3e-18, 0.0, 0.0, 0.0, 0.0, 0 };          // cm^2

#elif (SX_SOURCES==31)
// Starbench test (Bisbas et al. 2015)

double sxTS_num = 1;
double sxTS_positions[][3] = {
  {0.5, 0.5, 0.5}
};
double sxTS_ionizationRates[][SX_NFREQ] = {
  {  0.0, 0.0, 1e+49, 0.0, 0.0}   // phot/s
};
//double sxTS_energy[SX_NENERGY] = { (15.0-13.6)*ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 };  // erg
double sxTS_energy[SX_NENERGY] = { ELECTRONVOLT_IN_ERGS, 0, 0, 0, 0, 0 };              // erg
double sxTS_sigma[SX_NSIGMA] = { 0, 0, 6.3e-18, 0.0, 0.0, 0.0, 0.0, 0 };               // cm^2

#endif


