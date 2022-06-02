/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/parameters.c
 * \date        MM/YYYY
 * \author
 * \brief       Parses the parameter file
 * \details     This file contains the routine to parse the parameter file.
 *              Additionally the output list is also parsed.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

enum
{
  MAX_PARAMETERS = 300 /**< maximum number of parameters in the parameter file */
};

/*! \brief This function parses the parameter file.
 *
 *  Each parameter is defined by a keyword ("tag"), and can be either
 *  of type douple, int, or character string. Three arrays containing the name,
 *  type and address of the parameter are filled first. The routine then parses
 *  the parameter file and fills the referenced variables. The routine makes
 *  sure that each parameter appears exactly once in the parameter file,
 *  otherwise error messages are produced that complain about the missing
 *  parameters.
 *
 *  \param[in] fname The file name of the parameter file
 *
 *  \return void
 */
void read_parameter_file(const char *const fname)
{
  if(sizeof(int) != 4)
    mpi_terminate("Type `int' is not 32 bit on this platform. Stopping.");
  if(sizeof(long long) != 8)
    mpi_terminate("Type `long long' is not 64 bit on this platform. Stopping.");
  if(sizeof(float) != 4)
    mpi_terminate("Type `float' is not 32 bit on this platform. Stopping.");
  if(sizeof(double) != 8)
    mpi_terminate("Type `double' is not 64 bit on this platform. Stopping.");

  char buf[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200], buf1[MAXLEN_PARAM_TAG + 200], buf2[MAXLEN_PARAM_VALUE + 200],
      buf3[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 500];
  enum PARAM_TYPE id[MAX_PARAMETERS];
  void *addr[MAX_PARAMETERS];
  char tag[MAX_PARAMETERS][MAXLEN_PARAM_TAG];
  int param_handled[MAX_PARAMETERS];

  for(int i = 0; i < MAX_PARAMETERS; i++)
    param_handled[i] = 0;

  /* read parameter file on process 0 */
  int nt = 0;
  if(ThisTask == 0)
    {
      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LimitUBelowThisDensity");
      addr[nt] = &All.LimitUBelowThisDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LimitUBelowCertainDensityToThisValue");
      addr[nt] = &All.LimitUBelowCertainDensityToThisValue;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MultipleDomains");
      addr[nt] = &All.MultipleDomains;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "TopNodeFactor");
      addr[nt] = &All.TopNodeFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ActivePartFracForNewDomainDecomp");
      addr[nt] = &All.ActivePartFracForNewDomainDecomp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "GasSoftFactor");
      addr[nt] = &All.GasSoftFactor;
      id[nt++] = PARAM_REAL;

      for(int i = 0; i < NSOFTTYPES; i++)
        {
          sprintf(tag[nt], "SofteningComovingType%d", i);
          addr[nt] = &All.SofteningComoving[i];
          id[nt++] = PARAM_REAL;
        }

      for(int i = 0; i < NSOFTTYPES; i++)
        {
          sprintf(tag[nt], "SofteningMaxPhysType%d", i);
          addr[nt] = &All.SofteningMaxPhys[i];
          id[nt++] = PARAM_REAL;
        }

      for(int i = 0; i < NTYPES; i++)
        {
          sprintf(tag[nt], "SofteningTypeOfPartType%d", i);
          addr[nt] = &All.SofteningTypeOfPartType[i];
          id[nt++] = PARAM_INT;
        }

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinEgySpec");
      addr[nt] = &All.MinEgySpec;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinimumDensityOnStartUp");
      addr[nt] = &All.MinimumDensityOnStartUp;
      id[nt++] = PARAM_REAL;

#ifdef TOLERATE_WRITE_ERROR
      strcpy(tag[nt], "AlternativeOutputDir");
      addr[nt] = AlternativeOutputDir;
      id[nt++] = PARAM_STRING;
#endif

#ifdef REDUCE_FLUSH
      strcpy(tag[nt], "FlushCpuTimeDiff");
      addr[nt] = &All.FlushCpuTimeDiff;
      id[nt++] = PARAM_REAL;
#endif

#ifdef LEGACY_DISPLACEMENT_CONSTRAINT
      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SUBFIND
      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BECDM
      strcpy(tag[nt], "AxionMassEv");
      addr[nt] = &All.AxionMassEv;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DM_WINDTUNNEL
      strcpy(tag[nt], "DMWindtunnelInjectionRegion");
      addr[nt] = &All.DMWindtunnelInjectionRegion;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelSigmaVX");
      addr[nt] = &All.DMWindtunnelSigmaVX;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelSigmaVY");
      addr[nt] = &All.DMWindtunnelSigmaVY;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelSigmaVZ");
      addr[nt] = &All.DMWindtunnelSigmaVZ;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelVX");
      addr[nt] = &All.DMWindtunnelVX;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelVY");
      addr[nt] = &All.DMWindtunnelVY;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelVZ");
      addr[nt] = &All.DMWindtunnelVZ;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DMWindtunnelInjectionDensity");
      addr[nt] = &All.DMWindtunnelInjectionDensity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DM_WINDTUNNEL_STARS
      strcpy(tag[nt], "StarWindtunnelSigmaVX");
      addr[nt] = &All.StarWindtunnelSigmaVX;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "StarWindtunnelSigmaVY");
      addr[nt] = &All.StarWindtunnelSigmaVY;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "StarWindtunnelSigmaVZ");
      addr[nt] = &All.StarWindtunnelSigmaVZ;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "StarWindtunnelInjectionDensity");
      addr[nt] = &All.StarWindtunnelInjectionDensity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DM_WINDTUNNEL_EXTERNAL_SOURCE
      strcpy(tag[nt], "DMWindtunnelExternalSourceFile");
      addr[nt] = All.DMWindtunnelExternalSourceFile;
      id[nt++] = PARAM_STRING;
#endif

#ifdef WINDTUNNEL
      strcpy(tag[nt], "InjectionDensity");
      addr[nt] = &All.InjectionDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InjectionVelocity");
      addr[nt] = &All.InjectionVelocity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InjectionUtherm");
      addr[nt] = &All.InjectionUtherm;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InjectionRegion");
      addr[nt] = &All.InjectionRegion;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InjectionVolume");
      addr[nt] = &All.InjectionVolume;
      id[nt++] = PARAM_REAL;

#ifdef WINDTUNNEL_EXTERNAL_SOURCE
      strcpy(tag[nt], "WindTunnelExternalSourceFile");
      addr[nt] = All.WindTunnelExternalSourceFile;
      id[nt++] = PARAM_STRING;
#endif

#if defined(WINDTUNNEL_FIXVARIABLESININJECTIONREGION) && defined(MHD)
      strcpy(tag[nt], "InjectionBx_InGauss");
      addr[nt] = &All.InjectionBx_InGauss;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InjectionBy_InGauss");
      addr[nt] = &All.InjectionBy_InGauss;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InjectionBz_InGauss");
      addr[nt] = &All.InjectionBz_InGauss;
      id[nt++] = PARAM_REAL;
#endif

#ifdef WINDTUNNEL_READ_IN_BFIELD
      strcpy(tag[nt], "WindtunnelReadIn_InputFileName");
      addr[nt] = All.WindtunnelReadIn_InputFileName;
      id[nt++] = PARAM_STRING;
#endif

#endif /* WINDTUNNEL */

#ifdef CIRCUMSTELLAR
#ifdef CIRCUMSTELLAR_REFINEMENTS
      strcpy(tag[nt], "CircumstellarDerefinementDistance");
      addr[nt] = &All.CircumstellarDerefinementDistance;
      id[nt++] = PARAM_REAL;
#endif

#ifdef CIRCUMSTELLAR_IRRADIATION
      strcpy(tag[nt], "IrradiationTempScaling");
      addr[nt] = &All.IrradiationTempScaling;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef LOCALLY_ISOTHERM_DISK
      strcpy(tag[nt], "AspectRatio");
      addr[nt] = &All.AspectRatio;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CentralMass");
      addr[nt] = &All.CentralMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InnerRadius");
      addr[nt] = &All.InnerRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OuterRadius");
      addr[nt] = &All.OuterRadius;
      id[nt++] = PARAM_REAL;
#endif

#ifdef CIRCUMSTELLAR_WBOUNDARIES
      strcpy(tag[nt], "InnerRadius");
      addr[nt] = &All.inner_radius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OuterRadius");
      addr[nt] = &All.outer_radius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CircumstellarBoundaryDensity");
      addr[nt] = &All.CircumstellarBoundaryDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "EvanescentBoundaryStrength");
      addr[nt] = &All.EvanescentBoundaryStrength;
      id[nt++] = PARAM_REAL;
#endif

#ifdef CENTRAL_MASS_POTENTIAL
      strcpy(tag[nt], "CentralMass");
      addr[nt] = &All.CentralMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SofteningCentral");
      addr[nt] = &All.SofteningCentral;
      id[nt++] = PARAM_REAL;
#endif

#if defined(ACCRETE_ONTO_CENTRAL_POTENTIAL) && defined(CENTRAL_MASS_POTENTIAL)
      strcpy(tag[nt], "CentralAccretionRadius");
      addr[nt] = &All.CentralAccretionRadius;
      id[nt++] = PARAM_REAL;
#endif

#ifdef PERTURB_VELOCITIES
      strcpy(tag[nt], "VelocityPerturbation");
      addr[nt] = &All.VelocityPerturbation;
      id[nt++] = PARAM_REAL;
#endif

#ifdef STAR_PLANET_POTENTIAL
      strcpy(tag[nt], "MassRatio");
      addr[nt] = &All.MassRatio;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SofteningPlanet");
      addr[nt] = &All.SofteningPlanet;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PlanetGrowthTime");
      addr[nt] = &All.PlanetGrowthTime;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BINARY_POTENTIAL
      strcpy(tag[nt], "BinaryMassRatio");
      addr[nt] = &All.BinaryMassRatio;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BinarySoftening");
      addr[nt] = &All.BinarySoftening;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BinaryGrowthTime");
      addr[nt] = &All.BinaryGrowthTime;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BinaryEccentricity");
      addr[nt] = &All.BinaryEccentricity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BinaryBarycentricCoord");
      addr[nt] = &All.BinaryBarycentricCoord;
      id[nt++] = PARAM_INT;
#endif

#if defined(ISOTHERM_EQS) || defined(VS_TURB) || defined(AB_TURB)
      strcpy(tag[nt], "IsoSoundSpeed");
      addr[nt] = &All.IsoSoundSpeed;
      id[nt++] = PARAM_REAL;
#endif

#if(defined(VS_TURB) || defined(AB_TURB)) && defined(POWERSPEC_GRID)
      strcpy(tag[nt], "TimeBetTurbSpectrum");
      addr[nt] = &All.TimeBetTurbSpectrum;
      id[nt++] = PARAM_REAL;
#endif

#if defined(AB_TURB) && defined(AB_TURB_DECAYING)
      strcpy(tag[nt], "TimeTurbDecay");
      addr[nt] = &All.TimeTurbDecay;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DVR_RENDER
#if DVR_RENDER == 1
      strcpy(tag[nt], "DvrRenderTimeIntverallInGyr");
      addr[nt] = &All.DvrRenderTimeIntverallInGyr;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DvrPixelX");
      addr[nt] = &All.DvrPixelX;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "DvrPixelY");
      addr[nt] = &All.DvrPixelY;
      id[nt++] = PARAM_INT;
#endif
      strcpy(tag[nt], "DvrTauScaleFactor");
      addr[nt] = &All.DvrTauScaleFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DvrTauFloor");
      addr[nt] = &All.DvrTauFloor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DvrOutputDir");
      addr[nt] = &All.DvrOutputDir;
      id[nt++] = PARAM_STRING;
#endif

#ifdef ADAPTIVE_HYDRO_SOFTENING
      strcpy(tag[nt], "MinimumComovingHydroSoftening");
      addr[nt] = &All.MinimumComovingHydroSoftening;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AdaptiveHydroSofteningSpacing");
      addr[nt] = &All.AdaptiveHydroSofteningSpacing;
      id[nt++] = PARAM_REAL;
#endif

#ifdef NODEREFINE_BACKGROUND_GRID
      strcpy(tag[nt], "MeanVolume");
      addr[nt] = &All.MeanVolume;
      id[nt++] = PARAM_REAL;
#endif

#ifdef REFINEMENT_CGM
      strcpy(tag[nt], "TargetVolumeRelativeToSFRThreshold");
      addr[nt] = &All.TargetVolumeRelativeToSFRThreshold;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinMassForCGMRefinement");
      addr[nt] = &All.MinMassForCGMRefinement;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "FracRadiusForCGMRefinement");
      addr[nt] = &All.FracRadiusForCGMRefinement;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_NEW_CENTERING
      strcpy(tag[nt], "BlackHoleCenteringMassMultiplier");
      addr[nt] = &All.BlackHoleCenteringMassMultiplier;
      id[nt++] = PARAM_REAL;
#endif

#if defined(BLACK_HOLES) && defined(BH_FRICTION)
      strcpy(tag[nt], "BHFrictionCoefficient");
      addr[nt] = &All.BHFrictionCoefficient;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BHFrictionAvgTime");
      addr[nt] = &All.BHFrictionAvgTime;
      id[nt++] = PARAM_REAL;
#endif

#if defined(BLACK_HOLES) && defined(BH_SPIN_EVOLUTION)
      strcpy(tag[nt], "BHInitialSpin");
      addr[nt] = &All.BHInitialSpin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ShakuraSunyaevParameter");
      addr[nt] = &All.ShakuraSunyaevParameter;
      id[nt++] = PARAM_REAL;
#endif

#ifdef REFINEMENT_AROUND_BH
#ifdef REFINEMENT_AROUND_BH_FIXED
      strcpy(tag[nt], "RefBHRadius");
      addr[nt] = &All.RefBHRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RefBHMinCellRadius");
      addr[nt] = &All.RefBHMinCellRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RefBHMaxCellRadius");
      addr[nt] = &All.RefBHMaxCellRadius;
      id[nt++] = PARAM_REAL;
#else
      strcpy(tag[nt], "RefBHRadiusHSML");
      addr[nt] = &All.RefBHRadiusHSML;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RefBHMinCellRadiusRBondi");
      addr[nt] = &All.RefBHMinCellRadiusRBondi;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RefBHMaxCellRadiusHSML");
      addr[nt] = &All.RefBHMaxCellRadiusHSML;
      id[nt++] = PARAM_REAL;
#endif
      strcpy(tag[nt], "RefBHMinCellMass");
      addr[nt] = &All.RefBHMinCellMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RefBHLowerFactorC");
      addr[nt] = &All.RefBHLowerFactorC;
      id[nt++] = PARAM_REAL;
#endif

#ifdef REFINEMENT_AROUND_DM
      strcpy(tag[nt], "RefinementCellsPerSoftening");
      addr[nt] = &All.RefinementCellsPerSoftening;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
      strcpy(tag[nt], "BHBipolarTheta");
      addr[nt] = &All.BHBipolarTheta;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BHBipolarEfficiency");
      addr[nt] = &All.BHBipolarEfficiency;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BHBipolarColdFraction");
      addr[nt] = &All.BHBipolarColdFraction;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BHBipolarColdTemp");
      addr[nt] = &All.BHBipolarColdTemp;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MIN_METALLICITY_ON_STARTUP
      strcpy(tag[nt], "MinimumMetallicityOnStartUp");
      addr[nt] = &All.MinimumMetallicityOnStartUp;
      id[nt++] = PARAM_REAL;
#endif

#ifndef VORONOI_STATIC_MESH
#ifdef REGULARIZE_MESH_FACE_ANGLE
      strcpy(tag[nt], "CellMaxAngleFactor");
      addr[nt] = &All.CellMaxAngleFactor;
      id[nt++] = PARAM_REAL;
#else
      strcpy(tag[nt], "CellShapingFactor");
      addr[nt] = &All.CellShapingFactor;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "CellShapingSpeed");
      addr[nt] = &All.CellShapingSpeed;
      id[nt++] = PARAM_REAL;
#endif

#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES)
      strcpy(tag[nt], "PicXpixels");
      addr[nt] = &All.PicXpixels;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "PicYpixels");
      addr[nt] = &All.PicYpixels;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "PicXaxis");
      addr[nt] = &All.PicXaxis;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "PicYaxis");
      addr[nt] = &All.PicYaxis;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "PicZaxis");
      addr[nt] = &All.PicZaxis;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "PicXmin");
      addr[nt] = &All.PicXmin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PicXmax");
      addr[nt] = &All.PicXmax;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PicYmin");
      addr[nt] = &All.PicYmin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PicYmax");
      addr[nt] = &All.PicYmax;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PicZmin");
      addr[nt] = &All.PicZmin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PicZmax");
      addr[nt] = &All.PicZmax;
      id[nt++] = PARAM_REAL;
#endif

#ifdef VORONOI_FREQUENT_IMAGES
      strcpy(tag[nt], "TimeBetweenImages");
      addr[nt] = &All.TimeBetweenImages;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SUBBOX_SNAPSHOTS
      strcpy(tag[nt], "SubboxCoordinatesPath");
      addr[nt] = &All.SubboxCoordinatesPath;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "SubboxMinTime");
      addr[nt] = &All.SubboxMinTime;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SubboxMaxTime");
      addr[nt] = &All.SubboxMaxTime;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SubboxSyncModulo");
      addr[nt] = &All.SubboxSyncModulo;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SubboxNumFilesPerSnapshot");
      addr[nt] = &All.SubboxNumFilesPerSnapshot;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SubboxNumFilesWrittenInParallel");
      addr[nt] = &All.SubboxNumFilesWrittenInParallel;
      id[nt++] = PARAM_INT;
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_SINKS)
      strcpy(tag[nt], "CircumstellarSinkRadius");
      addr[nt] = &All.CircumstellarSinkRadius;
      id[nt++] = PARAM_REAL;
#endif

#if defined(FOF) && (defined(BLACK_HOLES) || defined(GFM_WINDS_VARIABLE) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_WINDS_LOCAL))
      strcpy(tag[nt], "TimeBetOnTheFlyFoF");
      addr[nt] = &All.TimeBetOnTheFlyFoF;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BLACK_HOLES

#ifdef BH_RELATIVE_NGB_DEVIATION
      strcpy(tag[nt], "DesNumNgbBlackHoleRelDeviationFactor");
      addr[nt] = &All.DesNumNgbBlackHoleRelDeviationFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_THERMALFEEDBACK_ACC
      strcpy(tag[nt], "BlackholeDeltaTemp");
      addr[nt] = &All.BlackholeDeltaTemp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BlackholeDeltaTime");
      addr[nt] = &All.BlackholeDeltaTime;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_BONDI_DENSITY
      strcpy(tag[nt], "BlackHoleAccretionSlope");
      addr[nt] = &All.BlackHoleAccretionSlope;
      id[nt++] = PARAM_REAL;
#else
      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt]           = &All.BlackHoleAccretionFactor;
      id[nt++]           = PARAM_REAL;
#endif

#ifdef BH_BONDI_DISK_VORTICITY
      strcpy(tag[nt], "DiskVorticityRadius");
      addr[nt] = &All.DiskVorticityRadius;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = PARAM_REAL;

#ifdef FOF
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MASSIVE_SEEDS
      strcpy(tag[nt], "DesNumNgbSeed");
      addr[nt] = &All.DesNumNgbSeed;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SeedMaxAccretionRadius");
      addr[nt] = &All.SeedMaxAccretionRadius;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "DesNumNgbBlackHole");
      addr[nt] = &All.DesNumNgbBlackHole;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "BlackHoleMaxAccretionRadius");
      addr[nt] = &All.BlackHoleMaxAccretionRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BlackHoleRadiativeEfficiency");
      addr[nt] = &All.BlackHoleRadiativeEfficiency;
      id[nt++] = PARAM_REAL;

#ifdef BH_NF_RADIO
      strcpy(tag[nt], "RadioModeMachnumber");
      addr[nt] = &All.RadioModeMachnumber;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioRelativeBubbleSize");
      addr[nt] = &All.RadioRelativeBubbleSize;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioRelativeBubbleEnergy");
      addr[nt] = &All.RadioRelativeBubbleEnergy;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioRelativeMaxDist");
      addr[nt] = &All.RadioRelativeMaxDist;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioModeMetallicityInSolar");
      addr[nt] = &All.RadioModeMetallicityInSolar;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_BASED_CGM_ZOOM
      strcpy(tag[nt], "CGM_MinRadius");
      addr[nt] = &All.CGM_MinRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CGM_MaxRadius");
      addr[nt] = &All.CGM_MaxRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "IGM_Radius");
      addr[nt] = &All.IGM_Radius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CGM_RefinementFactor");
      addr[nt] = &All.CGM_RefinementFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CGM_RadiiRedshift");
      addr[nt] = &All.CGM_RadiiRedshift;
      id[nt++] = PARAM_REAL;
#endif

#endif /* BLACK_HOLES */

#ifdef MRT

#ifdef MRT_SUBCYCLE
      strcpy(tag[nt], "RTNumSubCycles");
      addr[nt] = &All.RTNumSubCycles;
      id[nt++] = PARAM_INT;
#else
      All.RTNumSubCycles = 1.0;
#endif

#ifdef MRT_CHEM_SG
      strcpy(tag[nt], "UVKappa");
      addr[nt] = &All.UVKappa;
      id[nt++] = PARAM_REAL;
#ifdef MRT_IR
      strcpy(tag[nt], "IRKappa");
      addr[nt] = &All.IRKappa;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef MRT_RIEMANN_HLLE
      strcpy(tag[nt], "HLLEFile");
      addr[nt] = &All.HLLEFile;
      id[nt++] = PARAM_STRING;
#endif

#ifdef MRT_SINGLE_STAR
      strcpy(tag[nt], "StellarParamFile");
      addr[nt] = &All.StellarParamFile;
      id[nt++] = PARAM_STRING;
#endif

#ifdef MRT_IR_GRAIN_KAPPA
      strcpy(tag[nt], "GrainKappaPath");
      addr[nt] = &All.GrainKappaPath;
      id[nt++] = PARAM_STRING;
#endif
#ifdef MRT_BH
#ifdef MRT_BH_PULSED
      strcpy(tag[nt], "AGNPulseTime");
      addr[nt] = &All.AGNPulseTime;
      id[nt++] = PARAM_REAL;
#endif
#ifdef MRT_BH_BIPOLAR
      strcpy(tag[nt], "MRTBH_OpeningAngle");
      addr[nt] = &All.MRTBH_OpeningAngle;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "LogAGNLuminosity");
      addr[nt] = &All.LogAGNLuminosity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PhotonDesNumNgb");
      addr[nt] = &All.PhotonDesNumNgb;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "PhotonMaxNumNgbDeviation");
      addr[nt] = &All.PhotonMaxNumNgbDeviation;
      id[nt++] = PARAM_INT;
#endif

#if(defined(MRT_SOURCES) && defined(MRT_STARS)) || (defined(MRT_LOCAL_FEEDBACK))
      strcpy(tag[nt], "SpectrumTablePath");
      addr[nt] = &All.SpectrumTablePath;
      id[nt++] = PARAM_STRING;
#endif

#if defined(MRT_SOURCES) && defined(MRT_STARS)
      strcpy(tag[nt], "EscapeFraction");
      addr[nt] = &All.EscapeFraction;
      id[nt++] = PARAM_REAL;
#endif

#endif

#if defined(COOLING) && !defined(SIMPLE_COOLING) && !defined(GRACKLE) && !defined(CHIMES)
      strcpy(tag[nt], "TreecoolFile");
      addr[nt] = &All.TreecoolFile;
      id[nt++] = PARAM_STRING;
#ifdef UVB_START
      strcpy(tag[nt], "UVBStartRedshift");
      addr[nt] = &All.UVBStartRedshift;
      id[nt++] = PARAM_REAL;
#endif
#ifdef UVB_SELF_SHIELDING
      strcpy(tag[nt], "SelfShieldingFile");
      addr[nt] = &All.SelfShieldingFile;
      id[nt++] = PARAM_STRING;
#endif
#ifdef GFM_AGN_RADIATION
      strcpy(tag[nt], "TreecoolFileAGN");
      addr[nt] = &All.TreecoolFileAGN;
      id[nt++] = PARAM_STRING;
#endif
#ifdef RADCOOL
      strcpy(tag[nt], "TreecoolFileRAD");
      addr[nt] = &All.TreecoolFileRAD;
      id[nt++] = PARAM_STRING;
#endif
#endif

#ifdef GRACKLE
      strcpy(tag[nt], "GrackleOn");
      addr[nt] = &All.GrackleOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "GrackleRadiativeCooling");
      addr[nt] = &All.GrackleRadiativeCooling;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "GrackleMetalCooling");
      addr[nt] = &All.GrackleMetalCooling;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "GrackleUVB");
      addr[nt] = &All.GrackleUVB;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "GrackleSelfShieldingMethod");
      addr[nt] = &All.GrackleSelfShieldingMethod;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "GrackleDataFile");
      addr[nt] = &All.GrackleDataFile;
      id[nt++] = PARAM_STRING;
#ifndef METALS
      strcpy(tag[nt], "GrackleInitialMetallicity");
      addr[nt] = &All.GrackleInitialMetallicity;
      id[nt++] = PARAM_REAL;
#endif
#ifdef GRACKLE_PHOTOELECTRIC
      strcpy(tag[nt], "GracklePhotoelectricHeatingRate");
      addr[nt] = &All.GracklePhotoelectricHeatingRate;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef CHIMES
      strcpy(tag[nt], "Chimes_data_path");
      addr[nt] = All.ChimesDataPath;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "PhotoIon_table_path");
      addr[nt] = ChimesGlobalVars.PhotoIonTablePath;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "EqAbundance_table_path");
      addr[nt] = ChimesGlobalVars.EqAbundanceTablePath;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "Thermal_Evolution_On");
      addr[nt] = &All.ChimesThermEvolOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "Chemistry_eqm");
      addr[nt] = &All.ChimesForceEqOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "Reduction_On");
      addr[nt] = &ChimesGlobalVars.reductionOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "StaticMolCooling");
      addr[nt] = &ChimesGlobalVars.StaticMolCooling;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "CellSelfShieldingOn");
      addr[nt] = &ChimesGlobalVars.cellSelfShieldingOn;
      id[nt++] = PARAM_INT;

#if defined(CHIMES_JEANS_SHIELDING) || defined(CHIMES_SOBOLEV_SHIELDING)
      strcpy(tag[nt], "ShieldingLengthFactor");
      addr[nt] = &All.ChimesShieldingLengthFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxShieldingLength_kpc");
      addr[nt] = &All.ChimesMaxShieldingLength_kpc;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "Reduction_N_Ions_Low");
      addr[nt] = &ChimesGlobalVars.n_ions_low;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "Reduction_N_Ions_Med");
      addr[nt] = &ChimesGlobalVars.n_ions_med;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "Reduction_N_Ions_High");
      addr[nt] = &ChimesGlobalVars.n_ions_high;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "Grain_Temperature");
      addr[nt] = &ChimesGlobalVars.grain_temperature;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CrRate");
      addr[nt] = &All.ChimesCrRate;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "macrostep_tolerance");
      addr[nt] = &ChimesGlobalVars.time_tolerance;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "min_macrostep");
      addr[nt] = &ChimesGlobalVars.min_subcyclestep;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "max_mol_temperature");
      addr[nt] = &ChimesGlobalVars.T_mol;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "IsotropicPhotonDensity");
      addr[nt] = &All.ChimesIsotropicPhotonDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TemperatureEqm_Thresh");
      addr[nt] = &ChimesGlobalVars.T_EqThresh;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "relativeTolerance");
      addr[nt] = &ChimesGlobalVars.relativeTolerance;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "absoluteTolerance");
      addr[nt] = &ChimesGlobalVars.absoluteTolerance;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "thermalAbsoluteTolerance");
      addr[nt] = &ChimesGlobalVars.thermalAbsoluteTolerance;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "scale_metal_tolerances");
      addr[nt] = &ChimesGlobalVars.scale_metal_tolerances;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeCarbon");
      addr[nt] = &ChimesGlobalVars.element_included[0];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeNitrogen");
      addr[nt] = &ChimesGlobalVars.element_included[1];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeOxygen");
      addr[nt] = &ChimesGlobalVars.element_included[2];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeNeon");
      addr[nt] = &ChimesGlobalVars.element_included[3];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeMagnesium");
      addr[nt] = &ChimesGlobalVars.element_included[4];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeSilicon");
      addr[nt] = &ChimesGlobalVars.element_included[5];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeSulphur");
      addr[nt] = &ChimesGlobalVars.element_included[6];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeCalcium");
      addr[nt] = &ChimesGlobalVars.element_included[7];
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "IncludeIron");
      addr[nt] = &ChimesGlobalVars.element_included[8];
      id[nt++] = PARAM_INT;
#endif

#if defined(RADCOOL) && defined(COOLING)
      strcpy(tag[nt], "NewStarsOn");
      addr[nt] = &All.NewStarsOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "OldStarsOn");
      addr[nt] = &All.OldStarsOn;
      id[nt++] = PARAM_INT;

#ifdef RADCOOL_HOTHALO
      strcpy(tag[nt], "HotHaloOn");
      addr[nt] = &All.HotHaloOn;
      id[nt++] = PARAM_INT;
#endif

      strcpy(tag[nt], "SelfShieldingOn");
      addr[nt] = &All.SelfShieldingOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SelfShieldingDensity");
      addr[nt] = &All.SelfShieldingDensity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_AGN_RADIATION
      strcpy(tag[nt], "SelfShieldingDensity");
      addr[nt] = &All.SelfShieldingDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ObscurationFactor");
      addr[nt] = &All.ObscurationFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ObscurationSlope");
      addr[nt] = &All.ObscurationSlope;
      id[nt++] = PARAM_REAL;
#endif

#if defined(BH_ADIOS_WIND) || defined(BH_BUBBLES)
      strcpy(tag[nt], "RadioFeedbackFactor");
      addr[nt] = &All.RadioFeedbackFactor;
      id[nt++] = PARAM_REAL;
#endif

#if defined(BH_ADIOS_WIND)
#ifdef BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY
      strcpy(tag[nt], "RadioFeedbackMinDensityFactor");
      addr[nt] = &All.RadioFeedbackMinDensityFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_ADIOS_DENS_DEP_EFFICIANCY
      strcpy(tag[nt], "RadioFeedbackFactorPivotMass");
      addr[nt] = &All.RadioFeedbackFactorPivotMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioFeedbackFactorSlope");
      addr[nt] = &All.RadioFeedbackFactorSlope;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioFeedbackFactorRefDensityFactor");
      addr[nt] = &All.RadioFeedbackFactorRefDensityFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadioFeedbackFactorMaxEfficiency");
      addr[nt] = &All.RadioFeedbackFactorMaxEfficiency;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_ADIOS_RANDOMIZED
      strcpy(tag[nt], "RadioFeedbackReiorientationFactor");
      addr[nt] = &All.RadioFeedbackReiorientationFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_ADIOS_WIND_WITH_QUASARTHRESHOLD
      strcpy(tag[nt], "QuasarThreshold");
      addr[nt] = &All.QuasarThreshold;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef BH_BUBBLES
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BlackHoleRadioTriggeringFactor");
      addr[nt] = &All.BlackHoleRadioTriggeringFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DefaultICMDensity");
      addr[nt] = &All.DefaultICMDensity;
      id[nt++] = PARAM_REAL;

#ifdef UNIFIED_FEEDBACK
      strcpy(tag[nt], "RadioThreshold");
      addr[nt] = &All.RadioThreshold;
      id[nt++] = PARAM_REAL;
#endif

#ifdef BH_MAGNETIC_BUBBLES
      strcpy(tag[nt], "MagneticEnergyFraction");
      addr[nt] = &All.MagneticEnergyFraction;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef GFM_UVB_CORRECTIONS
      strcpy(tag[nt], "UV_HeII_threshold");
      addr[nt] = &All.UV_HeII_threshold;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UV_HeII_alpha");
      addr[nt] = &All.UV_HeII_alpha;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UV_HeII_beta");
      addr[nt] = &All.UV_HeII_alpha;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UV_HeatBoost");
      addr[nt] = &All.UV_HeatBoost;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WIND_ENERGY_METAL_DEPENDENCE
      strcpy(tag[nt], "WindEnergyReductionFactor");
      addr[nt] = &All.WindEnergyReductionFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "WindEnergyReductionMetallicity");
      addr[nt] = &All.WindEnergyReductionMetallicity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "WindEnergyReductionExponent");
      addr[nt] = &All.WindEnergyReductionExponent;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_STELLAR_EVOLUTION
#ifndef GFM_SINGLE_CELL_INJECTION
      strcpy(tag[nt], "DesNumNgbEnrichment");
      addr[nt] = &All.DesNumNgbEnrichment;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MaxNumNgbDeviationEnrichment");
      addr[nt] = &All.MaxNumNgbDeviationEnrichment;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "SNII_MinMass_Msun");
      addr[nt] = &All.SNII_MinMass_Msun;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SNII_MaxMass_Msun");
      addr[nt] = &All.SNII_MaxMass_Msun;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "IMF_MinMass_Msun");
      addr[nt] = &All.IMF_MinMass_Msun;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "IMF_MaxMass_Msun");
      addr[nt] = &All.IMF_MaxMass_Msun;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGB_MassTransferOn");
      addr[nt] = &All.AGB_MassTransferOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SNIa_Rate_Norm");
      addr[nt] = &All.SNIa_Rate_Norm;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SNIa_Rate_TAU");
      addr[nt] = &All.SNIa_Rate_TAU;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SNIa_MassTransferOn");
      addr[nt] = &All.SNIa_MassTransferOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SNII_MassTransferOn");
      addr[nt] = &All.SNII_MassTransferOn;
      id[nt++] = PARAM_INT;

#ifdef SMUGGLE_AGB_WINDS
      strcpy(tag[nt], "OB_MassTransferOn");
      addr[nt] = &All.OB_MassTransferOn;
      id[nt++] = PARAM_INT;
#endif

#ifdef GFM_RPROCESS
      strcpy(tag[nt], "NSNS_MassTransferOn");
      addr[nt] = &All.NSNS_MassTransferOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "NSNS_Rate_TAU");
      addr[nt] = &All.NSNS_Rate_TAU;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "NSNS_MassPerEvent");
      addr[nt] = &All.NSNS_MassPerEvent;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "NSNS_per_SNIa");
      addr[nt] = &All.NSNS_per_SNIa;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_SET_METALLICITY
      strcpy(tag[nt], "GasMetallicityInSolar");
      addr[nt] = &All.GasMetallicityInSolar;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_RPROCESS_CHANNELS
      for(int i = 0; i < GFM_RPROCESS_NSNS; i++)
        {
          sprintf(buf, "NSNS_MassPerEvent%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_MassPerEvent;
          id[nt++] = PARAM_REAL;

          sprintf(buf, "NSNS_RateNorm%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_RateNorm;
          id[nt++] = PARAM_REAL;

          sprintf(buf, "NSNS_RateTAU%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_RateTAU;
          id[nt++] = PARAM_REAL;

          sprintf(buf, "NSNS_PowerlawIndex%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_PowerlawIndex;
          id[nt++] = PARAM_REAL;

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
          sprintf(buf, "NSNS_KickAverage%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_KickAverage;
          id[nt++] = PARAM_REAL;

          sprintf(buf, "NSNS_KickSigma%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_KickSigma;
          id[nt++] = PARAM_REAL;
#endif
        }

#if GFM_RPROCESS_NSNS < GFM_RPROCESS_CHANNELS
      for(int i = GFM_RPROCESS_NSNS; i < GFM_RPROCESS_CHANNELS; i++)
        {
          sprintf(buf, "RPSN_MassPerEvent%d", i - GFM_RPROCESS_NSNS);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rpSN[i - GFM_RPROCESS_NSNS].RPSN_MassPerEvent;
          id[nt++] = PARAM_REAL;

          sprintf(buf, "RPSN_FractionPerSN%d", i - GFM_RPROCESS_NSNS);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rpSN[i - GFM_RPROCESS_NSNS].RPSN_FractionPerSN;
          id[nt++] = PARAM_REAL;
        }
#endif
#endif

      strcpy(tag[nt], "YieldTablePath");
      addr[nt] = &All.YieldTablePath;
      id[nt++] = PARAM_STRING;
#endif

#ifdef GFM_INJECT_B_FROM_SN
      strcpy(tag[nt], "SupernovaInjectedMagneticEnergyInErgs");
      addr[nt] = &All.SupernovaInjectedMagneticEnergyInErgs;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_DUST
      strcpy(tag[nt], "AGB_Dust_Delta_C");
      addr[nt] = &All.AGB_Dust_Delta_C;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGB_Dust_Delta_Metal");
      addr[nt] = &All.AGB_Dust_Delta_Metal;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SN_Dust_Delta_C");
      addr[nt] = &All.SN_Dust_Delta_C;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SN_Dust_Delta_Metal");
      addr[nt] = &All.SN_Dust_Delta_Metal;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Dust_Growth_Tau");
      addr[nt] = &All.Dust_Growth_Tau;
      id[nt++] = PARAM_REAL;

#ifndef GFM_DUST_MRN
      strcpy(tag[nt], "DustSingleGrainSize");
      addr[nt] = &All.DustSingleGrainSize;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_DUST_SPUTTERING
#if(GFM_DUST_SPUTTERING == 1)
      strcpy(tag[nt], "Dust_Sputter_Tau_Fac");
      addr[nt] = &All.Dust_Sputter_Tau_Fac;
      id[nt++] = PARAM_REAL;
#endif
#endif

#if(GFM_DUST_DESTMODE == 1)
      strcpy(tag[nt], "Dust_Destruction_Tau");
      addr[nt] = &All.Dust_Destruction_Tau;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef GFM_PREENRICH
#ifndef CHIMES_PREENRICH_AT_START
      strcpy(tag[nt], "PreEnrichTime");
      addr[nt] = &All.PreEnrichTime;
      id[nt++] = PARAM_INT;
#endif

      strcpy(tag[nt], "PreEnrichAbundanceFile");
      addr[nt] = &All.PreEnrichAbundanceFile;
      id[nt++] = PARAM_STRING;
#endif

#ifdef GFM_COOLING_METAL
      strcpy(tag[nt], "CoolingTablePath");
      addr[nt] = &All.CoolingTablePath;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "MinMetalTemp");
      addr[nt] = &All.MinMetalTemp;
      id[nt++] = PARAM_REAL;

#ifdef GFM_COOLING_METAL_START
      strcpy(tag[nt], "MetalCoolStartRedshift");
      addr[nt] = &All.MetalCoolStartRedshift;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
      strcpy(tag[nt], "PhotometricsTablePath");
      addr[nt] = &All.PhotometricsTablePath;
      id[nt++] = PARAM_STRING;
#endif

#if defined(REFINEMENT) || defined(GFM_STELLAR_EVOLUTION) || defined(BLACK_HOLES)
      strcpy(tag[nt], "ReferenceGasPartMass");
      addr[nt] = &All.ReferenceGasPartMass;
      id[nt++] = PARAM_REAL;
#endif

#if defined(REFINEMENT) && !defined(DG)
      strcpy(tag[nt], "TargetGasMassFactor");
      addr[nt] = &All.TargetGasMassFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RefinementCriterion");
      addr[nt] = &All.RefinementCriterion;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "DerefinementCriterion");
      addr[nt] = &All.DerefinementCriterion;
      id[nt++] = PARAM_INT;
#endif

#ifdef GMC_REFINEMENT
      strcpy(tag[nt], "GMCRefCellsPerJeansLength");
      addr[nt] = &All.GMCRefCellsPerJeansLength;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "GMCRefMinDensity");
      addr[nt] = &All.GMCRefMinDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "GMCRefMaxDensity");
      addr[nt] = &All.GMCRefMaxDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "GMCDerefMinDensity");
      addr[nt] = &All.GMCDerefMinDensity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef STICKYFLAGS
#if(REFLECTIVE_X == 2) || (REFLECTIVE_Y == 2) || (REFLECTIVE_Z == 2)
      strcpy(tag[nt], "StickyLayerMaxDist");
      addr[nt] = &All.StickyLayerMaxDist;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef GFM_STELLAR_FEEDBACK
      strcpy(tag[nt], "EnergyPerSNIa");
      addr[nt] = &All.EnergyPerSNIa;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGBWindVelocity");
      addr[nt] = &All.AGBWindVelocity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef USE_SFR

#ifdef QUICK_LYALPHA_LATETIMEONLY
      strcpy(tag[nt], "TimeOfCoolingStart");
      addr[nt] = &All.TimeOfCoolingStart;
      id[nt++] = PARAM_REAL;
#endif

#ifdef STEEPER_SFR_FOR_STARBURST
      strcpy(tag[nt], "StarburstPowerLawIndex");
      addr[nt] = &All.StarburstPowerLawIndex;
      id[nt++] = PARAM_REAL;
#endif

#if !defined(SMUGGLE_SFR) && !defined(LOCAL_FEEDBACK) && !defined(SFR_MCS)
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TemperatureThresh");
      addr[nt] = &All.TemperatureThresh;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = PARAM_REAL;

      /* the latter case is needed for the eEOS (and only there) */
#if !defined(GFM_STELLAR_EVOLUTION) || (defined(GFM_STELLAR_EVOLUTION) && defined(GFM_VARIABLE_IMF))
      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = PARAM_REAL;

#ifdef SOFTEREQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TempForSofterEQS");
      addr[nt] = &All.TempForSofterEQS;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MODIFIED_EOS
      strcpy(tag[nt], "FactorDensThresh");
      addr[nt] = &All.FactorDensThresh;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "FactorUthermAtThresh");
      addr[nt] = &All.FactorUthermAtThresh;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "FactorUthermJoin");
      addr[nt] = &All.FactorUthermJoin;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WINDS_LOCAL
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "WindEnergyIn1e51erg");
      addr[nt] = &All.WindEnergyIn1e51erg;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinWindVel");
      addr[nt] = &All.MinWindVel;
      id[nt++] = PARAM_REAL;
#endif

#if defined(GFM_CONST_IMF) && GFM_CONST_IMF == 1
      strcpy(tag[nt], "IMFslope");
      addr[nt] = &All.IMFslope;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WINDS_STRIPPING
      strcpy(tag[nt], "WindDumpFactor");
      addr[nt] = &All.WindDumpFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WINDS
#ifndef GFM_WINDS_VARIABLE
      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_STELLAR_EVOLUTION
      strcpy(tag[nt], "WindEnergyIn1e51erg");
      addr[nt] = &All.WindEnergyIn1e51erg;
      id[nt++] = PARAM_REAL;
#else
      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = PARAM_REAL;

#ifdef GFM_WINDS_THERMAL
      strcpy(tag[nt], "ThermalWindFactor");
      addr[nt] = &All.ThermalWindFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WINDS_THERMAL_NEWDEF
      strcpy(tag[nt], "ThermalWindFraction");
      addr[nt] = &All.ThermalWindFraction;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WINDS_VARIABLE
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "VariableWindSpecMomentum");
      addr[nt] = &All.VariableWindSpecMomentum;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinWindVel");
      addr[nt] = &All.MinWindVel;
      id[nt++] = PARAM_REAL;

#ifdef GFM_WINDS_MASSSCALING
      strcpy(tag[nt], "VariableWindMassScale");
      addr[nt] = &All.VariableWindMassScale;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GFM_WINDS_HUBBLESCALING
      strcpy(tag[nt], "WindSuppressionRedshift");
      addr[nt] = &All.WindSuppressionRedshift;
      id[nt++] = PARAM_REAL;
#endif

#if GFM_WINDS_VARIABLE == 0
      strcpy(tag[nt], "HaloConcentrationNorm");
      addr[nt] = &All.HaloConcentrationNorm;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "HaloConcentrationSlope");
      addr[nt] = &All.HaloConcentrationSlope;
      id[nt++] = PARAM_REAL;
#endif
#endif
#endif
#endif

#ifdef SMUGGLE_SFR
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DensThreshold");
      addr[nt] = &All.DensThreshold;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SfrEfficiency");
      addr[nt] = &All.SfrEfficiency;
      id[nt++] = PARAM_REAL;

#ifdef SMUGGLE_USE_POLYTROPIC_EQSTATE
      strcpy(tag[nt], "TemperatureAtThreshold");
      addr[nt] = &All.UthermThreshold;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef SMUGGLE_STAR_FEEDBACK
      strcpy(tag[nt], "FeedbackEfficiency");
      addr[nt] = &All.FeedbackEfficiency;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SMUGGLE_RADIATION_FEEDBACK
      strcpy(tag[nt], "DustOpacityRadiationFeedback");
      addr[nt] = &All.DustOpacityRadiationFeedback;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InputTimeHeatRadiationFeedback");
      addr[nt] = &All.InputTimeHeatRadiationFeedback;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InputTimeMomRadiationFeedback");
      addr[nt] = &All.InputTimeMomRadiationFeedback;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LumToMassRatioRadiationFeedback");
      addr[nt] = &All.LumToMassRatioRadiationFeedback;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RadiationFeedbackAvgPhotonEnergyineV");
      addr[nt] = &All.RadiationFeedbackAvgPhotonEnergyineV;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PhotoionizationGasTemp");
      addr[nt] = &All.PhotoionizationGasTemp;
      id[nt++] = PARAM_REAL;
#endif

#if defined(SMUGGLE_SN_COOLING_RADIUS_BOOST) || defined(SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION)
      strcpy(tag[nt], "FeedbackRadiusLimiterFactor");
      addr[nt] = &All.FeedbackRadiusLimiterFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SFR_MCS
#if SFR_MCS_SELECT_CRITERIA == 1
      strcpy(tag[nt], "SfrCritLength");
      addr[nt] = &All.SfrCritLength;
      id[nt++] = PARAM_REAL;
#elif(SFR_MCS_SELECT_CRITERIA == 2) || (SFR_MCS_SELECT_CRITERIA == 3)
      strcpy(tag[nt], "SfrCritJeansMassN");
      addr[nt] = &All.SfrCritJeansMassN;
      id[nt++] = PARAM_REAL;
#endif
#if(SFR_MCS_SELECT_CRITERIA == 0) || (SFR_MCS_SELECT_CRITERIA == 3)
      strcpy(tag[nt], "DensThreshold");
      addr[nt] = &All.DensThreshold;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "SfrEfficiency");
      addr[nt] = &All.SfrEfficiency;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = PARAM_REAL;

#if defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)
#ifdef SFR_MCS_DELAY
      strcpy(tag[nt], "TimeDelayFactor");
      addr[nt] = &All.TimeDelayFactor;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "MaxFBStarEvalTimestep");
      addr[nt] = &All.MaxFBStarEvalTimestep;
      id[nt++] = PARAM_REAL;

#ifndef IMF_SAMPLING_MCS
      strcpy(tag[nt], "MaxFBStarEvalTimestepCut");
      addr[nt] = &All.MaxFBStarEvalTimestepCut;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "MaxFBStarMassFac");
      addr[nt] = &All.MaxFBStarMassFac;
      id[nt++] = PARAM_REAL;
#endif  // defined(SN_MCS) || defined(HII_MCS) || defined(PE_MCS)

#if((defined(SN_MCS) && !defined(SN_MCS_SINGLE_INJECTION)) || (defined(HII_MCS) && !defined(HII_MCS_TEST)) || defined(PE_MCS)) && \
    !(defined(SN_MCS_INITIAL_DRIVING) || defined(IMF_SAMPLING_MCS))
      strcpy(tag[nt], "SB99TablesPath");
      addr[nt] = &All.SB99TablesPath;
      id[nt++] = PARAM_STRING;

#ifdef SB99_FIXED_Z
      strcpy(tag[nt], "SB99_metallicity_name");
      addr[nt] = &All.SB99_metallicity_name;
      id[nt++] = PARAM_STRING;
#endif
#endif

#ifdef SN_MCS
#ifndef IMF_SAMPLING_MCS
#ifndef SN_MCS_INITIAL_DRIVING
#ifdef SN_MCS_SINGLE_INJECTION
      strcpy(tag[nt], "SNDelay");
      addr[nt] = &All.SNDelay;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SNMassUnit");
      addr[nt] = &All.SNMassUnit;
      id[nt++] = PARAM_REAL;
#endif
#else
      strcpy(tag[nt], "DrivingZoneRadius");
      addr[nt] = &All.DrivingZoneRadius;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifndef SN_MCS_HOST_ONLY
      strcpy(tag[nt], "HostCellFeedbackFraction");
      addr[nt] = &All.HostCellFeedbackFraction;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "SNKineticRatio");
      addr[nt] = &All.SNKineticRatio;
      id[nt++] = PARAM_REAL;

#if !defined(IMF_SAMPLING_MCS) && !defined(SN_MCS_VARIABLE_EJECTA)
      strcpy(tag[nt], "SNMassReturn");
      addr[nt] = &All.SNMassReturn;
      id[nt++] = PARAM_REAL;
#ifdef METALS
      strcpy(tag[nt], "SNEjectaMetallicity");
      addr[nt] = &All.SNEjectaMetallicity;
      id[nt++] = PARAM_REAL;
#endif
#endif
#endif  // SN_MCS

#ifdef HII_MCS
      strcpy(tag[nt], "PhotoionizationGasTemp");
      addr[nt] = &All.PhotoionizationGasTemp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "R_Stromgren_Max");
      addr[nt] = &All.R_Stromgren_Max;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinimumPhotoRateFactor");
      addr[nt] = &All.MinimumPhotoRateFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "HiiAttenFac");
      addr[nt] = &All.HiiAttenFac;
      id[nt++] = PARAM_REAL;

#ifndef HII_MCS_EVERY_SYNCPOINT
      strcpy(tag[nt], "TimeBetweenHiiPlace");
      addr[nt] = &All.TimeBetweenHiiPlace;
      id[nt++] = PARAM_REAL;
#endif
#ifdef HII_MCS_DENSCUT
      strcpy(tag[nt], "HiiDensCut");
      addr[nt] = &All.HiiDensCut;
      id[nt++] = PARAM_REAL;
#endif
#ifdef HII_MCS_LR
      strcpy(tag[nt], "UVBEnergyDensHii");
      addr[nt] = &All.UVBEnergyDensHii;
      id[nt++] = PARAM_REAL;
#endif
#endif  // HII_MCS

#ifdef PE_MCS
      strcpy(tag[nt], "G_min");
      addr[nt] = &All.G_min;
      id[nt++] = PARAM_REAL;

#ifdef PE_MCS_FIXED_EPS
      strcpy(tag[nt], "PhotoelectricHeatingEps");
      addr[nt] = &All.PhotoelectricHeatingEps;
      id[nt++] = PARAM_REAL;
#endif
#endif  // PE_MCS

#ifdef IMF_SAMPLING_MCS
      strcpy(tag[nt], "MinimumIMFStarMass");
      addr[nt] = &All.MinimumIMFStarMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaximumIMFStarMass");
      addr[nt] = &All.MaximumIMFStarMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinimumImportantStellarMass");
      addr[nt] = &All.MinimumImportantStellarMass;
      id[nt++] = PARAM_REAL;
#ifdef SN_MCS
      strcpy(tag[nt], "SNStarMinMass");
      addr[nt] = &All.SNStarMinMass;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SNStarMaxMass");
      addr[nt] = &All.SNStarMaxMass;
      id[nt++] = PARAM_REAL;

#ifdef SN_MCS_PROMPT
      strcpy(tag[nt], "SNLifetime");
      addr[nt] = &All.SNLifetime;
      id[nt++] = PARAM_REAL;
#endif
#endif
#endif  // IMF_SAMPLING_MCS
#endif  // SFR_MCS

#ifdef TURB_APPROX_MCS
      strcpy(tag[nt], "MinTurbSpecEnergy");
      addr[nt] = &All.MinTurbSpecEnergy;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SMAUG_PRESSURE_FLOOR
      strcpy(tag[nt], "Polytrope_Tstar");
      addr[nt] = &All.Polytrope_Tstar;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Polytrope_nstar");
      addr[nt] = &All.Polytrope_nstar;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Polytrope_gstar");
      addr[nt] = &All.Polytrope_gstar;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DARKENERGY
#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = PARAM_STRING;
#endif
#endif

#ifdef RELAXOBJECT
      strcpy(tag[nt], "RelaxBaseFac");
      addr[nt] = &All.RelaxBaseFac;
      id[nt++] = PARAM_REAL;
#endif

#ifdef RELAXOBJECT_COOLING
      strcpy(tag[nt], "RelaxTemperature");
      addr[nt] = &All.RelaxTemperature;
      id[nt++] = PARAM_REAL;
#endif

#ifdef RELAXOBJECT_COOLING2
      strcpy(tag[nt], "TempCore");
      addr[nt] = &All.TempCore;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TempShell");
      addr[nt] = &All.TempShell;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MassCore");
      addr[nt] = &All.MassCore;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ShellBaseRadius");
      addr[nt] = &All.ShellBaseRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TempMin");
      addr[nt] = &All.TempMin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RTempMin");
      addr[nt] = &All.RTempMin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BackgroundTemp");
      addr[nt] = &All.BackgroundTemp;
      id[nt++] = PARAM_REAL;
#endif

#ifdef INSPIRAL
      strcpy(tag[nt], "InspiralVelocity");
      addr[nt] = &All.InspiralVelocity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef EOS_DEGENERATE
      strcpy(tag[nt], "EosTable");
      addr[nt] = All.EosTable;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "EosSpecies");
      addr[nt] = All.EosSpecies;
      id[nt++] = PARAM_STRING;
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(tag[nt], "NetworkRates");
      addr[nt] = All.NetworkRates;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "NetworkPartFunc");
      addr[nt] = All.NetworkPartFunc;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "NetworkMasses");
      addr[nt] = All.NetworkMasses;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "NetworkWeakrates");
      addr[nt] = All.NetworkWeakrates;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "NetworkTempThreshold");
      addr[nt] = &All.NetworkTempThreshold;
      id[nt++] = PARAM_REAL;

#ifdef NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE
      strcpy(tag[nt], "NuclearNetworkMaxCompositionChange");
      addr[nt] = &All.NuclearNetworkMaxCompositionChange;
      id[nt++] = PARAM_REAL;
#else
      strcpy(tag[nt], "NuclearNetworkMaxTempChange");
      addr[nt] = &All.NuclearNetworkMaxTempChange;
      id[nt++] = PARAM_REAL;
#endif

#ifdef NUCLEAR_NETWORK_TIMESTEP_LIMITER
      strcpy(tag[nt], "NuclearNetworkMaxEnergyDiff");
      addr[nt] = &All.NuclearNetworkMaxEnergyDiff;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef NETWORK_NSE
      strcpy(tag[nt], "NetworkNSEThreshold");
      addr[nt] = &All.NetworkNSEThreshold;
      id[nt++] = PARAM_REAL;
#endif

#ifdef EOS_OPAL
      strcpy(tag[nt], "EosOpalTable");
      addr[nt] = All.EosOpalTable;
      id[nt++] = PARAM_STRING;
#endif

#ifdef TRACER_TRAJECTORY
#ifdef TRACER_TRAJECTORY_GENERATE
      strcpy(tag[nt], "NumberOfTracersToGenerate");
      addr[nt] = &All.NumberOfTracersToGenerate;
      id[nt++] = PARAM_INT;
#else
      strcpy(tag[nt], "TracerInitFile");
      addr[nt] = All.TracerInitFile;
      id[nt++] = PARAM_STRING;
#endif

      strcpy(tag[nt], "TracerOutputFile");
      addr[nt] = All.TracerOutputFile;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "TracerOutputConfFile");
      addr[nt] = All.TracerOutputConfFile;
      id[nt++] = PARAM_STRING;
#endif

#ifdef MHD_CT
      strcpy(tag[nt], "CTMeanBx");
      addr[nt] = &All.CT_mean_Bx;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CTMeanBy");
      addr[nt] = &All.CT_mean_By;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CTMeanBz");
      addr[nt] = &All.CT_mean_Bz;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MHD_SEEDFIELD
      strcpy(tag[nt], "MHDSeedDir");
      addr[nt] = &All.B_dir;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MHDSeedValue");
      addr[nt] = &All.B_value;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MHD_SEEDPSPEC
      strcpy(tag[nt], "MHDSeedPSpecSlope");
      addr[nt] = &All.B_pspec_slope;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MHDSeedPSpecAmpl");
      addr[nt] = &All.B_pspec_ampl;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MHDSeedPSpecKCut");
      addr[nt] = &All.B_pspec_kcut;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MHDSeedPSpecHelical");
      addr[nt] = &All.B_pspec_helical;
      id[nt++] = PARAM_REAL;
#endif

#ifdef TGSET
      strcpy(tag[nt], "MaxTimeBinDiff");
      addr[nt] = &TGD.MaxTimeBinDiff;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SnapFreeFallTimeFac");
      addr[nt] = &TGD.SnapFreeFallTimeFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "NHInit");
      addr[nt] = &TGD.NHInit;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "NHTerm");
      addr[nt] = &TGD.NHTerm;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "StreamVel");
      addr[nt] = &TGD.StreamVel;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "JeansTemp");
      addr[nt] = &TGD.JeansTemp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "JeansDensLow");
      addr[nt] = &TGD.JeansDensLow;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "JeansDensHigh");
      addr[nt] = &TGD.JeansDensHigh;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "JeansNumberLow");
      addr[nt] = &TGD.JeansNumberLow;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "JeansNumberHigh");
      addr[nt] = &TGD.JeansNumberHigh;
      id[nt++] = PARAM_REAL;
#endif

#ifdef TGCHEM
      strcpy(tag[nt], "ChemMode");
      addr[nt] = &TGCD.ChemMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ChemIOMode");
      addr[nt] = &TGCD.ChemIOMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ChemH2Mode");
      addr[nt] = &TGCD.ChemH2Mode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ChemInitAbH2");
      addr[nt] = &TGCD.ChemInitAbH2;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ChemInitAbHII");
      addr[nt] = &TGCD.ChemInitAbHII;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ChemJ21");
      addr[nt] = &TGCD.ChemJ21;
      id[nt++] = PARAM_REAL;
#endif

#ifdef HEALRAY
      strcpy(tag[nt], "RaySourceFile");
      addr[nt] = &HRD.RaySourceFile;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "RayMode");
      addr[nt] = &HRD.RayMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayIOMode");
      addr[nt] = &HRD.RayIOMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayNumThreads");
      addr[nt] = &HRD.RayNumThreads;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayDebugMode");
      addr[nt] = &HRD.RayDebugMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayTimeFac");
      addr[nt] = &HRD.RayTimeFac;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayMultiFreqMode");
      addr[nt] = &HRD.RayMultiFreqMode;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayNumSources");
      addr[nt] = &HRD.RayNumSources;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayNumTrace");
      addr[nt] = &HRD.RayNumTrace;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayLvlInit");
      addr[nt] = &HRD.RayLvlInit;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayNumLines");
      addr[nt] = &HRD.RayNumLines;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RayNumFreq");
      addr[nt] = &HRD.RayNumFreq;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "RaySplitFac");
      addr[nt] = &HRD.RaySplitFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RayPartFrac");
      addr[nt] = &HRD.RayPartFrac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RayMinEgyFrac");
      addr[nt] = &HRD.RayMinEgyFrac;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SGCHEM
      strcpy(tag[nt], "SGChemConstInitAbundances");
      addr[nt] = &All.SGChemConstInitAbundances;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SGChemInitH2Abund");
      addr[nt] = &All.SGChemInitH2Abund;
      id[nt++] = PARAM_REAL;

      mpi_printf("READING chemical variable parameters\n");

      strcpy(tag[nt], "SGChemInitHPAbund");
      addr[nt] = &All.SGChemInitHPAbund;
      id[nt++] = PARAM_REAL;

#if CHEMISTRYNETWORK == 1
      strcpy(tag[nt], "SGChemInitDIIAbund");
      addr[nt] = &All.SGChemInitDIIAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitHDAbund");
      addr[nt] = &All.SGChemInitHDAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitHeIIIAbund");
      addr[nt] = &All.SGChemInitHeIIIAbund;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "SGChemInitCPAbund");
      addr[nt] = &All.SGChemInitCPAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitCOAbund");
      addr[nt] = &All.SGChemInitCOAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitCHxAbund");
      addr[nt] = &All.SGChemInitCHxAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitOHxAbund");
      addr[nt] = &All.SGChemInitOHxAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitHCOPAbund");
      addr[nt] = &All.SGChemInitHCOPAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitHePAbund");
      addr[nt] = &All.SGChemInitHePAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SGChemInitMPAbund");
      addr[nt] = &All.SGChemInitMPAbund;
      id[nt++] = PARAM_REAL;

#ifndef SGCHEM_VARIABLE_Z
      strcpy(tag[nt], "CarbAbund");
      addr[nt] = &All.CarbAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OxyAbund");
      addr[nt] = &All.OxyAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MAbund");
      addr[nt] = &All.MAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ZAtom");
      addr[nt] = &All.ZAtom;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "DeutAbund");
      addr[nt] = &All.DeutAbund;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InitDustTemp");
      addr[nt] = &All.InitDustTemp;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "UVFieldStrength");
      addr[nt] = &All.UVFieldStrength;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LWBGType");
      addr[nt] = &All.LWBGType;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "LWBGStartRedsh");
      addr[nt] = &All.LWBGStartRedsh;
      id[nt++] = PARAM_REAL;

#ifndef SGCHEM_VARIABLE_Z
      strcpy(tag[nt], "DustToGasRatio");
      addr[nt] = &All.DustToGasRatio;
      id[nt++] = PARAM_REAL;
#endif

#ifndef SGCHEM_VARIABLE_CRION
      strcpy(tag[nt], "CosmicRayIonRate");
      addr[nt] = &All.CosmicRayIonRate;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "InitRedshift");
      addr[nt] = &All.InitRedshift;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ExternalDustExtinction");
      addr[nt] = &All.ExternalDustExtinction;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "H2FormEx");
      addr[nt] = &All.H2FormEx;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "H2FormKin");
      addr[nt] = &All.H2FormKin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PhotoApprox");
      addr[nt] = &All.PhotoApprox;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ISRFOption");
      addr[nt] = &All.ISRFOption;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "AtomicCoolOption");
      addr[nt] = &All.AtomicCoolOption;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "H2OpacityOption");
      addr[nt] = &All.H2OpacityOption;
      id[nt++] = PARAM_INT;

#ifdef SGCHEM_ACCRETION_LUMINOSITY
      strcpy(tag[nt], "SGChemAccretionLuminosityOn");
      addr[nt] = &All.SGChemAccretionLuminosityOn;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SinkAccretionRateSmoothingMass");
      addr[nt] = &All.SinkAccretionRateSmoothingMass;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SGCHEM_TEMPERATURE_FLOOR
      strcpy(tag[nt], "SGChemTemperatureFloor");
      addr[nt] = &All.SGChemTemperatureFloor;
      id[nt++] = PARAM_REAL;
#endif

#endif /*SGCHEM*/

#ifdef GRADIENTREFINEMENT
      strcpy(tag[nt], "GradVelRefinement");
      addr[nt] = &All.GradVelRefinement;
      id[nt++] = PARAM_INT;
      strcpy(tag[nt], "GradBfldRefinement");
      addr[nt] = &All.GradBfldRefinement;
      id[nt++] = PARAM_INT;
#endif

#ifdef TREECOLV2
      strcpy(tag[nt], "TreeColMaxDistance");
      addr[nt] = &All.TreeColMaxDistance;
      id[nt++] = PARAM_REAL;
#endif

#ifdef TREECOLV2_VEL
      strcpy(tag[nt], "FracOverlap"); /* mean line overlap following Hartwig et al. 2015, ApJ, 799, 114 */
      addr[nt] = &All.FracOverlap;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SINKS
      strcpy(tag[nt], "NHThresh");
      addr[nt] = &SKD.NHThresh;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AccRad");
      addr[nt] = &SKD.AccRad;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SNE_FEEDBACK
      strcpy(tag[nt], "SNEInjectionCriterion");
      addr[nt] = &All.SNEInjectionCriterion;
      id[nt++] = PARAM_INT;
      strcpy(tag[nt], "SNESeed");
      addr[nt] = &All.SNESeed;
      id[nt++] = PARAM_INT;
      strcpy(tag[nt], "SNETargetMass");
      addr[nt] = &All.SNETargetMass;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "SNEPeriodInYears");
      addr[nt] = &All.SNEPeriodInYears;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "SNEMinimalParticleNumber");
      addr[nt] = &All.SNEMinimalParticleNumber;
      id[nt++] = PARAM_INT;

#ifdef CLUSTERED_SNE
      strcpy(tag[nt], "SNEClusterTfin");
      addr[nt] = &All.SNEClusterTfin;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "SNEClusterTinit");
      addr[nt] = &All.SNEClusterTinit;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "SNENumber");
      addr[nt] = &All.SNENumber;
      id[nt++] = PARAM_INT;
#endif

#ifdef INJECT_TRACER_INTO_SN
      strcpy(tag[nt], "SNETracersForEachSn");
      addr[nt] = &All.SNETracersForEachSn;
      id[nt++] = PARAM_INT;
      strcpy(tag[nt], "SNETracerBitmask");
      addr[nt] = &All.SNETracerBitmask;
      id[nt++] = PARAM_INT;
#endif

#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
      strcpy(tag[nt], "SNEScatterAroundSink");
      addr[nt] = &All.SNEScatterAroundSink;
      id[nt++] = PARAM_REAL;
#endif

#endif /*SNE_FEEDBACK*/

#ifdef REFINE_ONLY_WITH_TRACER
      strcpy(tag[nt], "MaxTracerVolume");
      addr[nt] = &All.MaxTracerVolume;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "MinTracerVolume");
      addr[nt] = &All.MinTracerVolume;
      id[nt++] = PARAM_REAL;
#endif

#ifdef ROTATING_HIGHRES_REGION
      strcpy(tag[nt], "Highres_x0");
      addr[nt] = &All.Highres_x0;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_y0");
      addr[nt] = &All.Highres_y0;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_deltaR");
      addr[nt] = &All.Highres_deltaR;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_deltaThetaxR");
      addr[nt] = &All.Highres_deltaThetaxR;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_deltaz");
      addr[nt] = &All.Highres_deltaz;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_vrot");
      addr[nt] = &All.Highres_vrot;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_t0");
      addr[nt] = &All.Highres_t0;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Highres_targetmass");
      addr[nt] = &All.Highres_targetmass;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SINK_PARTICLES
      strcpy(tag[nt], "SinkCreationDensityCodeUnits");
      addr[nt] = &SinkCreationDensityCodeUnits;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SinkFormationRadius");
      addr[nt] = &SinkFormationRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SinkEvolutionDumpRateYears");
      addr[nt] = &SinkEvolutionDumpRateYears;
      id[nt++] = PARAM_REAL;

#ifdef SINK_PARTICLES_FEEDBACK
      strcpy(tag[nt], "SINKStarFormationEfficiency");
      addr[nt] = &All.SINKStarFormationEfficiency;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxStellarMassPerSink");
      addr[nt] = &All.MaxStellarMassPerSink;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SINK_PHOTOION_FEEDBACK
      strcpy(tag[nt], "StromgrenTemperature");
      addr[nt] = &All.StromgrenTemperature;
      id[nt++] = PARAM_REAL;
#endif

#endif /*SINK_PARTICLES*/

#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
      strcpy(tag[nt], "DynamicViscosity");
      addr[nt] = &All.dyn_visc;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BulkViscosity");
      addr[nt] = &All.bulk_visc;
      id[nt++] = PARAM_REAL;
#else
#ifdef USE_KINEMATIC_VISCOSITY
      strcpy(tag[nt], "KinematicViscosity");
      addr[nt] = &All.KinematicViscosity;
      id[nt++] = PARAM_REAL;
#else
#ifdef ALPHA_VISCOSITY
      strcpy(tag[nt], "AlphaCoefficient");
      addr[nt] = &All.AlphaCoefficient;
      id[nt++] = PARAM_REAL;
#endif
#endif
#endif
#endif

#ifdef THERMAL_CONDUCTION
      strcpy(tag[nt], "ThermalConductivity");
      addr[nt] = &All.ThermalConductivity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef TRACER_DIFFUSION
      strcpy(tag[nt], "TracerDiffusivity");
      addr[nt] = &All.TracerDiffusivity;
      id[nt++] = PARAM_REAL;
#endif

#ifdef TRACER_MC
      strcpy(tag[nt], "TracerMCPerCell");
      addr[nt] = &All.TracerMCPerCell;
      id[nt++] = PARAM_INT;
#endif

#ifdef TRACER_PARTICLE
      strcpy(tag[nt], "MinimumTracerHsml");
      addr[nt] = &All.MinimumTracerHsml;
      id[nt++] = PARAM_REAL;
#endif

#ifdef GRAVITY_TABLE
      strcpy(tag[nt], "ExternalGravForcesFile");
      addr[nt] = All.ExternalGravForcesFile;
      id[nt++] = PARAM_STRING;
#endif

#ifdef SPECIAL_BOUNDARY
      strcpy(tag[nt], "BoundaryLayerScaleFactor");
      addr[nt] = &All.BoundaryLayerScaleFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SpecialBoundarySpeed");
      addr[nt] = &All.SpecialBoundarySpeed;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SpecialBoundaryMotion");
      addr[nt] = &All.SpecialBoundaryMotion;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "SpecialBoundaryType");
      addr[nt] = &All.SpecialBoundaryType;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "OutflowPressure");
      addr[nt] = &All.OutflowPressure;
      id[nt++] = PARAM_REAL;
#ifdef THERMAL_CONDUCTION
      strcpy(tag[nt], "BoundaryTemperature");
      addr[nt] = &All.BoundaryTemperature;
      id[nt++] = PARAM_REAL;
#endif
#endif

#if defined(COAXIAL_BOUNDARIES) && !defined(CIRCUMSTELLAR)
      strcpy(tag[nt], "InnerRadius");
      addr[nt] = &All.inner_radius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OuterRadius");
      addr[nt] = &All.outer_radius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "InnerAngVel");
      addr[nt] = &All.omega_in;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "OuterAngVel");
      addr[nt] = &All.omega_out;
      id[nt++] = PARAM_REAL;
#endif

#ifdef AMR
      strcpy(tag[nt], "MinRefLevel");
      addr[nt] = &All.MinRefLevel;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MaxRefLevel");
      addr[nt] = &All.MaxRefLevel;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MeshSmoothing");
      addr[nt] = &All.AMRMeshSmoothing;
      id[nt++] = PARAM_INT;
#endif

#ifdef AB_TURB
      strcpy(tag[nt], "ST_decay");
      addr[nt] = &All.StDecay;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_energy");
      addr[nt] = &All.StEnergy;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_DtFreq");
      addr[nt] = &All.StDtFreq;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_Kmin");
      addr[nt] = &All.StKmin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_Kmax");
      addr[nt] = &All.StKmax;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_SolWeight");
      addr[nt] = &All.StSolWeight;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_AmplFac");
      addr[nt] = &All.StAmplFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ST_Seed");
      addr[nt] = &All.StSeed;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "ST_SpectForm");
      addr[nt] = &All.StSpectForm;
      id[nt++] = PARAM_INT;
#endif

#ifdef PREHEATING
      strcpy(tag[nt], "TimePreheating");
      addr[nt] = &All.TimePreheating;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TempPreheating");
      addr[nt] = &All.TempPreheating;
      id[nt++] = PARAM_REAL;
#endif

#ifdef ADJ_BOX_POWERSPEC
      strcpy(tag[nt], "BoxWidth");
      addr[nt] = &All.BoxWidth;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BoxCenter_x");
      addr[nt] = &All.BoxCenter_x;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BoxCenter_y");
      addr[nt] = &All.BoxCenter_y;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "BoxCenter_z");
      addr[nt] = &All.BoxCenter_z;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "TransformSize");
      addr[nt] = &All.FourierGrid;
      id[nt++] = PARAM_INT;
#endif

#ifdef VEL_POWERSPEC_BOX
      strcpy(tag[nt], "PSCenterX");
      addr[nt] = &All.PSCenterX;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PSCenterY");
      addr[nt] = &All.PSCenterY;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PSCenterZ");
      addr[nt] = &All.PSCenterZ;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PSRadius");
      addr[nt] = &All.PSRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "PSMinRadius");
      addr[nt] = &All.PSMinRadius;
      id[nt++] = PARAM_REAL;
#endif

#ifdef ATOMIC_DM
      strcpy(tag[nt], "ADMProtonMassInkeV");
      addr[nt] = &All.ADMProtonMassInkeV;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ADMElectronMassInkeV");
      addr[nt] = &All.ADMElectronMassInkeV;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ADMFineStructure");
      addr[nt] = &All.ADMFineStructureConstant;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SIDM
      strcpy(tag[nt], "DtimeFac");
      addr[nt] = &All.DtimeFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DtimeFacLim");
      addr[nt] = &All.DtimeFacLim;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SIDMDesNumNgb");
      addr[nt] = &All.SIDMDesNumNgb;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SIDMMaxNumNgbDeviation");
      addr[nt] = &All.SIDMMaxNumNgbDeviation;
      id[nt++] = PARAM_REAL;
#ifdef SIDM_CONST_CROSS
      strcpy(tag[nt], "CrossSectionPerMass_in_cgs");
      addr[nt] = &All.CrossSectionPerMass_in_cgs;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef DUST_LIVE
      strcpy(tag[nt], "DesNumNgbDust");
      addr[nt] = &All.DesNumNgbDust;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "MaxNumNgbDeviationDust");
      addr[nt] = &All.MaxNumNgbDeviationDust;
      id[nt++] = PARAM_REAL;

#if !defined(DL_DRAG_SEMI_IMPLICIT) && !defined(DL_NODRAG)
      strcpy(tag[nt], "StoppingTimeFrac");
      addr[nt] = &All.StoppingTimeFrac;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DL_GRAIN_BINS
      strcpy(tag[nt], "MinGrainSize");
      addr[nt] = &All.MinGrainSize;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxGrainSize");
      addr[nt] = &All.MaxGrainSize;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "GrainDensity");
      addr[nt] = &All.GrainDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxBinFracMassChg");
      addr[nt] = &All.MaxBinFracMassChg;
      id[nt++] = PARAM_REAL;

#if defined(DL_SHATTERING) || defined(DL_COAGULATION) || defined(DL_PRODUCTION)
      strcpy(tag[nt], "GrainDataPath");
      addr[nt] = &All.GrainDataPath;
      id[nt++] = PARAM_STRING;
#endif

#ifdef DL_PRODUCTION
      strcpy(tag[nt], "DustTargetFrac");
      addr[nt] = &All.DustTargetFrac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "NumDustPerSpawn");
      addr[nt] = &All.NumDustPerSpawn;
      id[nt++] = PARAM_INT;

#ifdef DL_REFINEMENT
      strcpy(tag[nt], "DustMaxFrac");
      addr[nt] = &All.DustMaxFrac;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DL_DEREFINEMENT
      strcpy(tag[nt], "DustMinFrac");
      addr[nt] = &All.DustMinFrac;
      id[nt++] = PARAM_REAL;
#endif
#endif /* DL_PRODUCTION */

#ifdef DL_SUBCYCLE
      strcpy(tag[nt], "DustSubcycleFac");
      addr[nt] = &All.DustSubcycleFac;
      id[nt++] = PARAM_REAL;
#endif
#endif /* DL_GRAIN_BINS */

#ifdef DL_RADIATION
      strcpy(tag[nt], "GrainRPPath");
      addr[nt] = &All.GrainRPPath;
      id[nt++] = PARAM_STRING;
#endif
#endif /* DUST_LIVE */

#ifdef REFINEMENT_VOLUME_LIMIT
      strcpy(tag[nt], "MaxVolumeDiff");
      addr[nt] = &All.MaxVolumeDiff;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinVolume");
      addr[nt] = &All.MinVolume;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxVolume");
      addr[nt] = &All.MaxVolume;
      id[nt++] = PARAM_REAL;
#endif

#ifdef REFINEMENT_LIMIT_STARFORMING_GAS
      strcpy(tag[nt], "HighDensityMaxGasDerefinementFactor");
      addr[nt] = &All.HighDensityMaxGasDerefinementFactor;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "HighDensityMaxGasDerefinementFactor");
      addr[nt] = &All.HighDensityMaxGasDerefinementFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef REFINEMENT_BY_DENSITY
      strcpy(tag[nt], "MinimumDensityForRefinement");
      addr[nt] = &All.MinimumDensityForRefinement;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinimumVolumeForDensityRefinement");
      addr[nt] = &All.MinimumVolumeForDensityRefinement;
      id[nt++] = PARAM_REAL;
#endif

#ifdef TILE_ICS
      strcpy(tag[nt], "TileICsFactor");
      addr[nt] = &All.TileICsFactor;
      id[nt++] = PARAM_INT;
#endif

#ifdef OTVET
      strcpy(tag[nt], "IonizingLumPerSolarMass");
      addr[nt] = &All.IonizingLumPerSolarMass;
      id[nt++] = PARAM_REAL;
#ifdef OTVET_MULTI_FREQUENCY
      strcpy(tag[nt], "star_Teff");
      addr[nt] = &All.star_Teff;
      id[nt++] = PARAM_REAL;
#endif
#ifdef OTVET_SCATTER_SOURCE
      strcpy(tag[nt], "OtvetNumNgbSource");
      addr[nt] = &All.OtvetNumNgbSource;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "OtvetMaxNumNgbDeviationSource");
      addr[nt] = &All.OtvetMaxNumNgbDeviationSource;
      id[nt++] = PARAM_REAL;
#endif
#endif

#if defined(CONDUCTION) || defined(MONOTONE_CONDUCTION)
      strcpy(tag[nt], "ConductionEfficiency");
      addr[nt] = &All.ConductionEfficiency;
      id[nt++] = PARAM_REAL;
#endif

#ifdef CONDUCTION
      strcpy(tag[nt], "MaxSizeConductionStep");
      addr[nt] = &All.MaxSizeConductionStep;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MONOTONE_CONDUCTION
      strcpy(tag[nt], "MaxConductionSubCycles");
      addr[nt] = &All.MaxConductionSubCycles;
      id[nt++] = PARAM_INT;

#ifdef RESTRICT_KAPPA
      strcpy(tag[nt], "MaxDiffusivity");
      addr[nt] = &All.MaxDiffusivity;
      id[nt++] = PARAM_REAL;
#endif

#endif

#ifdef SMUGGLE_RADPRESS_OPT_THIN
      strcpy(tag[nt], "RadPressure_MaxAge");
      addr[nt] = &All.RadPressure_MaxAge;
      id[nt++] = PARAM_REAL;
#endif

#ifdef ADDBACKGROUNDGRID
      strcpy(tag[nt], "GridSize");
      addr[nt] = &All.GridSize;
      id[nt++] = PARAM_INT;
#endif

#if defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
      // standard parameters
      All.RayMemFac   = 1.5;
      All.MachMin     = 1.3;
      All.RayStepsMax = 5;
#endif

#ifdef SHOCK_FINDER_POST_PROCESSING
      strcpy(tag[nt], "RayMemFac");
      addr[nt] = &All.RayMemFac;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MachMin");
      addr[nt] = &All.MachMin;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "RayStepsMax");
      addr[nt] = &All.RayStepsMax;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "NumFilesPerOutput");
      addr[nt] = &All.NumFilesPerOutput;
      id[nt++] = PARAM_INT;

      strcpy(tag[nt], "OutputDirShockFinder");
      addr[nt] = All.OutputDirShockFinder;
      id[nt++] = PARAM_STRING;
#endif

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
#ifdef SKIP_BORDER
      strcpy(tag[nt], "DistToBorder");
      addr[nt] = &All.DistToBorder;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SHOCK_SUBVOLUME
      strcpy(tag[nt], "ShockSubvolumeNum");
      addr[nt] = &All.ShockSubvolumeNum;
      id[nt++] = PARAM_INT;
#endif
#endif

#ifdef DG
#if defined(CHARACTERISTIC_LIMITER) || defined(CONSERVED_LIMITER)
      strcpy(tag[nt], "DG_beta");
      addr[nt] = &All.DG_beta;
      id[nt++] = PARAM_REAL;
#endif

#ifdef DISCONTINUITY_DETECTION
      strcpy(tag[nt], "DG_alpha");
      addr[nt] = &All.DG_alpha;
      id[nt++] = PARAM_REAL;
#endif

#ifdef ANGLE_BOUND
      strcpy(tag[nt], "DG_min_angle");
      addr[nt] = &All.DG_min_angle;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MINMOD_B
      strcpy(tag[nt], "DG_M");
      addr[nt] = &All.DG_M;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MACHNUM_B
      strcpy(tag[nt], "DG_lim_mach_min");
      addr[nt] = &All.DG_lim_mach_min;
      id[nt++] = PARAM_REAL;
#endif

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
      strcpy(tag[nt], "DG_TargetSlope");
      addr[nt] = &All.DG_TargetSlope;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DG_SlopeRangeFactor");
      addr[nt] = &All.DG_SlopeRangeFactor;
      id[nt++] = PARAM_REAL;
#endif

#ifdef RK2
      All.DG_RK2_alpha = 1;  // use SSP RK2 Heun's method as a standard
#endif
#endif

#ifdef COSMIC_RAYS
      strcpy(tag[nt], "GammaCR");
      addr[nt] = &All.GammaCR;
      id[nt++] = PARAM_REAL;

#ifdef COSMIC_RAYS_SN_INJECTION
      strcpy(tag[nt], "CREnergyInputPerSolarMassOfStarFormation");
      addr[nt] = &All.CREnergyInputPerSolarMassOfStarFormation;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "DesNumNgbCRInjection");
      addr[nt] = &All.DesNumNgbCRInjection;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MaxNumNgbDeviationCRInjection");
      addr[nt] = &All.MaxNumNgbDeviationCRInjection;
      id[nt++] = PARAM_REAL;
#endif

      strcpy(tag[nt], "MinimumCREnergyDensity");
      addr[nt] = &All.MinimumCREnergyDensity;
      id[nt++] = PARAM_REAL;

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
      strcpy(tag[nt], "AccelerationEfficiency");
      addr[nt] = &All.AccelerationEfficiency;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CriticalMachnumber");
      addr[nt] = &All.CriticalMachnumber;
      id[nt++] = PARAM_REAL;
#endif

#ifdef COSMIC_RAYS_DIFFUSION
#ifdef COSMIC_RAYS_DIFFUSION_CONSTANT_TIMESTEP
      strcpy(tag[nt], "CRDiffusionTimestep");
      addr[nt] = &All.CRDiffusionTimestep;
      id[nt++] = PARAM_REAL;
#else
      strcpy(tag[nt], "CourantCRDiffusion");
      addr[nt] = &All.CourantCRDiffusion;
      id[nt++] = PARAM_REAL;
#endif
      strcpy(tag[nt], "CRDiffusionCoefficient");
      addr[nt] = &All.CR_Diffusion_Coefficient;
      id[nt++] = PARAM_REAL;

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_X
      strcpy(tag[nt], "CR_Boundary_X_Upper");
      addr[nt] = &All.CR_Boundary_X_Upper;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CR_Boundary_X_Lower");
      addr[nt] = &All.CR_Boundary_X_Lower;
      id[nt++] = PARAM_REAL;
#endif

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_Y
      strcpy(tag[nt], "CR_Boundary_Y_Upper");
      addr[nt] = &All.CR_Boundary_Y_Upper;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CR_Boundary_Y_Lower");
      addr[nt] = &All.CR_Boundary_Y_Lower;
      id[nt++] = PARAM_REAL;
#endif

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_Z
      strcpy(tag[nt], "CR_Boundary_Z_Upper");
      addr[nt] = &All.CR_Boundary_Z_Upper;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CR_Boundary_Z_Lower");
      addr[nt] = &All.CR_Boundary_Z_Lower;
      id[nt++] = PARAM_REAL;
#endif

#endif

#ifdef COSMIC_RAYS_STREAMING_EXPLICIT
      strcpy(tag[nt], "CR_Chi");
      addr[nt] = &All.CR_Chi;
      id[nt++] = PARAM_REAL;
#endif

#endif

#ifdef BRAGINSKII_VISCOSITY
#ifdef BRAGINSKII_VISCOSITY_CONSTANT_TIMESTEP
      // Manually set a constant time step. For testing purposes only.
      strcpy(tag[nt], "BragViscosityTimestep");
      addr[nt] = &All.BragViscosityTimestep;
      id[nt++] = PARAM_REAL;
#else
      // The Courant number. See equation 52 in the Braginskii code paper.
      strcpy(tag[nt], "BragViscosityCourant");
      addr[nt] = &All.BragViscosityCourant;
      id[nt++] = PARAM_REAL;
#endif
#ifndef BRAGINSKII_SPITZER
      // Use a constant viscosity coefficient. Mostly for testing purposes.
      strcpy(tag[nt], "BragViscosityCoefficient");
      addr[nt] = &All.BragViscosityCoefficient;
      id[nt++] = PARAM_REAL;
#endif
#ifdef BRAGINSKII_SPITZER
      // The Spitzer viscosity coefficient depends on density and temperature as ν_∥ ∝ T^(5/2)/ρ.
      // To avoid extremely small time steps, we set a maximum value for ν_∥.
      strcpy(tag[nt], "BragViscosityMaximumCoefficient");
      addr[nt] = &All.BragViscosityMaximumCoefficient;
      id[nt++] = PARAM_REAL;
#endif
#ifdef BRAGINSKII_RKL2_SUPER_TIME_STEPPING
      // See section 3.2 in the Braginskii code paper.
      All.BragViscosityRKL2Stages = 3;                // Initialize to the minimum number of stages allowed
      strcpy(tag[nt], "BragViscosityMaxRKL2Stages");  // The maximum number of stages allowed
      addr[nt] = &All.BragViscosityMaxRKL2Stages;
      id[nt++] = PARAM_INT;
#endif
#ifdef BRAGINSKII_VISCOSITY_SUBCYCLE
      All.BragViscositySubcycles = 1;
      strcpy(tag[nt], "BragViscosityMaxSubcycles");  // Maximum # of sub-cycling steps
      addr[nt] = &All.BragViscosityMaxSubcycles;
      id[nt++] = PARAM_INT;
#endif
#endif

#ifdef FLD
#ifdef FLD_TEST_BOUNDARY
      strcpy(tag[nt], "FLDFlux");
      addr[nt] = &All.fld_Flux;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "FLDGamma");
      addr[nt] = &All.fld_n_gamma;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "FLDDensity");
      addr[nt] = &All.fld_density;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "FLDEgySpec");
      addr[nt] = &All.fld_u;
      id[nt++] = PARAM_REAL;
#endif
#ifdef FLD_CONST_KAPPA
      strcpy(tag[nt], "Kappa_R");
      addr[nt] = &All.Kappa_R;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "Kappa_P");
      addr[nt] = &All.Kappa_P;
      id[nt++] = PARAM_REAL;
#endif
#ifdef FLD_MARSHAK
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef EXTERNALSHEARBOX
      strcpy(tag[nt], "ShearBoxSigma0");
      addr[nt] = &All.ShearBoxSigma0;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ShearBoxFg");
      addr[nt] = &All.ShearBoxFg;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "ShearBoxMu");
      addr[nt] = &All.ShearBoxMu;
      id[nt++] = PARAM_REAL;

#endif

#ifdef LOCAL_FEEDBACK
      strcpy(tag[nt], "LocalFeedbackSNEnergy");
      addr[nt] = &All.LocalFeedbackSNEnergy;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LocalFeedbackSNMassReturn");
      addr[nt] = &All.LocalFeedbackSNMassReturn;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LocalFeedbackSNRate");
      addr[nt] = &All.LocalFeedbackSNRate;
      id[nt++] = PARAM_REAL;

#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_SHOCK_ACCELERATION)
      strcpy(tag[nt], "LocalFeedbackCRInjectionFraction");
      addr[nt] = &All.LocalFeedbackCRInjectionFraction;
      id[nt++] = PARAM_REAL;
#endif

#ifdef LOCAL_FEEDBACK_PARTICLES
      strcpy(tag[nt], "LocalFeedbackSNTimeDelay");
      addr[nt] = &All.LocalFeedbackSNTimeDelay;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LocalFeedbackSNTimeSpread");
      addr[nt] = &All.LocalFeedbackSNTimeSpread;
      id[nt++] = PARAM_REAL;

#endif

#if !defined(EXTERNALSHEARBOX_KSRATE) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
      strcpy(tag[nt], "LocalFeedbackSFEff");
      addr[nt] = &All.LocalFeedbackSFEff;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "LocalFeedbackSFDenThresh");
      addr[nt] = &All.LocalFeedbackSFDenThresh;
      id[nt++] = PARAM_REAL;

#endif

#ifdef LOCAL_KINETIC
      strcpy(tag[nt], "LocalFeedbackKineticInjectionFraction");
      addr[nt] = &All.LocalFeedbackKineticInjectionFraction;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef LOCALIZED_SOFTENINGS
      strcpy(tag[nt], "LocalizedSofteningsParticleNumber");
      addr[nt] = &All.LocalizedSofteningsParticleNumber;
      id[nt++] = PARAM_INT;
#endif

#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION)
      strcpy(tag[nt], "OhmicDiffusionCoefficient");
      addr[nt] = &All.OhmicDiffusionCoefficient;
      id[nt++] = PARAM_REAL;
#endif
#ifdef AMBIPOLAR_DIFFUSION
      strcpy(tag[nt], "AmbipolarDiffusionCoefficient");
      addr[nt] = &All.AmbipolarDiffusionCoefficient;
      id[nt++] = PARAM_REAL;
#endif
#ifdef NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP
      strcpy(tag[nt], "NonidealMHDTimelimit");
      addr[nt] = &All.NonidealMHDTimelimit;
      id[nt++] = PARAM_REAL;
#endif
#endif

#ifdef SF_STELLAR_MASS_TO_GAS_MASS_RATIO
      strcpy(tag[nt], "StellarMassToGasMassRatio");
      addr[nt] = &All.StellarMassToGasMassRatio;
      id[nt++] = PARAM_REAL;
#endif

#ifdef AURIGA_MOVIE
      strcpy(tag[nt], "Auriga_Movie_CenterRadius");
      addr[nt] = &All.Auriga_Movie_CenterRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "Auriga_Movie_Directory");
      addr[nt] = All.Auriga_Movie_Directory;
      id[nt++] = PARAM_STRING;

      for(int ix = 0; ix < 3; ix++)
        for(int iy = 0; iy < 3; iy++)
          {
            sprintf(buf, "Auriga_Movie_Galaxy_Rotation%d%d", ix, iy);
            strcpy(tag[nt], buf);
            addr[nt] = &All.Auriga_Movie_Galaxy_Rotation[ix][iy];
            id[nt++] = PARAM_REAL;
          }

      strcpy(tag[nt], "Auriga_Movie_OutputListFilename");
      addr[nt] = All.Auriga_Movie_OutputListFilename;
      id[nt++] = PARAM_STRING;
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
      strcpy(tag[nt], "HighFreqStarsPath");
      addr[nt] = All.HighFreqStarsPath;
      id[nt++] = PARAM_STRING;
#endif

#ifdef HCOUTPUT
      strcpy(tag[nt], "HCOutput_RadialCut");
      addr[nt] = &All.HCOutput_RadialCut;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "HCOutput_CenterRadius");
      addr[nt] = &All.HCOutput_CenterRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "HCOutput_Directory");
      addr[nt] = All.HCOutput_Directory;
      id[nt++] = PARAM_STRING;

      strcpy(tag[nt], "HCOutput_OutputListFilename");
      addr[nt] = All.HCOutput_OutputListFilename;
      id[nt++] = PARAM_STRING;
#endif

#ifdef ONEDIMS_SPHERICAL
      strcpy(tag[nt], "CoreRadius");
      addr[nt] = &All.CoreRadius;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CoreMass");
      addr[nt] = &All.CoreMass;
      id[nt++] = PARAM_REAL;
#endif

#ifdef MODGRAV
      strcpy(tag[nt], "ModifiedGravityFile");
      addr[nt] = &All.ModifiedGravityFile;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "MaxAMRLevel");
      addr[nt] = &All.MaxAMRLevel;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "MinLevelTopLeaf");
      addr[nt] = &All.MinLevelTopLeaf;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "MaxLevelTopLeaf");
      addr[nt] = &All.MaxLevelTopLeaf;
      id[nt++] = PARAM_REAL;
      strcpy(tag[nt], "MaxLevelFullTree");
      addr[nt] = &All.MaxLevelFullTree;
      id[nt++] = PARAM_REAL;
#endif

#ifdef SGS_TURBULENCE
      strcpy(tag[nt], "MinimumSgsTSpecificEnergy");
      addr[nt] = &All.SgsTConst.MinimumSgsTSpecificEnergy;
      id[nt++] = PARAM_REAL;

#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
      strcpy(tag[nt], "EddyViscosityParameterCnu");
      addr[nt] = &All.SgsTConst.EddyViscosityParameterCnu;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "EddyViscosityParameterCepsilon");
      addr[nt] = &All.SgsTConst.EddyViscosityParameterCepsilon;
      id[nt++] = PARAM_REAL;
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */
#endif /* #ifdef SGS_TURBULENCE */

#ifdef SIMPLEX
#if SX_SOURCES == 10
      strcpy(tag[nt], "TestSrcFile");
      addr[nt] = All.sxTestSrcFile;
      id[nt++] = PARAM_STRING;
#endif

      strcpy(tag[nt], "UnitPhotons_per_s");
      addr[nt] = &All.UnitPhotons_per_s;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "MinNumPhotons");
      addr[nt] = &All.sxMinNumPhotons;
      id[nt++] = PARAM_REAL;
#endif

#ifdef AGB_WIND
      strcpy(tag[nt], "AGBWindDensity");
      addr[nt] = &All.AGBWindDensity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGBWindVelocity");
      addr[nt] = &All.AGBWindVelocity;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGBWindSpecificEnergy");
      addr[nt] = &All.AGBWindSpecificEnergy;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGBWindCenterX");
      addr[nt] = &All.AGBWindCenterX;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGBWindCenterY");
      addr[nt] = &All.AGBWindCenterY;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "AGBWindCenterZ");
      addr[nt] = &All.AGBWindCenterZ;
      id[nt++] = PARAM_REAL;
#endif

/* Stellar module */
#if defined(SOLAR_RADIATIVE_TRANSFER_DIFF) || defined(SOLAR_RADIATIVE_TRANSFER_EDD)
      strcpy(tag[nt], "VolumetricHeatingRate");
      addr[nt] = &All.VolumetricHeatingRate;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "VolumetricCoolingRate");
      addr[nt] = &All.VolumetricCoolingRate;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SurfaceRadiusInner");
      addr[nt] = &All.SurfaceRadiusInner;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "SurfaceRadiusOuter");
      addr[nt] = &All.SurfaceRadiusOuter;
      id[nt++] = PARAM_REAL;

      strcpy(tag[nt], "CoreRadius");
      addr[nt] = &All.CoreRadius;
      id[nt++] = PARAM_REAL;
#endif

      int errorflag = 0;
      FILE *fd      = fopen(fname, "r");
      if(fd)
        {
          sprintf(buf, "%s%s", fname, "-usedvalues");
          FILE *fdout = fopen(buf, "w");
          if(!fdout)
            {
              mpi_terminate("Error opening file '%s'", buf);
            }
          else
            {
              printf("Obtaining parameters from file '%s':\n\n", fname);
              while(!feof(fd))
                {
                  buf[0]       = 0;
                  char *status = fgets(buf, sizeof(buf), fd);
                  if(!status && ferror(fd))
                    terminate("error reading line from parameter file");

                  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                    continue;

                  if(buf1[0] == '%')
                    continue;

                  int j = -1;
                  for(int i = 0; i < nt; i++)
                    if(strcmp(buf1, tag[i]) == 0)
                      {
                        if(param_handled[i] == 0)
                          {
                            j                = i;
                            param_handled[i] = 1;
                            break;
                          }
                        else
                          {
                            j = -2;
                            break;
                          }
                      }

                  if(j >= 0)
                    {
                      switch(id[j])
                        {
                          case PARAM_REAL:
                            *((double *)addr[j]) = atof(buf2);
                            sprintf(buf3, "%%-%ds%%g\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((double *)addr[j]));
                            printf("        ");
                            printf(buf3, buf1, *((double *)addr[j]));
                            break;
                          case PARAM_STRING:
                            strcpy((char *)addr[j], buf2);
                            sprintf(buf3, "%%-%ds%%s\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, buf2);
                            printf("        ");
                            printf(buf3, buf1, buf2);
                            break;
                          case PARAM_INT:
                            *((int *)addr[j]) = atoi(buf2);
                            sprintf(buf3, "%%-%ds%%d\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((int *)addr[j]));
                            printf("        ");
                            printf(buf3, buf1, *((int *)addr[j]));
                            break;
                          default:
                            terminate("Incorrect parameter type");
                        }
                    }
                  else if(j == -2)
                    {
#ifdef ALLOWEXTRAPARAMS
                      warn("Tag '%s' ignored from file %s!", buf1, fname);
#else
                      mpi_printf("Error in file %s: Tag '%s' multiply defined.\n", fname, buf1);
                      errorflag = 1;
#endif
                    }
                  else
                    {
#ifdef ALLOWEXTRAPARAMS
                      warn("Tag '%s' ignored from file %s!", buf1, fname);
#else
                      mpi_printf("Error in file %s: Tag '%s' not allowed.\n", fname, buf1);
                      errorflag = 1;
#endif
                    }
                }
              fclose(fd);
              fclose(fdout);
              printf("\n");

              const int i = strlen(All.OutputDir);
              if(i > 0)
                if(All.OutputDir[i - 1] != '/')
                  strcat(All.OutputDir, "/");

              char *output_dir = All.OutputDir;
#ifdef SHOCK_FINDER_POST_PROCESSING
              output_dir = All.OutputDirShockFinder;
#endif
              mkdir(output_dir, MKDIR_MODE);
              file_path_sprintf(buf1, "%s-usedvalues", fname);
              file_path_sprintf(buf2, "%s/parameters-usedvalues", output_dir);
              sprintf(buf3, "cp %s %s", buf1, buf2);
              my_system(buf3);
            }
        }
      else
        {
          mpi_terminate("Parameter file %s not found.", fname);
        }

      for(int i = 0; i < nt; i++)
        if(param_handled[i] != 1)
          {
            mpi_printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
            errorflag = 1;
          }

      if(errorflag)
        mpi_terminate("Parameter file %s not valid.\n", fname);

      if(All.OutputListOn)
        {
          read_outputlist(All.OutputListFilename);
          mpi_printf("BEGRUN: Found %d times in output list.\n", All.OutputListLength);
        }
      else
        All.OutputListLength = 0;
    }

  All.NParameters = nt;

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef CHIMES
  MPI_Bcast(&ChimesGlobalVars, sizeof(struct globalVariables), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

#ifdef TOLERATE_WRITE_ERROR
  MPI_Bcast(AlternativeOutputDir, MAXLEN_PATH, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

#ifdef HOST_MEMORY_REPORTING
  check_maxmemsize_setting();
#endif

  mymalloc_init();

  Parameters      = (char(*)[MAXLEN_PARAM_TAG])mymalloc("Parameters", All.NParameters * MAXLEN_PARAM_TAG * sizeof(char));
  ParametersValue = (char(*)[MAXLEN_PARAM_VALUE])mymalloc("ParametersValue", All.NParameters * MAXLEN_PARAM_VALUE * sizeof(char));
  ParametersType  = (char *)mymalloc("ParamtersType", All.NParameters * sizeof(char));

  if(ThisTask == 0)
    {
      for(int i = 0; i < All.NParameters; i++)
        {
          snprintf(Parameters[i], MAXLEN_PARAM_TAG, "%s", tag[i]);
          ParametersType[i] = id[i];
          void *tmp         = ParametersValue[i];
          switch(id[i])
            {
              case PARAM_REAL:
                *((double *)tmp) = *((double *)addr[i]);
                break;
              case PARAM_STRING:
                snprintf((char *)tmp, MAXLEN_PARAM_VALUE, "%s", (char *)addr[i]);
                break;
              case PARAM_INT:
                tmp           = ParametersValue[i];
                *((int *)tmp) = *((int *)addr[i]);
                break;
            }
        }
    }

  MPI_Bcast(Parameters, sizeof(char) * All.NParameters * MAXLEN_PARAM_TAG, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ParametersValue, sizeof(char) * All.NParameters * MAXLEN_PARAM_VALUE, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ParametersType, sizeof(char) * All.NParameters, MPI_BYTE, 0, MPI_COMM_WORLD);

  for(int i = 0; i < NTYPES; i++)
    if(All.SofteningTypeOfPartType[i] >= NSOFTTYPES || All.SofteningTypeOfPartType[i] < 0)
      mpi_terminate("SofteningTypeOfPartType%d invalid (NSOFTTYPES=%d)", i, NSOFTTYPES);

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
        warn("NOTICE: Reducing requested NumFilesWrittenInParallel=%d to %d\n", All.NumFilesWrittenInParallel, NTask);
      All.NumFilesWrittenInParallel = NTask;
    }

  if(All.NumFilesWrittenInParallel == 0)
    {
      mpi_printf("NOTICE: All.NumFilesWrittenInParallel has been set to be equal to the number of processors\n");
      All.NumFilesWrittenInParallel = NTask;
    }

#ifdef SUBBOX_SNAPSHOTS
  if(All.SubboxNumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
        warn("NOTICE: Reducing requested SubboxNumFilesWrittenInParallel=%d to %d\n", All.SubboxNumFilesWrittenInParallel, NTask);
      All.SubboxNumFilesWrittenInParallel = NTask;
    }
  if(All.SubboxNumFilesWrittenInParallel > All.SubboxNumFilesPerSnapshot)
    {
      if(ThisTask == 0)
        warn("Reducing requested SubboxNumFilesWrittenInParallel=%d to be equal SubboxNumFilesPerSnapshot=%d",
             All.SubboxNumFilesWrittenInParallel, All.SubboxNumFilesPerSnapshot);
      All.SubboxNumFilesWrittenInParallel = All.SubboxNumFilesPerSnapshot;
    }
#endif

#ifndef GRAVITY_NOT_PERIODIC
  if(!All.PeriodicBoundariesOn)
    mpi_terminate(
        "Code was compiled with gravity periodic boundary conditions switched on.\n"
        "You must set `PeriodicBoundariesOn=1', or recompile the code.");
#else
  if(All.PeriodicBoundariesOn)
    mpi_terminate(
        "Code was compiled with gravity periodic boundary conditions switched off.\n"
        "You must set `PeriodicBoundariesOn=0', or recompile the code.");
#endif

#ifdef COOLING
  if(!All.CoolingOn)
    mpi_terminate(
        "Code was compiled with cooling switched on.\n"
        "You must set `CoolingOn=1', or recompile the code.");
#else
  if(All.CoolingOn)
    mpi_terminate(
        "Code was compiled with cooling switched off.\n"
        "You must set `CoolingOn=0', or recompile the code.");
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    mpi_terminate("The specified timestep criterion is not valid.");

#if NTYPES < 6
#error "NTYPES < 6 is not allowed."
#endif
#if NTYPES > 15
#error "NTYPES > 15 is not supported yet."
#endif
#if NTYPES > 8
  if(All.ICFormat == SNAP_FORMAT_GADGET || All.ICFormat == SNAP_FORMAT_GADGET_VARIANT)
    mpi_terminate("NTYPES > 8 is not allowed with ICFormat=%d, since the header block is limited to %zu bytes.", All.ICFormat,
                  sizeof(struct io_header));
#endif
#ifdef USE_SFR
  if(!All.StarformationOn)
    mpi_terminate(
        "Code was compiled with star formation switched on.\n"
        "You must set `StarformationOn=1', or recompile the code.");
#ifndef LOCAL_FEEDBACK
  if(!All.CoolingOn)
    mpi_terminate(
        "You try to use the code with star formation enabled,\nbut you did not switch on cooling.\n"
        "This mode is not supported.");
#endif
#else
  if(All.StarformationOn)
    mpi_terminate(
        "Code was compiled with star formation switched off.\n"
        "You must set `StarformationOn=0', or recompile the code.");
#endif

#if defined(TRACER_PARTICLE) && ((TRACER_PART_TMAX_TIME) || (TRACER_PART_TMAX_RHO)) && !(TRACER_PART_TMAX)
#error "Code was compiled with TRACER_PART_TMAX_TIME or TRACER_PART_TMAX_RHO.\nYou must compile with TRACER_PART_TMAX as well."
#endif

#if defined(TRACER_PART_NUM_FLUID_QUANTITIES) && !defined(TRACER_PARTICLE)
#error "Code was compiled with TRACER_PART_NUM_FLUID_QUANTITIES but without TRACER_PARTICLE.\nThis is not allowed."
#endif

#if defined(TRACER_MC_NUM_FLUID_QUANTITIES) && !defined(TRACER_MC)
#error "Code was compiled with TRACER_MC_NUM_FLUID_QUANTITIES but without TRACER_MC.\nThis is not allowed."
#endif

#ifdef TRACER_MC
#ifdef REFINEMENT_MERGE_PAIRS
#error \
    "Code was compiled with TRACER_MC and REFINEMENT_MERGE_PAIRS together.\n" "This is not supported yet."
#endif
#if((TRACER_MC_TMAX_TIME) || (TRACER_MC_TMAX_RHO)) && !(TRACER_MC_TMAX)
#error \
    "Code was compiled with TRACER_MC_TMAX_TIME or TRACER_MC_TMAX_RHO.\n" "You must compile with TRACER_MC_TMAX as well."
#endif
#if TRACER_MC_LAST_STAR_TIME && !(defined(GFM_STELLAR_EVOLUTION) || defined(GFM_WINDS))
#error \
    "Code was compiled with TRACER_MC and TRACER_MC_LAST_STAR_TIME but without GFM_STELLAR_EVOLUTION or GFM_WINDS.\nThere is no point doing that."
#endif
#if TRACER_MC_WIND_COUNTER && !defined(GFM_WINDS)
#error "Code was compiled with TRACER_MC_WIND_COUNTER but without GFM_WINDS.\nThis is not allowed."
#endif
#if defined(ADD_GROUP_PROPERTIES) || defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
#error "TRACER_MC compiled with ADD_GROUP_PROPERTIES or RECOMPUTE_POTENTIAL_IN_SNAPSHOT. More thought required in I/O."
#endif
#endif /* TRACER_MC */

#if defined(ENFORCE_JEANS_STABILITY_OF_CELLS) &&                                                              \
    (defined(ISOTHERM_EQS) || defined(LOCALLY_ISOTHERM_DISK) || defined(TGCHEM) || defined(EOS_DEGENERATE) || \
     (defined(USE_SFR) && !defined(SMUGGLE_SFR)) || defined(EOS_OPAL))
  if(ThisTask == 0)
    warn("Code was compiled with ENFORCE_JEANS_STABILITY_OF_CELLS together with another EOS. Please make sure you really want this.");
#endif

#if(defined(GFM_CONST_IMF) || defined(GFM_VARIABLE_IMF)) && !defined(GFM_STELLAR_EVOLUTION)
#error "Code was compiled with GFM_*_IMF, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed."
#endif

#if defined(GFM_DUST) && !defined(GFM_STELLAR_EVOLUTION)
#error "Code was compiled with GFM_DUST, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed."
#endif

#if defined(DL_GRAIN_BINS) && !defined(GFM_STELLAR_EVOLUTION)
#error "Code was compiled with DL_GRAIN_BINS, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed."
#endif

#if defined(GFM_CHEMTAGS) && !defined(GFM_STELLAR_EVOLUTION)
#error "Code was compiled with GFM_CHEMTAGS, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed."
#endif

#if defined(GFM_SPLITFE) && !defined(GFM_STELLAR_EVOLUTION)
#error "Code was compiled with GFM_SPLITFE, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed."
#endif

#if defined(GFM_SPLITFE) && !defined(GFM_CHEMTAGS)
#error "Code was compiled with GFM_SPLITFE, but not with GFM_CHEMTAGS.\nThis is not allowed."
#endif

#if defined(GFM_RPROCESS) && !defined(GFM_SPLITFE)
#error "Code was compiled with GFM_RPROCESS, but not with GFM_SPLITFE.\nThis is not allowed."
#endif

#if defined(GFM_RPROCESS) && !defined(GFM_CHEMTAGS)
#error "Code was compiled with GFM_RPROCESS, but not with GFM_CHEMTAGS.\nThis is not allowed."
#endif

#if defined(GFM_DUST) && defined(GFM_SPLITFE)
#error "Code was compiled with GFM_SPLITFE, and GFM_DUST.\nThis is not supported yet."
#endif

#if defined(GFM_CONST_IMF) && defined(GFM_VARIABLE_IMF)
#error "Code was compiled with GFM_CONST_IMF and GFM_VARIABLE_IMF together.\nThis is not allowed."
#endif

#if defined(GFM_WINDS_VARIABLE) && !defined(FOF)
#error "Code was compiled with GFM_WINDS_VARIABLE, but not with FOF.\nThis is not allowed."
#endif

#if(defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 1)) && !defined(FOF)
#error "Code was compiled with GFM_BIPOLAR_WINDS == 1, but not with FOF.\nThis is not allowed."
#endif

#if defined(GFM_WINDS_SAVE_PARTTYPE) && !defined(GFM_WINDS)
#error "GFM_WINDS_SAVE_PARTTYPE requires GFM_WINDS."
#endif

#if defined(GFM_DISCRETE_ENRICHMENT) && (!defined(GFM) || !defined(GFM_STELLAR_EVOLUTION))
#error "GFM_DISCRETE_ENRICHMENT requires both GFM and GFM_STELLAR_EVOLUTION."
#endif

#if defined(BH_BUBBLES) && !defined(BLACK_HOLES)
#error "Code was compiled with BH_BUBBLES, but not with BLACK_HOLES.\nThis is not allowed."
#endif

#if defined(REPOSITION_ON_POTMIN) && !defined(EVALPOTENTIAL)
#error "REPOSITION_ON_POTMIN requires EVALPOTENTIAL."
#endif

#if defined(BH_BASED_CGM_ZOOM) && (!defined(BLACK_HOLES) || !defined(REPOSITION_ON_POTMIN) || !defined(BH_DO_NOT_PREVENT_MERGERS))
#error "BH_BASED_CGM_ZOOM: Expect also BLACK_HOLES, REPOSITION_ON_POTMIN, and BH_DO_NOT_PREVENT_MERGERS."
#endif

#if defined(VORONOI_FREQUENT_IMAGES) && defined(IMAGES_FOREACHSNAPSHOT)
#error "VORONOI_FREQUENT_IMAGES and IMAGES_FOREACHSNAPSHOT can't be used together."
#endif

#if defined(METALS) && !defined(USE_SFR) && !defined(ADDBACKGROUNDGRID)
#error "Code was compiled with METALS, but not with USE_SFR.\nThis is not allowed."
#endif

#if defined(MIN_METALLICITY_ON_STARTUP) && !defined(METALS)
#error "Code was compiled with MIN_METALLICITY_ON_STARTUP but not METALS.\nThis is not allowed."
#endif

#if defined(TIMEDEPDE) && !defined(DARKENERGY)
#error "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\nThis is not allowed."
#endif

#if defined(OTVET) && defined(GFM_COOLING_METAL)
#error "Code was compiled with OTVET and GFM_COOLING_METAL together.\nThis is not supported yet."
#endif

#ifdef OTVET_FIXTIMESTEP
  if(All.MaxSizeTimestep != All.MinSizeTimestep)
    mpi_terminate("OTVET_FIXTIMESTEP requires All.MaxSizeTimestep == All.MinSizeTimestep. Fix it in parameterfile. Stopping...");
#endif

#if defined(OTVET_SCATTER_SOURCE) && defined(EDDINGTON_TENSOR_STARS) && defined(GFM)
#error "Code was compiled with OTVET_SCATTER_SOURCE for stars and GFM together.\nThis is not supported yet."
#endif

#if defined(SMUGGLE_RADIATION_FEEDBACK) && !defined(GFM_STELLAR_EVOLUTION)
#error \
    "Code was compiled with SMUGGLE_RADIATION_FEEDBACK that requires SMUGGLE_STAR_FEEDBACK and GFM_STELLAR_EVOLUTION active.\nPlease activate both."
#endif

#ifdef MODIFIED_EOS
  check_modified_eos_parameters();
#endif

/*
#ifdef REFINEMENT
  # If we restarted from a snapshot, and are not adiabatic, the initialisation of ReferenceGasPartMass will change
  if(RestartFlag == RESTART_SNAPSHOT && All.StarformationOn && All.ReferenceGasPartMass == 0)
    {
      terminate(
          "Restarting from a snapshot will recompute ReferenceGasPartMass but the result will be different if gas has been converted "
          "to stars already."
          "Set ReferenceGasPartMass in param.txt explicitely to the desired value. If you are sure you want to recompute "
          "ReferenceGasPartMass comment out this terminate statement.");
    }
#endif 
*/

#if defined(BLACK_HOLES) && (!defined(BH_RELATIVE_NGB_DEVIATION))
  if(All.DesNumNgbBlackHole == All.MaxNumNgbDeviation)
    mpi_terminate(
        "Code was compiled without BH_RELATIVE_NGB_DEVIATIONS and hence requires All.DesNumNgbBlackHole > All.MaxNumNgbDeviation.\n"
        "Fix choices in parameterfile or the code will crash. Stopping.");
#endif

#ifdef GFM_SINGLE_CELL_INJECTION
  All.DesNumNgbEnrichment          = 1;
  All.MaxNumNgbDeviationEnrichment = 0.1;
#endif
}

/*! \brief This function checks the consistency of the input parameters.
 *
 *  If you encounter some possible misuse and a corresponding error message
 *  that is hard to interpret, a check should be placed in this function with
 *  a terminate statement and a clear explanation why this does not work.
 *
 *  \return void
 */
void check_parameters(void)
{
  /* check whether time max is larger than max timestep */
  if(All.TimeMax - All.TimeBegin <= All.MaxSizeTimestep)
    mpi_terminate(
        "check_parameters: Your total runtime is smaller or equal than the maximum allowed timestep!\n"
        "Choose an appropriate value for MaxSizeTimestep < TimeMax - TimeBegin! (TimeBegin = %g, TimeMax = %g, MaxSizeTimestep = %g)",
        All.TimeBegin, All.TimeMax, All.MaxSizeTimestep);
  /* check if number of snapshot files to write is larger than number of tasks */
  if(NTask < All.NumFilesPerSnapshot)
    {
      warn(
          "Number of processors (%d) must be at least as large as All.NumFilesPerSnapshot (%d)! Reducing All.NumFilesPerSnapshot "
          "accordingly.\n",
          NTask, All.NumFilesPerSnapshot);
      All.NumFilesPerSnapshot = NTask;
    }
#ifdef SUBFIND
  /* check that DesLinkNgb is not larger than FOF_GROUP_MIN_LEN (otherwise the SUBFIND
   * neighbor search won’t work for small groups) */
  if(All.DesLinkNgb > FOF_GROUP_MIN_LEN)
    terminate("DesLinkNgb (%d) should not be larger than FOF_GROUP_MIN_LEN (%d)!", All.DesLinkNgb, FOF_GROUP_MIN_LEN);
#endif
}

/*! \brief This function reads a table with a list of desired output times.
 *
 *  The table does not have to be ordered in any way, but may not contain more
 *  than MAXLEN_OUTPUTLIST entries.
 *
 *  \param[in] fname The file name of the outputlist.
 *
 *  \return void
 */
void read_outputlist(const char *const fname)
{
  FILE *fd = fopen(fname, "r");
  if(!fd)
    terminate("Can't read output list in file '%s'", fname);
  All.OutputListLength = 0;
  while(1)
    {
      char buf[512];
      if(fgets(buf, sizeof(buf), fd) != buf)
        break;
      int flag;
      int count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);
      if(count == 1)
        flag = 1;
      if(count == 1 || count == 2)
        {
          if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
            terminate("too many entries in output list. You should increase MAXLEN_OUTPUTLIST=%d.", MAXLEN_OUTPUTLIST);
          All.OutputListFlag[All.OutputListLength] = flag;
          All.OutputListLength++;
        }
    }
  fclose(fd);
}
