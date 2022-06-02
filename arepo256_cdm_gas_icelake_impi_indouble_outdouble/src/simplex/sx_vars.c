/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \date        24/03/2021
 * \author      Simon May
 * \brief       Definitions of variables and structures in SPRAI
 * \details
 */

/* see sx_proto.h for documentation on the variables defined here */
#include "sx_proto.h"

struct sxRunData_struct *sxRunData;

struct sxSource_struct *sxSources;

struct sxCell_struct *sxCell;

int sxNumSPD;

struct sxABPP_struct *sxAPP, *sxBPP, *sxRealAPP, *sxRealBPP;
int sxMaxNumABPP;
int sxSwapABPP;
int sxNumReallocABPP;

int sxNumAPP, sxNumBPP;
int *sxSPDtoAPP, *sxSPDtoBPP;
int *sxRealSPDtoAPP, *sxRealSPDtoBPP;

int sxNumTotEDP;
int sxNumTotEDPD;

struct sxQPP_struct *sxQPP;
int sxMaxNumQPP;
int *sxEDPDtoQPP;
int sxNumReallocQPP;

struct sxEPP_struct *sxExportPP;
int sxNumExportPP, sxNumImportPP;
int *sxNumEDP, *sxOffEDP;
int *sxDPtoEDP;

int SXPID;
int SXID;

gsl_rng *sxRand;
