/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/tgset/tgset.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Special settings for primordial runs
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#define TGSET_JEANS_TABLE_SIZE 100

extern struct TGD_struct
{
  // Parameters
  int MaxTimeBinDiff;
  double SnapFreeFallTimeFac;
  double NHInit;
  double NHTerm;
  double StreamVel;
  double JeansTemp;
  double JeansDensLow;
  double JeansDensHigh;
  double JeansNumberLow;
  double JeansNumberHigh;

  // Other
  int MinTimeBin;
  int WriteLogs;
  int NHMaxTask;
  int NHMaxIdx;
  int DoJeansRef;
  double NHMax;
  double JeansLogNHTableRange;
  double JeansLogNHTable[TGSET_JEANS_TABLE_SIZE];
  double TargetJeansNumberTable[TGSET_JEANS_TABLE_SIZE];

} TGD;
