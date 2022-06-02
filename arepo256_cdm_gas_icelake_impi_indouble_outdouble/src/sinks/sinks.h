/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sinks/sinks.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles structures and constant definitions
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifdef SINKS

#define SKD_INIT_MAX_NUM_ACC 10
#define MAX_MASS_FAC 1.0
#define ACCRETION_SAFETY_FACTOR 0.02
#define MIN_TARGET_MASS_FACTOR_FOR_ACC 0.00001

extern struct ACC_struct
{
  int Task;
  int Index;
  int Active;
  double Pos[3];
  double Vel[3];
  double Mass;
#ifdef SGCHEM
  double Utherm;
  double Accel[3];
#endif
  int rcells;
} * ACC;

struct SINK_struct
{
  double Pos[3];
  double Mass;
};

extern struct SKD_struct
{
  int Task;
  int Index;
  int Flag;
  int NumAcc;
  int MaxNumAcc;
  int NumSinks;
  int TotNumSinks;
  int TimeBin;
  MyIDType ID;
  double NHThresh;
  double AccRad;
  double AccRad2;
  double DistFac;
  double NHFac;
  double NHMax;
  double Pos[3];
  double Mass;
  double CoMPos[3];
  double CoMVel[3];
  double Temp;

  // struct ACC_struct *ACC;
  struct SINK_struct *SINK;

} SKD;

extern struct SinksAux_struct
{
  int SinksAuxID;
  signed char MinTimeBin;
  signed char MaxTimeBin;
  double dMass;
  double MassNorm;
} * SinksAux;

extern int NSinks, TotNSinks;

#endif
