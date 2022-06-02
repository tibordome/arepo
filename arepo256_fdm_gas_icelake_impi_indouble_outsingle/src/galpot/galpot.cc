/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/barpot/barpot.cc
 * \date        09/2016
 * \author      Mattia Sormani, Robin Tress 
 * \brief       This code calculates the contribution of different galactic components to the potential.
 *
 * \details     
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "common.h"
#include "potential.h"
#include "../allvars.h"

ExpDiskPotential* Diskpot1;
ExpDiskPotential* Diskpot2;
DPLPotential* Halopot;
McMillanBulgePotential* Bulgepot;
SpiralArmPotential* Spiralpot;

#ifdef __cplusplus
extern "C" 
{
#endif
  void galpot_init(void)
  {
    double G = All.G;   

  /*Initializing Bulge potential*/
    double rminbulge = 0., rmaxbulge = 0.;
    int nrbulge = 2000, npolybulge = 6, ngaussbulge = 8;
    double rho0bulge = 9.84e7, a0bulge=0.75, acutbulge = 21.0;
    double alfbulge = 1.8, qbulge = 0.5;

    McMillanBulgeDensity   rhobulge(rho0bulge, alfbulge, a0bulge, acutbulge, qbulge, rminbulge, rmaxbulge);
    Bulgepot = new McMillanBulgePotential(G, rhobulge, nrbulge, npolybulge, ngaussbulge);

  /*Initializing exponential disk potential 
   *TODO: multipole expansion not ideal for disks. Consider switching to GalPot by Paul McMillan and Walter Dehnen */
    /*Thick stellar disk*/
    double rmind1 = 0., rmaxd1 = 0.;
    int nrd1 = 2000, npolyd1 = 20, ngaussd1 = 8;
    double Sigma0d1 = 1.83e6, Rd1 = 30.2, zd1 = 9.;

    ExpDiskDensity   rhodisk1(Sigma0d1, zd1, Rd1, rmind1, rmaxd1);
    Diskpot1 = new ExpDiskPotential(G, rhodisk1, nrd1, npolyd1, ngaussd1);

    /*Thin stellar disk*/
    double rmind2 = 0., rmaxd2 = 0.;
    int nrd2 = 2000, npolyd2 = 20, ngaussd2 = 8;
    double Sigma0d2 = 8.96e6, Rd2 = 25.0, zd2 = 3.0;

    ExpDiskDensity   rhodisk2(Sigma0d2, zd2, Rd2, rmind2, rmaxd2);
    Diskpot2 = new ExpDiskPotential(G, rhodisk2, nrd2, npolyd2, ngaussd2);

  /*Initializing Halo potential*/
    double rminh = 0., rmaxh = 0.;
    int nrh = 2000, npolyh = 6, ngaussh = 8;
    double rho0h = 8.54e3, ah = 196., alphah = 1.0, betah = 3.0;

    DPLDensity   rhohalo(rho0h, ah, alphah, betah, rminh, rmaxh);
    Halopot = new DPLPotential(G, rhohalo, nrh, npolyh, ngaussh);

  /*Initializing spiral perturbation*/
    double SNarms = 4.0; //numbers of arms 
    double Sphir = 6.3e-16*All.UnitTime_in_s; //angular velocity of the spiral arm pattern  
    double Salpha = 0.261799; //pitch angle 

    double SHz = 5.56916e20/All.UnitLength_in_cm; //scale-height of the stellar arm perturbation
    double Sp0 = 2.12889e-24*pow(All.UnitLength_in_cm,3)/All.UnitMass_in_g; //midplane arm density at radius Sr0
    double Sr0 = 2.47518e22/All.UnitLength_in_cm; //radius at which the midplane arm density is Sp0
    double Srs = 2.16578e22/All.UnitLength_in_cm; //radial scale length of the drop-off in density amplitude of the arms
    double St0 = 3.153e+15/All.UnitTime_in_s; //initial time

    Spiralpot = new SpiralArmPotential(Sr0, St0, SNarms, Sphir, Salpha, SHz, Sp0, Srs, G);

    return;
  }

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" 
{
#endif
  void galpot_dPhidx(double x, double y, double z, double time, double *dPhidx)
  {
    Vec3 X(x,y,z);

    dPhidx[0] = 0.; 
    dPhidx[1] = 0.; 
    dPhidx[2] = 0.;

  /*Bulge potential*/
    Vec3 dBulgePhidx = Bulgepot->dPhidx(X);
    
    dPhidx[0] += dBulgePhidx[0];
    dPhidx[1] += dBulgePhidx[1];
    dPhidx[2] += dBulgePhidx[2];

  /*Disk potential*/ 
    /*Thick disk*/
    Vec3 dDisk1Phidx = Diskpot1->dPhidx(X);

    dPhidx[0] += dDisk1Phidx[0];
    dPhidx[1] += dDisk1Phidx[1];
    dPhidx[2] += dDisk1Phidx[2];

    /*Thin disk*/
    Vec3 dDisk2Phidx = Diskpot2->dPhidx(X);

    dPhidx[0] += dDisk2Phidx[0];
    dPhidx[1] += dDisk2Phidx[1];
    dPhidx[2] += dDisk2Phidx[2];

  /*Halo potential*/
    Vec3 dDPLPhidx = Halopot->dPhidx(X);

    dPhidx[0] += dDPLPhidx[0];
    dPhidx[1] += dDPLPhidx[1];
    dPhidx[2] += dDPLPhidx[2];

  /*Spiral arm perturbation*/
    Vec3 dSpiralPhidx = Spiralpot->dPhidx(X, time);

    dPhidx[0] += dSpiralPhidx[0];
    dPhidx[1] += dSpiralPhidx[1];
    dPhidx[2] += dSpiralPhidx[2];

    /*    printf("TestPOT Bulge %f|%f|%f, Disc1 %g|%g|%g, Disc2 %g|%g|%g, Halo %g|%g|%g, Spiral %g|%g|%g \n",dBulgePhidx[0],dBulgePhidx[1],dBulgePhidx[2],dDisk1Phidx[0],dDisk1Phidx[1],dDisk1Phidx[2],dDisk2Phidx[0],dDisk2Phidx[1],dDisk2Phidx[2],dDPLPhidx[0],dDPLPhidx[1],dDPLPhidx[2],dSpiralPhidx[0],dSpiralPhidx[1],dSpiralPhidx[2] );*/

    return;
  }
#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
extern "C"
{
#endif
/*! \brief Returns the mass of stars in a given gas particle defined by the stellar density distribution 
 *         defined by galpot. This is used for instance by the SN routines to have a random SN 
 *         distribution with a profile following the stellar profile. These then would represent SNIa.
 *
 *         Make sure that this piece of code is consistent with the definition of the stellar profile 
 *         defined in the functions above. Follow the example witch is outcommented below to do this.
 */
  double galpot_find_stellar_mass_of_particle(double x, double y, double z, double time, double volume)
  {
//
//    Vec3 X(x,y,z);
//
//  /*Bar potential*/
//    double omega_b = -4.; //internal units
//    double t_grow = 1.5; //internal units
//
//    double costheta = cos(omega_b * time);
//    double sintheta = sin(omega_b * time);
//    Vec3 Xbar((x*costheta + y*sintheta), (-x*sintheta + y*costheta), z);
//
//    double bar_density;
//    if (time > t_grow)
//      bar_density = Barpot->rho(Xbar);
//    else
//      bar_density = (time / t_grow) * Barpot->rho(Xbar) + (1. - (time / t_grow)) * Barpot_axi->rho(Xbar);
//
//  /*Bulge potential*/
//    double bulge_density = Bulgepot->rho(X);
//
//  /*Disk potential*/
//    /*Thick disc*/
//    double disc1_density = Diskpot1->rho(X);
//
//    /*Thin disc*/
//    double disc2_density = Diskpot2->rho(X);
//
//    return volume * (bar_density + bulge_density + disc1_density + disc2_density);
    terminate("CAREFUL: the function galpot_find_stellar_mass_of_particle() is ill defined. Modify it so that it returnes the actual stellar mass consistent with your definition in galpot_init()");
    return 1;
  }
#ifdef __cplusplus
}
#endif

