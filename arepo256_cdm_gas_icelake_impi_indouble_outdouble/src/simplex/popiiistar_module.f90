! /*!
!  * \copyright   This file is part of the AREPO code developed by Volker Springel.
!  * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
!  * \copyright   and contributing authors.
!  *  
!  * \file        src/simplex/popiiistar_module.f90
!  * \date        8/2017
!  * \author      Sam Geen (STG)
!  * \brief       Read + interpolate over PopIII stellar properties 
!  * \details     
!  * 
!  * 
!  * \par Major modifications and contributions:
!  * 
!  * - 23.08.2017 Created (STG)
!  */

! /********************************************
!  *         Pop III Properties Tables        *
!  ********************************************/

MODULE popiiistar_module

! Gets values from popiii star tables
! Sam Geen, August 2017

! Use lookup table module for accretion, mass values
use lookup_table_module

implicit none

public

integer,parameter::dp=kind(1.0D0) ! real*8  

! Tables for metal cooling with and without flux
type(lookup_table) :: radius_table
type(lookup_table) :: temperature_table

! Unit conversions to cgs
real(dp),parameter::solar_mass_cgs=1.9891d33
real(dp),parameter::solar_mass_per_year_cgs=6.30321217d25
real(dp),parameter::solar_radius_cgs=6.957d10

END MODULE

!************************************************************************
! Sets up the tables, and then clears them (e.g. on program exit)


SUBROUTINE setup_popiii_tables
  use popiiistar_module
  character(len=128)              ::dir,filename
  !dir = "./"
  filename = TRIM('./popiii_grid_radius.dat')
  !write(*,*) "READING FILE", filename
  call setup_table(radius_table, filename)
  filename = TRIM('./popiii_grid_T.dat')
  !write(*,*) "READING FILE", filename
  call setup_table(temperature_table, filename)
  !write(*,*) "SX: Set up Pop III stellar evolution tables"
END SUBROUTINE setup_popiii_tables



SUBROUTINE popiii_lookup(accretion,mass,radius,temperature)
  ! Look up (radius, temperature) of a Pop III star based on its accretion rate and mass
  ! RETURNS
  ! 

  use popiiistar_module

  !integer,parameter::dp=kind(1.0D0) ! real*8  
  real(dp),intent(in)::accretion
  real(dp),intent(in)::mass
  real(dp),intent(out)::radius
  real(dp),intent(out)::temperature
  real(dp)::asolar,msolar
  
  ! Convert input values (already in cgs, but accretion is in log units)
  if (accretion.eq.0d0) then
     asolar=-20d0
  else
     asolar = log10(accretion)
  endif
  msolar = mass
  ! Read radius (output = solar radii)
  call find_value2(radius_table,     asolar,msolar,radius)
  ! Read temperature (output = log10(K))
  call find_value2(temperature_table,asolar,msolar,temperature)
  ! Convert output to linear cgs
  radius = 10d0**radius * solar_radius_cgs
  temperature = 10d0**temperature


END SUBROUTINE popiii_lookup
