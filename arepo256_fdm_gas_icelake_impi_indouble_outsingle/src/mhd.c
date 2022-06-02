#include "allvars.h"
#include "proto.h"

#ifdef MHD

#ifdef MHD_DEDNER
static void decay_psi(void);
static void dedner_source_terms(void);
#endif

static void do_mhd_source_terms(void);
#ifdef MHD_POWELL_SPLIT
static void do_mhd_powell_source_terms(void);
#endif

/*! \brief First half of the MHD source terms.
 *
 *  Before hydrodynamics timestep.
 *
 *  \return void
 */
void do_mhd_source_terms_first_half(void)
{
#ifdef MHD_DEDNER
  decay_psi();
#endif
#ifdef MHD_POWELL_SPLIT
  do_mhd_powell_source_terms();
#endif
  do_mhd_source_terms();
  update_primitive_variables();
}

/*! \brief Second half of the MHD source terms.
 *
 *  After hydrodynamics timestep.
 *
 *  \return void
 */
void do_mhd_source_terms_second_half(void)
{
  do_mhd_source_terms();
#ifdef MHD_POWELL_SPLIT
  do_mhd_powell_source_terms();
#endif
#ifdef MHD_DEDNER
  decay_psi();
#endif
  update_primitive_variables();
}

void do_mhd_source_terms(void)
{
  TIMER_START(CPU_MHD);

  if(All.ComovingIntegrationOn)
    {
      double atime    = All.Time;
      double hubble_a = hubble_function(atime);

      int idx, i;
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval /
                           hubble_a; /* half the timestep of the cell */
          SphP[i].Energy += dt_cell * 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) *
                            SphP[i].Volume * atime * hubble_a;
        }
    }

  TIMER_STOP(CPU_MHD);
}

#ifdef MHD_POWELL_SPLIT
void do_mhd_powell_source_terms(void)
{
  TIMER_START(CPU_MHD);

  if(All.ComovingIntegrationOn)
    terminate("do_mhd_powell_source_terms still lacks the cosmological factors");

  set_cosmo_factors_for_current_time();
  calculate_gradients();  // this updates divB

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval /
                       All.cf_hubble_a; /* half the timestep of the cell */
      dt_cell *= SphP[i].DivB;

      for(int k = 0; k < 3; k++)
        {
          double dMomentum = -SphP[i].B[k] * dt_cell;
          SphP[i].Momentum[k] += dMomentum;
          All.Powell_Momentum[k] += dMomentum;

          double dEnergy = -P[i].Vel[k] * SphP[i].B[k] * dt_cell;
          SphP[i].Energy += dEnergy;
          All.Powell_Energy += dEnergy;

          double dB = -P[i].Vel[k] * dt_cell;
          SphP[i].BConserved[k] += dB;
        }
    }

  TIMER_STOP(CPU_MHD);
}
#endif

#ifdef MHD_DEDNER
/* divB cleening is implemented as described in Dedner et. al 2002 */
void decay_psi(void)
{
  TIMER_START(CPU_MHD_DEDNER);

  double facmin = MAX_DOUBLE_NUMBER;
  double facmax = 0;

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef MHD_DEDNER_VARIABLE_SPEED
      double cf   = get_dedner_speed(i);
      double vpot = cf > 0. ? 2. * SphP[i].Psi / cf : 0.;
      double fac  = 0.5 * fmax(sqrt(cf * cf + vpot * vpot), 0.01 * All.DednerSpeed) / get_cell_radius(i);

      if(fac > facmax)
        facmax = fac;
      if(fac < facmin)
        facmin = fac;
#else
      double fac = 0.5 * All.DednerSpeed / get_cell_radius(i);
#endif
      double dt_cell =
          0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;

      SphP[i].PsiConserved *= exp(-dt_cell * fac);
      SphP[i].Psi *= exp(-dt_cell * fac);
    }

  double facminall, facmaxall;
  MPI_Allreduce(&facmax, &facmaxall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&facmin, &facminall, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  mpi_printf("MHD_DEDNER: decay timescale min=%g, max=%g\n", facminall, facmaxall);

  TIMER_STOP(CPU_MHD_DEDNER);
}

double get_dedner_speed(int i)
{
#ifdef MHD_DEDNER_VARIABLE_SPEED
  return SphP[i].DednerSpeed;
#else
  return All.DednerSpeed;
#endif
}
#endif

#ifdef MHD_DEDNER_VARIABLE_SPEED
static void compute_dedner_speeds_local(void)
{
  double DednerSpeedMin = MAX_DOUBLE_NUMBER;
  double DednerSpeedMax = 0;

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].DednerSpeed = get_sound_speed(i);
#ifdef VORONOI_STATIC_MESH
      SphP[i].DednerSpeed += sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) / All.cf_atime;
#endif

      if(SphP[i].DednerSpeed < DednerSpeedMin)
        DednerSpeedMin = SphP[i].DednerSpeed;
      if(SphP[i].DednerSpeed > DednerSpeedMax)
        DednerSpeedMax = SphP[i].DednerSpeed;
    }

  double DednerSpeedMaxAll, DednerSpeedMinAll;
  MPI_Allreduce(&DednerSpeedMax, &DednerSpeedMaxAll, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&DednerSpeedMin, &DednerSpeedMinAll, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  mpi_printf("MHD_DEDNER: local speed min=%g, max=%g\n", DednerSpeedMinAll, DednerSpeedMaxAll);
}
#endif

#ifdef MHD_DEDNER
void compute_dedner_speed(void)
{
  TIMER_START(CPU_MHD_DEDNER);

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    {
      double cmax = 0;

      for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          int i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          double csnd = get_sound_speed(i);

#ifdef VORONOI_STATIC_MESH
          csnd += sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]) / All.cf_atime;
#endif

          if(csnd > cmax)
            {
              cmax = csnd;
            }
        }

      double cmaxall;
      MPI_Allreduce(&cmax, &cmaxall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      All.DednerSpeed = cmaxall;

      if(All.DednerSpeed <= 0.)
        All.DednerSpeed = 1.0;

      mpi_printf("MHD_DEDNER: global speed=%g\n", All.DednerSpeed);
    }

#ifdef MHD_DEDNER_VARIABLE_SPEED
  compute_dedner_speeds_local();
#endif

  TIMER_STOP(CPU_MHD_DEDNER);
}
#endif

#endif
