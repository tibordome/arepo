/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/hii_mcs_anisotropic.c
 * \date        02/2019
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - 17.01.2022 Ported into current codebase
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(HII_MCS) && defined(HII_MCS_ANISO)

#define R_STROMGREN_FIRST_GUESS_FACTOR 0.9
#define R_STROMGREN_MINIMUM (1e-3 * PARSEC / All.UnitLength_in_cm)
#define CASE_B_RECOMBINATION_COEFF 2.56e-13  // cm^3 s^-1
#define RREC_HII_TOLERANCE_FACTOR 10
#define HII_MCS_CIE_FACTOR 1.05

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static int stromgren_evaluate(int target, int mode, int threadid);
static int stromgren_isactive(int auxid);
static int host_fetch_compare(const void *a, const void *b);

typedef struct
{
  int NConverged;
  int ConvergedPix[HII_MCS_N_PIX];
  int NumNgb[HII_MCS_N_PIX];
  int NumNgb_old[HII_MCS_N_PIX];
  MyFloat Rrec[HII_MCS_N_PIX];
  MyFloat Rrec_old[HII_MCS_N_PIX];
  MyFloat R_StromgrenPix[HII_MCS_N_PIX];
  MyFloat R_StromgrenPix_old[HII_MCS_N_PIX];
  MyFloat Left[HII_MCS_N_PIX];
  MyFloat Right[HII_MCS_N_PIX];
#ifdef HII_MCS_LR
  MyFloat Leftover;
#endif
} source_pixel_data;

static source_pixel_data *Source_Pixel;

static int *Sourcecell_pids;
static int NSources;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml[HII_MCS_N_PIX];
  MyFloat HsmlOld[HII_MCS_N_PIX];
  int ConvergedPix[HII_MCS_N_PIX];
  MyIDType SourceID;
  int Firstnode;
} data_in;

static data_in *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  int auxid  = SphP[i].SourceAuxID;
  in->Pos[0] = SphP[i].SourcePos[0];
  in->Pos[1] = SphP[i].SourcePos[1];
  in->Pos[2] = SphP[i].SourcePos[2];
  memcpy(in->Hsml, Source_Pixel[auxid].R_StromgrenPix, HII_MCS_N_PIX * sizeof(MyFloat));
  memcpy(in->HsmlOld, Source_Pixel[auxid].R_StromgrenPix_old, HII_MCS_N_PIX * sizeof(MyFloat));
  memcpy(in->ConvergedPix, Source_Pixel[auxid].ConvergedPix, HII_MCS_N_PIX * sizeof(int));
  in->SourceID = P[i].ID;

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Rrec[HII_MCS_N_PIX];
  int Ngb[HII_MCS_N_PIX];
} data_out;

static data_out *DataResult;

typedef struct
{
  double rate;
  double ratepos[3];
#ifdef HII_MCS_LR
  double luminosity;
#endif
#ifdef HII_MCS_RECORDS
  MyIDType StarID;
#endif
} star_to_cell_data;

static struct host_fetch
{
  int Task;
  int pid;
  double rate;
  double ratepos[3];
#ifdef HII_MCS_LR
  double luminosity;
#endif
#ifdef HII_MCS_RECORDS
  MyIDType StarID;
#endif
} * Host_fetch_out, *Host_fetch_in;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int pid, int mode)
{
  // n.b. we don't care about mode as we always increment
  int auxid = SphP[pid].SourceAuxID;
  for(int pix = 0; pix < HII_MCS_N_PIX; pix++)
    {
      Source_Pixel[auxid].NumNgb[pix] += out->Ngb[pix];  // Note that we always make adjustments relative to the last iteration,
      Source_Pixel[auxid].Rrec[pix] += out->Rrec[pix];   // hence += not =
    }
}

#define OMIT_GENERIC_COMM_PATTERN_FOR_GIVEN_PARTICLES
#include "../generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, TimeBinsHydro.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NSources)
          break;

        int pid = Sourcecell_pids[i];

        if(stromgren_isactive(i))
          stromgren_evaluate(pid, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        stromgren_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void place_hii_regions(void)
{
#ifdef HII_MCS_EVERY_SYNCPOINT
  place_hii_regions_internal();
#else
  double delta;

  if(All.ComovingIntegrationOn)
    delta = 1e3 * get_time_difference_in_Gyr(All.TimePrevHiiPlace, All.Time); /*Now in Myr*/
  else
    delta = (All.Time - All.TimePrevHiiPlace) / (SEC_PER_MEGAYEAR / All.UnitTime_in_s * All.HubbleParam); /*Now in Myr*/

  if(delta >= All.TimeBetweenHiiPlace)
    {
      place_hii_regions_internal();
      All.TimePrevHiiPlace = All.Time;
    }
#endif
}

void place_hii_regions_internal(void)
{
  TIMER_START(CPU_HII_FEEDBACK);

  int i, pix, pid, npleft, npleft_pix, iter = 0;
  int ion_stars_local;
  long long ntot;
  double host_rrec;
  double diff, diff_old, du;
  double t0, t1;
  MyFloat rstrom_temp;
  source_pixel_data *data_temp;

  mpi_printf("HII_MCS: Placing Hii regions\n");

  ion_stars_local = set_star_particle_ionization_properties();

  NSources         = 0;
  double s_hii_sum = 0;

  /*Reset ionisation tags and check for cells that should be ignored (e.g. temperature or density cuts)*/
  for(pid = 0; pid < NumGas; pid++)
    {
      if(SphP[pid].Utherm > (All.PhotoionizationEgySpec * HII_MCS_CIE_FACTOR)) /*Too hot*/
        {
          SphP[pid].StromgrenSourceID = HII_MCS_IGNORE_FLAG;
        }
#ifdef HII_MCS_DENSCUT
      else if((SphP[pid].Density * All.cf_a3inv) < All.HiiDensCut)
        {
          SphP[pid].StromgrenSourceID = HII_MCS_IGNORE_FLAG;
        }
#endif
      else
        SphP[pid].StromgrenSourceID = 0;
    }

  /*Only proceed with this stage if there are local sources*/
  if(ion_stars_local > 0)
    {
      /*First deal with host cell ionisation and potential source cells.*/
      Sourcecell_pids = (int *)mymalloc("Sourcecell_pids", ion_stars_local * sizeof(int));

      for(pid = 0; pid < NumGas; pid++)
        {
          if(SphP[pid].HostPhotonRate > 0)
            {
              s_hii_sum += SphP[pid].HostPhotonRate;
              if(SphP[pid].StromgrenSourceID == HII_MCS_IGNORE_FLAG) /*This cell is a source, but the cell itself is ignored */
                {
                  if(SphP[pid].R_Stromgren == 0) /*Initialise first guess*/
                    SphP[pid].R_Stromgren = All.R_Strom_prefactors *
                                            cbrt(SphP[pid].HostPhotonRate / (SphP[pid].Density * SphP[pid].Density)) * All.cf_atime;
                  SphP[pid].HostPhotonRate /=
                      HII_MCS_N_PIX; /*No contribution from host cell to recombination, so all photons enter pixels*/
                  SphP[pid].SourceAuxID     = NSources;
                  Sourcecell_pids[NSources] = pid;
                  NSources++;
                }
              else
                {
                  host_rrec = P[pid].Mass * SphP[pid].Density * All.Rrec_prefactors * All.cf_a3inv;
                  if(SphP[pid].HostPhotonRate > host_rrec)
                    {
                      SphP[pid].StromgrenSourceID = P[pid].ID; /*This cell has ionised itself*/
                      SphP[pid].HostPhotonRate -= host_rrec;   /*Remove contribution from host cell*/
                      SphP[pid].HostPhotonRate /= HII_MCS_N_PIX;
                      if(SphP[pid].HostPhotonRate >= host_rrec * All.MinimumPhotoRateFactor)
                        {
                          if(SphP[pid].R_Stromgren == 0)
                            SphP[pid].R_Stromgren = All.R_Strom_prefactors *
                                                    cbrt(SphP[pid].HostPhotonRate / (SphP[pid].Density * SphP[pid].Density)) *
                                                    All.cf_atime;
                          SphP[pid].SourceAuxID     = NSources;
                          Sourcecell_pids[NSources] = pid;
                          NSources++;
                        }
                      else
                        {
                          SphP[pid].HostPhotonRate =
                              0;  // If there are barely any photons left, don't bother carrying on beyond the host cell
                          SphP[pid].R_Stromgren = 0.0;
#ifdef HII_MCS_LR
                          SphP[pid].L_Hii = 0.0;
#endif
                        }
                    }
                  else  // Didn't manage to ionise host cell.
                    {
                      SphP[pid].StromgrenSourceID = 0;
                      SphP[pid].HostPhotonRate    = 0;
                      SphP[pid].R_Stromgren       = 0.0;
#ifdef HII_MCS_LR
                      SphP[pid].L_Hii = 0.0;
#endif
                    }
                }
            }
        }  // Loop over cells
    }      // Block only executed if local sources
  /* Now for those cells that have ionised themselves and still have photons to spare, we look for neighbours to ionise */

  int nsource_tot;
  MPI_Allreduce(&NSources, &nsource_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(nsource_tot > 0)  // If some task somewhere has sources
    {
      if(NSources > 0)  // If this task has sources
        Source_Pixel = (source_pixel_data *)mymalloc("Source_Pixel", NSources * sizeof(source_pixel_data));

      for(i = 0; i < NSources; i++)  // Note if no local sources, this block skipped
        {
          /* R_Stromgren is either left over from last time or calculated from host cell density if this is the first
          time for a star particle. We reduce initial guess for stromgren radius so ionisation front 'walks'
          outwards to deal with overlapping regions */
          rstrom_temp = SphP[Sourcecell_pids[i]].R_Stromgren * R_STROMGREN_FIRST_GUESS_FACTOR;
          if(rstrom_temp > (All.R_Stromgren_Max / All.cf_atime))
            rstrom_temp = All.R_Stromgren_Max * R_STROMGREN_FIRST_GUESS_FACTOR / All.cf_atime;

          Source_Pixel[i].NConverged = 0;

#ifdef HII_MCS_LR
          Source_Pixel[i].Leftover = 0;
#endif

          for(pix = 0; pix < HII_MCS_N_PIX; pix++)
            {
              Source_Pixel[i].ConvergedPix[pix]       = 0;
              Source_Pixel[i].Left[pix]               = 0;
              Source_Pixel[i].Right[pix]              = 0;
              Source_Pixel[i].NumNgb[pix]             = 0;
              Source_Pixel[i].NumNgb_old[pix]         = 0;
              Source_Pixel[i].Rrec[pix]               = 0;
              Source_Pixel[i].Rrec_old[pix]           = 0;
              Source_Pixel[i].R_StromgrenPix[pix]     = rstrom_temp;
              Source_Pixel[i].R_StromgrenPix_old[pix] = 0;
            }
        }

      generic_set_MaxNexport();

      /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
      do
        {
          t0 = second();

          generic_comm_pattern(NSources, kernel_local, kernel_imported);

          data_temp = Source_Pixel;
          for(i = 0, npleft = 0; i < NSources; i++)  // Note if no local sources, this block skipped
            {
              npleft_pix = 0;
              if(stromgren_isactive(i))
                {
                  pid = Sourcecell_pids[i];
                  for(pix = 0; pix < HII_MCS_N_PIX; pix++)
                    {
                      if(data_temp->ConvergedPix[pix] > 0)
                        continue;

                      data_temp->R_StromgrenPix_old[pix] = data_temp->R_StromgrenPix[pix];
                      diff                               = SphP[pid].HostPhotonRate - data_temp->Rrec[pix];
                      if(fabs(diff) > All.Rrec_Hii_tolerance)
                        {
                          /* need to redo this particle */
                          npleft_pix++;

                          if((abs(data_temp->NumNgb[pix] - data_temp->NumNgb_old[pix]) == 1) && data_temp->NumNgb_old[pix] > 0)
                            {
                              diff_old = SphP[pid].HostPhotonRate - data_temp->Rrec_old[pix];
                              if(((diff < 0) == (diff_old > 0)))  // sign_change
                                {
                                  /* this one should be ok */
                                  npleft_pix--;
                                  data_temp->ConvergedPix[pix] = 1; /* Mark as inactive */
                                  data_temp->NConverged += 1;
                                  continue;
                                }
                            }

                          if(data_temp->Left[pix] > 0 && data_temp->Right[pix] > 0)
                            {
                              if((data_temp->Right[pix] - data_temp->Left[pix]) <
                                 1.0e-3 * data_temp->Left[pix])  // Front has barely moved
                                {
                                  /* this one should be ok */
                                  npleft_pix--;
                                  data_temp->ConvergedPix[pix] = 1; /* Mark as inactive */
                                  data_temp->NConverged += 1;
                                  continue;
                                }
                            }

                          if(diff > All.Rrec_Hii_tolerance)  // Undershot
                            {
                              if(data_temp->R_StromgrenPix[pix] >= (All.R_Stromgren_Max / All.cf_atime))  // But already beyond maximum
                                {
                                  npleft_pix--;
                                  data_temp->ConvergedPix[i] = 1;
                                  data_temp->NConverged += 1;
#ifdef HII_MCS_LR
                                  data_temp->Leftover += diff;
#endif
                                  continue;
                                }
                              else  // prepare to guess larger next time
                                data_temp->Left[pix] = fmax(data_temp->R_StromgrenPix[pix], data_temp->Left[pix]);
                            }
                          else  // Overshot
                            {
                              if(data_temp->R_StromgrenPix[pix] <
                                 R_STROMGREN_MINIMUM)  // But already too small radius, avoid machine precision issues
                                {
                                  printf(
                                      "HII_MCS WARNING: Radius got too small, will be unconverged, so setting to 0 and continuing. "
                                      "ThisTask %d ID %d pix %d NumNgb %d HostPhotonRate %g Rrec %g diff %g R_Stromgren %g\n",
                                      ThisTask, P[pid].ID, pix, data_temp->NumNgb[pix], SphP[pid].HostPhotonRate, data_temp->Rrec[pix],
                                      diff, data_temp->R_StromgrenPix[pix]);
                                  data_temp->R_StromgrenPix[pix] = 0.0;
                                  npleft_pix--;
                                  data_temp->ConvergedPix[pix] = 1;
                                  data_temp->NConverged += 1;
                                  continue;
                                }
                              if(data_temp->Right[pix] != 0)
                                {
                                  if(data_temp->R_StromgrenPix[pix] < data_temp->Right[pix])
                                    data_temp->Right[pix] = data_temp->R_StromgrenPix[pix];
                                }
                              else
                                data_temp->Right[pix] = data_temp->R_StromgrenPix[pix];
                            }

                          if(iter >= MAXITER - 10)
                            {
                              printf(
                                  "iter %d pid=%d i=%d task=%d ID=%d pix %d R_Stromgren=%g Left=%g Right=%g Ngbs=%d Right-Left=%g "
                                  "Rrec=%g diff=%g tol=%g pos=(%g|%g|%g) sourcepos=(%g|%g|%g)\n",
                                  iter, pid, i, ThisTask, (int)P[pid].ID, pix, data_temp->R_StromgrenPix[pix], data_temp->Left[pix],
                                  data_temp->Right[pix], data_temp->NumNgb[pix], data_temp->Right[pix] - data_temp->Left[pix],
                                  data_temp->Rrec[pix], diff, All.Rrec_Hii_tolerance, P[pid].Pos[0], P[pid].Pos[1], P[pid].Pos[2],
                                  SphP[pid].SourcePos[0], SphP[pid].SourcePos[1], SphP[pid].SourcePos[2]);
                              myflush(stdout);
                            }

                          if(data_temp->Right[pix] > 0 && data_temp->Left[pix] > 0)
                            {
                              data_temp->R_StromgrenPix[pix] =
                                  pow(0.5 * (pow(data_temp->Left[pix], 3) + pow(data_temp->Right[pix], 3)), 1.0 / 3);
                            }
                          else
                            {
                              if(data_temp->Right[pix] == 0 && data_temp->Left[pix] == 0)
                                {
                                  terminate("should not occur");
                                }

                              if(data_temp->Right[pix] == 0 && data_temp->Left[pix] > 0)
                                {
                                  data_temp->R_StromgrenPix[pix] *= 1.1;
                                }

                              if(data_temp->Right[pix] > 0 && data_temp->Left[pix] == 0)
                                {
                                  data_temp->R_StromgrenPix[pix] /= 1.1;
                                }
                            }

                          data_temp->NumNgb_old[pix] = data_temp->NumNgb[pix];
                          data_temp->Rrec_old[pix]   = data_temp->Rrec[pix];
                        }
                      else
                        {
                          data_temp->ConvergedPix[pix] = 1; /* Mark as inactive */
                          data_temp->NConverged += 1;
                        }
                    }
                  if(npleft_pix > 0)
                    npleft++;
                }
              data_temp++;  // Point to next element of Source_Pixel
            }

          sumup_large_ints(1, &npleft, &ntot);

          t1 = second();

          if(ntot > 0)
            {
              iter++;

              if(iter > 0)
                mpi_printf("HII_MCS: R_Stromgren iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                           timediff(t0, t1));

              if(iter > MAXITER)
                terminate("failed to converge in neighbour iteration in place_hii_regions()\n");
            }
        }
      while(ntot > 0);

      for(i = 0; i < NSources; i++)  // Note if no local sources, this block skipped
        {
          rstrom_temp = 0;
          for(pix = 0; pix < HII_MCS_N_PIX; pix++)
            {
              rstrom_temp = fmin(rstrom_temp, Source_Pixel[i].R_StromgrenPix[pix]);
#ifdef HII_MCS_RECORDS
              SphP[Sourcecell_pids[i]].R_StromgrenArr[pix] = Source_Pixel[i].R_StromgrenPix[pix];
#endif
            }
          SphP[Sourcecell_pids[i]].R_Stromgren = rstrom_temp;
#ifdef HII_MCS_LR
          SphP[Sourcecell_pids[i]].L_Hii *= Source_Pixel[i].Leftover;  // Now in units of 1e49 ergs s^-1
#endif
        }

      if(NSources > 0)
        myfree(Source_Pixel);
    }  // nsource_tot > 0

  if(ion_stars_local > 0)
    myfree(Sourcecell_pids);

  /* Now adjust energies of gas cells and count ionized*/
  int nstrom, nstrom_tot;
  double mstrom, mstrom_tot;
  nstrom = 0;
  mstrom = 0;
  for(pid = 0; pid < NumGas; pid++)
    {
      if((SphP[pid].StromgrenSourceID > 0) && (SphP[pid].StromgrenSourceID != HII_MCS_IGNORE_FLAG))
        {
          nstrom++;
          mstrom += P[pid].Mass;
          du = All.PhotoionizationEgySpec - SphP[pid].Utherm;
          if(du > 0)
            {
              SphP[pid].Utherm += du;
              SphP[pid].Energy += P[pid].Mass * du * All.cf_atime * All.cf_atime;
              set_pressure_of_cell(pid);
            }
        }
    }

  int nactive_tot;
  MPI_Allreduce(&ion_stars_local, &nactive_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  double s_hii_mean = 0;
  nstrom_tot = mstrom_tot = 0;
  if(nactive_tot > 0)
    {
      MPI_Reduce(&nstrom, &nstrom_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&mstrom, &mstrom_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&s_hii_sum, &s_hii_mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      s_hii_mean /= nactive_tot;
    }
  mpi_printf("HII_MCS: %d sources ionized %d cells, ionized mass %g, s_hii_mean %g\n", nactive_tot, nstrom_tot, mstrom_tot,
             s_hii_mean);
#ifdef HII_MCS_LOG
  if(ThisTask == 0)
    {
      fprintf(FdHii, "%14e %14d %14d %14e %14e\n", All.Time, nactive_tot, nstrom_tot, mstrom_tot, s_hii_mean);
      myflush(FdHii);
    }
#endif

  /* collect some timing information */

  TIMER_STOP(CPU_HII_FEEDBACK);
}

/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int stromgren_evaluate(int target, int mode, int threadid)
{
  int j, n, pix;
  int nfound, numnodes, *firstnode, *converged;
  int numngb[HII_MCS_N_PIX];
  double rsearch, r;
  MyFloat *hpix, *holdpix;
  MyFloat rrec[HII_MCS_N_PIX];
  double dx[3];
  MyDouble *pos;

  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos               = target_data->Pos;
  MyIDType sourceid = target_data->SourceID;
  hpix              = target_data->Hsml;
  holdpix           = target_data->HsmlOld;
  converged         = target_data->ConvergedPix;

  rsearch = 0;
  for(pix = 0; pix < HII_MCS_N_PIX; pix++)
    {
      rsearch = fmax(rsearch, hpix[pix]);
      rsearch = fmax(rsearch, holdpix[pix]);
    }

  memset(numngb, 0, HII_MCS_N_PIX * sizeof(int));
  memset(rrec, 0, HII_MCS_N_PIX * sizeof(MyFloat));

  /*The front is increasing, flag everything within h as ionised*/
  nfound = ngb_treefind_variable_threads(pos, rsearch, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      if(P[j].Mass == 0 && P[j].ID == 0)  // Swallowed or elimated so ignore
        continue;

      if(SphP[j].StromgrenSourceID == HII_MCS_IGNORE_FLAG)
        continue;

      dx[0] = nearest_x(pos[0] - P[j].Pos[0]);
      dx[1] = nearest_y(pos[1] - P[j].Pos[1]);
      dx[2] = nearest_z(pos[2] - P[j].Pos[2]);

      vec2pix_ring(dx, &r, &pix);  // Find radius and healpix bin

      if(converged[pix] > 0)
        continue;

      if(r < hpix[pix])
        {
          if(SphP[j].StromgrenSourceID > 0)  // Already being ionised by a source, possibly this one
            continue;

          SphP[j].StromgrenSourceID = sourceid;  // Mark as ionized by this cell
          numngb[pix] += 1;
          rrec[pix] += All.Rrec_prefactors * P[j].Mass * SphP[j].Density * All.cf_a3inv - SphP[j].HostPhotonRate;
        }
      else if((r < holdpix[pix]) && (SphP[j].StromgrenSourceID == sourceid) &&
              (P[j].ID !=
               sourceid)) /*Reset the ionization state as we retract the front if this cell was being ionisised by this source*/
        {
          SphP[j].StromgrenSourceID = 0;  // Mark as no longer being ionised by this cell
          numngb[pix] -= 1;
          rrec[pix] -= All.Rrec_prefactors * P[j].Mass * SphP[j].Density * All.cf_a3inv + SphP[j].HostPhotonRate;
        }
    }

  for(pix = 0; pix < HII_MCS_N_PIX; pix++)
    {
      out.Ngb[pix]  = numngb[pix];
      out.Rrec[pix] = rrec[pix];
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

static int stromgren_isactive(int i)
{
  if(Source_Pixel[i].NConverged >= HII_MCS_N_PIX)
    return 0;

  if(SphP[Sourcecell_pids[i]].HostPhotonRate > 0)
    return 1;

  return 0;
}

void init_hii(void)
{
  double meanweight;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));  // Fully ionized

  All.PhotoionizationEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.PhotoionizationGasTemp;
  All.PhotoionizationEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  /* Recombination rate of one cell at density of RREC_HII_TOLERANCE_FACTOR cm^-3, in s^-1 */
  All.Rrec_Hii_tolerance = RREC_HII_TOLERANCE_FACTOR * CASE_B_RECOMBINATION_COEFF *
                           (All.TargetGasMass * All.UnitMass_in_g / All.HubbleParam) * HYDROGEN_MASSFRAC / PROTONMASS;
  All.Rrec_Hii_tolerance /= 1e49;  // Rescale to avoid issues if MyFloat = float.

  All.Rrec_prefactors =
      CASE_B_RECOMBINATION_COEFF * HYDROGEN_MASSFRAC * HYDROGEN_MASSFRAC / (PROTONMASS * PROTONMASS);  // cm^3 g^-2 s^-1
  All.Rrec_prefactors *= (All.UnitMass_in_g * All.UnitDensity_in_cgs);
  All.Rrec_prefactors *= All.HubbleParam;
  All.Rrec_prefactors /= 1e49;  // Rescale to avoid issues if MyFloat = float.
  All.R_Strom_prefactors =
      cbrt(3e49 * PROTONMASS * PROTONMASS /
           (4.0 * M_PI * HYDROGEN_MASSFRAC * HYDROGEN_MASSFRAC * CASE_B_RECOMBINATION_COEFF * All.UnitDensity_in_cgs *
            All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.HubbleParam * All.HubbleParam)) /
      (All.UnitLength_in_cm / All.HubbleParam);

  for(int i = 0; i < NumGas; i++)
    SphP[i].StromgrenSourceID = 0;

#ifndef HII_MCS_EVERY_SYNCPOINT
  if(All.ComovingIntegrationOn)
    All.TimePrevHiiPlace = 0.0;
  else
    All.TimePrevHiiPlace = -(All.Time - All.TimePrevHiiPlace) * (SEC_PER_MEGAYEAR / All.UnitTime_in_s * All.HubbleParam);
#endif
#ifdef HII_MCS_DENSCUT
  All.HiiDensCut *= ((PROTONMASS / HYDROGEN_MASSFRAC) / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam));
#endif

#ifdef HII_MCS_LR
  /* Factor to multiply luminosity / r^2 to get energy density in cgs. Factor of 1e49 undoes rescaling we carried out earlier */
  All.Factor_Hii = 1e49 / (4.0 * M_PI * CLIGHT * All.UnitLength_in_cm * All.UnitLength_in_cm / All.HubbleParam / All.HubbleParam);
#endif
}

int set_star_particle_ionization_properties(void)
{
  int i, j;
  int nexport, nimport, ngrp, recvTask;
  int pid, auxid, ion_stars, ion_stars_local;
  mesh_search_data *searchdata;
  star_to_cell_data *star_data;
  MPI_Status status;

  for(pid = 0; pid < NumGas; pid++)
    {
      SphP[pid].HostPhotonRate = 0.0;
      SphP[pid].SourcePos[0]   = 0.0;
      SphP[pid].SourcePos[1]   = 0.0;
      SphP[pid].SourcePos[2]   = 0.0;
#ifdef HII_MCS_LR
      SphP[pid].L_Hii = 0.0;
#endif
#ifdef HII_MCS_RECORDS
      SphP[pid].N_Sources            = 0;
      SphP[pid].StarIDArr[0]         = 0;
      SphP[pid].StarIDArr[1]         = 0;
      SphP[pid].StarIDArr[2]         = 0;
      SphP[pid].HiiRecombinationRate = P[pid].Mass * SphP[pid].Density * All.Rrec_prefactors * All.cf_a3inv;
#endif
    }

#ifdef HII_MCS_TEST
  /* Fixed photon rate for tests, to be replaced with dependent rate*/
  ion_stars = 0;
  for(auxid = 0; auxid < N_star; auxid++)
    {
      StarP[auxid].S_Hii = 0.25;  // 1e49 photons per second, quarter of O star.
      ion_stars++;
    }
#else
  ion_stars = update_photon_rates();
#endif

  /*Find host cells of active stars*/

  searchdata = (mesh_search_data *) mymalloc("searchdata", ion_stars * sizeof(mesh_search_data));
  star_data  = (star_to_cell_data *) mymalloc("star_data", ion_stars * sizeof(star_to_cell_data));

  for(auxid = 0, j = 0; auxid < N_star; auxid++)
    {
      if((StarP[auxid].S_Hii > 0))
        {
          pid                     = StarP[auxid].PID;
          searchdata[j].Pos[0]    = P[pid].Pos[0];
          searchdata[j].Pos[1]    = P[pid].Pos[1];
          searchdata[j].Pos[2]    = P[pid].Pos[2];
          star_data[j].rate       = StarP[auxid].S_Hii;
          star_data[j].ratepos[0] = P[pid].Pos[0] * star_data[j].rate;
          star_data[j].ratepos[1] = P[pid].Pos[1] * star_data[j].rate;
          star_data[j].ratepos[2] = P[pid].Pos[2] * star_data[j].rate;
#ifdef HII_MCS_LR
          star_data[j].luminosity = StarP[auxid].EnergyPerPhoton * star_data[j].rate;
#endif
#ifdef HII_MCS_RECORDS
          star_data[j].StarID = P[pid].ID;
#endif
          j++;
        }
    }

  assert(j == ion_stars);

  find_nearest_meshpoint_global(searchdata, ion_stars, 0, 0);

  /* Check how many particles need to be exported */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0, nexport = 0; i < ion_stars; i++)
    {
      if(searchdata[i].Task != ThisTask)
        {
          Send_count[searchdata[i].Task] += 1;
          nexport++;
        }
    }

  /* Send the export counts to other tasks and get receive counts */
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  Host_fetch_in  = (struct host_fetch *)mymalloc("Host_fetch_in", nimport * sizeof(struct host_fetch));
  Host_fetch_out = (struct host_fetch *)mymalloc("Host_fetch_out", nexport * sizeof(struct host_fetch));

  ion_stars_local = 0;

  for(i = 0, j = 0; i < ion_stars; i++)
    {
      if(searchdata[i].Task == ThisTask)
        {
          pid = searchdata[i].u.Index;
          SphP[pid].HostPhotonRate += star_data[i].rate;
          SphP[pid].SourcePos[0] += star_data[i].ratepos[0];
          SphP[pid].SourcePos[1] += star_data[i].ratepos[1];
          SphP[pid].SourcePos[2] += star_data[i].ratepos[2];
#ifdef HII_MCS_LR
          SphP[pid].L_Hii += star_data[i].luminosity;
#endif
#ifdef HII_MCS_RECORDS
          if(SphP[pid].N_Sources < 3)
            SphP[pid].StarIDArr[SphP[pid].N_Sources] = star_data[i].StarID;
          SphP[pid].N_Sources++;
#endif

          ion_stars_local++;
        }
      else
        {
          Host_fetch_out[j].Task       = searchdata[i].Task;
          Host_fetch_out[j].pid        = searchdata[i].u.Index;
          Host_fetch_out[j].rate       = star_data[i].rate;
          Host_fetch_out[j].ratepos[0] = star_data[i].ratepos[0];
          Host_fetch_out[j].ratepos[1] = star_data[i].ratepos[1];
          Host_fetch_out[j].ratepos[2] = star_data[i].ratepos[2];
#ifdef HII_MCS_LR
          Host_fetch_out[j].luminosity = star_data[i].luminosity;
#endif
#ifdef HII_MCS_RECORDS
          Host_fetch_out[j].StarID = star_data[i].StarID;
#endif
          j++;
        }
    }

  /* Get export buffer into correct order */
  qsort(Host_fetch_out, nexport, sizeof(struct host_fetch), host_fetch_compare);

  /* Exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&Host_fetch_out[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct host_fetch), MPI_BYTE,
                           recvTask, TAG_DENS_A, &Host_fetch_in[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct host_fetch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
        }
    }

  myfree(Host_fetch_out);

  /*Deal with imported particles */
  for(i = 0; i < nimport; i++)
    {
      pid = Host_fetch_in[i].pid;
      SphP[pid].HostPhotonRate += Host_fetch_in[i].rate;
      SphP[pid].SourcePos[0] += Host_fetch_in[i].ratepos[0];
      SphP[pid].SourcePos[1] += Host_fetch_in[i].ratepos[1];
      SphP[pid].SourcePos[2] += Host_fetch_in[i].ratepos[2];
#ifdef HII_MCS_LR
      SphP[pid].L_Hii += Host_fetch_in[i].luminosity;
#endif
#ifdef HII_MCS_RECORDS
      if(SphP[pid].N_Sources < 3)
        SphP[pid].StarIDArr[SphP[pid].N_Sources] = Host_fetch_in[i].StarID;
      SphP[pid].N_Sources++;
#endif
      ion_stars_local++;
    }

  myfree(Host_fetch_in);
  myfree(star_data);
  myfree(searchdata);

  for(pid = 0; pid < NumGas; pid++)
    {
      if(SphP[pid].HostPhotonRate > 0)
        {
          SphP[pid].SourcePos[0] /= SphP[pid].HostPhotonRate;
          SphP[pid].SourcePos[1] /= SphP[pid].HostPhotonRate;
          SphP[pid].SourcePos[2] /= SphP[pid].HostPhotonRate;
#ifdef HII_MCS_LR
          SphP[pid].L_Hii /= SphP[pid].HostPhotonRate;  // Temporarily this variable stores the mean photon energy in ergs
#endif
        }
    }

  return ion_stars_local;
}

/* Sorts by task (we don't care about index)*/
static int host_fetch_compare(const void *a, const void *b)
{
  if(((struct host_fetch *)a)->Task < (((struct host_fetch *)b)->Task))
    return -1;

  if(((struct host_fetch *)a)->Task > (((struct host_fetch *)b)->Task))
    return +1;

  return 0;
}

#ifdef HII_MCS_LR
double calculate_uv_background_boost_factor(int i)
{
  double boost = 1.0;

  if((SphP[i].StromgrenSourceID == 0) || (SphP[i].StromgrenSourceID == HII_MCS_IGNORE_FLAG))
    {
      boost = 1.0 + SphP[i].EnergyDensHii / All.UVBEnergyDensHii;
      assert(boost >= 1.0);  // Overly cautious, can be removed if confident EnergyDensHii is correct.
    }

  return boost;
}
#endif
#endif
