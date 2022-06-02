/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_gradients.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "allvars.h"
#include "proto.h"

#if(defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)) && !defined(GRADIENTS_GREEN_GAUSS) && !defined(ONEDIMS) && \
    !defined(AMR_GRADIENTS)

static double *minvalues, *maxvalues;

static void limit_gradients(void);

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
static void correct_for_reflective_boundaries(double *ValueOther, double Value, int type, unsigned int *image_flags);
#endif

#if defined(OUTPUT_DIVVEL) || defined(MHD) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES)
static void compute_divergences(void);
#endif

/*! \brief Adds row to another one in matrix equation.
 *
 *  Auxiliary routine to solve_matrix_problem.
 *
 *  \param[in, out] X Matrix.
 *  \param[in, out] y Vector.
 *  \param[in] source_row Index of row that should be added.
 *  \param[in] fac Factor by which row is multiplied before adding.
 *  \param[in] target_row Index of row to which to add source row.
 *
 *  \return void
 */
static inline void add_row(double X[NUMDIMS][NUMDIMS], double y[NUMDIMS], int source_row, double fac, int target_row)
{
  y[target_row] += fac * y[source_row];

  for(int i = 0; i < NUMDIMS; i++)
    {
      X[target_row][i] += fac * X[source_row][i];
    }
}

/*! \brief Solve a matrix problem X*grad = y.
 *
 *   Note that we know here that X is symmetric, and that we can pivot on the
 *   diagonal elements.
 *
 *  \param[in, out] x Matrix.
 *  \param[in, out] y Vector.
 *  \param[out] grad Gradient.
 *
 */
static void solve_matrix_problem(double X[NUMDIMS][NUMDIMS], double y[NUMDIMS], double grad[NUMDIMS])
{
#if NUMDIMS == 2
  int perm[NUMDIMS];

  if(fabs(X[0][0]) > fabs(X[1][1]))
    {
      perm[0] = 0;
      perm[1] = 1;
    }
  else
    {
      perm[0] = 1;
      perm[1] = 0;
    }

  add_row(X, y, perm[0], -X[perm[1]][perm[0]] / X[perm[0]][perm[0]], perm[1]);

  grad[perm[1]] = y[perm[1]] / X[perm[1]][perm[1]];
  grad[perm[0]] = (y[perm[0]] - X[perm[0]][perm[1]] * grad[perm[1]]) / X[perm[0]][perm[0]];

#else

  int perm[NUMDIMS];

  if(fabs(X[2][2]) > fabs(X[1][1]) && fabs(X[2][2]) > fabs(X[0][0]))
    {
      perm[0] = 2;
      perm[1] = 0;
      perm[2] = 1;
    }
  else if(fabs(X[1][1]) > fabs(X[0][0]))
    {
      perm[0] = 1;
      perm[1] = 0;
      perm[2] = 2;
    }
  else
    {
      perm[0] = 0;
      perm[1] = 1;
      perm[2] = 2;
    }

  add_row(X, y, perm[0], -X[perm[1]][perm[0]] / X[perm[0]][perm[0]], perm[1]);
  add_row(X, y, perm[0], -X[perm[2]][perm[0]] / X[perm[0]][perm[0]], perm[2]);

  if(fabs(X[perm[1]][perm[1]]) < fabs(X[perm[2]][perm[2]]))
    {
      int p   = perm[1];
      perm[1] = perm[2];
      perm[2] = p;
    }

  add_row(X, y, perm[1], -X[perm[2]][perm[1]] / X[perm[1]][perm[1]], perm[2]);

  grad[perm[2]] = y[perm[2]] / X[perm[2]][perm[2]];
  grad[perm[1]] = (y[perm[1]] - X[perm[1]][perm[2]] * grad[perm[2]]) / X[perm[1]][perm[1]];
  grad[perm[0]] = (y[perm[0]] - X[perm[0]][perm[1]] * grad[perm[1]] - X[perm[0]][perm[2]] * grad[perm[2]]) / X[perm[0]][perm[0]];

#endif
}

/*! \brief Loop through all active cells and calculate gradients.
 *
 *  \return void
 */
void calculate_gradients(void)
{
  TIMER_START(CPU_GRADIENTS);

  mpi_printf("VORONOI: Calculating Gradients...\n");

  minvalues = (double *)mymalloc("gradmin", NumGas * N_Grad * sizeof(double));
  maxvalues = (double *)mymalloc("gradmax", NumGas * N_Grad * sizeof(double));

  struct matrix_vec_data
  {
    double X[NUMDIMS][NUMDIMS]; /* input matrix */
    double y[NUMDIMS];          /* input vector */
    double grad[NUMDIMS];       /* output */
  } * mdata;

  mdata = (struct matrix_vec_data *)mymalloc("mdata", N_Grad * sizeof(struct matrix_vec_data));

  double *Value = (double *)mymalloc("Value", N_Grad * sizeof(double));

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef DEREFINE_GENTLY
      SphP[i].DoNotDerefFlag = 0;
#endif

      for(int j = 0; j < 3; j++)
        SphP[i].Grad.Center[j] = SphP[i].Center[j]; /* center of mass of cell at time of computation of gradient */

      for(int k = 0; k < N_Grad; k++)
        {
          minvalues[i * N_Grad + k] = +MAX_REAL_NUMBER;
          maxvalues[i * N_Grad + k] = -MAX_REAL_NUMBER;

          if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) ||
             (grad_elements[k].type == GRADIENT_TYPE_VELZ))
            {
              Value[k] = *(MyFloat *)(((char *)(&P[i])) + grad_elements[k].offset) / All.cf_atime;
            }
          else
            Value[k] = *(MyFloat *)(((char *)(&SphP[i])) + grad_elements[k].offset);
        }

      double *Center = SphP[i].Center;

      /* reset matrix and vector to 0 */
      memset(mdata, 0, N_Grad * sizeof(struct matrix_vec_data));

#ifdef REFLECTIVE_X
      int OutFlowX = 0;
#endif
#ifdef REFLECTIVE_Y
      int OutFlowY = 0;
#endif
#ifdef REFLECTIVE_Z
      int OutFlowZ = 0;
#endif

      int q = SphP[i].first_connection;

      while(q >= 0)
        {
          int dp       = DC[q].dp_index;
          int vf       = DC[q].vf_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              /* cell has been removed */
              q = DC[q].next;
              continue;
            }

          if(Mesh.VF[vf].area > 1e-10 * SphP[i].SurfaceArea)
            {
              double *CenterOther, Mirror[3];

              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;

#ifdef REFLECTIVE_X
              if((Mesh.DP[dp].image_flags & REFL_X_FLAGS) && (Mesh.DP[dp].image_flags & OUTFLOW_X))
                OutFlowX = 1;
#endif
#ifdef REFLECTIVE_Y
              if((Mesh.DP[dp].image_flags & REFL_Y_FLAGS) && (Mesh.DP[dp].image_flags & OUTFLOW_Y))
                OutFlowY = 1;
#endif
#ifdef REFLECTIVE_Z
              if((Mesh.DP[dp].image_flags & REFL_Z_FLAGS) && (Mesh.DP[dp].image_flags & OUTFLOW_Z))
                OutFlowZ = 1;
#endif

              if(Mesh.DP[dp].task == ThisTask)
                {
#ifndef VORONOI_STATIC_MESH
                  if(P[particle].Ti_Current != All.Ti_Current)
                    terminate("surprise! we don't expect this here anymore");
#endif

                  if(P[particle].ID == P[i].ID)
                    {
                      /* mirrored cell, we have to mirror the Center */

                      /* calculate normal vector of the interface */
                      double nx = Mesh.DP[dp].x - P[i].Pos[0];
                      double ny = Mesh.DP[dp].y - P[i].Pos[1];
                      double nz = Mesh.DP[dp].z - P[i].Pos[2];

                      /* perpendicular on the surface */
                      double nn = sqrt(nx * nx + ny * ny + nz * nz);
                      nx /= nn;
                      ny /= nn;
                      nz /= nn;
                      double fx = (Center[0] - Mesh.VF[vf].cx);
                      double fy = (Center[1] - Mesh.VF[vf].cy);
                      double fz = (Center[2] - Mesh.VF[vf].cz);
                      double ff = (fx * nx + fy * ny + fz * nz);

                      double px = Center[0] - ff * nx;
                      double py = Center[1] - ff * ny;
                      double pz = Center[2] - ff * nz;

                      Mirror[0]   = 2. * px - Center[0];
                      Mirror[1]   = 2. * py - Center[1];
                      Mirror[2]   = 2. * pz - Center[2];
                      CenterOther = Mirror;
                    }
                  else
                    CenterOther = SphP[particle].Center;
                }
              else
                CenterOther = PrimExch[particle].Center;

              double norm[3];
              norm[0] = nearest_x(CenterOther[0] - Center[0]);
              norm[1] = nearest_y(CenterOther[1] - Center[1]);
              norm[2] = nearest_z(CenterOther[2] - Center[2]);

              double dist    = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
              double distinv = 1.0 / dist;
              norm[0] *= distinv;
              norm[1] *= distinv;
              norm[2] *= distinv;

              double weight = Mesh.VF[vf].area;

              for(int k = 0; k < N_Grad; k++)
                {
                  double ValueOther;

                  if(Mesh.DP[dp].task == ThisTask)
                    {
                      if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) ||
                         (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                        {
                          ValueOther = *(MyFloat *)(((char *)(&P[particle])) + grad_elements[k].offset);
                        }
                      else
                        ValueOther = *(MyFloat *)(((char *)(&SphP[particle])) + grad_elements[k].offset);
                    }
                  else
                    {
                      ValueOther = *(MyFloat *)(((char *)(&PrimExch[particle])) + grad_elements[k].offset_exch);
                    }

                  if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) ||
                     (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                    {
                      ValueOther /= All.cf_atime;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
                      correct_for_reflective_boundaries(&ValueOther, Value[k], grad_elements[k].type, &Mesh.DP[dp].image_flags);
#endif
                      if(grad_elements[k].type == GRADIENT_TYPE_VELX)
                        ValueOther += norm[0] * dist * All.cf_atime * All.cf_Hrate;
                      else if(grad_elements[k].type == GRADIENT_TYPE_VELY)
                        ValueOther += norm[1] * dist * All.cf_atime * All.cf_Hrate;
                      else if(grad_elements[k].type == GRADIENT_TYPE_VELZ)
                        ValueOther += norm[2] * dist * All.cf_atime * All.cf_Hrate;
                    }

                  double fac = weight * (ValueOther - Value[k]) / dist;

                  for(int ia = 0; ia < NUMDIMS; ia++)
                    {
                      mdata[k].y[ia] += fac * norm[ia];

                      for(int ib = 0; ib < NUMDIMS; ib++)
                        mdata[k].X[ia][ib] += weight * norm[ia] * norm[ib];
                    }

                  if(ValueOther < minvalues[i * N_Grad + k])
                    minvalues[i * N_Grad + k] = ValueOther;

                  if(ValueOther > maxvalues[i * N_Grad + k])
                    maxvalues[i * N_Grad + k] = ValueOther;
                }
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      for(int k = 0; k < N_Grad; k++)
        {
          solve_matrix_problem(mdata[k].X, mdata[k].y, mdata[k].grad);

          MySingle *data = (MySingle *)(((char *)(&(SphP[i].Grad))) + grad_elements[k].offset_grad);
          for(int j = 0; j < NUMDIMS; j++)
            data[j] = mdata[k].grad[j];
          for(int j = NUMDIMS; j < 3; j++)
            data[j] = 0.;

#ifdef REFLECTIVE_X
          if(OutFlowX)
            data[0] = 0;
#endif
#ifdef REFLECTIVE_Y
          if(OutFlowY)
            data[1] = 0;
#endif
#ifdef REFLECTIVE_Z
          if(OutFlowZ)
            data[2] = 0;
#endif
        }
    }

  myfree(Value);
  myfree(mdata);

#ifdef DEREFINE_GENTLY
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef SINKS
      /* cells that are inside the sink accretion radius are eligible for derefinement */
      if(SphP[i].InAccrRadius)
        continue;
#endif

      if(SphP[i].Density > GENTLE_DEREFINE_FACTOR * minvalues[i * N_Grad + GDensity->idx])
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Density < maxvalues[i * N_Grad + GDensity->idx] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Utherm > GENTLE_DEREFINE_FACTOR * minvalues[i * N_Grad + GUtherm->idx])
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Utherm < maxvalues[i * N_Grad + GUtherm->idx] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[0] / All.cf_atime > GENTLE_DEREFINE_FACTOR * minvalues[i * N_Grad + GVelx->idx])
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[0] / All.cf_atime < maxvalues[i * N_Grad + GVelx->idx] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[1] / All.cf_atime > GENTLE_DEREFINE_FACTOR * minvalues[i * N_Grad + GVely->idx])
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[1] / All.cf_atime < maxvalues[i * N_Grad + GVely->idx] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[2] / All.cf_atime > GENTLE_DEREFINE_FACTOR * minvalues[i * N_Grad + GVelz->idx])
        SphP[i].DoNotDerefFlag = 1;

      if(P[i].Vel[2] / All.cf_atime < maxvalues[i * N_Grad + GVelz->idx] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Pressure > GENTLE_DEREFINE_FACTOR * minvalues[i * N_Grad + GPressure->idx])
        SphP[i].DoNotDerefFlag = 1;

      if(SphP[i].Pressure < maxvalues[i * N_Grad + GPressure->idx] / GENTLE_DEREFINE_FACTOR)
        SphP[i].DoNotDerefFlag = 1;
    }
#endif

#ifdef TVD_SLOPE_LIMITER
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      SphP[i].GradUl = SphP[i].Grad;
    }
#endif

#ifdef NUCLEAR_NETWORK
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(int j = 0; j < 3; j++)
        {
          SphP[i].Grad.drhoU[j]   = SphP[i].Grad.drho[j];
          SphP[i].Grad.dpressU[j] = SphP[i].Grad.dpress[j];
        }
    }
#endif

#ifdef MHD
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].CurlB[0] = SphP[i].Grad.dB[2][1] - SphP[i].Grad.dB[1][2];
      SphP[i].CurlB[1] = SphP[i].Grad.dB[0][2] - SphP[i].Grad.dB[2][0];
      SphP[i].CurlB[2] = SphP[i].Grad.dB[1][0] - SphP[i].Grad.dB[0][1];
    }
#endif

#ifdef TURB_APPROX_MCS
#ifdef TURB_APPROX_MCS_GRAD_UNLIM
  /* Save velocity gradients before slope limiting */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      memcpy(SphP[i].dvel_unlim, SphP[i].Grad.dvel, 9 * sizeof(MySingle));
    }
#endif  // TURB_APPROX_MCS_GRAD_UNLIM
#elif SFR_MCS_RATE_CRITERIA > 0
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      SphP[i].gradv_sq = 0;
      for(int j = 0; j < 3; j++)
        {
          for(int k = 0; k < 3; k++)
            SphP[i].gradv_sq += SphP[i].Grad.dvel[j][k] * SphP[i].Grad.dvel[j][k];
        }
    }
#endif

#ifndef UNLIMITED_GRADIENTS
  limit_gradients();
#endif

#ifdef REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
  /* compute magnitude of curl */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      double curlx = SphP[i].Grad.dvel[2][1] - SphP[i].Grad.dvel[1][2];
      double curly = SphP[i].Grad.dvel[0][2] - SphP[i].Grad.dvel[2][0];
      double curlz = SphP[i].Grad.dvel[1][0] - SphP[i].Grad.dvel[0][1];

      SphP[i].CurlVel = sqrt(curlx * curlx + curly * curly + curlz * curlz);
    }
#endif

#if defined(CALCULATE_QUANTITIES_IN_POSTPROCESS) || defined(AB_TURB)

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double curlx = SphP[i].Grad.dvel[2][1] - SphP[i].Grad.dvel[1][2];
      double curly = SphP[i].Grad.dvel[0][2] - SphP[i].Grad.dvel[2][0];
      double curlz = SphP[i].Grad.dvel[1][0] - SphP[i].Grad.dvel[0][1];

      SphP[i].Vorticity[0] = curlx;
      SphP[i].Vorticity[1] = curly;
      SphP[i].Vorticity[2] = curlz;
    }
#endif

  myfree(maxvalues);
  myfree(minvalues);

#if defined(OUTPUT_DIVVEL) || defined(MHD) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES) || \
    defined(SF_STELLAR_MASS_TO_GAS_MASS_RATIO)
  compute_divergences();
#endif

#ifdef NO_SCALAR_GRADIENTS
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(int k = 0; k < N_Scalar; k++)
        {
          SphP[i].Grad.dscalars[k][0] = 0;
          SphP[i].Grad.dscalars[k][1] = 0;
          SphP[i].Grad.dscalars[k][2] = 0;
        }
    }
#endif /* #ifdef NO_SCALAR_GRADIENTS */

  TIMER_STOP(CPU_GRADIENTS);
}

#if defined(OUTPUT_DIVVEL) || defined(MHD) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES) || \
    defined(SF_STELLAR_MASS_TO_GAS_MASS_RATIO)
void compute_divergences(void)
{
  mpi_printf("VORONOI: Computing divergences... \n");
#if defined(EVALPOTENTIAL) && defined(SINK_PARTICLES)
  mpi_printf("VORONOI: Also geting POTENTIAL PEAKS \n");
#endif

  exchange_primitive_variables_and_gradients();

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#if defined(OUTPUT_DIVVEL) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES)
      SphP[i].DivVel = 0;
#endif
#ifdef MHD
      SphP[i].DivB    = 0;
      SphP[i].DivBalt = 0;
#endif

#ifdef SF_STELLAR_MASS_TO_GAS_MASS_RATIO
      SphP[i].VelDisp = 0;
      int NNgb        = 0;
#endif

      double *CenterOther, Mirror[3];
#if defined(OUTPUT_DIVVEL) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES) || defined(SF_STELLAR_MASS_TO_GAS_MASS_RATIO)
      MyFloat *VelOther;
#endif
#ifdef MHD
      MyFloat *BOther, B[3], dB[3];
      struct grad_data *GradOther;
#endif
#if defined(EVALPOTENTIAL) && defined(SINK_PARTICLES)
      SphP[i].PotentialPeak = 1; /* starts off true */
      double PotenOther;
#endif

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          int dp       = DC[q].dp_index;
          int vf       = DC[q].vf_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              /* cell has been removed */
              q = DC[q].next;
              continue;
            }

          if(Mesh.VF[vf].area > 1e-10 * SphP[i].SurfaceArea)
            {
#ifdef MHD
              double dx = nearest_x(Mesh.VF[vf].cx - SphP[i].Center[0]);
              double dy = nearest_y(Mesh.VF[vf].cy - SphP[i].Center[1]);
              double dz = nearest_z(Mesh.VF[vf].cz - SphP[i].Center[2]);

              for(int j = 0; j < 3; j++)
                B[j] = SphP[i].B[j] + SphP[i].Grad.dB[j][0] * dx + SphP[i].Grad.dB[j][1] * dy + SphP[i].Grad.dB[j][2] * dz;
#endif

              if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                particle -= NumGas;

              if(Mesh.DP[dp].task == ThisTask)
                {
                  if(P[particle].ID == P[i].ID)
                    {
                      /* mirrored cell, we have to mirror the Center */
                      /* calculate normal vector of the interface */
                      double nx = Mesh.DP[dp].x - P[i].Pos[0];
                      double ny = Mesh.DP[dp].y - P[i].Pos[1];
                      double nz = Mesh.DP[dp].z - P[i].Pos[2];
                      /* perpendicular on the surface */
                      double nn = sqrt(nx * nx + ny * ny + nz * nz);
                      nx /= nn;
                      ny /= nn;
                      nz /= nn;
                      double fx   = (SphP[i].Center[0] - Mesh.VF[vf].cx);
                      double fy   = (SphP[i].Center[1] - Mesh.VF[vf].cy);
                      double fz   = (SphP[i].Center[2] - Mesh.VF[vf].cz);
                      double ff   = (fx * nx + fy * ny + fz * nz);
                      double px   = SphP[i].Center[0] - ff * nx;
                      double py   = SphP[i].Center[1] - ff * ny;
                      double pz   = SphP[i].Center[2] - ff * nz;
                      Mirror[0]   = 2. * px - SphP[i].Center[0];
                      Mirror[1]   = 2. * py - SphP[i].Center[1];
                      Mirror[2]   = 2. * pz - SphP[i].Center[2];
                      CenterOther = Mirror;
                    }
                  else
                    CenterOther = SphP[particle].Center;

#if defined(OUTPUT_DIVVEL) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES) || defined(SF_STELLAR_MASS_TO_GAS_MASS_RATIO)
                  VelOther = P[particle].Vel;
#endif
#ifdef MHD
                  GradOther = &SphP[particle].Grad;
                  BOther    = SphP[particle].B;
#endif
#if defined(EVALPOTENTIAL) && defined(SINK_PARTICLES)
                  PotenOther = P[particle].Potential;
#endif
                }
              else
                {
                  CenterOther = PrimExch[particle].Center;
#if defined(OUTPUT_DIVVEL) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES) || defined(SF_STELLAR_MASS_TO_GAS_MASS_RATIO)
                  VelOther = PrimExch[particle].VelGas;
#endif
#ifdef MHD
                  GradOther = &GradExch[particle];
                  BOther    = PrimExch[particle].B;
#endif
#if defined(EVALPOTENTIAL) && defined(SINK_PARTICLES)
                  PotenOther = PrimExch[particle].Potential;
#endif
                }

#ifdef MHD
              dx = nearest_x(Mesh.VF[vf].cx - CenterOther[0]);
              dy = nearest_y(Mesh.VF[vf].cy - CenterOther[1]);
              dz = nearest_z(Mesh.VF[vf].cz - CenterOther[2]);

              for(int j = 0; j < 3; j++)
                {
                  dB[j] = B[j] - (BOther[j] + GradOther->dB[j][0] * dx + GradOther->dB[j][1] * dy + GradOther->dB[j][2] * dz);
                  B[j]  = 0.5 * (B[j] + BOther[j] + GradOther->dB[j][0] * dx + GradOther->dB[j][1] * dy + GradOther->dB[j][2] * dz);
                }
#endif

              double norm[3];
              norm[0] = nearest_x(Mesh.DP[dp].x - P[i].Pos[0]);
              norm[1] = nearest_y(Mesh.DP[dp].y - P[i].Pos[1]);
              norm[2] = nearest_z(Mesh.DP[dp].z - P[i].Pos[2]);

              double dist = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
              norm[0] /= dist;
              norm[1] /= dist;
              norm[2] /= dist;

#if defined(OUTPUT_DIVVEL) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES)
              double Vel[3];
              for(int j = 0; j < 3; j++)
                Vel[j] = 0.5 * (P[i].Vel[j] + VelOther[j]);
              double nVel = Vel[0] * norm[0] + Vel[1] * norm[1] + Vel[2] * norm[2];
              SphP[i].DivVel += Mesh.VF[vf].area * nVel;
#endif
#ifdef MHD
              double nB = B[0] * norm[0] + B[1] * norm[1] + B[2] * norm[2];
              SphP[i].DivB += Mesh.VF[vf].area * nB;
              double ndB = dB[0] * norm[0] + dB[1] * norm[1] + dB[2] * norm[2];
              SphP[i].DivBalt += Mesh.VF[vf].area * ndB;
#endif
#if defined(EVALPOTENTIAL) && defined(SINK_PARTICLES)
              if(PotenOther < P[i].Potential)
                SphP[i].PotentialPeak = 0;
#endif
#ifdef SF_STELLAR_MASS_TO_GAS_MASS_RATIO
              for(int j = 0; j < 3; j++)
                {
                  SphP[i].VelDisp += (P[i].Vel[j] - VelOther[j]) * (P[i].Vel[j] - VelOther[j]);
                  NNgb++;
                }
#endif
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

#if defined(OUTPUT_DIVVEL) || defined(SGCHEM) || defined(TGCHEM) || defined(CHIMES)
      SphP[i].DivVel /= SphP[i].Volume;
#endif
#ifdef MHD
      SphP[i].DivB /= SphP[i].Volume;
      SphP[i].DivBalt /= SphP[i].Volume;
#endif
#ifdef SF_STELLAR_MASS_TO_GAS_MASS_RATIO
      SphP[i].VelDisp /= 3. * NNgb;
#endif
    }
}
#endif

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
void correct_for_reflective_boundaries(double *ValueOther, double Value, int type, unsigned int *image_flags)
{
#if defined(REFLECTIVE_X)
  if(type == GRADIENT_TYPE_VELX)
    {
      if((*image_flags & REFL_X_FLAGS) && !(*image_flags & OUTFLOW_X))
        *ValueOther *= -1;
      if((*image_flags & REFL_X_FLAGS) && (*image_flags) && (REFLECTIVE_X == 3))
        {
          *ValueOther *= -1;
          *ValueOther *= -1;
        }
      if((*image_flags & REFL_X_FLAGS) && (*image_flags & OUTFLOW_X))
        *ValueOther = Value;
    }

#endif
#if defined(REFLECTIVE_Y)
  if(type == GRADIENT_TYPE_VELY)
    {
      if((*image_flags & REFL_Y_FLAGS) && !(*image_flags & OUTFLOW_Y))
        *ValueOther *= -1;
      if((*image_flags & REFL_Y_FLAGS) && (*image_flags) && (REFLECTIVE_Y == 3))
        {
          *ValueOther *= -1;
          *ValueOther *= -1;
        }
      if((*image_flags & REFL_Y_FLAGS) && (*image_flags & OUTFLOW_Y))
        *ValueOther = Value;
    }

#endif
#if defined(REFLECTIVE_Z)
  if(type == GRADIENT_TYPE_VELZ)
    {
      if((*image_flags & REFL_Z_FLAGS) && !(*image_flags & OUTFLOW_Z))
        *ValueOther *= -1;
      if((*image_flags & REFL_Z_FLAGS) && (*image_flags) && (REFLECTIVE_Z == 3))
        {
          *ValueOther *= -1;
          *ValueOther *= -1;
        }
      if((*image_flags & REFL_Z_FLAGS) && (*image_flags & OUTFLOW_Z))
        *ValueOther = Value;
    }

#endif
}
#endif

/*! \brief Loops through mesh and limits associated gradients.
 *
 *  \return void
 */
void limit_gradients(void)
{
#if defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
  if(!SfVars.limit_gradients)
    {
      mpi_printf("VORONOI: Gradients are not limited!\n");
      return;
    }
#endif

  mpi_printf("VORONOI: Limiting gradients...\n");

  point *DP = Mesh.DP;
  face *VF  = Mesh.VF;

  for(int i = 0; i < Mesh.Nvf; i++)
    {
      if(DP[VF[i].p1].index < 0 || DP[VF[i].p2].index < 0)
        continue;
      for(int j = 0; j < 2; j++)
        {
          point *p;
#ifdef GRADIENT_LIMITER_DUFFELL
          point *pother;
#endif
          if(j == 0)
            {
              p = &DP[VF[i].p1];
#ifdef GRADIENT_LIMITER_DUFFELL
              pother = &DP[VF[i].p2];
#endif
            }
          else
            {
              p = &DP[VF[i].p2];
#ifdef GRADIENT_LIMITER_DUFFELL
              pother = &DP[VF[i].p1];
#endif
            }

          if(p->task == ThisTask && p->index >= 0 && p->index < NumGas)
            {
              int q = p->index;
              if(TimeBinSynchronized[P[q].TimeBinHydro])
                {
                  double d[3];
                  d[0] = nearest_x(VF[i].cx - SphP[q].Center[0]);
                  d[1] = nearest_y(VF[i].cy - SphP[q].Center[1]);
                  d[2] = nearest_z(VF[i].cz - SphP[q].Center[2]);
                  double value;
                  MySingle *data;
                  if(VF[i].area > 1.0e-10 * SphP[q].SurfaceArea)
                    {
                      for(int k = 0; k < N_Grad; k++)
                        {
                          if(grad_elements[k].type == GRADIENT_TYPE_VELX || grad_elements[k].type == GRADIENT_TYPE_VELY ||
                             grad_elements[k].type == GRADIENT_TYPE_VELZ)
                            {
                              value = *(MyFloat *)(((char *)(&P[q])) + grad_elements[k].offset);
                              value /= All.cf_atime;
                            }
                          else
                            value = *(MyFloat *)(((char *)(&SphP[q])) + grad_elements[k].offset);

                          data = (MySingle *)(((char *)(&(SphP[q].Grad))) + grad_elements[k].offset_grad);

#ifndef GRADIENT_LIMITER_DUFFELL
                          if(grad_elements[k].type != GRADIENT_TYPE_RTF)
                            limit_gradient(d, value, minvalues[q * N_Grad + k], maxvalues[q * N_Grad + k], data);
#else
                          double ValueOther;

                          int particle = pother->index;
                          if(particle >= NumGas && pother->task == ThisTask)
                            particle -= NumGas;

                          if(pother->task == ThisTask)
                            {
                              if((grad_elements[k].type == GRADIENT_TYPE_VELX) || (grad_elements[k].type == GRADIENT_TYPE_VELY) ||
                                 (grad_elements[k].type == GRADIENT_TYPE_VELZ))
                                {
                                  ValueOther = *(MyFloat *)(((char *)(&P[particle])) + grad_elements[k].offset);
                                }
                              else
                                ValueOther = *(MyFloat *)(((char *)(&SphP[particle])) + grad_elements[k].offset);
                            }
                          else
                            {
                              ValueOther = *(MyFloat *)(((char *)(&PrimExch[particle])) + grad_elements[k].offset_exch);
                            }

                          limit_gradient(d, value, ValueOther, ValueOther, data);
#endif
                        }
                    }
                }
            }
        }
    }

#ifndef DISABLE_VELOCITY_CSND_SLOPE_LIMITING
  for(int i = 0; i < Mesh.Nvf; i++)
    {
      if(DP[VF[i].p1].index < 0 || DP[VF[i].p2].index < 0)
        continue;
      for(int j = 0; j < 2; j++)
        {
          point *p;

          if(j == 0)
            {
              p = &DP[VF[i].p1];
            }
          else
            {
              p = &DP[VF[i].p2];
            }

          if(p->task == ThisTask && p->index >= 0 && p->index < NumGas)
            {
              int q = p->index;
              if(TimeBinSynchronized[P[q].TimeBinHydro])
                {
                  double d[3];
                  d[0] = nearest_x(VF[i].cx - SphP[q].Center[0]);
                  d[1] = nearest_y(VF[i].cy - SphP[q].Center[1]);
                  d[2] = nearest_z(VF[i].cz - SphP[q].Center[2]);

                  if(VF[i].area > 1.0e-10 * SphP[q].SurfaceArea)
                    {
                      /* let's now limit the overall size of the velocity gradient */
                      MySingle *grad_vx = (MySingle *)(((char *)(&(SphP[q].Grad))) + GVelx->offset_grad);
                      MySingle *grad_vy = (MySingle *)(((char *)(&(SphP[q].Grad))) + GVely->offset_grad);
                      MySingle *grad_vz = (MySingle *)(((char *)(&(SphP[q].Grad))) + GVelz->offset_grad);
                      limit_vel_gradient(d, grad_vx, grad_vy, grad_vz, get_sound_speed(q));
                    }
                }
            }
        }
    }
#endif
}

/*! \brief Limits velocity gradient.
 *
 *  Limit velocity change to the sound speed.
 *
 *  \param[in] d Direction vector.
 *  \param[in, out] grad_vx X-velocity gradient.
 *  \param[in, out] grad_vy Y-velocity gradient.
 *  \param[in, out] grad_vz Z-velocity gradient.
 *  \param[in] csnd sound speed.
 *
 *  \return void
 */
void limit_vel_gradient(double *d, MySingle *grad_vx, MySingle *grad_vy, MySingle *grad_vz, double csnd)
{
#define VEL_GRADIENT_LIMIT_FAC 1.0
  if(All.ComovingIntegrationOn)
    {
      grad_vx[0] -= All.cf_atime * All.cf_Hrate;
      grad_vy[1] -= All.cf_atime * All.cf_Hrate;
      grad_vz[2] -= All.cf_atime * All.cf_Hrate;
    }

  double dvx = fabs(grad_vx[0] * d[0] + grad_vx[1] * d[1] + grad_vx[2] * d[2]);
  double dvy = fabs(grad_vy[0] * d[0] + grad_vy[1] * d[1] + grad_vy[2] * d[2]);
  double dvz = fabs(grad_vz[0] * d[0] + grad_vz[1] * d[1] + grad_vz[2] * d[2]);
  double dv  = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

  if(dv > VEL_GRADIENT_LIMIT_FAC * csnd)
    {
      double fac = VEL_GRADIENT_LIMIT_FAC * csnd / dv;
      for(int i = 0; i < 3; i++)
        {
          grad_vx[i] *= fac;
          grad_vy[i] *= fac;
          grad_vz[i] *= fac;
        }
    }

  if(All.ComovingIntegrationOn)
    {
      grad_vx[0] += All.cf_atime * All.cf_Hrate;
      grad_vy[1] += All.cf_atime * All.cf_Hrate;
      grad_vz[2] += All.cf_atime * All.cf_Hrate;
    }
}

/*! \brief Limits gradients.
 *
 *  Slope limiter.
 *
 *  \param[in] d Direction vector.
 *  \param[in] phi Value.
 *  \param[in] min_phi Lower bound for value+gradient*dx.
 *  \param[in] max_phi Upper bound for value+gradient*dx.
 *  \param[in, out] dphi Gradient.
 *
 *  \return void
 */
void limit_gradient(double *d, double phi, double min_phi, double max_phi, MySingle *dphi)
{
  double dp = dphi[0] * d[0] + dphi[1] * d[1] + dphi[2] * d[2];

  if(dp > 0)
    {
      if(phi + dp > max_phi)
        {
          double fac;

          if(max_phi > phi)
            fac = (max_phi - phi) / dp;
          else
            fac = 0;
          if(fac < 0 || fac > 1)
            terminate("fac=%g\ndp=%g max_phi=%g phi=%g", fac, dp, max_phi, phi);
          dphi[0] *= fac;
          dphi[1] *= fac;
          dphi[2] *= fac;
        }
    }
  else if(dp < 0)
    {
      if(phi + dp < min_phi)
        {
          double fac;

          if(min_phi < phi)
            fac = (min_phi - phi) / dp;
          else
            fac = 0;
          if(fac < 0 || fac > 1)
            terminate("fac=%g\ndp=%g max_phi=%g phi=%g", fac, dp, max_phi, phi);
          dphi[0] *= fac;
          dphi[1] *= fac;
          dphi[2] *= fac;
        }
    }
}

#endif
