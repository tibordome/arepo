/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/forcetree_walk.c
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"
#if defined(RADCOOL) && !defined(GFM)
#error "Need to compile with GFM if RADCOOL option is chosen"
#endif

#ifdef RADCOOL_HOTHALO
static double effectivemass;
#endif

/*! \file forcetree_walk.c
 *  \brief Gravitational tree walk code
 *
 *  This file contains the various gravitational tree walks.
 */

/*! \brief variable for short-range lookup table
 *
 *  contains the factor needed for the short range
 *  contribution of the tree to the gravity force
 */

static float shortrange_table[NTAB + 1];

/*! \brief variable for short-range lookup table
 *
 *  contains the factor needed for the short range
 *  contribution of the tree to the potential energy
 */
static float shortrange_table_potential[NTAB + 1];

/*! \brief Initializes the short range table.
 *
 *  The short range table contains the complementary error function
 *  needed for the computation of the short range part of the gravity
 *  force/potential in case of the TreePM algorithm.
 *
 *  \return void
 */
void force_short_range_init(void)
{
  for(int i = 0; i <= NTAB; i++)
    {
      double u = ((RCUT / 2.0) / NTAB) * i;

      shortrange_table_potential[i] = -erfc(u); /* -r * g(r) */

      if(u > 0)
        shortrange_table[i] = (erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u) - 1.0) / (u * u); /* -g'(r) - 1/r^2 */
      else
        shortrange_table[i] = 0;
    }
}

/*! \brief This routine calculates the (short range) force contribution
 *   for a given particle in case the Tree(PM) algorithm is used.
 *
 *  In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access penalty (which reduces cache performance) incurred by the
 *  table.
 *
 *  Depending on the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 *
 *  \param[in] in Gravdata communicated into function.
 *  \param[in, out] out Gravdata communicated from function.
 *  \param[in] target Index of the particle to be processed.
 *  \param[in] mode 0: process local particle (phase 1), 1: process imported
 *             particle (phase 2).
 *  \param[in] thread_id Id of this thread.
 *  \param[in, out] firstnode First node involved in this algorithm.
 *  \param[in] measure_cost_flag Whether the cost of the tree walk should be
 *             measured.
 *
 *  \return Number of interactions processed for particle i.
 */
int force_treeevaluate(gravdata_in *in, gravdata_out *out, int target, int mode, int thread_id, int numnodes, int *firstnode,
                       int measure_cost_flag)
{
  struct NODE *nop = NULL;
#ifdef MULTIPLE_NODE_SOFTENING
  struct ExtNODE *extnop = 0;
#endif

  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
#ifdef EVALPOTENTIAL
  double pot = 0.0;
#endif

#ifdef TREECOLV2
  double Projection[NPIX]   = {0};
  double ProjectionH2[NPIX] = {0};
  double ProjectionCO[NPIX] = {0};
  double ProjectionC[NPIX]  = {0};
  double tcv2_gasmass, tcv2_h2mass, tcv2_comass, tcv2_cmass, node_size;
#ifdef TREECOLV2_DEBUG
  int idebug_tcv2 = 1;
#else
  int idebug_tcv2 = 0;
#endif
#ifdef TREECOLV2_VEL
  double vx_p = in->Vel[0]; /* velocities of particles */
  double vy_p = in->Vel[1];
  double vz_p = in->Vel[2];
  double vx_n, vy_n, vz_n;
  double v_th2 = in->Vth2;
#endif
#endif

#ifdef RADCOOL
  double youngstellarmass, oldstellarmass;
  double walkBirthTime, r2_ys, r2_os, r2_invfac_ys, r2_invfac_os;
  double dx_ys, dy_ys, dz_ys, dx_os, dy_os, dz_os;
  double walk_Phios, walk_Phins;
  walk_Phios = 0.0;
  walk_Phins = 0.0;
#ifdef RADCOOL_HOTHALO
  double T6gasmass, T7gasmass, T8gasmass;
  double log10temp, r2_T6, r2_T7, r2_T8, r2_invfac_T6, r2_invfac_T7, r2_invfac_T8;
  double dx_T6, dy_T6, dz_T6, dx_T7, dy_T7, dz_T7, dx_T8, dy_T8, dz_T8;
  double walk_PhiT6, walk_PhiT7, walk_PhiT8;
  walk_PhiT6 = 0.0;
  walk_PhiT7 = 0.0;
  walk_PhiT8 = 0.0;
#endif
#endif

#ifdef MODGRAV_EFF_MASS
  double modgrav_acc_x = 0;
  double modgrav_acc_y = 0;
  double modgrav_acc_z = 0;
#endif

#ifdef PE_MCS
  double lum_FUV;
  double r2_FUV, r2_invfac_FUV, dx_FUV, dy_FUV, dz_FUV;
  double walk_G_FUV = 0.0;
#endif

#ifdef HII_MCS_LR
  double lum_Hii;
  double r2_Hii, r2_invfac_Hii, dx_Hii, dy_Hii, dz_Hii;
  double walk_e_Hii = 0.0;
#endif

  int ninteractions = 0;

  double pos_x = in->Pos[0];
  double pos_y = in->Pos[1];
  double pos_z = in->Pos[2];
  double aold  = All.ErrTolForceAcc * in->OldAcc;
  double h_i   = All.ForceSoftening[in->SofteningType];

#ifdef PMGRID
  double rcut  = All.Rcut[0];
  double asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(pmforce_is_particle_high_res(in->Type, in->Pos))
    {
      rcut  = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif

  double rcut2     = rcut * rcut;
  double asmthinv  = 0.5 / asmth;
  double asmthinv2 = asmthinv * asmthinv;
  double asmthfac  = asmthinv * (NTAB / (RCUT / 2.0));
#endif

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == 0)
        no = Tree_MaxPart; /* root node */
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          double dx, dy, dz, r2, mass, hmax;

#ifdef MULTIPLE_NODE_SOFTENING
          int indi_flag1 = -1, indi_flag2 = 0;
#endif

          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              mass = P[no].Mass;

#ifdef TREECOLV2
              tcv2_gasmass = 0;
              tcv2_h2mass  = 0;
              tcv2_comass  = 0;
              tcv2_cmass   = 0;
              if(P[no].Type == 0)
                {
                  tcv2_gasmass = P[no].Mass;
#ifdef TREECOLV2_H2
                  tcv2_h2mass = P[no].Mass * 2.0 * SphP[no].TracAbund[IH2] / (1.0 + 4.0 * ABHE);
#endif
#ifdef TREECOLV2_CO
                  tcv2_comass = P[no].Mass * 28.0 * SphP[no].TracAbund[ICO] / (1.0 + 4.0 * ABHE);
#endif
#ifdef TREECOLV2_C
                  tcv2_cmass = P[no].Mass * 12.0 * SphP[no].TracAbund[ICATOM] / (1.0 + 4.0 * ABHE);
#endif
                  node_size = pow(SphP[no].Volume, 1. / 3.); /* Will be converted to angle later */
#ifdef TREECOLV2_NO_GAS_SELFGRAVITY
                  mass = 0.; /* point mass set to zero so that accel = 0*/
#endif
                }
#ifdef TREECOLV2_VEL
              vx_n = P[no].Vel[0];
              vy_n = P[no].Vel[1];
              vz_n = P[no].Vel[2];
#endif
#endif

#ifdef RADCOOL
              if(P[no].Type == 4)
                {
                  r2_ys         = r2;
                  r2_os         = r2;
                  walkBirthTime = get_time_difference_in_Gyr(STP(no).BirthTime, All.Time);
                  if(walkBirthTime <= TIMEON_NEWSTARS)
                    {
                      youngstellarmass = P[no].Mass;
                      oldstellarmass   = 0.0;
                    }
                  else if(walkBirthTime >= TIMEON_OLDSTARS)
                    {
                      oldstellarmass   = P[no].Mass;
                      youngstellarmass = 0.0;
                    }
                  else
                    {
                      oldstellarmass   = 0.0;
                      youngstellarmass = 0.0;
                    }
                }
              else
                {
                  oldstellarmass   = 0.0;
                  youngstellarmass = 0.0;
                  r2_ys            = r2;
                  r2_os            = r2;
                }
#ifdef RADCOOL_HOTHALO
              if(P[no].Type == 0)
                {
                  r2_T6     = r2;
                  r2_T7     = r2;
                  r2_T8     = r2;
                  log10temp = log10(calculate_HH_temperature(SphP[no].Utherm
#ifdef COOLING
                                                             ,
                                                             SphP[no].Ne
#endif
                                                             ));
                  if((log10temp > TLOGMIN6) && (log10temp <= TLOGMAX6))
                    {
                      T6gasmass = P[no].Mass * SphP[no].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                                  * (1.0 + T6METALBOOSTFACTOR * SphP[no].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN7) && (log10temp <= TLOGMAX7))
                    {
                      T7gasmass = P[no].Mass * SphP[no].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                                  * (1.0 + T7METALBOOSTFACTOR * SphP[no].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T8gasmass = 0.0;
                      T6gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN8) && (log10temp <= TLOGMAX8))
                    {
                      T8gasmass = P[no].Mass * SphP[no].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                                  * (1.0 + T8METALBOOSTFACTOR * SphP[no].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                    }
                  else
                    {
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                }
              else
                {
                  T6gasmass = 0.0;
                  T7gasmass = 0.0;
                  T8gasmass = 0.0;
                  r2_T6     = r2;
                  r2_T7     = r2;
                  r2_T8     = r2;
                }
#endif
#endif

#ifdef PE_MCS
              if(P[no].Type == 4)
                {
                  r2_FUV  = r2;
                  lum_FUV = STP(no).L_FUV;
                }
              else
                {
                  r2_FUV  = r2;
                  lum_FUV = 0.0;
                }
#endif

#ifdef HII_MCS_LR
              if(P[no].Type == 0)
                {
                  r2_Hii  = r2;
                  lum_Hii = SphP[no].L_Hii;
                }
              else
                {
                  r2_Hii  = r2;
                  lum_Hii = 0.0;
                }
#endif

              if(measure_cost_flag)
                Thread[thread_id].P_CostCount[no]++;

              double h_j = All.ForceSoftening[P[no].SofteningType];

              hmax = (h_j > h_i) ? h_j : h_i;

              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* we have an  internal node */
            {
              if(mode == 1)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              nop = &Nodes[no];

#ifdef MODGRAV_EFF_MASS
              if(All.PM_Ti_endstep == All.Ti_Current && nop->mg_bitflag.amr_node == 1)
                {
                  // nop->GravCost += 1;
                  modgrav_tree_fifth_force(pos_x, pos_y, pos_z, no, asmthfac, h_i, shortrange_table, &modgrav_acc_x, &modgrav_acc_y,
                                           &modgrav_acc_z
#ifdef EVALPOTENTIAL
                                           ,
                                           &pot
#endif
                  );
                }
#endif

              mass = nop->u.d.mass;
              dx   = GRAVITY_NEAREST_X(nop->u.d.s[0] - pos_x);
              dy   = GRAVITY_NEAREST_Y(nop->u.d.s[1] - pos_y);
              dz   = GRAVITY_NEAREST_Z(nop->u.d.s[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

#if defined(PMGRID) && !defined(TREECOLV2)
              if(r2 > rcut2)
                {
                  /* check whether we can stop walking along this branch */
                  double eff_dist = rcut + 0.5 * nop->len;

                  double dist = GRAVITY_NEAREST_X(nop->center[0] - pos_x);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }

                  dist = GRAVITY_NEAREST_Y(nop->center[1] - pos_y);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }

                  dist = GRAVITY_NEAREST_Z(nop->center[2] - pos_z);
                  if(dist < -eff_dist || dist > eff_dist)
                    {
                      no = nop->u.d.sibling;
                      continue;
                    }
                }
#endif

#ifdef TREECOLV2
#ifdef TREECOLV2_VEL
              vx_n = nop->Vel[0];
              vy_n = nop->Vel[1];
              vz_n = nop->Vel[2];
#endif
              tcv2_gasmass = nop->u.d.gasmass;
              if(isnan(tcv2_gasmass))
                printf("Node mass is NaN %d\n", no);
#ifdef TREECOLV2_H2
              tcv2_h2mass = nop->u.d.h2mass;
#endif
#ifdef TREECOLV2_CO
              tcv2_comass = nop->u.d.comass;
#endif
#ifdef TREECOLV2_C
              tcv2_cmass = nop->u.d.cmass;
#endif
              node_size = nop->len;
#ifdef TREECOLV2_NO_GAS_SELFGRAVITY
              mass -= nop->u.d.gasmass;
#endif
#endif

#ifdef RADCOOL

              dx_ys = GRAVITY_NEAREST_X(nop->young_stellar_s[0] - pos_x);
              dy_ys = GRAVITY_NEAREST_X(nop->young_stellar_s[1] - pos_y);
              dz_ys = GRAVITY_NEAREST_X(nop->young_stellar_s[2] - pos_z);

              dx_os = GRAVITY_NEAREST_X(nop->old_stellar_s[0] - pos_x);
              dy_os = GRAVITY_NEAREST_X(nop->old_stellar_s[1] - pos_y);
              dz_os = GRAVITY_NEAREST_X(nop->old_stellar_s[2] - pos_z);

              r2_ys = dx_ys * dx_ys + dy_ys * dy_ys + dz_ys * dz_ys;
              r2_os = dx_os * dx_os + dy_os * dy_os + dz_os * dz_os;

              youngstellarmass = nop->young_stellar_mass;
              oldstellarmass   = nop->old_stellar_mass;
#ifdef RADCOOL_HOTHALO
              dx_T6 = GRAVITY_NEAREST_X(nop->T6_gas_s[0] - pos_x);
              dy_T6 = GRAVITY_NEAREST_Y(nop->T6_gas_s[1] - pos_y);
              dz_T6 = GRAVITY_NEAREST_Z(nop->T6_gas_s[2] - pos_z);

              dx_T7 = GRAVITY_NEAREST_X(nop->T7_gas_s[0] - pos_x);
              dy_T7 = GRAVITY_NEAREST_Y(nop->T7_gas_s[1] - pos_y);
              dz_T7 = GRAVITY_NEAREST_Z(nop->T7_gas_s[2] - pos_z);

              dx_T8 = GRAVITY_NEAREST_X(nop->T8_gas_s[0] - pos_x);
              dy_T8 = GRAVITY_NEAREST_Y(nop->T8_gas_s[1] - pos_y);
              dz_T8 = GRAVITY_NEAREST_Z(nop->T8_gas_s[2] - pos_z);

              r2_T6 = dx_T6 * dx_T6 + dy_T6 * dy_T6 + dz_T6 * dz_T6;
              r2_T7 = dx_T7 * dx_T7 + dy_T7 * dy_T7 + dz_T7 * dz_T7;
              r2_T8 = dx_T8 * dx_T8 + dy_T8 * dy_T8 + dz_T8 * dz_T8;

              T6gasmass = nop->T6_gas_mass;
              T7gasmass = nop->T7_gas_mass;
              T8gasmass = nop->T8_gas_mass;
#endif
#endif

#ifdef PE_MCS
              dx_FUV = GRAVITY_NEAREST_X(nop->lum_FUV_s[0] - pos_x);
              dy_FUV = GRAVITY_NEAREST_X(nop->lum_FUV_s[1] - pos_y);
              dz_FUV = GRAVITY_NEAREST_X(nop->lum_FUV_s[2] - pos_z);

              r2_FUV = dx_FUV * dx_FUV + dy_FUV * dy_FUV + dz_FUV * dz_FUV;

              lum_FUV = nop->lum_FUV;
#endif

#ifdef HII_MCS_LR
              dx_Hii = GRAVITY_NEAREST_X(nop->lum_Hii_s[0] - pos_x);
              dy_Hii = GRAVITY_NEAREST_X(nop->lum_Hii_s[1] - pos_y);
              dz_Hii = GRAVITY_NEAREST_X(nop->lum_Hii_s[2] - pos_z);

              r2_Hii = dx_Hii * dx_Hii + dy_Hii * dy_Hii + dz_Hii * dz_Hii;

              lum_Hii = nop->lum_Hii;
#endif

              if(All.ErrTolTheta) /* check Barnes-Hut opening criterion */
                {
                  if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                }
              else /* check relative opening criterion */
                {
                  double len2 = nop->len * nop->len;

                  if(len2 > r2 * (1.2 * 1.2)) /* add a worst case protection */
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }

#ifdef ACTIVATE_MINIMUM_OPENING_ANGLE
                  if(mass * len2 > r2 * r2 * aold && len2 > r2 * (0.4 * 0.4))
#else
                  if(mass * len2 > r2 * r2 * aold)
#endif
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }

                  /* check in addition whether we lie inside the cell */

                  if(fabs(GRAVITY_NEAREST_X(nop->center[0] - pos_x)) < 0.60 * nop->len)
                    {
                      if(fabs(GRAVITY_NEAREST_Y(nop->center[1] - pos_y)) < 0.60 * nop->len)
                        {
                          if(fabs(GRAVITY_NEAREST_Z(nop->center[2] - pos_z)) < 0.60 * nop->len)
                            {
                              no = nop->u.d.nextnode;
                              continue;
                            }
                        }
                    }
                }

              double h_j = All.ForceSoftening[nop->u.d.maxsofttype];

              if(h_j > h_i)
                {
#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(nop->u.d.maxhydrosofttype != nop->u.d.minhydrosofttype)
                    if(ExtNodes[no].mass_per_type[0] > 0)
                      if(r2 < All.ForceSoftening[nop->u.d.maxhydrosofttype] * All.ForceSoftening[nop->u.d.maxhydrosofttype])
                        {
                          /* open cell */
                          no = nop->u.d.nextnode;
                          continue;
                        }
#endif
                  indi_flag1 = 0;
                  indi_flag2 = NSOFTTYPES;
#else
                  if(r2 < h_j * h_j)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
#endif
                  hmax = h_j;
                }
              else
                hmax = h_i;

                /* ok, node can be used */
#ifdef MULTIPLE_NODE_SOFTENING
              extnop = &ExtNodes[no];
#endif
              if(measure_cost_flag && mass)
                Thread[thread_id].Node_CostCount[no]++;

#ifdef MODGRAV_EFF_MASS
              if(All.PM_Ti_endstep == All.Ti_Current && nop->mg_bitflag.amr_node == 0)
                {
                  // nop->GravCost += 1;
                  modgrav_tree_fifth_force(pos_x, pos_y, pos_z, no, asmthfac, h_i, shortrange_table, &modgrav_acc_x, &modgrav_acc_y,
                                           &modgrav_acc_z
#ifdef EVALPOTENTIAL
                                           ,
                                           &pot
#endif
                  );
                }
#endif

              no = nop->u.d.sibling;
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              mass = Tree_Points[n].Mass;

              if(measure_cost_flag)
                Thread[thread_id].TreePoints_CostCount[n]++;

#ifdef TREECOLV2
              tcv2_gasmass = 0;
              tcv2_h2mass  = 0;
              tcv2_comass  = 0;
              tcv2_cmass   = 0;
              if(Tree_Points[n].Type == 0)
                {
                  tcv2_gasmass = mass;
#ifdef TREECOLV2_H2
                  tcv2_h2mass = mass * Tree_Points[n].TracAbund[IH2];
#endif
#ifdef TREECOLV2_CO
                  tcv2_comass = mass * Tree_Points[n].TracAbund[ICO];
#endif
#ifdef TREECOLV2_C
                  tcv2_cmass = mass * Tree_Points[n].TracAbund[ICATOM];
#endif
                  node_size = 2 * Tree_Points[n].Hsml; /* Will be converted to angle later */
#ifdef TREECOLV2_NO_GAS_SELFGRAVITY
                  mass = 0.; /* point mass set to zero so that accel = 0*/
#endif
                }
#ifdef TREECOLV2_VEL
              vx_n = Tree_Points[n].Vel[0];
              vy_n = Tree_Points[n].Vel[1];
              vz_n = Tree_Points[n].Vel[2];
#endif
#endif

#ifdef RADCOOL
              if((Tree_Points[n].Type) == 4)
                {
                  r2_ys         = r2;
                  r2_os         = r2;
                  walkBirthTime = get_time_difference_in_Gyr(Tree_Points[n].BirthTime, All.Time);
                  if(walkBirthTime <= TIMEON_NEWSTARS)
                    {
                      youngstellarmass = Tree_Points[n].Mass;
                      oldstellarmass   = 0.0;
                    }
                  else if(walkBirthTime >= TIMEON_OLDSTARS)
                    {
                      oldstellarmass   = Tree_Points[n].Mass;
                      youngstellarmass = 0.0;
                    }
                  else
                    {
                      youngstellarmass = 0.0;
                      oldstellarmass   = 0.0;
                    }
                }
              else
                {
                  youngstellarmass = 0.0;
                  oldstellarmass   = 0.0;
                  r2_ys            = r2;
                  r2_os            = r2;
                }
#ifdef RADCOOL_HOTHALO
              if((Tree_Points[n].Type) == 0)
                {
                  r2_T6     = r2;
                  r2_T7     = r2;
                  r2_T8     = r2;
                  log10temp = log10(calculate_HH_temperature(Tree_Points[n].Utherm
#ifdef COOLING
                                                             ,
                                                             Tree_Points[n].Ne
#endif
                                                             ));
                  if((log10temp > TLOGMIN6) && (log10temp <= TLOGMAX6))
                    {
                      T6gasmass = Tree_Points[n].Mass * Tree_Points[n].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                                  * (1.0 + T6METALBOOSTFACTOR * Tree_Points[n].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN7) && (log10temp <= TLOGMAX7))
                    {
                      T7gasmass = Tree_Points[n].Mass * Tree_Points[n].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                                  * (1.0 + T7METALBOOSTFACTOR * Tree_Points[n].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T8gasmass = 0.0;
                      T6gasmass = 0.0;
                    }
                  else if((log10temp > TLOGMIN8) && (log10temp <= TLOGMAX8))
                    {
                      T8gasmass = Tree_Points[n].Mass * Tree_Points[n].Density
#ifdef RADCOOL_HOTHALO_METAl_BOOST
                                  * (1.0 + T8METALBOOSTFACTOR * Tree_Points[n].Metallicity / GFM_SOLAR_METALLICITY)
#endif
                          ;
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                    }
                  else
                    {
                      T6gasmass = 0.0;
                      T7gasmass = 0.0;
                      T8gasmass = 0.0;
                    }
                }
              else
                {
                  T6gasmass = 0.0;
                  T7gasmass = 0.0;
                  T8gasmass = 0.0;
                  r2_T6     = r2;
                  r2_T7     = r2;
                  r2_T8     = r2;
                }
#endif
#endif

#ifdef PE_MCS
              if((Tree_Points[n].Type) == 4)
                {
                  r2_FUV  = r2;
                  lum_FUV = Tree_Points[n].lum_FUV;
                }
              else
                {
                  r2_FUV  = r2;
                  lum_FUV = 0.0;
                }
#endif

#ifdef HII_MCS_LR
              if((Tree_Points[n].Type) == 0)
                {
                  r2_Hii  = r2;
                  lum_Hii = Tree_Points[n].lum_Hii;
                }
              else
                {
                  r2_Hii  = r2;
                  lum_Hii = 0.0;
                }
#endif

              double h_j = All.ForceSoftening[Tree_Points[n].SofteningType];

              hmax = (h_j > h_i) ? h_j : h_i;

              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
            {
              if(mode == 0)
                {
                  tree_treefind_export_node_threads(no, target, thread_id);
                }

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }

#ifdef TREECOLV2
          /* Do TreeCol projection */
          treecolv2_add_node_contribution(dx, dy, dz
#ifdef TREECOLV2_VEL
                                          ,
                                          vx_p, vy_p, vz_p, vx_n, vy_n, vz_n, v_th2
#endif
                                          ,
                                          node_size, tcv2_gasmass, Projection
#ifdef TREECOLV2_H2
                                          ,
                                          tcv2_h2mass, ProjectionH2
#endif
#ifdef TREECOLV2_CO
                                          ,
                                          tcv2_comass, ProjectionCO
#endif
#ifdef TREECOLV2_C
                                          ,
                                          tcv2_cmass, ProjectionC
#endif
                                          ,
                                          idebug_tcv2);
#endif

          /* now evaluate the multipole moment */
          if(mass)
            {
              double r = sqrt(r2);

#ifdef PMGRID
              double tabentry = asmthfac * r;
              int tabindex    = (int)tabentry;

              if(tabindex < NTAB)
                {
                  double tabweight    = tabentry - tabindex;
                  double factor_force = (1.0 - tabweight) * shortrange_table[tabindex] + tabweight * shortrange_table[tabindex + 1];
#ifdef EVALPOTENTIAL
                  double factor_pot =
                      (1.0 - tabweight) * shortrange_table_potential[tabindex] + tabweight * shortrange_table_potential[tabindex + 1];
#endif
#endif

#ifdef MULTIPLE_NODE_SOFTENING
                  for(int type = indi_flag1; type < indi_flag2; type++)
                    {
                      if(type >= 0)
                        {
                          mass = extnop->mass_per_type[type];
                          double h_j;
#ifdef ADAPTIVE_HYDRO_SOFTENING
                          if(type == 0)
                            h_j = All.ForceSoftening[nop->u.d.maxhydrosofttype];
                          else
#endif
                            h_j = All.ForceSoftening[type];

                          hmax = (h_j > h_i) ? h_j : h_i;
                        }

                      if(mass)
                        {
#endif
                          double fac;
#ifdef EVALPOTENTIAL
                          double wp;
#endif

                          if(r >= hmax)
                            {
                              double rinv  = 1.0 / r;
                              double rinv3 = rinv * rinv * rinv;
#ifdef PMGRID
                              fac = rinv3 + rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                              wp = rinv * factor_pot; /* wp   = -g(r)    */
#endif
#else
                  fac = rinv3;
#ifdef EVALPOTENTIAL
                  wp  = -rinv;
#endif
#endif
                            }
                          else
                            {
                              double h_inv  = 1.0 / hmax;
                              double h3_inv = h_inv * h_inv * h_inv;
                              double u      = r * h_inv;

                              if(u < 0.5)
                                {
                                  double u2 = u * u;
                                  fac       = h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
#ifdef EVALPOTENTIAL
                                  wp = h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
#endif
                                }
                              else
                                {
                                  double u2 = u * u;
                                  double u3 = u2 * u;
                                  fac       = h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
#ifdef EVALPOTENTIAL
                                  wp = h_inv * (SOFTFAC13 + SOFTFAC14 / u +
                                                u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
#endif
                                }

#ifdef PMGRID
                              if(r > 0)
                                {
                                  double rinv = 1.0 / r;
                                  fac += rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                                  wp += rinv * (factor_pot + 1.0); /* wp   = -g(r)    */
#endif
                                }
#endif
                            }

#ifdef RADCOOL
#ifdef MULTIPLE_NODE_SOFTENING
                          if(type == indi_flag1)
#endif
                            {
                              r2_invfac_ys = 1.0 / (r2_ys + hmax * hmax);
                              r2_invfac_os = 1.0 / (r2_os + hmax * hmax);
                              walk_Phins += youngstellarmass * r2_invfac_ys;
                              walk_Phios += oldstellarmass * r2_invfac_os;
#ifdef RADCOOL_HOTHALO
                              r2_invfac_T6 = 1.0 / (r2_T6 + hmax * hmax);
                              r2_invfac_T7 = 1.0 / (r2_T7 + hmax * hmax);
                              r2_invfac_T8 = 1.0 / (r2_T8 + hmax * hmax);
                              walk_PhiT6 += T6gasmass * r2_invfac_T6;
                              walk_PhiT7 += T6gasmass * r2_invfac_T7;
                              walk_PhiT8 += T6gasmass * r2_invfac_T8;
#endif
                            }
#endif

#ifdef PE_MCS
#ifdef MULTIPLE_NODE_SOFTENING
                          if(type == indi_flag1)
#endif
                            {
                              r2_invfac_FUV = 1.0 / (r2_FUV + hmax * hmax);
                              walk_G_FUV += lum_FUV * r2_invfac_FUV;
                            }
#endif

#ifdef HII_MCS_LR
#ifdef MULTIPLE_NODE_SOFTENING
                          if(type == indi_flag1)
#endif
                            {
                              r2_invfac_Hii = 1.0 / (r2_Hii + hmax * hmax);
                              walk_e_Hii += lum_Hii * r2_invfac_Hii;
                            }
#endif

#ifdef EVALPOTENTIAL
                          pot += mass * wp;
#endif
                          fac *= mass;

                          acc_x += dx * fac;
                          acc_y += dy * fac;
                          acc_z += dz * fac;

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
                          double fcorr[3];
                          ewald_corr(dx, dy, dz, fcorr);
                          acc_x += mass * fcorr[0];
                          acc_y += mass * fcorr[1];
                          acc_z += mass * fcorr[2];
#ifdef EVALPOTENTIAL
                          pot += mass * ewald_pot_corr(dx, dy, dz);
#endif
#endif

#ifdef MULTIPLE_NODE_SOFTENING
                        }
                    }
#endif
                  ninteractions++;
#ifdef PMGRID
                }
#endif
            }
        }
    }

  out->Acc[0] = acc_x;
  out->Acc[1] = acc_y;
  out->Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
  out->Potential = pot;
#endif
#ifdef NO_GRAVITY_TYPE
  if(in->Type == NO_GRAVITY_TYPE)
    {
      out->Acc[0] = 0.0;
      out->Acc[1] = 0.0;
      out->Acc[2] = 0.0;
#ifdef EVALPOTENTIAL
      out->Potential = 0.0;
#endif
    }
#endif
#ifdef OUTPUTGRAVINTERACTIONS
  out->GravInteractions = ninteractions;
#endif
#ifdef RADCOOL
  out->Phios = walk_Phios;
  out->Phins = walk_Phins;
#ifdef RADCOOL_HOTHALO
  out->PhiT6 = walk_PhiT6;
  out->PhiT7 = walk_PhiT7;
  out->PhiT8 = walk_PhiT8;
#endif
#endif
#ifdef TREECOLV2
  for(int j = 0; j < NPIX; j++)
    {
      out->Projection[j] = Projection[j];
#ifdef TREECOLV2_H2
      out->ProjectionH2[j] = ProjectionH2[j];
#endif
#ifdef TREECOLV2_CO
      out->ProjectionCO[j] = ProjectionCO[j];
#endif
#ifdef TREECOLV2_C
      out->ProjectionC[j] = ProjectionC[j];
#endif
    }
#endif
#ifdef MODGRAV
  out->ModgravAcc[0] = modgrav_acc_x;
  out->ModgravAcc[1] = modgrav_acc_y;
  out->ModgravAcc[2] = modgrav_acc_z;
#endif

#ifdef PE_MCS
  out->G_FUV = walk_G_FUV * All.Factor_FUV * All.cf_a2inv; /* Now energy density in Habing Units */
#endif

#ifdef HII_MCS_LR
  out->EnergyDensHii = walk_e_Hii * All.Factor_Hii * All.cf_a2inv; /* Now energy density in Habing Units */
#endif

  return ninteractions;
}

/*! \brief Prepares node to be exported.
 *
 *  \param[in] no Index of node.
 *  \param[in] i Index of particle.
 *  \param[in] thread_id ID of thread.
 *
 *  \return 0
 */
int tree_treefind_export_node_threads(int no, int i, int thread_id)
{
  /* The task indicated by the pseudoparticle node */
  int task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != i)
    {
      Thread[thread_id].Exportflag[task]     = i;
      int nexp                               = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = i;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp                      = Thread[thread_id].NexportNodes++;
  nexp                          = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task           = task;
  nodelist[nexp].Index          = i;
  nodelist[nexp].Node           = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}

#ifdef ALLOW_DIRECT_SUMMATION
/*! \brief Kernel of direct summation force calculation.
 *
 *  \param[in] target Index of particle in import array.
 *  \param[in] result_idx Index in result array.
 *  \param[in] nimport number of imported particles.
 *
 *  \return void
 */
void force_evaluate_direct(int target, int result_idx, int nimport)
{
  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;
#ifdef EVALPOTENTIAL
  double pot = 0.0;
#endif

  double pos_x = DirectDataAll[target].Pos[0];
  double pos_y = DirectDataAll[target].Pos[1];
  double pos_z = DirectDataAll[target].Pos[2];
  double h_i   = All.ForceSoftening[DirectDataAll[target].SofteningType];

#ifdef PMGRID
  double asmth = All.Asmth[0];
#if defined(PLACEHIGHRESREGION)
  int ptype_i = DirectDataAll[target].Type;
  if(pmforce_is_particle_high_res(ptype_i, DirectDataAll[target].Pos))
    asmth = All.Asmth[1];
#endif
  double asmthinv  = 0.5 / asmth;
  double asmthinv2 = asmthinv * asmthinv;
  double asmthfac  = asmthinv * (NTAB / (RCUT / 2.0));
#endif

  for(int j = 0; j < nimport; j++)
    {
      double h_j = All.ForceSoftening[DirectDataAll[j].SofteningType];

      double hmax = (h_j > h_i) ? h_j : h_i;

      double dx = GRAVITY_NEAREST_X(DirectDataAll[j].Pos[0] - pos_x);
      double dy = GRAVITY_NEAREST_Y(DirectDataAll[j].Pos[1] - pos_y);
      double dz = GRAVITY_NEAREST_Z(DirectDataAll[j].Pos[2] - pos_z);

      double r2 = dx * dx + dy * dy + dz * dz;

      double mass = DirectDataAll[j].Mass;

      /* now evaluate the force component */

      double r = sqrt(r2);

#ifdef PMGRID
      double tabentry = asmthfac * r;
      int tabindex    = (int)tabentry;

      if(tabindex < NTAB)
        {
          double tabweight    = tabentry - tabindex;
          double factor_force = (1.0 - tabweight) * shortrange_table[tabindex] + tabweight * shortrange_table[tabindex + 1];
#ifdef EVALPOTENTIAL
          double factor_pot =
              (1.0 - tabweight) * shortrange_table_potential[tabindex] + tabweight * shortrange_table_potential[tabindex + 1];
#endif
#endif

          double fac;
#ifdef EVALPOTENTIAL
          double wp;
#endif

          if(r >= hmax)
            {
              double rinv  = 1.0 / r;
              double rinv3 = rinv * rinv * rinv;
#ifdef PMGRID
              fac = rinv3 + rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
              wp = rinv * factor_pot; /* wp   = -g(r)    */
#endif
#else
          fac = rinv3;
#ifdef EVALPOTENTIAL
          wp  = -rinv;
#endif
#endif
            }
          else
            {
              double h_inv  = 1.0 / hmax;
              double h3_inv = h_inv * h_inv * h_inv;
              double u      = r * h_inv;

              if(u < 0.5)
                {
                  double u2 = u * u;
                  fac       = h3_inv * (SOFTFAC1 + u2 * (SOFTFAC2 * u + SOFTFAC3));
#ifdef EVALPOTENTIAL
                  wp = h_inv * (SOFTFAC4 + u2 * (SOFTFAC5 + u2 * (SOFTFAC6 * u + SOFTFAC7)));
#endif
                }
              else
                {
                  double u2 = u * u;
                  double u3 = u2 * u;
                  fac       = h3_inv * (SOFTFAC8 + SOFTFAC9 * u + SOFTFAC10 * u2 + SOFTFAC11 * u3 + SOFTFAC12 / u3);
#ifdef EVALPOTENTIAL
                  wp = h_inv * (SOFTFAC13 + SOFTFAC14 / u + u2 * (SOFTFAC1 + u * (SOFTFAC15 + u * (SOFTFAC16 + SOFTFAC17 * u))));
#endif
                }
#ifdef PMGRID
              if(r > 0)
                {
                  double rinv = 1.0 / r;
                  fac += rinv * factor_force * asmthinv2; /* fac  = -g'(r)/r */
#ifdef EVALPOTENTIAL
                  wp += rinv * (factor_pot + 1.0); /* wp   = -g(r)    */
#endif
                }
#endif
            }

#ifdef EVALPOTENTIAL
          pot += mass * wp;
#endif
          fac *= mass;

          acc_x += dx * fac;
          acc_y += dy * fac;
          acc_z += dz * fac;

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && !defined(ONEDIMS_SPHERICAL)
          {
            double fcorr[3];
            ewald_corr(dx, dy, dz, fcorr);
            acc_x += mass * fcorr[0];
            acc_y += mass * fcorr[1];
            acc_z += mass * fcorr[2];
#ifdef EVALPOTENTIAL
            pot += mass * ewald_pot_corr(dx, dy, dz);
#endif
          }
#endif

#ifdef PMGRID
        }
#endif
    }

  DirectAccOut[result_idx].Acc[0] = acc_x;
  DirectAccOut[result_idx].Acc[1] = acc_y;
  DirectAccOut[result_idx].Acc[2] = acc_z;
#ifdef EVALPOTENTIAL
  DirectAccOut[result_idx].Potential = pot;
#endif
}
#endif
