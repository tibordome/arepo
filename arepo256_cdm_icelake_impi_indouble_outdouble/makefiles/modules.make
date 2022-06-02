# AREPO Makefile
#   see documentation/getting_started.md
#
# Add source and header files for new modules below.
# Add libraries for new modules in makefiles/modules-lib.make.

ifeq (VORONOI, $(findstring VORONOI, $(CONFIGVARS)))
ifeq (TWODIMS, $(findstring TWODIMS, $(CONFIGVARS)))
  OBJS    += voronoi_2d.o
endif
endif

ifeq (NEW_FFT, $(findstring NEW_FFT, $(CONFIGVARS)))
  OBJS    += my_fft/my_fft.o
  SUBDIRS += my_fft
endif

ifeq (MYIBARRIER, $(findstring MYIBARRIER, $(CONFIGVARS)))
  OBJS    += mpi_utils/myIBarrier.o
  INCL    += mpi_utils/myIBarrier.h
endif

ifeq (MHD, $(findstring MHD, $(CONFIGVARS)))
  OBJS    += mhd.o
endif

ifeq (ADDBACKGROUNDGRID, $(findstring ADDBACKGROUNDGRID, $(CONFIGVARS)))
  OBJS    += add_backgroundgrid/add_bggrid.o add_backgroundgrid/calc_weights.o add_backgroundgrid/distribute.o
  INCL    += add_backgroundgrid/add_bggrid.h
  SUBDIRS += add_backgroundgrid
endif

ifeq (DVR_RENDER, $(findstring DVR_RENDER, $(CONFIGVARS)))
  OBJS    += dvr_render/dvr_render.o
  INCL    += dvr_render/dvr_render.h
  SUBDIRS += dvr_render
endif

ifeq (COOLING, $(findstring COOLING, $(CONFIGVARS)))
  OBJS    += cooling/cooling.o cooling/simple_cooling.o
  INCL    += cooling/cooling_vars.h cooling/cooling_proto.h
  SUBDIRS += cooling
endif

ifeq (RT_ENABLE, $(findstring RT_ENABLE, $(CONFIGVARS)))
  OBJS    += rt/rt_CGmethod.o rt/rt_gradients.o rt/rt_voronoi.o  \
             rt/rt_stellar_sources.o rt/rt_advect.o rt/rt_cooling.o rt/rt_chem.o \
             rt/rt_inject_photons_sfr.o rt/rt_inject_photons.o\
             rt/rt_optical_depth.o rt/pix2vec_ring.o
  SUBDIRS += rt
endif

ifeq (BLACK_HOLES, $(findstring BLACK_HOLES, $(CONFIGVARS)))
  OBJS    += blackhole/blackhole_neighbors.o blackhole/blackhole_mergers.o blackhole/blackhole_bubbles.o blackhole/blackhole_adios_wind.o \
             blackhole/blackhole_swallowgas.o blackhole/blackhole.o blackhole/blackhole_density.o blackhole/blackhole_adios_wind_randomized.o \
             blackhole/blackhole_disk_vorticity.o blackhole/blackhole_bubbles_nf.o blackhole/blackhole_friction.o \
             blackhole/blackhole_refinement.o blackhole/blackhole_mdot.o  blackhole/blackhole_feedback.o  blackhole/blackhole_centering.o \
             blackhole/blackhole_bipolar.o blackhole/blackhole_spin.o
  INCL    += blackhole/blackhole_proto.h
  SUBDIRS += blackhole
endif

ifeq (FOF, $(findstring FOF, $(CONFIGVARS)))
  OBJS    += fof/fof.o fof/fof_vars.o fof/fof_distribute.o fof/fof_findgroups.o fof/fof_nearest.o fof/fof_io.o fof/fof_sort_kernels.o \
             fof/fof_fuzz.o fof/fof_gfm.o fof/fof_bh.o fof/fof_spinmeasurement.o fof/fof_massiveseeds.o
  INCL    += fof/fof.h
  SUBDIRS += fof
endif

ifeq (GFM_AGN_RADIATION, $(findstring GFM_AGN_RADIATION, $(CONFIGVARS)))
  OBJS    += GFM/agn_radiation.o
  SUBDIRS += GFM
endif

ifeq (GFM_STELLAR_PHOTOMETRICS, $(findstring GFM_STELLAR_PHOTOMETRICS, $(CONFIGVARS)))
  OBJS    += GFM/stellar_photometrics.o GFM/stellar_photometrics_vars.o
  INCL    += GFM/stellar_photometrics_vars.h GFM/stellar_photometrics_proto.h
  SUBDIRS += GFM
endif

ifeq (GFM_STELLAR_EVOLUTION, $(findstring GFM_STELLAR_EVOLUTION, $(CONFIGVARS)))
  OBJS    += GFM/stellar_evolution_init.o GFM/stellar_evolution_main.o GFM/stellar_evolution_evolve.o \
             GFM/stellar_evolution_util.o GFM/stellar_evolution_vars.o GFM/stellar_density.o          \
             GFM/stellar_evolution_logs.o GFM/stellar_evolution_dust.o GFM/stellar_evolution_enrich.o
  INCL    += GFM/stellar_evolution_proto.h GFM/stellar_evolution_vars.h
  SUBDIRS += GFM
endif

ifeq (GFM_COOLING_METAL, $(findstring GFM_COOLING_METAL, $(CONFIGVARS)))
  OBJS    += GFM/cooling_metal.o GFM/cooling_metal_vars.o
  INCL    += GFM/cooling_metal_proto.h GFM/cooling_metal_vars.h
  SUBDIRS += GFM
endif

ifeq (GFM_DUST_COOLING, $(findstring GFM_DUST_COOLING, $(CONFIGVARS)))
  OBJS    += GFM/cooling_dust.o
  INCL    += GFM/cooling_dust_proto.h
  SUBDIRS += GFM
endif

ifeq (GFM_WINDS, $(findstring GFM_WINDS, $(CONFIGVARS)))
  OBJS    += GFM/winds.o GFM/winds_variable.o GFM/winds_vars.o  GFM/winds_local.o GFM/winds_recouple.o GFM/winds_findcells.o
  INCL    += GFM/winds_proto.h
  SUBDIRS += GFM
endif

ifeq (GFM, $(findstring GFM, $(CONFIGVARS)))
  OBJS    += GFM/helper.o
  INCL    += GFM/helper_proto.h
  SUBDIRS += GFM
endif

inc_stellar_feedback =
ifeq (GFM_STELLAR_FEEDBACK, $(findstring GFM_STELLAR_FEEDBACK, $(CONFIGVARS)))
  inc_stellar_feedback = yes
endif
ifeq (GFM_WINDS_LOCAL, $(findstring GFM_WINDS_LOCAL, $(CONFIGVARS)))
  inc_stellar_feedback = yes
endif

ifdef inc_stellar_feedback
  OBJS    += GFM/stellar_feedback.o GFM/stellar_feedback_kernels.o
  INCL    += GFM/stellar_feedback_kernels.h
  SUBDIRS += GFM
endif

ifeq (SUBFIND, $(findstring SUBFIND, $(CONFIGVARS)))
  OBJS    += subfind/subfind.o subfind/subfind_vars.o  subfind/subfind_serial.o  subfind/subfind_coll_tree.o \
             subfind/subfind_properties.o subfind/subfind_so.o  subfind/subfind_distribute.o \
             subfind/subfind_collective.o subfind/subfind_findlinkngb.o subfind/subfind_nearesttwo.o \
             subfind/subfind_loctree.o subfind/subfind_coll_domain.o  subfind/subfind_coll_treewalk.o \
             subfind/subfind_io.o subfind/subfind_sort_kernels.o subfind/subfind_mark_cgm.o \
             subfind/subfind_reprocess.o subfind/subfind_so_potegy.o
## subdir subfind and header subfind.h and file subfind/subfind_density.o already default
endif

ifeq (AMR, $(findstring AMR, $(CONFIGVARS)))
  OBJS    += amr/amr.o amr/amr_mesh.o amr/amr_refinement.o amr/amr_exchange.o amr/amr_ngb.o amr/amr_refinement_criterion.o amr/amr_validate.o amr/amr_remap.o amr/amr_gradients.o amr/amr_generate_gas.o amr/amr_update_nodes.o amr/amr_walk.o
  ifeq (TWODIMS, $(findstring TWODIMS, $(CONFIGVARS)))
    OBJS  += amr/amr_2d.o
  endif
  INCL    += amr/amr.h amr/amr_proto.h
  SUBDIRS += amr
endif

ifeq (DG, $(filter DG, $(CONFIGVARS)))
  OBJS    += dg/dg_vars.o dg/dg_legendre.o dg/dg_core.o dg/dg_debug.o dg/dg_set_get.o dg/dg_limiter.o dg/dg_fluxes.o dg/dg_recompute.o dg/dg_projection.o dg/dg_time_integration.o dg/dg_limiter_special.o
  OBJS    += dg/dg_io.o dg/dg_refinement.o
  INCL    += dg/dg_vars.h dg/dg_proto.h dg/dg_defines.h dg/dg_core_inline.h
  SUBDIRS += dg
  EXEC = Tenet
endif

ifeq (DG_TEST_PROBLEM, $(filter DG_TEST_PROBLEM, $(CONFIGVARS)))
  OBJS    += dg/dg_test_problems.o dg/cell_projection.o
  INCL    += dg/dg_test_problems.h
  SUBDIRS += dg
endif

ifeq (VS_TURB, $(findstring VS_TURB, $(CONFIGVARS)))
  OBJS    += turb/turb_driving.o turb/turb_powerspectra.o
  SUBDIRS += turb
endif

ifeq (AB_TURB, $(findstring AB_TURB, $(CONFIGVARS)))
  OBJS    += turb/ab_turb.o  turb/turb_driving.o turb/turb_powerspectra.o
  SUBDIRS += turb
endif

ifeq (MHD_SEEDPSPEC, $(findstring MHD_SEEDPSPEC, $(CONFIGVARS)))
  OBJS    += constrained_transport/mhd_seedpspec.o
  SUBDIRS += constrained_transport
endif

ifeq (REGULARIZE_MESH_LLOYD, $(findstring REGULARIZE_MESH_LLOYD, $(CONFIGVARS)))
  OBJS    += constrained_transport/regularize_mesh_more_options.o
  SUBDIRS += constrained_transport
endif

ifeq (REGULARIZE_MESH_SMOOTH, $(findstring REGULARIZE_MESH_SMOOTH, $(CONFIGVARS)))
  OBJS    += constrained_transport/regularize_mesh_more_options.o
  SUBDIRS += constrained_transport
endif

ifeq (TGSET, $(findstring TGSET, $(CONFIGVARS)))
  OBJS    += tgset/tgset.o
  INCL    += tgset/tgset.h tgset/tgset_proto.h
  SUBDIRS += tgset
endif

ifeq (TGCHEM, $(findstring TGCHEM, $(CONFIGVARS)))
  OBJS    += tgchem/tgchem.o tgchem/tgchem_init.o tgchem/tgchem_rates.o tgchem/tgchem_step.o tgchem/tgchem_utils.o
  INCL    += tgchem/tgchem.h tgchem/tgchem_proto.h
  SUBDIRS += tgchem
endif

ifeq (HEALRAY, $(findstring HEALRAY, $(CONFIGVARS)))
  OBJS    += healray/healray.o healray/healray_comm.o healray/healray_experimental.o healray/healray_finish.o healray/healray_init.o healray/healray_rayout.o healray/healray_sources.o
  INCL    += healray/healray.h healray/healray_proto.h
  SUBDIRS += healray
endif

ifeq (SGCHEM, $(findstring SGCHEM, $(CONFIGVARS)))
  OBJS    += SGChem/calc_shield.o SGChem/cheminmo.o SGChem/cma-setup.o SGChem/cma-util.o \
            SGChem/const_rates.o SGChem/cool_func.o SGChem/cool_util.o \
            SGChem/coolinmo.o SGChem/dvode.o SGChem/evolve_abundances.o SGChem/init_chemistry_parameters.o \
            SGChem/jac.o SGChem/photoinit_ism.o SGChem/photoinit_lowZ.o SGChem/rate_eq_simple.o SGChem/rate_eq_primordial.o \
            SGChem/rate_eq_nl99.o SGChem/rate_eq_mrt.o SGChem/sgchem.o SGChem/spline.o SGChem/validate_rates.o SGChem/set_local_abundances.o SGChem/lwbg.o \
            SGChem/calc_temp.o SGChem/compute_heating.o SGChem/compute_md_kappa.o SGChem/rate_eq_gong17.o SGChem/rate_eq_primordial_noH2.o
  INCL    += SGChem/cma.h SGChem/cool.h SGChem/f2c.h SGChem/fs_data.h SGChem/h2heat.h  SGChem/isrf.h \
            SGChem/mol_data.h SGChem/non_eq.h SGChem/sgchem_def.h SGChem/sgchem_proto.h SGChem/shield_data.h \
            SGChem/gamma.h SGChem/kappa-planck.h
  SUBDIRS += SGChem
endif

ifeq (TREECOLV2, $(findstring TREECOLV2, $(CONFIGVARS)))
  OBJS    += TreeColV2/treecolv2_utils.o
  INCL    += TreeColV2/treecolv2_proto.h TreeColV2/treecolv2.h
  SUBDIRS += TreeColV2
endif

ifeq (SNE_FEEDBACK, $(findstring SNE_FEEDBACK, $(CONFIGVARS)))
  OBJS    += sne/sne.o sne/sne_injection_criteria.o sne/sne_utility.o sne/sne_time_stepping.o
  INCL    += sne/sne.h sne/sne_proto.h
  SUBDIRS += sne
endif

ifeq (SINK_PARTICLES, $(findstring SINK_PARTICLES, $(CONFIGVARS)))
  OBJS    += sink_particles/accrete_onto_sink_particles.o sink_particles/dump_sink_particle_info.o \
             sink_particles/init_sink_particles.o  sink_particles/sink_particles.o \
             sink_particles/create_sink_particles.o sink_particles/get_all_sink_particle_info.o
  INCL    += sink_particles/proto_sink_particles.h sink_particles/sink_particles.h
  SUBDIRS += sink_particles
endif

ifeq (SINK_PARTICLES_FEEDBACK, $(findstring SINK_PARTICLES_FEEDBACK, $(CONFIGVARS)))
  OBJS    += sink_particles/sink_feedback.o
endif

ifeq (SINK_MERGERS, $(findstring SINK_MERGERS, $(CONFIGVARS)))
  OBJS    += sink_particles/perform_mergers.o
endif

ifeq (SINK_PHOTOION_FEEDBACK, $(findstring SINK_PHOTOION_FEEDBACK, $(CONFIGVARS)))
  OBJS    += sink_particles/sink_photoion.o
endif

ifeq (SINKS, $(findstring SINKS, $(CONFIGVARS)))
  OBJS    += sinks/sinks.o sinks/sinks_accrete.o sinks/sinks_aux.o sinks/sinks_create.o sinks/sinks_init.o \
             sinks/sinks_dmass.o sinks/sinks_neighbors.o sinks/sinks_mergers.o
  INCL    += sinks/sinks.h sinks/sinks_proto.h
  SUBDIRS += sinks
endif

ifeq (SIDM, $(findstring SIDM, $(CONFIGVARS)))
  OBJS    += sidm/sidm_hsml.o sidm/sidm_ngb.o sidm/sidm_scatter.o sidm/sidm_cross.o sidm/sidm.o sidm/sidm_vars.o sidm/sidm_scatter_process.o
  INCL    += sidm/sidm_vars.h sidm/sidm_proto.h
  SUBDIRS += sidm
endif

ifeq (ADJ_BOX_POWERSPEC, $(findstring ADJ_BOX_POWERSPEC, $(CONFIGVARS)))
  OBJS    += power_spec/adj_box_powerspec.o
  INCL    += power_spec/adj_box_powerspec_proto.h
  SUBDIRS += power_spec
endif

ifeq (SMUGGLE_STAR_FEEDBACK,$(findstring SMUGGLE_STAR_FEEDBACK,$(CONFIGVARS)))
OBJS    += GFM/stellar_feedback_kernels.o GFM/stellar_feedback.o SMUGGLE/stellar_feedback_util.o SMUGGLE/smuggle_feedback_kernels.o
INCL    += GFM/stellar_feedback_kernels.h SMUGGLE/stellar_feedback_proto.h SMUGGLE/smuggle_feedback_kernels.h

ifeq (SMUGGLE_OUTPUT_STELLAR_FEEDBACK,$(findstring SMUGGLE_OUTPUT_STELLAR_FEEDBACK,$(CONFIGVARS)))
OBJS    += SMUGGLE/stellar_feedback_logs.o
endif

SUBDIRS += SMUGGLE GFM
endif

ifeq (SMUGGLE_MOLEC_COOLING,$(findstring SMUGGLE_MOLEC_COOLING,$(CONFIGVARS)))
OBJS    += SMUGGLE/cooling_molecules.o
INCL    += SMUGGLE/cooling_molecules_proto.h
SUBDIRS += SMUGGLE
endif

ifeq (SMUGGLE_DUST_HEATING_COOLING,$(findstring SMUGGLE_DUST_HEATING_COOLING,$(CONFIGVARS)))
OBJS    += SMUGGLE/heating_cooling_dust.o
INCL    += SMUGGLE/heating_cooling_dust_proto.h
SUBDIRS += SMUGGLE
endif

ifeq (SMUGGLE_SFR, $(findstring SMUGGLE_SFR, $(CONFIGVARS)))
  OBJS    += SMUGGLE/sfr.o
  SUBDIRS += SMUGGLE
endif

ifeq (LOCAL_FEEDBACK, $(findstring LOCAL_FEEDBACK, $(CONFIGVARS)))
  OBJS    += local_feedback/compute_sfr.o local_feedback/inject_feedback.o local_feedback/create_star_particles.o local_feedback/RT_inject_photons.o
  SUBDIRS += local_feedback
endif

ifeq (SMUGGLE_COMPUTE_SFR_FROM_H2, $(findstring SMUGGLE_COMPUTE_SFR_FROM_H2, $(CONFIGVARS)))
  OBJS    += SMUGGLE/H2_frac.o
  INCL    += SMUGGLE/H2_frac_proto.h
  SUBDIRS += SMUGGLE
endif

ifeq (SMUGGLE_TEST_SFR, $(findstring SMUGGLE_TEST_SFR, $(CONFIGVARS)))
  OBJS    += SMUGGLE/test_sfr.o
  SUBDIRS += SMUGGLE
endif

ifeq (OTVET, $(findstring OTVET, $(CONFIGVARS)))
  OBJS    += OTVET/do_otvet.o OTVET/otvet_eddington.o OTVET/otvet_star_lum.o OTVET/otvet_CGmethod.o OTVET/otvet_chem_ps2009.o OTVET/otvet_chem_ps2011.o OTVET/otvet_cooling.o
  INCL    += OTVET/otvet_proto.h
  SUBDIRS += OTVET
endif

ifeq (SMUGGLE_RADPRESS_OPT_THIN, $(findstring SMUGGLE_RADPRESS_OPT_THIN, $(CONFIGVARS)))
  OBJS    += SMUGGLE/radpressthin.o
  INCL    += SMUGGLE/radpressthin_proto.h
  SUBDIRS += SMUGGLE
endif

ifeq (SMUGGLE_RADPRESS_OPT_THICK, $(findstring SMUGGLE_RADPRESS_OPT_THICK, $(CONFIGVARS)))
  OBJS    += SMUGGLE/radpressthick.o
  INCL    += SMUGGLE/radpressthick_proto.h
  SUBDIRS += SMUGGLE
endif

ifeq (SMUGGLE_RADIATION_FEEDBACK, $(findstring SMUGGLE_RADIATION_FEEDBACK, $(CONFIGVARS)))
  OBJS    += SMUGGLE/radiation_stellar_feedback_find_cells.o SMUGGLE/radiation_stellar_feedback.o SMUGGLE/radiation_stellar_feedback_util.o
  INCL    += SMUGGLE/radiation_stellar_feedback_proto.h
  SUBDIRS += SMUGGLE GFM
endif

ifeq (TEST_COOLING_METAL, $(findstring TEST_COOLING_METAL, $(CONFIGVARS)))
  OBJS    += Test_Cooling_Metal/test_cooling_metal.o
  INCL    += Test_Cooling_Metal/test_cooling_metal_proto.h
  SUBDIRS += Test_Cooling_Metal
endif

ifeq (NUCLEAR_NETWORK, $(findstring NUCLEAR_NETWORK, $(CONFIGVARS)))
  OBJS    += network/network.o network/network_solver.o network/utilities.o network/integrate.o network/network_nse.o network/integrate_mpi.o
  INCL    += network/network.h network/network_solver.h network/utilities.h network/integrate.h network/network_nse.h
  SUBDIRS += network
endif

ifeq (HCOUTPUT, $(findstring HCOUTPUT, $(CONFIGVARS)))
  OBJS    += highcadoutput/hcoutput.o
  INCL    += highcadoutput/hcoutput.h
  SUBDIRS += highcadoutput
endif

ifeq (MRT, $(findstring MRT, $(CONFIGVARS)))
  OBJS    += MRT/RT.o MRT/RT_set_VET.o MRT/RT_cooling.o MRT/RT_chem_ps2011.o MRT/RT_chem_ps2009.o MRT/RT_riemann_rosunov.o \
             MRT/RT_IR.o MRT/RT_setup.o MRT/RT_riemann_HLLE.o MRT/RT_comoving.o MRT/RT_add_sources_stellar.o \
             MRT/RT_add_sources_blackholes.o MRT/RT_add_radiation_stellar.o MRT/RT_add_radiation_blackholes.o \
             MRT/RT_source_utils.o MRT/RT_finite_volume_solver.o MRT/RT_gradients.o MRT/RT_gradients_lsf.o MRT/RT_exchange.o \
             MRT/RT_run.o MRT/RT_update_primitive_variables.o MRT/RT_spectra.o MRT/RT_riemann_rosunov_new.o MRT/RT_riemann_HLLE_new.o \
             MRT/RT_thermochem.o MRT/RT_init.o MRT/RT_sgchem.o MRT/RT_set_kappa.o MRT/RT_radiation_pressure.o
  INCL    += MRT/RT.h MRT/RT_proto.h
  SUBDIRS += MRT
endif

ifeq (CONDUCTION, $(findstring CONDUCTION, $(CONFIGVARS)))
  OBJS    += conduction.o
endif

ifeq (MONOTONE_CONDUCTION, $(findstring MONOTONE_CONDUCTION, $(CONFIGVARS)))  # Conduction module
  ifeq (IMPLICIT_TI, $(findstring IMPLICIT_TI, $(CONFIGVARS)))
    ifeq (SEMI_IMPLICIT_TI, $(findstring SEMI_IMPLICIT_TI, $(CONFIGVARS)))
      OBJS += conduction/conduction_semi_HYPRE.o
    else
      OBJS += conduction/conduction_HYPRE.o
    endif
  else
    OBJS   += conduction/conduction.o
  endif
  ifeq (MULTIPLE_TIME_STEPPING, $(findstring MULTIPLE_TIME_STEPPING, $(CONFIGVARS)))
    OBJS   += conduction/conduction_multiple_timestep.o
  endif
  INCL     += conduction/conduction.h
  SUBDIRS  += conduction
endif

ifeq (CALCULATE_QUANTITIES_IN_POSTPROCESS, $(findstring CALCULATE_QUANTITIES_IN_POSTPROCESS, $(CONFIGVARS)))
  OBJS    += calculate_quantities_in_postprocess/calculate_quantities.o
  INCL    += calculate_quantities_in_postprocess/calculate_quantities.h
  SUBDIRS += calculate_quantities_in_postprocess
endif

ifeq (CIRCUMSTELLAR, $(findstring CIRCUMSTELLAR, $(CONFIGVARS)))
  OBJS    += circumstellar/circumstellar.o circumstellar/circumstellar_sinks.o circumstellar/circumstellar_gravity.o circumstellar/circumstellar_kepler.o
  INCL    += circumstellar/circumstellar_proto.h
  SUBDIRS += circumstellar
endif

ifeq (ATOMIC_DM, $(findstring ATOMIC_DM, $(CONFIGVARS)))
  OBJS    += atomic_dm/cooling_atomic_dm.o
  INCL    += atomic_dm/cooling_atomic_dm_proto.h atomic_dm/cooling_atomic_dm_vars.h
  SUBDIRS += atomic_dm
endif

ifeq (FLD, $(filter FLD, $(CONFIGVARS)))
  OBJS    += fld/fld.o fld/fld_MGmethod.o fld/fld_HYPRE.o fld/fld_HYPRE_IJ.o
  INCL    += fld/fld_proto.h fld/fld.h
  SUBDIRS += fld
endif

ifeq (CUDA, $(findstring CUDA, $(CONFIGVARS)))
  OBJS    += cuda_util.o
  INCL    += cuda_util.h
endif

ifeq (DM_WINDTUNNEL, $(findstring DM_WINDTUNNEL, $(CONFIGVARS)))
  OBJS    += dmwindtunnel/dmwindtunnel.o
  INCL    += dmwindtunnel/dmwindtunnel.h
  SUBDIRS += dmwindtunnel
endif

ifeq (SHOCK_FINDER, $(findstring SHOCK_FINDER, $(CONFIGVARS)))
  OBJS    += shock_finder/shock_finder.o shock_finder/shock_finder_ryu.o shock_finder/shock_finder_skillman.o shock_finder/shock_finder_arepo.o
  INCL    += shock_finder/shock_finder.h shock_finder/shock_finder_fields.h shock_finder/shock_finder_rays.h
  SUBDIRS += shock_finder
endif

ifeq (MHD_CT, $(findstring MHD_CT, $(CONFIGVARS)))
  OBJS    += constrained_transport/constrained_transport.o
  INCL    += constrained_transport/constrained_transport.h
  SUBDIRS += constrained_transport
endif

ifeq (BECDM, $(findstring BECDM, $(CONFIGVARS)))
  OBJS    += becdm/becdm.o
  INCL    += becdm/becdm.h
  SUBDIRS += becdm
endif

ifeq (COSMIC_RAYS, $(findstring COSMIC_RAYS, $(CONFIGVARS)))
  OBJS    += cosmic_rays/cosmic_rays.o
  INCL    += cosmic_rays/cosmic_rays.h
  SUBDIRS += cosmic_rays
endif

ifeq (COSMIC_RAYS_SN_INJECTION, $(findstring COSMIC_RAYS_SN_INJECTION, $(CONFIGVARS)))
  OBJS    += cosmic_rays/cosmic_rays_injection.o cosmic_rays/cosmic_rays_find_ngbs.o
endif

ifeq (COSMIC_RAYS_DIFFUSION, $(findstring COSMIC_RAYS_DIFFUSION, $(CONFIGVARS)))
  OBJS    += cosmic_rays/cosmic_rays_diffusion.o cosmic_rays/cosmic_rays_set_diffusion_coefficients.o
endif

ifeq (COSMIC_RAYS_STREAMING, $(findstring COSMIC_RAYS_STREAMING, $(CONFIGVARS)))
  OBJS    += cosmic_rays/cosmic_rays_streaming.o
endif

ifeq (COSMIC_RAYS_SHOCK_ACCELERATION, $(findstring COSMIC_RAYS_SHOCK_ACCELERATION, $(CONFIGVARS)))
  OBJS    += cosmic_rays/cosmic_rays_shock_acceleration.o
endif

ifeq (BRAGINSKII_VISCOSITY,$(findstring BRAGINSKII_VISCOSITY,$(CONFIGVARS)))
  OBJS    += braginskii_viscosity/braginskii_viscosity.o
  INCL    += braginskii_viscosity/braginskii_viscosity.h
  SUBDIRS += braginskii_viscosity
endif

ifeq (NON_IDEAL_MHD, $(findstring NON_IDEAL_MHD, $(CONFIGVARS)))
  OBJS    += mhd_nonideal/ohmic_diffusion.o
  OBJS    += mhd_nonideal/ambipolar_diffusion.o
  INCL    += mhd_nonideal/proto_non_ideal_mhd.h
  SUBDIRS += mhd_nonideal
endif

ifeq (IMPLICIT_OHMIC_DIFFUSION, $(findstring IMPLICIT_OHMIC_DIFFUSION, $(CONFIGVARS)))
  ifeq (MHD_CT, $(findstring MHD_CT, $(CONFIGVARS)))
    ifeq (OHM_CRANK_NICHOLSON, $(findstring OHM_CRANK_NICHOLSON, $(CONFIGVARS)))
      OBJS += mhd_nonideal/implicit_ohmic_diffusion_ct_cn.o
    else
      OBJS += mhd_nonideal/implicit_ohmic_diffusion_ct.o
    endif
  else
    ifeq (OHM_CRANK_NICHOLSON, $(findstring OHM_CRANK_NICHOLSON, $(CONFIGVARS)))
      OBJS += mhd_nonideal/implicit_ohmic_diffusion_cn.o
      OBJS += mhd_nonideal/implicit_ohmic_diffusion_cn_nonperiodic.o
    else
      OBJS += mhd_nonideal/implicit_ohmic_diffusion.o
      OBJS += mhd_nonideal/implicit_ohmic_diffusion_nonperiodic.o
    endif
  endif
  INCL     += mhd_nonideal/proto_non_ideal_mhd.h
  SUBDIRS  += mhd_nonideal
endif

ifeq (EOS_OPAL, $(findstring EOS_OPAL, $(CONFIGVARS)))
  OBJS    += opal_eos.o
  INCL    += opal_eos.h
endif

ifeq (AURIGA_MOVIE, $(findstring AURIGA_MOVIE, $(CONFIGVARS)))
  OBJS    += auriga_movie/movie.o auriga_movie/movie_util.o auriga_movie/movie_density.o auriga_movie/movie_projection.o
  INCL    += auriga_movie/movie.h
  SUBDIRS += auriga_movie
endif

ifeq (SPECIAL_RELATIVITY, $(findstring SPECIAL_RELATIVITY, $(CONFIGVARS)))
  OBJS    += special_relativity.o
  ifeq (SPECIAL_RELATIVITY_HLLC, $(findstring SPECIAL_RELATIVITY_HLLC, $(CONFIGVARS)))
    OBJS  += riemann_hllc_special_relativity.o
  else
    OBJS  += riemann_hlle_special_relativity.o
  endif
endif

ifeq (GENERAL_RELATIVITY, $(findstring GENERAL_RELATIVITY, $(CONFIGVARS)))
  OBJS    += general_relativity.o
  ifeq (SPECIAL_RELATIVITY_HLLC, $(findstring SPECIAL_RELATIVITY_HLLC, $(CONFIGVARS)))
    OBJS  += riemann_hllc_general_relativity.o
  else
    OBJS  += riemann_hlle_general_relativity.o
  endif
endif

ifeq (GRACKLE, $(findstring GRACKLE, $(CONFIGVARS)))
  OBJS    += grackle/grackle.o
  INCL    += grackle/grackle_def.h
  SUBDIRS += grackle
endif

ifeq (GRAVITY_TABLE, $(findstring GRAVITY_TABLE, $(CONFIGVARS)))
  OBJS    += gravity_table.o
endif

ifeq (MODGRAV, $(findstring MODGRAV, $(CONFIGVARS)))
  OBJS    += modgrav/modgrav_forcetree.o modgrav/modgrav_pm_periodic.o modgrav/modgrav_pm_nonperiodic.o
  INCL    += modgrav/modgrav_forcetree.h modgrav/modgrav_pm.h
  SUBDIRS += modgrav
endif

ifeq (SIMPLEX, $(findstring SIMPLEX, $(CONFIGVARS)))
  OBJS    += simplex/sx_vars.o simplex/sx_evolve.o simplex/sx_math.o simplex/sx_debug.o # general files
  OBJS    += simplex/sx_photon.o simplex/sx_photon_transport.o simplex/sx_photon_exchange.o simplex/sx_photon_sources.o
  OBJS    += simplex/sx_chem_sgchem.o simplex/sx_chem_sgchem_int.o       # SGChem chemistry
  OBJS    += simplex/sx_chem_fiby.o                                      # FiBY chemistry
  OBJS    += simplex/lookup_table_module.o simplex/popiiistar_module.o   # PopIII star module
  INCL    += simplex/sx_proto.h simplex/sx_direction_bins.h simplex/sx_test_sources.h simplex/sx_api.h
  SUBDIRS += simplex
endif

ifeq (DUST_LIVE, $(findstring DUST_LIVE, $(CONFIGVARS)))
  OBJS    += dust_live/drag_backreaction.o dust_live/drag_kernels.o dust_live/drag_kicks.o dust_live/dust_density.o dust_live/dust_derefinement.o dust_live/dust_init.o dust_live/dust_neighbors.o dust_live/dust_production.o dust_live/dust_timestep.o dust_live/dust_transfer.o dust_live/dust_transfer_kernels.o dust_live/dust_util.o dust_live/dust_vars.o dust_live/dust_winds.o dust_live/grain_sizes.o dust_live/radiation_pressure.o dust_live/radiation_pressure_kernels.o dust_live/radiation_thermal_kernels.o
  INCL    += dust_live/dust_proto.h dust_live/dust_vars.h
  SUBDIRS += dust_live
endif

ifeq (BAROTROPIC, $(findstring BAROTROPIC, $(CONFIGVARS)))
  OBJS    += barotropic/barotropic.o
  SUBDIRS += barotropic
endif

ifeq (PERTURB_VELOCITIES, $(findstring PERTURB_VELOCITIES, $(CONFIGVARS)))
  OBJS    += perturb_velocities.o
endif

ifeq (OPACITIES, $(findstring OPACITIES, $(CONFIGVARS)))
  OBJS    += opacities/opacities_combined.o opacities/xztrin21.o
  INCL    += opacities/opacities_combined.h
  SUBDIRS += opacities
endif

ifeq (SGS_TURBULENCE, $(findstring SGS_TURBULENCE, $(CONFIGVARS)))
  OBJS    += sgs_turbulence/sgs_turbulence.o sgs_turbulence/sgs_turbulence_stress_tensor.o sgs_turbulence/sgs_turbulence_turbulent_production.o sgs_turbulence/sgs_turbulence_eddy_viscosity_closure.o sgs_turbulence/sgs_turbulence_viscous_dissipation.o
  INCL    += sgs_turbulence/sgs_turbulence.h
  SUBDIRS += sgs_turbulence
endif

ifeq (TURBULENT_METALDIFFUSION, $(findstring TURBULENT_METALDIFFUSION, $(CONFIGVARS)))
  OBJS    += mm.o
endif

ifeq (CHIMES, $(findstring CHIMES, $(CONFIGVARS)))
  OBJS    += chimes/chimes.o chimes/cooling.o chimes/init_chimes.o chimes/init_chimes_parallel.o chimes/interpol.o chimes/optimise.o chimes/rate_coefficients.o chimes/rate_equations.o chimes/set_rates.o cooling/chimes_cooling.o
  INCL    += chimes/allvars.h chimes/proto.h
  SUBDIRS += cooling chimes
endif

ifeq (SFR_MCS, $(findstring SFR_MCS, $(CONFIGVARS)))
  OBJS    += sfr_mcs/sfr_mcs.o sfr_mcs/sfr_mcs_cooling.o sfr_mcs/feedback_mcs.o sfr_mcs/feedback_mcs_utils.o sfr_mcs/imf_sampling_mcs.o sfr_mcs/turb_approx_mcs.o
  ifeq (SN_MCS, $(findstring SN_MCS, $(CONFIGVARS)))
    OBJS    += sfr_mcs/sn_mcs.o sfr_mcs/sn_mcs_utils.o sfr_mcs/sn_mcs_inject.o
  endif
  ifeq (HII_MCS,$(findstring HII_MCS,$(CONFIGVARS)))
    OBJS    += sfr_mcs/hii_mcs.o sfr_mcs/hii_mcs_anisotropic.o sfr_mcs/hii_mcs_healpix_utils.o
  endif
  ifeq (PE_MCS,$(findstring PE_MCS,$(CONFIGVARS)))
    OBJS    += sfr_mcs/pe_mcs.o sfr_mcs/estimate_local_dust_column_mcs.o
  endif
  INCL    += sfr_mcs/sfr_mcs_vars.h sfr_mcs/sfr_mcs_proto.h
  SUBDIRS += sfr_mcs
endif

ifeq (GALPOT, $(findstring GALPOT, $(CONFIGVARS)))
  OBJS    += galpot/galpot.o galpot/potential.o
  INCL    += galpot/galpot.h galpot/common.h galpot/potential.h
  SUBDIRS += galpot
endif


ifeq (BIERMANN_BATTERY,$(findstring BIERMANN_BATTERY,$(CONFIGVARS)))
  OBJS    += magnetic_batteries.o
else
  ifeq (DURRIVE_BATTERY,$(findstring DURRIVE_BATTERY,$(CONFIGVARS)))
    OBJS    += magnetic_batteries.o
  endif
endif

ifeq (SOLAR, $(findstring SOLAR, $(CONFIGVARS)))
  OBJS    += solar/solar.o
  INCL    += solar/solar.h
  SUBDIRS += solar
  ifeq (SOLAR_RADIATIVE_TRANSFER_DIFF, $(findstring SOLAR_RADIATIVE_TRANSFER_DIFF, $(CONFIGVARS)))
    OBJS  += solar/radiative_transfer.o
  endif
  ifeq (SOLAR_RADIATIVE_TRANSFER_EDD, $(findstring SOLAR_RADIATIVE_TRANSFER_EDD, $(CONFIGVARS)))
    OBJS  += solar/radiative_transfer.o solar/radiative_transfer_edd.o
  endif
endif
