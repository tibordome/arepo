/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/main.c
 * \date        MM/YYYY
 * \author
 * \brief       Start of the program
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

static void print_command_line_help(void);
static void print_command_line_arguments(int argc, char **argv, const char *parameter_file, int restart_flag, int restart_snap_num);
static void endrun_with_status(int status);

/*! \brief The entry point of the program.
 *
 *  This function initializes the MPI communication packages, and sets
 *  CPU time counters to 0. Then begrun1() is called, which sets up
 *  the simulation. Then either IC's or restart files are loaded. In
 *  case of IC's init() is called which prepares the IC's for the run.
 *  A call to begrun2() finishes the initialization. Finally, run() is
 *  started, the main simulation loop, which iterates over the timesteps.
 *
 *  \param[in] argc Argument count from command line.
 *  \param[in] argv Argument vector from command line.
 *
 *  \return status of exit; 0 for normal exit.
 */
int main(int argc, char **argv)
{
#ifdef IMPOSE_PINNING
  detect_topology();
  get_core_set();
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

#ifdef CHIMES_PTHREADS
  /* Create a communicator containing MPI tasks on the same node. */
  /* To determine which node a given task is on, it would be
   * simplest to use MPI_COMM_TYPE_SHARED. However, this was only
   * introduced in MPI-3, and it looks like it is incompatible
   * with older MPI implementations. Instead, convert the host name
   * to an integer, which can be passed to MPI_Comm_split() to
   * identify the node. */
  char node_name[MPI_MAX_PROCESSOR_NAME];
  int node_name_length, i;
  int node_colour = 0;
  MPI_Get_processor_name(node_name, &node_name_length);
  for(i = 0; i < node_name_length; i++)
    node_colour += (i + 1) * node_name[i];

  MPI_Comm_split(MPI_COMM_WORLD, node_colour, 0, &node_comm);
  MPI_Comm_rank(node_comm, &ThisTask_node);
  MPI_Comm_size(node_comm, &NTask_node);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* output a welcome message */
  hello();

  /* initialize CPU time/wallclock time measurement */
  init_cpu_log();

  determine_compute_nodes();

#ifdef IMPOSE_PINNING
  /* pin the MPI ranks to the available core set */
  pin_to_core_set();
#endif

#if NUM_THREADS > 1
  char *p = getenv("OMP_NUM_THREADS");
  int num_threads;
  if(p)
    num_threads = atoi(p);
  else
    num_threads = 1;

  if(num_threads != NUM_THREADS)
    terminate(
        "You have activated NUM_THREADS, but the value of the environment\n variable OMP_NUM_THREADS=%d is not equal to "
        "NUM_THREADS=%d!",
        num_threads, NUM_THREADS);

  omp_set_num_threads(NUM_THREADS);
  omp_set_dynamic(0); /* explicitly turn off dynamic threads */
  mpi_printf("\nUsing %d OpenMP threads\n", omp_get_max_threads());
#ifdef IMPOSE_PINNING
  set_pinning_openmp_threads();
  report_pinning_openmp_threads();
#endif
#else
  mpi_printf("\nPINNING: We are not using OpenMP.\n\n");
#ifdef IMPOSE_PINNING
  report_pinning();
#endif
#endif

#ifdef HOST_MEMORY_REPORTING
  mpi_report_committable_memory();
#endif

  Argc = argc;
  Argv = argv;

  for(PTask = 0; NTask > (1 << PTask); PTask++)
    ;

  begrun0();

  if(Argc < 2)
    {
      mpi_printf("\nParameters are missing.\n");
      print_command_line_help();
      endrun_with_status(EXIT_FAILURE);
    }

  strcpy(ParameterFile, Argv[1]);

  if(Argc >= 3)
    RestartFlag = atoi(Argv[2]);
  else
    RestartFlag = RESTART_IC;

  if(Argc >= 4)
    RestartSnapNum = atoi(Argv[3]);
  else
    RestartSnapNum = -1;

  print_command_line_arguments(Argc, Argv, ParameterFile, RestartFlag, RestartSnapNum);

  /* Do minimal validation of arguments here rather than in random places in the code */
  if((RestartFlag == RESTART_FOF_SUBFIND || RestartFlag == RESTART_SLICE || RestartFlag == RESTART_PROJECTION ||
      RestartFlag == RESTART_SNAP_CONVERSION || RestartFlag == RESTART_GAS_VELOCITY_POWER_SPECTRUM ||
      RestartFlag == RESTART_PROJECTION_AXIS || RestartFlag == RESTART_TRACER_POWER_SPECTRA ||
      RestartFlag == RESTART_CORRELATION_FUNCTION || RestartFlag == RESTART_VORONOI_MESH || RestartFlag == RESTART_SHOCK_FINDER ||
      RestartFlag == RESTART_VORONOI_MESH_SLICE || RestartFlag == RESTART_GRADIENTS || RestartFlag == RESTART_RECALC_POTENTIAL ||
      RestartFlag == RESTART_CALC_ADDITIONAL || RestartFlag == RESTART_AURIGA_MOVIE || RestartFlag == RESTART_SIMPLEX ||
      RestartFlag == RESTART_CALC_VORONOI_DM_DENSITY) &&
     RestartSnapNum < 0)
    {
      mpi_terminate("Need to give the snapshot number");
    }
  if(RestartFlag == RESTART_POWER_SPECTRUM && Argc < 5)
    mpi_terminate("Need to give a particle type bitmask");

#ifndef RECOMPUTE_POTENTIAL_IN_SNAPSHOT
  if(RestartFlag == RESTART_RECALC_POTENTIAL)
    mpi_terminate("Need RECOMPUTE_POTENTIAL_IN_SNAPSHOT for this option");
#endif

#ifndef COMPUTE_VORONOI_DM_DENSITY_IN_POSTPROC
  if(RestartFlag == RESTART_CALC_VORONOI_DM_DENSITY)
    mpi_terminate("Need COMPUTE_VORONOI_DM_DENSITY_IN_POSTPROC for this option");
#endif

#if defined(GFM_STELLAR_EVOLUTION) && GFM_STELLAR_EVOLUTION == 2
  /* call only the stellar evolution test routine */
  test_stellar_evolution();
  /* clean up & finalize MPI */
  MPI_Finalize();
  return EXIT_SUCCESS;
#endif

#ifdef TEST_COOLING_METAL
  test_cooling_function(); /* call only the metal cooling test routine */
  test_cooling();
  MPI_Finalize(); /* clean up & finalize MPI */
  return EXIT_SUCCESS;
#endif

#ifdef SMUGGLE_TEST_SFR
  test_sfr();
  MPI_Finalize(); /* clean up & finalize MPI */
  return EXIT_SUCCESS;
#endif

#ifdef RUNNING_SAFETY_FILE
  /* do not run if 'running' safety file exists */
  int runningflag = 0;
  if(ThisTask == 0)
    {
      FILE *fd;
      char runningfname[MAXLEN_PATH];

      sprintf(runningfname, "./running");
      /* Is the running file present? If yes, interrupt the run. */
      if((fd = fopen(runningfname, "r")))
        {
          fclose(fd);
          printf("running-file detected. stopping.\n");
          runningflag = 1;
        }
    }
  MPI_Bcast(&runningflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(runningflag)
    {
      MPI_Finalize(); /* do not call endrun() */
      return EXIT_SUCCESS;
    }
  else
    {
      /* touch a running safety file */
      if(ThisTask == 0)
        {
          FILE *fd;
          char runningfname[MAXLEN_PATH];

          sprintf(runningfname, "./running");
          if((fd = fopen(runningfname, "w")))
            {
              fclose(fd);
              printf("touching a running-file: %s \n", runningfname);
            }
          else
            terminate("could not touch a running-file: %s\n", runningfname);
        }
    }
#endif

  /* set-up run  */
  begrun1();

  /* see if we are loading a restart file or an IC file */
  if(RestartFlag == RESTART_RESTART)
    loadrestart();
  else
    {
      /* We're reading an IC file. Is it a snapshot or really an IC? */
      char fname[MAXLEN_PATH];

      if(RestartFlag >= RESTART_SNAPSHOT && RestartSnapNum >= 0)
        {
          if(All.NumFilesPerSnapshot > 1)
            file_path_sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase, RestartSnapNum);
          else
            file_path_sprintf(fname, "%s/%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);

#ifdef SUBBOX_SNAPSHOTS
          if(RestartFlag == RESTART_FOF_SUBFIND || RestartFlag == RESTART_PROJECTION || RestartFlag == RESTART_PROJECTION_CAMERA ||
             RestartFlag == RESTART_SHOCK_FINDER)
            {
              if(RestartFlag == RESTART_FOF_SUBFIND && Argc < 5)
                mpi_terminate("subbox missing: %d", Argc);
              if(RestartFlag == RESTART_PROJECTION && Argc < 16)
                mpi_terminate("subbox missing: %d", Argc);
              if(RestartFlag == RESTART_PROJECTION_CAMERA && Argc != 8)
                {
                  mpi_terminate("subbox missing: %d", Argc);
                }
              else
                {
                  if(RestartFlag == RESTART_FOF_SUBFIND)
                    {
                      All.SubboxNumber              = atoi(Argv[4]);
                      All.NumFilesPerSnapshot       = All.SubboxNumFilesPerSnapshot;
                      All.NumFilesWrittenInParallel = All.SubboxNumFilesWrittenInParallel;
                    }
                  if(RestartFlag == RESTART_PROJECTION)
                    All.SubboxNumber = atoi(Argv[15]);
                  if(RestartFlag == RESTART_PROJECTION_CAMERA)
                    All.SubboxNumber = atoi(Argv[7]);
                  if(RestartFlag == RESTART_SHOCK_FINDER)
                    {
                      if(Argv[4] == 0)
                        mpi_terminate("<SubboxNum> is missing!");
                      All.SubboxNumber = atoi(Argv[4]);
                    }
                }
              if(All.SubboxNumFilesPerSnapshot > 1)
                file_path_sprintf(fname, "%s/snapdir_subbox%d_%03d/%s_subbox%d_%03d", All.OutputDir, All.SubboxNumber, RestartSnapNum,
                                  All.SnapshotFileBase, All.SubboxNumber, RestartSnapNum);
              else
                file_path_sprintf(fname, "%s/%s_subbox%d_%03d", All.OutputDir, All.SnapshotFileBase, All.SubboxNumber, RestartSnapNum);
            }
#endif
        }
      else
        strcpy(fname, All.InitCondFile);

#ifdef SHOCK_FINDER_POST_PROCESSING
      file_path_sprintf(All.OutputDir, "%s/", All.OutputDirShockFinder);
#endif

#ifdef POWERSPECTRUM_IN_POSTPROCESSING
      if(RestartFlag == RESTART_POWER_SPECTRUM)
        {
          int typemask = atoi(Argv[4]);

#ifdef POWERSPECTRUM_IN_POSTPROCESSING_ICS
          strcpy(All.InputFileName, All.InitCondFile);
#else
          strcpy(All.InputFileName, fname);
#endif

          int typeflag[NTYPES];
          for(int type = 0; type < NTYPES; type++)
            typeflag[type] = typemask & (0x1 << type);

          calculate_power_spectra_and_ntot(RestartSnapNum, typeflag);

          endrun();

          return EXIT_SUCCESS;
        }
      else
        {
          mpi_terminate("POWERSPECTRUM_IN_POSTPROCESSING only allows RestartFlag %d", RESTART_POWER_SPECTRUM);
        }
#endif

      /* now we can load the file */
      if(RestartFlag == RESTART_CORRELATION_FUNCTION || RestartFlag == RESTART_POWER_SPECTRUM)
        {
          int typemask = atoi(Argv[4]);
          read_ic(fname, typemask);
        }
      else /* readTypes=0x01: just load the gas particles. readTypes=LOAD_TYPES: load all particles */
        {
#ifdef READ_DM_AS_GAS
          read_ic(fname,
                  (RestartFlag == RESTART_SLICE || RestartFlag == RESTART_PROJECTION || RestartFlag == RESTART_PROJECTION_CAMERA ||
                   RestartFlag == RESTART_VORONOI_MESH || RestartFlag == RESTART_VORONOI_MESH_SLICE)
                      ? 0x02
                      : LOAD_TYPES);
#else
          read_ic(fname, (RestartFlag == RESTART_SLICE || RestartFlag == RESTART_PROJECTION ||
                          RestartFlag == RESTART_PROJECTION_CAMERA || RestartFlag == RESTART_VORONOI_MESH ||
                          RestartFlag == RESTART_SHOCK_FINDER || RestartFlag == RESTART_VORONOI_MESH_SLICE)
                             ? 0x01
                             : LOAD_TYPES);
#endif
        }

#ifdef CHIMES_INITIALISE_IN_EQM
      int partIndex;
      if(RestartFlag != RESTART_RESTART)
        {
          for(partIndex = 0; partIndex < NumGas; partIndex++)
            {
              SphP[partIndex].ChimesGasVars.index = partIndex;
              chimes_set_pointers(partIndex);
            }
        }
#endif /* CHIMES_INITIALISE_IN_EQM */

      /* If we are supposed to just convert the file, write and exit here. */
      if(RestartFlag == RESTART_SNAP_CONVERSION)
        {
#ifdef GFM_COOLING_METAL
          read_cooling_tables_current_time();
#endif
          /* important for proper functioning of FOF+SUBFIND */
          if(All.ComovingIntegrationOn) /* change to new velocity variable */
            {
              int i, j;
              for(i = 0; i < NumPart; i++)
                for(j = 0; j < 3; j++)
                  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
            }
          set_softenings();
          All.TopNodeAllocFactor = 0.08;
          All.TreeAllocFactor    = 0.7;
          All.NgbTreeAllocFactor = 0.7;

          strncat(All.SnapshotFileBase, "_converted", MAXLEN_PATH - strlen(All.SnapshotFileBase) - 1);
          mpi_printf("Start writing file %s\nRestartSnapNum %d\n", All.SnapshotFileBase, RestartSnapNum);
          savepositions(RestartSnapNum, 0);
          endrun();
        }

      /* init returns a status code, where a value of >= 0 means that endrun() should be called. */
      const int status = init();

      if(status >= 0)
        {
          if(status > 0)
            mpi_printf("init() returned with %d\n", status);
          endrun_with_status(status);
        }
    }

  begrun2();

  /* main simulation loop */
  run();

  /* clean up & finalize MPI */
  endrun();

  return EXIT_SUCCESS;
}

/*! \brief Print information on how to use the command line interface.
 *
 *  \return void
 */
static void print_command_line_help(void)
{
  if(ThisTask == 0)
    {
      printf("Call with <ParameterFile> [<RestartFlag>] [<RestartSnapNum>] [<SpecialOptions>]\n");
      printf("\n");
      printf("   RestartFlag    Action\n");
      printf("      %2d          Read initial conditions and start simulation\n", RESTART_IC);
      printf("      %2d          Read restart files and resume simulation\n", RESTART_RESTART);
      printf("      %2d          Restart from specified snapshot dump and resume simulation\n", RESTART_SNAPSHOT);
      printf(
          "      %2d          "
          "Run FOF and optionally SUBFIND: [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n",
          RESTART_FOF_SUBFIND);
      printf(
          "      %2d          Make an image slice:"
          "    <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> <xmin> <xmax> <ymin>"
          " <ymax> <zval>\n",
          RESTART_SLICE);
      printf(
          "      %2d          Make a projected image:"
          " <SnapNum> <pixelsX> <pixelsY> <axisX> <axisY> <axisZ> <xmin> <xmax>"
          " <ymin> <ymax> <zmin> <zmax> [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n",
          RESTART_PROJECTION);
      printf(
          "      %2d          "
          "Convert snapshot file to different format [input=ICFormat  output=SnapFormat]"
          "   NOTE: derived quantities have round-off errors!\n",
          RESTART_SNAP_CONVERSION);
      printf("      %2d          Calculate a velocity power spectrum for the gas cells\n", RESTART_GAS_VELOCITY_POWER_SPECTRUM);
      printf(
          "      %2d          Make a grid projection: <SnapNum>"
          " <pixelsX> <pixelsY> <pixelsZ> 0 0  <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>\n",
          RESTART_PROJECTION_GRID_RAYTRACING);
      printf(
          "      %2d          Make a projection along an arbitrary axis:"
          " <SnapNum> <pixelsX> <pixelsY> <centerX> <centerY> <centerZ>"
          " <dirX> <dirY> <dirZ> <boxX> <boxY> <boxZ>\n",
          RESTART_PROJECTION_AXIS);
      printf(
          "      %2d          Make a perspective camera projection:"
          " <SnapNum> <pixelsX> <pixelsY> <filename of camera file>"
          " [<SubboxSnapNum> for SUBBOX_SNAPSHOTS]\n",
          RESTART_PROJECTION_CAMERA);
      printf("      %2d          Calculate power spectra of various quantities for TRACER_PARTICLEs\n", RESTART_TRACER_POWER_SPECTRA);
      printf(
          "      %2d          Calculate two-point correlation function:"
          " <SnapNum> <parttype bitmask> [output path]\n",
          RESTART_CORRELATION_FUNCTION);
      printf(
          "      %2d          Calculate power spectrum:"
          " <SnapNum> <parttype bitmask> [output path]\n",
          RESTART_POWER_SPECTRUM);
      printf("      %2d          Write out the Voronoi mesh: <SnapNum>\n", RESTART_VORONOI_MESH);
      printf(
          "      %2d          Run the post-processing shock finder:"
          " <SnapNum> [<SubboxNum> for SUBBOX_SNAPSHOTS]\n",
          RESTART_SHOCK_FINDER);
      printf(
          "      %2d          Write out a two-dimensional slice of the Voronoi mesh:"
          " <SnapNum> <center_x> <center_y> <center_z> <normal_x> <normal_y> <normal_z>\n",
          RESTART_VORONOI_MESH_SLICE);
      printf("      %2d          Write out snapshot dump with measured gradients\n", RESTART_GRADIENTS);
      printf(
          "      %2d          Recalculate gravitational potential values for"
          " specified snaphot dump: <snapnum>\n",
          RESTART_RECALC_POTENTIAL);
      printf(
          "      %2d          Calculate additional quantities from a snapshot dump:"
          " <snapnum>\n",
          RESTART_CALC_ADDITIONAL);
      printf(
          "      %2d          Render Auriga movie frame from a snapshot dump:"
          " <snapnum>\n",
          RESTART_AURIGA_MOVIE);
      printf(
          "      %2d          Run SimpleX RT in post-processing mode on a snapshot:"
          " <snapnum> <runTime> <numSteps>\n",
          RESTART_SIMPLEX);
      printf("      %2d          Compute Voronoi densities for DM particles for specified snapshot dump: <snapnum>\n",
             RESTART_CALC_VORONOI_DM_DENSITY);
      printf("\n");
    }
}

/*! \brief Print information about the command line arguments.
 *
 *  \return void
 */
static void print_command_line_arguments(const int argc, char **argv, const char *const parameter_file, const int restart_flag,
                                         const int restart_snap_num)
{
  if(ThisTask == 0)
    {
      mpi_printf("Command line arguments:\n  ");
      for(int i = 0; i < argc; i++)
        mpi_printf("<%s> ", argv[i]);
      mpi_printf("\nParameterFile='%s', RestartFlag=%d", parameter_file, restart_flag);
      if(restart_snap_num != -1)
        mpi_printf(", RestartSnapNum=%d", restart_snap_num);
      mpi_printf("\n\n");
    }
}

/*! \brief This function ends the simulation in case of no error.
 *
 *  This method has to be called by all processes. It should be used only
 *  if the simulation ends without errors.
 *  Otherwise, terminate() should be used instead.
 *
 *  \return void
 */
void endrun(void) { endrun_with_status(EXIT_SUCCESS); }

/*! \brief Like endrun(), but with an argument to set the exit status for exit().
 *
 *  \return void
 */
static void endrun_with_status(const int status)
{
  mpi_printf("Code run for %f seconds!\n", timediff(StartOfRun, second()));
  mpi_printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

#ifdef HAVE_HDF5
  /* The HDF5 library will sometimes register an atexit() handler that calls
   * the error handler. This will cause an error message if the versions of the
   * functions H5Eset_auto() (called in AREPO) and H5Eget_auto() (used in the
   * exit handler) don’t match (e.g. H5Eset_auto1() and H5Eget_auto2()).
   * Hence, unset the HDF5 error handler here. */
  H5Eset_auto(NULL, NULL);
#endif

#ifdef CUDA
  cuda_finish();
#endif

#ifdef SNE_FEEDBACK
  sne_destroy();
#endif

#ifdef RUNNING_SAFETY_FILE
  /* simulation has not reached the final time */
  if(All.Ti_Current < TIMEBASE)
    {
      char running_fname[MAXLEN_PATH], running_done_fname[MAXLEN_PATH];
      sprintf(running_fname, "./running");
      sprintf(running_done_fname, "./running_done");
      rename(running_fname, running_done_fname);
      mpi_printf("moved ./running file to ./running_done, job can now restart.\n");
    }
  else
    mpi_printf("leaving ./running file in place since run is complete to prevent any restarts.\n");
#endif

  MPI_Finalize();
  exit(status);
}
