#!/bin/sh
# shell script for code testing

#################################
## perform predefined examples ##
#################################

# Number of cores to compile and run the problem on
# (Some examples will override this setting, e.g. 1d test problems can only
#  run on 1 task!)
NUMBER_OF_TASKS=8
NUMBER_OF_COMPILERS=8
PLOT=False  # create plots for runs? True/False
PYTHON=python3

## choose your tests
TESTS=''
## available 1d test cases
TESTS+='wave_1d '
TESTS+='shocktube/shocktube_1d '
TESTS+='shocktube/shocktube_sod_1d '
TESTS+='interacting_blastwaves_1d '
TESTS+='polytrope_spherical_1d '
TESTS+='mhd_shocktube_1d '
TESTS+='alfven_wave_1d '

## available 2d test cases
TESTS+='Gresho_2d '
TESTS+='Noh_2d/Noh_2d '
TESTS+='Noh_2d/Noh_refinement_2d '
TESTS+='Yee/Yee_2d '
TESTS+='Yee/Yee_Cartesian_2d '
TESTS+='khi_hydro_2d '
TESTS+='khi_mhd_2d '
TESTS+='brag_2d '
TESTS+='fastwave/fastwave_2d '
TESTS+='fastwave/fastwave_local_timestepping_2d '
TESTS+='fastwave/fastwave_subcycling_2d '
TESTS+='fastwave/fastwave_RKL2_timestepping_2d '
TESTS+='fastwave/fastwave_RKL2_plus_local_timestepping_2d '
TESTS+='rising_bubble_2d '
TESTS+='current_sheet_2d '
TESTS+='OrszagTang/OrszagTang_2d '
TESTS+='OrszagTang/OrszagTang_Boost_2d '
TESTS+='OrszagTang/OrszagTang_LabFrameExtrapolation_2d '

## 3d test cases
TESTS+='Noh_3d/Noh_3d '
TESTS+='Noh_3d/Noh_refinement_3d '
TESTS+='cosmo_box/gravity_only_3d '
TESTS+='cosmo_box/star_formation_3d '
#TESTS+='cosmo_zoom_gravity_only_3d '
TESTS+='isolated_galaxy_collisionless_3d '
TESTS+='galaxy_merger_star_formation_3d '
TESTS+='brag_decay_3d '

## AMR test cases
TESTS+='AMR/shocktube_2d '

## compilations tests
TESTS+='sanity_checks/compilation_basic '
TESTS+='sanity_checks/compilation_examples/cosmo_zoom_gravity_only_3d_parent '


# clean up
rm -rf run/

# set variables
export NUMBER_OF_TASKS
export NUMBER_OF_COMPILERS
DIR_BASE='examples'
RUNDIR_BASE='run/examples'

# loop over all tests
for TEST in $TESTS; do
  DIR="$DIR_BASE/$TEST"
  RUNDIR="$RUNDIR_BASE/$TEST"

  # create run directory
  mkdir -p "$RUNDIR/"

  # copy setup to run directory
  cp -RL "$DIR/." "$RUNDIR/"

  # create ICs in run directory
  echo "$DIR/"
  "$PYTHON" "$DIR/create.py" "$RUNDIR/"
  return_value="$?"      # get return value
  if [ "$return_value" != 0 ]; then    # check return value
    printf 'ERROR: test.sh:\t%s\t %s create.py failed!\n' "$DIR/" "$PYTHON"
    exit $return_value
  fi

  # compile Arepo
  make -j "$NUMBER_OF_COMPILERS" CONFIG="$RUNDIR/Config.sh" BUILD_DIR="$RUNDIR/build" EXEC="$RUNDIR/Arepo"
  return_value="$?"      # get return value
  if [ "$return_value" != 0 ]; then    # check return value
    printf 'ERROR: test.sh:\t%s\t make failed!\n' "$DIR/"
    exit $return_value
  fi

  # do not run executable if the test is compilation-only
  if ! grep -q '^\s*COMPILE_ONLY\s*=\s*True\b' "$DIR/create.py"; then
    # check if the test wants a specific number of MPI tasks
    num_tasks_line=$(grep -o '^\s*numTasksMPI\s*=\s*[1-9][0-9]*\b' \
      "$DIR/create.py")
    num_tasks_override=$(echo "$num_tasks_line" | grep -o '[1-9][0-9]*')
    if [ -n "$num_tasks_override" ]; then
      NTASK="$num_tasks_override"
    else
      NTASK="$NUMBER_OF_TASKS"
    fi
    # change to RUNDIR in subshell and execute test simulation
    (cd "$RUNDIR/" && mpiexec -n "$NTASK" ./Arepo param.txt)
    return_value="$?"    # get return value
    if [ "$return_value" != 0 ]; then  # check return value
      printf 'ERROR: test.sh:\t%s\t execution failed!\n' "$DIR/"
      exit $return_value
    fi

    # check result in example directory; this also creates some check plots
    $PYTHON "$DIR/check.py" "$RUNDIR/" "$PLOT"
    return_value="$?"    # get return value
    if [ "$return_value" != 0 ]; then  # check return value
      echo 'ERROR: test.sh: test failed!'
      exit $return_value
    fi
  fi

  if [ "$return_value" == 0 ]; then    # check return value
    printf 'test.sh:\t%s\t test passed!\n\n' "$DIR/"
  fi

  # clean up
  if [ -z "$PLOT" ] || [ "$PLOT" == 0 ]; then
    echo 'cleaning up...'
    rm -rf run/
  fi

  echo
done

echo
echo
echo 'Tests'
echo $TESTS
echo 'passed!'

exit $return_value
