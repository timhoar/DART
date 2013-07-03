#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Script to manage the compilation of all components for this model;
# executes a known "perfect model" experiment using an existing
# observation sequence file (obs_seq.in) and initial conditions appropriate 
# for both 'perfect_model_obs' (perfect_ics) and 'filter' (filter_ics).
# There are enough initial conditions for 80 ensemble members in filter.
# Use ens_size = 81 and it WILL bomb. Guaranteed.
# The 'input.nml' file controls all facets of this execution.
#
# 'create_obs_sequence' and 'create_fixed_network_sequence' were used to
# create the observation sequence file 'obs_seq.in' - this defines 
# what/where/when we want observations. This script does not run these 
# programs - intentionally. 
#
# 'perfect_model_obs' results in a True_State.nc file that contains 
# the true state, and obs_seq.out - a file that contains the "observations"
# that will be assimilated by 'filter'.
#
# 'filter' results in three files (at least): Prior_Diag.nc - the state 
# of all ensemble members prior to the assimilation (i.e. the forecast), 
# Posterior_Diag.nc - the state of all ensemble members after the 
# assimilation (i.e. the analysis), and obs_seq.final - the ensemble 
# members' estimate of what the observations should have been.
#
# Once 'perfect_model_obs' has advanced the model and harvested the 
# observations for the assimilation experiment, 'filter' may be run 
# over and over by simply changing the namelist parameters in input.nml.
#
# The result of each assimilation can be explored in model-space with
# matlab scripts that directly read the netCDF output, or in observation-space.
# 'obs_diag' is a program that will create observation-space diagnostics
# for any result of 'filter' and results in a couple data files that can
# be explored with yet more matlab scripts.
#
#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the 
# resulting source file is used by all the remaining programs, 
# so this MUST be run first.
#----------------------------------------------------------------------
#
# If you get a ton of compile ERRORS (not warnings) read on ...
#
# Since a lot of the code is inherited from wrf, it comes with a .F 
# extension even though it is free-format. This makes it necessary to
# compile with flags that force interpretation of free-format.
# One way around this is to rename the files and the references in 
# the path_names_xxxxx files - which is done by ChangeExtensions.csh
# in the PBL_1d/shell_scripts directory, or try to hunt down the right
# flag for the compiler you are using. Your choice.
#
# The code also relies on the autopromotion flag that coerces all 
# 'real' variables to be 8byte real variables.
#
# Intel     -free -r8
# gfortran  -ffree-form -fdefault-real-8
# pathscale -freeform -r8
# pgi       -Mfree -Mr8
# absoft    -ffree  (see the mkmf.template for absoft for more on r8)
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set MODEL = "PBL_1d"

@ n = 1

echo
echo
echo "---------------------------------------------------------------"
echo "${MODEL} build number ${n} is preprocess"

csh  mkmf_preprocess
make || exit $n

./preprocess || exit 99

#----------------------------------------------------------------------
# Build all the single-threaded targets
#----------------------------------------------------------------------
echo
echo "Building this model requires the real*8 override flag be added to the"
echo "default mkmf.template rules. If the following compile fails read the"
echo "comments in the quickbuild.csh script for more help."
echo

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   switch ( $TARGET )
   case mkmf_preprocess:
      breaksw
   default:
      @ n = $n + 1
      echo
      echo "---------------------------------------------------"
      echo "${MODEL} build number ${n} is ${PROG}" 
      \rm -f ${PROG}
      csh $TARGET || exit $n
      make        || exit $n
      breaksw
   endsw
end

\rm -f *.o *.mod input.nml*_default

if ( $#argv == 1 && "$1" == "-mpi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script now compiling MPI parallel versions of the DART programs."
else if ( $#argv == 1 && "$1" == "-nompi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script is exiting without building the MPI version of the DART programs."
  exit 0
else
  echo ""
  echo "Success: All DART programs compiled."
  echo "Script is exiting before building the MPI version of the DART programs."
  echo "Run the quickbuild.csh script with a -mpi argument or"
  echo "edit the quickbuild.csh script and remove the exit line"
  echo "to compile with MPI to run in parallel on multiple cpus."
  echo ""
  exit 0
endif

#----------------------------------------------------------------------
# to enable an MPI parallel version of filter for this model, 
# call this script with the -mpi argument, or if you are going to build
# with MPI all the time, remove or comment out the entire section above.
#----------------------------------------------------------------------

\rm -f filter wakeup_filter

@ n = $n + 1
echo
echo "---------------------------------------------------"
echo "build number $n is mkmf_filter"
csh   mkmf_filter -mpi
make

if ($status != 0) then
   echo
   echo "If this died in mpi_utilities_mod, see code comment"
   echo "in mpi_utilities_mod.f90 starting with 'BUILD TIP' "
   echo
   exit $n
endif

@ n = $n + 1
echo
echo "---------------------------------------------------"
echo "build number $n is mkmf_wakeup_filter"
csh  mkmf_wakeup_filter -mpi
make || exit $n

\rm -f *.o *.mod input.nml*_default

echo
echo 'time to run filter here:'
echo ' for lsf run "bsub < runme_filter"'
echo ' for pbs run "qsub runme_filter"'
echo ' for lam-mpi run "lamboot" once, then "runme_filter"'
echo ' for mpich run "mpd" once, then "runme_filter"'

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
