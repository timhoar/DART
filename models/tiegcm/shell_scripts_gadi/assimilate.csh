#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# script to perform a single assimilation. filter does not advance the model.
# The observation sequence file is chopped into the exact timeframe needed for
# the assimilation.
#
# This script is designed to be submitted as a batch job but may be run from 
# the command line (as a single thread) to check for file motion, etc.
# If running interactively, please comment out the part that actually runs filter.
#
#-------------------------------------------------------------------------------
#
#PBS -P n23
#PBS -l walltime=0:05:00
#PBS -l wd
#PBS -l ncpus=16
#PBS -l mem=320000GB
#PBS -N assim
#PBS -m ae
#PBS -M g.bowden@adfa.edu.au
#PBS -j oe

#===============================================================================

if ($?PBS_O_HOST) then

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $HOSTNAME
   setenv MPI_RUN_CMD mpirun

   # MP_DEBUG_NOTIMEOUT may alleviate MPI timeouts that may occur under
   # certain task geometries. It is NOT a good idea to use it in general. 
   # setenv MP_DEBUG_NOTIMEOUT yes

else

   #----------------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #----------------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     tiegcm_filter
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI_RUN_CMD ''

endif

#-------------------------------------------------------------------------------
# Just an echo of job attributes
#-------------------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $MYHOST"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#-------------------------------------------------------------------------------
# STEP 1: Set the environment (modules, variables, etc.) for this experiment.
#-------------------------------------------------------------------------------

# This string gets replaced by stage_experiement when it get copies into place.
# We need to be in the right directory before we can source DART_params.csh
cd CENTRALDIRSTRING

if ( -e DART_params.csh ) then
   source DART_params.csh
else
   echo "ERROR: resource file 'DART_params.csh' not found."
   echo "       need one in "`pwd`
   exit 1
endif

# ------------------------------------------------------------------------------
# STEP 2: Clear out previous results to avoid using stale information
# ------------------------------------------------------------------------------

${REMOVE} obs_seq.out obs_seq.final
${REMOVE} dart_log.out dart_log.nml

\unlink tiegcm_restart_p.nc tiegcm_s.nc tiegcm.nml 

${LINK} tiegcm_restart_p.nc.0001 tiegcm_restart_p.nc   || exit 3
${LINK} tiegcm_s.nc.0001         tiegcm_s.nc           || exit 3
${LINK} tiegcm.nml.0001          tiegcm.nml            || exit 3

# ------------------------------------------------------------------------------
# STEP 3: Get the right observation sequence file for this date/time
# ------------------------------------------------------------------------------

set mtime = `ncdump -v mtime tiegcm_restart_p.nc | tail -n 2 | head -n 1`

set YEAR   = 2018
set DOY    = `echo $mtime[1] | sed -e "s/,//"`
set HOUR   = `echo $mtime[2] | sed -e "s/,//"`
set MINUTE = `echo $mtime[3] | sed -e "s/,//"`

echo "MODEL TIME read as $YEAR $DOY $HOUR $MINUTE"

set TIMESTAMP = `printf %02d%02d%02d $HOUR $MINUTE $SECOND`

# Get the appropriate set of observations for this timeframe. 

if ( -e $OBSDIR/obs_seq.$TIMESTAMP.out ) then
  ${LINK} $OBSDIR/obs_seq.$TIMESTAMP.out obs_seq.out || exit 3
else
   # could run obs_sequence_tool to cut out just the obs of interest
   echo "ERROR: no observations for $YEAR $DOY $HOUR $MINUTE"
   echo "ERROR: expecting a file: $OBSDIR/obs_seq.$TIMESTAMP.out"
   exit 3
endif

# ------------------------------------------------------------------------------
# STEP 4: [OPTIONAL] prepare for DART INFLATION 
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(:) /= 0 AND inf_initial_from_restart = .TRUE.
#
# The strategy is to use the LATEST (last) inflation file, which is
# always available in the contents of the pointer file.
# ------------------------------------------------------------------------------
#
# Only supporting prior state-space inflation for now.

# If the inflation pointer file does not exist, we need to make one.
# This should only happen on the first assimilation cycle.

if (! -e priorinf_pointer_file.txt ) then
   if ( -e prior_inflate_ics ) then # run fill_inflation_restart.f90 
      # ./fill_inflation_restart
      $MOVE prior_inflate_ics prior_inflate_ics.$TIMESTAMP
      echo "prior_inflate_ics.$TIMESTAMP" >! priorinf_pointer_file.txt
   else
      echo "ERROR: inflation staging failed."
      exit 4
   endif
endif

# Link the inflation file (output_priorinf_[mean,sd].DATE.nc) from
# the previous assimilation cycle to be used as input for the current
# assimilation. These are simply read, so linking is sufficient.

set inflation_file = `head -n 1 priorinf_pointer_file.txt`

${LINK} -f $inflation_file prior_inflate_ics

# ------------------------------------------------------------------------------
# STEP 5: assimilate the ensemble of new states.
#
# The stage_experiment script linked all the TIE-GCM files in the instance
# directories to the expected names required by filter. Should be good to go.
# ------------------------------------------------------------------------------

${MPI_RUN_CMD} ./filter || exit 3

if ( ! -e obs_seq.final ) then
   echo "ERROR: filter did not finish successfully - no obs_seq.final"
   exit 5
endif

# ------------------------------------------------------------------------------
# STEP 6: Tag the output with the valid time of the model state.
# Update the inflation pointer file with the most recent file.
# ------------------------------------------------------------------------------

${MOVE} -v obs_seq.final     obs_seq.final.$TIMESTAMP
${MOVE} -v namelist_update namelist_update.$TIMESTAMP

if ( -e prior_inflate_restart ) then
   ${MOVE} -v prior_inflate_restart prior_inflate_restart.$TIMESTAMP
   echo "prior_inflate_restart.$TIMESTAMP" >! priorinf_pointer_file.txt
endif

# Must copy each filter_restart.nnnn to the right place, update link
# so that 'dart_to_model' gets what it needs.
foreach RESTART ( filter_restart.???? )
   set instance = $RESTART:e
   set INSTANCE_DIRECTORY = `printf "instance_%04d" $instance`
   cd $INSTANCE_DIRECTORY
   ${MOVE} -v ../$RESTART  $RESTART.$TIMESTAMP
   \unlink dart_restart
   ${LINK} $RESTART.$TIMESTAMP dart_restart
   cd ..
end

# ------------------------------------------------------------------------------
# STEP 7: Finish.
# ------------------------------------------------------------------------------

echo
echo "filter completed successfully at "`date`

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "These are the files in the run directory at completion:"
ls -lrt

exit 0

