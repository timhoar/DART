#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
#-------------------------------------------------------------------------------
# This block of directives constitutes the preamble for the PBS queuing system
#
# the normal way to submit to the queue is:    qsub run_filter.csh
#
#PBS -P n23
#PBS -l walltime=0:05:00
#PBS -l wd
#PBS -l ncpus=16
#PBS -l mem=320000GB
#PBS -N tiegcm_advance
#PBS -m ae
#PBS -M g.bowden@adfa.edu.au
#PBS -j oe
#PBS -J 1-ENSEMBLESIZESTRING

#===============================================================================

if ($?PBS_O_HOST) then

   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv INSTANCE    $PBS_ARRAY_INDEX
   set    NODELIST = `cat "${PBS_NODEFILE}"`
   set    nproc    = `cat ${PBS_NODEFILE} | wc -l`

else
   echo "This script must be run with a queueing system that supports job arrays."
   exit 1
endif

#-------------------------------------------------------------------------------
# Just an echo of job attributes
#-------------------------------------------------------------------------------

echo
echo "# Advancing TIE-GCM instance $INSTANCE started at "`date`
echo
echo "${JOBNAME} ($JOBID) working directory "`pwd`
echo "${JOBNAME} ($JOBID) running on ${NODELIST}"
echo "${JOBNAME} ($JOBID) started at "`date`
echo

#-------------------------------------------------------------------------------
# STEP 1: Set the environment (modules, variables, etc.) for this experiment.
#-------------------------------------------------------------------------------

# This string gets replaced by stage_experiment when it get copies into place.
# We need to be in the right directory before we can source DART_params.csh
cd CENTRALDIRSTRING

if ( -e DART_params.csh ) then
   source DART_params.csh
else
   echo "ERROR: resource file 'DART_params.csh' not found."
   echo "       need one in "`pwd`
   exit 1
endif

#-------------------------------------------------------------------------------
# STEP 2: Move to the directory specifi for this instance, prepare.
#-------------------------------------------------------------------------------

set INSTANCE_DIRECTORY = `printf instance_%04d $INSTANCE`
cd ${INSTANCE_DIRECTORY}

# Make sure we have a clean logfile for this entire advance
# and set the filename for the DART prior for this member.
set     logfile = `printf log_advance.%04d.txt $INSTANCE`
set output_file = `printf filter_ics.%04d $INSTANCE`

rm -f $logfile
rm -f input.nml
rm -f tiegcm.nml.updated
rm -f dart_ics

# Ensure that the input.nml has the required value for
# dart_to_model_nml:advance_time_present for this context.

sed -e "/advance_time_present /c\ advance_time_present = .TRUE." \
       ../input.nml >! input.nml || exit 2

set mpirun = `cat ../mpirun.command`
set RUN_CMD = "$mpirun -np $nproc"

#-------------------------------------------------------------------------------
# STEP 3: Convert the DART output file to form needed by model.
# Overwrite the appropriate variables of a TIEGCM netCDF restart file.
# The DART output file (namelist_update) has the 'advance_to' time 
# which must be communicated to the model ... through the tiegcm namelist
#-------------------------------------------------------------------------------

ls -l dart_restart        >>&   $logfile || exit 3
ls -l tiegcm_s.nc         >>&   $logfile || exit 3
ls -l tiegcm_restart_p.nc >>&   $logfile || exit 3
ls -l tiegcm.nml          >>&   $logfile || exit 3

# dart_to_model also generates 'namelist_update'
../dart_to_model >>& $logfile || exit 3

# update tiegcm namelist variables by grabbing the values from namelist_update  
# and then overwriting whatever is in the tiegcm.nml
# There is a danger that you match multiple things ... F107 and F107A,
# SOURCE_START and START, for example ... so try to grep whitespace too ...

set start_year   = `grep " START_YEAR "   namelist_update`
set start_day    = `grep " START_DAY "    namelist_update`
set source_start = `grep " SOURCE_START " namelist_update`
set start        = `grep " START "        namelist_update`
set secstart     = `grep " SECSTART "     namelist_update`
set stop         = `grep " STOP "         namelist_update`
set secstop      = `grep " SECSTOP "      namelist_update`
set hist         = `grep " HIST "         namelist_update`
set sechist      = `grep " SECHIST "      namelist_update`
set f107         = `grep " F107 "         namelist_update`

# FIXME TOMOKO ... do we want F107 and F107A to be identical
#
# the way to think about the following sed syntax is this:
# / SearchStringWithWhiteSpaceToMakeUnique  /c\ the_new_contents_of_the_line 

sed -e "/ START_YEAR /c\ ${start_year}" \
    -e "/ START_DAY /c\ ${start_day}" \
    -e "/ SOURCE_START /c\ ${source_start}" \
    -e "/ START /c\ ${start}" \
    -e "/ STOP /c\ ${stop}" \
    -e "/ HIST /c\ ${hist}" \
    -e "/ SECSTART /c\ ${secstart}" \
    -e "/ SECSTOP /c\ ${secstop}" \
    -e "/ SECHIST /c\ ${sechist}" \
    -e "/ F107 /c\ ${f107}" \
    tiegcm.nml >! tiegcm.nml.updated

if ( -e tiegcm.nml.updated ) then
   echo "tiegcm.nml updated with new start/stop time for ensemble member $INSTANCE"
   mv -v tiegcm.nml tiegcm.nml.original
   mv -v tiegcm.nml.updated tiegcm.nml
else
   echo "ERROR tiegcm.nml did not update correctly for ensemble member $INSTANCE."
   exit 3
endif

#-------------------------------------------------------------------------------
# STEP 4: Run the model
#-------------------------------------------------------------------------------

echo "ensemble member $INSTANCE : before tiegcm" >>  $logfile
ncdump -v mtime tiegcm_restart_p.nc              >>& $logfile
echo "Starting tiegcm at "`date`                 >>  $logfile

${RUN_CMD} ../tiegcm tiegcm.nml >>& $logfile

grep -q "NORMAL EXIT" $logfile
set tiegcm_status = $status

if ($tiegcm_status != 0) then
   echo "ERROR: tiegcm model advance failed."
   echo "ERROR: check $INSTANCE_DIRECTORY/$logfile"
   exit 4
endif

echo "ensemble member $INSTANCE : after tiegcm" >>  $logfile
ncdump -v mtime tiegcm_restart_p.nc             >>& $logfile

# grab the bits needed for unique timestamp
set YEAR = $start_year[3]
set  DOY = $stop[3]
set HOUR = $stop[5]
set TIMESTAMP = ${YEAR}.${DOY}.${HOUR}

#-------------------------------------------------------------------------------
# STEP 5: Convert the model output to form needed by DART
# AT this point, the model has updated the information in tiegcm_restart_p.nc
# We need to get that information back into the DART state vector.
#
# The updated information needs to be moved into CENTRALDIR in
# preparation for the next cycle.
#-------------------------------------------------------------------------------

echo "Starting model_to_dart at "`date`  >>  $logfile
../model_to_dart                         >>& $logfile

if ( -e dart_ics ) then
   ${COPY} dart_ics   dart_priors.${TIMESTAMP} || exit 5
   ${MOVE} dart_ics   ../$output_file          || exit 5
   echo
   echo "Finished advance_tiegcm.csh for member $INSTANCE at "`date`

else
   echo
   echo "ERROR: model_to_dart failed for member $INSTANCE"
   exit 5
endif

exit 0
