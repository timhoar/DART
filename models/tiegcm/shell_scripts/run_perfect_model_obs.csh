#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to generate observations and a TRUE state.
#
# Unlike the more complex job.csh, this script only processes a single 
# observation file.  Still fairly complex; requires a raft of
# data files and most of them are in hardcoded locations.
#
# This script is designed to be run from the command line (as a single thread)
# and should only take a few seconds to a minute to complete, depending on
# the filesystem performance and data file size.
#
# The script moves the necessary files to the current directory - in DART
# nomenclature, this will be called CENTRALDIR. 
# After everything is confirmed to have been assembled, it is possible
# to edit the data, data.cal, and input.nml files for the specifics of 
# the experiment; as well as allow final configuration of a 'nodelist' file.
#
# Once the 'table is set', all that remains is to start/submit the 
# 'runme_filter' script. That script will spawn 'filter' as a 
# parallel job on the appropriate nodes; each of these tasks will 
# call a separate advance_model.csh when necessary.
#
# The central directory is where the scripts reside and where script and 
# program I/O are expected to happen.
#-----------------------------------------------------------------------------
#
#BSUB -J tiegcm_perfect
#BSUB -o tiegcm_perfect.%J.log
#BSUB -P P3507xxxx
#BSUB -q regular
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_HOSTS) then

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     tiegcm_perfect
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST

endif

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $MYHOST"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------

setenv CENTRALDIR /glade/scratch/${user}/DART/${JOBNAME}/job_${JOBID}

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

set myname = $0          # this is the name of this script

# some systems don't like the -v option to any of the following 

set OSTYPE = `uname -s` 
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -v'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
# DARTDIR      The location of the DART tiegcm model directory
# TIEGCMDIR    The location of the TIEGCM executable
# EXPERIMENT   The location of the initial ensemble of TIEGCM files
#-----------------------------------------------------------------------------

set    DARTDIR = /glade/u/home/${USER}/work/DART/tiegcm/models/tiegcm
set  TIEGCMDIR = /glade/u/home/${USER}/work/DART/tiegcm/models/tiegcm/src
set EXPERIMENT = /glade/p/work/${USER}/initial_ensemble

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the tiegcm executable, control files, and data files.
# The tiegcm initial conditions are in the next block.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/perfect_model_obs          . || exit 1
${COPY} ${DARTDIR}/work/dart_to_model              . || exit 1
${COPY} ${DARTDIR}/work/model_to_dart              . || exit 1
${COPY} ${DARTDIR}/work/input.nml                  . || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh . || exit 1
${COPY} ${DARTDIR}/work/obs_seq.in                 . || exit 1

${COPY} ${TIEGCMDIR}/tiegcm-nompi             tiegcm || exit 1

${COPY} ${TIEGCMDIR}/tiegcm_restart_p.nc           . || exit 1
${COPY} ${TIEGCMDIR}/tiegcm_s.nc                   . || exit 1
${COPY} ${TIEGCMDIR}/tiegcm.nml  tiegcm.nml.original || exit 1

#-----------------------------------------------------------------------------
# Remove all the comments that follow (;) symbol from tiegcm.nml namelist file
# That is a non-standard syntax for fortran namelists.
#-----------------------------------------------------------------------------

grep -v "^;" tiegcm.nml.original >! tiegcm.nml  || exit 1

#-----------------------------------------------------------------------------
# Convert a TIEGCM file 'tiegcm_restart.nc' to a DART ics file 'perfect_ics'
# 
# REQUIREMENTS: for input.nml
# model_nml            : tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
# model_nml            : tiegcm_secondary_file_name = 'tiegcm_s.nc'
# model_nml            : tiegcm_namelist_file_name  = 'tiegcm.nml'
# model_to_dart_nml    : file_out                   = 'dart_ics'
#-----------------------------------------------------------------------------
# perfect_model_obs_nml: async                      = 2
# perfect_model_obs_nml: adv_ens_command            = 'advance_model.csh'
# perfect_model_obs_nml: start_from_restart         = .TRUE.
# perfect_model_obs_nml: restart_in_file_name       = 'perfect_ics'
#-----------------------------------------------------------------------------
# dart_to_model_nml    : file_in                    = 'dart_restart'
# dart_to_model_nml    : file_namelist_out          = 'namelist_update'

./model_to_dart              || exit 2
${MOVE} dart_ics perfect_ics || exit 2

#-----------------------------------------------------------------------------
# Run perfect_model_obs ... harvest the observations to populate obs_seq.out
# model_mod expects a generic name // advance_model.csh expects a filename
# with the ensemble member ID tacked on - must provide both.
#-----------------------------------------------------------------------------

ln -sf tiegcm_restart_p.nc tiegcm_restart_p.nc.0001 || exit 3
ln -sf tiegcm_s.nc         tiegcm_s.nc.0001         || exit 3
ln -sf tiegcm.nml          tiegcm.nml.0001          || exit 3

./perfect_model_obs || exit 3

#-----------------------------------------------------------------------------
# At this point, all the restart,diagnostic files are in the run/CENTRALDIR.
# You may want to move them to someplace more 'permanent'.
#
# TJH: At this point, the output files have pretty 'generic' names.
# The files could be archived with the assimilation date in their name.
#-----------------------------------------------------------------------------

# ${MOVE} tiegcm_s.nc.0001           ${EXPERIMENT}/perfect/tiegcm_s.nc
# ${MOVE} tiegcm_restart_p.nc.0001   ${EXPERIMENT}/perfect/tiegcm_restart_p.nc
# ${MOVE} tiegcm.nml                 ${EXPERIMENT}/perfect
# ${MOVE} obs_seq.out                ${EXPERIMENT}/perfect
# ${MOVE} True_State.nc              ${EXPERIMENT}/perfect

# ${MOVE} tiegcm_out_1               ${EXPERIMENT}/perfect/tiegcm_out
# ${MOVE} dart_log.out               ${EXPERIMENT}/perfect
# ${MOVE} dart_log.nml               ${EXPERIMENT}/perfect
# Good style dictates that you save the scripts so you can see what worked.

# ${COPY} input.nml                  ${EXPERIMENT}/DART
# ${COPY} *.csh                      ${EXPERIMENT}/DART
# ${COPY} $myname                    ${EXPERIMENT}/DART

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "These are the files in the run directory at completion:"
ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

