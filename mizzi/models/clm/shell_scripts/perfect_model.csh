#!/bin/tcsh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# The FORCE options are not optional.
# the VERBOSE options are useful for debugging.
set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.clm2_${ensemble_member}.r.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `head -1 ../rpointer.lnd`
set FILE = $FILE:t
set FILE = $FILE:r
set MODEL_DATE_EXT = `echo $FILE:e`
set MODEL_DATE_STR = `echo $FILE:e | sed -e "s#-# #g"`
set MODEL_DATE = `echo $MODEL_DATE_STR`
@ MODEL_YEAR    = $MODEL_DATE[1]
@ MODEL_MONTH   = $MODEL_DATE[2]
@ MODEL_DAY     = $MODEL_DATE[3]
@ MODEL_SECONDS = $MODEL_DATE[4]

echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_SECONDS"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = ${HOME}/DART/models/clm/work
set DARTDIR = /glade/home/thoar/clm/models/clm/work
set  OBSDIR = /glade/scratch/thoar/CLM_PMO_leafc

#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART
# Either get them from the CCSM 'run' directory or some stock repository
# The grid files are absolute paths ... so they need not move.
#-------------------------------------------------------------------------

foreach FILE ( input.nml perfect_model_obs clm_to_dart dart_to_clm )

   ${COPY} ${DARTDIR}/${FILE} . || exit 1

end

#-------------------------------------------------------------------------
# Block 1: convert 1 clm restart file to a DART initial conditions file.
# At the end of the block, we have a DART restart file  perfect_ics
# that came from the pointer file ../rpointer.lnd_0001
#
# DART namelist settings appropriate/required:
# &perfect_model_obs_nml:  restart_in_file_name    = 'perfect_ics'
# &clm_to_dart_nml:        clm_to_dart_output_file = 'dart.ud'
#-------------------------------------------------------------------------

set member = 1

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set POINTER_FILENAME = rpointer.lnd
   set MODEL_RESTART_FILENAME = `head -1 ../../${POINTER_FILENAME}`
   set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s/\.r\./\.h0\./"`
   ${LINK} ../../$MODEL_RESTART_FILENAME clm_restart.nc
   ${LINK} ../../$MODEL_HISTORY_FILENAME clm_history.nc
   ${COPY} ../input.nml                  .

   # patch the CLM restart files to ensure they have the proper
   # _FillValue and missing_value attributes.
   ncatted -O -a    _FillValue,frac_sno,o,d,1.0e+36   clm_restart.nc
   ncatted -O -a missing_value,frac_sno,o,d,1.0e+36   clm_restart.nc
   ncatted -O -a    _FillValue,DZSNO,o,d,1.0e+36      clm_restart.nc
   ncatted -O -a missing_value,DZSNO,o,d,1.0e+36      clm_restart.nc
   ncatted -O -a    _FillValue,H2OSOI_LIQ,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a missing_value,H2OSOI_LIQ,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a    _FillValue,H2OSOI_ICE,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a missing_value,H2OSOI_ICE,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a    _FillValue,T_SOISNO,o,d,1.0e+36   clm_restart.nc
   ncatted -O -a missing_value,T_SOISNO,o,d,1.0e+36   clm_restart.nc

   echo "starting clm_to_dart for member ${member} at "`date`
   ../clm_to_dart >! output.${member}.clm_to_dart
   echo "finished clm_to_dart for member ${member} at "`date`

   mv dart.ud ../perfect_ics
   cd ..

#-------------------------------------------------------------------------
# Block 2: Advance the model and harvest the synthetic observations.
# Will result in a single file : 'perfect_restart' which we don't need
# for a perfect model experiment with CESM.
#
# DART namelist settings required:
# &perfect_model_obs_nml:           async                  = 0,
# &perfect_model_obs_nml:           adv_ens_command        = "./no_model_advance.csh",
# &perfect_model_obs_nml:           restart_in_file_name   = 'perfect_ics'
# &perfect_model_obs_nml:           restart_out_file_name  = 'perfect_restart'
# &perfect_model_obs_nml:           obs_sequence_in_name   = 'obs_seq.in'
# &perfect_model_obs_nml:           obs_sequence_out_name  = 'obs_seq.out'
# &perfect_model_obs_nml:           init_time_days         = -1,
# &perfect_model_obs_nml:           init_time_seconds      = -1,
# &perfect_model_obs_nml:           first_obs_days         = -1,
# &perfect_model_obs_nml:           first_obs_seconds      = -1,
# &perfect_model_obs_nml:           last_obs_days          = -1,
# &perfect_model_obs_nml:           last_obs_seconds       = -1,
#
#-------------------------------------------------------------------------

# clm always needs a clm_restart.nc, and a clm_history.nc to start.

set MODEL_RESTART_FILENAME = `head -1 ../rpointer.lnd`
set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s/\.r\./\.h0\./"`

${LINK} ../$MODEL_RESTART_FILENAME clm_restart.nc
${LINK} ../$MODEL_HISTORY_FILENAME clm_history.nc

# Determine proper observation sequence file.

set OBSFNAME = `printf obs_seq.in.0Z.%04d%02d%02d ${MODEL_YEAR} ${MODEL_MONTH} ${MODEL_DAY}`
set OBS_FILE = ${OBSDIR}/${MODEL_YEAR}/${OBSFNAME}

${LINK} ${OBS_FILE}   obs_seq.in

./perfect_model_obs || exit 2

${MOVE} True_State.nc    ../True_State.${MODEL_DATE_EXT}.nc
${MOVE} obs_seq.out      ../obs_seq.${MODEL_DATE_EXT}.out
${MOVE} dart_log.out     ../dart_log.${MODEL_DATE_EXT}.out

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

