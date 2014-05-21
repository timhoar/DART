#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by the model
# 3) runs the model
# 4) converts the model output to input expected by DART
#
# The error code from the script reflects which block it failed.
#
# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process

echo "starting advance_model"

set process = $1
set num_states = $2
set control_file = $3

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to run tiegcm.
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get the data files needed to run tiegcm. One directory up is 'CENTRALDIR'

cp ../input.nml .

# Ensure that the input.nml has the required value for
# dart_to_model_nml:advance_time_present for this context.

echo '1'                      >! ex_commands
echo '/dart_to_model_nml'     >> ex_commands
echo '/advance_time_present'  >> ex_commands
echo ':s/\.false\./\.true\./' >> ex_commands
echo ':wq'                    >> ex_commands

( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -n $ensemble_member_line ../$control_file | tail -n 1`
   set input_file      = `head -n $input_file_line      ../$control_file | tail -n 1`
   set output_file     = `head -n $output_file_line     ../$control_file | tail -n 1`

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by model.
   # Overwrite the appropriate variables of a TIEGCM netCDF restart file.
   # The DART output file (namelist_update) has the 'advance_to' time 
   # which must be communicated to the model ... through the tiegcm namelist
   #----------------------------------------------------------------------

   # The EXPECTED input DART 'output' file name is 'dart_restart'

   set tiesecond   = `printf "tiegcm_s.nc.%04d"         $ensemble_member`
   set tierestart  = `printf "tiegcm_restart_p.nc.%04d" $ensemble_member`
   set tieinp      = `printf "tiegcm.nml.%04d"          $ensemble_member`

   cp -p   ../$input_file dart_restart        || exit 2
   cp -p   ../$tiesecond  tiegcm_s.nc         || exit 2
   cp -p   ../$tierestart tiegcm_restart_p.nc || exit 2
   cp -p   ../$tieinp     tiegcm.nml          || exit 2

#  echo "ensemble member $ensemble_member : before dart_to_model"
#  ncdump -v mtime tiegcm_restart.nc

   ../dart_to_model || exit 2  # dart_to_model generates namelist_update

   # update tiegcm namelist variables   

   set start_year   = " START_YEAR = "`  head -n 1 namelist_update | tail -n 1`
   set start_day    = " START_DAY = "`   head -n 2 namelist_update | tail -n 1`
   set source_start = " SOURCE_START = "`head -n 3 namelist_update | tail -n 1`
   set start        = " START = "`       head -n 3 namelist_update | tail -n 1`
   set secstart     = " SECSTART = "`    head -n 3 namelist_update | tail -n 1`
   set stop         = " STOP = "`        head -n 4 namelist_update | tail -n 1`
   set secstop      = " SECSTOP = "`     head -n 4 namelist_update | tail -n 1`
   set hist         = " HIST = "`        head -n 5 namelist_update | tail -n 1`
   set sechist      = " SECHIST = "`     head -n 5 namelist_update | tail -n 1`
   set save         = " SAVE = "`        head -n 5 namelist_update | tail -n 1`
   set secsave      = " SECSAVE = "`     head -n 5 namelist_update | tail -n 1`
   set f107         = " F107 ="`         head -n 6 namelist_update | tail -n 1`

# FIXME TOMOKO - can I remove the SAVE and SECSAVE ...
# According to Alex, these are deprecated from most versions of tiegcm
#  -e 's/'"`grep 'SAVE'         ../tiegcm.nml | head -n 1`"'/'"$save"'/' \
#  -e 's/'"`grep 'SECSAVE'      ../tiegcm.nml`"'/'"$secsave"'/' \

   # There is a danger that you match multiple things ... F107 and F107A,
   # for example ... so try to grep whitespace too ...
   # FIXME TOMOKO ... do we want F107 and F107A to be identical

   sed \
   -e 's/'"`grep ' START_YEAR '   ../tiegcm.nml`"'/'"$start_year"'/' \
   -e 's/'"`grep ' START_DAY '    ../tiegcm.nml`"'/'"$start_day"'/' \
   -e 's/'"`grep ' SOURCE_START ' ../tiegcm.nml`"'/'"$source_start"'/' \
   -e 's/'"`grep ' START '        ../tiegcm.nml | head -n 4 | tail -n 1`"'/'"$start"'/' \
   -e 's/'"`grep ' STOP '         ../tiegcm.nml | head -n 1`"'/'"$stop"'/' \
   -e 's/'"`grep ' HIST '         ../tiegcm.nml | head -n 1`"'/'"$hist"'/' \
   -e 's/'"`grep ' SECSTART '     ../tiegcm.nml`"'/'"$secstart"'/' \
   -e 's/'"`grep ' SECSTOP '      ../tiegcm.nml`"'/'"$secstop"'/' \
   -e 's/'"`grep ' SECHIST '      ../tiegcm.nml`"'/'"$sechist"'/' \
   -e 's/'"`grep ' F107 '         ../tiegcm.nml | head -n 1`"'/'"$f107"'/' \
   ../tiegcm.nml >! tiegcm.nml

   #----------------------------------------------------------------------
   # Block 3: Run the model
   #----------------------------------------------------------------------

#  echo "ensemble member $ensemble_member : before tiegcm"
#  ncdump -v mtime tiegcm_restart.nc

# mpi
#  setenv TARGET_CPU_LIST "-1" 
#  mpirun.lsf ../tiegcm < tiegcm.nml |& tee tiegcm_out_$ensemble_member

# without mpi
   ../tiegcm < tiegcm.nml >& tiegcm_out_$ensemble_member || exit 3

#  echo "ensemble member $ensemble_member : after tiegcm"
#  ncdump -v mtime tiegcm_restart.nc

   ls -lrt

   #----------------------------------------------------------------------
   # Block 4: Convert the model output to form needed by DART
   # AT this point, the model has updated the information in tiegcm_restart_p.nc
   # We need to get that information back into the DART state vector.
   #
   # model_to_dart expects the tiegcm input file     to be 'tiegcm_restart_p.nc'
   # model_to_dart expects the tiegcm secondary file to be 'tiegcm_s.nc'
   # model_to_dart writes out a file with the default name 'dart_ics'
   #
   # The updated information needs to be moved into CENTRALDIR in
   # preparation for the next cycle.
   #----------------------------------------------------------------------

   ../model_to_dart || exit 4

   mv dart_ics            ../$output_file || exit 4
   mv tiegcm_s.nc         ../$tiesecond   || exit 4
   mv tiegcm_restart_p.nc ../$tierestart  || exit 4
   mv tiegcm.nml          ../$tieinp      || exit 4
   mv tiegcm_out_*        ../.            || exit 4

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line      = $input_file_line + 3
   @ output_file_line     = $output_file_line + 3
end

# Change back to original directory 
cd ..

# After you are assured this script works as expected, you can actually 
# remove the temporary directory. For now ... leave this commented OUT.
#\rm -rf $temp_dir

# Remove the filter_control file to signal completion
\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

