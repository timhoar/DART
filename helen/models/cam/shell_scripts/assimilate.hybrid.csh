#!/bin/csh -f
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# This block is an attempt to localize all the machine-specific 
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "starting assimilate script at "`date`

switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      # The FORCE options are not optional.
      # the VERBOSE options are useful for debugging.
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

      set BASEOBSDIR = /glade/proj3/image/Observations/ACARS
      set DARTDIR    = ${HOME}/svn/DART/dev
      set LAUNCHCMD  = mpirun.lsf

   breaksw
   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
      set DARTDIR    = ${HOME}/devel
      set LAUNCHCMD  = "aprun -n $NTASKS"
   breaksw
endsw 

set ensemble_size = ${NINST_ATM}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
echo "temp_dir is $temp_dir"
mkdir -p $temp_dir
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `ls -1t ../*.cam_0001.i.* | head -n 1`
set FILE = $FILE:t
set FILE = $FILE:r
set MYCASE = `echo $FILE | sed -e "s#\..*##"`
set MODEL_DATE_EXT = `echo $FILE:e`
set MODEL_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set MODEL_YEAR     = $MODEL_DATE[1]
set MODEL_MONTH    = $MODEL_DATE[2]
set MODEL_DAY      = $MODEL_DATE[3]
set MODEL_SECONDS  = $MODEL_DATE[4]
set MODEL_HOUR     = `echo $MODEL_DATE[4] / 3600 | bc`

echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_SECONDS (seconds)"
echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_HOUR (hours)"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DART_OBS_DIR = ${MODEL_YEAR}${MODEL_MONTH}_6H
set OBSDIR       = ${BASEOBSDIR}/${DART_OBS_DIR}

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/${FILE} not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/${FILE} not found ... ERROR"
   exit 1
endif

# Modify the DART input.nml such that
# the DART ensemble size matches the CESM number of instances
# WARNING: the output files contain ALL enemble members ==> BIG

ex input.nml <<ex_end
g;ens_size ;s;= .*;= $ensemble_size;
g;num_output_state_members ;s;= .*;= $ensemble_size;
g;num_output_obs_members ;s;= .*;= $ensemble_size;
wq
ex_end

#=========================================================================
# Block 2: Stage the files needed for SAMPLING ERROR CORRECTION
#
# Each ensemble size has its own (static) file which does not need to be archived.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

set  MYSTRING = `grep sampling_error_correction input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set SECSTRING = `echo $MYSTRING[2] | tr 'A-Z' 'a-z'`

if ( $SECSTRING == true ) then
   set SAMP_ERR_FILE = ${DARTDIR}/system_simulation/final_full_precomputed_tables/final_full.${ensemble_size}
   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo "ERROR: no sampling error correction file for this ensemble size."
      echo "ERROR: looking for ${SAMP_ERR_FILE}"
      exit 2
   endif
else
   echo "Sampling Error Correction not requested for this assimilation."
endif

#=========================================================================
# Block 3: DART INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
# inf_out_file_name           = 'prior_inflate_restart', 'post_inflate_restart',
# inf_diag_file_name          = 'prior_inflate_diag',    'post_inflate_diag',
#
# NOTICE: the archiving scripts more or less require the names of these
# files to be as listed above. When being archived, the filenames get a
# unique extension (describing the assimilation time) appended to them.
#
# The inflation file is essentially a duplicate of the model state ...
# it is slaved to a specific geometry. The initial files are created
# offline with values of unity. For the purpose of this script, they are
# thought to be the output of a previous assimilation, so they should be
# named something like prior_inflate_restart.YYYY-MM-DD-SSSSS
#
# The first inflation file can be created with 'fill_inflation_restart'
# which can be built in the usual DART manner.
#
# The strategy is to use the LATEST inflation file from the CESM 'rundir'.
# After an assimilation, the new inflation values/files will be moved to
# the CESM rundir to be used for subsequent assimilations. If the short-term
# archiver has worked correctly, only the LATEST files will available. Of
# course, it is not required to have short-term archiving turned on, so ...
#=========================================================================

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep inf_initial_from_restart input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr [A-Z] [a-z]`
set  POSTE_TF = `echo $MYSTRING[3] | tr [A-Z] [a-z]`

# its a little tricky to remove both styles of quotes from the string.

set  MYSTRING = `grep inf_in_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_IFNAME = $MYSTRING[2]
set  POSTE_INF_IFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_out_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_OFNAME = $MYSTRING[2]
set  POSTE_INF_OFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_diag_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_DIAG = $MYSTRING[2]
set  POSTE_INF_DIAG = $MYSTRING[3]

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == ".false.") then
      echo "ERROR: inf_flavor(1) = $PRIOR_INF, yet inf_initial_from_restart = $PRIOR_TF"
      echo "ERROR: fix input.nml to reflect whether you want prior inflation or not."
      exit 3
   endif

   # Look for the output from the previous assimilation
   (ls -rt1 ../${PRIOR_INF_OFNAME}.* | tail -1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   # If one exists, use it as input for this assimilation
   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest ${PRIOR_INF_IFNAME}
   else
      echo "ERROR: Requested prior inflation but specified no incoming prior inflation file."
      echo "ERROR: expected something like ../${PRIOR_INF_OFNAME}.YYYY-MM-DD-SSSSS"
      exit 4
   endif
else
   echo "Prior Inflation not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ( $POSTE_INF > 0 ) then

   if ($POSTE_TF == ".false.") then
      echo "ERROR: inf_flavor(2) = $POSTE_INF, yet inf_initial_from_restart = $POSTE_TF"
      echo "ERROR: fix input.nml to reflect whether you want posterior inflation or not."
      exit 5
   endif

   # Look for the output from the previous assimilation
   (ls -rt1 ../${POSTE_INF_OFNAME}.* | tail -1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   # If one exists, use it as input for this assimilation
   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest ${POSTE_INF_IFNAME}
   else
      echo "ERROR: Requested POSTERIOR inflation but specified no incoming POSTERIOR inflation file."
      echo "ERROR: expected something like ../${POSTE_INF_OFNAME}.YYYY-MM-DD-SSSSS"
      exit 6
   endif
else
   echo "Posterior Inflation not requested for this assimilation."
endif

#=========================================================================
#
# Block 4: Convert CAM restart files to DART initial condition files.
# cam_to_dart is serial code, we can do all of these at the same time
# as long as we can have unique namelists for all of them.
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ic_old'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &cam_to_dart_nml:      cam_to_dart_output_file = 'dart_ics',
#
#=========================================================================

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   # Turns out the .h0. files are timestamped with the START of the 
   # run, which is *not* MODEL_DATE_EXT ...  I just link to a whatever 
   # is convenient (since the info is static).

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set ATM_INITIAL_FILENAME = `printf ../../${MYCASE}.cam_%04d.i.${MODEL_DATE_EXT}.nc ${member}`
   set ATM_HISTORY_FILENAME = `ls -1t ../../${MYCASE}.cam*.h0.* | head -n 1`
   set DART_IC_FILE = `printf ../filter_ic_old.%04d ${member}`

   ${LINK} $ATM_INITIAL_FILENAME caminput.nc
   ${LINK} $ATM_HISTORY_FILENAME cam_phis.nc
   ${LINK} $DART_IC_FILE         dart_ics

   cp ../input.nml .

   echo "starting cam_to_dart for member ${member} at "`date`
   ${EXEROOT}/cam_to_dart >! output.${member}.cam_to_dart &
   echo "finished cam_to_dart for member ${member} at "`date`

   cd ..

   @ member++
end

wait

if ($status != 0) then
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'cam_to_dart' ... ERROR"
   exit 7
endif

#=========================================================================
# Block 5: Actually run the assimilation.
# Will result in a set of files : 'filter_restart.xxxx'
#
# DART namelist settings required:
# &filter_nml:           async                  = 0,
# &filter_nml:           adv_ens_command        = "./no_model_advance.csh",
# &filter_nml:           restart_in_file_name   = 'filter_ic_old'
# &filter_nml:           restart_out_file_name  = 'filter_ic_new'
# &filter_nml:           obs_sequence_in_name   = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name  = 'obs_seq.final'
# &filter_nml:           init_time_days         = -1,
# &filter_nml:           init_time_seconds      = -1,
# &filter_nml:           first_obs_days         = -1,
# &filter_nml:           first_obs_seconds      = -1,
# &filter_nml:           last_obs_days          = -1,
# &filter_nml:           last_obs_seconds       = -1,
# &ensemble_manager_nml: single_restart_file_in = .false.
#
#=========================================================================

# CAM:static_init_model() always needs a caminput.nc and a cam_phis.nc
# for geometry information, etc.

set ATM_INITIAL_FILENAME = ../${MYCASE}.cam_0001.i.${MODEL_DATE_EXT}.nc
set ATM_HISTORY_FILENAME = `ls -1t ../${MYCASE}.cam_0001.h0.* | head -n 1`

${LINK} $ATM_INITIAL_FILENAME caminput.nc
${LINK} $ATM_HISTORY_FILENAME cam_phis.nc

# Determine proper observation sequence file.

set OBSFNAME = `printf obs_seq${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY}%02d ${MODEL_HOUR}`
set OBS_FILE = ${OBSDIR}/${OBSFNAME}

${LINK} ${OBS_FILE} obs_seq.out

echo "assimilate:starting filter at "`date`
${LAUNCHCMD} ${EXEROOT}/filter || exit 7
echo "assimilate:finished filter at "`date`

${MOVE} Prior_Diag.nc      ../Prior_Diag.${MODEL_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../Posterior_Diag.${MODEL_DATE_EXT}.nc
${MOVE} obs_seq.final      ../obs_seq.${MODEL_DATE_EXT}.final
${MOVE} dart_log.out       ../dart_log.${MODEL_DATE_EXT}.out

# Accomodate any possible inflation files
# 1) rename file to reflect current date
# 2) move to CENTRALDIR so the DART INFLATION BLOCK works next time and
#    that they can get archived.

foreach FILE ( ${PRIOR_INF_OFNAME} ${POSTE_INF_OFNAME} ${PRIOR_INF_DIAG} ${POSTE_INF_DIAG} )
   if ( -e ${FILE} ) then
      ${MOVE} ${FILE} ../${FILE}.${MODEL_DATE_EXT}
   else
      echo "No ${FILE} for ${MODEL_DATE_EXT}"
   endif
end

#=========================================================================
# Block 6: Update the cam restart files. 
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_ic_new'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_cam_nml:      dart_to_cam_input_file = 'temp_ic',
# &dart_to_cam_nml:      advance_time_present   = .false.
# &atm_in_xxxx:ncdata = 'cam_initial_x.nc'
#=========================================================================

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory, which already exists
   # and has the required input files remaining from 'Block 4'

   cd member_${member}

   set DART_RESTART_FILE = `printf ../filter_ic_new.%04d ${member}`
   ${LINK} $DART_RESTART_FILE temp_ic

   echo "starting dart_to_cam for member ${member} at "`date`
   ${EXEROOT}/dart_to_cam >! output.${member}.dart_to_cam &
   echo "finished dart_to_cam for member ${member} at "`date`

   cd ..

   @ member++
end

wait

if ($status != 0) then
   echo "ERROR ... DART died in 'dart_to_cam' ... ERROR"
   echo "ERROR ... DART died in 'dart_to_cam' ... ERROR"
   exit 8
endif

#-------------------------------------------------------------------------
# Now that everything is staged, we have to communicate the current
# model time to the drv_in&seq_timemgr_inparm namelist
# which is built from CASEROOT/user_nl_drv by the *.run script
#-------------------------------------------------------------------------
set mydir = `pwd`
cd ${CASEROOT}
./xmlchange -file env_run.xml  -id START_TOD     -val ${MODEL_SECONDS}
./xmlchange -file env_conf.xml -id RUN_STARTDATE -val ${MODEL_YEAR}-${MODEL_MONTH}-${MODEL_DAY}
./xmlchange -file env_conf.xml -id RUN_REFDATE   -val ${MODEL_YEAR}-${MODEL_MONTH}-${MODEL_DAY}
./xmlchange -file env_conf.xml -id RUN_REFTOD    -val ${MODEL_SECONDS}
${COPY} env_conf.xml LockedFiles/env_conf.xml.locked
cd $mydir

# we (DART) do not need these files, and CESM does not need them either
# to continue a run.  if we remove them here they do not get moved to
# the short-term archiver.
${REMOVE} ../*.rs.*
${REMOVE} ../*.rh0.*
${REMOVE} ../*.rs1.*
${REMOVE} ../PET*ESMF_LogFile

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------
echo "finished assimilate script at "`date`

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

