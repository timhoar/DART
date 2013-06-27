#!/bin/csh -f
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# The FORCE options are not optional.
# the VERBOSE options are useful for debugging.
set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fv --preserve=timestamps'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'
set LAUNCHCMD = mpirun.lsf

set ensemble_size = ${NINST_LND}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
echo "temp_dir is $temp_dir"
mkdir -p $temp_dir
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.clm2_${ensemble_member}.r.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `head -1 ../rpointer.lnd_0001`
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
# Get observation sequence file ... or die right away. Cannot specify -f on
# the link command and still check status.
#-----------------------------------------------------------------------------

set DART_OBS_DIR = ${MODEL_YEAR}${MODEL_MONTH}
set OBSDIR = ${WORK}/DART/observations/snow/work/obs_seqs
set OBSFNAME = obs_seq.0Z.${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY}
set OBS_FILE = ${OBSDIR}/${OBSFNAME}

\rm -f              obs_seq.out
\ln -vs ${OBS_FILE} obs_seq.out

set lnstat = $status
if ($lnstat != 0) then
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... ln died with status $lnstat"
   exit 1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
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
# The sampling error correction is a lookup table. The tables are stored
# in the DART distribution. Each ensemble size has its own (static) file 
# which does not need to be archived.  It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

set DARTDIR = ${WORK}/DART/models/clm/work

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
# Block 4: Convert CLM restart files to DART initial condition files.
# clm_to_dart is serial code, we can do all of these at the same time
# as long as we can have unique namelists for all of them.
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &clm_to_dart_nml:      clm_to_dart_output_file = 'dart_ics',
#
#=========================================================================

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml'

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set LND_RESTART_FILENAME = `printf ../../${MYCASE}.clm2_%04d.r.${MODEL_DATE_EXT}.nc  ${member}`
   set LND_HISTORY_FILENAME = `printf ../../${MYCASE}.clm2_%04d.h0.${MODEL_DATE_EXT}.nc ${member}`
   set     DART_IC_FILENAME = `printf ../filter_ics.%04d ${member}`

   ${LINK} $LND_RESTART_FILENAME clm_restart.nc
   ${LINK} $LND_HISTORY_FILENAME clm_history.nc
   ${LINK}     $DART_IC_FILENAME dart_ics

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

   cp ../input.nml .

   echo "starting clm_to_dart for member ${member} at "`date`
   ${EXEROOT}/clm_to_dart >! output.${member}.clm_to_dart &
   echo "finished clm_to_dart for member ${member} at "`date`

   cd ..

   @ member++
end

wait

if ($status != 0) then
   echo "ERROR ... DART died in 'clm_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'clm_to_dart' ... ERROR"
   exit 7
endif

#=========================================================================
# Block 5: Actually run the assimilation.
# Will result in a set of files : 'filter_restart.xxxx'
#
# DART namelist settings required:
# &filter_nml:           async                  = 0,
# &filter_nml:           adv_ens_command        = "./no_model_advance.csh",
# &filter_nml:           restart_in_file_name   = 'filter_ics'
# &filter_nml:           restart_out_file_name  = 'filter_restart'
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

# clm always needs a clm_restart.nc, clm_history.nc for geometry information, etc.

set LND_RESTART_FILENAME = ../${MYCASE}.clm2_0001.r.${MODEL_DATE_EXT}.nc
set LND_HISTORY_FILENAME = ../${MYCASE}.clm2_0001.h0.${MODEL_DATE_EXT}.nc

${LINK} $LND_RESTART_FILENAME clm_restart.nc
${LINK} $LND_HISTORY_FILENAME clm_history.nc

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
# Block 6: Update the clm restart files.
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_clm_nml:      dart_to_clm_input_file = 'dart_restart',
# &dart_to_clm_nml:      advance_time_present   = .false.
#=========================================================================

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory, which already exists
   # and has the required input files remaining from 'Block 4'

   cd member_${member}

   set DART_RESTART_FILE = `printf ../filter_restart.%04d ${member}`
   ${LINK} $DART_RESTART_FILE dart_restart

   echo "starting dart_to_clm for member ${member} at "`date`
   ${EXEROOT}/dart_to_clm >! output.${member}.dart_to_clm &
   echo "finished dart_to_clm for member ${member} at "`date`

   cd ..

   @ member++
end

wait

if ($status != 0) then
   echo "ERROR ... DART died in 'dart_to_clm' ... ERROR"
   echo "ERROR ... DART died in 'dart_to_clm' ... ERROR"
   exit 8
endif

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

#\rm -f ../$CASE.*.rh0.* 
#\rm -f ../$CASE.*.rs1.* 
#\rm -f ../PET*.ESMF_LogFile 

ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

