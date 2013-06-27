#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# ---------------------
# Purpose
# ---------------------
#
# The real purpose of this set of notes is to record what is needed to configure
# and build a CESM instance that has CAM, CLM, and CICE as active components
# in a multi-instance configuration over a single data ocean ... etc.
# Despite looking like a script, it might best be used as a set of notes.
# ---------------------
# How to set the script
# ---------------------
# -- Copy this script into your directory
# -- Choose a case name (by changing "setenv case" ) and save the script as $case.csh
# -- Set the case options at the top of the script
# -- If you have source mods, the script assumes that your mods are in: mods_$case.
#    So, if you have source mods, create a subdirectory mods_$case that contains your mods.
#    If you don t have any source mods, the script creates an empty subdirectory mods_$case.
# -- If you have namelist mods, you need to add them to the namelist template: user_nl_cam
#    Set your namelist variables over there (without modifying the syntax to create user_nl_cam
# -- Now, you are ready to go. Save your script and submit your run with the command: ./$case.csh
#    The script creates your case, configure, compile and submit your job.
# -- The script also creates a subdirectory (nml_$case) that contains your namelists.
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime settings,
# you need to delete everything and start the run from scratch.
#
# ./${CASENAME}.*.clean_build
# ./configure -cleanall
#
# ====================================================================
# === IMPORTANT modifications to the distribution code before ANYTHING 
# ====================================================================

# had to edit the following to remove the LSB_PJL_... word too long error
# cesm1_1_beta08/scripts/ccsm_utils/Machines/mkbatch.bluefire
# as long as OMP_NUM_THREADS == 1 ... the default is fine.

# Ensure that RUN_REFTOD has been added to 
# cesm1_1_beta08/scripts/ccsm_utils/Case.template/config_definition.xml

# Change cesm1_1_beta08/models/atm/cam/bld/cam.cpl7.template
# ncdata  = '${RUN_REFCASE}.cam.i.${RUN_REFDATE}-00000.nc'
# to
# ncdata  = '\${RUN_REFCASE}.cam_\${atm_inst_string}.i.\${RUN_REFDATE}-\${RUN_REFTOD}.nc'

# Change cesm1_1_beta08/models/ice/cice/bld/cice.cpl7.template
# set ice_ic = ${RUN_REFCASE}.cice.r.${RUN_REFDATE}-00000.nc
# to
# set ice_ic = \${RUN_REFCASE}.cice_\${ice_inst_string}.r.\${RUN_REFDATE}-\${RUN_REFTOD}.nc

# Change cesm1_1_beta08/models/lnd/clm/bld/clm.cpl7.template
# There's more to it than this for the beta08 distribution.
# xxdiff the DART development branch models/cam/shell_scripts/clm.cpl7.template and ...

# ====================================================================
# ====  Set case options
# ====================================================================

setenv case                 hybrid1
setenv compset              F_2000
setenv ccsmtag              cesm1_1_beta08
setenv resolution           f09_f09
setenv num_instances        4
setenv reuse_existing_case  true

# --------------------------------------------------------------------
# define machines and directories
# --------------------------------------------------------------------
#
# mach            computer name
# cesm_datadir    location of public CESM data files
# cesm_public     location of public CESM code distributions
# caseroot        your (future) cesm case directory
# rundir          (future) run-time directory
# archdir         (future) short-term archive directory
# ccsmroot        location of the cesm code base
# DARTdir         location of DART executables, scripts and input

setenv mach         bluefire
setenv cesm_datadir /glade/proj3/cseg/inputdata
setenv cesm_public  /glade/proj3/cseg
setenv caseroot     /glade/user/${USER}/cases/${case}
setenv rundir       /glade/scratch/${USER}/${case}
setenv archdir      /glade/scratch/${USER}/archive/${case}

setenv ccsmroot     ${cesm_public}/collections/${ccsmtag}
setenv ccsmroot     /glade/home/thoar/${ccsmtag}
setenv DARTdir      /glade/home/thoar/svn/DART/dev

# ======================
# configure settings
# ======================

setenv run_refdate 2008-10-31
setenv run_reftod  00000

setenv sst_dataset ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010

# ======================
# runtime settings
# ======================

setenv resubmit      0
setenv stop_n        6
setenv stop_option   nhours

# ======================
# job settings
# ======================

setenv proj         93300315
setenv timewall     1:30
setenv queue        lrg_regular

# ====================================================================
# Create the case.
# For list of the cases: ./create_newcase -list
# ====================================================================

# if reuse_existing_case is false and the directory does not exist, ...

if ("${reuse_existing_case}" == "false") then
   echo "removing old files from ${caseroot} and ${rundir}"
   \rm -fr ${caseroot}
   \rm -fr ${rundir}
   ${ccsmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset} -skip_rundb

   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit 1
   endif
else
   cd ${caseroot}
   ./configure  -cleannamelist
endif

# ====================================================================
# Configure the case.
# ====================================================================

cd ${caseroot}

./xmlchange -file env_build.xml    -id EXEROOT        -val ${rundir}
./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
#./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val ${nancy_scratch}/esmf-mpi

set num_tasks_per_instance = 16
set nthreads = 1
@ total_nt = $num_instances * $num_tasks_per_instance

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $total_nt
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_ATM -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $total_nt
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_LND -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $total_nt
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_ICE -val $num_instances

./xmlchange -file env_conf.xml -id RUN_TYPE                -val hybrid
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_refdate
./xmlchange -file env_conf.xml -id RUN_REFDATE             -val $run_refdate
./xmlchange -file env_conf.xml -id RUN_REFTOD              -val $run_reftod
./xmlchange -file env_conf.xml -id RUN_REFCASE             -val ${case}
./xmlchange -file env_conf.xml -id GET_REFCASE             -val FALSE
./xmlchange -file env_conf.xml -id BRNCH_RETAIN_CASENAME   -val TRUE
./xmlchange -file env_conf.xml -id SSTICE_DATA_FILENAME    -val $sst_dataset
./xmlchange -file env_conf.xml -id SSTICE_YEAR_ALIGN       -val $year_start
./xmlchange -file env_conf.xml -id SSTICE_YEAR_START       -val $year_start
./xmlchange -file env_conf.xml -id SSTICE_YEAR_END         -val $year_end
./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-rtm off'

# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics.

./xmlchange -file env_run.xml      -id CONTINUE_RUN     -val FALSE
./xmlchange -file env_run.xml      -id RESUBMIT         -val $resubmit
./xmlchange -file env_run.xml      -id STOP_OPTION      -val $stop_option
./xmlchange -file env_run.xml      -id STOP_N           -val $stop_n
./xmlchange -file env_run.xml      -id CALENDAR         -val GREGORIAN

# Substantial archiving changes exist in the Tools/st_archive.sh script.
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S                     -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val FALSE 
./xmlchange -file env_run.xml -id DOUT_L_HTAR                -val FALSE

# ====================================================================
# Create namelist template: user_nl_cam, user_nl_clm
# ====================================================================

cd ${caseroot}

cat <<EOF >! user_nl_cam
&camexp
 inithist             = 'ENDOFRUN'
 div24del2flag        = 4
 empty_htapes         = .true.
 fincl1               = 'PHIS:I'
 nhtfrq               = -$stop_n
 iradae               = -$stop_n
 aerodep_flx_datapath = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero'
 aerodep_flx_file     = 'aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'
 aerodep_flx_cycle_yr = 2000
 aerodep_flx_type     = 'CYCLICAL'
/
EOF

# at least one of these caused problems when running with the full cesm:
# empty_htapes                 = .true.
# nhtfrq                       = -12
# faerdep ... perhaps only used if CN modeling on ... 

cat <<EOF >! user_nl_clm
&clmexp
  fatmgrid = '${cesm_datadir}/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc'
/
EOF

# at least one of these caused problems when running with the full cesm:
#  hist_nhtfrq = -12
#  hist_empty_htapes = .true.

# ====================================================================
# Update source files if need be
# ====================================================================

\cp -rf ~thoar/${ccsmtag}/SourceMods/* ${caseroot}/SourceMods/
if ( $status == 0) then
   echo "FYI - Local Source Modifications used for this case:"
   ls -lr ${caseroot}/SourceMods/*
else
   echo "FYI - No SourceMods for this case"
endif

# ====================================================================
# Configure
# ====================================================================

cd ${caseroot}

./configure -case

if ( $status != 0 ) then
   echo "ERROR: Case could not be configured."
   exit 2
endif

# ====================================================================
# Stage a copy of the DART assimilate.csh script HERE
# ====================================================================

cd ${caseroot}

\mv Tools/st_archive.sh Tools/st_archive.sh.org
\cp -f ${DARTdir}/models/cam/shell_scripts/st_archive_hybrid.sh Tools/st_archive.sh

# TJH Until I get it working, use a link ... 
ln -sf ${DARTdir}/models/cam/shell_scripts/assimilate.hybrid.csh assimilate.csh

# ====================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText"
# ====================================================================

cd ${caseroot}

echo ''
echo 'Adding the call to assimilate.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform an assimilation with DART.
# -------------------------------------------------------------------------

set CplLogFile = `ls -1t cpl.log* | head -1`
if ($CplLogFile == "") then
 echo 'ERROR: Model did not complete - no cpl.log file present - exiting'
 echo 'ERROR: Assimilation will not be attempted.'
 exit -4
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
  ${CASEROOT}/assimilate.csh

  if ( $status == 0 ) then
     echo "`date` -- DART HAS FINISHED"
  else
     echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
     exit -5
  endif
endif

# END OF DART BLOCK
# -------------------------------------------------------------------------

"EndOfText"

# Now that the "here" document is created, 
# determine WHERE to insert it -- ONLY IF it is not already there.

grep --line-number "ABANDON HOPE" ${case}.run
if ( $status > 1 ) then

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ orglen = `cat ${case}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $orglen - $keep 

   mv ${case}.run ${case}.run.orig

   head -$keep      ${case}.run.orig >! ${case}.run
   cat              add_to_run.txt   >> ${case}.run
   tail -$lastlines ${case}.run.orig >> ${case}.run

endif

# ====================================================================
# We are trying to make all resubmits to be continuation ... 
# ====================================================================

cd ${caseroot}/Tools

echo ''
echo 'Require all resubmits to be hybrid starts, which means'
echo 'CONTINUE_RUN should be FALSE in Tools/ccsm_postrun.csh'
echo ''

ex ccsm_postrun.csh <<ex_end
/use COMP_RUN_BARRIERS as surrogate for timing run logical/
/CONTINUE_RUN/
s;TRUE;FALSE;
wq
ex_end

# ====================================================================
# build
# ====================================================================

cd ${caseroot}

echo ''
echo 'Building the case'
echo ''

./$case.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit 3
endif

# ====================================================================
# Stage the restarts now that the run directory exists
# ====================================================================

# 20081031 ... /ptmp/thoar/restarts
# 20080801 ... /glade/proj3/DART/raeder/FV1deg_4.0/Exp1/obs_0000
set stagedir = /glade/proj3/DART/raeder/FV1deg_4.0/Exp1/obs_0000
set stagedir = /ptmp/thoar/restarts

echo ''
echo "Staging the restarts from {$stagedir}"
echo ''

@ n = 1
while ($n <= $num_instances)

   echo "Staging restarts for instance $n of $num_instances"

   set ATMFILE = ${stagedir}/CAM/caminput_${n}.nc
   set LNDFILE = ${stagedir}/CLM/clminput_${n}.nc
   set ICEFILE = ${stagedir}/ICE/iceinput_${n}.nc

   # Must decode the valid time of each file and reconstitute
   # a name of the appropriate form. 

   set DATESTR = `ncdump -v date    ${ATMFILE} | grep "date ="` 
   set  TODSTR = `ncdump -v datesec ${ATMFILE} | grep "datesec ="` 
   set DATESTR = `echo $DATESTR | sed -e "s#=# #g"`
   set  TODSTR = `echo  $TODSTR | sed -e "s#=# #g"`
   set  MODEL_YEAR = `echo "$DATESTR[2] / 10000" | bc`
   set   REMAINDER = `echo "$DATESTR[2] - $MODEL_YEAR*10000" | bc`
   set MODEL_MONTH = `echo "$REMAINDER / 100" | bc`
   set   MODEL_DAY = `echo "$REMAINDER - $MODEL_MONTH*100" | bc`
   set   MODEL_TOD = `echo $TODSTR[2] | bc`

   # sanity check that these are for the same time as atm ...

   set CLMDATESTR = `ncdump -v timemgr_rst_curr_ymd ${LNDFILE} | grep "timemgr_rst_curr_ymd ="` 
   set  CLMTODSTR = `ncdump -v timemgr_rst_curr_tod ${LNDFILE} | grep "timemgr_rst_curr_tod ="` 
   set CLMDATESTR = `echo $CLMDATESTR | sed -e "s#=# #g"`
   set  CLMTODSTR = `echo  $CLMTODSTR | sed -e "s#=# #g"`

   if ( $CLMDATESTR[2] != $DATESTR[2] ) then
      echo "Times of CAM and CLM input files do not match." 
      echo "$ATMFILE has $DATESTR[2]"
      echo "$LNDFILE has $CLMDATESTR[2]"
      echo "$ICEFILE has no time information"
   endif

   # create new filenames

   set NEWATMFILE = `printf "${case}.cam_%04d.i.%04d-%02d-%02d-%05d.nc"  $n $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_TOD`
   set NEWLNDFILE = `printf "${case}.clm2_%04d.r.%04d-%02d-%02d-%05d.nc" $n $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_TOD`
   set NEWICEFILE = `printf "${case}.cice_%04d.r.%04d-%02d-%02d-%05d.nc" $n $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_TOD`

   cp --preserve=timestamps ${ATMFILE} ${rundir}/run/${NEWATMFILE}
   cp --preserve=timestamps ${LNDFILE} ${rundir}/run/${NEWLNDFILE}
   cp --preserve=timestamps ${ICEFILE} ${rundir}/run/${NEWICEFILE}

 @ n++
end

echo 'If inflation is being used ... '
echo "must stage a ${rundir}/[prior,pos]_inflate_restart.YYYY-MM-DD-SSSSS"

# ====================================================================
# Edit the run script to reflect project, queue, and wallclock
# ====================================================================

set PROJ=`grep BSUB $case.run | grep -e '-P' `
sed s/$PROJ[3]/$proj/ < $case.run >! temp
/bin/mv temp  $case.run

set TIMEWALL=`grep BSUB $case.run | grep -e '-W' `
sed s/$TIMEWALL[3]/$timewall/ < $case.run >! temp
/bin/mv temp  $case.run

set QUEUE=`grep BSUB $case.run | grep -e '-q' `
sed s/$QUEUE[3]/$queue/ < $case.run >! temp
/bin/mv temp  $case.run

chmod 0744 $case.run

# ====================================================================
# Submit job
# ====================================================================

set MYSTRING = `grep "set DARTDIR" assimilate.csh`
set DARTDIR = $MYSTRING[4]

echo ''
echo 'case is ready to submit after you check the'
echo "DART settings in ${DARTDIR}/input.nml"
echo 'After you check them,'
echo "cd into ${caseroot} and run: ./$case.submit"
echo ''

