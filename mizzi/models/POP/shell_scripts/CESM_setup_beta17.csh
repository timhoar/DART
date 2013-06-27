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
# This script is designed to configure and build a multi-instance CESM model
# that has CAM, CLM, and CICE as active components over a single data ocean,
# and will use DART to assimilate observations at regular intervals.
# This script does not build DART.

# ---------------------
# How to set up the script
# ---------------------
# -- Either edit and run this script in the $DART/models/clm/shell_scripts
#    directory where it first appears,
#    or copy it to somewhere that it will be preserved and run it there.
#    It will create a 'case' directory, where the model will be built,
#    and an execution directory, where each forecast and assimilation will
#    take place.  The short term archiver will use a third directory for
#    storage of model output until it can be moved to long term storage (HPSS)
# -- Examin the whole script to identify things to change for your experiments.
# -- Provide any initial files needed by your run:
#       inflation
#       sampling error correction
#       CAM/CLM/CICE initial ensemble
#       ...
# -- Run this script.
# -- Edit the DART input.nml that appears in the $caseroot directory.
# -- Submit the job using $caseroot/${case}.${mach}.submit
#      ($mach may not be needed for cesm releases after cesm1_1_beta04)
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime
# settings, you need to delete everything and start the run from scratch.
#
# ./${CASENAME}.*.clean_build
# ./configure -cleanall
#
# ====================================================================
# ====  Set case options
# ====================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members
# reuse_existing_case:
#    false; Remove $caseroot and $exeroot and rebuild
#    true;  configure -cleannamelist

setenv case                 pop_beta17
setenv compset              CIAF
setenv cesmtag              cesm1_1_beta17
setenv resolution           T62_gx1v6
setenv num_instances        4
setenv reuse_existing_case  false

# ====================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1_beta04 on bluefire, MUST use 'thoar' value provided.
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/dev/...
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                 caseroot will be deleted if reuse_existing_case is false (below)
#                 So this script, and other useful things should be kept elsewhere.
# exeroot         (Future) Run-time directory; scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# ====================================================================

setenv mach         bluefire
setenv cesm_datadir /glade/proj3/cseg/inputdata
setenv cesmroot     /glade/proj3/cseg/collections/$cesmtag

setenv DARTroot     /glade/home/${USER}/svn/DART/dev

setenv caseroot     /glade/user/${USER}/cases/${case}
setenv exeroot      /ptmp/${USER}/${case}
setenv archdir      /ptmp/${USER}/archive/${case}

# ====================================================================
# configure settings
# ====================================================================

setenv refyear 2000
setenv stream_year_first 2000
setenv stream_year_last  2006
setenv stream_year_align 2000
setenv refmon  01
setenv refday  01
setenv run_refdate $refyear-$refmon-$refday
setenv run_reftod  00000

# ====================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
#               Changing stop_option requires changes to user_nl_clm below.
# stop_n        Number of time units in the forecast
# ====================================================================

setenv resubmit      4
setenv stop_option   ndays
setenv stop_n        1

# ====================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#
# TJH How many f19_f19 CLM instances can fit on 1 'regular' node?
# ====================================================================

setenv proj         93300315
setenv timewall     0:29
setenv queue        premium

# ====================================================================
# set these standard commands based on the machine you are running on.
# ====================================================================

switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      # The FORCE options are not optional.
      # the VERBOSE options are useful for debugging.
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

   breaksw
   default:
      # NERSC "hopper"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

   breaksw
endsw


# ====================================================================
# Create the case.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ====================================================================

if ("${reuse_existing_case}" == "false") then
   echo "removing old files from ${caseroot} and ${exeroot}"
   ${REMOVE} ${caseroot}
   ${REMOVE} ${exeroot}
   ${cesmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset}

   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit 1
   endif
else
   cd ${caseroot}
   ./configure -cleannamelist
   ./configure -cleanmach

endif

# ====================================================================
# Configure the case.
# ====================================================================

cd ${caseroot}

./xmlchange -file env_build.xml -id CALENDAR      -val GREGORIAN
./xmlchange -file env_build.xml -id EXEROOT       -val ${exeroot}

# with I compset ... CLM & DATM get num_instances, all else ninst == 1
#
# I really dont understand the logic of task layout.
# "everything on 320 pes except land on 960.
# cpl rootpe=0
# lnd rootpe=320
# ice rootpe=320
# ocn rootpe=640
# atm rootpe=960" -- Tony Craig

if ($num_instances == 4) then
   # Tiny case - can fit on one node 64 tasks.
   # 1/8 cpl,  1/8 cice, 1/8 ocn, 5/8 lnd.
   @ total_nt = 64
   @ atm_pes  = $total_nt
   @ cpl_pes  = $total_nt / 8
   @ cice_pes = $total_nt / 8
   @  ocn_pes = $total_nt / 8
   @ lnd_pes  = $total_nt - ($cpl_pes + $cice_pes + $ocn_pes)
else if ($num_instances == 40) then
   @ total_nt = 640
   # 10 nodes
   # 2/8 cpl,  1/8 cice, 1/8 ocn, 5/8 lnd.
   @ atm_pes  = $total_nt
   @ cpl_pes  = $total_nt / 4
   @ cice_pes = $total_nt / 8
   @  ocn_pes = $total_nt / 8
   @ lnd_pes  = $total_nt - ($cpl_pes + $cice_pes + $ocn_pes)
else
   echo "Unsupported number of instances ... dying a horrible death."
   exit 4
endif

# first pe of each component (counted from 0)
@  atm_rootpe = 0
@  cpl_rootpe = 0

@ cice_rootpe = $cpl_rootpe  + $cpl_pes
@  ocn_rootpe = $cice_rootpe + $cice_pes
@  lnd_rootpe = $ocn_rootpe + $ocn_pes

echo "check pe counting ..."
@ last_pe = $atm_rootpe + $atm_pes
echo "last pe => $last_pe =?= $total_nt <= total_nt"

echo "task layout"
echo "[$atm_rootpe .......................... ATM ..................... $total_nt]"
echo "[$cpl_rootpe ... CPL ... $cice_rootpe ... ICE ... $ocn_rootpe ... OCN ... $lnd_rootpe ... LND ... $total_nt]"
echo ""
echo "ATM gets $atm_pes"
echo "CPL gets $cpl_pes"
echo "ICE gets $cice_pes"
echo "OCN gets $ocn_pes"
echo "LND gets $lnd_pes"
echo ""

./xmlchange -file env_pesetup.xml -id NTASKS_ATM -val $atm_pes
./xmlchange -file env_pesetup.xml -id NTHRDS_ATM -val 1
./xmlchange -file env_pesetup.xml -id ROOTPE_ATM -val $atm_rootpe
./xmlchange -file env_pesetup.xml -id  NINST_ATM -val $num_instances

./xmlchange -file env_pesetup.xml -id NTASKS_LND -val $lnd_pes
./xmlchange -file env_pesetup.xml -id NTHRDS_LND -val 1
./xmlchange -file env_pesetup.xml -id ROOTPE_LND -val $lnd_rootpe
./xmlchange -file env_pesetup.xml -id  NINST_LND -val 1

./xmlchange -file env_pesetup.xml -id NTASKS_ICE -val $cice_pes
./xmlchange -file env_pesetup.xml -id NTHRDS_ICE -val 1
./xmlchange -file env_pesetup.xml -id ROOTPE_ICE -val $cice_rootpe
./xmlchange -file env_pesetup.xml -id  NINST_ICE -val 1

./xmlchange -file env_pesetup.xml -id NTASKS_OCN -val $ocn_pes
./xmlchange -file env_pesetup.xml -id NTHRDS_OCN -val 1
./xmlchange -file env_pesetup.xml -id ROOTPE_OCN -val $ocn_rootpe
./xmlchange -file env_pesetup.xml -id  NINST_OCN -val $num_instances

./xmlchange -file env_pesetup.xml -id NTASKS_CPL -val $cpl_pes
./xmlchange -file env_pesetup.xml -id NTHRDS_CPL -val 1
./xmlchange -file env_pesetup.xml -id ROOTPE_CPL -val $cpl_rootpe

./xmlchange -file env_pesetup.xml -id NTASKS_GLC -val $ocn_pes
./xmlchange -file env_pesetup.xml -id NTHRDS_GLC -val 1
./xmlchange -file env_pesetup.xml -id ROOTPE_GLC -val $ocn_rootpe
./xmlchange -file env_pesetup.xml -id  NINST_GLC -val 1

./xmlchange -file env_run.xml -id RUN_TYPE                   -val branch
./xmlchange -file env_run.xml -id RUN_STARTDATE              -val $run_refdate
./xmlchange -file env_run.xml -id RUN_REFCASE                -val c.da48
./xmlchange -file env_run.xml -id RUN_REFDATE                -val $run_refdate
./xmlchange -file env_run.xml -id RUN_REFTOD                 -val $run_reftod
./xmlchange -file env_run.xml -id BRNCH_RETAIN_CASENAME      -val FALSE
./xmlchange -file env_run.xml -id GET_REFCASE                -val FALSE

./xmlchange -file env_run.xml -id STOP_OPTION                -val $stop_option
./xmlchange -file env_run.xml -id STOP_N                     -val $stop_n
./xmlchange -file env_run.xml -id CONTINUE_RUN               -val FALSE
./xmlchange -file env_run.xml -id RESUBMIT                   -val $resubmit

./xmlchange -file env_run.xml -id DOUT_S                     -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val FALSE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MSROOT              -val "csm/${case}"
./xmlchange -file env_run.xml -id DOUT_L_HTAR                -val FALSE

./xmlchange -file env_run.xml -id DATM_MODE                  -val CPLHIST3HrWx
./xmlchange -file env_run.xml -id DATM_CPLHIST_CASE          -val $case
./xmlchange -file env_run.xml -id DATM_CPLHIST_YR_ALIGN      -val $refyear
./xmlchange -file env_run.xml -id DATM_CPLHIST_YR_START      -val $refyear
./xmlchange -file env_run.xml -id DATM_CPLHIST_YR_END        -val $refyear

./setup

if ( $status != 0 ) then
   echo "ERROR: Case could not be set up."
   exit 2
endif

# ====================================================================
# Modify namelist templates for each instance.
#
# hist_empty_htapes = .true.     suppresses the creation of all history files
# hist_fincl1 = 'TG',            except the first one, which will have one variable
# hist_nhtfrq = -$stop_n,        create one every $stop_n HOURS
# hist_mfilt  =  1,              with precisely one day in it
# hist_avgflag_pertape = 'I'     use instantaneous values - no average
#
# The fincl2 history tape has the half-hourly flux tower observations.
# The observation operators in obs_def_tower_mod.f90
# are going to read from the .h1. history file for these values.
# ====================================================================

@ inst = 1
while ($inst <= $num_instances)

   set instance = `printf %04d $inst`

   set fname = "user_nl_datm_$instance"
   echo "dtlimit = 1.5, 1.5, 1.5"                    >> $fname
   echo "fillalgo = 'nn', 'nn', 'nn'"                >> $fname
   echo "fillmask = 'nomask','nomask','nomask'"      >> $fname
   echo "mapalgo = 'bilinear','bilinear','bilinear'" >> $fname
   echo "mapmask = 'nomask','nomask','nomask'"       >> $fname
   echo "streams = 'datm.streams.txt.CPLHIST3HrWx.Solar_"$instance             $stream_year_align $stream_year_first $stream_year_last"'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.Precip_"$instance            $stream_year_align $stream_year_first $stream_year_last"'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.nonSolarNonPrecip_"$instance $stream_year_align $stream_year_first $stream_year_last"'" >> $fname
   echo "taxmode = 'cycle','cycle','cycle'"          >> $fname
   echo "tintalgo = 'coszen','coszen','coszen'"      >> $fname
   echo "dataMode = 'CPLHIST'"      >> $fname

   @ inst ++
end

set fname = "user_nl_dice"
echo "streams = 'dice.streams.txt.SSMI_IAF" 1948 1948 2004"'" >> $fname
echo "dataMode = 'SSTDATA'"                                   >> $fname

set fname = "user_nl_dlnd"
echo "streams = 'dlnd.streams.txt.rof.diatren_iaf_rx1" 1 1 1"'" >> $fname
echo "dataMode = 'DIATREN_IAF_RX1'"                             >> $fname

# ====================================================================
# to create custom streamfiles ...
# ====================================================================

./preview_namelists

foreach FILE (CaseDocs/*streams*)
   set FNAME = $FILE:t
   cp -v $FILE user_$FNAME
end

echo "Change each of the user_ streams files"
echo "Change each of the user_ streams files"
echo "Change each of the user_ streams files"
echo "Change each of the user_ streams files"

exit

./preview_namelists

# ====================================================================
# Update source files if need be
#    Ideally, using DART will not require any source mods.
#    Until then, this script accesses source mods from a hard-wired location below.
#    Those may eventually be packaged into the DART repository.
#    If you have additional source mods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ====================================================================

# this one needs a recursive copy to get all the files in the subdirs
${COPY} -r  ~thoar/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
if ( $status == 0) then
   echo "FYI - Local Source Modifications used for this case:"
   ls -lr ${caseroot}/SourceMods/*
else
   echo "FYI - No SourceMods for this case"
endif

# ====================================================================
# Stage the required parts of DART in the caseroot directory.
# ====================================================================

# if ("${reuse_existing_case}" == "false") then
#    ${MOVE} Tools/st_archive.sh Tools/st_archive.sh.orig
# endif
# ${COPY} ${DARTroot}/models/POP/shell_scripts/st_archive.sh Tools/
# ${COPY} ${DARTroot}/models/POP/shell_scripts/datm.buildnml.csh Buildconf/
# 
# ${COPY} ${DARTroot}/models/POP/shell_scripts/assimilate.csh  assimilate.csh
# ${COPY} ${DARTroot}/models/POP/work/input.nml                .

# ====================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText",
# not 'EndOfText' or EndOfText.
# ====================================================================

echo ''
echo 'Adding the call to assimilate.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform an assimilation with DART.
# -------------------------------------------------------------------------

set CplLogFile = `ls -1t cpl.log* | head -n 1`
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

grep "ABANDON HOPE" ${case}.run
set STATUSCHECK = $status

if ( ${STATUSCHECK} == 0 ) then
   echo "DART block already present in ${case}.run"
else if ( ${STATUSCHECK} == 1 ) then

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ origlen = `cat ${case}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $origlen - $keep

   mv ${case}.run ${case}.run.orig

   head -n $keep      ${case}.run.orig >! ${case}.run
   cat                add_to_run.txt   >> ${case}.run
   tail -n $lastlines ${case}.run.orig >> ${case}.run

endif

# ====================================================================
# Edit the run script to reflect project, queue, and wallclock
# ====================================================================

echo ''
echo 'Updating the run script to set project number, wallclock, and queue.'
echo ''

set PROJ=`grep BSUB $case.run | grep -e '-P' `
sed s/$PROJ[3]/$proj/ < $case.run >! temp
${MOVE} temp  $case.run

set TIMEWALL=`grep BSUB $case.run | grep -e '-W' `
sed s/$TIMEWALL[3]/$timewall/ < $case.run >! temp
${MOVE} temp  $case.run

set QUEUE=`grep BSUB $case.run | grep -e '-q' `
sed s/$QUEUE[3]/$queue/ < $case.run >! temp
${MOVE} temp  $case.run

chmod 0744 $case.run

# ====================================================================
# build
# ====================================================================

echo ''
echo 'Building the case'
echo ''

./$case.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit 3
endif

# ====================================================================
# Stage the required parts of DART in the execution root directory,
# now that EXEROOT exists.
# ====================================================================

foreach FILE ( filter clm_to_dart dart_to_clm )
   ${COPY} ${DARTroot}/models/clm/work/${FILE} ${exeroot}/
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/clm/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/clm/work/${FILE} not copied to ${exeroot}"
      exit 3
   endif
end

# ====================================================================
# Stage the restarts now that the run directory exists
# ====================================================================
#
# obs sequences files: /ptmp/yfzhang/Obs_seqs

# 20000501 ... /ptmp/afox/MD_40_PME/run
set stagedir = /ptmp/afox/MD_40_PME/run

# 20021101 ... /ptmp/yfzhang/inputdata_cam/lnd/clm2/initdata
# set stagedir = /ptmp/yfzhang/inputdata_cam/lnd/clm2/initdata

echo ''
echo "Copying the restart files from ${stagedir}"
echo ''

@ n = 1
while ($n <= $num_instances)

   echo "Staging restarts for instance $n of $num_instances"

#  set LANDFILE = `printf ${stagedir}/init1998.clm2_%04d.r.2002-11-01-00000.nc $n`
   set LANDFILE = `printf ${stagedir}/MD_40_PME.clm2_%04d.r.2000-01-31-00000.nc $n`
   set LND_RESTART_FILENAME = `printf "${case}.clm2_%04d.r.%04d-%02d-%02d-%05d.nc" $n $refyear $refmon $refday $run_reftod`

   ${COPY} ${LANDFILE} ${exeroot}/run/${LND_RESTART_FILENAME}

 @ n++
end

# ====================================================================
# What to do next
# ====================================================================

echo ''
echo 'Time to check the case.'
echo "cd into ${caseroot}"
echo 'Modify what you like in input.nml, make sure the observation directory'
echo 'names set in assimilate.csh match those on your system, and submit'
echo 'the CESM job by running:'
echo "./$case.submit"
echo ''

