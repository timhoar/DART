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
# -- If you have namelist mods, you need to add them to the namelist template: user_nl_cam_$case
#    Set your namelist variables over there (without modifying the syntax to create user_nl_cam_$case
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

# ====================================================================
# ====  Set case options
# ====================================================================

# ======================
# case settings
# ======================

# beta05 is also available, but i'm still using 04 for now
# because i have a fixed copy in my home directory.
setenv ccsmtag        cesm1_1_beta04
setenv case           F_PMO_3
setenv compset        F_2000
setenv resolution     f09_f09
setenv num_instances  1
setenv coldbuild      true

# had to edit the following to remove the LSB_PJL_... word too long error
# /glade/home/thoar/cesm1_1_beta04/scripts/ccsm_utils/Machines/mkbatch.bluefire
# as long as OMP_NUM_THREADS == 1 ... the default is fine.

# ================================
# define machines and directories
# ================================

# my home, and temp/scratch space
setenv nancy_home $HOME
setenv nancy_scratch $SCRATCH

setenv mach hopp2                                       ;# machine - must match cesm names
setenv DARTdir ${nancy_home}/devel

setenv cesm_datadir /project/projectdirs/ccsm1/inputdata

setenv ccsmroot ${nancy_home}/${ccsmtag}                   ;# this is where the cesm code lives
#setenv ccsmroot ${cesm_public}/collections/${ccsmtag}  ;# this is where the cesm code lives

setenv caseroot ${nancy_home}/cases/$case                  ;## your cesm case directory
setenv lockroot ${nancy_home}/locked_cases                 ;## locked cases
setenv rundir   ${nancy_scratch}/${case}
setenv archdir  ${nancy_scratch}/archive

# for data files not in the public cems area
setenv nancy_datadir ${nancy_scratch}/cesm_datafiles

# ======================
# clear out previous builds
# ======================

if ("${coldbuild}" == "true") then
   echo "removing old files from ${caseroot} and ${rundir}"
   \rm -fr ${caseroot}
   \rm -fr ${rundir}
endif

# ======================
# configure settings
# ======================

setenv run_startdate 2008-08-01

setenv sst_dataset ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010

# ======================
# runtime settings
# ======================

setenv resubmit      0
setenv stop_n        12
setenv stop_option   nhours

# ======================
# job settings
# ======================

setenv timewall     00:40:00
setenv queue        premium

# ======================
# namelist variables
# ======================
# Create namelist templates that get copied once the case has been configured.

setenv this_dir `pwd`

cat <<EOF >! user_nl_cam_${case}
&camexp
 inithist                     = 'ENDOFRUN'
 div24del2flag                = 4
 aerodep_flx_datapath         = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero'
 aerodep_flx_file             = 'aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'
 aerodep_flx_cycle_yr         = 2000
 aerodep_flx_type             = 'CYCLICAL'
 iradae                       = -12
/
EOF

# kevin had these in his cam namelist but they seem to cause
# problems with the restart files with cesm.
#  empty_htapes                 = .true.
#  nhtfrq                       = -12

cat <<EOF >! user_nl_clm_${case}
&clmexp
  fatmgrid = '${cesm_datadir}/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc'
  faerdep  = '${nancy_datadir}/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
  outnc_large_files = .true.
/
EOF

# where it would be in the cesm datadir; not on hopper:
#  faerdep  = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'

# kevin had these in his clm namelist but they seem to cause
# problems with the restart files with cesm.
#  hist_nhtfrq = -12
#  hist_empty_htapes = .true.

# ====================================================================
# ====  End of case options
# ====================================================================

# ====================================================================
# Create a new case, configure, and build.
# For list of the cases: ./create_newcase -list
# ====================================================================

if ("${coldbuild}" == "true") then
   ${ccsmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset} -skip_rundb
   
   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit 1
   endif
endif

cd ${caseroot}

./xmlchange -file env_build.xml    -id EXEROOT        -val ${rundir}
./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val ${nancy_scratch}/esmf-mpi

set num_tasks_per_instance = 12
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

./xmlchange -file env_conf.xml -id RUN_TYPE                -val startup
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_startdate
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_FILENAME   -val $sst_dataset
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_START -val $year_start
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_END   -val $year_end

#./xmlchange -file env_conf.xml     -id CLM_CONFIG_OPTS  -val '-rtm off'

./xmlchange -file env_run.xml      -id RESUBMIT         -val $resubmit
./xmlchange -file env_run.xml      -id STOP_OPTION      -val $stop_option
./xmlchange -file env_run.xml      -id STOP_N           -val $stop_n
./xmlchange -file env_run.xml      -id CALENDAR         -val GREGORIAN

# Substantial archiving changes exist in the Tools/st_archive.sh script.
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S                     -val FALSE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val FALSE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val FALSE
./xmlchange -file env_run.xml -id DOUT_L_HTAR                -val FALSE

echo ALL ARCHIVING DISABLED

# ====================================================================
# Create namelist template: user_nl_cam
# ====================================================================

\mv -f ${this_dir}/user_nl_cam_{$case} ${caseroot}/user_nl_cam
\mv -f ${this_dir}/user_nl_clm_{$case} ${caseroot}/user_nl_clm

# ====================================================================
# Update source files if need be
# ====================================================================
cp ~/cesm_mods/seq*F90           $caseroot/SourceMods/src.drv
cp ~/cesm_mods/hist*F90          $caseroot/SourceMods/src.clm
cp ~/cesm_mods/ccsm_comp_mod.F90 $caseroot/SourceMods/src.drv

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
# Stage a copy of the DART pmo.csh script HERE
# ====================================================================

cd ${caseroot}

\mv Tools/st_archive.sh Tools/st_archive.sh.org
\cp -f ${DARTdir}/models/cam/shell_scripts/st_archive.sh Tools/st_archive.sh
# only needed for beta04 - fixed in more recent versions
\cp -f ${ccsmroot}/scripts/ccsm_utils/Tools/lt_archive.csh Tools/lt_archive.csh

\cp -f ${DARTdir}/models/cam/shell_scripts/pmo.ned.csh pmo.csh

# ====================================================================
# update the namelists and scripts as needed
# ====================================================================

echo ''
echo 'Editing the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo ''

cd ${caseroot}/Buildconf

cp -f  cam.buildnml.csh  cam.buildnml.csh.org
cp -f cice.buildnml.csh cice.buildnml.csh.org
cp -f  clm.buildnml.csh  clm.buildnml.csh.org

# The CAM buildnml script only needs changing in one place.

ex cam.buildnml.csh <<ex_end
/cam_inparm/
/ncdata/
s;= '.*';= "cam_initial_1.nc";
wq
ex_end

# The CICE buildnml script only needs changing in one place.

ex cice.buildnml.csh <<ex_end
/setup_nml/
/ice_ic/
s;= '.*';= "ice_restart_1.nc";
wq
ex_end

# The CLM buildnml script needs changing in 1 place.

ex clm.buildnml.csh <<ex_end
/lnd_in/
/finidat/
s;= '.*';= "clm_restart_1.nc";
wq
ex_end

chmod 0755 clm.buildnml.csh

echo 'if you do a clean namelist, you must repeat the previous namelist update section.'
echo ' ' 

# ====================================================================
# The *.run script must be modified to call the DART perfect model script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText"
# ====================================================================

cd ${caseroot}

echo ''
echo 'Adding the call to pmo.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform a perfect model run with DART.
# -------------------------------------------------------------------------

set CplLogFile = `ls -1t cpl.log* | head -1`
if ($CplLogFile == "") then
 echo 'ERROR: Model did not complete - no cpl.log file present - exiting'
 echo 'ERROR: Perfect model run will not be attempted.'
 exit -4
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
  ${CASEROOT}/pmo.csh

  if ( $status == 0 ) then
     echo "`date` -- DART HAS FINISHED"
  else
     echo "`date` -- DART ERROR - ABANDON HOPE"
     exit -5
  endif
endif

# END OF DART BLOCK
# -------------------------------------------------------------------------

"EndOfText"

# Now that the "here" document is created, 
# determine WHERE to insert it.

set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.${mach}.run`
set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

@ orglen = `cat ${case}.${mach}.run | wc -l`
@ keep = $MYSTRING[1]
@ lastlines = $orglen - $keep 

mv ${case}.${mach}.run ${case}.${mach}.run.orig

head -$keep      ${case}.${mach}.run.orig >! ${case}.${mach}.run
cat              add_to_run.txt           >> ${case}.${mach}.run
tail -$lastlines ${case}.${mach}.run.orig >> ${case}.${mach}.run

chmod 0755 ${case}.${mach}.run

# ====================================================================
# IMPORTANT: All resubmits must be coldstarts.
# Change Tools/ccsm_postrun.csh line 83 to CONTINUE_RUN -val FALSE'
# ====================================================================

echo 'Changing Tools/ccsm_postrun.csh:83 to CONTINUE_RUN -val FALSE'
echo 'so all the resubmits are coldstarts.'

cd ${caseroot}/Tools
cp ccsm_postrun.csh ccsm_postrun.orig
sed -e '/CONTINUE_RUN/s/TRUE/FALSE/' ccsm_postrun.orig >! ccsm_postrun.csh
chmod 0755 ccsm_postrun.csh 

# if run fails, you need to back up, or you want to check what
# the current run time is, look here:
#   Buildconf/cpl.buildnml.csh: start_ymd, start_tod

# ====================================================================
# build
# ====================================================================

cd ${caseroot}

echo ''
echo 'Building the case'
echo ''

./$case.$mach.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit 3
endif

# ====================================================================
# Stage the restarts now that the run directory exists
# ====================================================================

set stagedir = ${nancy_scratch}/ned_datafiles

echo ''
echo "Staging the restarts from {$stagedir}"
echo ''

@ n = 1
while ($n <= $num_instances)

   echo "Staging restarts for instance $n of $num_instances"
   cp ${stagedir}/CAM/caminput_${n}.nc ${rundir}/run/cam_initial_${n}.nc
   cp ${stagedir}/CLM/clminput_${n}.nc ${rundir}/run/clm_restart_${n}.nc
   cp ${stagedir}/ICE/iceinput_${n}.nc ${rundir}/run/ice_restart_${n}.nc

 @ n++
end

echo 'If inflation is being used ... '
echo "must stage a ${rundir}/[prior,pos]_inflate_restart.YYYY-MM-DD-SSSSS"

# ====================================================================
# Edit the run script to reflect project, queue, and wallclock
# ====================================================================


if ($?proj) then
  set PROJ=`grep '^#PBS ' $case.$mach.run | grep -e '-P' `
  sed -e s/$PROJ[3]/$proj/ < $case.$mach.run >! temp
  /bin/mv temp  $case.$mach.run
endif

if ($?timewall) then
  set TIMEWALL=`grep '^#PBS ' $case.$mach.run | grep walltime `
  sed -e /"${TIMEWALL}"/s/=.\$/=$timewall/ < $case.$mach.run >! temp
  /bin/mv temp  $case.$mach.run
endif

if ($?queue) then
  set QUEUE=`grep '^#PBS ' $case.$mach.run | grep -e '-q' `
  sed -e s/$QUEUE[3]/$queue/ < $case.$mach.run >! temp
  /bin/mv temp  $case.$mach.run
endif
chmod 0744 $case.$mach.run

# ====================================================================
# Submit job
# ====================================================================

set MYSTRING = `grep "set DARTDIR" pmo.csh`
set DARTDIR = $MYSTRING[4]

echo ''
echo 'case is ready to submit after you check the'
echo "DART settings in ${DARTDIR}/input.nml"
echo 'After you check them,'
echo "cd into ${caseroot} and run: ./$case.$mach.submit"
echo ''

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
