#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
##------------------------------------------------------------------------------
## This block of directives constitutes the preamble for the PBS queuing system
##
## the normal way to submit to the queue is:    qsub run_filter.csh
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error
## -o <arg>  filename for standard out
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok
##                     and calgary, there is no way to 'share' the processors
##                     on the node with another job, so you might as well use
##                     them both. (ppn == Processors Per Node)
##
#PBS -P n23
#PBS -l walltime=3:00:00
#PBS -l wd
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -N run_tiegcm_assim_gadi
#PBS -m ae
#PBS -M g.bowden@adfa.edu.au
#PBS -j oe

#===============================================================================

# Load required modules
module load openmpi
module load netcdf/4.7.3
module load hdf5

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

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     tiegcm_filter
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI_RUN_CMD ''

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

#setenv CENTRALDIR /scratch/n23/gwb112/swm_project/DART/dart_tiegcm/${JOBNAME}/job_${JOBID}
cd ../../..
setenv CENTRALDIR `pwd`/${JOBNAME}/job_${JOBID}

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
      setenv   LINK 'ln -s'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      setenv   LINK 'ln -s'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -v'
      setenv   MOVE 'mv -fv'
      setenv   LINK 'ln -s'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) CENTRALDIR == $CENTRALDIR"

# Set variables containing various directory names where we will GET things
# DARTDIR      The location of the DART tiegcm model directory
# TIEGCMDIR    The location of the TIEGCM executable
# ENSEMBLEDIR  The location of the initial ensemble of TIEGCM files
# EXPERIMENT   The (safe) location for the results of this run.
#-----------------------------------------------------------------------------

set    DARTDIR = /scratch/n23/gwb112/swm_project/DART/dart_tiegcm/models/tiegcm
set  TIEGCMDIR = /scratch/n23/gwb112/swm_project/TIEGCM
set TIEGCMDATA = ${TIEGCMDIR}/tiegcm_res5.0_data
#set EXPERIMENT = /scratch/n23/gwb112/swm_project/gold_test/initial
#set ENSEMBLEDIR = /scratch/n23/gwb112/swm_project/gold_test/initial
#set EXPERIMENT = /scratch/n23/gwb112/swm_project/gold_test_restart/initial
#set ENSEMBLEDIR = /scratch/n23/gwb112/swm_project/gold_test_restart/initial
#set EXPERIMENT = /scratch/n23/gwb112/swm_project/gold_test_large/initial
#set ENSEMBLEDIR = /scratch/n23/gwb112/swm_project/gold_test_large/initial
#set EXPERIMENT = /scratch/n23/gwb112/swm_project/DART/dart_tiegcm/run_tiegcm_assim_gadi/job_14647484.gadi-pbs
#set ENSEMBLEDIR = /scratch/n23/gwb112/swm_project/DART/dart_tiegcm/run_tiegcm_assim_gadi/job_14647484.gadi-pbs
set EXPERIMENT = /scratch/n23/gwb112/swm_project/DART/dart_tiegcm/run_tiegcm_assim_gadi/job_15037203.gadi-pbs 
set ENSEMBLEDIR = /scratch/n23/gwb112/swm_project/DART/dart_tiegcm/run_tiegcm_assim_gadi/job_15037203.gadi-pbs

# Need to set TGCMDATA environment variable
setenv TGCMDATA $TIEGCMDATA

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the tiegcm executable, control files, and data files.
# The tiegcm initial conditions are in the next block.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/filter                       . || exit 1
${COPY} ${DARTDIR}/work/dart_to_model                . || exit 1
${COPY} ${DARTDIR}/work/model_to_dart                . || exit 1
${COPY} ${DARTDIR}/work/input.nml   input.nml.original || exit 1
${COPY} ${DARTDIR}/shell_scripts_gadi/advance_model.csh       . || exit 1
${COPY} ${DARTDIR}/work/wakeup_filter                . || exit 1

#${COPY} ${EXPERIMENT}/observation/obs_seq.out        . || exit 1
${COPY} ${EXPERIMENT}/obs_seq.out            . || exit 1
${COPY} ${EXPERIMENT}/tiegcm_restart_p.nc    . || exit 1
${COPY} ${EXPERIMENT}/tiegcm_s.nc            . || exit 1
${COPY} ${EXPERIMENT}/tiegcm.nml   tiegcm.nml.original || exit 1
${COPY} ${EXPERIMENT}/gpi*.nc.*              . || exit 1
${COPY} ${EXPERIMENT}/imf*.nc.*              . || exit 1

${COPY} ${TIEGCMDIR}/tiegcm.exec/tiegcm2.0      tiegcm || exit 1
${COPY} ${TIEGCMDIR}/tiegcm.exec/machines.ini        . || exit 1
${COPY} ${TIEGCMDIR}/tiegcm.exec/mpirun.command      . || exit 1

${COPY} ${TIEGCMDATA}/gswm*                          . || exit 1
${COPY} ${TIEGCMDATA}/wei05sc.nc                     . || exit 1

#-----------------------------------------------------------------------------
# Put all of the DART initial conditions files and all of the TIEGCM files
# in the CENTRALDIR - preserving the ensemble member ID for each filename.
# The advance_model.csh script will copy the appropriate files for each 
# ensemble member into the model advance directory.
# These files may be linked to CENTRALDIR since they get copied to the
# model advance directory. 
#
# REQUIREMENTS: for input.nml
# model_nml            : tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
# model_nml            : tiegcm_secondary_file_name = 'tiegcm_s.nc'
# model_nml            : tiegcm_namelist_file_name  = 'tiegcm.nml'
# model_to_dart_nml    : file_out                   = 'dart_ics'
#-----------------------------------------------------------------------------
# ensemble_manager_nml : single_restart_file_in     = .false.
# filter_nml           : async                      = 4
# filter_nml           : adv_ens_command            = './advance_model.csh'
# filter_nml           : start_from_restart         = .TRUE.
# filter_nml           : restart_in_file_name       = 'filter_ics'
# filter_nml           : restart_out_file_name      = 'filter_restart'
#-----------------------------------------------------------------------------
# dart_to_model_nml    : file_in                    = 'dart_restart'
# dart_to_model_nml    : file_namelist_out          = 'namelist_update'

sed -e "/ tiegcm_restart_file_name /c\ tiegcm_restart_file_name = 'tiegcm_restart_p.nc'" \
    -e "/ tiegcm_secondary_file_name /c\ tiegcm_secondary_file_name = 'tiegcm_s.nc'" \
    -e "/ tiegcm_namelist_file_name /c\ tiegcm_namelist_file_name = 'tiegcm.nml'" \
    -e "/ file_out /c\ file_out = 'dart_ics'" \
    -e "/ single_restart_file_in /c\ single_restart_file_in = .FALSE." \
    -e "/ async /c\ async = 4" \
    -e "/ adv_ens_command /c\ adv_ens_command = './advance_model.csh'" \
    -e "/ start_from_restart /c\ start_from_restart = .TRUE." \
    -e "/ restart_in_file_name /c\ restart_in_file_name = 'filter_ics'" \
    -e "/ restart_out_file_name /c\ restart_out_file_name = 'filter_restart'" \
    -e "/ file_in /c\ file_in = 'dart_restart'" \
    -e "/ file_namelist_out /c\ file_namelist_out = 'namelist_update'" \
    input.nml.original >! input.nml  || exit 2

# Determine the number of ensemble members from input.nml,
# it may exist in more than one place.
# Parse out the filter_nml string and see which 
# one is immediately after it ...

if ( ! -e input.nml ) then
   echo "ERROR - input.nml does not exist in local directory."
   echo "ERROR - input.nml needed to determine number of ensemble members."
   exit 1
endif

set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set NUM_ENS = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`

# FIXME ... read the async value from input.nml and set parallel_model accordingly.
# if async=2, e.g. you are going to run './modelxxx', single process
# (or possibly 'mpirun -np 1 ./modelxxx'), so each processor advances
# one ensemble independently of the others, leave this as false.
#
# if async=4, e.g. all the processors advance each modelxxx in turn with
# mpirun -np 64 modelxxx (or whatever) for as many ensembles as you have,
# set this to "true"

set parallel_model = "true"

# A common strategy for the beginning is to check for the existence of
# some variables that get set by the different queuing mechanisms.
# This way, we know which queuing mechanism we are working with,
# and can set 'queue-independent' variables for use for the remainder
# of the script.

if ($?LS_SUBCWD) then

    # LSF has a list of processors already in a variable (LSB_HOSTS)
    # alias submit 'bsub < \!*'
    echo "LSF - using mpirun.lsf for execution"
    setenv MPICMD mpirun.lsf

else if ($?PBS_O_WORKDIR) then

    # PBS has a list of processors in a file whose name is (PBS_NODEFILE)
    # alias submit 'qsub \!*'
    echo "PBS - using mpirun for execution"
    setenv MPICMD mpirun

else

    # If you have a linux cluster with no queuing software, use this
    # section. The list of computational nodes is given to the mpirun
    # command and it assigns them as they appear in the file. In some
    # cases it seems to be necessary to wrap the command in a small
    # script that changes to the current directory before running.

    echo "running with no queueing system"

    # before running this script, do this once. the syntax is
    # node name : how many tasks you can run on it
    #setenv MYNODEFILE ~/nodelist
    #echo "node7:2" >! $MYNODEFILE
    #echo "node5:2" >> $MYNODEFILE
    #echo "node3:2" >> $MYNODEFILE
    #echo "node1:2" >> $MYNODEFILE

#   one possibility
    setenv NUM_PROCS `cat nodelist-pgi | wc -l`
    set MPIRUN = /apps/openmpi/4.0.2/bin/mpirun
    set MPICMD = $MPIRUN -np $NUM_PROCS -nolocal -machinefile nodelist-pgi

    echo "MPICMD = ${MPICMD}"

endif

@ instance = 1
while ( $instance <= $NUM_ENS )

  set darticname  = `printf "filter_ics.%04d"          $instance`
  set tiesecond   = `printf "tiegcm_s.nc.%04d"         $instance`
  set tierestart  = `printf "tiegcm_restart_p.nc.%04d" $instance`
  set tieinp      = `printf "tiegcm.nml.%04d"          $instance`

  ${COPY} ${ENSEMBLEDIR}/$tiesecond  .                   || exit 2
  ${COPY} ${ENSEMBLEDIR}/$tierestart .                   || exit 2
  ${COPY} ${ENSEMBLEDIR}/$tieinp     tiegcm.nml.original || exit 2

  # Ensure that the tiegcm.nml for all the ensemble members is identical
  # in all the ways that matter. This will result in a miniumum of changes
  # in the advance_model.csh script. This script REQUIRES that there is a  
  # SINGLE tiegcm_restart_p.nc. Just keep appending all the timesteps to
  # the same file. If you need to subset the large file, use the NCO
  # operators. for example    ncks -d time,20,30 tiegcm_restart_p.nc bob.nc 
  # If you need more than 300 timesteps in the file, increase it here.
  
  sed -e 's/;.*//' -e '/^$/ d' \
      -e "/ MXHIST_PRIM /c\ MXHIST_PRIM = 300" \
      -e "/ MXHIST_SECH /c\ MXHIST_SECH = 300" \
      -e "/ SOURCE /c\ SOURCE = 'tiegcm_restart_p.nc'" \
      -e "/ OUTPUT /c\ OUTPUT = 'tiegcm_restart_p.nc'" \
      -e "/ SECOUT /c\ SECOUT = 'tiegcm_s.nc'"         \
      tiegcm.nml.original >! $tieinp  || exit 2

  # If an existing ensemble of filter_ics.#### exist, use it.
  # If not, generate one. Be aware - even if they exist, they may
  # not have the same variable set as your current input.nml
  # If that is the case, you will have to generate your own set anyway.
  # If you get an error from aread_state_restart(), this is likely the case.

  if (  -e  ${ENSEMBLEDIR}/initial/$darticname.GENERATE ) then
     ${REMOVE} $darticname
     ${LINK} ${ENSEMBLEDIR}/initial/$darticname . || exit 2
  else
     # We must convert a tiegcm_restart_p.nc file to a dart_ics file
     # for each ensemble member. So - momentarily, we must
     # create links to the static filenames expected by model_to_dart

     ${REMOVE} tiegcm_restart_p.nc tiegcm_s.nc tiegcm.nml 

     ${LINK} $tierestart tiegcm_restart_p.nc   || exit 2
     ${LINK} $tiesecond  tiegcm_s.nc           || exit 2
     ${LINK} $tieinp     tiegcm.nml            || exit 2

     ./model_to_dart || exit 2

     if (-e dart_ics ) then
        ${MOVE} dart_ics $darticname
     else
        echo "ERROR: File conversion from $tierestart to $darticname failed."
        echo "ERROR: File conversion from $tierestart to $darticname failed."
        echo "ERROR: File conversion from $tierestart to $darticname failed."
        exit 2
     endif
  endif

  @ instance++
end

#-------------------------------------------------------------------------------
# Everything below this separator should not need to be modified if everything
# above the separator is set correctly.
#-------------------------------------------------------------------------------

if ( "$parallel_model" == "false" ) then

   # each filter task advances the ensembles, each running on 1 proc.

   ${MPICMD} ./filter

else

   # filter runs in parallel until time to do a model advance,
   # and then this script starts up the modelxxx jobs, each one
   # running in parallel. then it runs wakeup_filter to wake
   # up filter so it can continue. The communication happens through
   # 'named pipes' created by the mkfifo command.

   \rm -f model_to_filter.lock filter_to_model.lock
   mkfifo model_to_filter.lock filter_to_model.lock

   set filterhome = ~/.filter$$
   if ( ! -e $filterhome) mkdir $filterhome

   # start filter and immediately return control back to this script

   (setenv HOME $filterhome; ${MPICMD} ./filter) &

   while ( -e filter_to_model.lock )

      set todo=`cat < filter_to_model.lock`
      echo "todo received, value = ${todo}"

      if ( "${todo}" == "finished" ) then
         echo "main script: filter done."
         wait
         break

      else if ( "${todo}" == "advance" ) then

         # FIXME : in input.nml, the advance model command must
         # have -np N with N equal to the number of processors this job is using.

         echo "calling model advance now:"
         ./advance_model.csh 0 ${NUM_ENS} filter_control00000 || exit 9

         echo "restarting filter."
         ${MPICMD} ./wakeup_filter

      else
          echo "main script: unexpected value received."
          break
      endif
   end

   echo "filter finished, removing pipes."
   \rm -f model_to_filter.lock filter_to_model.lock

   if ( -d $filterhome) \rmdir $filterhome

endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

