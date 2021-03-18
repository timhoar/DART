#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# The intent of this script is to create a single 'CENTRALDIR' directory that
# has sub-directories for each of the TIE-GCM model instances. All instances
# should be assumed to run simultaneously, so each directory should be entirely
# independent. If there are readonly resources, those can be linked instead
# of being copied.
#
# This script should only need to be run once per experiment.
# The intended workflow is:
# 0) stage the experiment (by running this script)
# 1) configure and run (interactively) the XXXXX script to submit a series of
#    dependent jobs.
#    Job 1) a job array that will dispatch N jobs where each job will advance
#           a single ensemble member
#    Job 2) If Job1 terminates normally, Job2 will run the assimilation.
#    Job 3) same as Job 1 ...
#    Job 4) same as Job 2 ...
#    ...
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
#PBS -l walltime=0:10:00
#PBS -l wd
#PBS -l ncpus=1
#PBS -l mem=32GB
#PBS -N stage_files
#PBS -m ae
#PBS -M g.bowden@adfa.edu.au
#PBS -j oe

#===============================================================================

if ($?PBS_O_HOST) then

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $HOSTNAME

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     stage_files
   setenv JOBID       1234
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

# CENTRALDIR   Where the experiment will be run
# ENSEMBLEDIR  The location of the initial ensemble of TIEGCM files
# TIEGCMDIR    The location of the TIEGCM executable
# TGCMDATA     The location of the TIEGCM support data
# DARTDIR      The location of the DART tiegcm model directory
# OBSDIR       The location of the input observation sequence files

setenv  CENTRALDIR /glade/scratch/thoar/tiegcm_experiments/job_${JOBID}
setenv ENSEMBLEDIR /glade/work/thoar/tiegcm_ensemble
setenv   TIEGCMDIR /glade/work/thoar/git/DART_development/models/tiegcm/tiegcm
setenv    TGCMDATA ${TIEGCMDIR}/tiegcm_res5.0_data
setenv     DARTDIR /glade/work/thoar/git/DART_development/models/tiegcm
setenv      OBSDIR /glade/work/thoar/tiegcm_ensemble

mkdir -p ${CENTRALDIR}

# Archive a copy of the script in the experiment directory
set myname = $0
cp ${myname} ${CENTRALDIR}/${myname}.original

cd ${CENTRALDIR}

#----------------------------------------------------------------------
# Make a 'resource file' that has all the variables for this experiment
#----------------------------------------------------------------------

cat << EndOfText >! DART_params.csh
#!/bin/csh
# Resource file for a DART experiment

# Load required odules
# module load openmpi
# module load netcdf/4.7.3
# module load hdf5

setenv REMOVE 'rm -rf'
setenv   COPY 'cp -p'
setenv   MOVE 'mv -f'
setenv   LINK 'ln -s'

setenv  CENTRALDIR ${CENTRALDIR}
setenv   TIEGCMDIR ${TIEGCMDIR}
setenv    TGCMDATA ${TGCMDATA}
setenv     DARTDIR ${DARTDIR}
setenv ENSEMBLEDIR ${ENSEMBLEDIR}

echo "${JOBNAME} ($JOBID) CENTRALDIR == \${CENTRALDIR}"

EndOfText

# make the resource file executable
chmod 755 DART_params.csh

# We need some of those resources, this makes them available.
source DART_params.csh

#-----------------------------------------------------------------------------
# Populate the CENTRALDIR with everything that is common to all instances.
# Get the DART executables, scripts, and input files
# Get the tiegcm executable, control files, and data files.
# The individual ensemble member directories will be populated in a loop below.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/filter                                   . || exit 1
${COPY} ${DARTDIR}/work/dart_to_model                            . || exit 1
${COPY} ${DARTDIR}/work/model_to_dart                            . || exit 1
${COPY} ${DARTDIR}/work/advance_time                             . || exit 1
${COPY} ${DARTDIR}/work/fill_inflation_restart                   . || exit 1
${COPY} ${DARTDIR}/work/obs_sequence_tool                        . || exit 1
${COPY} ${DARTDIR}/work/input.nml                                . || exit 1
${COPY} ${DARTDIR}/shell_scripts/submit_multiple_cycles.csh      . || exit 1

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

# Need to convey the location of the experiment to the following scripts.
# Would not need to, but for a quirk of PBS.

echo  "s#CENTRALDIRSTRING#$CENTRALDIR#" 

sed -e "s#CENTRALDIRSTRING#$CENTRALDIR#" \
    ${DARTDIR}/shell_scripts/assimilate.csh     >! assimilate.csh     || exit 1
sed -e "s#CENTRALDIRSTRING#$CENTRALDIR#" \
    -e "s#ENSEMBLESIZESTRING#$NUM_ENS#" \
    ${DARTDIR}/shell_scripts/advance_tiegcm.csh >! advance_tiegcm.csh || exit 1

# Make these scripte executable
chmod 755 assimilate.csh advance_tiegcm.csh

${COPY} ${OBSDIR}/obs_seq.out                         obs_seq.input || exit 1
${LINK} ${ENSEMBLEDIR}/tiegcm_restart_p.nc                        . || exit 1
${LINK} ${ENSEMBLEDIR}/tiegcm_s.nc                                . || exit 1
${COPY} ${ENSEMBLEDIR}/tiegcm.nml               tiegcm.nml.original || exit 1

${COPY} ${TIEGCMDIR}/tiegcm.exec/tiegcm2.0                   tiegcm || exit 1
${COPY} ${TIEGCMDIR}/tiegcm.exec/machines.ini                     . || exit 1
${COPY} ${TIEGCMDIR}/tiegcm.exec/mpirun.command                   . || exit 1

${COPY} ${TGCMDATA}/gswm*                                         . || exit 1
${COPY} ${TGCMDATA}/wei05sc.nc                                    . || exit 1

#-----------------------------------------------------------------------------
# Make a unique directory for each model instance ... populate that
# directory with everything a TIE-GCM instance needs to run. If there
# are readonly resources, they can be 'linked' into the directory.

@ instance = 1
while ( $instance <= $NUM_ENS )

   echo "Staging ensemble member $instance of $NUM_ENS at "`date`

   set member = `printf "%04d" $instance`
   set rundir = instance_$member
   mkdir $rundir
   cd $rundir

   cp ../input.nml .

   set darticname  = `printf "filter_ics.%04d"          $instance`
   set tiesecond   = `printf "tiegcm_s.nc.%04d"         $instance`
   set tierestart  = `printf "tiegcm_restart_p.nc.%04d" $instance`

   ${COPY} -v ${ENSEMBLEDIR}/$tiesecond  tiegcm_s.nc         || exit 2
   ${COPY} -v ${ENSEMBLEDIR}/$tierestart tiegcm_restart_p.nc || exit 2

   # The gpi and imf files can get copied directly into each ensemble
   # member rundir but retain the 'simple' name without the ensemble suffix.
   # The ':r:t' retrieves just the simple name (no directory, no member).
   set WHOLE_GPINAME = ${ENSEMBLEDIR}/gpi_*.nc.$member
   set WHOLE_IMFNAME = ${ENSEMBLEDIR}/imf_*.nc.$member
   set GPINAME = $WHOLE_GPINAME:r:t
   set IMFNAME = $WHOLE_IMFNAME:r:t

   ${COPY} -v ${WHOLE_GPINAME}  ${GPINAME} || exit 2
   ${COPY} -v ${WHOLE_IMFNAME}  ${IMFNAME} || exit 2

   # Ensure that the tiegcm.nml for all the ensemble members is identical
   # in all the ways that matter.

   sed -e 's/;.*//' -e '/^$/ d' \
       -e "/ MXHIST_PRIM /c\ MXHIST_PRIM = 300" \
       -e "/ MXHIST_SECH /c\ MXHIST_SECH = 300" \
       -e "/ SOURCE /c\ SOURCE = 'tiegcm_restart_p.nc'" \
       -e "/ OUTPUT /c\ OUTPUT = 'tiegcm_restart_p.nc'" \
       -e "/ SECOUT /c\ SECOUT = 'tiegcm_s.nc'"         \
       ${ENSEMBLEDIR}/tiegcm.nml >! tiegcm.nml  || exit 2

   # If an existing ensemble of filter_ics.#### exist, use it.
   # If not, generate one. Be aware - even if they exist, they may
   # not have the same variable set as your current input.nml
   # If that is the case, you will have to generate your own set anyway.
   # If you get an error from aread_state_restart(), this is likely the case.

   if (  -e  ${ENSEMBLEDIR}/initial/$darticname.GENERATE ) then
      ${REMOVE} $darticname
      ${COPY} ${ENSEMBLEDIR}/initial/$darticname filter_ics || exit 2
   else
      # We must convert a tiegcm_restart_p.nc file to a filter_ics file
      # for each ensemble member.

      ../model_to_dart || exit 2

      if (! -e dart_ics ) then
         echo "ERROR: File conversion from $tierestart to $darticname failed."
         echo "ERROR: File conversion from $tierestart to $darticname failed."
         echo "ERROR: File conversion from $tierestart to $darticname failed."
         exit 2
      else
         # Rename as expected input for filter ... 
         ${MOVE} dart_ics ../$darticname
      endif
   endif

   cd ..

   @ instance++

end

echo "Finished staging the experiment at "`date`

#-------------------------------------------------------------------------------
# Create some instructions for what do do next.

cat << EndOfText >! README.txt

# REQUIREMENTS: for input.nml
#-----------------------------------------------------------------------------
# obs_sequence_tool_nml: filename_seq               = 'obs_seq.input'
# obs_sequence_tool_nml: filename_out               = 'obs_seq.out'
#-----------------------------------------------------------------------------
# model_nml            : tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
# model_nml            : tiegcm_secondary_file_name = 'tiegcm_s.nc'
# model_nml            : tiegcm_namelist_file_name  = 'tiegcm.nml'
# model_to_dart_nml    : file_out                   = 'dart_ics'
#-----------------------------------------------------------------------------
# ensemble_manager_nml : single_restart_file_in     = .false.
# filter_nml           : async                      = 0
# filter_nml           : adv_ens_command            = './no_advance_model.csh'
# filter_nml           : start_from_restart         = .TRUE.
# filter_nml           : restart_in_file_name       = 'filter_ics'
# filter_nml           : restart_out_file_name      = 'filter_restart'
#-----------------------------------------------------------------------------
# dart_to_model_nml    : file_in                    = 'dart_restart'
# dart_to_model_nml    : file_namelist_out          = 'namelist_update'

1) cd ${CENTRALDIR}
2) check namelists
3) ./submit_multiple_cycles 3

EndOfText

cat README.txt

exit 0
