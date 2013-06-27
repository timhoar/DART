#!/bin/csh

set RUN_DIR = /net/scud-ib/s1/arellano/DART/rundir_gabi
set SRC_DIR = /cu1/arellano/DARTvK/DART/models/wrf/work	# filter execution files
set SCRIPT_DIR = /net/scud-ib/s1/arellano/DART/rundir_C		# scripts
set ICBC_OUT_DIR = /net/scud-ib/s1/arellano/DART/datadir/icbc/icbc2/icbc3/icbc4/output
set MEGAN_FILE_DIR = /net/scud-ib/s1/arellano/TEMP/EMIS/code/MEGAN/gabi
set WRF_RUN_DIR = /net/scud-ib/s1/arellano/DART/rundir_gabi/WRF_RUN

set expn    = June08_Gabi
set nens    = 40
set start_year=2008
set start_month=6
set start_day=1
set start_hour=0

set mon0=`printf %02d $start_month`
set day0=`printf %02d $start_day`
set  hr0=`printf %02d $start_hour`
set curhr=`echo ${start_year}${mon0}${day0}${hr0} 0 | ./advance_time`
set  year0=`echo $curhr | cut -b1-4`
set month0=`echo $curhr | cut -b5-6`
set   day0=`echo $curhr | cut -b7-8`
set  hour0=`echo $curhr | cut -b9-10`

# get gregday and gregsec
set g=(`echo ${year0}${month0}${day0}${hour0} 0 -g | ./advance_time`)
set gregday=$g[1]
set gregsec=$g[2]

# change directory
cd ${RUN_DIR}

# prepare scripts
#cp ${SCRIPT_DIR}/advance_model.csh .
#cp ${SCRIPT_DIR}/advance_wrf.csh .
#cp ${SCRIPT_DIR}/run_filter_mpiwrf.csh .
#cp ${SCRIPT_DIR}/step09_param.csh .
#cp ${SCRIPT_DIR}/filter.template.lsf .
#cp ${SCRIPT_DIR}/advance_model.template.lsf .
#cp ${SCRIPT_DIR}/save_files.lsf .
#cp ${SCRIPT_DIR}/input.nml.template.0 .
#cp ${SCRIPT_DIR}/input.nml.template.1 .
#cp ${SCRIPT_DIR}/namelist.input .

# other executables?
#ln -sf ${SRC_DIR}/convertdate .
#ln -sf ${SRC_DIR}/restart_file_tool .
#ln -sf ${SRC_DIR}/update_wrf_bc .
#
# assumes that LBCs and WRFINPUTs are in ./WRF directory

#
## Set up for the first cycle
set cycle_0_dir = $RUN_DIR/${expn}/Cycle_0
if(! -d $cycle_0_dir) mkdir -p $cycle_0_dir

# loop
set ie = 1
while ( $ie <= $nens )
   set ensstring = `echo $ie + 10000 | bc | cut -b2-5`
   ln -sf $ICBC_OUT_DIR/filter_ic.${ensstring} ${cycle_0_dir}/filter_ic_old.${ensstring}
   ln -sf $ICBC_OUT_DIR/wrfinput_d01_${gregday}_${gregsec}_${ie} ${cycle_0_dir}/wrfinput_d01_${ie}
   @ ie ++
end
cp ./WRF/wrfinput_d01 .

# link biochem files --this is sloppy for now
# we preprocess MEGAN (although it should be static) after running
# icbc_real.ksh --since we need the wrfinput_d01s
#ln -sf ${MEGAN_FILE_DIR}/wrfbiochemi_d01_* $WRF_RUN_DIR/

echo 'done prep'
