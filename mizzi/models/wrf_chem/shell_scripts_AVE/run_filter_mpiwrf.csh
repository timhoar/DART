#!/bin/csh -x
###############################################################################
#
#  run_filter_mpiwrf.csh
#
# Experiment settings 
set EXP_NAME = AVE_ORIG_June08
set HSI_NAME = /MIZZI/${EXP_NAME}
set NUM_MEMBERS = 3
set DOMAIN = 1
set CYCLE_PERIOD = 6
set CYCLE = 49
set NUM_CYCLES = 49
set DATE_START = 2008061300
set DATE_END = 2008061306
#
# Task to run
# Need RUN_COPY_FILES for cold start only
set RUN_COPY_FILES = false
# Need RUN_SETUP for cold and warm start.
set RUN_SETUP = true
set RUN_FILTER = false
set RUN_ADVANCE_MODEL = true
#
# Set DART inflation files name
set RESTART_IN_FILE_NAME = 'filter_ic_old'
set RESTART_OUT_FILE_NAME = 'filter_ic_new'
set INF_IN_FILE_NAME_PRIOR = 'prior_inflate_ic_old'
set INF_IN_FILE_NAME_POST = 'post_inflate_ics'
set INF_OUT_FILE_NAME_PRIOR = 'prior_inflate_ic_new'
set INF_OUT_FILE_NAME_POST = 'post_inflate_restart'
set INF_DIAG_FILE_NAME_PRIOR = 'prior_inflate_diag'
set INF_DIAG_FILE_NAME_POST = 'post_inflate_diag'
set RESTART_TOOL_IN_FILE_NAME = 'filter_ic_new'
set RESTART_TOOL_OUT_FILE_NAME = 'assim_model_state_ic'
#
# Job submission settings
set PROJ_NUMBER = 25000079
set JOB_TIME_FILTER = 02:20
set JOB_TIME_WRFCHEM = 02:20
set NUM_TASKS = 64
set TASKS_PER_NODE = 32
set JOB_CLASS = regular
# 
# Code version settings
set DART_VER = DART_DEVEL
set WRFCHEM_VER = WRFCHEMv3.4_dmpar
set WRF_VER = WRFv3.4_dmpar
set WRFDA_VER = WRFDAv3.4_dmpar
#
# Independent path settings
set TRUNK_DIR = /glade/home/mizzi/TRUNK   
set DAT_DIR = /ptmp/mizzi/AVE_ORIG_DATA
set DAT_TEST_DIR = /ptmp/mizzi/AVE_TEST_DATA
#
# Dependent path settings
set DART_DIR = ${TRUNK_DIR}/${DART_VER}
set WRF_DIR = ${TRUNK_DIR}/${WRF_VER}
set WRFCHEM_DIR = ${TRUNK_DIR}/${WRFCHEM_VER}
set WRFDA_DIR = ${TRUNK_DIR}/${WRFDA_VER}/var
set RUN_DIR = /ptmp/mizzi/${EXP_NAME}
set WRF_RUN_DIR = ${RUN_DIR}/WRF_RUN
set OBS_DIR = ${DAT_DIR}/OBS_SETS
set WRF_DAT_DIR = ${DAT_DIR}/WRF
#
# Set commands
set REMOVE = '/bin/rm -rf'
set COPY = 'cp -pf'
set MOVE = 'mv -f'
unalias cd
unalias ls
#
# BEGIN RUN_COPY_FILES CODE BLOCK
if (${RUN_COPY_FILES} == true) then
#
# Create ${RUN_DIR}
   if (! -d ${RUN_DIR}) then
      mkdir -p ${RUN_DIR}
      cd ${RUN_DIR}
   else
      cd ${RUN_DIR}
   endif
#
# Copy DART, WRF, WRFCHEM, and WRFDA executables to ${RUN_DIR}
   cp ${DART_DIR}/models/wrf/work/dart_to_wrf ./.
   cp ${DART_DIR}/models/wrf/work/wrf_to_dart ./.
   cp ${DART_DIR}/models/wrf/work/update_wrf_bc ./.
   cp ${DART_DIR}/models/wrf/work/pert_wrf_bc ./.
   cp ${DART_DIR}/models/wrf/work/filter ./.
   cp ${DART_DIR}/models/wrf/work/wakeup_filter ./.
   cp ${DART_DIR}/models/wrf/work/obs_diag ./.
   cp ${DART_DIR}/models/wrf/work/obs_seq_to_netcdf ./.
   cp ${DART_DIR}/models/wrf/work/restart_file_tool ./.
   cp ${DART_DIR}/models/wrf/work/advance_time ./.
   cp ${DART_DIR}/models/wrf/work/convertdate ./.
   cp ${DART_DIR}/models/wrf_chem/templates/input.nml.template ./.
   cp ${DART_DIR}/models/wrf_chem/templates/input.nml.template input.nml
   cp ${DART_DIR}/models/wrf_chem/templates/namelist.input.template ./.
   cp ${DART_DIR}/models/wrf_chem/templates/namelist.input.template ./namelist.input
   cp ${DART_DIR}/models/wrf_chem/templates/advance_model.template.lsf ./.
   cp ${DART_DIR}/models/wrf_chem/templates/wrfinput_d01 ./. 
   if(! -e input.nml.template || ! -e wrfinput_d01 || ! -e namelist.input.template) then
      echo ERROR: CHECK WHETHER input.nml.template, wrfinput_d01, AND namelist.input EXIST.
      exit
   endif
#
# Copy DART, WRF, WRFCHEM, and WRFDA scripts to ${RUN_DIR}
   cp ${DART_DIR}/models/wrf_chem/shell_scripts/advance_model.csh ./.
#
# Create ${WRF_RUN_DIR} - This is needed for interfacing with advance_model.csh 
   if (! -d ${WRF_RUN_DIR}) then
      mkdir -p ${WRF_RUN_DIR}
      cd ${WRF_RUN_DIR}
   else
      cd ${WRF_RUN_DIR}
   endif
#
# Copy WRF, WRFCHEM, and WRFDA files to ${WRF_RUN_DIR} 
   cp ${WRFDA_DIR}/build/da_wrfvar.exe ./.
   cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ./.
   cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ../.
   cp ${WRFCHEM_DIR}/test/em_real/CAM_ABS_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/CAM_AEROPT_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/GENPARM.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/LANDUSE.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/RRTMG_LW_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/RRTMG_SW_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/RRTM_DATA ./.
   cp ${WRFCHEM_DIR}/test/em_real/SOILPARM.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/URBPARM.TBL ./.
   cp ${WRFCHEM_DIR}/test/em_real/VEGPARM.TBL ./.
#
   cp ${DAT_TEST_DIR}/chem_static/clim_p_trop.nc ./.
   cp ${DAT_TEST_DIR}/chem_static/exo_coldens_d01 ./.
   cp ${DAT_TEST_DIR}/chem_static/hist_io_mods_moz_d01 ./.
   cp ${DAT_TEST_DIR}/chem_static/ubvals_b40.20th.track1_1996-2005.nc ./.
   cp ${DAT_TEST_DIR}/chem_static/wrf_season_wes_usgs_d01.nc ./.
   cp ${DAT_TEST_DIR}/chem_static/wrfbiochemi_d01_2008-06-13_00:00:00 ./.
   cp ${DAT_TEST_DIR}/chem_static/wrfbiochemi_d01_2008-06-13_06:00:00 ./.
   cp ${DAT_TEST_DIR}/chem_static/20080613/wrfchemi_d01_2008-06-13_**:00:00 ./.
   cp ${DAT_TEST_DIR}/chem_static/20080613/wrffirechemi_d01_2008-06-13_**:00:00 ./.
endif
# END RUN_COPY_FILES CODE BLOCK   
#
# BEGIN SETUP CODE BLOCK
if (${RUN_SETUP} == true) then
   cd ${RUN_DIR}
   if (! -d ${RUN_DIR}/WRF) then 
      mkdir -p ${RUN_DIR}/WRF
   endif
#
# APM: copy in wrfbdy files - this is temp fix
   cd ${RUN_DIR}/WRF
   rm wrfbdy*
   set iens = 1
   while ( $iens <= ${NUM_MEMBERS} )
      set imem=${iens}
      if( ${iens} < 1000 ) then
         set imem=0${iens}
      endif
      if( ${iens} < 100 ) then
         set imem=00${iens}
      endif
      if( ${iens} < 10 ) then
         set imem=000${iens}
      endif
      set jmem=${iens}
      if( ${iens} < 100 ) then
         set jmem=0${iens}
      endif
      if( ${iens} < 10 ) then
         set jmem=00${iens}
      endif
      cp ${DAT_TEST_DIR}/wpb_rc/2008061300/wrfbdy_d01_2008-06-13_00:00:00.e${jmem} wrfbdy_148817_21600_${iens}
      @ iens ++
   end
   cd ${RUN_DIR}
# APM: begin temp copies
#   rm -rf advance_model.template.lsf
#   rm -rf advance_model.csh
#   rm -rf namelist.input
#   cp ${DART_DIR}/models/wrf_chem/shell_scripts/advance_model.template.lsf ./.
#   cp ${DART_DIR}/models/wrf_chem/shell_scripts/advance_model.csh ./.
#   cp ${DART_DIR}/models/wrf_chem/namelist/namelist.input.template ./namelist.input
#   cp ${DART_DIR}/models/wrf_chem/work/input.nml.template ./.
#   cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ./.
   cd ${WRF_RUN_DIR}
   cp ${DAT_TEST_DIR}/chem_static/clim_p_trop.nc ./.
   cp ${DAT_TEST_DIR}/chem_static/exo_coldens_d01 ./.
   cp ${DAT_TEST_DIR}/chem_static/hist_io_mods_moz_d01 ./.
   cp ${DAT_TEST_DIR}/chem_static/ubvals_b40.20th.track1_1996-2005.nc ./.
   cp ${DAT_TEST_DIR}/chem_static/wrf_season_wes_usgs_d01.nc ./.
   cp ${DAT_TEST_DIR}/chem_static/wrfbiochemi_d01_2008-06-13_00:00:00 ./.
   cp ${DAT_TEST_DIR}/chem_static/wrfbiochemi_d01_2008-06-13_06:00:00 ./.
   cp ${DAT_TEST_DIR}/chem_static/20080613/wrfchemi_d01_2008-06-13_**:00:00 ./.
   cp ${DAT_TEST_DIR}/chem_static/20080613/wrffirechemi_d01_2008-06-13_**:00:00 ./.
   cd ${RUN_DIR}
# APM: end temp copies
#
# Set script parameters
   set expn = ${EXP_NAME}
   set nens = ${NUM_MEMBERS}
   set ndom = ${DOMAIN}
   set cutoff = 0.1
   set adaptive_inf = true
   set vcoord_local = 3
   set obs_pre = obs_seq
#
# Set HSI parameters
   set sav_mss = false
   if($sav_mss == true) then
      set  mdir = ${EXP_NAME}
      set  mssdir = $mdir
      set  MSRCP = 'msrcp -pe 3650 -proj '${PROJ_NUMBER}
   endif
#
# Set cycling and timing parameters
   set icyc_start = 1
   set icyc = ${CYCLE}
#
   set YY = `echo $DATE_START | cut -c1-4`
   set MM = `echo $DATE_START | cut -c5-6`
   set DD = `echo $DATE_START | cut -c7-8`
   set HH = `echo $DATE_START | cut -c9-10`
   set init_dat = ${YY}-${MM}-${DD}_${HH}:00:00
#
   set YY = `echo $DATE_END | cut -c1-4`
   set MM = `echo $DATE_END | cut -c5-6`
   set DD = `echo $DATE_END | cut -c7-8`
   set HH = `echo $DATE_END | cut -c9-10`
   set last_dat = ${YY}-${MM}-${DD}_${HH}:00:00
#
   set init_grg = `echo $init_dat 0 -g | ./advance_time`
   set last_grg = `echo $last_dat 0 -g | ./advance_time`
   set intv_scs = `expr ${CYCLE_PERIOD} \* 3600`
   set intv_hrs = `expr $intv_scs \/ 3600`
   set diff_day = `expr $last_grg[1] \- $init_grg[1]`
   set diff_scs = `expr $last_grg[2] \- $init_grg[2]`
   set diff_tot = `expr $diff_day \* 86400 \+ $diff_scs`
   set n_cycles = `expr $diff_tot \/ $intv_scs`
   set ncyc = `expr $icyc \+ $n_cycles \- 1`
#
   set obs_time_days = $init_grg[1]
   set obs_time_secs = $init_grg[2]
   set cycle_restart = 0
   echo Initial obs time: $obs_time_days $obs_time_secs
endif
# END RUN_SETUP CODE BLOCK
#
# Begin cycling loop
while ( $icyc <= $ncyc )
   set iens =  1
   if($icyc == $icyc_start) then
      set prev_cyc = `expr $icyc \- 1`
      set old_dir = $DAT_DIR/$expn/Cycle_${prev_cyc}
   else
      set prev_cyc = `expr $icyc \- 1`
      set old_dir = $DAT_DIR/$expn/Cycle_${prev_cyc}
   endif
#
# BEGIN RUN_FILTER CODE BLOCK
   if (${RUN_FILTER} == true) then
      ${REMOVE} script.sed
      ${REMOVE} script1.sed
      ${REMOVE} script2.sed
#
# Prepare input.nml
      if ( $obs_time_secs == 0 ) then
         set obs_time_days1 = `expr $obs_time_days \- 1`
         set obs_time_secs1 = `expr 86400 \- $intv_scs \/ 2 \+ 1`
      else
         set obs_time_days1 = $obs_time_days
         set obs_time_secs1 = `expr $obs_time_secs \- $intv_scs \/ 2 \+ 1`
      endif
      set obs_time_secs2 = `expr $obs_time_secs \+ $intv_scs \/ 2`
#
      cat >! script.sed << EOF
      /ens_size /c\
      ens_size                 = $nens,
      /first_obs_days /c\
      first_obs_days           = $obs_time_days1,
      /first_obs_seconds /c\
      first_obs_seconds        = $obs_time_secs1,
      /last_obs_days /c\
      last_obs_days            = $obs_time_days,
      /last_obs_seconds /c\
      last_obs_seconds         = $obs_time_secs2,
      /cutoff  /c\
      cutoff                   = $cutoff,
      /num_domains /c\
      num_domains              = $ndom,
      /assimilation_period_seconds/c\
      assimilation_period_seconds = $intv_scs,
      /vert_localization_coord /c\
      vert_localization_coord     = $vcoord_local,
EOF
#
      if( $adaptive_inf == true ) then
         cat >! script1.sed << EOF
         /inf_flavor/c\
         inf_flavor = 2, 0,
EOF
         if($icyc == $icyc_start) then
            cat >! script2.sed << EOF
            /inf_initial_from_restart/c\
            inf_initial_from_restart = .false., .false.,
            /inf_sd_initial_from_restart/c\
            inf_sd_initial_from_restart = .false., .false.,
            /inf_in_file_name/c\
            inf_in_file_name = '${INF_IN_FILE_NAME_PRIOR}', '${INF_IN_FILE_NAME_POST}'
            /inf_out_file_name/c\
            inf_out_file_name = '${INF_OUT_FILE_NAME_PRIOR}', '${INF_OUT_FILE_NAME_POST}'
            /inf_diag_file_name/c\
            inf_diag_file_name = '${INF_DIAG_FILE_NAME_PRIOR}', '${INF_DIAG_FILE_NAME_POST}'
            /restart_in_file_name/c\
            restart_in_file_name = '${RESTART_IN_FILE_NAME}'
            /restart_out_file_name/c\
            restart_out_file_name = '${RESTART_OUT_FILE_NAME}'
EOF
         else
            cat >! script2.sed << EOF
            /inf_initial_from_restart/c\
            inf_initial_from_restart = .true., .true.,
            /inf_sd_initial_from_restart/c\
            inf_sd_initial_from_restart = .true., .true.,
            /inf_in_file_name/c\
            inf_in_file_name = '${INF_IN_FILE_NAME_PRIOR}', '${INF_IN_FILE_NAME_POST}'
            /inf_out_file_name/c\
            inf_out_file_name = '${INF_OUT_FILE_NAME_PRIOR}', '${INF_OUT_FILE_NAME_POST}'
            /inf_diag_file_name/c\
            inf_diag_file_name = '${INF_DIAG_FILE_NAME_PRIOR}', '${INF_DIAG_FILE_NAME_POST}'
            /restart_in_file_name/c\
            restart_in_file_name = '${RESTART_IN_FILE_NAME}'
            /restart_out_file_name/c\
            restart_out_file_name = '${RESTART_OUT_FILE_NAME}'
EOF
         endif
         cat script1.sed >> script.sed
         cat script2.sed >> script.sed
      else	
         cat >! script1.sed << EOF
         /inf_flavor /c\
         inf_flavor = 0, 0,
EOF
         cat script1.sed >> script.sed
      endif
#
      if ( -e input.nml) ${REMOVE} input.nml
      sed -f script.sed input.nml.template >! input.nml
#
# Run filter
      if( -e ${RESTART_OUT_FILE_NAME}.0001) ${REMOVE} ${RESTART_OUT_FILE_NAME}.*
      if($icyc == $icyc_start) then
         while ( $iens <= $nens )
            set imem=${iens}
            if( ${iens} < 1000 ) then
               set imem=0${iens}
            endif
            if( ${iens} < 100 ) then
               set imem=00${iens}
            endif
            if( ${iens} < 10 ) then
               set imem=000${iens}
            endif
            ln -sf ${old_dir}/${RESTART_IN_FILE_NAME}.$imem .
            @ iens ++
         end
      else   
         if (! -d $old_dir || ! -e ${old_dir}/${INF_IN_FILE_NAME_PRIOR}) then
            echo We need ${INF_OUT_FILE_NAME} in $old_dir for an adaptive inflation. Stop.
            exit
         else
            while ( $iens <= $nens )
               set imem=${iens}
               if( ${iens} < 1000 ) then
                  set imem=0${iens}
               endif
               if( ${iens} < 100 ) then
                  set imem=00${iens}
               endif
               if( ${iens} < 10 ) then
                  set imem=000${iens}
               endif
               ln -sf ${old_dir}/${RESTART_IN_FILE_NAME}.$imem .
               @ iens ++
            end
            ln -sf ${old_dir}/${INF_IN_FILE_NAME_PRIOR} ./. 
         endif
      endif
#
      if(! -e ${RESTART_IN_FILE_NAME}.0001) then
         echo We need filter_ics to start filter. Stop.
         exit
      endif
#
# Update the observation file for each cycle. 
      set datetmp = `./convertdate | tail -1 |cut -d: -f2` << EOF
      2
      $obs_time_days $obs_time_secs
EOF
#
      set date = `echo $datetmp[1]$datetmp[2]$datetmp[3]$datetmp[4]`
      if(! -e $OBS_DIR/${obs_pre}${date}) then
         echo $OBS_DIR/${obs_pre}${date} for obs_seq.out does not exist. Stop.
         exit
      else
         ln -sf $OBS_DIR/${obs_pre}${date} obs_seq.out
      endif
#
# Job name for filter 
      set JOB_NAME = ${expn}.${icyc}
      cat >! filter.sed << EOF
      s#RUN_DIR#${RUN_DIR}#g
      s#PROJ_NUMBER#${PROJ_NUMBER}#g
      s#NUM_TASKS#${NUM_TASKS}#g
      s#TASKS_PER_NODE#${TASKS_PER_NODE}#g
      s#JOB_NAME#${JOB_NAME}#g
      s#JOB_TIME#${JOB_TIME_FILTER}#g
      s#JOB_CLASS#${JOB_CLASS}#g
EOF
#
      cp ${DART_DIR}/models/wrf_chem/templates/filter.template.lsf .
      sed -f filter.sed filter.template.lsf >! filter.lsf
#
# Submit filter job
      bsub -K < filter.lsf
#
# Check errors in filter.
      if ( ! -e filter_done ) then
         if ( -e filter_started ) then
            if ( -e filter_system_error ) then
               ${REMOVE} filter_system_error
               echo "System failure in filter. We try it again."
#               bsub -K medium filter.lsf
               exit
            else
               echo "Filter was not normally finished. Exiting."
               ${REMOVE} filter_started
               exit
            endif
         endif
      endif
#
      ${REMOVE} filter_started filter_done
      echo Filter is done for cycle ${icyc}: $obs_time_days $obs_time_secs
      ${MOVE} input.nml input.nml.filter.${icyc}
   endif
# END RUN_FILTER CODE BLOCK
#
# BEGIN RUN_ADVANCE_MODEL CODE BLOCK   
   if(${RUN_ADVANCE_MODEL} == true) then
# 
# Update time for advance_model and next cycle
      @ obs_time_secs += $intv_scs
      if ($obs_time_secs >= 86400) then
         @ obs_time_secs -= 86400
         @ obs_time_days += 1
      endif
#
# Convert the analysis to advance model for each member
#      ${REMOVE} assim_model_state_ic.*
      cat >! input.nml.ic.restart_file_tool << EOF
      input_file_name              = ${RESTART_TOOL_IN_FILE_NAME},
      output_file_name             = ${RESTART_TOOL_OUT_FILE_NAME},,
      ens_size                     = $nens,
      single_restart_file_in       = .false.,
      single_restart_file_out      = .false.,
      write_binary_restart_files   = .false.,
      overwrite_data_time          = .false.,
      new_data_days                = -1,
      new_data_secs                = -1,
      input_is_model_advance_file  = .false.,
      output_is_model_advance_file = .true.,
      overwrite_advance_time       = .true.,
      new_advance_days             = $obs_time_days,
      new_advance_secs             = $obs_time_secs 
EOF
      sed -e '/^   input_file_name/,/^   new_advance_secs/d' \
      -e '/&restart_file_tool_nml/r ./input.nml.ic.restart_file_tool' input.nml.filter.${icyc} >! input.nml 
#
# AFAJ: Modify restart tool utility
# Since it takes forever, instead of using restart tool, we just insert the newdate
#
      set iens_filter_ic_new = 1
      rm -f add_newdate
      printf "%6s%11s\n" $obs_time_secs $obs_time_days >! add_newdate
      while ( $iens_filter_ic_new <= $nens ) 
         set ensstring = `echo $iens_filter_ic_new + 10000 | bc | cut -b2-5`
         cat add_newdate filter_ic_new.${ensstring} >! assim_model_state_ic.${ensstring}
         @ iens_filter_ic_new ++
      end
#
      set n_ics = `ls -1 assim_model_state_ic.* | wc -l`
#      if( $n_ics != $nens ) then
#         echo Not enough assim_model_state_ic files: $n_ics. Stop.
#         exit
#      endif
#
# Run advance_model
      echo Advance models for $nens members now...
      ${REMOVE} assim_model_state_ud.*
#
      set n = 1
      while ( $n <= $nens)
         set job_ensemble = GCycle${icyc}_${n}
         cat >! advance.sed << EOF
         s#RUN_DIR#${RUN_DIR}#g
         s#JOB_NAME#${job_ensemble}#g
         s#PROJ_NUMBER#${PROJ_NUMBER}#g
         s#ENS_MEM#${n}#g
         s%WRFINPUT_DIR%${old_dir}%g
         s%ICYCLE%${cycle_restart}%g
         s#NUM_TASKS#${NUM_TASKS}#g
         s#TASKS_PER_NODE#${TASKS_PER_NODE}#g
         s#JOB_TIME#${JOB_TIME_WRFCHEM}#g
         s#JOB_CLASS#${JOB_CLASS}#g
EOF
#
         if (-e advance_model.lsf) $REMOVE advance_model.lsf
         sed -f advance.sed advance_model.template.lsf >! advance_model.lsf
# APM: temp mods
#chmod +x advance_model.lsf
#./advance_model.lsf
#exit
         bsub -K < advance_model.lsf
         cp -f advance_model.lsf $RUN_DIR/advance_temp${n}/
         @ n++
      end
exit
#
# Save model output
      set out_dir = $RUN_DIR/${expn}/Cycle_${icyc}
      if(! -d $out_dir) mkdir -p $out_dir
      echo Saving the output files now...
      ls -lrt >! ${out_dir}/list
      ${COPY} input.nml.filter.${icyc} ${out_dir}/
      ${COPY} input.nml ${out_dir}/input.nml.ic.restart_file_tool.${icyc}
      ${MOVE} wrfinput_d0*_* ${out_dir}/
      ${MOVE} wrf.out* ${out_dir}/
      ${MOVE} dart_log.out ${out_dir}/
      ${MOVE} dart_log.nml ${out_dir}/
#
# We can back up these files in MSS to save some space on bluefire.
# AFAJ deleted filter_ics here since it came from old_dir anyway
# the out_dir shoudl contain the following:
# Prior/Posterior_Diag.nc , obs_seq.final, filter_ic_old corresponding
# to new wrfinput, prior_inf_ic_old corresponding to prior_inf_ic_new
      foreach FILE ( Prior_Diag.nc Posterior_Diag.nc obs_seq.final )
#         echo gzipping $FILE...
#         gzip -f $FILE
         if ( $sav_mss == true ) then
            set md = $mssdir/Cycle_${icyc}
            echo ${MSRCP} $FILE mss:$md/$FILE 
            if($FILE == obs_seq.final) then
               ${MSRCP} $FILE mss:$md/$FILE 
               ${MOVE} obs_seq.final ${out_dir}/
            else
               ${MSRCP} -srcdelete $FILE mss:$md/$FILE 
            endif
            if ( ! $status == 0 ) then
               echo "Failed msrcping $FILE onto $md/"
               echo ${MOVE} $FILE ${out_dir}/ instead.
               ${MOVE} $FILE ${out_dir}/
#               touch BOMBED
            endif
         else
            echo ${MOVE} $FILE ${out_dir}/
            ${MOVE} $FILE ${out_dir}/
            if ( ! $status == 0 ) then
               echo "Failed moving $FILE onto $out_dir/"
               touch BOMBED
            endif
         endif
      end
#
# Save for adaptive inflation for next cycle.
      if( $adaptive_inf == true ) then
         if ( -e $inf_fout && ! -z $inf_fout ) then
#
# need to ask if save_to mss? AFAJ comment for now
#            echo ${MSRCP} $inf_fout mss:$md/$inf_fout
#            ${MSRCP} $inf_fout mss:$md/$inf_fout
# rename inflation AFAJ
           echo ${MOVE}  $inf_fout ${out_dir}/${inf_fin}
           ${MOVE}  $inf_fout ${out_dir}/${inf_fin}
           if ( ! $status == 0 ) then
              echo "Failed moving $inf_fout to ${out_dir}"
              touch BOMBED
           endif
         endif
      endif
# APM unkn endif
#
# A fatal error has occurred.  STOP
      if ( -e BOMBED ) then
         echo "FATAL SYSTEM ERROR"
         touch ABORT_FILTER_RUN
         ${REMOVE} BOMBED
         exit
      endif
#
# For PBS only (AFAJ)
# Check if all members are done advancing model.
      set is_all_done = `qstat | grep GCycle | wc -l`
#      set is_all_done = `ls $RUN_DIR | grep filter_control | wc -l`
      while ( $is_all_done > 0 )
#         echo Wait until wrf runs are completed...
         sleep 20
         set is_all_done = `qstat | grep GCycle | wc -l`
#         set is_all_done = `ls $RUN_DIR | grep filter_control | wc -l`
      end
      echo $is_all_done
#
###########################################################
# Error checking for wrf runs
###########################################################
      cd $RUN_DIR
      set n_prior = `ls -1 assim_model_state_ud.* | wc -l`
      set fn_list = blown_${obs_time_days}_${obs_time_secs}.out
      while ($n_prior != $nens)
         echo There are only $n_prior assim_model_state_ud files.
         set n = 1
         while ( $n <= $nens )
            set iens = `echo $n + 10000 | bc | cut -b2-5`
            if(! -e assim_model_state_ud.${iens}) then
               if( -e Failed_member$n ) then
               echo We already reran this member. Need to check it manually. Stop.
               exit -1
            else
               echo $icyc > Failed_member$n
               cp $RUN_DIR/advance_temp${n}/advance_model.lsf .
               qsub -q short advance_model.lsf
            endif
            echo Member $n failed in cycle $icyc. 
            set logfiles = `ls ./logs/GCycle${icyc}_${n}.*.log`
            set     nlog = `ls ./logs/GCycle${icyc}_${n}.*.log | wc -l`
            set  logfile = $logfiles[$nlog]	# Checking the latest log file.
            set chk_error = `grep lsftmpdir $logfile | wc -l`
#            endif
            @ n++
         end
#
# For PBS only (AFAJ)
         set is_it_done = `qstat | grep GCycle${icyc} | wc -l`
#         set is_all_done = `ls $RUN_DIR | grep filter_control | wc -l`
         while ( $is_it_done > 0 )
            sleep 20
            set is_it_done = `qstat | grep GCycle${icyc} | wc -l`
#            set is_all_done = `ls $RUN_DIR | grep filter_control | wc -l`
         end 
         ls -l $RUN_DIR/advance_temp*/GCycle${icyc}*.log
         mv wrfinput_d01_* $RUN_DIR/${expn}/Cycle_${icyc}/
         set n_prior = `ls -1 assim_model_state_ud.* | wc -l`
      end
      echo Model_advance is successfully done for cycle $icyc. 

###########################################################
# 4. Concatenate priors for all members to run filter
###########################################################
      echo start to clean up
      echo delete filter_ic_old.* 
      if(-e filter_ic_old.0001) ${REMOVE} filter_ic_old.*
      ${COPY} input.nml ${out_dir}/input.nml.ud.restart_file_tool.${icyc}
#     move assim_model_state_ud to filter_ic_old files
      set iens_filter_ic = 1
      while ( $iens_filter_ic <= $nens )
         set ensstring = `echo $iens_filter_ic + 10000 | bc | cut -b2-5`
         ${MOVE} assim_model_state_ud.${ensstring} ${out_dir}/filter_ic_old.${ensstring}
         ${MOVE} filter_ic_new.${ensstring} ${out_dir}/posterior_state_vector.${ensstring}
         @ iens_filter_ic ++
      end
#
# clean up AFAJ
     ${REMOVE} assim_model_state_ic.*
     ${REMOVE} *.sed
  endif
# END RUN_ADVANCE_MODEL CODE BLOCK
#   
  @ icyc++
end
