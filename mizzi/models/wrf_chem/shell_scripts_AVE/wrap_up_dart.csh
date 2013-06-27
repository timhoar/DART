#!bin/csh

set REMOVE = '/bin/rm -rf'
set   COPY = 'cp -pf'
set   MOVE = 'mv -f'

set expn = June08_Gabi          # experiment name - match with $mssdir in advance_model.csh
set nens = 40                   # ens_size
set ndom = 1  
set RUN_DIR = `pwd`

set icyc = 38
set out_dir = $RUN_DIR/${expn}/Cycle_${icyc}
###########################################################
  # 4. Concatenate priors for all members to run filter
  ###########################################################
  echo start to clean up
  echo delete filter_ic_old.*
  if(-e filter_ic_old.0001) ${REMOVE} filter_ic_old.*

  ${COPY} input.nml ${out_dir}/input.nml.ud.restart_file_tool.${icyc}
  #  move assim_model_state_ud to filter_ic_old files
     set iens_filter_ic = 1
     while ( $iens_filter_ic <= $nens )
        set ensstring = `echo $iens_filter_ic + 10000 | bc | cut -b2-5`
        ${MOVE} assim_model_state_ud.${ensstring} ${out_dir}/filter_ic_old.${ensstring}
        ${MOVE} filter_ic_new.${ensstring} ${out_dir}/posterior_state_vector.${ensstring}
        @ iens_filter_ic ++
     end

  # clean up AFAJ
  ${REMOVE} assim_model_state_ic.*
  ${REMOVE} *.sed


