#!/bin/ksh

#-- This is just a wrapper to launch num_states instances
#-- of model advance, here filter is serial but launches
#-- num_states number of model advancement process concurrently

process=$1
num_states=$2
control_file=$3

rm -f assim_model_state_ud.*

if [ ${process} == "0" ] ; then
   source ./env.lsf
   mpirun.lsf ./adv_model_par.csh ${num_states} ${control_file}
fi

\rm -rf $control_file

exit 0
