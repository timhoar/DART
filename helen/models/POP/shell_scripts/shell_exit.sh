#!/bin/sh

# LSF does not reliably return an exit code from
# a serial section of a batch script, only the return
# from the last call to mpirun.lsf.  so to make a
# batch script exit, call this:
# 
#  setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
#  setenv EXITCODE -1
#  mpirun.lsf shell_exit.sh
#

if [[ $# -gt 0 ]]; then
   EXITCODE=$1
fi

if [[ "$EXITCODE" == "" ]]; then
   EXITCODE=0
fi

echo exiting with status code $EXITCODE
exit $EXITCODE
   
