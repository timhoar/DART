#! /bin/csh -f
#
# This file is one I raided from Andy. He has no knowledge that I copied it,
# and I have no knowledge that it actually works. TJH 15 Feb 2012

set casename = "clm_tim"
set mach     = bluefire
set res      = 1.9x2.5_1.9x2.5
set compset  = ICN

set andydir = /glade/home/afox/cesm1_1/scripts
set timdir = /glade/home/thoar/svn/DART/clm/models/clm/shell_scripts

# TJH ... not for me ...
# cd ${scripdir}
# echo $casename
# rm -rf $casename

cd /glade/user/thoar/cases/

~thoar/cesm1_1_beta08/scripts/create_newcase -case $casename \
                 -mach $mach \
                 -res $res \
                 -compset $compset\
		 -skip_rundb

cd $casename

./xmlchange -file env_conf.xml  -id DATM_MODE         -val CPLHIST3HrWx
./xmlchange -file env_conf.xml  -id DATM_CPL_CASE     -val $casename
./xmlchange -file env_conf.xml  -id DATM_CPL_YR_START -val 2000
./xmlchange -file env_conf.xml  -id DATM_CPL_YR_END   -val 2000
./xmlchange -file env_conf.xml  -id DATM_CPL_YR_ALIGN -val 2000
./xmlchange -file env_conf.xml  -id RUN_STARTDATE     -val 2000-05-01

./xmlchange -file env_mach_pes.xml -id NINST_LND  -val 4
./xmlchange -file env_mach_pes.xml -id NINST_ATM  -val 4
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 64
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 64

./xmlchange -file env_run.xml -id RESUBMIT     -val 0
./xmlchange -file env_run.xml -id STOP_OPTION  -val ndays
./xmlchange -file env_run.xml -id STOP_N       -val 1

./configure -case

cp  ${andydir}/PME_assim/Buildconf/datm.buildnml.csh Buildconf/datm.buildnml.csh.Andy
cp  ${andydir}/PME_assim/Buildconf/clm.buildnml.csh  Buildconf/clm.buildnml.csh

cp  ${timdir}/datm.buildnml.csh                      Buildconf/datm.buildnml.csh

# xxdiff \
# /glade/home/afox/cesm1_1/scripts/PME_assim/SourceMods/src.datm/datm_comp_mod.F90 \
# /glade/home/afox/cesm1_1/models/atm/datm/datm_comp_mod.F90

[thoar@mirage1 SourceMods]$ cd src.datm/
[thoar@mirage1 src.datm]$ vi datm_comp_mod.F90
[thoar@mirage1 src.datm]$ pwd

cp -r ~thoar/cesm1_1_beta08/SourceMods .

./xmlchange -file env_build.xml -id DEBUG  -val TRUE

./*.build

exit

sed '/echo "`date` -- CSM EXECUTION HAS FINISHED"/ r ../call_dart.txt' < $casename.$mach.run >temm
mv temm $casename.$mach.run

sed 's/0:50/1:10/' < $casename.$mach.run >temm
mv temm $casename.$mach.run
