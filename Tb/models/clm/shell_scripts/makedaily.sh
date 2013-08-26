#!/bin/bash
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# split a yearly file into "daily" files which start at 00:00Z
# the previous day and end at 23:59Z on the day that matches the
# day in the filename.

# set the first and last days to be split.
# depending on the window and the input file, 
# the data from outside these bounds may be needed.

let start_year=2001
let start_month=01
let start_day=01

let end_year=2001
let end_month=1
let end_day=2

# end of things you should have to set in this script if you are
# content to have 'daily' files with observations +/- 12 hours from
# date in the filename.

# convert the start and stop times to gregorian days, so we can compute
# total number of days including rolling over month and year boundaries.
# do the end time first so we can use the same values to set the 
# initial day while we are doing the total day calculation.

# make sure there is an initial input.nml for advance_time
# input.nml gets overwritten in the subsequent loop.
cp -f  ../work/input.nml.template input.nml || exit -1
ln -sf ../work/clm_history.nc             . || exit -2
ln -sf ../work/clm_restart.nc             . || exit -3

# advance_time (with the -g flag) outputs 2 integers:
# days seconds
year1=`echo  $start_year  | bc`
month1=`echo $start_month | bc`
day1=`echo   $start_day   | bc`
dart_1=`printf %04d%02d%02d%02d $year1 $month1 $day1 0`

year2=`echo  $end_year  | bc`
month2=`echo $end_month | bc`
day2=`echo   $end_day   | bc`
dart_N=`printf %04d%02d%02d%02d $year2 $month2 $day2 0`

# these outputs from advance time (with the -g flag) are
# 2 integers: gregorian_day_number seconds
# and since we don't set hours, minutes, or seconds, the second
# number is always 0 and uninteresting for us.
mon2=`printf %02d $end_month`
day2=`printf %02d $end_day`
end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "last  day,seconds is day ${end_d[0]} ${end_d[1]}"

mon2=`printf %02d $start_month`
day2=`printf %02d $start_day`
start_d=(`echo ${start_year}${mon2}${day2}00 0 -g | ../work/advance_time`)
echo "first day,seconds is day ${start_d[0]} ${start_d[1]}"

# how many total days are going to be split (for the loop counter)
# (pull out the first of the 2 numbers which are output from advance_time)
let totaldays=${end_d[0]}-${start_d[0]}+1

# form some strings for logging.
# time_one    .... the first time in the file
# time_end    .... the last  time in the file
# filetime    .... the time in the file NAME ... usually the center.
# with no -g option advance_time returns strings in the format YYYYMMDDHH

time_one=(`echo ${start_year}${mon2}${day2}00 -12h | ../work/advance_time`)
filetime=(`echo ${start_year}${mon2}${day2}00    0 | ../work/advance_time`)
time_end=(`echo ${start_year}${mon2}${day2}00 +12h | ../work/advance_time`)

# loop over each day
let d=1
while (( d <= totaldays)) ; do

  echo "subsetting $d of $totaldays ..."
  #echo $filetime $time_end

  # string for first time in the file
   pyear=${time_one:0:4}
  pmonth=${time_one:4:2}
    pday=${time_one:6:2}
   phour=${time_one:8:2}

  # string for time in the file NAME
   cyear=${filetime:0:4}
  cmonth=${filetime:4:2}
    cday=${filetime:6:2}
   chour=${filetime:8:2}

  # string for last time in the file
   nyear=${time_end:0:4}
  nmonth=${time_end:4:2}
    nday=${time_end:6:2}
   nhour=${time_end:8:2}

  # compute the equivalent DART timestamps here - seconds and days.
  g=(`echo ${cyear}${cmonth}${cday}${chour} -12h -g | ../work/advance_time`)
  dart0d=${g[0]}
  darts0=${g[1]}
  let dart0s=${darts0}+1

  g=(`echo ${cyear}${cmonth}${cday}${chour}    0 -g | ../work/advance_time`)
  dart1d=${g[0]}
  darts1=${g[1]}
  let dart1s=${darts1}

  g=(`echo ${cyear}${cmonth}${cday}${chour} +12h -g | ../work/advance_time`)
  dart2d=${g[0]}
  darts2=${g[1]}
  let dart2s=${darts2}

  echo prev $pyear $pmonth $pday $phour which is dart $dart0d $dart0s
  echo curr $cyear $cmonth $cday $chour which is dart $dart1d $dart1s
  echo next $nyear $nmonth $nday $nhour which is dart $dart2d $dart2s

  # I have annual files  ...
  # I'll need to revisit this when I wrap over year boundaries ... TJH

  sed -e "s/YYYY/${cyear}/g"    \
      -e "s/MM/${cmonth}/g"     \
      -e "s/DD/${cday}/g"       \
      -e "s/SSSSS/${chour}/g"   \
      -e "s/DART0D/${dart0d}/g" \
      -e "s/DART0S/${dart0s}/g" \
      -e "s/DART2D/${dart2d}/g" \
      -e "s/DART2S/${dart2s}/g" < ../work/input.nml.template > input.nml

  # make sure output dir exists
  if [[ ! -d ../${cyear}${cmonth} ]] ; then
     mkdir ../${cyear}${cmonth}
  fi

  # do the extract here
  ../work/obs_sequence_tool

  # advance the day; the output is YYYYMMDD00
  time_one=(`echo ${pyear}${pmonth}${pday}${phour} +1d | ../work/advance_time`)
  filetime=(`echo ${cyear}${cmonth}${cday}${chour} +1d | ../work/advance_time`)
  time_end=(`echo ${nyear}${nmonth}${nday}${nhour} +1d | ../work/advance_time`)
  echo "next set of times are: $time_one $filetime $time_end"

  # advance the loop counter
  let d=d+1

done

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

