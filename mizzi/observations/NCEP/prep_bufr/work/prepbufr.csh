#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
#**************************************************************
# MODIFIED BY Arthur Mizzi SO daily=no PROCESSES ONLY SINGLE TIME LEVEL
#**************************************************************
#
#--------------------------------------------------------------
# DESCRIPTION:
#
#  This script is used to generate daily (3:01Z to 3:00Z of next day) decoded 
#  NCEP reanalysis PREPBUFR text/ascii data.
#
# there are two ways to run this script - either submit it to a batch
# system (LSF job commands are included), or this script can be called
# as a command line executable.  see the sections below for what options
# to set at the top of the script before it is run.
#
# LSF batch system settings, which will certainly need to be customized
# for your system, e.g. -P project charge code, different queue name, etc.
#BSUB -o prepbufr.out
#BSUB -e prepbufr.err
#BSUB -J prepbufr
#BSUB -q regular
#BSUB -W 1:00
#BSUB -P XXXXXXXX
#BSUB -n 1
# 
# to run this from the command line, see below for a choice of whether
# to invoke this script with 4 args: year, month, start/end day, 
# or hardcode the dates into the script and run it by name with no args.
#
#--------------------------------------------------------------
# USER SET PARAMETERS

# if daily is 'yes', 4 6-hour files will be processed and a single, 1-day
# output file will be created.  if daily is set to anything else, then
# each 6-hour input file will be converted to a single 6-hour output file.
# this script still processes a day's worth of files at a time even if
# daily is 'no' - it just makes 4 individual files per day.
#
# SCRIPT MODIFIED SO daily=no PROCESSES ONE TIME LEVEL ONLY

set    daily = no

# if daily is 'no' and zeroZ is 'yes', input files at 0Z will be translated 
# into output files also marked 0Z.  otherwise, they will be named with the 
# previous day number and 24Z (chose here based on what script will be 
# processing these files next.  the 'create_real_obs' script assumes 
# filenames with the pattern 6Z,12Z,18Z,24Z, so 'no' is right for it.)
# this variable is ignored completely if daily is 'yes'.

set zeroZ = no

# if convert is 'yes', then the big-endian BUFR files will be converted
# to little-endian files before processing. this is needed if you are running
# on a machine that uses Intel chips (e.g. linux clusters, altix, pcs, etc).
# it is not needed for ibm power systems.  any value other than 'yes' will
# skip the convert step.

set  convert = yes 

# if block is 'yes', then the cword program will be run to convert an
# unblocked file into a blocked one.  this is not required for recent
# prepbufr files, but older ones may require it.

set block = no

# starting year, month, day, and ending day.  this script does not allow
# you to do more than a single month at a time, but does handle the last
# day of the month, leap day in feb, and the last day of the year correctly.
# this version of the conversion tool takes up to 3 hours of observations
# from the day *following* the end day, so you must have at least the 6Z
# file from one day beyond the last day for this to finish ok (to get obs
# at exactly 3Z from the next file).

# if this is being called by another program, and passing in the dates
# as command line arguments (or if you just prefer to call this with args), 
# then set commline to 'yes'.  otherwise, set the year, month, and days 
# directly in the variables below.

set commline = no

if ($commline == 'yes') then
  set year     = $argv[1]
  set month    = $argv[2]
  set beginday = $argv[3]
  set endday   = $argv[4]
else
  set year     = ${DT_YYYY}
  set month    = ${DT_MM}
  set beginday = ${DT_DD}
  set endday   = ${DT_DD}
endif

# directory where the BUFR files are located.  the script assumes the
# files will be located in subdirectories by month, with the names following
# the pattern YYYYMM, and then inside the subdirectories, the files are
# named by the pattern 'prepqmYYMMDDHH'.  for example, if the dir below
# is the default ../data, then the 6Z file for nov 27th, 2010 would be:
#  ../data/201011/prepqm10112706
# but the conventions for names of prepqm files have changed over the years,
# so if the prepqm files do *not* follow this pattern, you will have to edit
# the BUFR_in variable in the script below to match the filenames you have.
# the setting of BUFR_out matches what the 'create_real_obs' script expects
# to have as input, but if you want to generate a different set of names
# you can change it below as well.
# there are several shell variables in the loop you can use to construct
# alternate names:
# 'year' is 4 digits; 'yy' is 2.
# 'mm', 'dd', and 'hh' are 0 padded so they are always 2 digits.
# 'oyear', 'omm', 'odd', 'ohh' are the original date, if the day, month, 
# and/or year have rolled over.

#set BUFR_dir = ../data
set BUFR_dir = ${ASIM_DIR}/prep_bufr

# directory where DART prepbufr programs are located, relative
# to the directory where this script is going to be executed.
set DART_exec_dir = ${DART_DIR}/observations/NCEP/prep_bufr/exe

# END USER SET PARAMETERS
#--------------------------------------------------------------

      # the prepqm input files.  match general naming pattern with the
      # data time encoded in the filename.  if the pattern of the filename
      # is different (2 digit year vs 4, extra fixed text in the name, etc)
      # fix the BUFR_in line below to match what you have.  if the file is
      # gzipped, you can leave it and this program will unzip it before
      # processing it.

      set BUFR_in = ${BUFR_dir}/${YYYY}${MM}/prepqm${DT_YY}${MM}${DD}${HH}

      if ( -e ${BUFR_in} ) then
         echo "copying ${BUFR_in} into prepqm.in"
         rm -f prepqm.in
         cp -f ${BUFR_in} prepqm.in
      else if ( -e ${BUFR_in}.gz ) then
         echo "unzipping ${BUFR_in}.gz into prepqm.in"
         rm -f prepqm.in
         gunzip -c -f ${BUFR_in}.gz >! prepqm.in
      else
         echo "MISSING INPUT FILE: cannot find either"
         echo ${BUFR_in}
         echo   or 
         echo ${BUFR_in}.gz
         echo "Script will abort now."
         exit -1
      endif

      # blocking
      if ($block == 'yes') then
         echo "blocking prepqm.in"
         mv -f prepqm.in prepqm.unblocked
         echo 'block' >! in
         echo 'prepqm.unblocked' >> in
         echo 'prepqm.blocked' >> in
         ${DART_exec_dir}/cword.x < in
         mv -f prepqm.blocked prepqm.in
         rm -f prepqm.unblocked in
      endif

      # byte swapping
      if ($convert == 'yes') then
         echo "byteswapping bigendian to littleendian prepqm.in"
         mv -f prepqm.in prepqm.bigendian
         ${DART_exec_dir}/grabbufr.x prepqm.bigendian prepqm.littleendian
         mv -f prepqm.littleendian prepqm.in
         rm -f prepqm.bigendian
      endif

      ${DART_exec_dir}/prepbufr.x
      set BUFR_out = ${BUFR_dir}/${YYYY}${MM}/temp_obs.${YYYY}${MM}${DD}${HH}
      if (${DT_HH} == 0) then
         if ($zeroZ == 'no') then
            # if 0Z, output named with previous day and 24Z
            set l_date=${PAST_DATE}
            set l_yyyy=${PAST_YYYY}
            set l_mm=${PAST_MM}
            set l_dd=${PAST_DD}
            set l_hh=24
            set BUFR_out = ${BUFR_dir}/${YYYY}${MM}/temp_obs.${l_yyyy}${l_mm}${l_dd}${l_hh}
         endif
      endif
      echo "moving output to ${BUFR_out}"
      mv -fv prepqm.out ${BUFR_out}

   rm -f prepqm.in

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

