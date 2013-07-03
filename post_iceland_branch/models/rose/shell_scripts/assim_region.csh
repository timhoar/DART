#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# Standard script for use in assimilation applications
# where regions are assimilated by separate executables.
# Can be used with most low-order models, the bgrid model 
# and most large models.

# This script copies the necessary files into the temporary directory
# and then executes the fortran program assim_region.

set  CENTRALDIR = $1
set     element = $2
set    temp_dir = $3

set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'

# Standard script for use in assimilation applications
# where regions are to be assimilated by separate executables.

# Create a clean temporary directory and go there
${REMOVE} $temp_dir
mkdir -p  $temp_dir
cd        $temp_dir

# Copy the initial condition file to the temp directory
${COPY} ${CENTRALDIR}/filter_assim_region__in$element filter_assim_region_in

# Copy the DART namelist to the temp directory
${COPY} ${CENTRALDIR}/input.nml .
${COPY} ${CENTRALDIR}/rose.nml rose.nml

# Copy the assim_region executable to the temporary directory
${COPY} ${CENTRALDIR}/assim_region .

# Copy the observation sequence file to the temporary directory
${COPY} ${CENTRALDIR}/filter_assim_obs_seq .


echo ls $temp_dir for element $element
echo junk >! element$element
ls -lRt 

./assim_region                      >! rose_reg_temp

echo "after executing assim_region"
echo "after executing assim_region" >> rose_reg_temp
echo "element $element"             >> rose_reg_temp
ls -lt                              >> rose_reg_temp
cat rose_reg_temp          >> ${CENTRALDIR}/rose_reg_temp$element

mv filter_assim_region_out ${CENTRALDIR}/filter_assim_region_out$element

# Change back to working directory and get rid of temporary directory
cd ${CENTRALDIR}
#${REMOVE} $temp_dir