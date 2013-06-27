#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: compile_pert_sounding.csh 5074 2011-07-15 17:06:58Z thoar $


# Compiles pert_sounding.f90 for ocotillo

#ifort -c pert_sounding_module.f90
#ifort pert_sounding.f90 -o pert_sounding pert_sounding_module.o
#\rm ./pert_sounding_mod.mod



# Compiles pert_sounding.f90 for bluefire

xlf -c pert_sounding_module.f90
xlf pert_sounding.f90 -o pert_sounding pert_sounding_module.o
rm ./pert_sounding_mod.mod

exit $status

# <next few lines under version control, do not edit>
# $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/experiments/Radar/IC/sounding_perturbation/compile_pert_sounding.csh $
# $Revision: 5074 $
# $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $

