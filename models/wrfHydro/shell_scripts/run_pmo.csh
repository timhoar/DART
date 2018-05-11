#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This is an example script for how to run a perfect_model experiment.
# All the required files have been staged in CENTRALDIR and the namelists
# have been customized for the experiment. Since perfect_model_obs is a single-
# threaded program ... and we are only running one instance of NOAH ... we are
# running this on the command line.

#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Check to make sure all the required files have been staged in CENTRALDIR
#==============================================================================

##./setup_filter.csh forceCopyAll##
@ BAIL = $?

./perfect_model_obs

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
