#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: start_over.csh 5074 2011-07-15 17:06:58Z thoar $

\rm -r advance_temp*
\rm assim_model_state*
\rm filter_control*
\rm dart_log*
\rm filter.*.log
\rm Prior_Diag.nc
\rm Posterior_Diag.nc
\rm blown*.out
\rm obs_seq.final
\rm finished_grid*
\rm last_p*
\rm refl_obs*.txt
\rm finished.g*
\rm WRFOUT/*

exit $status

# <next few lines under version control, do not edit>
# $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/experiments/Radar/start_over.csh $
# $Revision: 5074 $
# $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $
