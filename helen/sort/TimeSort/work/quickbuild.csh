#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: quickbuild.csh 5074 2011-07-15 17:06:58Z thoar $
#
# This script compiles all executables in this directory.


set MODEL = "timing"

@ n = 1

foreach TARGET ( mkmf_* ) 

  set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

  echo "${MODEL} build number ${n} is preprocess"

  csh $TARGET || exit $n
  make        || $n

end

rm -f *.o *.mod

exit 0
# <next few lines under version control, do not edit>
# $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/development/models/wrf/work/quickbuild.csh $
# $Revision: 5074 $
# $Date: 2011-07-15 11:06:58 -0600 (Fri, 15 Jul 2011) $

