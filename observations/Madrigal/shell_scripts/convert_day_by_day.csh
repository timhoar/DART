#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Convert all the observation sequence files in the downloads directory.

ls -1 ../downloads/gps* >! obstemp
sort obstemp | uniq >! file_list.txt

# Simply loop over them, converting them one at a time
# The output filenames automatically get appended with the YYYY_MM_DD
# of the observations.

set nfile = `cat file_list.txt | wc -l`
@ ifile = 1
while ( $ifile <= $nfile )

  set fname = `head -n $ifile file_list.txt | tail -n 1`

  sed -e "s#INPUT_FILENAME#${fname}#" input.nml.template >! input.nml

  echo ""
  echo "converting ${fname}"

  ../work/TEC_text_to_obs  || exit 1

  @ ifile += 1
 
end

# if all goes well, clean up
\rm -f obstemp file_list.txt input.nml dart_log.*

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

