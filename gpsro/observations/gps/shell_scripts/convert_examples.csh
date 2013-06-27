#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# example script for special purpose use.  see the 'do_convert.csh'
# script for one that loops over any number of days, including rolling
# over month and year boundaries correctly.
#
#
# calls the gpsro_to_obsseq script with 6 args:
#
#  - the date in YYYYMMDD format
#  - the processing directory location, relative to the 'work' dir.
#  - whether to download the data automatically from the cosmic web site.
#     'yes' will do the download, 'no' assumes the data is already downloaded
#     and on local disk.  (downloading data requires signing up for a 
#     username/password to access the site, and then setting the username 
#     and password in the cosmic_to_obsseq script before running it.)
#  - whether to convert the data.  set to 'yes' to make obs_seq files (the
#     usual use of this script). 'no' if just downloading or just cleaning up.
#  - whether to delete the data automatically from the local disk after the
#     conversion is done.  valid values are 'yes' or 'no'.
#  - a text file containing the list of satellites to convert data from.
#     valid strings, one per line, are currently:  
#       cosmic sacc ncofs grace tsx metopa champ 
#     alternatives for the realtime obs for these three satellites:
#       cosmicrt saccrt ncofsrt
#

# examples of common use follow.  

# download only:
./gpsro_to_obsseq.csh 20071001 ../gpsro yes no no satlist
./gpsro_to_obsseq.csh 20071002 ../gpsro yes no no satlist

# convert only.  assume all data already downloaded:
./gpsro_to_obsseq.csh 20071001 ../gpsro no yes no satlist
./gpsro_to_obsseq.csh 20071002 ../gpsro no yes no satlist

# download and convert, not removing files:
./gpsro_to_obsseq.csh 20071001 ../gpsro yes yes no satlist
./gpsro_to_obsseq.csh 20071002 ../gpsro yes yes no satlist

# clean up only after verifying conversion worked:
./gpsro_to_obsseq.csh 20071001 ../gpsro no no yes satlist
./gpsro_to_obsseq.csh 20071002 ../gpsro no no yes satlist

# download, convert, and clean up all in one go:
./gpsro_to_obsseq.csh 20071001 ../gpsro yes yes yes satlist
./gpsro_to_obsseq.csh 20071002 ../gpsro yes yes yes satlist

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

