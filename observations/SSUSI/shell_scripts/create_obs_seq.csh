#!/bin/csh
#
# A script to produce multiple obs_seqs from the SSUSI files and combine them together 
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$


# The ON2_UNCERTAINTY variable in the netcdf files have IEEE NaN values,
# but none of the required metadata to interpret them correctly.
# These 2 lines will add the required attributes so that NaNs are replaced with
# a fill value that can be queried and checked for.
# Since the ON2_UNCERTAINTY is a standard deviation, it is enough to make it negative

module load nco

# The intent is to convert all the files in the data directory, one-by-one
# As we do that, we can create a list of those files and feed that list
# into obs_sequence_tool to consolidate them into a single observation sequence file.

\rm file_list.txt

# ensure the file names are as expected

cp ../work/input.nml input.nml

sed -i "s#input_netcdf_file .*#input_netcdf_file = 'netcdf_input.nc'#" input.nml
sed -i "s#output_obs_file .*#output_obs_file = 'obs_seq.out'#" input.nml
sed -i "s#filename_seq_list .*#filename_seq_list = 'file_list.txt'#" input.nml
sed -i "s#filename_seq .*#filename_seq = ''#" input.nml

set i = 1
foreach fname (`ls ../data/*.NC`)  

    # ensure the netCDF metadata is correct

    ncatted -a _FillValue,ON2_UNCERTAINTY,o,f,NaN    $fname
    ncatted -a _FillValue,ON2_UNCERTAINTY,m,f,-1.0   $fname

    ln -sf $fname netcdf_input.nc
  
    ../work/convert_f16_edr_dsk || exit 1

    # create a nice sequential output file name
    set ofname = `printf obs_seq.out_%04d $i`
   
    mv obs_seq.out $ofname || exit 2

    echo $ofname >> file_list.txt
   
    endif
    @ i++
end   

../work/obs_sequence_tool || exit 3

echo "Wrote to obs_seq.processed"

\rm -f file_list.txt dart_log.*

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

