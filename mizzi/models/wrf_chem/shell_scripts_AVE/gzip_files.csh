#!/bin/csh

set ncyc = 30
set icyc = 1
set rootdir = /s1/arellano/DART/rundir_gabi/June08_Gabi

while ( $icyc <= $ncyc ) 
  cd ${rootdir}/Cycle_${icyc}
  gzip filter*
  @ icyc ++
end
