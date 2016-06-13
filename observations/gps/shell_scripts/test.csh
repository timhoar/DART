#!/bin/csh
set i = 0
while ($i == 0)
./convert_cosmic_gps_cdf
echo 'continuous? (y:0/n:1)'
set i = $< 
end


