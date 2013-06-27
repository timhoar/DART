#!/bin/csh
# generate multiple days of gps observations

# args are date, working directory location, and whether to download
# data automatically from the cosmic web site (downloading data
# requires signing up for a username/password to access the site,
# and then setting env vars before running the script.)

setenv cosmic_user xxx
setenv cosmic_pw   yyy

./cosmic_to_obsseq.csh 20061101 .. yes
./cosmic_to_obsseq.csh 20061102 .. yes
./cosmic_to_obsseq.csh 20061103 .. yes
./cosmic_to_obsseq.csh 20061104 .. yes
./cosmic_to_obsseq.csh 20061105 .. yes
./cosmic_to_obsseq.csh 20061106 .. yes
./cosmic_to_obsseq.csh 20061107 .. yes
