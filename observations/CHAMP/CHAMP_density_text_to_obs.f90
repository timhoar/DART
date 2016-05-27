! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program CHAMP_density_text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CHAMP_density_text_to_obs - reads fixed-format ASCII files from
!   http://sisko.colorado.edu/sutton/data/ver2.2/champ/density/2002/ascii/
!   work/Density_3deg_02_335.ascii is an example of an input file
!
!   created  29 Mar 2010 Nancy Collins NCAR/IMAGe
!   modified 15 Aug 2012 Alexey Morozov (Univ. of Michigan)
!   modified 25 May 2016 Tim Hoar NCAR/IMAGe
!
!   Since the GRACE and CHAMP data available from Erik Sutton are in the
!   same format, this converter now has namelist options to specify the
!   observation type of the output. Appending to an existing file 
!   (inserting, actually) is also namelist-controlled.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, PI, DEG2RAD, obstypelength

use     utilities_mod, only : initialize_utilities, finalize_utilities, to_upper, &
                              open_file, close_file, find_namelist_in_file, &
                              check_namelist_read, error_handler, E_MSG, E_ERR

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_date, &
                              operator(-), operator(+), operator(>=)

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use      obs_kind_mod, only : CHAMP_NEUTRAL_DENSITY, GRACEA_NEUTRAL_DENSITY, GRACEB_NEUTRAL_DENSITY, SAT_RHO

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! things which can/should be in the text_to_obs_nml

character(len=256) :: text_input_file            = 'Density_3deg_05_020.ascii'
character(len=256) :: obs_out_file               = 'obs_seq.out'
character(len=obstypelength) :: observation_type = 'CHAMP_NEUTRAL_DENSITY'
logical            :: append_to_existing_file    = .false.
logical            :: debug                      = .true.

namelist /CHAMP_density_text_to_obs_nml/  &
     text_input_file, &
     obs_out_file,    &
     observation_type,  &
     append_to_existing_file,  &
     debug

character(len=512)  :: string1, string2
character(len=1024) :: input_line !162 is the nominal width of CHAMPdens2.2 data records,
                                  ! but the second line, the descriptor, has 731 characters.

integer :: oday, osec, rcio, iunit
integer :: year, day, second
integer :: num_copies, num_qc, max_obs, linenum
integer :: observation_type_int

logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc
real(r8) :: lat, lon, vert

real(r8) :: second_r !CHAMP seconds are reals instead of ints

!variables to be discarded (only needed so that the read line works)
integer  :: ignore_i
real     :: ignore_r

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! start of executable code

call initialize_utilities('CHAMP_density_text_to_obs')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'CHAMP_density_text_to_obs_nml', iunit)
read(iunit, nml = CHAMP_density_text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'CHAMP_density_text_to_obs_nml')

call set_observation_type()

! time setup
call set_calendar_type(GREGORIAN)

! open input text file
iunit = open_file(text_input_file, 'formatted', 'read')
write(string1,*) 'opened input file "' // trim(text_input_file) // '"'
write(string2,*) 'converting them as observation type '//trim(observation_type)

call error_handler(E_MSG, 'CHAMP_density_text_to_obs', string1, &
           source, revision, revdate, text2=string2)

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 100000
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! If you want to append to existing files (e.g. you have a lot of
! small text files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files
! once they are in DART obs_seq format.

if (append_to_existing_file) then
   inquire(file=obs_out_file, exist=file_exist)
   if ( file_exist ) then
     write(string1,*)'..  inserting into "'//trim(obs_out_file)//'"'
     call error_handler(E_MSG,'CHAMP_density_text_to_obs',string1)
     call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
   endif
endif

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

! The first  line is the version and origin information.
! The second line is a description of the columns and units.
! As long as these are constant, we can skip them.
! column 01 * Two-digit Year (years)
! column 02 * Day of the Year (days)
! column 03 * Second of the Day (GPS time,sec)
! column 04 * Center Latitude of 3-degree Bin (deg)
! column 05 * Satellite Geodetic Latitude (deg)
! column 06 * Satellite Longitude (deg)
! column 07 * Satellite Height (km)
! column 08 * Satellite Local Time (hours)
! column 09 * Satellite Quasi-Dipole Latitude (deg)
! column 10 * Satellite Magnetic Longitude (deg)
! column 11 * Satellite Magnetic Local Time (hours)
! column 12 * Neutral Density (kg/m^3)
! column 13 * Neutral Density Normalized to 400km using NRLMSISe00
! column 14 * Neutral Density Normalized to 410km using NRLMSISe00
! column 15 * NRLMSISe00 Neutral Density at Satellite Height
! column 15 * Uncertainty in Neutral Density (kg/m^3)
! column 17 * Number of Data Points in Current Averaging Bin
! column 18 * Number of Points in Current Averaging Bin that Required Interpolation
! column 19 * Average Coefficient of Drag Used in Current Averaging Bin

read(iunit,"(A)") input_line
read(iunit,"(A)") input_line

linenum = 2

obsloop: do    ! no end limit - have the loop break when input ends

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! read the whole line into a buffer and parse it later
   read(iunit, "(A)", iostat=rcio) input_line

   if (rcio < 0) then
      write(string1,*) trim(text_input_file)//' had ', linenum-2,' observations.'
      call error_handler(E_MSG,'CHAMP_density_text_to_obs',string1, &
                         source, revision, revdate)
      exit obsloop
   endif

   if (rcio /= 0) then
      write(string1,*) 'got bad read code (', rcio,') on line ',linenum
      write(string2,*) 'of ',trim(text_input_file)
      call error_handler(E_ERR,'CHAMP_density_text_to_obs',string1, &
                   source, revision, revdate, text2=string2)
   endif

   linenum = linenum + 1

   ! here is a line from sisko.colorado.edu/sutton/data/ver2.2/champ/density/2002/ascii/,
   !data format is:
   !+ 1)year(2I), 2)day(3I), 3)second(8.3F), 4)round(lat), 5)lat(d,-90 90), 6)lon(d,-180 180), 7)alt(km),
   !+ 8)LT, 9)Mlat, 10)Mlon, 11)MLT, 12)Rho(Density!), 13)MSISRho400, 14)MSISRho410, 15)MSISRhoSat
   !+ 16)Rho uncertainty (I'm guessing std deviation from units: kg/m^3) 17)points averaged over
   !+ 18)points needing interpolation 19)coeff of drag averaged over bin

   read(input_line, *, iostat=rcio) &
        year, day, second_r, ignore_i, lat, lon, vert, &
        ignore_r, ignore_r, ignore_r, ignore_r, temp, ignore_r, ignore_r, ignore_r, &
        terr, ignore_i, &
        ignore_i, ignore_r

   if (rcio /= 0) then
      write(string1,*) 'unable to parse line ',linenum
      write(string2,*) 'of ',trim(text_input_file)
      call error_handler(E_ERR,'CHAMP_density_text_to_obs',string1, &
                      source, revision, revdate, text2=string2)
   endif

   vert = vert * 1000.0_r8 ! DART needs alt in meters, CHAMP has km

   if (debug) print *, 'this observation located at lat, lon, vert = ', lat, lon, vert

   ! if lon comes in between -180 and 180, use these lines instead:
   if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   if ( lon <   0.0_r8 ) lon = lon + 360.0_r8 ! changes into 0-360

   ! put date into a dart time format

   year      = 2000 + year !because year in file is (2I) - 2 digits
   comp_day0 = set_date(year, 1, 1, 0, 0, 0)  ! always Jan 1 of whatever year.
   second    = nint(second_r)
   time_obs  = comp_day0 + set_time(second, day-1)

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   if (debug) call print_date(time_obs, 'this obs time is')

   ! make an obs derived type, and then add it to the sequence

   call create_3d_obs(lat, lon, vert, VERTISHEIGHT, temp, &
                      observation_type_int, terr, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   if (debug) print *, 'added '//trim(observation_type)//' obs to output seq'

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   !if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
else
   call error_handler(E_MSG,'CHAMP_density_text_to_obs','no observations in sequence', &
                      source, revision, revdate)
endif

! end of main program
call finalize_utilities()

contains

!-----------------------------------------------------------------------
!> 

subroutine set_observation_type

! sets the global variable 'observation_type_int' based on the character
! string namelist input 'observation_type'

character(len=obstypelength) :: observation_string

! must create local copy because to_upper works in-place
observation_string = observation_type
call to_upper(observation_string)

if     (trim(observation_string) ==  'CHAMP_NEUTRAL_DENSITY') then
              observation_type_int =  CHAMP_NEUTRAL_DENSITY
elseif (trim(observation_string) == 'GRACEA_NEUTRAL_DENSITY') then
              observation_type_int = GRACEA_NEUTRAL_DENSITY
elseif (trim(observation_string) == 'GRACEB_NEUTRAL_DENSITY') then
              observation_type_int = GRACEB_NEUTRAL_DENSITY
elseif (trim(observation_string) == 'SAT_RHO') then
              observation_type_int = SAT_RHO
else 
   write(string1,*)'Unable to interpret observation string "'//trim(observation_type)//'"'
   write(string2,*)'valid strings are "CHAMP_NEUTRAL_DENSITY", "GRACEA_NEUTRAL_DENSITY", "GRACEB_NEUTRAL_DENSITY", "SAT_RHO"'
   call error_handler(E_ERR, 'CHAMP_density_text_to_obs', string1, &
                      source, revision, revdate, text2=string2)
endif  

end subroutine set_observation_type



end program CHAMP_density_text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
