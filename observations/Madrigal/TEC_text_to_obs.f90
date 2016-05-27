! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use a
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program TEC_text_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   TEC_text_to_obs - reads fixed-format ASCII files from
!   httpsiskocoloradoedu/sutton/data/ver2.2/champ/density/2002/ascii/
!   work/Density_3deg_02_335.ascii is an example of an input file
!
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

use      obs_kind_mod, only : GPS_VTEC_EXTRAP

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! namelist variables

character(len=256) :: text_input_file = 'gps050122g.002.txt'
character(len=256) :: obs_out_base    = 'obs_seq.out'

namelist /TEC_text_to_obs_nml/  &
     text_input_file, &
     obs_out_base

! everbody else

character(len=obstypelength) :: observation_type = 'GPS_VTEC_EXTRAP'

character(len=256)  :: obs_out_file
character(len=512)  :: string1, string2
character(len=1024) :: input_line

integer :: oday, osec, rcio, iunit
integer :: iloc, num_missing
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs, linenum
integer :: observation_type_int

logical  :: first_obs

real(r8) :: tec, tec_std, qc
real(r8) :: lat, lon

!The error variance was chosen by Alex Chartier ... this includes
! observation error and error of representativeness.
real(r8), PARAMETER :: observation_error_variance = 25.0_r8

!variables to be discarded (only needed so that the read line works)
integer  :: ignore_i

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

! start of executable code

call initialize_utilities('TEC_text_to_obs')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'TEC_text_to_obs_nml', iunit)
read(iunit, nml = TEC_text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'TEC_text_to_obs_nml')

! Not really needed, but already here when imported from CHAMP ...
call set_observation_type()

! time setup
call set_calendar_type(GREGORIAN)

! open input text file
iunit = open_file(text_input_file, 'formatted', 'read')
write(string1,*) 'opened input file "' // trim(text_input_file) // '"'
write(string2,*) 'converting them as observation type '//trim(observation_type)

call error_handler(E_MSG, 'TEC_text_to_obs', string1, &
           source, revision, revdate, text2=string2)

! each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 2000000
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

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.

qc = 0.0_r8

! The first line is a description of the columns and units.
! As long as these are constant, we can skip them.
! column 01 * YEAR
! column 02 * MONTH
! column 03 * DAY
! column 04 * HOUR
! column 05 * MIN
! column 06 * SEC
! column 07 * UT1_UNIX   (skipping)
! column 08 * UT2_UNIX   (skipping)
! column 09 * RECNO      (skipping)
! column 10 * GDLAT Geodetic latitude of measurement - Units: deg
! column 11 * GLON Geographic longitude of measurement - Units: deg
! column 12 * TEC Vertically integrated electron density - Units: tec
! column 13 * DTEC Error in Vertically integrated electron density - Units: tec
!
! 1 TEC unit = 1E16 electrons per square meter

read(iunit,"(A)") input_line

linenum = 1
num_missing = 0

obsloop: do    ! no end limit - have the loop break when input ends

   if (mod(linenum,100000) == 0) write(*,*)'Processing line ',linenum 

   ! read in a line from the text file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)

   ! read the whole line into a buffer and parse it later
   read(iunit, "(A)", iostat=rcio) input_line

   if (rcio < 0) then
      write(string1,*) trim(text_input_file)//' had ', linenum-1,' data lines.'
      write(string2,*) 'of those, ',num_missing,' had "missing" values.'
      call error_handler(E_MSG,'TEC_text_to_obs',string1, &
                 source, revision, revdate, text2=string2)
      exit obsloop
   endif

   if (rcio /= 0) then
      write(string1,*) 'got bad read code (', rcio,') on line ',linenum
      write(string2,*) 'of ',trim(text_input_file)
      call error_handler(E_ERR,'TEC_text_to_obs',string1, &
                 source, revision, revdate, text2=string2)
   endif

   linenum = linenum + 1

   ! Some of the lines have 'missing' (as a character string). If this 
   ! is present - count them up, skip the line and cycle to the next line.

   iloc = index(input_line,'missing')
   if (iloc /= 0) then
      num_missing = num_missing + 1
      cycle obsloop
   endif

   ! assumed format:
   !  1    2     3   4    5   6     7        8       9    10    11   12  13
   ! YEAR MONTH DAY HOUR MIN SEC UT1_UNIX UT2_UNIX RECNO GDLAT GLON TEC DTEC

   read(input_line, *, iostat=rcio) &
        year, month, day, hour, minute, second, ignore_i, ignore_i, ignore_i, lat, lon, tec, tec_std

   if (rcio /= 0) then
      write(string1,*) 'unable to parse line ',linenum
      write(string2,*) 'of ',trim(text_input_file)
      call error_handler(E_ERR,'TEC_text_to_obs',string1, &
                 source, revision, revdate, text2=string2)
   endif

   ! if lon comes in between -180 and 180, use these lines instead:
   if ( lat >   90.0_r8 ) lat =  90.0_r8
   if ( lat <  -90.0_r8 ) lat = -90.0_r8
   if ( lon <    0.0_r8 ) lon = lon + 360.0_r8
   if ( lon >= 360.0_r8 ) lon = 0.0_r8

   ! put date into a dart time format

   time_obs = set_date(year, month, day, hour, minute, second)

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

   ! make an obs derived type, and then add it to the sequence
   ! putting the observation at 350 km allows for vertical localization
   ! DART location units for VERTISHEIGHT are meters.

   call create_3d_obs(lat, lon, 350000.0_r8, VERTISHEIGHT, tec, &
               observation_type_int, observation_error_variance, oday, osec, qc, obs)
   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
! Since the files are in daily chunks, I am going to append the year/month/day
! to the output file base. Whatever is in the year/month/day variables should
! be correct.

write(string1,'(''.'',i4.4,''_'',i2.2,''_'',i2.2)')year,month,day
write(obs_out_file,'(A)')trim(obs_out_base)//trim(string1)

if ( get_num_obs(obs_seq) > 0 ) then
   write(string1, *)'obs_count = ', get_num_obs(obs_seq)
   call error_handler(E_MSG, 'TEC_text_to_obs', string1)
   call write_obs_seq(obs_seq, obs_out_file)
else
   call error_handler(E_MSG,'TEC_text_to_obs','no observations in sequence', &
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

if     (trim(observation_string) == 'GPS_VTEC_EXTRAP') then
              observation_type_int = GPS_VTEC_EXTRAP
else 
   write(string1,*)'Unable to interpret observation string "'//trim(observation_type)//'"'
   write(string2,*)'valid string is "GPS_VTEC_EXTRAP"'
   call error_handler(E_ERR, 'TEC_text_to_obs', string1, &
              source, revision, revdate, text2=string2)
endif  

end subroutine set_observation_type


end program TEC_text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
