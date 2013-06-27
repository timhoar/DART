! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_track_obs_interactive

! Author: Yongsheng Chen
! Date: 2005/12/09
! Revision: 1.0

use types_mod,        only : r8, rad2deg
use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use time_manager_mod, only : time_type, operator(>), operator(<), operator(>=), &
                             operator(/=), set_date, set_calendar_type, get_time, &
                             get_date, set_time, GREGORIAN
use     location_mod, only : location_type, set_location
use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, set_obs_def, &
                             write_obs_seq, assignment(=), &
                             static_init_obs_sequence, destroy_obs_sequence 
use    utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, timestamp, E_ERR, E_MSG, logfileunit, &
                             get_unit, open_file
use     obs_kind_mod, only : KIND_VORTEX_LAT, KIND_VORTEX_LON, &
                             KIND_VORTEX_PMIN, KIND_VORTEX_WMAX, &
                             VORTEX_LAT, VORTEX_LON, &
                             VORTEX_PMIN, VORTEX_WMAX


implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: DART/models/WRF/CYS_NEW/create_track_obs_interactive.f90,v $", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: track_obs_seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time, stime, etime

integer :: calendar_type = GREGORIAN

character(len = 129) :: output_file_name='obs_seq.out'
integer :: iunit, i, end_of_file
integer :: year, month, day, hour ,minute, second
integer :: start_month, start_day, start_hour
integer :: end_month, end_day, end_hour
integer :: select_obs_range = 0, in_obs_range = 0
integer :: select_obs = 0
!  if select_obs = 1, select output observations using the logicals:
!  out_obs_pos   out_obs_pmin   out_obs_wmax
logical :: out_obs_pos  = .true., out_obs_pmin  = .true., out_obs_wmax  = .false.

integer :: max_num_obs = 800000, num_copies = 1, num_qc = 1
integer :: num_obs, which_vert, more_obs

real (r8) :: vloc

!  assume constant obs error (std)
real (r8) :: obs_pos_err = 0.3, &
             obs_pmin_err = 200.0, &
             obs_wmax_err = 5.0
!            obs_pmin_err = 3000.0, &    !make it large 30 hPa to offset strong vortex

character (len=3) :: adv 
character (len=9) :: obstime 
character (len=20) :: stat 
real (r8) :: clat, clon, wmax, pmin, obs_value(1), clon360
character (len=120) :: track_fmt_str = '(a3,f7.2,1x,f7.2,1x,I2,1x,I2,1x,I2,2x,f4.0,1x,f5.0,a20)'

!   eg:   '15A  17.00  -82.20 10/19/06Z  130   901 HURRICANE-4'

character(len = 129) :: date_str, storm_id_str, item_name_str
character(len = 129) :: obs_file_name = 'track.dat'
character(len = 129) :: copy_meta_data, qc_meta_data

integer :: track_format


call initialize_utilities('create_track_obs_interactive')
call register_module(source,revision,revdate)

print*,'----------- Interactive Input Track Observations -------------'
print*,'DART Observation sequence filename: '
read(*,'(a)') output_file_name

call set_calendar_type(calendar_type)

! Initialize the obs_sequence module ...
call static_init_obs_sequence()

num_copies = 1
num_qc = 1
copy_meta_data = 'Track observation'    !!!! must have word "observation"
qc_meta_data = 'Track QC Index'
call init_obs_sequence(track_obs_seq,  num_copies, num_qc, max_num_obs)

do i = 1, num_copies
   call set_copy_meta_data(track_obs_seq, i, copy_meta_data)
end do
do i = 1, num_qc
   call set_qc_meta_data(track_obs_seq, i, qc_meta_data)
end do

call init_obs(obs,num_copies,num_qc)
call init_obs(prev_obs,num_copies,num_qc)

num_obs = 0
more_obs = 1

print*,'Track info format: '
print*,' 0 = free format:  yyyy mm dd hh lat lon pmin'
print*,' 1 = Unisys     :  '//track_fmt_str
print*,' other = free format for now'
read(*,*) track_format

if (track_format .eq. 1) then
   print*,'Read Unisys Track Info'
   print*,'year= '
   read(*,*) year
endif

do while (more_obs == 1)

   if (track_format .eq. 1) then
      read(*,track_fmt_str) adv, clat, clon, month, day, hour, wmax, pmin, stat
      print*, adv, clat, clon, month, day, hour, wmax, pmin, stat
   else
      print*,'... yyyy mm dd hh  lat (deg)  lon (deg) pmin (hPa,-1 for missing) ...'
      read(*,*) year, month, day, hour, clat, clon, pmin
   endif

   minute = 0
   second = 0
   time = set_date(year,month,day,hour,minute,second)
   call set_obs_def_time(obs_def, time)
   
   ! Insert obs into DART obs sequence
   
   if ( out_obs_pos ) then          !output position
      clon360 = clon
      if(clon360 < 0.0_r8) clon360 = clon360 + 360.0_r8
      which_vert = -1     !-1 for surface, but what will be the value for vloc?
                          ! 1 for model level
      vloc=1.0
      location = set_location(clon360,clat,vloc,which_vert)
      call set_obs_def_location(obs_def, location)

      num_obs = num_obs + 1
      obs_value(1) = clat
      call set_obs_def_kind(obs_def, VORTEX_LAT)
      call set_obs_def_error_variance(obs_def, obs_pos_err*obs_pos_err)
      call set_obs_def(obs, obs_def)
      call set_obs_values(obs, obs_value)
!     if ( num_obs == 1 ) then
         call insert_obs_in_seq(track_obs_seq, obs)
!     else
!        call insert_obs_in_seq(track_obs_seq, obs, prev_obs)
!     endif
      prev_obs = obs

      num_obs = num_obs + 1
      obs_value(1) = clon
      call set_obs_def_kind(obs_def, VORTEX_LON)
      call set_obs_def_error_variance(obs_def, obs_pos_err*obs_pos_err)
      call set_obs_def(obs, obs_def)
      call set_obs_values(obs, obs_value)
!     call insert_obs_in_seq(track_obs_seq, obs, prev_obs)
         call insert_obs_in_seq(track_obs_seq, obs)
      prev_obs = obs

   endif

   if ( out_obs_pmin ) then        !output minimum pressure

      if ( pmin > 0) then 
      num_obs = num_obs + 1
      obs_value(1) = pmin * 100._r8  !convert to Pa
      call set_obs_def_kind(obs_def, VORTEX_PMIN)
      call set_obs_def_error_variance(obs_def, obs_pmin_err*obs_pmin_err)
      call set_obs_def(obs, obs_def)
      call set_obs_values(obs, obs_value)
!     if ( num_obs == 1 ) then
         call insert_obs_in_seq(track_obs_seq, obs)
!     else
!        call insert_obs_in_seq(track_obs_seq, obs, prev_obs)
!     endif
      prev_obs = obs
      endif  ! if (pmin>0)
   endif

   if ( out_obs_wmax) then         !output maximum wind speed
      print*,'... Intensity (Wind Max) in m/s :'
      read(*,*) wmax
      num_obs = num_obs + 1
!     obs_value(1) = wmax * 0.5144_r8  !convert knot to m/s
      obs_value(1) = wmax
      call set_obs_def_kind(obs_def, VORTEX_WMAX)
      call set_obs_def_error_variance(obs_def, obs_wmax_err*obs_wmax_err)
      call set_obs_def(obs, obs_def)
      call set_obs_values(obs, obs_value)
!     if ( num_obs == 1 ) then
         call insert_obs_in_seq(track_obs_seq, obs)
!     else
!        call insert_obs_in_seq(track_obs_seq, obs, prev_obs)
!     endif
      prev_obs = obs
   endif

   print*,' ---- more observations? (1=yes, 0=no) ------'
   read(*,*) more_obs

enddo

25 print*,'End of track data'

call write_obs_seq(track_obs_seq, output_file_name)

print*,'before destroy_obs_sequence'
! release the memory of the seq.
call destroy_obs_sequence(track_obs_seq)
print*,'after destroy_obs_sequence'

call timestamp(source,revision,revdate,'end') ! close the log file.

goto 88

! Error messages
99 continue
call error_handler(E_ERR,'create_track_obs', &
     'Error when reading header of the track data file', source, revision, revdate)
goto 88
77 continue
call error_handler(E_ERR,'create_track_obs', &
     'Error when reading track data', source, revision, revdate)
goto 88

88 continue

end program create_track_obs_interactive
