! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program create_sat_wind_obs

! Author: Yongsheng Chen
! Date: 2005/12/09
! Revision: 1.0

use types_mod,        only : r8, rad2deg, deg2rad, missing_r8
use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use time_manager_mod, only : time_type, operator(>), operator(<), &
                             operator(>=), operator(<=), operator(/=), &
                             set_date, set_calendar_type, get_time, &
                             get_date, set_time, GREGORIAN, operator(+)
use     location_mod, only : location_type, set_location
use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, set_obs_def, &
                             write_obs_seq, assignment(=), &
                             static_init_obs_sequence, destroy_obs_sequence 
use    utilities_mod, only : initialize_utilities, register_module, &
                             error_handler, timestamp, E_ERR, E_MSG, logfileunit, &
                             get_unit, open_file, find_namelist_in_file, check_namelist_read
use     obs_kind_mod, only : SAT_U_WIND_COMPONENT, SAT_V_WIND_COMPONENT

use netcdf

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source: DART/models/WRF/CYS_NEW/create_sat_wind_obs.f90,v $", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: sat_wind_obs_seq
type(obs_type)          :: obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time, time0, stime, etime

integer :: calendar_type = GREGORIAN

character(len = 129) :: output_file_name=''

integer :: max_num_obs = 800000, num_copies = 1, num_qc = 1
integer :: num_obs, which_vert

!  assume constant obs error (std)
real (r8) :: obs_err = 5.0   ! 5 m/s error

real (r8) :: obs_value(1)

character(len = 129) :: obs_file_name=''
character(len = 129) :: copy_meta_data, qc_meta_data

integer           :: dimid
integer           :: i

logical, parameter :: debug = .false.
integer            :: mode, io, iunit, var_id, ncid

integer            :: recnum
real (r8), allocatable, dimension(:) :: oblat, oblon, pressure, winddir, windspd
integer, allocatable, dimension(:)   :: validtime
real (r8)          :: oblat_valid_range(2), oblon_valid_range(2), pressure_valid_range(2), &
                              winddir_valid_range(2), windspd_valid_range(2)
real (r8)          :: uu, vv
logical            :: valid_rec
real (r8)          :: oblon360

integer            :: start_days, start_seconds      !bounding box for time
integer            :: end_days, end_seconds
logical            :: is_limit_start_time = .false., & 
                      is_limit_end_time = .false.
real (r8)          :: latmin=-90., latmax=90., lonmin=0., lonmax=360.  !bounding box for region
integer            :: rec_interval = 1

namelist /create_sat_wind_obs_nml/ latmin, latmax, lonmin, lonmax, rec_interval, &
                                   start_days,start_seconds,end_days,end_seconds, &
                                   obs_file_name, output_file_name

! --------------------------------------------------------------------------------

call initialize_utilities('create_sat_wind_obs')
call register_module(source,revision,revdate)

call find_namelist_in_file("input.nml", "create_sat_wind_obs_nml", iunit)
read(iunit, nml = create_sat_wind_obs_nml , iostat = io)
call check_namelist_read(iunit, io, "create_sat_wind_obs_nml")
call error_handler(E_MSG,'create_sat_wind_obs','create_sat_wind_obs_nml values are',' ',' ',' ')
write(logfileunit, nml=create_sat_wind_obs_nml)
write(     *     , nml=create_sat_wind_obs_nml)

rec_interval = max(1,rec_interval)                  !minimum 1 -- write out every record 
if ( lonmin < 0.0_r8 ) lonmin = lonmin + 360.0_r8   !DART use longitude 0 to 360
if ( lonmax < 0.0_r8 ) lonmax = lonmax + 360.0_r8

if(obs_file_name .eq. '') then
  print*,'Satellite Wind Observation input filename: '
  read(*,'(a)') obs_file_name
endif
if(output_file_name .eq. '') then
  print*,'DART Observation sequence filename: '
  read(*,'(a)') output_file_name
endif

call set_calendar_type(calendar_type)

! Initialize the obs_sequence module ...
call static_init_obs_sequence()

num_copies = 1
num_qc = 1
copy_meta_data = 'Satellite wind observation'    !!!! must have word "observation"
qc_meta_data = 'Satellite wind QC Index'
call init_obs_sequence(sat_wind_obs_seq,  num_copies, num_qc, max_num_obs)

do i = 1, num_copies
   call set_copy_meta_data(sat_wind_obs_seq, i, copy_meta_data)
end do
do i = 1, num_qc
   call set_qc_meta_data(sat_wind_obs_seq, i, qc_meta_data)
end do

call init_obs(obs,num_copies,num_qc)

! Open obs data file, obtain record length
mode = NF90_NOWRITE
call check ( nf90_open(obs_file_name, mode, ncid) )
call check ( nf90_inq_dimid(ncid, "recNum", dimid) )
call check ( nf90_inquire_dimension(ncid, dimid, len=recnum) )

if (recnum == 0) then
   print*, 'Empty record, stop here!'
   stop
endif


! Allocate space and retrieve observations, 
allocate(oblat(recnum),oblon(recnum),pressure(recnum),validtime(recnum),winddir(recnum),windspd(recnum))

call check ( nf90_inq_varid(ncid, "obLat", var_id) )
call check ( nf90_get_att(ncid, var_id, "valid_range", oblat_valid_range) )
call check ( nf90_get_var(ncid, var_id, oblat) )

call check ( nf90_inq_varid(ncid, "obLon", var_id) )
call check ( nf90_get_att(ncid, var_id, "valid_range", oblon_valid_range) )
call check ( nf90_get_var(ncid, var_id, oblon) )

call check ( nf90_inq_varid(ncid, "pressure", var_id) )
call check ( nf90_get_att(ncid, var_id, "valid_range", pressure_valid_range) )
call check ( nf90_get_var(ncid, var_id, pressure) )

call check ( nf90_inq_varid(ncid, "validTime", var_id) )  
call check ( nf90_get_var(ncid, var_id, validtime) )

call check ( nf90_inq_varid(ncid, "windDir", var_id) )
call check ( nf90_get_att(ncid, var_id, "valid_range", winddir_valid_range) )
call check ( nf90_get_var(ncid, var_id, winddir) )

call check ( nf90_inq_varid(ncid, "windSpd", var_id) )
call check ( nf90_get_att(ncid, var_id, "valid_range", windspd_valid_range) )
call check ( nf90_get_var(ncid, var_id, windspd) )


! get time for the record reference starting time (1970-1-1 00:00:0.0)
!     and  for time bounding box
time0 = set_date(1970,1,1,0,0,0)
if (start_seconds >= 0 .and. start_days >= 0 ) then
  is_limit_start_time = .true.
  stime = set_time(start_seconds, start_days)
endif
if (end_seconds   >=0  .and. end_days   >= 0 ) then
  is_limit_end_time = .true.
  etime = set_time(  end_seconds,   end_days)
endif

num_obs = 0
! Process observations

do i=1,recnum, rec_interval

   oblon360 = oblon(i)
   if(oblon360 < 0.0_r8 ) oblon360 = oblon360+360.0_r8

   !validTime is seconds since the reference time 1970-1-1 0:0:0
   time = time0 + set_time(validtime(i))

   !Validate the observation, if not valid, skip
   !also validate if the observation output the bounding box
   valid_rec = .true.
   if ( oblat(i) < oblat_valid_range(1) .or. oblat(i) > oblat_valid_range(2) .or. &
        oblon(i) < oblon_valid_range(1) .or. oblon(i) > oblon_valid_range(2) .or. &
        pressure(i) < pressure_valid_range(1) .or. pressure(i) > pressure_valid_range(2) .or. &
        winddir(i) < winddir_valid_range(1) .or. winddir(i) > winddir_valid_range(2) .or. &
        windspd(i) < windspd_valid_range(1) .or. windspd(i) > windspd_valid_range(2) .or. &
        oblat(i) < latmin .or. oblat(i) > latmax .or. oblon360 < lonmin .or. oblon360 > lonmax  .or. &
        ( is_limit_start_time .and. time<stime) .or. (is_limit_end_time .and. time>etime) ) &
      valid_rec = .false.


!  ----------------------------------------------------------------
!  if block to eject observation outside the selected range
!  ----------------------------------------------------------------
   if( valid_rec ) then

     call set_obs_def_time(obs_def, time)
   
     ! Dart longitude from 0 to 360
     if(oblon(i) < 0.0_r8) oblon(i) = oblon(i) + 360.0_r8

     which_vert = 2     !-1 for surface, but what will be the value for vloc?
                        ! 1 for model level
                        ! 2 for pressure
                        ! 3 for height
     location = set_location(oblon(i),oblat(i),pressure(i),which_vert)
     call set_obs_def_location(obs_def, location)

     !convert wind speed and direction to zonal and meridional components
     uu = missing_r8
     vv = missing_r8
     if ( windspd(i) /= 0.0 .and. abs(winddir(i)) < 1.e5 ) then
        winddir(i) = winddir(i)*DEG2RAD
        uu = -windspd(i)*sin(winddir(i))
        vv = -windspd(i)*cos(winddir(i))
     endif

     ! Insert obs into DART obs sequence
   
     obs_value(1) = uu

     call set_obs_def_kind(obs_def, SAT_U_WIND_COMPONENT)
     call set_obs_def_error_variance(obs_def, obs_err*obs_err)
     call set_obs_def(obs, obs_def)
     call set_obs_values(obs, obs_value)
     call insert_obs_in_seq(sat_wind_obs_seq, obs)
     num_obs = num_obs+1

     obs_value(1) = vv
     call set_obs_def_kind(obs_def, SAT_V_WIND_COMPONENT)
     call set_obs_def_error_variance(obs_def, obs_err*obs_err)
     call set_obs_def(obs, obs_def)
     call set_obs_values(obs, obs_value)
     call insert_obs_in_seq(sat_wind_obs_seq, obs)
     num_obs = num_obs + 1

   endif
!  --------------------------------------------------------
!  end if block to eject observation outside the selected range
!  --------------------------------------------------------

enddo

25 print*,'End of Sat Wind data'

print*,'Total number of observations in the sequence :', num_obs

call write_obs_seq(sat_wind_obs_seq, output_file_name)

! release the memory of the seq.
call destroy_obs_sequence(sat_wind_obs_seq)

call timestamp(source,revision,revdate,'end') ! close the log file.

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'create_sat_wind_obs', &
         trim(nf90_strerror(istatus)), source, revision, revdate)
  end subroutine check

end program create_sat_wind_obs
