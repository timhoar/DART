! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_model

!----------------------------------------------------------------------
! purpose: interface between TIEGCM and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into TIEGCM fields.
!         Replace those fields on the TIEGCM restart file with the new values,
!
!         Replace 'mtime' variable in the TIEGCM restart file
!         with the 'valid time' of the DART state vector.
!         Write out updated namelist variables (e.g., model_time, adv_to_time)
!         in a file called 'namelist_update'
!
!         The dart_to_model_nml namelist setting for advance_time_present
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, E_ERR, E_MSG, &
                             error_handler, finalize_utilities, do_output, &
                             find_namelist_in_file, check_namelist_read,   &
                             logfileunit
use        model_mod, only : get_model_size, static_init_model, get_f107_value, &
                             get_restart_file_name, update_TIEGCM_restart
use  assim_model_mod, only : aread_state_restart, open_restart_read, close_restart
use time_manager_mod, only : time_type, get_time, get_date, set_calendar_type, &
                             print_time, print_date, set_date, set_time, &
                             operator(*),  operator(+), operator(-), &
                             operator(>),  operator(<), operator(/), &
                             operator(/=), operator(<=)

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: file_in              = 'dart_restart'
character(len=128) :: file_namelist_out    = 'namelist_update'
logical            :: advance_time_present = .false.

namelist /dart_to_model_nml/ file_in, file_namelist_out, advance_time_present

!-----------------------------------------------------------------------
! global storage
!-----------------------------------------------------------------------

character(len=128)     :: file_name
type(time_type)        :: model_time, adv_to_time, jan1, tbase, target_time
real(r8), allocatable  :: x_state(:)
integer                :: iunit, io, file_unit, x_size

!======================================================================

call initialize_utilities(progname='dart_to_model', output_flag=.true.)

call find_namelist_in_file("input.nml", "dart_to_model_nml", iunit)
read(iunit, nml = dart_to_model_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_model_nml")

call static_init_model()        ! reads input.nml, etc., sets the table
x_size = get_model_size()       ! now that we know how big state vector is ...
allocate(x_state(x_size))       ! allocate space for the (empty) state vector

! Open the DART model state ...
! (possibly) Read in the time to which TIEGCM must advance.
! Read in the valid time for the model state
! Read in state vector from DART

file_unit = open_restart_read(trim(file_in))

if ( advance_time_present ) then
   call aread_state_restart(model_time, x_state, file_unit, adv_to_time)
else
   call aread_state_restart(model_time, x_state, file_unit)
endif
call close_restart(file_unit)

if (do_output()) then
    call print_time(model_time,'time of restart file '//trim(file_in))
    call print_date(model_time,'date of restart file '//trim(file_in))
endif

if (do_output() .and. advance_time_present) then
   call print_time(adv_to_time,'advance-to-time for  '//trim(file_in))
   call print_date(adv_to_time,'advance-to-date for  '//trim(file_in))
endif

! write fields to the binary TIEGCM restart file

file_name = get_restart_file_name()
call update_TIEGCM_restart(x_state, trim(file_name), model_time)

if ( advance_time_present ) then
   ! update TIEGCM namelist variables used in advance_model.csh
   call write_tiegcm_time_control(trim(file_namelist_out), model_time, adv_to_time)
endif

deallocate(x_state)

write(     *     ,*)''
write(logfileunit,*)''
call error_handler(E_MSG,'dart_to_model','Finished successfully.',source,revision,revdate)
call finalize_utilities()


!======================================================================
contains
!======================================================================


subroutine write_tiegcm_time_control(filename, model_time, adv_to_time)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: model_time
type(time_type),  intent(in) :: adv_to_time

integer :: model_doy, model_hour, model_minute    ! current tiegcm mtime (doy,hour,mitute)
integer :: adv_to_doy, adv_to_hour, adv_to_minute ! advance tiegcm mtime
integer :: target_doy, target_hour, target_minute ! forecast time interval in mtime format
integer :: utsec, year, month, day, sec, model_year

real(r8) :: f10_7

! Write updated TIEGCM namelist variables to a text file.
! It is up to advance_model.csh to update the TIEGCM namelist.

file_unit = get_unit()
open(unit = file_unit, file = trim(filename))

call get_date(model_time, model_year, month, day, model_hour, model_minute, sec)
jan1  = set_date(model_year,1,1)
tbase = model_time - jan1               ! total time since the start of the year
call get_time(tbase, utsec, model_doy)
model_doy = model_doy + 1               ! add jan1 back in

call get_date(adv_to_time, year, month, day, adv_to_hour, adv_to_minute, sec)
jan1  = set_date(year,1,1)
tbase = adv_to_time - jan1              ! total time since the start of the year
call get_time(tbase, utsec, adv_to_doy)
adv_to_doy = adv_to_doy + 1             ! add jan1 back in

! Calculate number of hours to advance tiegcm
target_time = adv_to_time - model_time
call get_time(target_time, sec, day)
target_doy    = day
target_hour   = sec/3600
target_minute = (sec - target_hour*3600)/60

!START_YEAR
write(file_unit, *, iostat = io )  model_year
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write model_year to '//trim(filename), &
         source,revision,revdate)
endif

!START_DAY
write(file_unit, *, iostat = io )  model_doy
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write model_day to '//trim(filename), &
         source,revision,revdate)
endif

!SOURCE_START, START, SECSTART
write(file_unit, *, iostat = io )  model_doy,',',model_hour,',',model_minute
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write (day,hour,minute) to '//trim(filename), &
         source,revision,revdate)
endif

!STOP, SECSTOP
write(file_unit, *, iostat = io )  adv_to_doy,',',adv_to_hour,',',adv_to_minute
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write adv_to_time (day,hour,minute) to '//trim(filename), &
         source,revision,revdate)
endif

!HIST, SAVE, SECHIST, SECSAVE
write(file_unit, *, iostat = io )  target_doy,',',target_hour,',',target_minute
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write target_time (day,hour,minute) to '//trim(filename), &
         source,revision,revdate)
endif

!F107
f10_7 = get_f107_value(x_state)

write(file_unit, *, iostat = io ) f10_7 
if (io /= 0 )then
   call error_handler(E_ERR,'dart_to_model:','cannot write f107 to '//trim(filename), &
        source,revision,revdate)
endif

close(file_unit)

end subroutine write_tiegcm_time_control

end program dart_to_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
