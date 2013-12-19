! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program model_to_dart

!----------------------------------------------------------------------
! purpose: interface between TIEGCM and DART
!
! method: Read TIEGCM restart file (netCDF format).
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, finalize_utilities, &
                             error_handler, E_MSG, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : model_type, static_init_model, get_model_size, &
                             init_model_instance, read_TIEGCM_restart,      &
                             read_TIEGCM_secondary, read_TIEGCM_namelist, &
                             prog_var_to_vector, get_restart_file_name, &
                             get_aux_file_name, get_namelist_file_name 
!                            clamp_bounded_variables
use  assim_model_mod, only : open_restart_write, awrite_state_restart, close_restart
use time_manager_mod, only : time_type, print_date, print_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: file_out = 'dart_ics'

namelist /model_to_dart_nml/ file_out

!-----------------------------------------------------------------------
! global storage
!-----------------------------------------------------------------------

character(len=128) :: file_name1
character(len=128) :: file_name2
character(len=128) :: file_namelist
type(model_type)   :: var ! custome storage for native TIEGCM state
type(time_type)    :: model_time
integer            :: iunit, io, file_unit, x_size
real(r8), allocatable :: x_state(:)

!=======================================================================
! Start the program.
!=======================================================================

call initialize_utilities(progname='model_to_dart', output_flag=.true.)

! read the namelist to get output file name
call find_namelist_in_file("input.nml", "model_to_dart_nml", iunit)
read(iunit, nml = model_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "model_to_dart_nml") ! closes, too.

! static_init_model sets the geometry, model size, etc.
call static_init_model()

x_size = get_model_size()
allocate(x_state(x_size))

! Allocate an empty instance of the TIEGCM model type for storage
! This is needed because read_TIEGCM_restart() requires it.  
call init_model_instance(var)

! Read the TIEGCM state variables into var and set the model_time
! to reflect the valid time of the TIEGCM state.

file_name1    =  get_restart_file_name()
file_name2    =      get_aux_file_name()
file_namelist = get_namelist_file_name()

call read_TIEGCM_restart(  trim(file_name1),    var, model_time)
call read_TIEGCM_secondary(trim(file_name2),    var)
call read_TIEGCM_namelist( trim(file_namelist), var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! clamp the bounded variables before writing to TIEGCM restart file.
! call clamp_bounded_variables()

! write out state vector in "proprietary" format
file_unit = open_restart_write(trim(file_out))
call awrite_state_restart(model_time, x_state, file_unit)
call close_restart(file_unit)

! write a little summary
call print_date(model_time, str='model_to_dart: tiegcm date')
call print_time(model_time, str='model_to_dart: DART   time')

call error_handler(E_MSG,'model_to_dart','finished successfully.',source,revision,revdate)
call finalize_utilities()

end program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
