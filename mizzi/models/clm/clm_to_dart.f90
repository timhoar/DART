! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program clm_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between clm and DART
!
! method: Read clm "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The clm filename is read from the clm_in namelist
!         <edit clm_to_dart_output_file in input.nml:clm_to_dart_nml>
!         clm_to_dart
!
! author: Tim Hoar 12 July 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, restart_file_to_sv, &
                             get_clm_restart_filename
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: clm_to_dart_output_file  = 'dart_ics'

namelist /clm_to_dart_nml/ clm_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: clm_restart_filename

!======================================================================

call initialize_utilities(progname='clm_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "clm_to_dart_nml", iunit)
read(iunit, nml = clm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "clm_to_dart_nml") ! closes, too.

call get_clm_restart_filename( clm_restart_filename )

write(*,*)
write(*,'(''clm_to_dart:converting clm restart file '',A, &
      &'' to DART file '',A)') &
       trim(clm_restart_filename), trim(clm_to_dart_output_file)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call restart_file_to_sv(clm_restart_filename, statevector, model_time) 

iunit = open_restart_write(clm_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

call print_date(model_time, str='clm_to_dart:clm  model date')
call print_time(model_time, str='clm_to_dart:DART model time')
call timestamp(string1=source, pos='end')

end program clm_to_dart

