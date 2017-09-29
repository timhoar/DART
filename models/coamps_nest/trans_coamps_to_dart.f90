! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program trans_coamps_to_dart

! trans_coamps_to_dart
! --------------------
! This program pulls pieces out of the large COAMPS restart file,
! then assembles them into a state vector that can be used by DART.
! This includes two pieces of information - the current time and
! the actual state
! --------------------

use coamps_translate_mod, only : initialize_translator,         &
                                 generate_coamps_filenames,     &
                                 open_coamps_files,             &
                                 coamps_read_all_fields,        &
                                 set_dart_current_time,         &
                                 open_dart_file, dart_write,    &
                                 finalize_translator,           &
                                 record_hdf_varnames

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! The translation module uses internal flags for whether it's
! reading or writing - these are just aliases so it's clearer
! what the calls are saying
logical, parameter :: READING_COAMPS = .false.
logical, parameter :: WRITING_DART   = .true.

call initialize_translator()

call generate_coamps_filenames(READING_COAMPS)
call open_coamps_files(READING_COAMPS)
call record_hdf_varnames()

call coamps_read_all_fields()

call set_dart_current_time()

call open_dart_file(WRITING_DART)
call dart_write()

call finalize_translator()

end program trans_coamps_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
