! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!! NOTE:  this version is the start of a utility to print out
!!  information about obs_seq files - the number of each type
!!  of obs, the times, etc.  you can get more info by running
!!  the obs_diag program, but for a quick review of a file
!!  for sanitys sake, this might help.

program query_obs_seq

! <next few lines under version control, do not edit>
! $URL$
! $Id: $
! $Revision$
! $Date: $

use        types_mod, only : r8
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit
use time_manager_mod, only : time_type, operator(>), print_time, set_time
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, assignment(=), &
                             get_num_obs, get_first_obs, get_last_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence, delete_seq_head, &
                             delete_seq_tail, get_num_key_range, set_copy_meta_data, &
                             get_obs_def
use      obs_def_mod, only : obs_def_type, get_obs_kind, get_obs_def_time
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name

implicit none

! <next few lines under version control, do not edit>
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date: $"

type(obs_sequence_type) :: seq_in
type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in, num_copies_in, num_qc_in
integer                 :: num_processed, iunit, io, i, j
integer                 :: max_num_obs, file_id, remaining_obs_count
integer                 :: first_seq
integer                 :: this_obs_kind
character(len = 129)    :: read_format, meta_data
logical                 :: pre_I_format, all_gone
logical                 :: trim_first, trim_last
character(len = 129)    :: msgstring

!----------------------------------------------------------------
! Namelist input with default values

! max_obs_kinds is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds
integer, parameter :: max_obs_input_types = max_obs_kinds
character(len = 32) :: obs_input_types(max_obs_input_types)
integer :: type_count(max_obs_input_types)
integer :: num_obs_input_types
logical :: restrict_by_obs_type

character(len = 129) :: filename_seq = 'obs_seq.out'

! Time of first and last observations to be counted from obs_sequence
! If negative, these are not used
integer  :: first_obs_days    = -1
integer  :: first_obs_seconds = -1
integer  :: last_obs_days     = -1
integer  :: last_obs_seconds  = -1

type(time_type) :: first_obs_time, last_obs_time


namelist /query_obs_seq_nml/ filename_seq,  &
         first_obs_days, first_obs_seconds, &
         last_obs_days, last_obs_seconds,   &
         obs_input_types


!----------------------------------------------------------------
! Start of the executable code.
!----------------------------------------------------------------

call query_obs_seq_modules_used()

! Initialize input obs_types
do i = 1, max_obs_input_types
   obs_input_types(i) = ""
   type_count(i) = 0
enddo

! Read the namelist entry
call find_namelist_in_file("input.nml", "query_obs_seq_nml", iunit)
read(iunit, nml = query_obs_seq_nml, iostat = io)
call check_namelist_read(iunit, io, "query_obs_seq_nml")

! FIXME: do we really want to do this?  at least not on stdout?
! Record the namelist values used for the run ...
write(nmlfileunit, nml=query_obs_seq_nml)
write(     *     , nml=query_obs_seq_nml)

! See if the user is restricting the obs types to be queried, and set up
! the values if so.
num_obs_input_types = 0
do i = 1, num_obs_input_types
   if ( len(obs_input_types(i)) .eq. 0 .or. obs_input_types(i) .eq. "" ) exit
   num_obs_input_types = i
enddo
if (num_obs_input_types == 0) then
   restrict_by_obs_type = .false.
else
   restrict_by_obs_type = .true.
endif

! Read header information for the sequence

! check to see if we are going to trim the sequence by time
if(first_obs_seconds >= 0 .or. first_obs_days >= 0) then
   first_obs_time = set_time(first_obs_seconds, first_obs_days)
   trim_first = .true.
else
   trim_first = .false.
endif
if(last_obs_seconds >= 0 .or. last_obs_days >= 0) then
   last_obs_time = set_time(last_obs_seconds, last_obs_days)
   trim_last = .true.
else
   trim_last = .false.
endif
if (trim_first .and. trim_last) then
   if (first_obs_time > last_obs_time) then
      call error_handler(E_ERR,'query_obs_seq', 'first time cannot be later than last time', &
                         source,revision,revdate)
   endif
endif

! TWO PASS algorithm; open each file, trim it if requested, and count
! the number of actual observations.  then the arrays can be created
! with the correct sizes, and various counts can be accumulated the
! second time through the file.

! pass 1:

! count up the number of observations we are going to eventually have.
! if all the observations in a file are not part of the linked list, the
! number of valid observations might be much smaller than the total size in 
! the header.  it is slower, but go ahead and read in the entire sequence
! and count up the real number of obs - trim_seq will do the count even if
! it is not trimming in time.

call read_obs_seq_header(filename_seq, num_copies_in, num_qc_in, size_seq_in, &
   max_num_obs, file_id, read_format, pre_I_format, close_the_file = .true.)

call read_obs_seq(filename_seq, 0, 0, 0, seq_in)
call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time,   &
              filename_seq, .true., remaining_obs_count)
call destroy_obs_sequence(seq_in)

! no valid obs found to begin with?
if (size_seq_in == 0) then
   msgstring = 'Input file is empty'
   call error_handler(E_ERR,'query_obs_seq',msgstring,source,revision,revdate)
endif

! no valid obs found after time trim?
if (remaining_obs_count == 0) then
   msgstring = 'All obs outside the first/last times'
   call error_handler(E_ERR,'query_obs_seq',msgstring,source,revision,revdate)
endif

! pass 2:

! Initialize individual observation variables 
call init_obs(     obs, num_copies_in, num_qc_in)
call init_obs(next_obs, num_copies_in, num_qc_in)

! FIXME: initialize counts here?
num_processed = 0

! Read obs seq and start counting up

write(msgstring,*) 'Starting to process input sequence file ', trim(filename_seq)
call error_handler(E_MSG,'query_obs_seq',msgstring,source,revision,revdate)

call read_obs_seq(filename_seq, 0, 0, 0, seq_in)

! If you get here, there better be observations in this file which
! are going to be used (the process_file flag wouldn't be set otherwise.)
call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time,   &
              filename_seq, .false., remaining_obs_count)

! This would be an error at this point.
if(remaining_obs_count == 0) then
   call destroy_obs_sequence(seq_in) 
   write(msgstring, *) 'Internal error trying to process file ', trim(filename_seq)
   call error_handler(E_ERR,'query_obs_seq',msgstring,source,revision,revdate)
endif

call print_metadata(seq_in, filename_seq)

size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

!-------------------------------------------------------------
! Start to process obs from sequence_in
!--------------------------------------------------------------
num_processed = 1
is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(adjustl(filename_seq))
   call error_handler(E_MSG,'query_obs_seq', msgstring, source, revision, revdate)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), 'First timestamp: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_kind = get_obs_kind(this_obs_def)
   type_count(this_obs_kind) = type_count(this_obs_kind) + 1
!   print *, 'obs kind index = ', this_obs_kind
!   print *, 'obs name = ', get_obs_kind_name(this_obs_kind)

!   if (mod(num_processed,1000) == 0) then
!      print*, 'processed number ',num_processed,' of ',size_seq_in
!   endif

   num_processed = num_processed + 1

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), 'Last timestamp: ')
   endif

enddo ObsLoop


call destroy_obs_sequence(seq_in)


print *, 'Number of obs processed  :          ', size_seq_in
print *, '---------------------------------------------------------'
do i = 1, max_obs_input_types
   if (type_count(i) > 0) print '(a32,i8,a)', trim(get_obs_kind_name(i)), &
                                   type_count(i), ' obs'
enddo

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

call timestamp(source,revision,revdate,'end')

contains


!=====================================================================


subroutine query_obs_seq_modules_used()

! Initialize modules used that require it
call initialize_utilities('query_obs_seq')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

end subroutine query_obs_seq_modules_used


subroutine print_metadata(seq1, fname1)
!
! print out the metadata strings, trimmed
!

type(obs_sequence_type), intent(in) :: seq1
character(len=*), optional :: fname1

integer :: num_copies , num_qc, i
character(len=129) :: str1
character(len=255) :: msgstring1

num_copies = get_num_copies(seq1)
num_qc     = get_num_qc(    seq1)

if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'query_obs_seq', msgstring1, source, revision, revdate)
endif

MetaDataLoop : do i=1, num_copies
   str1 = trim(adjustl(get_copy_meta_data(seq1,i)))

   write(msgstring1,*)'metadata ',trim(adjustl(str1))
   call error_handler(E_MSG, 'query_obs_seq', msgstring1, source, revision, revdate)

enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = trim(adjustl(get_qc_meta_data(seq1,i)))

   write(msgstring1,*)'qc metadata ', trim(adjustl(str1))
   call error_handler(E_MSG, 'query_obs_seq', msgstring1, source, revision, revdate)
   
enddo QCMetaData

end subroutine print_metadata

! pass in an already opened sequence and a start/end time.  this routine
! really trims the observations out of the sequence, and returns a count
! of how many remain.
subroutine trim_seq(seq, trim_first, first_time, trim_last, last_time, seqfilename, &
                    print_msg, remaining_obs_count)
 type(obs_sequence_type), intent(inout) :: seq
 logical, intent(in)                    :: trim_first, trim_last
 type(time_type), intent(in)            :: first_time, last_time
 character(len = *), intent(in)         :: seqfilename
 logical, intent(in)                    :: print_msg
 integer, intent(out)                   :: remaining_obs_count

   ! Need to find first obs with appropriate time, delete all earlier ones
   if(trim_first) then
      call delete_seq_head(first_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are before first_obs_days:first_obs_seconds'
            call error_handler(E_MSG,'query_obs_seq',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! Also get rid of observations past the last_obs_time if requested
   if(trim_last) then
      call delete_seq_tail(last_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are after last_obs_days:last_obs_seconds'
            call error_handler(E_MSG,'query_obs_seq',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   remaining_obs_count = get_num_key_range(seq)

end subroutine trim_seq




end program query_obs_seq
