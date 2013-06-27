! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_reseed

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : r4, r8, digits12, i8
use utilities_mod, only : register_module, error_handler, E_ERR, &
                          initialize_utilities, finalize_utilities, &
                          logfileunit, nmlfileunit, timestamp,  &
                          find_namelist_in_file, check_namelist_read, &
                          open_file, close_file, do_nml_file, do_nml_term
use time_manager_mod, only : time_type, operator(+), set_time, get_time, &
                             set_calendar_type, print_time, print_date
use random_seq_mod, only : random_seq_type, init_random_seq, &
                           random_uniform, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type (time_type) :: state_time

integer :: io
integer :: unitnum = 66, iunit
character(len=128) :: fname
integer :: reseed, hours

character(len=64) :: calendar_name = 'GREGORIAN'
integer :: nsamples = 1000000
integer :: start_day = 148866
integer :: start_sec = 0

namelist /test_reseed_nml/  calendar_name, nsamples,  &
start_day, start_sec

! --------------------------------------------

call initialize_utilities('test_reseed')
call register_module(source,revision,revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "test_reseed_nml", iunit)
read(iunit, nml = test_reseed_nml, iostat = io)
call check_namelist_read(iunit, io, "test_reseed_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_reseed_nml)
if (do_nml_term()) write(     *     , nml=test_reseed_nml)

! test no calendar like low order models use, also test
! the gregorian cal with a date close to current.

call set_calendar_type(calendar_name)
state_time = set_time(start_sec, start_day)

print *, ' '
call print_time(state_time, 'setting start time')
if (calendar_name == 'GREGORIAN') call print_date(state_time, 'is start date')

! -----

fname = 'uniform_baseline'
open(unit=unitnum, file=fname, action='write')
call test1(nsamples, .true., state_time, 66)
close(unitnum)

fname = 'gaussian_baseline'
open(unit=unitnum, file=fname, action='write')
call test1(nsamples, .false., state_time, 66)
close(unitnum)

! -----

hours = 1
reseed = 1

call dohourstest(reseed, hours)

! -----

hours = 1
reseed = 10

call dohourstest(reseed, hours)

! -----

hours = 6
reseed = 10

call dohourstest(reseed, hours)

! -----

hours = 12
reseed = 10

call dohourstest(reseed, hours)

! -----

hours = 6
reseed = 100

call dohourstest(reseed, hours)

! -----

hours = 12
reseed = 100

call dohourstest(reseed, hours)

! -----

hours = 24
reseed = 100

call dohourstest(reseed, hours)

! -----

hours = 24
reseed = 1000

call dohourstest(reseed, hours)


! -----

call finalize_utilities()

contains

!-------------------------------------------------------------------------

subroutine test1(nreps, unif, start_time, unit)
! output uniform distribution.  seeds gen with a single time
 integer,               intent(in)    :: nreps, unit
 logical,               intent(in)    :: unif
 type(time_type),       intent(in)    :: start_time

integer :: i, seed
type (random_seq_type) :: seq1

seed = generate_seed(start_time)
call init_random_seq(seq1, seed)

do i=1, nreps
   if (unif) then
      write(unit, *) random_uniform(seq1)
   else
      write(unit, *) random_gaussian(seq1, 0.0_r8, 1.0_r8)
   endif 
end do

end subroutine test1

!-------------------------------------------------------------------------

subroutine test2(nreps, per, start_time, delta_time, unif, unit)
! output uniform distribution.  generates nrep numbers, reseeding
! the generator each 'per' times.
 integer,               intent(in)    :: nreps, per, unit
 type(time_type),       intent(in)    :: start_time, delta_time
 logical,               intent(in)    :: unif

integer :: i, j, nextseed
type(time_type) :: t
type(random_seq_type) :: seq1

t = start_time
nextseed = generate_seed(t)
call init_random_seq(seq1, nextseed)

j = 1
do i=1, nreps
   if (unif) then
      write(unit, *) random_uniform(seq1)
   else
      write(unit, *) random_gaussian(seq1, 0.0_r8, 1.0_r8)
   endif
   if (j == per) then
      t = t + delta_time
      nextseed = generate_seed(t)
      call init_random_seq(seq1, nextseed)
      j = 1
   else
      j = j + 1
   endif
end do

end subroutine test2

!-------------------------------------------------------------------------

subroutine test3

real(r8), allocatable :: history(:)
real(r8) :: next_val
integer :: i, j, nextseed
type(time_type) :: t, base_time, state_time, delta_time, delta_time2
type(random_seq_type) :: seq1


base_time = set_time(0, 148000)
state_time = base_time
delta_time = set_time(1)
delta_time2 = set_time(0, 1)

print *, ' '
call print_time(state_time, 'setting start time')
call print_time(delta_time, 'delta time')
call print_time(delta_time2, 'delta time2')
print *, ' '

! generate the first number from 100 consecutive seeds
allocate(history(10000))
do i=1, 10000
   nextseed = generate_seed(state_time)

   call init_random_seq(seq1, nextseed)
   history(i) = random_uniform(seq1)

   state_time = state_time + delta_time
enddo

! now generate 100 values from each seed and see if
! they overlap in any of the previously generated vals.
state_time = base_time + delta_time2
do i=1, 1000
   nextseed = generate_seed(state_time)

   call init_random_seq(seq1, nextseed)
   next_val = random_uniform(seq1)

   do j=1, 10000-1
      if (history(j) == next_val) then
         print *, 'found match, ', i, j, history(j), next_val
         next_val = random_uniform(seq1)
         print *, 'next val from ran ', next_val, ' next in seq is ', history(j+1)
      endif
   enddo

   state_time = state_time + delta_time2

enddo

deallocate(history)

end subroutine test3

!-------------------------------------------------------------------------

subroutine dohourstest(reseed, hours)
 integer, intent(in) :: reseed, hours

type(time_type) :: delta_time

write(fname, '(A,I4.4,A,I2.2,A)') 'uniform_', reseed, 'dt_', hours, 'h'
open(unit=unitnum, file=fname, action='write')
delta_time = set_time(hours*3600, 0)
call test2(nsamples, reseed, state_time, delta_time, .true., unitnum)
close(unitnum)

write(fname, '(A,I4.4,A,I2.2,A)') 'gaussian_', reseed, 'dt_', hours, 'h'
open(unit=unitnum, file=fname, action='write')
delta_time = set_time(hours*3600, 0)
call test2(nsamples, reseed, state_time, delta_time, .false., unitnum)
close(unitnum)

end subroutine dohourstest

!-------------------------------------------------------------------------

function generate_seed(state_time)
! use the state time to set the seed for the (repeatable) random sequence 

type(time_type), intent(in) :: state_time
integer                     :: generate_seed

integer     :: days,seconds,bigint,bigtwo
integer(i8) :: bigtime,bigone
integer(i8), parameter :: secs_day = 86400_i8

call get_time(state_time, seconds, days)

! option 1: add days & secs, old gen needs negative seed,
! new one doesn't care. 
!generate_seed = days + seconds
!return

! option 2: compute total seconds in an i8 then use the lower
! 32 bits (which is all the init routine will take to stay
! portable across different platforms). 

! if generate_seed was i8, or if this routine was combined with
! the real init routine, we could avoid the bitwise and below.
! it's going to get repeated in the set seed routine.
generate_seed = iand((secs_day * days) + seconds, z'FFFFFFFF')
return

! option 3: works with old generator but maybe overkill?
! for the old generator, it has to be < 0 and apparently > something
! like 50,000 more than the largest negative int.

!bigtime       = int(days,i8)*100000_i8 + int(seconds,i8)
!bigint        = huge(seconds)        ! biggest 32bit integer
!bigone        = int(bigint,i8)       ! coerce to 64bit integer
!bigtwo        = int(mod(bigtime,bigone)) ! modulo arith on 64bit integers, then 32bit
!generate_seed = -1 * int(bigtwo)     ! coerce back to 32bit integer

! arb choices if the seed is going to cause an error
!if (generate_seed >= 0) generate_seed = -1234
!if (generate_seed < -2147432706) generate_seed = -4376
!return

! option 4: positive number instead of negative
!generate_seed = bigtwo
!return

! not reached anymore:
!write(*,*)'TJH DEBUG generate_seed ',bigtime,bigint,generate_seed

end function generate_seed

!-------------------------------------------------------------------------

end program test_reseed
