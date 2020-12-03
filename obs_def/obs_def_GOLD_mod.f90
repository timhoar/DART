! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! This module supports the observation types from the GOLD instruments.  
! https://gold.cs.ucf.edu/data/

! BEGIN DART PREPROCESS KIND LIST
! GOLD_TEMPERATURE,          KIND_TEMPERATURE,           COMMON_CODE
! GOLD_NEMAX,                KIND_NEMAX_DISK
! GOLD_ON2COLUMN,            KIND_ON2_DISK
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_GOLD_mod, only : get_expected_nemax
!  use obs_def_GOLD_mod, only : get_expected_on2
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(GOLD_NEMAX) 
!      call get_expected_nemax(state, location, obs_val, istatus)
! case(GOLD_ON2COLUMN) 
!      call get_expected_on2(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(GOLD_NEMAX) 
!      continue
! case(GOLD_ON2COLUMN) 
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(GOLD_NEMAX) 
!      continue
! case(GOLD_ON2COLUMN) 
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(GOLD_NEMAX) 
!      continue
! case(GOLD_ON2COLUMN) 
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_GOLD_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, get_location, set_location, &
                             VERTISHEIGHT, VERTISLEVEL, VERTISUNDEF
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_NEMAX_DISK, &
                             KIND_ON2_DISK, &
                             KIND_TEMPERATURE, &
                             KIND_ELECTRON_DENSITY, &
                             KIND_GEOMETRIC_HEIGHT, &
                             KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                             KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                             KIND_PRESSURE, &
                             KIND_TEMPERATURE

implicit none
private
public :: get_expected_nemax, &
          get_expected_on2

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_GOLD_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

real(r8), PARAMETER       :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
real(r8), PARAMETER       :: boltzmann_constant = 1.38064852E-23_r8 ! [J/K]
logical, save             :: module_initialized = .false.

contains

!-----------------------------------------------------------------------------

subroutine initialize_module

  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module

!-----------------------------------------------------------------------------
!> Given DART state vector and a location, 
!> retrieve ionospheric maximum electron density in column [/m3] 
!> The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_nemax(x, location, obs_val, istatus)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

integer  :: nAlts, iAlt
real(r8), allocatable :: ALT(:), IDensityS_ie(:) 
real(r8) :: loc_vals(3)
real(r8) :: nemax
type(location_type) :: probe

if ( .not. module_initialized ) call initialize_module

istatus = 36 !initially bad return code
obs_val = MISSING_R8

! something larger than the expected number of vert levels in the model
! TIE-GCM should have either 29 or 57 vertical levels
allocate(ALT(500), IDensityS_ie(500))

ALT          = 0.0_r8
IDensityS_ie = 0.0_r8

loc_vals = get_location(location)

nAlts = 0
LEVELS: do iAlt=1, size(ALT)+1
   ! loop over levels.  if we get to one more than the allocated array size,
   ! this model must have more levels than we expected.  increase array sizes,
   ! recompile, and try again.
   if (iAlt > size(ALT)) then
      call error_handler(E_ERR, 'get_expected_nemax', 'more than 500 levels in model', &
           source, revision, revdate, &
           text2='increase ALT, IDensityS_ie array sizes in code and recompile')
   endif
   ! At each altitude interpolate the 2D IDensityS_ie to the lon-lat where data 
   ! point is located after this loop we will get a column centered at data 
   ! point's lon-lat and at all model altitudes
   probe = set_location(loc_vals(1), loc_vals(2), real(iAlt, r8), VERTISLEVEL) !probe is where we have data

   call interpolate(x, probe, KIND_ELECTRON_DENSITY, IDensityS_ie(iAlt), istatus) 
   if (istatus /= 0) exit LEVELS

   nAlts = nAlts+1
enddo LEVELS

if (nAlts == 0) return

nemax = maxval(IDensityS_ie(1:nAlts))
obs_val = nemax

end subroutine get_expected_nemax


!-----------------------------------------------------------------------------
!> Given DART state vector and a location, compute column abundance of O
!> relative to N2 down to N2 reference depth of 1e17 [/cm2]
!>    
!> The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_on2(x, location, obs_val, istatus)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

integer  :: nAlts, iAlt
real(r8), allocatable :: ALT(:)
real(r8) :: loc_vals(3)
real(r8) :: o1_ratio, o2_ratio, press, temp, nd
real(r8), allocatable :: o1_number(:), n2_number(:)
real(r8) :: o1_sum, n2_sum
real(r8) :: o1_inc, n2_inc
type(location_type) :: probe

real(r8), PARAMETER       :: n2_sum_ref = 1.0E21

if ( .not. module_initialized ) call initialize_module

istatus = 36 !initially bad return code
obs_val = MISSING_R8

! something larger than the expected number of vert levels in the model
! TIE-GCM should have either 29 or 57 vertical levels
allocate(ALT(500), o1_number(500), n2_number(500))
ALT       = 0.0_r8
o1_number = 0.0_r8
n2_number = 0.0_r8

loc_vals = get_location(location)

nAlts = 0
LEVELS: do iAlt=1, size(ALT)+1
   ! loop over levels.  if we get to one more than the allocated array size,
   ! this model must have more levels than we expected.  increase array sizes,
   ! recompile, and try again.
   if (iAlt > size(ALT)) then
      call error_handler(E_ERR, 'get_expected_on2', 'more than 500 levels in model', &
           source, revision, revdate, &
           text2='increase ALT, IDensityS_ie array sizes in code and recompile')
   endif
   ! At each altitude interpolate the 2D IDensityS_ie to the lon-lat where data 
   ! point is located after this loop we will get a column centered at data 
   ! point's lon-lat and at all model altitudes
   probe = set_location(loc_vals(1), loc_vals(2), real(iAlt, r8), VERTISLEVEL) !probe is where we have data
   call interpolate(x, probe, KIND_ATOMIC_OXYGEN_MIXING_RATIO, o1_ratio, istatus) 
   if (istatus /= 0) exit LEVELS
   call interpolate(x, probe, KIND_MOLEC_OXYGEN_MIXING_RATIO, o2_ratio, istatus) 
   if (istatus /= 0) exit LEVELS
   call interpolate(x, probe, KIND_GEOMETRIC_HEIGHT, ALT(iAlt), istatus)
   if (istatus /= 0) exit LEVELS
   call interpolate(x, probe, KIND_PRESSURE, press, istatus) 
   if (istatus /= 0) exit LEVELS
   call interpolate(x, probe, KIND_TEMPERATURE, temp, istatus) 
   if (istatus /= 0) exit LEVELS
   nd = press / (boltzmann_constant * temp)
   o1_number(iAlt) = o1_ratio * nd
   n2_number(iAlt) = (1 - o1_ratio - o2_ratio) * nd
   nAlts = nAlts+1
enddo LEVELS

if (nAlts == 0) return

n2_sum = 0.0
o1_sum = 0.0

TRAPZ: do iAlt = nAlts-1,1,-1 !approximate the integral over the altitude as a sum of trapezoids
   !area of a trapezoid: A = (h2-h1) * (f2+f1)/2
   n2_inc = ( ALT(iAlt+1)-ALT(iAlt) )  * ( n2_number(iAlt+1)+n2_number(iAlt) ) /2.0_r8
   o1_inc = ( ALT(iAlt+1)-ALT(iAlt) )  * ( o1_number(iAlt+1)+o1_number(iAlt) ) /2.0_r8
   if ((n2_sum + n2_inc) .lt. n2_sum_ref) then
      n2_sum = n2_sum + n2_inc
      o1_sum = o1_sum + o1_inc
   else
      o1_sum = o1_sum + (o1_inc * (n2_sum_ref - n2_sum) / n2_inc)
      n2_sum = n2_sum_ref
      exit TRAPZ
   end if
enddo TRAPZ

! TJH is it possible that n2_sum (or n2_inc) is zero ...

obs_val = o1_sum / n2_sum

end subroutine get_expected_on2

    
end module obs_def_GOLD_mod
! END DART PREPROCESS MODULE CODE      

