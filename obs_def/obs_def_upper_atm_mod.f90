! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Fortran has a limit of 32 characters for variable names. Hence,
! each column can be at most 32 characters wide.
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! BEGIN DART PREPROCESS KIND LIST
! SAT_TEMPERATURE,                 KIND_TEMPERATURE,                COMMON_CODE
! SAT_TEMPERATURE_ELECTRON,        KIND_TEMPERATURE_ELECTRON,       COMMON_CODE
! SAT_TEMPERATURE_ION,             KIND_TEMPERATURE_ION,            COMMON_CODE
! SAT_DENSITY_NEUTRAL_O3P,         KIND_DENSITY_NEUTRAL_O3P,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_O2,          KIND_DENSITY_NEUTRAL_O2,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2,          KIND_DENSITY_NEUTRAL_N2,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_N4S,         KIND_DENSITY_NEUTRAL_N4S,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_NO,          KIND_DENSITY_NEUTRAL_NO,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2D,         KIND_DENSITY_NEUTRAL_N2D,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2P,         KIND_DENSITY_NEUTRAL_N2P,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_H,           KIND_DENSITY_NEUTRAL_H,          COMMON_CODE
! SAT_DENSITY_NEUTRAL_HE,          KIND_DENSITY_NEUTRAL_HE,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_CO2,         KIND_DENSITY_NEUTRAL_CO2,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_O1D,         KIND_DENSITY_NEUTRAL_O1D,        COMMON_CODE
! SAT_DENSITY_ION_O4SP,            KIND_DENSITY_ION_O4SP,           COMMON_CODE
! SAT_DENSITY_ION_O2P,             KIND_DENSITY_ION_O2P,            COMMON_CODE
! SAT_DENSITY_ION_N2P,             KIND_DENSITY_ION_N2P,            COMMON_CODE
! SAT_DENSITY_ION_NP,              KIND_DENSITY_ION_NP,             COMMON_CODE
! SAT_DENSITY_ION_NOP,             KIND_DENSITY_ION_NOP,            COMMON_CODE
! SAT_DENSITY_ION_O2DP,            KIND_DENSITY_ION_O2DP,           COMMON_CODE
! SAT_DENSITY_ION_O2PP,            KIND_DENSITY_ION_O2PP,           COMMON_CODE
! SAT_DENSITY_ION_HP,              KIND_DENSITY_ION_HP,             COMMON_CODE
! SAT_DENSITY_ION_HEP,             KIND_DENSITY_ION_HEP,            COMMON_CODE
! SAT_DENSITY_ION_E,               KIND_DENSITY_ION_E,              COMMON_CODE
! SAT_VELOCITY_U,                  KIND_VELOCITY_U,                 COMMON_CODE
! SAT_VELOCITY_V,                  KIND_VELOCITY_V,                 COMMON_CODE
! SAT_VELOCITY_W,                  KIND_VELOCITY_W,                 COMMON_CODE
! SAT_VELOCITY_U_ION,              KIND_VELOCITY_U_ION,             COMMON_CODE
! SAT_VELOCITY_V_ION,              KIND_VELOCITY_V_ION,             COMMON_CODE
! SAT_VELOCITY_W_ION,              KIND_VELOCITY_W_ION,             COMMON_CODE
! SAT_VELOCITY_VERTICAL_O3P,       KIND_VELOCITY_VERTICAL_O3P,      COMMON_CODE
! SAT_VELOCITY_VERTICAL_O2,        KIND_VELOCITY_VERTICAL_O2,       COMMON_CODE
! SAT_VELOCITY_VERTICAL_N2,        KIND_VELOCITY_VERTICAL_N2,       COMMON_CODE
! SAT_VELOCITY_VERTICAL_N4S,       KIND_VELOCITY_VERTICAL_N4S,      COMMON_CODE
! SAT_VELOCITY_VERTICAL_NO,        KIND_VELOCITY_VERTICAL_NO,       COMMON_CODE
! SAT_F107,                        KIND_1D_PARAMETER,               COMMON_CODE
! SAT_RHO,                         KIND_DENSITY,                    COMMON_CODE
! GPS_PROFILE,                     KIND_ELECTRON_DENSITY,           COMMON_CODE
! SSUSI_O_N2_RATIO,                KIND_O_N2_COLUMN_DENSITY_RATIO
! GND_GPS_VTEC,                    KIND_GND_GPS_VTEC
! GPS_VTEC_EXTRAP,                 KIND_GND_GPS_VTEC
! CHAMP_DENSITY,                   KIND_DENSITY
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_upper_atm_mod, only : get_expected_upper_atm_density
!  use obs_def_upper_atm_mod, only : get_expected_gnd_gps_vtec
!  use obs_def_upper_atm_mod, only : get_expected_gps_vtec_extrap
!  use obs_def_upper_atm_mod, only : get_expected_O_N2_ratio
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(CHAMP_DENSITY) 
!      call get_expected_upper_atm_density(state, location, obs_val, istatus)
! case(GND_GPS_VTEC)
!      call get_expected_gnd_gps_vtec(state, location, obs_val, istatus)
! case(GPS_VTEC_EXTRAP)
!      call get_expected_gps_vtec_extrap(state, location, obs_val, istatus)
! case(SSUSI_O_N2_RATIO)
!      call get_expected_O_N2_ratio(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_upper_atm_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, get_location, set_location, &
                             VERTISHEIGHT, VERTISLEVEL
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                             KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                             KIND_TEMPERATURE, &
                             KIND_PRESSURE, &
                             KIND_DENSITY_ION_E, &
                             KIND_GND_GPS_VTEC, &
                             KIND_GEOPOTENTIAL_HEIGHT, &
                             KIND_GEOMETRIC_HEIGHT, &
                             KIND_O_N2_COLUMN_DENSITY_RATIO

implicit none
private
public :: get_expected_upper_atm_density, &
          get_expected_gnd_gps_vtec, &
          get_expected_gps_vtec_extrap, &
          get_expected_O_N2_ratio

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

real(r8), PARAMETER :: N2_molar_mass = 28.0_r8
real(r8), PARAMETER :: O_molar_mass  = 16.0_r8
real(r8), PARAMETER :: O2_molar_mass = 32.0_r8
real(r8), PARAMETER :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
integer,  PARAMETER :: MAXLEVELS = 100 ! more than max levels expected in the model 
integer,       save :: nlevels = 0     ! actual number of levels in the model

character(len=512) :: string1, string2, string3

contains


subroutine initialize_module
!-----------------------------------------------------------------------------
call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module


subroutine get_expected_upper_atm_density(x, location, obs_val, istatus)
!-----------------------------------------------------------------------------
!Given DART state vector and a location, 
!it computes thermospheric neutral density [Kg/m3] 
!The istatus variable should be returned as 0 unless there is a problem
!
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
real(r8)                        :: mmro1, mmro2 ! mass mixing ratio 
real(r8)                        :: pressure, temperature 

if ( .not. module_initialized ) call initialize_module

call interpolate(x, location, KIND_ATOMIC_OXYGEN_MIXING_RATIO, mmro1, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif
call interpolate(x, location, KIND_MOLEC_OXYGEN_MIXING_RATIO, mmro2, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif
call interpolate(x, location, KIND_PRESSURE, pressure, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif
call interpolate(x, location, KIND_TEMPERATURE, temperature, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif

!density [Kg/m3] =  pressure [N/m2] * M [g/mol] / temperature [K] / R [N*m/K/kmol] 
obs_val          =  pressure &
                 /(mmro1/16.0_r8+mmro2/32.0_r8+(1-mmro1-mmro2)/28.0_r8) &
                 /temperature/universal_gas_constant 

end subroutine get_expected_upper_atm_density


subroutine get_expected_gnd_gps_vtec(state_vector, location, obs_val, istatus)
!-----------------------------------------------------------------------------
!Given DART state vector and a location, 
!it computes ground GPS vertical total electron content
!The istatus variable should be returned as 0 unless there is a problem
!
real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted total electron content that would be in the
! integrated column from an instrument looking straight down at the tangent point.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

integer  :: nAlts, iAlt
real(r8), allocatable :: ALT(:), IDensityS_ie(:) 
real(r8) :: loc_vals(3)
real(r8) :: tec
type(location_type) :: probe

if ( .not. module_initialized ) call initialize_module

istatus = 36 !initially bad return code
obs_val = MISSING_R8

! something larger than the expected number of vert levels in the model
allocate(ALT(MAXLEVELS), IDensityS_ie(MAXLEVELS))

loc_vals = get_location(location)

nAlts = 0
LEVELS: do iAlt=1, size(ALT)+1
   ! loop over levels.  if we get to one more than the allocated array size,
   ! this model must have more levels than we expected.  increase array sizes,
   ! recompile, and try again.
   if (iAlt > size(ALT)) then
      write(string1,'(''more than '',i4,'' levels in the model.'')') MAXLEVELS
      string2='increase MAXLEVELS in obs_def_upper_atm_mod.f90, rerun preprocess and recompile.'
      string3='increase ALT, IDensityS_ie array sizes in code and recompile'
      call error_handler(E_ERR, 'get_expected_gnd_gps_vtec', string1, &
           source, revision, revdate, text2=string2, text3=string3)
   endif

   ! At each altitude interpolate the 2D IDensityS_ie to the lon-lat where data 
   ! point is located. After this loop we will have a column centered at the data 
   ! point's lon-lat and at all model altitudes.
   probe = set_location(loc_vals(1), loc_vals(2), real(iAlt, r8), VERTISLEVEL) !probe is where we have data 
   call interpolate(state_vector, probe, KIND_DENSITY_ION_E, IDensityS_ie(iAlt), istatus) 
   if (istatus /= 0) exit LEVELS
   call interpolate(state_vector, probe, KIND_GEOPOTENTIAL_HEIGHT, ALT(iAlt), istatus) 
   if (istatus /= 0) exit LEVELS
   nAlts = nAlts+1
enddo LEVELS

if (nAlts == 0) return

tec=0.0_r8 !start with zero for the summation

do iAlt = 1, nAlts-1 !approximate the integral over the altitude as a sum of trapezoids
   !area of a trapezoid: A = (h2-h1) * (f2+f1)/2
   tec = tec + ( ALT(iAlt+1)-ALT(iAlt) )  * ( IDensityS_ie(iAlt+1)+IDensityS_ie(iAlt) ) /2.0_r8
enddo
obs_val = tec * 10.0**(-16) !units of TEC are "10^16" #electron/m^2 instead of just "1" #electron/m^2

deallocate(ALT, IDensityS_ie)

! Good return code. 
istatus = 0

end subroutine get_expected_gnd_gps_vtec



subroutine get_expected_gps_vtec_extrap(state_vector, location, obs_val, istatus)
!-----------------------------------------------------------------------------
! Given DART state vector and a location, 
! compute ground GPS vertical total electron content including an estimate of
! the contribution from above the model (the 'extra'polated part)
! The istatus variable should be returned as 0 unless there is a problem

real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

if ( .not. module_initialized ) call initialize_module

istatus = 1
obs_val = MISSING_R8

call error_handler(E_ERR, 'get_expected_gps_vtec_extrap', 'routine not written', &
           source, revision, revdate, &
           text2='routine in obs_def/obs_def_upper_atm_mod.f90')

! FIXME this should replace the tiegcm/model_mod:create_vtec() routine

end subroutine get_expected_gps_vtec_extrap




subroutine get_expected_O_N2_ratio(state_vector, location, obs_val, istatus)
!-----------------------------------------------------------------------------
! 
! First, find the number of levels in the model.
! Then, loop down through the levels to create a top-down vertical profile.
!       As we do that, we accumulate the amount of N2 and O, stopping when
!       the N2 reaches 10^21 M^-2. This will probably mean only using part
!       of the 'last' layer.

real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat, loc_value
type(location_type) :: loc

real(r8), parameter :: Max_N2_column_density = 1.0E22_r8

real(r8) :: voxel_number_density
real(r8) :: mmrn2
real(r8) :: N2_total
real(r8) :: O_total

real(r8) :: O_mmr(MAXLEVELS)
real(r8) :: O2_mmr(MAXLEVELS)
real(r8) :: pressure(MAXLEVELS)
real(r8) :: temperature(MAXLEVELS)
real(r8) :: heights(MAXLEVELS)
real(r8) :: thickness(MAXLEVELS)
real(r8) :: O_integrated
real(r8) :: N2_integrated

real(r8), allocatable :: N2_mmr(:)
real(r8), allocatable :: mbar(:)
real(r8), allocatable :: N2_number_density(:)
real(r8), allocatable :: O_number_density(:)
real(r8), allocatable :: O_number_density_middles(:)
real(r8), allocatable :: N2_number_density_middles(:)

integer  :: ilayer

if ( .not. module_initialized ) call initialize_module

istatus = 1
obs_val = MISSING_R8

call error_handler(E_MSG, 'get_expected_O_N2_ratio', 'routine not tested', &
           source, revision, revdate, &
           text2='routine in obs_def/obs_def_upper_atm_mod.f90')

if ( .not. module_initialized ) call initialize_module

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

nlevels   = determine_num_layers(location, state_vector)

write(*,*)'FIXME get_expected_O_N2_ratio has ',nlevels,'levels.'

! Fill arrays to get the entire column.

FILLARRAYS : do ilayer = 1,nlevels

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_vector,loc,KIND_GEOMETRIC_HEIGHT, heights(ilayer),istatus)
   if (istatus /= 0) then
      obs_val = MISSING_R8
      return
   endif

   call interpolate(state_vector, loc, KIND_PRESSURE, pressure(ilayer), istatus)
   if (istatus /= 0) then
      obs_val = MISSING_R8
      return
   endif

   call interpolate(state_vector, loc, KIND_TEMPERATURE, temperature(ilayer), istatus)
   if (istatus /= 0) then
      obs_val = MISSING_R8
      return
   endif

   call interpolate(state_vector, loc, KIND_ATOMIC_OXYGEN_MIXING_RATIO, O_mmr(ilayer), istatus)
   if (istatus /= 0) then
      obs_val = MISSING_R8
      return
   endif

   call interpolate(state_vector, loc, KIND_MOLEC_OXYGEN_MIXING_RATIO, O2_mmr(ilayer), istatus)
   if (istatus /= 0) then
      obs_val = MISSING_R8
      return
   endif

   write(*,'(''height,p,T,O1,O2 '',i2,1x,5(1x,f16.7))'), ilayer,heights(ilayer), pressure(ilayer), &
             temperature(ilayer), O_mmr(ilayer), O2_mmr(ilayer)

enddo FILLARRAYS

! calculate what we can using array notation

allocate(N2_mmr(nlevels), mbar(nlevels), O_number_density(nlevels), N2_number_density(nlevels))
allocate( O_number_density_middles(nlevels), N2_number_density_middles(nlevels))

N2_mmr = 1.0_r8 - O_mmr(1:nlevels) - O2_mmr(1:nlevels)
  mbar = 1.0_r8/( O2_mmr(1:nlevels)/O2_molar_mass + O_mmr(1:nlevels)/O_molar_mass + N2_mmr/N2_molar_mass )

 O_number_density =  O_mmr(1:nlevels) * pressure(1:nlevels) * mbar /  O_molar_mass
N2_number_density = N2_mmr(1:nlevels) * pressure(1:nlevels) * mbar / N2_molar_mass

thickness(1:nlevels-1) = heights(2:nlevels) - heights(1:nlevels-1)

 O_number_density_middles(1:nlevels-1) = &
       ( O_number_density(2:nlevels) +  O_number_density(1:nlevels-1))/2.0_r8
N2_number_density_middles(1:nlevels-1) = &
       (N2_number_density(2:nlevels) + N2_number_density(1:nlevels-1))/2.0_r8

N2_total = 0.0_r8
 O_total = 0.0_r8

TOPDOWN : do ilayer = nlevels,1,-1

   ! integrate over layer thickness
   O_integrated  =  O_number_density_middles(ilayer) * thickness(ilayer)
   N2_integrated = N2_number_density_middles(ilayer) * thickness(ilayer)

   ! if necessary, only store part of the final layer so as not to overshoot 10^21 m^-2      

   if ((N2_total + N2_integrated) >= Max_N2_column_density) then
       N2_total = Max_N2_column_density
      O_total   = O_total + O_integrated * Max_N2_column_density / (N2_total + N2_integrated)    
      exit TOPDOWN
   else
       N2_total = N2_total + N2_integrated
        O_total =  O_total +  O_integrated
   endif

enddo TOPDOWN

obs_val = O_total / N2_total
istatus = 0

deallocate(N2_mmr, mbar, O_number_density, N2_number_density)
deallocate( O_number_density_middles, N2_number_density_middles)

end subroutine get_expected_O_N2_ratio



function determine_num_layers(location, state_vector)

type(location_type), intent(in) :: location
real(r8),            intent(in) :: state_vector(:)
integer                         :: determine_num_layers

type(location_type) :: loc
real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat, loc_value
integer  :: ilayer, istatus

determine_num_layers = nlevels   ! nlevels comes from the module

! If nlevels is set already, we're done.
if (determine_num_layers > 0) return

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

! Determine the number of layers by interpolating
! lots of layers to get the height at that layer. When the
! interpolation fails, we have over-reached the number of layers in the model.
! the nlevels variable has the 'save' attribute and is avaiable to
! the entire module

determine_num_layers = 0

COUNTLEVELS : do ilayer = 1,MAXLEVELS
   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)
   call interpolate(state_vector,loc,KIND_GEOMETRIC_HEIGHT,loc_value,istatus)
   if (istatus /= 0) exit COUNTLEVELS
   determine_num_layers = determine_num_layers + 1
enddo COUNTLEVELS
   
if ((determine_num_layers == MAXLEVELS) .or. (determine_num_layers == 0)) then
   write(string1,*) 'FAILED to determine number of layers in model.'
   call error_handler(E_ERR,'determine_num_layers',string1,source,revision,revdate)
endif

end function determine_num_layers



end module obs_def_upper_atm_mod
! END DART PREPROCESS MODULE CODE      

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
