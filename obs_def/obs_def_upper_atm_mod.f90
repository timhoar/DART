! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

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
! GPS_PROFILE,                     KIND_ELECTRON_DENSITY,           COMMON_CODE
! GPS_VTEC_EXTRAP,                 KIND_VERTICAL_TEC,               COMMON_CODE
! COSMIC_ELECTRON_DENSITY,         KIND_ELECTRON_DENSITY,           COMMON_CODE
! SAT_RHO,                         KIND_DENSITY
! CHAMP_DENSITY,                   KIND_DENSITY
! GND_GPS_VTEC,		           KIND_GND_GPS_VTEC
! GROUND_SLANT_TEC,                KIND_SLANT_TEC
! SSUSI_O_N2_RATIO,                KIND_O_N2_COLUMN_DENSITY_RATIO
! MIDAS_TEC,                       KIND_VERTICAL_TEC
! END DART PREPROCESS KIND LIST

! NOTE: GPS_VTEC_EXTRAP can come from the COMMON_CODE because the tiegcm model_mod.f90
! creates a KIND_VERTICAL_TEC on-demand as it reads in the model state. Since not
! all models can be expected to do this, there is a stub routine
! get_expected_gps_vtec_extrap() that will need to be written. At that time, the
! DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF block will need to be extended and the
! COMMON_CODE removed from the table above.  TJH 17 May 2016

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_upper_atm_mod, only : get_expected_upper_atm_density, &
!                                    get_expected_gnd_gps_vtec, &
!                                    get_expected_gps_vtec_extrap, &
!                                    get_expected_O_N2_ratio, &
!                                    get_expected_slant_tec, &
!                                    read_slant_tec_metadata, &
!                                    write_slant_tec_metadata, &
!                                    interactive_slant_tec_metadata, &
!                                    get_expected_vtec, &
!                                    get_expected_O_N2_ratio
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! NOTE:
! CHAMP_DENSITY  can be created with observations/CHAMP/CHAMP_density_text_to_obs

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(SAT_RHO, CHAMP_DENSITY)
!      call get_expected_upper_atm_density(state, location, obs_val, istatus)
! case(GND_GPS_VTEC)
!      call get_expected_gnd_gps_vtec(state, location, obs_val, istatus)
! case(SSUSI_O_N2_RATIO)
!      call get_expected_O_N2_ratio(state, location, obs_val, istatus)
! case(GROUND_SLANT_TEC)
!      call get_expected_slant_tec(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(SAT_RHO) 
!      continue
! case(CHAMP_DENSITY)
!      continue
! case(MIDAS_TEC) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! case(GROUND_SLANT_TEC)
!      call read_slant_tec_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(SAT_RHO) 
!      continue
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! case(GROUND_SLANT_TEC)
!      call write_slant_tec_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(SAT_RHO) 
!      continue
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! case(GROUND_SLANT_TEC)
!      call interactive_slant_tec_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_upper_atm_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             ascii_file_format
use     location_mod, only : location_type, get_location, set_location, &
                             VERTISHEIGHT, VERTISLEVEL, interactive_location, &
                             read_location, write_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                             KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                             KIND_TEMPERATURE, &
                             KIND_PRESSURE, &
                             KIND_DENSITY, &
                             KIND_DENSITY_ION_E, &
                             KIND_GND_GPS_VTEC, &
                             KIND_GEOMETRIC_HEIGHT, &
                             KIND_O_N2_COLUMN_DENSITY_RATIO

implicit none
private
public :: get_expected_upper_atm_density, &
          get_expected_gnd_gps_vtec, &
          get_expected_gps_vtec_extrap, &
          get_expected_O_N2_ratio, &
          get_expected_slant_tec, &
          read_slant_tec_metadata, &
          write_slant_tec_metadata, &
          interactive_slant_tec_metadata, &
          get_expected_vtec

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_upper_atm_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.

! Derived type for slant tec observatons.
! Contains auxiliary information stored with each obs of this type.
! The information is needed to compute the forward operator.

type slant_tec_type
   private
   type(location_type) :: receiver_location
   type(location_type) :: transmitter_location
end type slant_tec_type
type(slant_tec_type), allocatable :: slant_tec_metadata(:)

integer :: MAX_slant_tec_key = 1000
integer :: slant_tec_counter = 0 ! Cumulative index into slant_tec_metadata array

integer, PARAMETER :: labellength = 5
character(len=labellength) :: RECEIVERSTRING = 'recvr'
character(len=labellength) :: TRANSMITTERSTRING = 'trans'

! constants required to convert to same units as observations

real(r8), PARAMETER :: N2_molar_mass = 28.0_r8 ! [g/mol]
real(r8), PARAMETER ::  O_molar_mass = 16.0_r8 ! [g/mol]
real(r8), PARAMETER :: O2_molar_mass = 32.0_r8 ! [g/mol]
real(r8), PARAMETER :: N2_molar_mass_kg = N2_molar_mass/1000.0_r8 ! [kg/mol]
real(r8), PARAMETER ::  O_molar_mass_kg =  O_molar_mass/1000.0_r8 ! [kg/mol]
real(r8), PARAMETER :: O2_molar_mass_kg = O2_molar_mass/1000.0_r8 ! [kg/mol]
real(r8), PARAMETER :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! [m^2 kg / s^2 / K]
integer,  PARAMETER :: MAXLEVELS = 100 ! more than max levels expected in the model

character(len=512) :: string1, string2, string3

contains


!-----------------------------------------------------------------------
!>

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

allocate(slant_tec_metadata(MAX_slant_tec_key))
slant_tec_counter = 0

end subroutine initialize_module


!-----------------------------------------------------------------------
!> Given DART state vector and a location,
!> it computes thermospheric neutral density [Kg/m3]
!> The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_upper_atm_density(x, location, obs_val, istatus)

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
real(r8)                        :: mmro1, mmro2 ! mass mixing ratio 
real(r8)                        :: mass_reciprocal, pressure, temperature 

if ( .not. module_initialized ) call initialize_module

! Some models (i.e. GITM) have density as part of the state.
! If it is available, just return it. If density is not state,
! then we need to create it from its constituents.

call interpolate(x, location, KIND_DENSITY, obs_val, istatus)
if (istatus == 0) return

! This part was implemented for TIEGCM. Check the units for use with
! other models.

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

! density [Kg/m3] =  pressure [N/m2] * M [g/mol] / temperature [K] / R [N*m/K/kmol] 
! where M is the mean molar mass 
! 1/M = sum(wi/Mi) where wi are mass mixing fractions and Mi are individual molar masses

mass_reciprocal = mmro1/O_molar_mass + mmro2/O2_molar_mass + &
                  (1.0_r8-mmro1-mmro2)/N2_molar_mass

obs_val = pressure / mass_reciprocal / temperature / universal_gas_constant 

end subroutine get_expected_upper_atm_density


!-----------------------------------------------------------------------
!> Given DART state vector and a location,
!> it computes ground GPS vertical total electron content
!> The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_gnd_gps_vtec(state_vector, location, obs_val, istatus)


real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer  :: nAlts, iAlt
real(r8), allocatable :: ALT(:), IDensityS_ie(:) 
real(r8) :: loc_vals(3)
real(r8) :: tec
type(location_type) :: probe

! Given a location and the state vector from one of the ensemble members,
! compute the model-predicted total electron content that would be in the
! integrated column from an instrument looking straight down at the tangent point.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

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
   call interpolate(state_vector, probe, KIND_GEOMETRIC_HEIGHT,    ALT(iAlt), istatus)
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


!-----------------------------------------------------------------------
!> Given DART state vector and a location,
!> compute ground GPS vertical total electron content including an estimate of
!> the contribution from above the model (the 'extrap'olated part)
!> The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_gps_vtec_extrap(state_vector, location, obs_val, istatus)

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


!-----------------------------------------------------------------------
!>

subroutine get_expected_O_N2_ratio(state_vector, location, obs_val, istatus)

real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat
type(location_type) :: loc

real(r8), parameter :: Max_N2_column_density = 1.0E21_r8

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
real(r8), allocatable :: total_number_density(:)
real(r8), allocatable :: O_number_density(:)

integer :: ilayer, nlevels, nilevels
integer :: vstatus(4)
real(r8) :: layerfraction

! First, find the number of levels in the model.
! Then, loop down through the levels to create a top-down vertical profile.
!       As we do that, we accumulate the amount of N2 and O, stopping when
!       the N2 reaches 10^21 M^-2. This will probably mean only using part
!       of the 'last' layer.

if ( .not. module_initialized ) call initialize_module

istatus = 1
obs_val = MISSING_R8

call error_handler(E_MSG, 'get_expected_O_N2_ratio', 'routine not tested', &
           source, revision, revdate, &
           text2='routine in obs_def/obs_def_upper_atm_mod.f90', &
           text3='test and inform the DART development team. Thanks -- Tim.')

if ( .not. module_initialized ) call initialize_module

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

! some variables are defined on interface layers

nilevels = 0
heights = 0.0_r8

FILLINTERFACES : do ilayer = 1,MAXLEVELS

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_vector, loc, KIND_GEOMETRIC_HEIGHT, heights(ilayer),istatus)
   if (istatus /= 0) exit FILLINTERFACES

   nilevels = nilevels + 1

enddo FILLINTERFACES

if (nilevels == 0) return

thickness = 0.0_r8
thickness(1:nilevels-1) = heights(2:nilevels) - heights(1:nilevels-1)

! Some variables are defined on midpoints of the layers

nlevels = 0

FILLMIDPOINTS : do ilayer = 1,MAXLEVELS

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_vector, loc, KIND_PRESSURE, &
                    pressure(ilayer), vstatus(1))

   call interpolate(state_vector, loc, KIND_TEMPERATURE, &
                    temperature(ilayer), vstatus(2))

   call interpolate(state_vector, loc, KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                    O_mmr(ilayer), vstatus(3))

   call interpolate(state_vector, loc, KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                    O2_mmr(ilayer), vstatus(4))

   if (any(vstatus /= 0)) exit FILLMIDPOINTS

   nlevels = nlevels + 1

enddo FILLMIDPOINTS

if (nlevels == 0) return

! Check to make sure we have more interfaces than layers.

if (nilevels /= (nlevels+1)) then
   write(string1,*)'Require there to be 1 more interfaces than midpoints.'
   write(string2,*)'Found ',nilevels,' interface layers.'
   write(string3,*)'Found ',nlevels,' midpoint layers.'
   call error_handler(E_MSG,'get_expected_O_N2_ratio', string1, &
              source, revision, revdate, text2=string2,text3=string3)
   return
endif

! calculate what we can using array notation

allocate(N2_mmr(nlevels), mbar(nlevels), total_number_density(nlevels), &
              O_number_density(nlevels),    N2_number_density(nlevels))

N2_mmr = 1.0_r8 - O_mmr(1:nlevels) - O2_mmr(1:nlevels)
  mbar = 1.0_r8/( O2_mmr(1:nlevels)/O2_molar_mass_kg + &
                   O_mmr(1:nlevels)/ O_molar_mass_kg + &
                  N2_mmr(1:nlevels)/N2_molar_mass_kg )

! O_mmr and N2_mmr defined at midpoints, heights defined at interfaces, so the
! calculated thicknesses apply directly to the O and N2 densities.

total_number_density = pressure(1:nlevels) / (k_constant * temperature(1:nlevels))

 O_number_density =  O_mmr(1:nlevels) * mbar /  O_molar_mass_kg * total_number_density
N2_number_density = N2_mmr(1:nlevels) * mbar / N2_molar_mass_kg * total_number_density

if ( 1 == 2 ) then ! DEBUG BLOCK NOT IN USE
   write(*,*)
   do ilayer = nlevels,1,-1
      write(*,*)'DEBUG level, thickness ...',ilayer, thickness(ilayer), &
            O_number_density(ilayer), N2_number_density(ilayer), &
                 temperature(ilayer), total_number_density(ilayer), &
                 O2_mmr(ilayer), O_mmr(ilayer), N2_mmr(ilayer), mbar(ilayer)
   enddo
   write(*,*)
endif

N2_total = 0.0_r8
 O_total = 0.0_r8

TOPDOWN : do ilayer = nlevels-1, 1, -1

   if (ilayer == 1) then
      write(string1,*)'Integrated all the way down to the surface.'
      write(string2,*)'Still do not have ',Max_N2_column_density,' nitrogen molecules per m^2'
      call error_handler(E_MSG,'get_expected_O_N2_ratio', string1, &
                 source, revision, revdate, text2=string2)
      istatus = 2
      return
   endif

   ! integrate over layer thickness
   O_integrated  =  O_number_density(ilayer) * thickness(ilayer)
   N2_integrated = N2_number_density(ilayer) * thickness(ilayer)

   if ((N2_total+N2_integrated) >= Max_N2_column_density) then
      ! only store part of the final layer so as not to overshoot 10^21 m^-2
      ! Let y2 == N2_total, y = Max_N2_column_density, y1 = N2_total + N2_integrated
      ! the layer fraction is (y - y2)/(y1-y2)
      ! (Max_N2_column_density - N2_total)/(N2_total + N2_integrated - N2_total)
      layerfraction = (Max_N2_column_density - N2_total) / N2_integrated
      N2_total = N2_total + N2_integrated*layerfraction
       O_total =  O_total +  O_integrated*layerfraction
      exit TOPDOWN
   else
      N2_total = N2_total + N2_integrated
       O_total =  O_total +  O_integrated
   endif

enddo TOPDOWN

obs_val = O_total / N2_total
istatus = 0

deallocate(N2_mmr, mbar, total_number_density, O_number_density, N2_number_density)

end subroutine get_expected_O_N2_ratio


!-----------------------------------------------------------------------
!> Given DART state vector and a location,
!> compute ground GPS vertical total electron content including an estimate of
!> the contribution from above the model (the 'extra'polated part)
!> The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_slant_tec(state_vector, location, key, obs_val, istatus)

real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: key              ! into module metadata
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

type(location_type) :: recvr_location
type(location_type) :: trans_location

call error_handler(E_ERR, 'get_expected_slant_tec', 'method not tested', &
           source, revision, revdate, text2='consider this "pre-alpha"')

if ( .not. module_initialized ) call initialize_module

istatus = 1
obs_val = MISSING_R8

call get_slant_tec_metadata(key, recvr_location, trans_location )

end subroutine get_expected_slant_tec


!-----------------------------------------------------------------------
!> read_slant_tec_metadata reads the metadata for slant tec observations
!> from an observation sequence file and fills a local array.

subroutine read_slant_tec_metadata(key, obsID, ifile, fform)

integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
logical           :: is_asciifile
integer           :: ierr
character(len=labellength) :: header
type(location_type) :: recvr_location, trans_location

call error_handler(E_ERR, 'read_slant_tec_metadata', 'method not tested', &
           source, revision, revdate, text2='consider this "pre-alpha"')

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

write(string2,*)'observation #',obsID

! Read the receiver location after confirming we are in sync

if ( is_asciifile ) then
   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_slant_tec_metadata','ASCII receiver header',string2)
else
   read(ifile, iostat=ierr) header
   call  check_iostat(ierr,'read_slant_tec_metadata','binary receiver header',string2)
endif

if (trim(header) /= trim(RECEIVERSTRING)) then
    write(string1,*)"Expected slant tec receiver header ["//trim(RECEIVERSTRING),&
                    &"] in input file, got ["//header//"]"
    call error_handler(E_ERR, 'read_slant_tec_metadata', string1, &
               source, revision, revdate, text2=string2)
endif

recvr_location = read_location(ifile, fform)

! Read the transmitter location after confirming we are in sync

if ( is_asciifile ) then
   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_slant_tec_metadata','ASCII transmitter header',string2)
else
   read(ifile, iostat=ierr) header
   call  check_iostat(ierr,'read_slant_tec_metadata','binary transmitter header',string2)
endif

if (trim(header) /= trim(TRANSMITTERSTRING)) then
    write(string1,*)"Expected slant tec transmitter header ["//trim(TRANSMITTERSTRING),&
                    &"] in input file, got ["//header//"]"
    call error_handler(E_ERR, 'read_slant_tec_metadata', string1, &
               source, revision, revdate, text2=string2)
endif

trans_location = read_location(ifile, fform)

! Store the metadata in module storage and record the new length of the metadata arrays.
call set_slant_tec_metadata(key, recvr_location, trans_location)

end subroutine read_slant_tec_metadata


!-----------------------------------------------------------------------
!> writes the metadata for slant tec observations.

subroutine write_slant_tec_metadata(key, ifile, fform)

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform

logical  :: is_asciifile
type(location_type) :: recvr_location, trans_location

call error_handler(E_ERR, 'write_slant_tec_metadata', 'method not tested', &
           source, revision, revdate, text2='consider this "pre-alpha"')

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

call get_slant_tec_metadata(key, recvr_location, trans_location)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   write(ifile, *) trim(RECEIVERSTRING)
else
   write(ifile   ) trim(RECEIVERSTRING)
endif
call write_location(ifile, recvr_location, fform)

if (is_asciifile) then
   write(ifile, *) trim(TRANSMITTERSTRING)
else
   write(ifile   ) trim(TRANSMITTERSTRING)
endif
call write_location(ifile, trans_location, fform)

end subroutine write_slant_tec_metadata


!-----------------------------------------------------------------------
!>

subroutine interactive_slant_tec_metadata(key)

integer, intent(out) :: key

type(location_type) :: recvr_location
type(location_type) :: trans_location

call error_handler(E_ERR, 'interactive_slant_tec_metadata', 'method not tested', &
           source, revision, revdate, text2='consider this "pre-alpha"')

if ( .not. module_initialized ) call initialize_module

! Prompt for input for the required metadata

write(*,*)'Enter location of receiver'
call interactive_location(recvr_location)

write(*,*)'Enter location of transmitter'
call interactive_location(trans_location)

call set_slant_tec_metadata(key, recvr_location, trans_location)

end subroutine interactive_slant_tec_metadata


!-----------------------------------------------------------------------
!> Common code to increment the current key count, and set the private
!> contents of this observation's auxiliary data.

subroutine set_slant_tec_metadata(key, recvr_location, trans_location )

integer,             intent(out) :: key
type(location_type), intent(in)  :: recvr_location
type(location_type), intent(in)  :: trans_location

if ( .not. module_initialized ) call initialize_module

slant_tec_counter = slant_tec_counter + 1

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(slant_tec_counter,'set_slant_tec_metadata')

key = slant_tec_counter ! now that we know it is legal

slant_tec_metadata(key)%receiver_location    = recvr_location
slant_tec_metadata(key)%transmitter_location = trans_location

end subroutine set_slant_tec_metadata


!-----------------------------------------------------------------------
!> Common code to increment the current key count, and set the private
!> contents of this observation's auxiliary data.

subroutine get_slant_tec_metadata(key, recvr_location, trans_location )


integer,             intent(in) :: key
type(location_type), intent(out)  :: recvr_location
type(location_type), intent(out)  :: trans_location

if ( .not. module_initialized ) call initialize_module

! Make sure the new key is within the length of the metadata arrays.
call key_within_range(key,'get_slant_tec_metadata')

recvr_location = slant_tec_metadata(key)%receiver_location
trans_location = slant_tec_metadata(key)%transmitter_location

end subroutine get_slant_tec_metadata


!-----------------------------------------------------------------------
!>

subroutine check_iostat(istat, routine, varname, msgstring)

integer,          intent(in) :: istat
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: msgstring

if ( istat /= 0 ) then
   write(string1,*)'istat should be 0 but is ',istat,' for ['//varname//']'
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
end if

end subroutine check_iostat


!-----------------------------------------------------------------------
!> Make sure we are addressing within the metadata arrays

subroutine key_within_range(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

! fine -- no problem.
if ((key > 0) .and. (key <= slant_tec_counter)) return

! Bad news. Tell the user.
write(string1, *) 'key (',key,') not within known range ( 1,', slant_tec_counter,')'
call error_handler(E_ERR,routine,string1,source,revision,revdate)

end subroutine key_within_range


!-----------------------------------------------------------------------
!> If the allocatable metadata arrays are not big enough ... try again

subroutine grow_metadata(key, routine)

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

integer :: orglength
type(slant_tec_type), allocatable :: safe_metadata(:)

! fine -- no problem.
if ((key > 0) .and. (key <= MAX_slant_tec_key)) return

orglength         = MAX_slant_tec_key
MAX_slant_tec_key = 2 * orglength

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*MAX_slant_tec_key) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

! News. Tell the user we are increasing storage.
write(string1, *) 'key (',key,') exceeds MAX_slant_tec_key (',orglength,')'
write(string2, *) 'Increasing MAX_slant_tec_key to ',MAX_slant_tec_key
call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

allocate(safe_metadata(orglength))
safe_metadata(:) = slant_tec_metadata(:)

deallocate(slant_tec_metadata)
  allocate(slant_tec_metadata(MAX_slant_tec_key))

slant_tec_metadata(1:orglength) = safe_metadata(:)

deallocate(safe_metadata)

end subroutine grow_metadata



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
real(r8) :: loc_lon, loc_lat
type(location_type) :: loc

real(r8), parameter :: Max_N2_column_density = 1.0E21_r8

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
real(r8), allocatable :: total_number_density(:)
real(r8), allocatable :: O_number_density(:)

real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
integer :: ilayer, nlevels, nilevels
integer :: vstatus(4)
real(r8) :: layerfraction

if ( .not. module_initialized ) call initialize_module

istatus = 1
obs_val = MISSING_R8

call error_handler(E_ERR, 'get_expected_O_N2_ratio', 'routine not tested', &
           source, revision, revdate, &
           text2='routine in obs_def/obs_def_upper_atm_mod.f90', &
           text3='test and inform the DART development team. Thanks -- Tim.')

if ( .not. module_initialized ) call initialize_module

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

! some variables are defined on interface layers

nilevels = 0
heights = 0.0_r8

FILLINTERFACES : do ilayer = 1,MAXLEVELS

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_vector, loc, KIND_GEOMETRIC_HEIGHT, heights(ilayer),istatus)
   if (istatus /= 0) exit FILLINTERFACES

   nilevels = nilevels + 1

enddo FILLINTERFACES

if (nilevels == 0) return

thickness = 0.0_r8
thickness(1:nilevels-1) = heights(2:nilevels) - heights(1:nilevels-1)

! Some variables are defined on midpoints of the layers

nlevels = 0

FILLMIDPOINTS : do ilayer = 1,MAXLEVELS

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_vector, loc, KIND_PRESSURE, &
                    pressure(ilayer), vstatus(1))

   call interpolate(state_vector, loc, KIND_TEMPERATURE, &
                    temperature(ilayer), vstatus(2))

   call interpolate(state_vector, loc, KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                    O_mmr(ilayer), vstatus(3))

   call interpolate(state_vector, loc, KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                    O2_mmr(ilayer), vstatus(4))

   if (any(vstatus /= 0)) exit FILLMIDPOINTS

   nlevels = nlevels + 1

enddo FILLMIDPOINTS

if (nlevels == 0) return

! Check to make sure we have more interfaces than layers.

if (nilevels /= (nlevels+1)) then
   write(string1,*)'Require there to be 1 more interfaces than midpoints.'
   write(string2,*)'Found ',nilevels,' interface layers.'
   write(string3,*)'Found ',nlevels,' midpoint layers.'
   call error_handler(E_MSG,'get_expected_O_N2_ratio', string1, &
              source, revision, revdate, text2=string2,text3=string3)
   return
endif

! calculate what we can using array notation

allocate(N2_mmr(nlevels), mbar(nlevels), total_number_density(nlevels), &
              O_number_density(nlevels),    N2_number_density(nlevels))

N2_mmr = 1.0_r8 - O_mmr(1:nlevels) - O2_mmr(1:nlevels)
  mbar = 1.0_r8/( O2_mmr(1:nlevels)/O2_molar_mass + &
                   O_mmr(1:nlevels)/ O_molar_mass + &
                  N2_mmr(1:nlevels)/N2_molar_mass )

! O_mmr and N2_mmr defined at midpoints, heights defined at interfaces, so the
! calculated thicknesses apply directly to the O and N2 densities.

total_number_density = pressure(1:nlevels) / (k_constant * temperature(1:nlevels))

 O_number_density =  O_mmr(1:nlevels) * mbar /  O_molar_mass * total_number_density 
N2_number_density = N2_mmr(1:nlevels) * mbar / N2_molar_mass * total_number_density

if ( 1 == 2 ) then ! DEBUG BLOCK NOT IN USE
   write(*,*)
   do ilayer = nlevels,1,-1
      write(*,*)'DEBUG level, thickness ...',ilayer, thickness(ilayer), &
            O_number_density(ilayer), N2_number_density(ilayer), &
                 temperature(ilayer), total_number_density(ilayer), &
                 O2_mmr(ilayer), O_mmr(ilayer), N2_mmr(ilayer), mbar(ilayer)
   enddo
   write(*,*)
endif

N2_total = 0.0_r8
 O_total = 0.0_r8

TOPDOWN : do ilayer = nlevels,1,-1

   if (ilayer == 1) then
      write(string1,*)'Integrated all the way down to the surface.'
      write(string2,*)'Still do not have ',Max_N2_column_density,' nitrogen molecules per m^2'
      call error_handler(E_MSG,'get_expected_O_N2_ratio', string1, &
                 source, revision, revdate, text2=string2)
      istatus = 2
      return
   endif

   ! integrate over layer thickness
   O_integrated  =  O_number_density(ilayer) * thickness(ilayer)
   N2_integrated = N2_number_density(ilayer) * thickness(ilayer)

   if ((N2_total+N2_integrated) >= Max_N2_column_density) then
      ! only store part of the final layer so as not to overshoot 10^21 m^-2
      ! Let y2 == N2_total, y = Max_N2_column_density, y1 = N2_total + N2_integrated
      ! the layer fraction is (y - y2)/(y1-y2) 
      ! (Max_N2_column_density - N2_total)/(N2_total + N2_integrated - N2_total)
      layerfraction = (Max_N2_column_density - N2_total) / N2_integrated
      N2_total = N2_total + N2_integrated*layerfraction
       O_total =  O_total +  O_integrated*layerfraction
      exit TOPDOWN
   else
      N2_total = N2_total + N2_integrated
       O_total =  O_total +  O_integrated
   endif

enddo TOPDOWN

obs_val = O_total / N2_total
istatus = 0

deallocate(N2_mmr, mbar, total_number_density, O_number_density, N2_number_density)

end subroutine get_expected_O_N2_ratio


end module obs_def_upper_atm_mod
! END DART PREPROCESS MODULE CODE

