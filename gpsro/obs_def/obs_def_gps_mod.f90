! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! Observation-specific routines for 'Refractivity' and 'Excess Phase' 
! values obtained from low-orbit satellite occultation instruments. 
! The fundamental observation from these instruments is the delay time
! in arrival of a signal as it passes through the atmosphere as a pair of
! satellites rise and set relative to each other.  The delay time is caused
! by the bending of the signal as it goes through the atmosphere (it
! traverses a longer path); the amount of bending/delay is a function of 
! the pressure, temperature, and moisture in the air.
!
! Refractivity is a Level 2 product that can be read directly from the
! input files and is a unitless quantity.  Excess phase is a line integral
! computed along the straight distance between the sending and receiving 
! satellites, with units of meters (the distance is multiplied by the
! local refractivity and not divided out at the end).
!
! To assimilate either of these observation types the model must be able
! to return values for the Temperature, Pressure, and Specific Humidity
! at any given location.  The conversion to refractivity or excess phase
! is done in the code in this file.  Temperature needs to be in kelvin,
! pressure in pascals, and specific humidity in kg/kg.
!
! In the DART obs_seq files, the refractivity does not need any additional
! metadata information; the refractivity is computed at the tangent point
! of the ray with the simplifying assumption that all bending takes place
! at this point.  However, the excess phase computation requires the
! local refractivity be computed at multiple points along the ray and
! summed, so those obervations need the geometry of the ray, the step-size
! for the integral, and the height at which to stop the integration.
! There is a default maximum number of excess phase observations across
! all obs_seq files being used (might be a large number for the obs_diag
! program) which can be increased by setting a larger value in the namelist.


! BEGIN DART PREPROCESS KIND LIST
! TEMPERATURE,          KIND_TEMPERATURE,              COMMON_CODE
! SPECIFIC_HUMIDITY,    KIND_SPECIFIC_HUMIDITY,        COMMON_CODE
! PRESSURE,             KIND_PRESSURE,                 COMMON_CODE
! GPS_RO_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY
! GPS_RO_EXCESSPHASE,   KIND_OCCULTATION_EXCESSPHASE
! END DART PREPROCESS KIND LIST

! alternatives, if someone wanted to separate the obs by source.
! the converter program that creates obs_seq files would need
! to use these same types when creating the obs (the scripts
! know which satellite data is being used, but the converter
! executable has no knowledge at this point about the source).
 !! COSMIC_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !!  CHAMP_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !!  GRACE_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !!   SACC_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !!  CNOFS_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !!    TSX_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !! METOPA_REFRACTIVITY,  KIND_OCCULTATION_REFRACTIVITY !!
 !!  COSMIC_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!
 !!   CHAMP_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!
 !!   GRACE_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!
 !!    SACC_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!
 !!   CNOFS_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!
 !!     TSX_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!
 !!  METOPA_EXCESSPHASE,  KIND_OCCULTATION_EXCESSPHASE  !!

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_gps_mod, only : get_expected_excess_phase, &
!                              interactive_excess_phase,  &
!                              read_excess_phase,         &
!                              write_excess_phase,        &
!                              get_expected_occult_refract
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(GPS_RO_EXCESSPHASE)
!            call get_expected_excess_phase(state, location, obs_def%key, obs_val, istatus)
!         case(GPS_RO_REFRACTIVITY)
!            call get_expected_occult_refract(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(GPS_RO_EXCESSPHASE)
!            call read_excess_phase(obs_def%key, ifile, fform)
!         case(GPS_RO_REFRACTIVITY)
!            continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(GPS_RO_EXCESSPHASE)
!            call write_excess_phase(obs_def%key, ifile, fform)
!         case(GPS_RO_REFRACTIVITY)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(GPS_RO_EXCESSPHASE)
!            call interactive_excess_phase(obs_def%key)
!         case(GPS_RO_REFRACTIVITY)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_gps_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             file_exist, open_file, close_file, nmlfileunit, &
                             check_namelist_read, find_namelist_in_file, &
                             do_output, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, &
                             write_location, read_location, vert_is_height, &
                             VERTISHEIGHT
use time_manager_mod, only : time_type, read_time, write_time, &
                             set_time, set_time_missing, interactive_time
use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY,     &
                             KIND_PRESSURE, KIND_OCCULTATION_REFRACTIVITY, &
                             KIND_OCCULTATION_EXCESSPHASE
   

implicit none
private

public :: get_expected_excess_phase,   &
          read_excess_phase,           &
          write_excess_phase,          &
          interactive_excess_phase,    &
          set_excess_phase,            &
          get_excess_phase,            &
          get_expected_occult_refract, &
          cartesian2geodetic,          &
          geodetic2cartesian,          &
          point_refractivity

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$"
character(len=128), parameter :: &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len=129) :: msgstring
integer            :: keycount     ! module global storage

! Storage for the special information required for excess phase
! computation.  It needs to integrate, at run time, along the
! occulation ray (assumed to be a straight line), so the direction,
! step-length, and max height need to be stored with the obs.
!
! The point computation for refractivity doesn't need any additional
! metadata beyond the standard obs information.

integer :: max_xs_phase_obs = 100000
type xs_phase_type                       ! all distances here in kilometers, 
   private 
   real(r8)         :: ray_direction(3)  ! cartesian unit vector, x/y/z
   real(r8)         :: radius_curvature  ! rfict read from input profile
   real(r8)         :: step_size         ! integration step length (km, was m)
   real(r8)         :: ray_stop_height   ! when integration endpts exceed, quit
end type xs_phase_type

type(xs_phase_type), allocatable :: xs_data(:)

namelist /obs_def_gps_nml/ max_xs_phase_obs

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine initialize_module()

! Initialize xs phase private key counter, allocate xs phase metadata arrays.

integer :: rc, iunit

call register_module(source, revision, revdate)
module_initialized = .true.

! global count of all xs phase observations from any input file
keycount = 0

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_gps_nml", iunit)
read(iunit, nml = obs_def_gps_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_gps_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_gps_nml)
if (do_nml_term()) write(     *     , nml=obs_def_gps_nml)

if (max_xs_phase_obs > 0) then
   allocate(xs_data(max_xs_phase_obs), stat = rc)
   if (rc /= 0) then
      write(msgstring, *) 'initial allocation failed for xs phase observation data,', &
                          'itemcount = ', max_xs_phase_obs
      call error_handler(E_ERR,'initialize_module', msgstring, &
                         source, revision, revdate)
   endif
endif

end subroutine initialize_module

!------------------------------------------------------------------------------

subroutine set_excess_phase(xskey, rayx, rayy, rayz, &
                            radius_curvature, step_size, ray_stop_height)
 
! increment key and set all private data for this observation

integer,  intent(out) :: xskey
real(r8), intent(in)  :: rayx, rayy, rayz
real(r8), intent(in)  :: radius_curvature, step_size, ray_stop_height

if ( .not. module_initialized ) call initialize_module

keycount = keycount + 1
xskey = keycount

if(xskey > max_xs_phase_obs) then
   write(msgstring, *) 'key (',xskey,') exceeds max_xs_phase_obs (',max_xs_phase_obs,')'
   call error_handler(E_ERR,'read_gpsro_ref', msgstring, &
                      source, revision, revdate)
endif

xs_data(xskey)%ray_direction(1) = rayx
xs_data(xskey)%ray_direction(2) = rayy
xs_data(xskey)%ray_direction(3) = rayz

xs_data(xskey)%radius_curvature = radius_curvature
xs_data(xskey)%step_size        = step_size
xs_data(xskey)%ray_stop_height  = ray_stop_height

end subroutine set_excess_phase

!------------------------------------------------------------------------------

subroutine get_excess_phase(xskey, height, longitude, latitude, &
                            radius_curvature, step_size, ray_stop_height)
 
! return all private data for this observation

integer,  intent(in)  :: xskey
real(r8), intent(out) :: height, longitude, latitude
real(r8), intent(out) :: radius_curvature, step_size, ray_stop_height


if ( .not. module_initialized ) call initialize_module

if (xskey < 1 .or. xskey > keycount) then
   write(msgstring, *) 'key (',xskey,') out of valid range (1<=key<=',keycount,')'
   call error_handler(E_ERR,'get_gpsro_ref', msgstring, &
                      source, revision, revdate)
endif

height    = xs_data(xskey)%ray_direction(1)
latitude  = xs_data(xskey)%ray_direction(2)
longitude = xs_data(xskey)%ray_direction(3)

radius_curvature = xs_data(xskey)%radius_curvature
step_size        = xs_data(xskey)%step_size
ray_stop_height  = xs_data(xskey)%ray_stop_height

end subroutine get_excess_phase



 subroutine write_excess_phase(xskey, ifile, fform)
!------------------------------------------------------------------------------
!

integer,          intent(in)           :: xskey, ifile
character(len=*), intent(in), optional :: fform

integer  :: ii

if ( .not. module_initialized ) call initialize_module

! Write out the obs_def key for this observation
if (ascii_file_format(fform)) then
   write(ifile, *) xs_data(xskey)%radius_curvature, xs_data(xskey)%step_size, &
                   xs_data(xskey)%ray_stop_height, &
                  (xs_data(xskey)%ray_direction(ii), ii=1, 3)
else
   write(ifile) xs_data(xskey)%radius_curvature, xs_data(xskey)%step_size, &
                xs_data(xskey)%ray_stop_height, &
               (xs_data(xskey)%ray_direction(ii), ii=1, 3)
endif

end subroutine write_excess_phase



 subroutine read_excess_phase(xskey, ifile, fform)
!------------------------------------------------------------------------------
!
! Every excess phase observation has its own (metadata) xskey.
! When you read multiple observation sequence files it is necessary 
! to track the total number of metadata xskeys read not just the number 
! in the current file.
! 

integer,          intent(out)          :: xskey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

integer :: keyin    ! the metadata key in the input obs sequence

real(r8) :: nx, ny, nz, rfict0, ds, htop

if ( .not. module_initialized ) call initialize_module

if (ascii_file_format(fform)) then
   read(ifile, *) rfict0, ds, htop, nx, ny, nz
else
   read(ifile)    rfict0, ds, htop, nx, ny, nz
endif


! increment key and set all private data for this observation
call set_excess_phase(xskey, nx, ny, nz, rfict0, ds, htop)

end subroutine read_excess_phase


subroutine interactive_excess_phase(gpskey)
!----------------------------------------------------------------------
!
! Interactively prompt for the info needed to create a gps excess phase 
! observation.  Increments the key number and returns it.

 integer, intent(out) :: gpskey

real(r8) :: longitude, latitude
real(r8) :: radius_curvature, step_size, ray_stop_height
real(r8) :: nx, ny, nz, azimuth


if ( .not. module_initialized ) call initialize_module

!Now interactively obtain reflectivity type information
! valid choices are local or non-local

write(*, *)
write(*, *) 'Beginning to inquire information on excess phase type.'
write(*, *)

10 continue
   write(*, *)
   write(*, *) 'Enter longitude, latitude of tangent point'
   write(*, *) ' (lon 0,360 deg, lat -90,90 deg)'
   write(*, *)
   read(*,*) longitude, latitude
   if ((longitude <   0.0_r8 .or. longitude > 360.0_r8) .or. &
       (latitude  < -90.0_r8 .or. latitude  >  90.0_r8)) then
      write(*, *) 'error: values out of range.  enter new values'
      goto 10
   endif

20 continue
   write(*, *)
   write(*, *) 'Enter angle of occultation plane from north'
   write(*, *) ' (between -180,180 deg)'
   write(*, *)
   read(*,*) azimuth
   if (azimuth < -180.0_r8 .or. azimuth > 180.0_r8) then
      write(*, *) 'error: value out of range.  enter new value'
      goto 20
   endif

   ! FIXME:  i have no idea what valid values are for the following
   !  item so i cannot add any error checking or guidance for the user.
   !  i'm assuming it must be positive but beyond that... ?

30 continue
   write(*, *)
   write(*, *) 'Enter local curvature radius (in meters)'
   write(*, *)
   read(*,*) radius_curvature
   if (radius_curvature < 0.0_r8) then
      write(*, *) 'error: value out of range.  enter new value'
      goto 30
   endif

40 continue
   write(*, *)
   write(*, *) 'Enter integration step size along tangent ray (in km)'
   write(*, *) ' (5 km is a common default)'
   write(*, *)
   read(*,*) step_size
   if (step_size <= 0.0_r8) then
      write(*, *) 'error: value out of range.  enter new value'
      goto 40
   endif

50 continue
   write(*, *)
   write(*, *) 'Enter max height where ray integration stops (in km)'
   write(*, *) ' (between 10 to 15 km is a common default.'
   write(*, *) '  profiles include values up to 60 km)'
   write(*, *)
   read(*,*) ray_stop_height
   if ((ray_stop_height <= 0.0_r8) .or. (ray_stop_height > 60.0_r8)) then
      write(*, *) 'error: value out of range.  enter new value'
      goto 50
   endif

! convert lat/lon and azimuth plane into a unit vector along the ray
call tanvec01(longitude, latitude, azimuth, nx, ny, nz)

! increment key and set all private data for this observation
! (the lon, lat, and observation height will be entered as the obs location
! which is prompted for as part of the default obs information.)
call set_excess_phase(gpskey, nx, ny, nz, radius_curvature, step_size, ray_stop_height)

write(*, *)
write(*, *) 'End of specialized section for gps observation data.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_excess_phase

!------------------------------------------------------------------------------
!
! Purpose: Calculate GPS RO local refractivity or non_local (integrated) 
!          refractivity (excess phase, Sergey Sokolovskiy et al., 2005)
!------------------------------------------------------------------------------
!
! inputs:
!    state_vector:    DART state vector
!
! output parameters:
!    ro_ref: modeled local refractivity (N-1)*1.0e6 or non_local 
!            refractivity (excess phase, m)
!            (according to the input data parameter subset)
!    istatus:  =0 normal; =1 outside of domain.
!------------------------------------------------------------------------------
!  Author: Hui Liu 
!  Version 1.1: June 15, 2004: Initial version CAM
!
!  Version 1.2: July 29, 2005: revised for new obs_def and WRF
!------------------------------------------------------------------------------
subroutine get_expected_excess_phase(state_vector, location, gpskey, ro_ref, istatus)
 real(r8),            intent(in)  :: state_vector(:)
 type(location_type), intent(in)  :: location
 integer,             intent(in)  :: gpskey
 real(r8),            intent(out) :: ro_ref
 integer,             intent(out) :: istatus

! local variables
real(r8) :: nx, ny, nz       ! unit tangent direction of ray at perigee
real(r8) :: xo, yo, zo       ! perigee location in Cartesian coordinate

real(r8) :: ref_perigee, ref00, ref1, ref2, dist_to_perigee
real(r8) :: phase
real(r8) :: xx, yy, zz, height1, lat1, lon1, delta_phase1, delta_phase2

integer  :: iter, istatus0
real(r8) :: lon, lat, height, obsloc(3)

if ( .not. module_initialized ) call initialize_module

if ( .not. vert_is_height(location)) then
   write(msgstring, *) 'vertical location must be height; gps obs key ', gpskey
   call error_handler(E_ERR,'get_expected_gpsro_ref', msgstring, &
                      source, revision, revdate)
endif

obsloc   = get_location(location)


lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90
height   = obsloc(3)                       ! (m)

! calculate refractivity at perigee

call point_refractivity(state_vector, location, height, lat, lon, ref_perigee, istatus0)
! if istatus > 0, the interpolation failed and we should return failure now.
if(istatus0 > 0) then
   istatus = istatus0
   ro_ref = missing_r8
   return
endif

! non_local refractivity(excess phase delay)

! Initialization
phase = 0.0_r8  
dist_to_perigee =  0.0_r8   ! distance to perigee from a point of the ray

nx = xs_data(gpskey)%ray_direction(1)
ny = xs_data(gpskey)%ray_direction(2)
nz = xs_data(gpskey)%ray_direction(3)

! convert location of the perigee from geodetic to Cartesian coordinate

call geodetic2cartesian (height, lat, lon, xo, yo, zo, xs_data(gpskey)%radius_curvature )

! currently, use a straight line passing the perigee point as ray model.
! later, more sophisticated ray models can be used.
!
! Start the horizontal integrate of the model refractivity along a 
! straight line path in cartesian coordinate
!
! (x-xo)/a = (y-yo)/b = (z-zo)/c,  (a,b,c) is the line direction

ref1 = ref_perigee
ref2 = ref_perigee

iter = 0
do 

   iter = iter + 1
   dist_to_perigee = dist_to_perigee + xs_data(gpskey)%step_size
   
   !  integrate to one direction of the ray for one step
   xx = xo + dist_to_perigee * nx
   yy = yo + dist_to_perigee * ny
   zz = zo + dist_to_perigee * nz
  
   ! convert the location of the point to geodetic coordinates 
   ! height(m), lat, lon(deg)
   
   call cartesian2geodetic(xx, yy, zz, height1, lat1, lon1, xs_data(gpskey)%radius_curvature )  
   if (height1 >= xs_data(gpskey)%ray_stop_height) exit
   
   ! get the refractivity at this ray point(ref00)
   call point_refractivity(state_vector, location, height1, lat1, lon1, ref00, istatus0)
   ! when any point of the ray is problematic, return failure
   if(istatus0 > 0) then
     istatus = istatus0
     ro_ref = missing_r8
     return
   endif
   
   ! get the excess phase due to this ray interval
   delta_phase1 = (ref1 + ref00) * xs_data(gpskey)%step_size * 0.5_r8
   
   ! save the refractivity for integration of next ray interval
   ref1 = ref00
   
   ! integrate to the other direction of the ray
   xx = xo - dist_to_perigee * nx
   yy = yo - dist_to_perigee * ny 
   zz = zo - dist_to_perigee * nz
  
   call cartesian2geodetic (xx, yy, zz, height1, lat1, lon1, xs_data(gpskey)%radius_curvature )  
   
   ! get the refractivity at this ray point(ref00)
   call point_refractivity(state_vector, location, height1, lat1, lon1, ref00, istatus0)
   ! when any point of the ray is problematic, return failure
   if(istatus0 > 0) then
     istatus = istatus0
     ro_ref = missing_r8
     return
   endif
   
   ! get the excess phase due to this ray interval
   delta_phase2 = (ref2 + ref00) * xs_data(gpskey)%step_size * 0.5_r8
   
   ! save the refractivity for integration of next ray interval
   ref2 = ref00
   
   phase = phase + delta_phase1 + delta_phase2
   ! print*, 'phase= ',  phase, delta_phase1, delta_phase2

end do

! finish the integration of the excess phase along the ray

ro_ref = phase    ! in m

! print*, 'xx = ', lon, lat, height, ro_ref

! if the original height was too high, for example.  do not return a
! negative or 0 excess phase or refractivity.
if (ro_ref == missing_r8 .or. ro_ref <= 0.0_r8) then
   istatus = 5
   ro_ref = missing_r8
   return
endif

! ended ok, return non-local excess phase accumulated value
istatus = 0

end subroutine get_expected_excess_phase


!------------------------------------------------------------------------------
!
! Purpose: Calculate GPS RO local refractivity
!------------------------------------------------------------------------------
!
! inputs:
!    state_vector:    DART state vector
!
! output parameters:
!    ro_ref: modeled local refractivity (N-1)*1.0e6 
!    istatus:  =0 normal; =1 outside of domain or other error.
!------------------------------------------------------------------------------
!  Author: Hui Liu 
!  Version 1.1: June 15, 2004: Initial version CAM
!
!  Version 1.2: July 29, 2005: revised for new obs_def and WRF
! 
!  updated, nancy collins  26 july 2011
!  separated code for local/nonlocal operators.  local operator needs no
!  additional metadata in obs_seq files.
!------------------------------------------------------------------------------
subroutine get_expected_occult_refract(state_vector, location, ro_ref, istatus)
 real(r8),            intent(in)  :: state_vector(:)
 type(location_type), intent(in)  :: location
 real(r8),            intent(out) :: ro_ref
 integer,             intent(out) :: istatus

! local variables
integer  :: istatus0
real(r8) :: ref_perigee, lon, lat, height, obsloc(3)

if ( .not. module_initialized ) call initialize_module

if ( .not. vert_is_height(location)) then
   write(msgstring, *) 'vertical location for gps obs must be height'
   call error_handler(E_ERR,'get_expected_occult_refract', msgstring, &
                      source, revision, revdate)
endif

obsloc   = get_location(location)

lon      = obsloc(1)                       ! degree: 0 to 360
lat      = obsloc(2)                       ! degree: -90 to 90
height   = obsloc(3)                       ! (m)

! calculate refractivity at perigee

call point_refractivity(state_vector, location, height, lat, lon, ref_perigee, istatus0)
! if istatus > 0, the interpolation failed and we should return failure now.
if(istatus0 > 0) then
   istatus = istatus0
   ro_ref = missing_r8
   return
endif

! if the original height was too high, for example.  do not return a
! negative or 0 excess phase or refractivity.
if (ro_ref == missing_r8 .or. ro_ref <= 0.0_r8) then
   istatus = 5
   ro_ref = missing_r8
   return
endif

! local refractivity - fix units and return.
ro_ref = ref_perigee * 1.0e6      ! in (N-1)*1.0e6

! ended ok, return local refractivity 
istatus = 0

end subroutine get_expected_occult_refract



!------------------------------------------------------------------------------
!
! Calculate local refractivity at any GPS ray point (height, lat, lon)
!
! inputs:
!    height, lat, lon:  GPS observation location (units: m, degree)
!
! output:
!    ref00: modeled local refractivity at ray point(unit: N-1, ~1.0e-4 to e-6)
!
!------------------------------------------------------------------------------
subroutine point_refractivity(state_vector, location, height, lat, lon, ref00, istatus0)
 real(r8), intent(in)  :: state_vector(:)
 real(r8), intent(in)  :: lon, lat, height
 real(r8), intent(out) :: ref00
 integer,  intent(out) :: istatus0

real(r8), parameter::  rd = 287.05_r8, rv = 461.51_r8, c1 = 77.6d-6 , &
                       c2 = 3.73d-1,  rdorv = rd/rv
real(r8) :: lon2, t, q, p, tv, ew
type(location_type) :: location, location2
integer :: which_vert

if ( .not. module_initialized ) call initialize_module

! for integration of GPS ray path beyond the wraparound point
lon2 = lon
if(lon > 360.0_r8 ) lon2 = lon - 360.0_r8
if(lon <   0.0_r8 ) lon2 = lon + 360.0_r8

which_vert = VERTISHEIGHT
location2 = set_location(lon2, lat, height,  which_vert)

! set return values assuming failure, so we can simply return if any
! of the interpolation calls below fail.
istatus0 = 3
ref00 = missing_r8

call interpolate(state_vector, location2,  KIND_TEMPERATURE,       t, istatus0)
if (istatus0 > 0) return
call interpolate(state_vector, location2,  KIND_SPECIFIC_HUMIDITY, q, istatus0)
if (istatus0 > 0) return
call interpolate(state_vector, location2,  KIND_PRESSURE,          p, istatus0)
if (istatus0 > 0) return

!  required variable units for calculation of GPS refractivity
!   t :  Kelvin, from top to bottom
!   q :  kg/kg, from top to bottom
!   p :  mb

p     = p * 0.01_r8      ! to mb

tv    = t * (1.0_r8+(rv/rd - 1.0_r8)*q)         ! virtual temperature
ew    = q * p/(rdorv + (1.0_r8-rdorv)*q )
ref00 = c1*p/t + c2*ew/(t**2)              ! (N-1)

! now we have succeeded, set istatus to good
istatus0 = 0

end subroutine point_refractivity


!------------------------------------------------------------------------------
!
!  Converts geodetical coordinates to cartesian with a reference sphere
!------------------------------------------------------------------------------
!  input parameters:
!   s - geodetical coordinates
!        (height (m), latitude (degree), longitude (degree))
!                     -90 to 90           0 to 360
!  output parameters:
!   x - cartesian coordinates (m) connected with the earth(x, y, z-coordinate)
!------------------------------------------------------------------------------
subroutine geodetic2cartesian (s1, s2, s3, x1, x2, x3, rfict0) 
 real(r8), intent(in)  :: s1, s2, s3, rfict0    ! units: m
 real(r8), intent(out) :: x1, x2 ,x3

real(r8) :: g3, g4

if ( .not. module_initialized ) call initialize_module

g3 = s1 + rfict0
g4 = g3 * cos(s2*DEG2RAD) 
x1 = g4 * cos(s3*DEG2RAD)
x2 = g4 * sin(s3*DEG2RAD)
x3 = g3 * sin(s2*DEG2RAD)

end subroutine geodetic2cartesian


!------------------------------------------------------------------------------
!
!  Converts cartesian coordinates to geodetical.
!
!   input parameters:
!        x - cartesian coordinates (x, y, z-coordinate, unit: m)
!
!   output parameters:
!        s - geodetical coordinates
!            (height (m), latitude (deg), longitude (deg))
!                          -90 to 90         0 to 360
!------------------------------------------------------------------------------
subroutine cartesian2geodetic (x1, x2, x3, s1, s2, s3, rfict0)
 real(r8), intent(in)  :: x1, x2, x3, rfict0
 real(r8), intent(out) :: s1, s2, s3

real(r8), parameter :: crcl  = 2.0_r8 * PI, &
                       crcl2 = 4.0_r8 * PI

real(r8) :: rho, sphi, azmth

if ( .not. module_initialized ) call initialize_module

rho   = sqrt (x1**2 + x2**2 + x3**2 ) 
sphi  = x3/rho
s1    = rho - rfict0
s2    = asin (sphi) 
azmth = atan2 (x2, x1)
s3    = mod((azmth + crcl2), crcl)

s2    = s2 * RAD2DEG
s3    = s3 * RAD2DEG

end  subroutine cartesian2geodetic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   tanvec01 - subroutine that computes the unit vector tangent of the
!              ray at the perigee.
!
!    lon0  - longitude of the tangent point
!    lat0  - latitude of the tangent point
!    azim0 - angle between occultation plane from north
!    uz    - x component of tangent vector at tangent point of array
!    uy    - y component of tangent vector at tangent point of array
!    uz    - z component of tangent vector at tangent point of array
!
!     created Hui Liu, NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tanvec01(lon0, lat0, azim0, ux, uy, uz)
 real(r8), intent(in)  :: lon0, lat0, azim0
 real(r8), intent(out) :: ux, uy, uz

real(r8) :: zz0(3), lon, lat, azim, rtp(3), rnm(3), rno(3), uon(3), vlen0

zz0(1) = 0.0_r8  ;  zz0(2) = 0.0_r8  ;  zz0(3) = 1.0_r8
lon = lon0 * deg2rad  ;  lat = lat0 * deg2rad  ;  azim = azim0 * deg2rad

rtp(1) = cos(lat) * cos(lon)
rtp(2) = cos(lat) * sin(lon)
rtp(3) = sin(lat)

!  compute unit vector normal to merdion plane through tangent point
call vprod(rtp, zz0, rnm)
vlen0 = sqrt(rnm(1)*rnm(1) + rnm(2)*rnm(2) + rnm(3)*rnm(3))
rnm(:) = rnm(:) / vlen0

!  compute unit vector toward north from perigee point
call vprod(rnm, rtp, rno)
vlen0 = sqrt(rno(1)*rno(1) + rno(2)*rno(2) + rno(3)*rno(3))
rno(:) = rno(:) / vlen0

!  rotate the vector rno around rtp for a single azim to get tangent vector
call spin(rno, rtp, azim, uon)
vlen0 = sqrt(uon(1)*uon(1) + uon(2)*uon(2) + uon(3)*uon(3))
ux = uon(1) / vlen0  ;  uy = uon(2) / vlen0  ;  uz = uon(3) / vlen0

return
end subroutine tanvec01

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   spin - subroutine that rotates vector v1 around vs clockwise by 
!          by a specified angle.
!
!    v1 - vector to rotate
!    vs - vector to rotate about
!     a - angle to rotate v1 around
!    v2 - output vector after rotation
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spin(v1, vs, a, v2)
 real(r8), intent(in)  :: v1(3), vs(3), a
 real(r8), intent(out) :: v2(3)

real(r8) :: vsabs, vsn(3), a1, a2, a3, s(3,3) 

! Calculation of the unit vector for the rotation
vsabs  = sqrt(vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))
vsn(:) = vs(:) / vsabs

! Calculation the rotation matrix
a1 = cos(a)  ;  a2 = 1.0_r8 - a1  ;  a3 = sin(a)
s(1,1) = a2 * vsn(1) * vsn(1) + a1
s(1,2) = a2 * vsn(1) * vsn(2) - a3 * vsn(3)
s(1,3) = a2 * vsn(1) * vsn(3) + a3 * vsn(2)
s(2,1) = a2 * vsn(2) * vsn(1) + a3 * vsn(3)
s(2,2) = a2 * vsn(2) * vsn(2) + a1
s(2,3) = a2 * vsn(2) * vsn(3) - a3 * vsn(1)
s(3,1) = a2 * vsn(3) * vsn(1) - a3 * vsn(2)
s(3,2) = a2 * vsn(3) * vsn(2) + a3 * vsn(1)
s(3,3) = a2 * vsn(3) * vsn(3) + a1

!  Compute the rotated vector
v2(:) = s(:,1) * v1(1) + s(:,2) * v1(2) + s(:,3) * v1(3)

end subroutine spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   vprod - subroutine that computes the vector product of two vectors
!
!    x - first vector
!    y - second vector
!    z - vector product
!
!     created Hui Liu NCAR/MMM
!     modified June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vprod(x, y, z)
 real(r8), intent(in)  :: x(3), y(3)
 real(r8), intent(out) :: z(3)

z(1) = x(2)*y(3) - x(3)*y(2)
z(2) = x(3)*y(1) - x(1)*y(3)
z(3) = x(1)*y(2) - x(2)*y(1)

end subroutine vprod


end module obs_def_gps_mod

! END DART PREPROCESS MODULE CODE

