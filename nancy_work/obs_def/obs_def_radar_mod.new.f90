! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------------
!  DART radar observation module.  Does ...
!
! < FIXME: need input from david here.  want about 10 lines of general info >
!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS KIND LIST
! DOPPLER_RADIAL_VELOCITY, KIND_VELOCITY
! RADAR_REFLECTIVITY, KIND_RADAR_REFLECTIVITY
! END DART PREPROCESS KIND LIST
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_radar_mod, only : write_radial_vel, read_radial_vel, &
!                            interactive_radial_vel, get_expected_radial_vel, &
!                            write_radar_ref, read_radar_ref, &
!                            interactive_radar_ref, get_expected_radar_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(DOPPLER_RADIAL_VELOCITY)
!            call get_expected_radial_vel(state, location, obs_def%key, obs_val, istatus)
!         case(RADAR_REFLECTIVITY)
!            call get_expected_radar_ref(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call read_radial_vel(obs_def%key, ifile, fileformat)
!      case(RADAR_REFLECTIVITY)
!         call read_radar_ref(obs_val, obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call write_radial_vel(obs_def%key, ifile, fileformat)
!      case(RADAR_REFLECTIVITY)
!         call write_radar_ref(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(DOPPLER_RADIAL_VELOCITY)
!         call interactive_radial_vel(obs_def%key)
!      case(RADAR_REFLECTIVITY)
!         call interactive_radar_ref(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_radar_mod

use        types_mod, only : r8, missing_r8, PI, DEG2RAD
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file, &
                             nmlfileunit, do_output
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_TEMPERATURE, KIND_VERTICAL_VELOCITY, &
                             KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, &
                             KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO

implicit none
private

public :: write_radar_ref, read_radar_ref, interactive_radar_ref, &
          get_expected_radar_ref, get_obs_def_radar_ref, &
          write_radial_vel, read_radial_vel, interactive_radial_vel, &
          get_expected_radial_vel, get_obs_def_radial_vel

! FIXME: set_radial_vel removed from public list.  could be added back
! if needed by any other programs.

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical :: module_initialized = .false.

! Derived type for radial velocity.  Contains auxiliary information stored
! with each obs of this type; used to compute the forward operator.
type radial_vel_type
   private
   type(location_type) :: rad_loc
   real(r8)            :: direction(3)
   real(r8)            :: nyquistvel
end type radial_vel_type

! Derived type for reflectivity obs.  Currently no auxiliary information
! is neded.

! Max number of these obs supported in a single run.  Can be overridden
! by specifying a larger value in the namelist.
! integer :: max_radar_ref_obs = 1000000
!define here:
!  :: globally_defined_reflectivity_metadata_variable(max_radar_ref_obs)

integer :: velkeycount = 0 ! cumulative index into rad. velocity metadata
integer :: refkeycount = 0 ! cumulative index into rad. reflectivity metadata

! for error message content
character(len=129) :: msgstring

! FIXME: is this accurate?
! constants which depend on the microphysics scheme used.  should be set
! in the namelist to match the case being run.

! FIXME: add description, and print this out in a msg at run time
! so the values are in the log for documentation of this run.
real(r8), parameter :: dief = 0.224_r8

! FIXME: add description
real(r8), parameter :: n0r = 8.0e6_r8, n0g = 4.0e4_r8, n0s = 3.0e6_r8
real(r8), parameter :: rho_r = 1000.0_r8, rho_g = 917.0_r8, rho_s = 100.0_r8

!-------------------------------------------------------------
! Namelist with default values
! 
! Obsolete: convert_to_dbz and dbz_threshold
!  convert_to_dbz and dbz_threshold have both been removed from the namelist.
!  Values will always be converted to dBZ, and threshold was only used to 
!  ensure the log() call never saw a real 0.0_r8.  Please remove these 
!  values from your namelist to avoid a run-time error. 
!
! FIXME: removing these values is a user-visible, non-backwards compatibility
!  issue.   When this file is committed remove this comment but ensure that
!  the dart main web page is updated to point out this change.
!
! FIXME: this table was with the dbz_threshold comments (which were removed), 
!        but it seems useful to keep around someplace anyway:
!
!    Some useful values: Z = 3.163  ->   5.0 dBZ
!                        Z = 2.512  ->   4.0 dBZ
!                        Z = 1.259  ->   1.0 dBZ
!                        Z = 1.0233 ->   0.1 dBZ
!                        Z = 0.1    -> -10.0 dBZ
!                        Z = 0.01   -> -20.0 dBZ
!                        Z = 0.001  -> -30.0 dBZ
!
!
! There are two replicated sets of 3 namelist values below.
! In each case, there is 1 logical value and 2 numeric values.
! If the logical is false, the numeric values are ignored.
! If true, then the 2 numeric values are:
!  1) a threshold which determines if this value is to be clamped, and 
!  2) what it should be clamped to.
! These are separate to allow, for example, the option of setting all values
! below -20 dBZ to -40 dBZ.
!
! The next 3 namelist items apply to the incoming observation values.  
! They are in the namelist so they can be changed at runtime, instead 
! of set only when the observation file is originally generated. 
!
! apply_ref_limit_to_obs:    
!   Logical.  If .TRUE., replace reflectivities below "reflectivity_limit_obs" 
!   with "lowest_reflectivity_obs". If .FALSE., ignore "reflectivity_limit"
!   and "lowest_reflectivity_obs".  
!
! reflectivity_limit_obs:    
!   Observed reflectivity below this value is set to "lowest_reflectivity_obs". 
!   Units are dBZ.
!
! lowest_reflectivity_obs:   
!   If reflectivity < reflectivity_limit_obs, reflectivity is set to this value.
!   The default value of 'missing' is useful for perfect model experiments.
!   Suggested options: lowest_reflectivity_obs = missing_r8 or 
!                      lowest_reflectivity_obs = any real value.
!   Units are dBZ.
!
! FIXME:  it is unclear if missing_r8 is actually a good suggestion here
!         since if the QC is good, this value (-888888.0) might actually be
!         assimilated as a valid observation.
!
! The next 3 namelist items apply to the forward operator values (the returned
! value from each ensemble member predicting what the observation value should
! be given the current state in this particular member).
!
! apply_ref_limit_to_fwdop:  
!   Similar to "apply_ref_limit_to_obs" but applied to the forward operator.
!
! reflectivity_limit_fwdop:  
!   Similar to "reflectivity_limit_obs" for forward operator.
!
! lowest_reflectivity_fwdop: 
!   Similar to "lowest_reflectivity_obs" for forward operator.
!
! FIXME:  The discussion is still open on whether it is a good idea to allow 
! different threshold values for the observation vs forward operator.  There
! was general agreement that 1) it seemed useful to allow the threshold to be
! different from the lowest value, and 2) also good to have a separate logical
! flag for the observation vs the forward operators.  
! But it was less clear if it was a good thing to allow users to specify
! *different* numeric values between the observations and forward operators,
! since what is observed and what is computed are supposed to be from matching
! distributions.  Verify that this is wanted, and if so, remove this comment.
! Otherwise, combine the 2 thresholds and lowest values into 1 set.
!
! max_radial_vel_obs:
!  Integer value.  Maximum number of observations of this type to support at
!  run time.  This is combined total of all obs_seq files, for example if 
!  running the obs_diag program which potentially opens multiple obs_seq.final
!  files, or the obs merge program which can also open multiple obs files.
!  Default is 1,000,000 observations.
! FIXME: this needs to be added to the .html file

logical  :: apply_ref_limit_to_obs    = .false.
real(r8) :: reflectivity_limit_obs    = 0.0_r8
real(r8) :: lowest_reflectivity_obs   = missing_r8    ! FIXME?
logical  :: apply_ref_limit_to_fwdop  = .false.
real(r8) :: reflectivity_limit_fwdop  = 0.0_r8
real(r8) :: lowest_reflectivity_fwdop = missing_r8    ! FIXME?
integer  :: max_radial_vel_obs = 1000000  ! max number of this type of obs.


namelist /obs_def_radar_mod_nml/ apply_ref_limit_to_obs,    &
                                 reflectivity_limit_obs,    &
                                 lowest_reflectivity_obs,   &
                                 apply_ref_limit_to_fwdop,  &
                                 reflectivity_limit_fwdop,  &
                                 lowest_reflectivity_fwdop, &
                                 max_radial_vel_obs


! module global storage for auxiliary obs data
type(radial_vel_type), allocatable :: radial_vel_data(:)

contains

!----------------------------------------------------------------------

subroutine initialize_module

implicit none

integer :: iunit, io, rc

call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_def_radar_mod_nml", iunit)
read(iunit, nml = obs_def_radar_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_radar_mod_nml")

! Record the namelist values used for the run ... 
if (do_output()) write(nmlfileunit, nml=obs_def_radar_mod_nml)
if (do_output()) write(     *     , nml=obs_def_radar_mod_nml)

! Allocate space for the auxiliary information associated with each obs
! This code must be placed after reading the namelist, so the user can
! increase or decrease the number of obs supported and use more or less
! memory at run time.
allocate(radial_vel_data(max_radial_vel_obs), stat = rc)
if (rc /= 0) then            
   write(msgstring, *) 'initial allocation failed for radial vel obs data,', &
                       'itemcount = ', max_radial_vel_obs
   call error_handler(E_ERR,'initialize_module', msgstring, &
                      source, revision, revdate)
endif                        

! log the values used for the constants
call print_constants()

module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine write_radial_vel(velkey, ifile, fform)

integer,          intent(in)           :: velkey, ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      call write_location(ifile,  radial_vel_data(velkey)%rad_loc, fileformat) 
      call write_orientation(ifile, radial_vel_data(velkey)%direction(:), fileformat) 
      call write_nyquistvel(ifile, radial_vel_data(velkey)%nyquistvel, fileformat) 
      ! Write out the obs_def velkey for this observation 
      write(ifile) velkey
   CASE DEFAULT
      ! Write the 5 character identifier for verbose formatted output
      write(ifile, 11)
11    format('platform')
      call write_location(ifile,  radial_vel_data(velkey)%rad_loc, fileformat) 
      call write_orientation(ifile, radial_vel_data(velkey)%direction(:), fileformat) 
      call write_nyquistvel(ifile, radial_vel_data(velkey)%nyquistvel, fileformat) 
      ! Write out the obs_def velkey for this observation 
      write(ifile, *) velkey
END SELECT

end subroutine write_radial_vel

!----------------------------------------------------------------------

subroutine get_obs_def_radial_vel(velkey, radarloc, beamdir, velnyquist)

integer, intent(in)              :: velkey
type(location_type), intent(out) :: radarloc
real(r8), intent(out)            :: beamdir(3)
real(r8), intent(out)            :: velnyquist

radarloc   = radial_vel_data(velkey)%rad_loc
beamdir    = radial_vel_data(velkey)%direction(:)
velnyquist = radial_vel_data(velkey)%nyquistvel

end subroutine get_obs_def_radial_vel

!----------------------------------------------------------------------

subroutine read_radial_vel(velkey, ifile, fform)

! Main read subroutine for the obs kind radial velocity metadata.

! location: Refers to the lat/lon/height/vertical coordinate option
!           for the radar. Its type is defined through location_type in
!           location_mod. Uses the read_location function in
!           location_mod.
!
! orientation: This is a 3-element array specific to the radar module:
!              orientation(1) = sin(azimuth)*cos(elevation)
!              orientation(2) = cos(azimuth)*cos(elevation)
!              orientation(3) = sin(elevation)

integer,          intent(out)          :: velkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=8)    :: header
character(len=32)   :: fileformat
real(r8)            :: orientation(3)
type(location_type) :: location
real(r8)            :: nyquistv

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      location = read_location(ifile, fileformat)
      orientation = read_orientation(ifile, fileformat)
      nyquistv = read_nyquistvel(ifile, fileformat)
      ! Read in the velkey for this particular observation
      read(ifile) velkey
   CASE DEFAULT
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT='(a8)') header
      if(header /= 'platform') then
         call error_handler(E_ERR,'obs_def_radar_mod:read_radial_vel', &
              'Expected location header "platform" in input file', &
              source, revision, revdate)
      endif
      location = read_location(ifile, fileformat)
      orientation = read_orientation(ifile, fileformat)
      nyquistv = read_nyquistvel(ifile, fileformat)
      ! Read in the velkey for this particular observation
      read(ifile, *) velkey
END SELECT

call set_radial_vel(velkey, location, orientation, nyquistv)

end subroutine read_radial_vel

!----------------------------------------------------------------------

subroutine set_radial_vel(velkey, rad_location, rad_orientation, rad_nyquistv)

integer,             intent(out) :: velkey
real(r8),            intent(in)  :: rad_orientation(3)
type(location_type), intent(in)  :: rad_location
real(r8),            intent(in)  :: rad_nyquistv

if ( .not. module_initialized ) call initialize_module

! Do the increment here: this code can now be called from anywhere
! and it will have a consistently incremented key count.

! the total velocity metadata key count from all sequences
velkeycount = velkeycount + 1    
velkey = velkeycount             ! copied to the output variable

!Make sure enough space is allocated
if(velkey > max_radial_vel_obs) then
   write(*, *) 'velkey (',velkey,') exceeds max_radial_vel_obs (',max_radial_vel_obs,')'
   call error_handler(E_ERR,'read_radial_vel:set_radial_vel', &
              'Increase max_radial_vel_obs.', source, revision, revdate)
endif

radial_vel_data(velkey)%rad_loc      = rad_location
radial_vel_data(velkey)%direction(:) = rad_orientation
radial_vel_data(velkey)%nyquistvel   = rad_nyquistv

end subroutine set_radial_vel

!----------------------------------------------------------------------

subroutine write_radar_ref(refkey, ifile, fform)

integer,          intent(in)           :: refkey, ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! At this point, this is empty as there is no metadata for radar reflectivity
! in the current obs sequence format. In the future, if more metadata pieces
! are added, this part should be modified.

end subroutine write_radar_ref

!----------------------------------------------------------------------

subroutine get_obs_def_radar_ref(refkey, is_this_dbz)

! Deprecated routine; values are always in dbz now (and not Z)
! Still here; will just always return true.

integer, intent(in)  :: refkey
logical, intent(out) :: is_this_dbz
 
is_this_dbz = .true.
 
end subroutine get_obs_def_radar_ref

!----------------------------------------------------------------------

subroutine read_radar_ref(obsvalue, refkey, ifile, fform)

! Main read subroutine for the obs kind radar reflectivity metadata.

! reftype: Denotes whether reflectivity obs units are in Z (reftype = .false.)
!                                              or are in dBZ (reftype = .true.)

real(r8),         intent(inout)        :: obsvalue
integer,          intent(out)          :: refkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32)   :: fileformat
!logical             :: reftype

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! At this point, no metadata is available for radar reflectivity. For this
! reason, internal counter "refkey" has also no consequence. In the
! future, if more metadata pieces are added, this part should be modified.
refkeycount = refkeycount + 1    ! the total reflectivity metadata key count from all sequences
refkey = refkeycount             ! copied to the output variable

! All computations are now done in dbz.  There used to be an option to work
! directly in Z, but that has been removed.
! FIXME: verify what happens if the QC is good but the 'lowest' value is
! in fact missing_r8.  This seems like it could cause problems.
if (obsvalue /= missing_r8) then
   if ( (apply_ref_limit_to_obs) .and. (obsvalue < reflectivity_limit_obs) ) then
      obsvalue = lowest_reflectivity_obs
   endif
endif

! Again, since there is no metadata for reflectivity to track, we don't
! need the following call at the moment; thus the
! following commented line and the subroutine "set_radar_ref". In the future,
! if more metadata pieces are added, this part should be modified and
! "set_radar_ref" should be uncommented.

!call set_radar_ref(refkey,some_reflectivity_metadata)

end subroutine read_radar_ref

!----------------------------------------------------------------------
!
!subroutine set_radar_ref(refkey, some_reflectivity_metadata)
!
!integer, intent(in)     :: refkey
!define here, intent(in) :: some_reflectivity_metadata
!
!if ( .not. module_initialized ) call initialize_module
!
!!Make sure enough space is allocated
!if(refkey > max_radar_ref_obs) then
!   write(*, *) 'refkey (',refkey,') exceeds max_radar_ref_obs (',max_radar_ref_obs,')'
!   call error_handler(E_ERR,'obs_def_radar_mod:read_radar_ref:set_radar_ref', &
!              'Increase max_radar_ref_obs.', source, revision, revdate)
!endif
!
!globally_defined_reflectivity_metadata_variable(refkey) = some_reflectivity_metadata
!
!end subroutine set_radar_ref

!----------------------------------------------------------------------

subroutine interactive_radial_vel(velkey)

! Interactively reads in location and orientation components of the obs
! kind radial velocity. See read_radial_vel for more information on these
! two components.

! Uses interactive_location of location_mod and the local subroutines
! interactive_orientation and set_radial_vel.

! velkey is internally incremented and only counts the index for this
! specialized observation kind. It is written in the obs file after 
! the radar location and beam direction.

implicit none

integer, intent(out) :: velkey

real(r8)             :: orientation(3)
type(location_type)  :: location
real(r8)             :: nyquistv

if ( .not. module_initialized ) call initialize_module

!Interactively obtain radar location and beam direction information
!Note: Obs location is additionally inquired by the standard DART module
!      "interactive_obs_def". No check is performed here whether radar
!      location, obs location, and beam direction match. It is the user's
!      responsibility to make sure that they do.


write(*, *)
write(*, *) 'Beginning to inquire information on radar location.'
write(*, *)
write(*, *) 'WARNING!! Make sure that you select 3 (height) for the vertical'
write(*, *) 'coordinate option and enter height in units gpm.  This location'
write(*, *) 'question will be repeated again later for this same observation,'
write(*, *) 'and you must enter the same longitude, latitude, and height'
write(*, *) 'as you enter here.'
write(*, *)
! FIXME:  this is the location of the center of the radar?  then this comment
! is wrong - these locations cannot be the same.  this one is the center
! of the radar cone; the second location is the reflection position (the actual
! observation location).

call interactive_location(location)

write(*, *)
write(*, *) 'Beginning to inquire information on radar beam direction.'
write(*, *)

call interactive_orientation(orientation)

write(*, *)
write(*, *) 'Beginning to inquire information on radar Nyquist velocity.'
write(*, *)

call interactive_nyquistvel(nyquistv)

call set_radial_vel(velkey, location, orientation, nyquistv)

write(*, *)
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_radial_vel

!----------------------------------------------------------------------

subroutine interactive_radar_ref(refkey)

! Interactively (sort of) reads in reftype component of the obs
! kind radar reflectivity. See read_radar_ref for more information on these
! two components.

! Uses the local subroutine interactive_reftype.

! refkey is internally incremented and only counts the index for this specialized
! observation kind. It is written in the obs file after the reflectivity type.

implicit none

integer, intent(out) :: refkey

if ( .not. module_initialized ) call initialize_module

! generate a unique local key number
refkeycount = refkeycount + 1
refkey = refkeycount

!Make sure enough space is allocated
if(refkey >= max_radar_ref_obs) then
   write(*, *) 'refkey (',refkey,') exceeds max_radar_ref_obs (',max_radar_ref_obs,')'
   call error_handler(E_ERR,'obs_def_radar_mod:interactive_radar_ref', &
              'Increase max_radar_ref_obs.', source, revision, revdate)
endif

!Now interactively obtain reflectivity type information
!Note: Reflectivity type will be actually read from the namelist variable
!      convert_to_dbz. Thus, for reftype, the code is not really interactive.
!      Nevertheless, this was chosen to be this way to be consistent with the
!      rest of the code and make future modifications more straightforward.

write(*, *)
write(*, *) 'Beginning to inquire information on reflectivity type.'
write(*, *) 'This will be read from the namelist variable convert_to_dbz.'
write(*, *)

! If reflectivity metadata is added in the future, insert here a call to a
! subroutine to read that metadata interactively (and write the subroutine
! that goes with it!)


write(*, *)
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_radar_ref

!----------------------------------------------------------------------

subroutine interactive_orientation(orientation)

implicit none

real(r8), intent(out) :: orientation(3)

write(*, *) 'Input first component: sin(azimuth)*cos(elevation)'
read(*, *)  orientation(1)

do while((orientation(1) > 1.0_r8) .or. (orientation(1) < -1.0_r8))
   write(*, *) 'Input value is not allowed, please enter a value -1.0 to 1.0'
   read(*, *) orientation(1)
end do

write(*, *) 'Input second component: cos(azimuth)*cos(elevation)'
read(*, *)  orientation(2)

do while((orientation(2) > 1.0_r8) .or. (orientation(2) < -1.0_r8))
   write(*, *) 'Input value is not allowed, please enter a value -1.0 to 1.0'
   read(*, *) orientation(2)
end do

write(*, *) 'Input third component: sin(elevation)'
read(*, *)  orientation(3)

do while((orientation(3) > 1.0_r8) .or. (orientation(3) < 0.0_r8))
   write(*, *) 'Input value is not allowed, please enter a value 0.0 to 1.0'
   read(*, *) orientation(3)
end do

end subroutine interactive_orientation

!----------------------------------------------------------------------

subroutine interactive_nyquistvel(nyquistv)

implicit none

real(r8), intent(out) :: nyquistv

write(*, *) 'Input Nyquist velocity for this obs point'
read(*, *)  nyquistv

end subroutine interactive_nyquistvel

!----------------------------------------------------------------------

subroutine get_expected_radial_vel(state_vector, location, velkey, vr, istatus)

! This is the main forward operator routine for radar Doppler velocity.
! Given an ob DART location, computes model-predicted radial velocity at
! the same location.

! In addition to u,v,w, from which and the direction of the radar beam is
! calculated the along-beam component of the 3-d wind vector, we also need
! the reflectivity value to compute the terminal velocity of the hydrometeor,
! so that w can be adjusted for it. See also get_expected_radar_ref.

! Reference: Lin et al., 1983 (J. Climate Appl.Meteor., 1065-1092)
! Note that the reflectivity-weighted mean terminal velocities are used here.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: velkey
real(r8),            intent(out) :: vr
integer,             intent(out) :: istatus

! FIXME: gravity is currently hardcoded here.  We should find a way to let
! the model have input if it uses a slightly different value for G, or if it
! is working on a different planet with an entirely different set of constants.
! Question:  how much impact on the results does changing G have?
real(r8), parameter :: gravity  = 9.81_r8     ! FIXME: dart val, not wrf
real(r8), parameter :: a        = 8.42e20_r8
real(r8), parameter :: b        = 0.8_r8
real(r8), parameter :: c        = 4.84e18_r8
real(r8), parameter :: d        = 0.25_r8
real(r8), parameter :: CD       = 0.6_r8
real(r8), parameter :: rhos0    = 1.0_r8
real(r8), parameter :: e        = 4.0_r8*gravity*rho_g/(3.0_r8*CD)
real(r8), parameter :: f        = 0.5_r8
real(r8), parameter :: gam7b    = 3376.92_r8
real(r8), parameter :: gam7d    = 1155.38_r8
real(r8), parameter :: gam7f    = 1871.25_r8
real(r8), parameter :: powr     = (7.0_r8 + b)/4.0_r8
real(r8), parameter :: pows     = (7.0_r8 + d)/4.0_r8
real(r8), parameter :: powg_dry = (7.0_r8 + f)/4.0_r8
real(r8), parameter :: powg_wet = 1.7875_r8

real(r8) :: ar
real(r8) :: as_wet
real(r8) :: as_dry
real(r8) :: ag_dry

real(r8) :: u, v, w, qr, qg, qs, alpha, wt, rho, temp, ref
real(r8) :: precip_r, precip_s, precip_g

ar       = n0r*a*gam7b / (PI*rho_r*n0r)**powr
as_wet   = n0s*c*gam7d / (PI*rho_s*n0s)**pows
as_dry   = dief*((rho_s/rho_r)**2.0_r8)*as_wet
ag_dry   = 1.0e18_r8*dief*((rho_g/rho_r)**2.0_r8)* &
           n0g*gam7f/(PI*rho_g*n0g)**powg_dry

if ( .not. module_initialized ) call initialize_module

call interpolate(state_vector, location, KIND_U_WIND_COMPONENT, u, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_V_WIND_COMPONENT, v, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_VERTICAL_VELOCITY, w, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO, qr, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   vr = missing_r8
   return
endif

precip_r = rho * qr
precip_s = rho * qs
precip_g = rho * qg
alpha    = sqrt(rhos0/rho)

wt = 0.0_r8

! Computing reflectivity-weighted terminal velocity wt
! RAIN
if (qr >= 1.0e-6_r8) then
   wt = alpha * ar * (precip_r**powr)
endif
! HAIL/GRAUPEL
if (qg >= 1.0e-6_r8) then
   wt = wt + sqrt(e/rho)*ag_dry*(precip_g**powg_dry)
   ! Note that we are assuming dry graupel/hail surface
   ! Thus, no temperature check for hail/graupel (see snow for comparison)
   ! ag_wet = ((7.2e20_r8)**0.95_r8)*gam7f/(720.0_r8*(n0g**0.8375_r8))
   ! wt     = wt + sqrt(e/rho)*ag_wet*(((rho * qg)/(PI*rho_g))**powg_wet)
endif
! SNOW
if (qs >= 1.0e-6_r8) then
   if ( temp < 273.15_r8 ) then
      wt = wt + alpha*as_dry*(precip_s**pows)
   else
      wt = wt + alpha*as_wet*(precip_s**pows)
   endif
endif

if (wt > 0.0_r8) then
   call get_reflectivity(qr, qg, qs, rho, temp, ref)
   wt = wt/ref
endif

vr = radial_vel_data(velkey)%direction(1) * u +    &
     radial_vel_data(velkey)%direction(2) * v +    &
     radial_vel_data(velkey)%direction(3) * (w-wt)

end subroutine get_expected_radial_vel

!----------------------------------------------------------------------

subroutine get_expected_radar_ref(state_vector, location, ref, istatus)
!
! This is the main forward operator routine for radar reflectivity.
! Given an ob DART location, computes model-predicted radar reflectivity at
! the same location.

! Radar reflectivity = 10 * log_10( radar reflectivity factor) in dBZ.

! If apply_ref_limit_to_fwdop, reflectivities below
! reflectivity_limit_fwdop are set to lowest_reflectivity_fwdop.

! Subroutine "interpolate" ultimately calls model_mod routine "model_interpolate"
! to get model values of qr, qg, qs, rho, and temp at the ob location. Then
! the routine "get_reflectivity" is called to compute the radar reflectivity
! factor value, Z, that corresponds to the hydrometeor and thermodynamic values.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: ref
integer,             intent(out) :: istatus

real(r8) :: qr, qg, qs, rho, temp

if ( .not. module_initialized ) call initialize_module

!FIXME: discussion of whether this code should first try to ask
! the model to interpolate reflectivity directly, and if it cannot,
! then resort to asking for the following 5 fields. 

call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO , qr , istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   ref = missing_r8
   return
endif

call get_reflectivity(qr, qg, qs, rho, temp, ref)

! Always convert to dbz.  Make sure value is always slightly positive.
ref = 10.0_r8 * log10(max(tiny(ref), ref))

if ( (apply_ref_limit_to_fwdop) .and. (ref < reflectivity_limit_fwdop) .and. &
     (ref /= missing_r8) ) then
   ref = lowest_reflectivity_fwdop
endif

if (ref == missing_r8) then
   istatus = 1
endif

end subroutine get_expected_radar_ref

!----------------------------------------------------------------------

subroutine get_reflectivity(qr, qg, qs, rho, temp, ref)
!
! Computes "radar reflectivity factor" (Z) in mm^6 m^-3
!
! References: Ferrier, 1994 (JAS, 249-280)
!             Smith et al., 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)

! According to Smith (1984), there are two choices for the dielectric
! factor (dief), depending on how the snow particle sizes are specified.
! If melted raindrop diameters are used, then the factor is 0.224.  If
! equivalent ice sphere diameters are used, then the factor is 0.189.

real(r8), intent(in)  :: qr, qg, qs, rho, temp
real(r8), intent(out) :: ref

real(r8) :: precip
real(r8) :: ar 
real(r8) :: ag_dry
real(r8) :: as_wet
real(r8) :: as_dry


ar     = 7.2e20_r8  /  (((PI*rho_r)**1.75_r8)*(n0r**0.75_r8))
ag_dry = dief*((rho_g/rho_r)**2.0_r8)*7.2e20_r8/ &
           (((PI*rho_g)**1.75_r8)*(n0g**0.75_r8))
as_wet = 7.2e20_r8/(((PI*rho_s)**1.75_r8)*(n0s**0.75_r8))
as_dry = dief*((rho_s/rho_r)**2.0_r8)*as_wet 

if ( .not. module_initialized ) call initialize_module

ref = 0.0_r8

! RAIN
if ( qr >= 1.0e-6_r8 ) then
   precip = rho * qr
   ref = ref + ar * (precip**1.75_r8)
endif

! HAIL / GRAUPEL
if ( qg >= 1.0e-6_r8 ) then
    precip = rho * qg
    ref = ref + ag_dry * (precip**1.75_r8)
    ! Note that we assume dry surface for hail/graupel
    ! Thus, no temperature check for hail/graupel (see snow for comparison)
    ! ag_wet = (7.2e20_r8/(((PI*rho_g)**1.75_r8)*(n0g**0.75_r8)))**0.95_r8
    ! ref    = ref + ag_wet * (precip**1.6625_r8)
endif

! SNOW
if ( qs >= 1.0e-6_r8 ) then
   precip = rho * qs
   if ( temp < 273.15_r8 ) then
      ref = ref + as_dry * (precip**1.75_r8)
   else
      ref = ref + as_wet * (precip**1.75_r8)
   endif
endif

end subroutine get_reflectivity

!----------------------------------------------------------------------

subroutine write_orientation(ifile, orientation, fform)

! Writes orientation to obs file.
! This variable is only associated with the obs kind radial velocity.
! See read_radial_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: orientation(3)
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) orientation(1), orientation(2), orientation(3)
   CASE DEFAULT
      write(ifile, '(''dir3d'')' ) 
      write(ifile, *) orientation(1), orientation(2), orientation(3)
END SELECT

end subroutine write_orientation

!----------------------------------------------------------------------

subroutine write_nyquistvel(ifile, nyquistv, fform)

! Writes Nyquist velocity to obs file.
! This variable is only associated with the obs kind radial velocity.
! See read_radial_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8),                   intent(in) :: nyquistv
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) nyquistv
   CASE DEFAULT
      write(ifile, *) nyquistv
END SELECT

end subroutine write_nyquistvel

!----------------------------------------------------------------------

subroutine write_reftype(ifile, fform)

! Writes reflectivity type to obs file.
! This variable is only associated with the obs kind radar reflectivity.
! See read_radar_ref for additional discussion.

! Note that when writing to the obs file, we already know that this value
! should be equal to the namelist variable convert_to_dbz.

integer,                    intent(in) :: ifile
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) convert_to_dbz
   CASE DEFAULT
      write(ifile, *) convert_to_dbz
END SELECT

end subroutine write_reftype

!----------------------------------------------------------------------

function read_orientation(ifile, fform)

! Reads orientation from file that was written by write_orientation.
! This variable is only associated with the obs kind radial velocity.
! See read_radial_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8)                               :: read_orientation(3)
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_orientation(1), read_orientation(2), read_orientation(3)
   CASE DEFAULT
      read(ifile, '(a5)' ) header

      if(header /= 'dir3d') then
         write(errstring,*)'Expected orientation header "dir3d" in input file, got ', header
         call error_handler(E_ERR, 'read_orientation', errstring, source, revision, revdate)
      endif
! Now read the orientation data value
      read(ifile, *) read_orientation(1), read_orientation(2), read_orientation(3)
END SELECT

end function read_orientation

!----------------------------------------------------------------------

function read_nyquistvel(ifile, fform)

! Reads Nyquist velocity from file that was written by write_nyquistvel.
! This variable is only associated with the obs kind radial velocity.
! See read_radial_vel for additional discussion.

integer,                    intent(in) :: ifile
real(r8)                               :: read_nyquistvel
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_nyquistvel
   CASE DEFAULT
      read(ifile, *) read_nyquistvel
END SELECT

end function read_nyquistvel

!----------------------------------------------------------------------

function read_reftype(ifile, fform)

! Reads reflectivity type from file that was written by write_reftype.
! This variable is only associated with the obs kind radar reflectivity.
! See read_radar_ref for additional discussion.

integer,                    intent(in) :: ifile
logical                                :: read_reftype
character(len=*), intent(in), optional :: fform

character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_reftype
   CASE DEFAULT
      read(ifile, *) read_reftype
END SELECT

end function read_reftype

!----------------------------------------------------------------------
subroutine print_constants()

! Log the constants set in the code.
! Prints to both the log file and standard output.

write(msgstring, *) 'constants used in this module:'
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  dief = ', dief
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  n0r = ', n0r
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  n0g = ', n0g
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  n0s = ', n0s
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  rho_r = ', rho_r
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  rho_g = ', rho_g
call error_handler(E_MSG,'', msgstring, '', '', '')

write(msgstring, *) '  rho_s = ', rho_s
call error_handler(E_MSG,'', msgstring, '', '', '')

end subroutine print_constants

end module obs_def_radar_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
