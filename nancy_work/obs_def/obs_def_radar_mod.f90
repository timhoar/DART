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
! PRECIPITATION_FALL_RATE, KIND_WEIGHTED_VELOCITY
! END DART PREPROCESS KIND LIST
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_radar_mod, only : write_radial_vel, read_radial_vel, &
!                            interactive_radial_vel, get_expected_radial_vel, &
!                            read_radar_ref, get_expected_radar_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(DOPPLER_RADIAL_VELOCITY)
!     call get_expected_radial_vel(state, location, obs_def%key, obs_val, istatus)
!  case(RADAR_REFLECTIVITY)
!     call get_expected_radar_ref(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(DOPPLER_RADIAL_VELOCITY)
!      call read_radial_vel(obs_def%key, ifile, fileformat)
!   case(RADAR_REFLECTIVITY)
!      call read_radar_ref(obs_val, obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(DOPPLER_RADIAL_VELOCITY)
!      call write_radial_vel(obs_def%key, ifile, fileformat)
!   case(RADAR_REFLECTIVITY)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(DOPPLER_RADIAL_VELOCITY)
!      call interactive_radial_vel(obs_def%key)
!   case(RADAR_REFLECTIVITY)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_radar_mod

use        types_mod, only : r8, missing_r8, PI, deg2rad
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             check_namelist_read, find_namelist_in_file, &
                             nmlfileunit, do_output
use     location_mod, only : location_type, write_location, read_location, &
                             interactive_location, get_location
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_TEMPERATURE, KIND_VERTICAL_VELOCITY, &
                             KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, &
                             KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO, &
                             KIND_WEIGHTED_VELOCITY, KIND_RADAR_REFLECTIVITY

implicit none
private

public :: read_radar_ref, get_expected_radar_ref,                         &
          write_radial_vel, read_radial_vel, interactive_radial_vel,      &
          get_expected_radial_vel, get_obs_def_radial_vel, set_radial_vel

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical :: module_initialized = .false.

! Derived type for radial velocity.  Contains auxiliary information stored
! with each obs of this type; used to compute the forward operator.
! See more extensive comments in the interactive_radial_vel() routine for
! expected units, etc.  Technically, the radar location is unused in the
! forward operators currently in this code, but it may be useful for post
! processing or diagnostics, especially if multiple radar locations are
! in the same file.  
type radial_vel_type
   private
   type(location_type) :: radar_location      ! location of radar 
   real(r8)            :: beam_direction(3)   ! direction of beam
   real(r8)            :: nyquist_velocity    ! nyquist velocity
end type radial_vel_type

! Cumulative index into radial velocity metadata array
integer :: velkeycount = 0 

! For error message content
character(len=129) :: msgstring

! Values which are initialized at run time so some can be changed by
! namelist.  After initialization, treated as parameters (values not changed).
! FIXME: short line for each parm to say what it is
real(r8) :: param_gravity    ! gravity
real(r8) :: param_a          !
real(r8) :: param_b          !
real(r8) :: param_c          !
real(r8) :: param_d          !
real(r8) :: param_CD         !
real(r8) :: param_rhos0      !
real(r8) :: param_e          !
real(r8) :: param_f          !

real(r8) :: param_gam7b      !
real(r8) :: param_gam7d      !
real(r8) :: param_gam7f      !

real(r8) :: param_powr       !
real(r8) :: param_pows       !
real(r8) :: param_powg_dry   !
real(r8) :: param_powg_wet   !

real(r8) :: param_ar_v       !
real(r8) :: param_as_wet_v   !
real(r8) :: param_as_dry_v   !
real(r8) :: param_ag_wet_v   !
real(r8) :: param_ag_dry_v   !

real(r8) :: param_ar_f       !
real(r8) :: param_as_wet_f   !
real(r8) :: param_as_dry_f   !
real(r8) :: param_ag_wet_f   !
real(r8) :: param_ag_dry_f   !

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
! There are two replicated sets of 3 namelist values below.
! In each case, there is 1 logical value and 2 numeric values.
! If the logical is false, the numeric values are ignored.
! If true, then the 2 numeric values are:
!  1) a threshold which determines if this value is to be clamped, and 
!  2) what it should be clamped to.
! These are separate to allow, for example, the option of setting all values
! below -20 dBZ to -40 dBZ.
!
! If the observation value or the forward operator already has a value of
! missing_r8 it is assumed either the istatus is marked as failed (for the
! forward operator) or that the QC (quality control) flag is set to not
! assimilate this observation and that value is left unchanged regardless of
! the setting on the apply_ref_limit flag.  Note however that it is not a
! good idea to reset a good but small observation value to missing_r8 -- do
! not use it as the lowest_reflectivity setting.
!
! The next 3 namelist items apply to the incoming observation values.  
! They are in the namelist so they can be changed at runtime, instead 
! of set only when the observation file is originally generated. 
!
! apply_ref_limit_to_obs:    
!   Logical.  If .TRUE., replace reflectivities below "reflectivity_limit_obs" 
!   with "lowest_reflectivity_obs". If .FALSE., ignore "reflectivity_limit"
!   and "lowest_reflectivity_obs".  Defaults to .FALSE.
!
! reflectivity_limit_obs:    
!   Observed reflectivity below this value is set to "lowest_reflectivity_obs". 
!   Units are dBZ.  Defaults to -10.0.
!
! lowest_reflectivity_obs:   
!   If reflectivity < reflectivity_limit_obs, reflectivity is set to this value.
!   Units are dBZ.  Defaults to -10.0.
!
! The next 3 namelist items apply to the forward operator values (the returned
! value from each ensemble member predicting what the observation value should
! be given the current state in this particular member).
!
! apply_ref_limit_to_fwd_op:  
!   Same as "apply_ref_limit_to_obs" but applied to the forward operator.
!
! reflectivity_limit_fwd_op:  
!   Same as "reflectivity_limit_obs" for forward operator.
!
! lowest_reflectivity_fwd_op: 
!   Same as "lowest_reflectivity_obs" for forward operator.
!
!
! max_radial_vel_obs:
!  Integer value.  Maximum number of observations of this type to support at
!  run time.  This is combined total of all obs_seq files, for example the
!  observation diagnostic program potentially opens multiple obs_seq.final
!  files, or the obs merge program can also open multiple obs files.
!  Default is 1,000,000 observations.
!
! dielectric_factor:
!  According to Smith (1984), there are two choices for the dielectric
!  factor, depending on how the snow particle sizes are specified.
!  If melted raindrop diameters are used, then the factor is 0.224.  If
!  equivalent ice sphere diameters are used, then the factor is 0.189.
!  So there are two values in the paper and possibly others.   
!  FIXME: what is the preferred default value here?
!
! n0_rain, n0_graupel, n0_snow:
!  Intercept parameter for size distribution of each constituent.
!  FIXME: is this right?  and what are defaults?
!
! rho_rain, rho_graupel, rho_snow:
!  Density of each hydrometeor type.  FIXME: defaults?
!
! use_wet_graupel:
!  Logical.  If .true. include a contribution in the precipitation fall speed
!  for wet graupel.  Defaults to .false.   FIXME: is this right?  ok?
!
! FIXME: everything from max_radial_vel_obs down needs to be added to 
! the .html doc file

logical  :: apply_ref_limit_to_obs     = .false.
real(r8) :: reflectivity_limit_obs     = -10.0_r8
real(r8) :: lowest_reflectivity_obs    = -10.0_r8
logical  :: apply_ref_limit_to_fwd_op  = .false.
real(r8) :: reflectivity_limit_fwd_op  = -10.0_r8
real(r8) :: lowest_reflectivity_fwd_op = -10.0_r8
integer  :: max_radial_vel_obs         = 1000000  ! max of this type of obs
logical  :: use_wet_graupel            = .FALSE.

! Constants which depend on the microphysics scheme used.  Should be set
! in the namelist to match the case being run.

! FIXME: add descriptions here
real(r8) :: dielectric_factor =  0.224_r8
real(r8) :: n0_rain           =  8.0e6_r8
real(r8) :: n0_graupel        =  4.0e4_r8
real(r8) :: n0_snow           =  3.0e6_r8
real(r8) :: rho_rain          = 1000.0_r8
real(r8) :: rho_graupel       =  917.0_r8
real(r8) :: rho_snow          =  100.0_r8



! FIXME: update .nml file to match
namelist /obs_def_radar_mod_nml/ apply_ref_limit_to_obs,     &
                                 reflectivity_limit_obs,     &
                                 lowest_reflectivity_obs,    &
                                 apply_ref_limit_to_fwd_op,  &
                                 reflectivity_limit_fwd_op,  &
                                 lowest_reflectivity_fwd_op, &
                                 max_radial_vel_obs,         &
                                 use_wet_graupel,            &
                                 dielectric_factor,          &
                                 n0_rain,                    &
                                 n0_graupel,                 &
                                 n0_snow,                    &
                                 rho_rain,                   &
                                 rho_graupel,                &
                                 rho_snow


! Module global storage for auxiliary obs data, allocated in init routine
type(radial_vel_type), allocatable :: radial_vel_data(:)

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module

! Called once to set values and allocate space

integer :: iunit, io, rc

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_radar_mod_nml", iunit)
read(iunit, nml = obs_def_radar_mod_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_radar_mod_nml")

! Record the namelist values used for the run ... 
if (do_output()) write(nmlfileunit, nml=obs_def_radar_mod_nml)
if (do_output()) write(     *     , nml=obs_def_radar_mod_nml)

! Consistency warning; print a message if the thresholds and lower values
! are going to be used and are different.
call check_namelist_limits(apply_ref_limit_to_obs,     &
                           reflectivity_limit_obs,     &
                           lowest_reflectivity_obs,    &
                           apply_ref_limit_to_fwd_op,  &
                           reflectivity_limit_fwd_op,  &
                           lowest_reflectivity_fwd_op)

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

! Set the module global values that do not change during the run.
! This code uses some values which are set in the namelist, so this call
! *must* happen after the namelist read above.
call initialize_constants()

! Log the values used for the constants.
call print_constants()

end subroutine initialize_module

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Radial velocity section
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine read_radial_vel(velkey, ifile, fform)

! Main read subroutine for the radial velocity observation auxiliary data.

integer,          intent(out)          :: velkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=8)    :: header
logical             :: is_asciifile
type(location_type) :: radar_location
real(r8)            :: beam_direction(3)
real(r8)            :: nyquist_velocity
integer             :: oldkey

! location: Refers to the lat/lon/height/vertical coordinate option
!           for the radar. Its type is defined through location_type in
!           location_mod. Uses the read_location function in
!           location_mod.
!
! beam_direction: This is a 3-element array specific to the radar module:
!              beam_direction(1) = sin(azimuth)*cos(elevation)
!              beam_direction(2) = cos(azimuth)*cos(elevation)
!              beam_direction(3) = sin(elevation)

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
      ! Read the character identifier for verbose formatted output
      read(ifile, FMT="(a8)") header
      if(header /= 'platform') then
         call error_handler(E_ERR,'read_radial_vel', &
              "Expected location header 'platform' in input file", &
              source, revision, revdate)
      endif
endif

! read_location is a dart library routine that expects an optional string
! arg, but the other two read routines are local to this module and we can
! tell them exactly what format to be reading because we already know it.
radar_location   = read_location        (ifile, fform)
beam_direction   = read_beam_direction  (ifile, is_asciifile)
nyquist_velocity = read_nyquist_velocity(ifile, is_asciifile)

! Read in the velkey for this particular observation, however, it will
! be discarded and a new, unique key will be generated in the set routine.
if (is_asciifile) then
   read(ifile, *) oldkey
else
   read(ifile) oldkey
endif

! Generate new unique radial velocity observation key, and set the contents
! of the private defined type.
call set_radial_vel(velkey, radar_location, beam_direction, nyquist_velocity)

end subroutine read_radial_vel

!----------------------------------------------------------------------

subroutine write_radial_vel(velkey, ifile, fform)

! Write radial velocity auxiliary information to the obs_seq file.

integer,          intent(in)           :: velkey, ifile
character(len=*), intent(in), optional :: fform

logical             :: is_asciifile
type(location_type) :: radar_location
real(r8)            :: beam_direction(3)
real(r8)            :: nyquist_velocity

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   ! Write the 5 character identifier for verbose formatted output
   write(ifile, "('platform')")
endif

! Extract the values for this key and call the appropriate write routines.
radar_location    = radial_vel_data(velkey)%radar_location
beam_direction(:) = radial_vel_data(velkey)%beam_direction(:)
nyquist_velocity  = radial_vel_data(velkey)%nyquist_velocity

! write_location routine is part of the dart library and wants the optional
! format string argument.  The other two routines are local to this module, 
! and we have already figured out if it is a unformatted/binary file or 
! formatted/ascii, so go ahead and pass that info directly down to the routines.
call         write_location(ifile, radar_location,    fform) 
call   write_beam_direction(ifile, beam_direction(:), is_asciifile) 
call write_nyquist_velocity(ifile, nyquist_velocity,  is_asciifile) 

! Write out the velkey used for this run, however this will be discarded
! when this observation is read in and a new key will be generated.
! (It may be useful for correlating error messages or identifying particular
! observations so we are leaving it as part of the aux data.)
if (is_asciifile) then
   write(ifile, *) velkey
else
   write(ifile) velkey
endif

end subroutine write_radial_vel

!----------------------------------------------------------------------

function read_beam_direction(ifile, is_asciiformat)

! Reads beam_direction from file that was written by write_beam_direction.
! See read_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_beam_direction(3)

character(len=5)   :: header
real(r8)           :: beam_direction(3)

if ( .not. module_initialized ) call initialize_module


if (is_asciiformat) then
   read(ifile, "(a5)" ) header

   if(header /= 'dir3d') then
      write(msgstring,*)"Expected beam_direction header 'dir3d' in input file, got ", header
      call error_handler(E_ERR, 'read_beam_direction', msgstring, source, revision, revdate)
   endif
   ! Now read the beam_direction data value into temporaries
   read(ifile, *) beam_direction(1), beam_direction(2), beam_direction(3)
else
   ! No header label, just the binary direction values.
   read(ifile)    beam_direction(1), beam_direction(2), beam_direction(3)
endif

! Check for illegal values
if (minval(beam_direction) < -1.0_r8 .or. maxval(beam_direction) > 1.0_r8) then
   write(msgstring,*) "beam_direction value must be between -1 and 1, got: ", &
                       beam_direction(1), beam_direction(2), beam_direction(3)
   call error_handler(E_ERR, 'read_beam_direction', msgstring, &
                      source, revision, revdate)
endif

! set function return value
read_beam_direction(:) = beam_direction(:)

end function read_beam_direction

!----------------------------------------------------------------------

subroutine write_beam_direction(ifile, beam_direction, is_asciiformat)

! Writes beam_direction to obs file.

integer,  intent(in) :: ifile
real(r8), intent(in) :: beam_direction(3)
logical,  intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, "('dir3d')" ) 
   write(ifile, *) beam_direction(1), beam_direction(2), beam_direction(3)
else
   write(ifile)    beam_direction(1), beam_direction(2), beam_direction(3)
endif

end subroutine write_beam_direction

!----------------------------------------------------------------------

function read_nyquist_velocity(ifile, is_asciiformat)

! Reads Nyquist velocity from file that was written by write_nyquist_velocity.
! See read_radial_vel for additional discussion.

integer, intent(in) :: ifile
logical, intent(in) :: is_asciiformat
real(r8)            :: read_nyquist_velocity

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   read(ifile, *) read_nyquist_velocity
else
   read(ifile)    read_nyquist_velocity
endif

! Check for illegal values; must be non-negative.
if (read_nyquist_velocity < 0.0_r8) then
   write(msgstring,*) "bad value for nyquist velocity: ", read_nyquist_velocity
   call error_handler(E_ERR, 'read_nyquist_velocity', msgstring, &
                      source, revision, revdate)
endif
end function read_nyquist_velocity

!----------------------------------------------------------------------

subroutine write_nyquist_velocity(ifile, nyquist_velocity, is_asciiformat)

! Writes Nyquist velocity to obs file.

integer,  intent(in) :: ifile
real(r8), intent(in) :: nyquist_velocity
logical,  intent(in) :: is_asciiformat

if ( .not. module_initialized ) call initialize_module

if (is_asciiformat) then
   write(ifile, *) nyquist_velocity
else
   write(ifile)    nyquist_velocity
endif

end subroutine write_nyquist_velocity

!----------------------------------------------------------------------

subroutine get_obs_def_radial_vel(velkey, radar_location, beam_direction, &
                                  nyquist_velocity)

! Return the auxiliary contents of a given radial velocity observation

integer,             intent(in)  :: velkey
type(location_type), intent(out) :: radar_location
real(r8),            intent(out) :: beam_direction(3)
real(r8),            intent(out) :: nyquist_velocity

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

radar_location    = radial_vel_data(velkey)%radar_location
beam_direction    = radial_vel_data(velkey)%beam_direction(:)
nyquist_velocity  = radial_vel_data(velkey)%nyquist_velocity

end subroutine get_obs_def_radial_vel

!----------------------------------------------------------------------

subroutine set_radial_vel(velkey, radar_location, beam_direction, nyquist_velocity)

! Common code to increment the current key count, and set the private
! contents of this observation's auxiliary data.

integer,             intent(out) :: velkey
type(location_type), intent(in)  :: radar_location
real(r8),            intent(in)  :: beam_direction(3)
real(r8),            intent(in)  :: nyquist_velocity

if ( .not. module_initialized ) call initialize_module

! The total velocity metadata key count from all sequences
velkeycount = velkeycount + 1    
velkey = velkeycount             ! set the return value

! Simple error check on key number before accessing the array
! This errors out if too key value now too large.
call velkey_out_of_range(velkey)

radial_vel_data(velkey)%radar_location    = radar_location
radial_vel_data(velkey)%beam_direction(:) = beam_direction(:)
radial_vel_data(velkey)%nyquist_velocity  = nyquist_velocity

end subroutine set_radial_vel

!----------------------------------------------------------------------

subroutine interactive_radial_vel(velkey)

! Interactively reads in auxiliary information for a radial velocity obs.

integer, intent(out) :: velkey

! Uses interactive_location of DART location_mod, plus the local subroutines
! interactive_beam_direction and set_radial_vel.
! See read_radial_vel for more information.

! velkey is internally incremented in the set routine, and only counts 
! the index for this specialized observation kind. 

type(location_type)  :: location
real(r8)             :: beam_direction(3)
real(r8)             :: nyquist_velocity

if ( .not. module_initialized ) call initialize_module

!Note: Obs location will subsequently be read in by the standard DART module
! "interactive_obs_def". No check is performed here whether radar location,
! obs location, and beam direction are self-consistent. It is the user's
! responsibility to make sure that they are.  This set of information
! does overspecify the problem slightly, but it is expensive to compute
! the true beam direction because of bending of the beam and earth curvature.


write(*, *)
write(*, *) 'Beginning to inquire for information on radar location.'
write(*, *)
write(*, *) 'WARNING!! Make sure that you select 3 (height) for the vertical'
write(*, *) 'coordinate option and enter height in geopotential height (gpm).'
write(*, *) 'This location is where the radar source is located.  The later'
write(*, *) 'location question will be asking about where the observation'
write(*, *) 'itself is located.'
write(*, *)

call interactive_location(location)

write(*, *)
write(*, *) 'Beginning to inquire for information on radar beam direction.'
write(*, *)

call interactive_beam_direction(beam_direction)

write(*, *)
write(*, *) 'Beginning to inquire for information on radar Nyquist velocity.'
write(*, *)

call interactive_nyquist_velocity(nyquist_velocity)


call set_radial_vel(velkey, location, beam_direction, nyquist_velocity)

write(*, *)
write(*, *) 'End of specialized section for radial velocity.'
write(*, *) 'You will now have to enter the regular obs information.'
write(*, *)

end subroutine interactive_radial_vel

!----------------------------------------------------------------------

subroutine interactive_beam_direction(beam_direction)

! Prompt for beam direction information in azimuth/elevation degrees.

real(r8), intent(out) :: beam_direction(3)

real(r8) :: az, el

! FIXME: backwards incompatibility.  when this file is checked in, make sure
! to update the web page about the interface changing.

az = -1.0
do while (az < 0.0 .or. az > 360.0) 
   write(*, *) 'Input the beam direction azimuth in degrees (0 <= az <= 360):'
   read(*, *) az
end do

el = -1.0
do while (el < 0.0 .or. el > 90.0)
   write(*, *) 'Input the beam direction elevation in degrees (0 <= el <= 90):'
   read(*, *) el
end do

! Convert to radians and compute the actual values stored with the observation.
az = az * deg2rad
el = el * deg2rad

beam_direction(1) = sin(az) * cos(el)
beam_direction(2) = cos(az) * cos(el)
beam_direction(3) = sin(el)

end subroutine interactive_beam_direction

!----------------------------------------------------------------------

subroutine interactive_nyquist_velocity(nyquist_velocity)

! Prompt for Nyquist velocity

real(r8), intent(out) :: nyquist_velocity

nyquist_velocity = -1.0

do while (nyquist_velocity < 0.0)
   write(*, *) 'Input Nyquist velocity for this obs point in m/sec'
   write(*, *) '(Typical values are around 20-40, must be >= 0):'
   read(*, *)  nyquist_velocity
end do

end subroutine interactive_nyquist_velocity

!----------------------------------------------------------------------

subroutine get_expected_radial_vel(state_vector, location, velkey, &
                                   radial_vel, istatus)

! This is the main forward operator routine for radar Doppler velocity.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: velkey
real(r8),            intent(out) :: radial_vel
integer,             intent(out) :: istatus


! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted radial velocity that would be observed
! at that location.  The value is returned in 'radial_vel'. 
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

! FIXME: here is original text, which does not quite parse as-is:
! In addition to u,v,w, from which and the direction of the radar beam is
! calculated the along-beam component of the 3-d wind vector, we also need
! the reflectivity value to compute the terminal velocity of the hydrometeor,
! so that w can be adjusted for it. See also get_expected_radar_ref.

! FIXME: here is proposed revised text.  is it accurate?  better?  be honest.
! The along-beam component of the 3-d wind vector is computed from the
! u,v,w fields plus the beam_direction().  The reflectivity value is also
! needed to compute the terminal velocity of the hydrometeor, so the
! w value can be adjusted for it.  

! Reference: Lin et al., 1983 (J. Climate Appl.Meteor., 1065-1092)
! Note that the reflectivity-weighted mean terminal velocities are used here.


real(r8) :: u, v, w, qr, qg, qs, alpha, rho, temp, ref, precip_fall_speed
real(r8) :: precip_r, precip_s, precip_g
real(r8) :: debug_location(3)
logical  :: debug = .false.   ! set to .true. to enable debug printout

if ( .not. module_initialized ) call initialize_module

! Simple error check on key number before accessing the array
call velkey_out_of_range(velkey)

call interpolate(state_vector, location, KIND_U_WIND_COMPONENT, u, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_V_WIND_COMPONENT, v, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_VERTICAL_VELOCITY, w, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO, qr, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, qg, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, qs, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_DENSITY, rho, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif
call interpolate(state_vector, location, KIND_TEMPERATURE, temp, istatus)
if (istatus /= 0) then
   radial_vel = missing_r8
   return
endif

! Based on the values interpolated above, compute the radial velocity
precip_r = rho * qr
precip_s = rho * qs
precip_g = rho * qg
alpha    = sqrt(param_rhos0/rho)

! FIXME: is this the right value to be interpolating, or is it (w-wt) directly?
! FIXME: and what is the actual underlying kind?

! Allow the model to return the terminal velocity directly; otherwise
! compute it based on the number of hydrometeors and the reflectivity.

call interpolate(state_vector, location, KIND_WEIGHTED_VELOCITY, &
                 precip_fall_speed, istatus)
if (istatus /= 0) then

   precip_fall_speed = 0.0_r8
   
   ! Computing reflectivity-weighted precipitation terminal velocity
   ! RAIN
   if (qr >= 1.0e-6_r8) then
      precip_fall_speed = alpha * param_ar_v * (precip_r**param_powr)
   endif
   ! HAIL/GRAUPEL
   if (qg >= 1.0e-6_r8) then
      if (.not. use_wet_graupel .or. temp < 273.15_r8) then
         precip_fall_speed = precip_fall_speed + sqrt(param_e/rho) * &
                             param_ag_dry_v * (precip_g**param_powg_dry)
      else
         precip_fall_speed = precip_fall_speed + sqrt(param_e/rho) * &
                 param_ag_wet_v * ((precip_g/(PI*rho_graupel))**param_powg_wet)
      endif
   endif
   ! SNOW
   if (qs >= 1.0e-6_r8) then
      if ( temp < 273.15_r8 ) then
         precip_fall_speed = precip_fall_speed + alpha * param_as_dry_v * &
                             (precip_s**param_pows)
      else
         precip_fall_speed = precip_fall_speed + alpha * param_as_wet_v * &
                             (precip_s**param_pows)
      endif
   endif
   
   if (precip_fall_speed > 0.0_r8) then
      call get_reflectivity(qr, qg, qs, rho, temp, ref)
      precip_fall_speed = precip_fall_speed/ref
   endif
endif

radial_vel = radial_vel_data(velkey)%beam_direction(1) * u +    &
             radial_vel_data(velkey)%beam_direction(2) * v +    &
             radial_vel_data(velkey)%beam_direction(3) * (w-precip_fall_speed)

! Good return code.  Reset possible istatus error from trying to compute
! weighted fall speed directly.
istatus = 0

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'radial velocity key: ', velkey
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated u: ', u
   print *, 'interpolated v: ', v
   print *, 'interpolated w: ', w
   print *, 'interpolated qr: ', qr
   print *, 'interpolated qg: ', qg
   print *, 'interpolated qs: ', qs
   print *, 'interpolated rho: ', rho
   print *, 'interpolated temp: ', temp
   print *, 'interp or derived fall speed: ', precip_fall_speed
   print *, 'final radial_vel: ', radial_vel
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_radial_vel

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Radar reflectivity
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine read_radar_ref(obsvalue, refkey, ifile, fform)

! Main read subroutine for the obs kind radar reflectivity metadata.

real(r8),         intent(inout)        :: obsvalue
integer,          intent(out)          :: refkey
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform


if ( .not. module_initialized ) call initialize_module

refkey = 0

if ((apply_ref_limit_to_obs) .and. &
    (obsvalue < reflectivity_limit_obs) .and. (obsvalue /= missing_r8)) then
   obsvalue = lowest_reflectivity_obs
endif

end subroutine read_radar_ref

!----------------------------------------------------------------------

subroutine get_expected_radar_ref(state_vector, location, ref, istatus)
 
! The main forward operator routine for radar reflectivity.

real(r8),            intent(in)  :: state_vector(:)
type(location_type), intent(in)  :: location
real(r8),            intent(out) :: ref
integer,             intent(out) :: istatus

! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted radar reflectivity that would be observed
! at that location.

! Radar reflectivity = 10 * log_10( radar reflectivity factor) in dBZ.

! If apply_ref_limit_to_fwd_op, reflectivity values which are below
! reflectivity_limit_fwd_op will be set to lowest_reflectivity_fwd_op.

! "interpolate()" ultimately calls model_mod routine "model_interpolate()"
! to get model values of qr, qg, qs, rho, and temp at the ob location. Then
! the routine "get_reflectivity()" is called to compute the radar reflectivity
! factor value, Z, that corresponds to the hydrometeor and thermodynamic values.

real(r8) :: qr, qg, qs, rho, temp
real(r8) :: debug_location(3)
logical  :: debug = .false.  ! set to .true. to enable debug printout

if ( .not. module_initialized ) call initialize_module

! Start with known values before calling interpolate routines.
qr   = 0.0_r8
qg   = 0.0_r8
qs   = 0.0_r8
rho  = 0.0_r8
temp = 0.0_r8

! If the model can return radar reflectivity data directly, give it a chance
! to do so.  Otherwise, compute the various fields individually and then do
! the computation here.

call interpolate(state_vector, location, KIND_RADAR_REFLECTIVITY, ref, istatus)
if (istatus /= 0) then

   call interpolate(state_vector, location, KIND_RAINWATER_MIXING_RATIO, &
                    qr, istatus)
   if (istatus /= 0) then
      ref = missing_r8
      return
   endif
   
   call interpolate(state_vector, location, KIND_GRAUPEL_MIXING_RATIO, &
                    qg, istatus)
   if (istatus /= 0) then
      ref = missing_r8
      return
   endif
   
   call interpolate(state_vector, location, KIND_SNOW_MIXING_RATIO, &
                    qs, istatus)
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
   
endif

! Always convert to dbz.  Make sure value is always slightly positive.
! (tiny() is a fortran intrinsic function that is > 0 by a very small amount.)
ref = 10.0_r8 * log10(max(tiny(ref), ref))

if ((apply_ref_limit_to_fwd_op) .and. &
    (ref < reflectivity_limit_fwd_op) .and. (ref /= missing_r8)) then
   ref = lowest_reflectivity_fwd_op
endif

! Do not return a missing data value with a successful return code.
if (ref == missing_r8 .and. istatus == 0) then
   istatus = 1
endif

if (debug) then
   debug_location = get_location(location)
   print *
   print *, 'obs location (deg): ', debug_location(1),         &
                                    debug_location(2),         debug_location(3)
   print *, 'obs location (rad): ', debug_location(1)*deg2rad, &
                                    debug_location(2)*deg2rad, debug_location(3)
   print *, 'interpolated qr: ', qr
   print *, 'interpolated qg: ', qg
   print *, 'interpolated qs: ', qs
   print *, 'interpolated rho: ', rho
   print *, 'interpolated temp: ', temp
   print *, 'final reflectivity: ', ref
   print *, 'istatus: ', istatus
endif

end subroutine get_expected_radar_ref

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Helper routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine get_reflectivity(qr, qg, qs, rho, temp, ref)
 
! Computes "radar reflectivity factor" (Z) in mm^6 m^-3

real(r8), intent(in)  :: qr, qg, qs, rho, temp
real(r8), intent(out) :: ref

!
! References: Ferrier, 1994 (JAS, 249-280)
!             Smith et al., 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)

real(r8) :: precip

if ( .not. module_initialized ) call initialize_module

ref = 0.0_r8

! RAIN
if ( qr >= 1.0e-6_r8 ) then
   precip = rho * qr
   ref = ref + param_ar_f * (precip**1.75_r8)
endif

! HAIL / GRAUPEL
if ( qg >= 1.0e-6_r8 ) then
    precip = rho * qg
    if (.not. use_wet_graupel .or. temp < 273.15_r8) then
       ref = ref + param_ag_dry_f * (precip**1.75_r8)
    else
       ref = ref + param_ag_wet_f * (precip**1.6625_r8)
    endif
endif

! SNOW
if ( qs >= 1.0e-6_r8 ) then
   precip = rho * qs
   if ( temp < 273.15_r8 ) then
      ref = ref + param_as_dry_f * (precip**1.75_r8)
   else
      ref = ref + param_as_wet_f * (precip**1.75_r8)
   endif
endif

end subroutine get_reflectivity

!----------------------------------------------------------------------

subroutine initialize_constants()

! Initialize module global constants. 

! IMPORTANT: Uses namelist values, so this routine cannot be called until
! after the namelist has been read.
!

! FIXME: gravity is currently hardcoded here.  We should find a way to let
! the model have input if it uses a slightly different value for G, or if it
! is working on a different planet with an entirely different set of constants.
! Question:  how much impact on the results does changing G have?

! References: Ferrier, 1994 (JAS, 249-280)
!             Smith et al., 1984 (JCAM, 1258-1260)
!             Smith at al., 1975 (JAM, 1156-1165)


param_gravity  = 9.81_r8        ! In general, this must match model value.
param_a        = 8.42e20_r8
param_b        = 0.8_r8
param_c        = 4.84e18_r8
param_d        = 0.25_r8
param_CD       = 0.6_r8
param_rhos0    = 1.0_r8
param_e        = 4.0_r8*param_gravity*rho_graupel/(3.0_r8*param_CD)
param_f        = 0.5_r8

param_gam7b    = 3376.92_r8
param_gam7d    = 1155.38_r8
param_gam7f    = 1871.25_r8
param_powr     = (7.0_r8 + param_b)/4.0_r8
param_pows     = (7.0_r8 + param_d)/4.0_r8
param_powg_dry = (7.0_r8 + param_f)/4.0_r8
param_powg_wet = 1.7875_r8


! There are two different computations for the following constants.
! _v is the radial velocity version, _f is the reflectivity one.

! radial velocity code:
param_ar_v     = n0_rain * param_a * param_gam7b / &
                 (PI * rho_rain * n0_rain)**param_powr

param_as_wet_v = n0_snow * param_c * param_gam7d / &
                 (PI * rho_snow * n0_snow)**param_pows

param_as_dry_v = dielectric_factor * ((rho_snow/rho_rain)**2.0_r8) * &
                 param_as_wet_v

param_ag_wet_v = ((7.2e20_r8)**0.95_r8) * param_gam7f / &
                 (720.0_r8 * (n0_graupel**0.8375_r8))

param_ag_dry_v = 1.0e18_r8 * dielectric_factor *                               &
                 ((rho_graupel/rho_rain)**2.0_r8) * n0_graupel * param_gam7f / &
                 (PI * rho_graupel * n0_graupel)**param_powg_dry

! reflectivity code:
param_ar_f     = 7.2e20_r8 / (((PI*rho_rain)**1.75_r8)*(n0_rain**0.75_r8))

param_as_wet_f = 7.2e20_r8 / (((PI*rho_snow)**1.75_r8)*(n0_snow**0.75_r8))

param_as_dry_f = dielectric_factor * ((rho_snow/rho_rain)**2.0_r8) * & 
                 param_as_wet_f

param_ag_wet_f = (7.2e20_r8/(((PI*rho_graupel)**1.75_r8) * &
                 (n0_graupel**0.75_r8)))**0.95_r8

param_ag_dry_f = dielectric_factor * ((rho_graupel/rho_rain)**2.0_r8) * &
                 7.2e20_r8 / (((PI*rho_graupel)**1.75_r8)*(n0_graupel**0.75_r8))



end subroutine initialize_constants

!----------------------------------------------------------------------

subroutine print_constants()

! Log the constants set in the code.
! Prints to both the log file and standard output.

! The values in this list which are also in the namelist will have their
! values written by the write(nml=) code, but this routine includes all
! the fixed constants so they are written in one place, both to standard 
! output and the log file. Using the correct values is critical to doing 
! the appropriate computation, so some duplication is probably a good thing.

write(msgstring, *) 'Constants used in the obs_def_radar module:'
call error_handler(E_MSG,'', msgstring, '', '', '')

call pr_con(dielectric_factor , "dielectric_factor" )
call pr_con(n0_rain           , "n0_rain"           )
call pr_con(n0_graupel        , "n0_graupel"        )
call pr_con(n0_snow           , "n0_snow"           )
call pr_con(rho_rain          , "rho_rain"          )
call pr_con(rho_graupel       , "rho_graupel"       )
call pr_con(rho_snow          , "rho_snow"          )
call pr_con(param_gravity     , "param_gravity"     )
call pr_con(param_a           , "param_a"           )
call pr_con(param_b           , "param_b"           )
call pr_con(param_c           , "param_c"           )
call pr_con(param_d           , "param_d"           )
call pr_con(param_CD          , "param_CD"          )
call pr_con(param_rhos0       , "param_rhos0"       )
call pr_con(param_e           , "param_e"           )
call pr_con(param_f           , "param_f"           )
call pr_con(param_gam7b       , "param_gam7b"       )
call pr_con(param_gam7d       , "param_gam7d"       )
call pr_con(param_gam7f       , "param_gam7f"       )
call pr_con(param_powr        , "param_powr"        )
call pr_con(param_pows        , "param_pows"        )
call pr_con(param_powg_dry    , "param_powg_dry"    )
call pr_con(param_powg_wet    , "param_powg_wet"    )
call pr_con(param_ar_v        , "param_ar_v"        )
call pr_con(param_as_wet_v    , "param_as_wet_v"    )
call pr_con(param_as_dry_v    , "param_as_dry_v"    )
call pr_con(param_ag_wet_v    , "param_ag_wet_v"    )
call pr_con(param_ag_dry_v    , "param_ag_dry_v"    )
call pr_con(param_ar_f        , "param_ar_f"        )
call pr_con(param_as_wet_f    , "param_as_wet_f"    )
call pr_con(param_as_dry_f    , "param_as_dry_f"    )
call pr_con(param_ag_wet_f    , "param_ag_wet_f"    )
call pr_con(param_ag_dry_f    , "param_ag_dry_f"    )

end subroutine print_constants

!----------------------------------------------------------------------

subroutine pr_con(c_val, c_str)

! Utility routine to print a string and value

real(r8),         intent(in) :: c_val
character(len=*), intent(in) :: c_str

write(msgstring, "(A30,A,ES28.8)") c_str, " = ", c_val
call error_handler(E_MSG,'', msgstring, '', '', '')

end subroutine pr_con

!----------------------------------------------------------------------

function ascii_file_format(fform)

! Common routine for determining input file format.

character(len=*), intent(in), optional :: fform
logical                                :: ascii_file_format

! Returns .true. if file is formatted/ascii, .false. if unformatted/binary
! Defaults (if fform not specified) to formatted/ascii.

if ( .not. module_initialized ) call initialize_module

! Default to formatted/ascii.
if ( .not. present(fform)) then
   ascii_file_format = .true.
   return 
endif

SELECT CASE (trim(adjustl(fform)))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      ascii_file_format = .false.
   CASE DEFAULT
      ascii_file_format = .true.
END SELECT


end function ascii_file_format

!----------------------------------------------------------------------

subroutine velkey_out_of_range(velkey)

! Range check velkey and trigger a fatal error if larger than allocated array.

integer, intent(in) :: velkey

! fine -- no problem.
if (velkey <= max_radial_vel_obs) return

! Bad news.  Tell the user.
write(msgstring, *) 'velkey (',velkey,') exceeds max_radial_vel_obs (', &
                     max_radial_vel_obs,')'
call error_handler(E_MSG,'set_radial_vel', msgstring, '', '', '')
call error_handler(E_ERR,'set_radial_vel', &
                   'Increase max_radial_vel_obs in namelist', &
                   source, revision, revdate)

end subroutine velkey_out_of_range

!----------------------------------------------------------------------

subroutine check_namelist_limits(apply_ref_limit_to_obs, &
   reflectivity_limit_obs, lowest_reflectivity_obs, apply_ref_limit_to_fwd_op,& 
   reflectivity_limit_fwd_op, lowest_reflectivity_fwd_op)

! Consistency warning; print a message if the thresholds and lower values
! are going to be used and are different.

logical,  intent(in) :: apply_ref_limit_to_obs
real(r8), intent(in) :: reflectivity_limit_obs
real(r8), intent(in) :: lowest_reflectivity_obs
logical,  intent(in) :: apply_ref_limit_to_fwd_op
real(r8), intent(in) :: reflectivity_limit_fwd_op
real(r8), intent(in) :: lowest_reflectivity_fwd_op

! FIXME: the point is to gently remind the user if they are setting different
! limits on the actual observation values and the forward operator, since
! that may be what they intend, but it isn't something they should be doing
! by mistake.  but we don't want to be annoying.  sanity check this code
! carefully to be sure it is at least accurate.

! if neither limit is being enforced, return silently.
if (.not. apply_ref_limit_to_obs .and. .not. apply_ref_limit_to_fwd_op) return

! if both are on, and the cutoffs and clamp-to values are the same, fine also.
if (apply_ref_limit_to_obs .and. apply_ref_limit_to_fwd_op) then
   if ((reflectivity_limit_obs  == reflectivity_limit_fwd_op) .and. &
       (lowest_reflectivity_obs == lowest_reflectivity_fwd_op)) return
endif

! ok, if we got this far, either one of the limits is not on, and/or
! the limits or clamp-to values do not match.  don't make this an error,
! but do print something to the log to point out they aren't the same.
if (apply_ref_limit_to_obs) then
   write(msgstring, *) 'reflectivity obs values below ', &
      reflectivity_limit_obs, 'will be set to', lowest_reflectivity_obs
else
   write(msgstring, *) 'reflectivity obs values will be processed unchanged'
endif
call error_handler(E_MSG,'check_namelist_limits', msgstring, '', '', '')

if (apply_ref_limit_to_fwd_op) then
   write(msgstring, *) 'reflectivity forward operator values below ', &
      reflectivity_limit_fwd_op, 'will be set to', lowest_reflectivity_fwd_op
else
   write(msgstring, *) 'reflectivity forward operator values will be processed unchanged'
endif
call error_handler(E_MSG,'check_namelist_limits', msgstring, '', '', '')

end subroutine check_namelist_limits

!----------------------------------------------------------------------

end module obs_def_radar_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
