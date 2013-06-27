! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Banner for checking maximum lengths (i.e. 32)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! BEGIN DART PREPROCESS KIND LIST
! AMSRE_BRIGHNTESS_T,             KIND_BRIGHTNESS_TEMPERATURE
! END DART PREPROCESS KIND LIST

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_brightnessT_mod, only : get_brightness_temperature, &
!                                      read_amsre_metadata, &
!                                      write_amsre_metadata, &
!                                      interactive_amsre_metadata
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(AMSRE_BRIGHNTESS_T)
!     call get_brightness_temperature(state, state_time, ens_index, location, obs_time, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!  case(AMSRE_BRIGHNTESS_T)
!     call read_amsre_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!  case(AMSRE_BRIGHNTESS_T)
!     call write_amsre_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!  case(AMSRE_BRIGHNTESS_T)
!     call interactive_amsre_metadata(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_brightnessT_mod

! This is the prototype forward observation operator code for brightness temperatures.
! The DART model state for CLM generally does not have all the pieces necessary
! to apply the forward observation operator directly, so this code gets what it
! needs from a CLM history file. These history files have become complicated now
! that CLM is trying to support unstructured grids. Sometimes the variables of
! interest are shaped XXX(time, lat, lon), sometimes XXX(time, lndgrid).
! 'single column' runs may appear as either lat=lon=1 or lndgrid=1
!
use        types_mod, only : r4, r8, digits12, MISSING_R8, MISSING_I, PI, deg2rad
use     location_mod, only : location_type, get_location, get_dist, &
                             set_location, VERTISUNDEF
use time_manager_mod, only : time_type, get_date, set_date, print_date, print_time, &
                             get_time, set_time, operator(-), operator(/=)
use    utilities_mod, only : register_module, E_ERR, E_MSG, error_handler,     &
                             check_namelist_read, find_namelist_in_file,       &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             nc_check, file_exist, is_longitude_between,       &
                             ascii_file_format
use typesizes
use netcdf

implicit none
private

public ::  set_amsre_metadata, &
           get_amsre_metadata, &
          read_amsre_metadata, &
         write_amsre_metadata, &
   interactive_amsre_metadata, &
   get_brightness_temperature

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical            :: module_initialized = .false.
logical            :: unstructured = .false.
character(len=129) :: string1, string2, string3
integer            :: nlon, nlat, ntime, ens_size
type(time_type)    :: initialization_time

character(len=129), allocatable, dimension(:) :: fname
integer,            allocatable, dimension(:) :: ncid
real(r8),           allocatable, dimension(:) :: lon, lat, area
real(digits12),     allocatable, dimension(:) :: rtime

real(r8), parameter :: RAD2KM = 40030.0_r8/(2.0_r8 * PI) ! (mean radius of earth ~6371km)

!----------------------------------------------------------------------
! Metadata for AMSR-E observations.
!----------------------------------------------------------------------

type amsre_metadata
   private
   real(r8)  :: frequency    ! gHz
   real(r8)  :: footprint    ! km^2  'support' of the observation
   character :: polarization ! vertical or horizontal
   integer   :: landcovercode
end type amsre_metadata

type(amsre_metadata), allocatable, dimension(:) :: observation_metadata
type(amsre_metadata) :: missing_metadata

character(len=8), parameter :: AMSRESTRING = 'amsr-e'
character(len=8), parameter ::  IGBPSTRING = 'igbp'

integer :: MAXamsrekey = 24*366  ! one year of hourly data - to start
integer ::    amsrekey = 0       ! useful length of metadata arrays

!----------------------------------------------------------------------
! namelist items
!----------------------------------------------------------------------

character(len=256) :: casename = 'clm_dart'
logical            :: debug = .false.
integer            :: hist_nhtfrq = -24 
! CLM variable hist_nhtfrq ... controls how often to write out the history files.
! Negative value means the output frequency is the absolute value (in hours).

namelist /obs_def_brightnessT_nml/ casename, debug, hist_nhtfrq

!----------------------------------------------------------------------
! This function name will be a problem if cosmos and amsrE are needed 
! at the same time. refactor into single function in utilities_mod?
!----------------------------------------------------------------------

!interface interactive
!   module procedure interactive_real
!   module procedure interactive_int
!end interface

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of public routines for obs_def_brightnessT_mod
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine set_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

! Fill the module storage metadata for a particular observation.

integer,   intent(out) :: key
real(r8),  intent(in)  :: frequency, footprint
character, intent(in)  :: polarization
integer,   intent(in)  :: landcovercode

if ( .not. module_initialized ) call initialize_module

amsrekey = amsrekey + 1  ! increase module storage used counter

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(amsrekey,'set_amsre_metadata')

key = amsrekey ! now that we know its legal

observation_metadata(key)%frequency     = frequency
observation_metadata(key)%footprint     = footprint
observation_metadata(key)%polarization  = polarization
observation_metadata(key)%landcovercode = landcovercode

end subroutine set_amsre_metadata


!======================================================================


subroutine get_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

! Query the metadata in module storage for a particular observation.
! This can be useful for post-processing routines, etc.

integer,   intent(in)  :: key
real(r8),  intent(out) :: frequency, footprint
character, intent(out) :: polarization
integer,   intent(out) :: landcovercode

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_amsre_metadata')

frequency     = observation_metadata(key)%frequency
footprint     = observation_metadata(key)%footprint
polarization  = observation_metadata(key)%polarization
landcovercode = observation_metadata(key)%landcovercode

end subroutine get_amsre_metadata


!======================================================================


 subroutine read_amsre_metadata(key,       obsID, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_amsre_metadata(obs_def%key, key, ifile, fform)
!
! This routine reads the metadata for neutron intensity observations.
!
integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
integer           :: ierr
character(len=9)  :: header

real(r8)          :: frequency, footprint
character         :: polarization
integer           :: landcovercode

if ( .not. module_initialized ) call initialize_module

write(string2,*)'observation #',obsID

if ( ascii_file_format(fform) ) then

   read(ifile, *, iostat=ierr) header
   call check_iostat(ierr,'read_amsre_metadata','header',string2)
   if (trim(header) /= trim(AMSRESTRING)) then
      write(string1,*)"Expected Tb header ["//AMSRESTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, *, iostat=ierr) frequency, footprint, polarization
   call check_iostat(ierr,'read_amsre_metadata','freq foot polar',string2)

   read(ifile, *, iostat=ierr) header
   if (trim(header) /= trim(IGBPSTRING)) then
      write(string1,*)"Expected IGBP header ["//IGBPSTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, *, iostat=ierr) landcovercode
   call check_iostat(ierr,'read_amsre_metadata','landcovercode',string2)

else

   read(ifile, iostat=ierr) header
   call check_iostat(ierr,'read_amsre_metadata','header',string2)
   if (trim(header) /= trim(AMSRESTRING)) then
      write(string1,*)"Expected Tb header ["//AMSRESTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, iostat=ierr) frequency, footprint, polarization
   call check_iostat(ierr,'read_amsre_metadata','freq foot polar',string2)

   read(ifile, iostat=ierr) header
   if (trim(header) /= trim(IGBPSTRING)) then
      write(string1,*)"Expected IGBP header ["//IGBPSTRING//"] in input file, got ["//header//"]"
      call error_handler(E_ERR,'read_amsre_metadata',string1,source,revision,revdate,text2=string2)
   endif
   read(ifile, iostat=ierr) landcovercode
   call check_iostat(ierr,'read_amsre_metadata','landcovercode',string2)

endif

! Store the metadata in module storage and record the new length of the metadata arrays.
call set_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

! The new 'key' is returned.

end subroutine read_amsre_metadata


!======================================================================


 subroutine write_amsre_metadata(key, ifile, fform)
!----------------------------------------------------------------------
! writes the metadata for AMSR-E brightness temperature observations.

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform

real(r8)          :: frequency, footprint
character         :: polarization
integer           :: landcovercode

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

call get_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

if ( ascii_file_format(fform)) then
   write(ifile, *) trim(AMSRESTRING)
   write(ifile, '(f8.2,2x,f8.2,2x,A)') frequency, footprint, polarization
   write(ifile, *) trim(IGBPSTRING)
   write(ifile, '(i8)') landcovercode
else
   write(ifile) trim(AMSRESTRING)
   write(ifile) frequency, footprint, polarization
   write(ifile) trim(IGBPSTRING)
   write(ifile) landcovercode
endif

end subroutine write_amsre_metadata


!======================================================================
! The AMSR-E Level-2A product (AE_L2A) contains brightness temperatures at 
! 6.9 GHz, 10.7 GHz, 18.7 GHz, 23.8 GHz, 36.5 GHz, and 89.0 GHz. 
! Data are resampled to be spatially consistent, and therefore are available 
! at a variety of resolutions that correspond to the footprint sizes of the 
! observations such as 56 km, 38 km, 21 km, 12 km, and 5.4 km, respectively. 
! Each swath is packaged with associated geolocation fields. 
! Data are stored in Hierarchical Data Format - Earth Observing System (HDF-EOS) 
! format and are available from 1 June 2002 to 4 October 2011 via FTP.

subroutine interactive_amsre_metadata( key )

integer, intent(out) :: key

real(r8)  :: frequency, footprint
character :: polarization
integer   :: landcovercode

if ( .not. module_initialized ) call initialize_module

! Prompt for input for the required metadata

frequency     = interactive_real('"frequency"    [GHz]',minvalue = 0.0_r8)
footprint     = interactive_real('"footprint"    [km]' ,minvalue = 0.0_r8)
polarization  = interactive_char('"polarization" [H,V]')
landcovercode = interactive_int( '"IGBP land cover code" [??]',minvalue = 0)

call set_amsre_metadata(key, frequency, footprint, polarization, landcovercode)

end subroutine interactive_amsre_metadata


!======================================================================


subroutine get_brightness_temperature(state, state_time, ens_index, location, obs_time, obs_val, istatus)

! This is THE forward observation operator. Given the state and a location, return the value

real(r8),            intent(in)  :: state(:)   ! DART state vector
type(time_type),     intent(in)  :: state_time ! valid time of DART state 
integer,             intent(in)  :: ens_index  ! Ensemble member number
type(location_type), intent(in)  :: location   ! target location (i.e. obs)
type(time_type),     intent(in)  :: obs_time   ! time of observation
real(r8),            intent(out) :: obs_val    ! model estimate of observation value
integer,             intent(out) :: istatus    ! status of the calculation

istatus = 1
obs_val = MISSING_R8

! a wee bit of work left ...

end subroutine get_brightness_temperature


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of internal routines for obs_def_brightnessT_mod
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine initialize_module

! Called once to set values and allocate space, open all the CLM files
! that have the observations, etc.

integer :: iunit, io, i

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_brightnessT_nml", iunit)
read(iunit, nml = obs_def_brightnessT_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_brightnessT_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_def_brightnessT_nml)
if (do_nml_term()) write(     *     , nml=obs_def_brightnessT_nml)

! Allocate the module array to store the key information.

missing_metadata%frequency     = MISSING_R8
missing_metadata%footprint     = MISSING_R8
missing_metadata%polarization  = 'x'
missing_metadata%landcovercode = MISSING_I

allocate( observation_metadata(MAXamsrekey) )
observation_metadata(:) = missing_metadata

end subroutine initialize_module


!======================================================================


subroutine GetDimensions(ncid, fname)

! Harvest information from the first observation file.
! The SingleColumMode files have
!        float lat(lndgrid) ;
!        float lon(lndgrid) ;
!        float area(lndgrid) ;
! while the 2D files have
!        float lat(lat) ;
!        float lon(lon) ;
!        float area(lat, lon) ;

integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname

! integer, intent(out) :: nlon, nlat, ntime ... module variables

integer, dimension(NF90_MAX_VAR_DIMS) :: londimids, latdimids
integer :: lonvarid, lonndims
integer :: latvarid, latndims
integer :: dimid

call nc_check(nf90_inq_varid(ncid,  'lon', lonvarid), &
              'obs_def_brightnessT.GetDimensions','inq_varid lon '//trim(fname))
call nc_check(nf90_inq_varid(ncid,  'lat', latvarid), &
              'obs_def_brightnessT.GetDimensions','inq_varid lat '//trim(fname))

call nc_check(nf90_inquire_variable(ncid, lonvarid, ndims=lonndims, dimids=londimids),&
              'obs_def_brightnessT.GetDimensions','inquire lon '//trim(fname))
call nc_check(nf90_inquire_variable(ncid, latvarid, ndims=latndims, dimids=latdimids),&
              'obs_def_brightnessT.GetDimensions','inquire lat '//trim(fname))

if ( (lonndims /= 1) .or. (latndims /= 1) ) then
   write(string1,*) 'Require "lon" and "lat" variables to be 1D. They are ', &
                     lonndims, latndims
   call error_handler(E_ERR,'obs_def_brightnessT.GetDimensions',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, londimids(1), len=nlon), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension lon '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, latdimids(1), len=nlat), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension lat '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'time', dimid), &
              'obs_def_brightnessT.GetDimensions','inq_dimid time '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=ntime), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension time '//trim(fname))

if ( nf90_inq_dimid(ncid, 'lndgrid', dimid) == NF90_NOERR) then
   unstructured = .true.
   if ((nlon /= 1) .or. (nlat /= 1)) then
      string1 = 'unstructured grids with more than a single gridcell are not supported.'
      call error_handler(E_ERR,'obs_def_brightnessT.GetDimensions',string1,source,revision,revdate)
   endif
endif

end subroutine GetDimensions


!======================================================================


subroutine get_scalar_from_history(varstring, state_time, ens_index, location, &
                                   obs_time, obs_val, istatus)

character(len=*),    intent(in)  :: varstring
type(time_type),     intent(in)  :: state_time
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

if ( unstructured ) then
   call get_scalar_from_2Dhistory(varstring,ens_index,location,obs_time,obs_val,istatus)
else
   call get_scalar_from_3Dhistory(varstring,ens_index,location,obs_time,obs_val,istatus)
endif

end subroutine get_scalar_from_history


!======================================================================


subroutine get_scalar_from_3Dhistory(varstring, ens_index, location, obs_time, &
                                     obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use
!
! The requirement is that the history file variable is a 3D variable shaped similarly:
!
! float NEP(time, lat, lon) ;
!          NEP:long_name = "net ecosystem production, blah, blah, blah" ;
!          NEP:units = "gC/m^2/s" ;
!          NEP:cell_methods = "time: mean" ;
!          NEP:_FillValue = 1.e+36f ;
!          NEP:missing_value = 1.e+36f ;

character(len=*),    intent(in)  :: varstring
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer,  dimension(NF90_MAX_VAR_DIMS) :: dimids
real(r8), dimension(3) :: loc
integer,  dimension(3) :: ncstart, nccount
integer,  dimension(1) :: loninds, latinds, timeinds
integer                :: gridloni, gridlatj, timei
integer                :: varid, xtype, ndims, natts, dimlen
integer                :: io1, io2, second, day
real(r8)               :: loc_lon, loc_lat, radius, distance
real(r4), dimension(1) :: hyperslab
real(r4)               :: spvalR4
real(r8)               :: scale_factor, add_offset
real(digits12)         :: otime
character(len=NF90_MAX_NAME+20)      :: strshort

type(location_type) :: gridloc

obs_val = MISSING_R8
istatus = 1

!----------------------------------------------------------------------
! if observation is outside region encompassed in the history file - fail
loc      = get_location(location) ! loc is in DEGREES
loc_lon  = loc(1)
loc_lat  = loc(2)

if ( (nlon==1) .and. (nlat==1) ) then

   ! Defining the region if running in an unstructured grid is tricky.
   ! Have lat, lon, and the area of the gridcell which we assume to be basically square.
   ! The square root of the area defines the length of the edge of the gridcell.
   ! Half the hypotenuse defines the radius of a circle. Any ob within
   ! that radius is close enough.

   gridloc   = set_location(lon(1),lat(1), 0.0_r8, VERTISUNDEF)
   distance  = get_dist(gridloc, location, no_vert = .TRUE.) * RAD2KM ! planet earth
   radius    = sqrt(2.0_r8 * area(1))/2.0_r8

   if (debug .and. do_output()) then
      write(string1,*)'    observation lon, lat is ',loc_lon, loc_lat
      write(string2,*)'gridcell    lon, lat is ',lon(1),lat(1)
      write(string3,*)'area,radius is ',area(1),radius,' distance ',distance
      call error_handler(E_MSG, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
                 string1, source, revision, revdate, text2=string2, text3=string3)
   endif

   if ( distance > radius ) return

else
   if ( .not. is_longitude_between(loc_lon, lon(1), lon(nlon), doradians=.FALSE.)) return
   if ((loc_lat < lat(1)) .or. (loc_lat > lat(nlat))) return
endif

!----------------------------------------------------------------------
! Now that we know the observation operator is possible, continue ...

write(strshort,'(''ens_index '',i4,1x,A)')ens_index,trim(varstring)

if (ens_index > ens_size) then
   write(string1,*)'Known to have ',ens_size,'ensemble members for observation operator.'
   write(string2,*)'asking to use operator for ensemble member ',ens_index
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate, text2=string2)
endif

!----------------------------------------------------------------------
! bombproofing ... make sure the netcdf file is open.

call nc_check(nf90_inquire(ncid(ens_index)), &
              'obs_def_brightnessT.get_scalar_from_3Dhistory', 'inquire '//trim(strshort))

! bombproofing ... make sure the variable is the shape and size we expect

call nc_check(nf90_inq_varid(ncid(ens_index), trim(varstring), varid), &
        'obs_def_brightnessT.get_scalar_from_3Dhistory', 'inq_varid '//trim(strshort))
call nc_check(nf90_inquire_variable(ncid(ens_index), varid, xtype=xtype, ndims=ndims, &
        dimids=dimids, natts=natts), &
        'obs_def_brightnessT.get_scalar_from_3Dhistory','inquire variable '//trim(strshort))

if (ndims /= 3) then
   write(string1,*)trim(varstring),' is supposed to have 3 dimensions, it has',ndims
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! If the variable is not a NF90_FLOAT, then the assumptions for processing
! the missing_value, _FillValue, etc., may not be correct.
if (xtype /= NF90_FLOAT) then
   write(string1,*)trim(varstring),' is supposed to be a 32 bit real. xtype = ', &
                   NF90_FLOAT,' it is ',xtype
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 1 is longitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(1), len=dimlen), &
        'obs_def_brightnessT.get_scalar_from_3Dhistory', 'inquire_dimension 1 '//trim(strshort))
if (dimlen /= nlon) then
   write(string1,*)'LON has length',nlon,trim(varstring),' has ',dimlen,'longitudes.'
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 2 is latitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(2), len=dimlen), &
        'obs_def_brightnessT.get_scalar_from_3Dhistory', 'inquire_dimension 2 '//trim(strshort))
if (dimlen /= nlat) then
   write(string1,*)'LAT has length',nlat,trim(varstring),' has ',dimlen,'latitudes.'
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 3 is time
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(3), len=dimlen), &
        'obs_def_brightnessT.get_scalar_from_3Dhistory', 'inquire_dimension 3'//trim(strshort))
if (dimlen /= ntime) then
   write(string1,*)'TIME has length',ntime,trim(varstring),' has ',dimlen,'times.'
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

!----------------------------------------------------------------------
! Find the grid cell and timestep of interest
! FIXME ... since the history file contents are for the previous 30m,
! perhaps the closest time is not the best approximation.
! Get the individual locations values

call get_time(obs_time, second, day)
otime    = real(day,digits12) + real(second,digits12)/86400.0_digits12

latinds  = minloc(abs(lat - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(lon - loc_lon))   ! these return 'arrays' ...
timeinds = minloc(abs(rtime - otime))   ! these return 'arrays' ...

gridlatj = latinds(1)
gridloni = loninds(1)
timei    = timeinds(1)

if (debug .and. do_output()) then
   write(*,*)'obs_def_brightnessT.get_scalar_from_3Dhistory:targetlon, lon, lon index is ', &
                                           loc_lon,lon(gridloni),gridloni
   write(*,*)'obs_def_brightnessT.get_scalar_from_3Dhistory:targetlat, lat, lat index is ', &
                                           loc_lat,lat(gridlatj),gridlatj
   write(*,*)'obs_def_brightnessT.get_scalar_from_3Dhistory:  targetT,   T,   T index is ', &
                                           otime,rtime(timei),timei
endif

if ( abs(otime - rtime(timei)) > 30*60 ) then
   if (debug .and. do_output()) then
      write(*,*)'obs_def_brightnessT.get_scalar_from_3Dhistory: no close time ... skipping observation'
      call print_time(obs_time,'obs_def_brightnessT.get_scalar_from_3Dhistory:observation time')
      call print_date(obs_time,'obs_def_brightnessT.get_scalar_from_3Dhistory:observation date')
   endif
   istatus = 2
   return
endif

!----------------------------------------------------------------------
! Grab exactly the scalar we want.

ncstart = (/ gridloni, gridlatj, timei /)
nccount = (/        1,        1,     1 /)

call nc_check(nf90_get_var(ncid(ens_index),varid,hyperslab,start=ncstart,count=nccount), &
     'obs_def_brightnessT.get_scalar_from_3Dhistory', 'get_var')

obs_val = hyperslab(1)

!----------------------------------------------------------------------
! Apply any netCDF attributes ...

io1 = nf90_get_att(ncid(ens_index), varid, '_FillValue' , spvalR4)
if ((io1 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io2 = nf90_get_att(ncid(ens_index), varid, 'missing_value' , spvalR4)
if ((io2 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io1 = nf90_get_att(ncid(ens_index), varid, 'scale_factor', scale_factor)
io2 = nf90_get_att(ncid(ens_index), varid, 'add_offset'  , add_offset)

if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor + add_offset
elseif (io1 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor
elseif (io2 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val + add_offset
endif

if (obs_val /= MISSING_R8) istatus = 0

end subroutine get_scalar_from_3Dhistory


!======================================================================


subroutine get_scalar_from_2Dhistory(varstring, ens_index, location, obs_time, &
                                     obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use
!
! The requirement is that the history file variable is a 2D variable shaped similarly:
!
! float NEP(time, lndgrid) ;
!          NEP:long_name = "net ecosystem production, blah, blah, blah" ;
!          NEP:units = "gC/m^2/s" ;
!          NEP:cell_methods = "time: mean" ;
!          NEP:_FillValue = 1.e+36f ;
!          NEP:missing_value = 1.e+36f ;
!
! Just because it is 2D does not mean it is a single column,
! although single columns are all that is really supported right now.

character(len=*),    intent(in)  :: varstring
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer,  dimension(NF90_MAX_VAR_DIMS) :: dimids
real(r8), dimension(3) :: loc
integer,  dimension(2) :: ncstart, nccount
integer,  dimension(1) :: timeinds
integer                :: gridij, timei
integer                :: varid, xtype, ndims, natts, dimlen
integer                :: io1, io2, second, day
real(r8)               :: loc_lon, loc_lat, radius, distance
real(r4), dimension(1) :: hyperslab
real(r4)               :: spvalR4
real(r8)               :: scale_factor, add_offset
real(digits12)         :: otime
character(len=NF90_MAX_NAME+20)      :: strshort

type(location_type) :: gridloc

obs_val = MISSING_R8
istatus = 1

!----------------------------------------------------------------------
! if observation is outside region encompassed in the history file - fail
loc      = get_location(location) ! loc is in DEGREES
loc_lon  = loc(1)
loc_lat  = loc(2)

! Defining the region if running in an unstructured grid is tricky.
! We have lat, lon, and the area of the gridcell which we assume to be basically square.
! The square root of the area defines the length of the edge of the gridcell.
! Half the hypotenuse defines the radius of a circle. Any ob within
! that radius is close enough.

! TJH FIXME This does not work with unstructured grid.
! latinds  = minloc(abs(lat - loc_lat))   ! these return 'arrays' ...
! loninds  = minloc(abs(lon - loc_lon))   ! these return 'arrays' ...
! gridij = the closest location

! The "1" in the following reflect the fact that only a single gridcell
! is currently supported in the unstructured grid configuration.
! (see GetDimensions())

gridij   = 1
gridloc  = set_location(lon(gridij),lat(gridij), 0.0_r8, VERTISUNDEF)
distance = get_dist(gridloc, location, no_vert = .TRUE.) * RAD2KM
radius   = sqrt(2.0_r8 * area(gridij))/2.0_r8

if (debug .and. do_output()) then
   write(string1,*)'    observation lon, lat is ',loc_lon, loc_lat
   write(string2,*)'gridcell    lon, lat is ',lon(gridij),lat(gridij)
   write(string3,*)'area,radius is ',area(gridij),radius,' distance ',distance
   call error_handler(E_MSG, 'obs_def_brightnessT.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate, text2=string2, text3=string3)
endif

if ( distance > radius ) return

!----------------------------------------------------------------------
! Now that we know the observation operator is possible, continue ...

write(strshort,'(''ens_index '',i4,1x,A)')ens_index,trim(varstring)

if (ens_index > ens_size) then
   write(string1,*)'believed to have ',ens_size,'ensemble members for observation operator.'
   write(string2,*)'asking to use operator for ensemble member ',ens_index
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate, text2=string2)
endif

!----------------------------------------------------------------------
! bombproofing ... make sure the netcdf file is open.

call nc_check(nf90_inquire(ncid(ens_index)), &
              'obs_def_brightnessT.get_scalar_from_2Dhistory', 'inquire '//trim(strshort))

! bombproofing ... make sure the variable is the shape and size we expect

call nc_check(nf90_inq_varid(ncid(ens_index), trim(varstring), varid), &
        'obs_def_brightnessT.get_scalar_from_2Dhistory', 'inq_varid '//trim(strshort))
call nc_check(nf90_inquire_variable(ncid(ens_index), varid, xtype=xtype, ndims=ndims, &
        dimids=dimids, natts=natts), &
        'obs_def_brightnessT.get_scalar_from_2Dhistory','inquire variable '//trim(strshort))

if (ndims /= 2) then
   write(string1,*)trim(varstring),' is supposed to have 2 dimensions, it has',ndims
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

! If the variable is not a NF90_FLOAT, then the assumptions for processing
! the missing_value, _FillValue, etc., may not be correct.
if (xtype /= NF90_FLOAT) then
   write(string1,*)trim(varstring),' is supposed to be a 32 bit real. xtype = ', &
                   NF90_FLOAT,' it is ',xtype
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 1 is spatial
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(1), len=dimlen), &
        'obs_def_brightnessT.get_scalar_from_2Dhistory', 'inquire_dimension 1 '//trim(strshort))
if (dimlen /= nlon) then
   write(string1,*)'LON has length',nlon,trim(varstring),' has ',dimlen,'longitudes.'
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 2 is time
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(2), len=dimlen), &
        'obs_def_brightnessT.get_scalar_from_2Dhistory', 'inquire_dimension 2 '//trim(strshort))
if (dimlen /= ntime) then
   write(string1,*)'TIME has length',ntime,trim(varstring),' has ',dimlen,'times.'
   call error_handler(E_ERR, 'obs_def_brightnessT.get_scalar_from_2Dhistory', &
              string1, source, revision, revdate)
endif

!----------------------------------------------------------------------
! Find the timestep of interest
! FIXME ... since the history file contents are for the previous 30m,
! perhaps the closest time is not the best approximation.

call get_time(obs_time, second, day)
otime    = real(day,digits12) + real(second,digits12)/86400.0_digits12
timeinds = minloc(abs(rtime - otime))   ! these return 'arrays' ...
timei    = timeinds(1)

if (debug .and. do_output()) then
   write(*,*)'obs_def_brightnessT.get_scalar_from_2Dhistory:  targetT,   T,   T index is ', &
                                           otime,rtime(timei),timei
endif

if ( abs(otime - rtime(timei)) > 30*60 ) then
   if (debug .and. do_output()) then
      write(*,*)'obs_def_brightnessT.get_scalar_from_2Dhistory: no close time ... skipping observation'
      call print_time(obs_time,'obs_def_brightnessT.get_scalar_from_2Dhistory:observation time')
      call print_date(obs_time,'obs_def_brightnessT.get_scalar_from_2Dhistory:observation date')
   endif
   istatus = 2
   return
endif

!----------------------------------------------------------------------
! Grab exactly the scalar we want.

ncstart = (/ gridij, timei /)
nccount = (/      1,     1 /)

call nc_check(nf90_get_var(ncid(ens_index),varid,hyperslab,start=ncstart,count=nccount), &
     'obs_def_brightnessT.get_scalar_from_2Dhistory', 'get_var')

obs_val = hyperslab(1)

!----------------------------------------------------------------------
! Apply any netCDF attributes ...

io1 = nf90_get_att(ncid(ens_index), varid, '_FillValue' , spvalR4)
if ((io1 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io2 = nf90_get_att(ncid(ens_index), varid, 'missing_value' , spvalR4)
if ((io2 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io1 = nf90_get_att(ncid(ens_index), varid, 'scale_factor', scale_factor)
io2 = nf90_get_att(ncid(ens_index), varid, 'add_offset'  , add_offset)

if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor + add_offset
elseif (io1 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor
elseif (io2 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val + add_offset
endif

if (obs_val /= MISSING_R8) istatus = 0

end subroutine get_scalar_from_2Dhistory


!======================================================================


!======================================================================


function interactive_real(str1,minvalue,maxvalue)
real(r8)                       :: interactive_real
character(len=*),   intent(in) :: str1
real(r8), optional, intent(in) :: minvalue
real(r8), optional, intent(in) :: maxvalue

integer :: i

interactive_real = MISSING_R8

! Prompt with a minimum amount of error checking

if     (present(minvalue) .and. present(maxvalue)) then

   interactive_real = minvalue - 1.0_r8
   MINMAXLOOP : do i = 1,10
      if ((interactive_real >= minvalue) .and. (interactive_real <= maxvalue)) exit MINMAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive_real = minvalue - 1.0_r8
   MINLOOP : do i=1,10
      if (interactive_real >= minvalue) exit MINLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive_real = maxvalue + 1.0_r8
   MAXLOOP : do i=1,10
      if (interactive_real <= maxvalue) exit MAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive_real
endif

end function interactive_real


!======================================================================


function interactive_char(str1)
character                      :: interactive_char
character(len=*),   intent(in) :: str1

interactive_char = ''

write(*, *) 'Enter '//str1
read( *, *) interactive_char

end function interactive_char


!======================================================================


function interactive_int(str1,minvalue,maxvalue)
integer                       :: interactive_int
character(len=*),  intent(in) :: str1
integer, optional, intent(in) :: minvalue
integer, optional, intent(in) :: maxvalue

integer :: i

interactive_int = MISSING_I

! Prompt with a minimum amount of error checking

if     (present(minvalue) .and. present(maxvalue)) then

   interactive_int = minvalue - 1
   MINMAXLOOP : do i = 1,10
      if ((interactive_int >= minvalue) .and. (interactive_int <= maxvalue)) exit MINMAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive_int = minvalue - 1
   MINLOOP : do i=1,10
      if (interactive_int >= minvalue) exit MINLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive_int = maxvalue + 1
   MAXLOOP : do i=1,10
      if (interactive_int <= maxvalue) exit MAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive_int
endif

end function interactive_int


!======================================================================


subroutine key_within_range(key, routine)

! Make sure we are addressing within the metadata arrays

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

! fine -- no problem.
if ((key > 0) .and. (key <= amsrekey)) return

! Bad news. Tell the user.
write(string1, *) 'key (',key,') not within known range ( 1,', amsrekey,')'
call error_handler(E_ERR,routine,string1,source,revision,revdate)

end subroutine key_within_range


!======================================================================


subroutine grow_metadata(key, routine)
!----------------------------------------------------------------------
! If the allocatable metadata arrays are not big enough ... try again

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

integer :: orglength
type(amsre_metadata), allocatable, dimension(:) :: safe_metadata

! fine -- no problem.
if ((key > 0) .and. (key <= MAXamsrekey)) return

orglength   =     MAXamsrekey
MAXamsrekey = 2 * orglength

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*MAXamsrekey) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

! News. Tell the user we are increasing storage.
write(string1, *) 'key (',key,') exceeds Nmax_amsre_Tb (',orglength,')'
write(string2, *) 'Increasing Nmax_amsre_Tb to ',MAXamsrekey
call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

allocate(safe_metadata(orglength))
safe_metadata(:) = observation_metadata(:)

deallocate(observation_metadata)
  allocate(observation_metadata(MAXamsrekey))

observation_metadata(1:orglength)             = safe_metadata(:)
observation_metadata(orglength+1:MAXamsrekey) = missing_metadata

deallocate(safe_metadata)

end subroutine grow_metadata


!======================================================================


subroutine check_iostat(istat, routine, varname, msgstring)

integer,          intent(in) :: istat
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: msgstring

if ( istat /= 0 ) then
   write(string1,*)'istat should be 0 but is ',istat,' for '//varname
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
end if

end subroutine check_iostat


!======================================================================


subroutine test_block
! Defining the region if running in a single column is tricky.
! We have lat, lon, and the area of the gridcell which we assume to be basically square.
! The square root of the area defines the length of the edge of the gridcell which
! can then be interpreted as the diameter of a circle. Any observation within this
! distance is close enough.

real(r8) :: x1, x2, y1, y2, d1, d2, d3, d4, radius, distance
type(location_type) :: gridloc, testloc

gridloc  = set_location(0.0_r8, 0.0_r8, 0.0_r8, VERTISUNDEF)
testloc  = set_location(1.0_r8, 0.0_r8, 0.0_r8, VERTISUNDEF)
distance = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
write(*,*)'TJH DEBUG: 1 degree at the equator has distance ',distance,' km.'

x1 = 286.8750               ! gridcell longitude is 287.50
x2 = 288.1250
y1 = 42.40837669372559      ! gridcell latitude is 42.8795814514160
y2 = 43.35078620910645

! compute the distance along the top of the grid cell
gridloc = set_location(x1, y1, 0.0_r8, VERTISUNDEF)
testloc = set_location(x2, y1, 0.0_r8, VERTISUNDEF)
d1      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
gridloc = set_location(x2, y2, 0.0_r8, VERTISUNDEF)
d2      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
testloc = set_location(x1, y2, 0.0_r8, VERTISUNDEF)
d3      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM
gridloc = set_location(x1, y1, 0.0_r8, VERTISUNDEF)
d4      = get_dist(gridloc, testloc, no_vert = .TRUE.) * RAD2KM

write(*,*)
write(*,*)'lengths         are ',d1,d2,d3,d4
write(*,*)'rectangular area is ',((d1+d3)/2.0_r8)*((d2+d4)/2.0_r8)

d1 = sqrt(area(1))
d2 = sqrt(2.0_r8 * d1**2)
d3 = d2 / 2.0_r8
radius = sqrt(2.0_r8 * area(1))/2.0_r8
write(*,*)'length of one side of a square is x = sqrt(area) = ',d1
write(*,*)'diagonal is                sqrt(x**2 + x**2)     = ',d2
write(*,*)'radius   is                sqrt(x**2 + x**2)/2.0 = ',d3
write(*,*)'-or-                       sqrt(area + area)/2.0 = ',radius

write(*,*)'TJH DEBUG: Radius,area is ',distance,PI*distance**2, &
          ' gridcell area is ',area(1)

end subroutine test_block


!======================================================================


end module obs_def_brightnessT_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

This come from the history file

float WATSAT(levgrnd, lat, lon) ;
      WATSAT:long_name = "saturated soil water content (porosity)" ;
      WATSAT:units = "mm3/mm3" ;
      WATSAT:_FillValue = 1.e+36f ;
      WATSAT:missing_value = 1.e+36f ;

All these come from the restart file

double snw_rds(column, levsno) ;
       snw_rds:long_name = "snow layer effective radius" ;
       snw_rds:units = "um" ;

double H2OSOI_LIQ(column, levtot) ;
       H2OSOI_LIQ:long_name = "liquid water" ;
       H2OSOI_LIQ:units = "kg/m2" ;

double H2OSOI_ICE(column, levtot) ;
       H2OSOI_ICE:long_name = "ice lens" ;
       H2OSOI_ICE:units = "kg/m2" ;

SWE is H2OSOI_LIQ + H2OSOI_ICE, basically.
density is SWE / DZSNO

double T_GRND(column) ;
       T_GRND:long_name = "ground temperature" ;
       T_GRND:units = "K" ;

double T_SOISNO(column, levtot) ;
       T_SOISNO:long_name = "soil-snow temperature" ;
       T_SOISNO:units = "K" ;

double DZSNO(column, levsno) ;
       DZSNO:long_name = "snow layer thickness" ;
       DZSNO:units = "m" ;

double ZSNO(column, levsno) ;
       ZSNO:long_name = "snow layer depth" ;
       ZSNO:units = "m" ;

double ZISNO(column, levsno) ;
       ZISNO:long_name = "snow interface depth" ;
       ZISNO:units = "m" ;


