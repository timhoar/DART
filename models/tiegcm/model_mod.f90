! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!-------------------------------------------------------------------------------
!
! Interface for HAO-TIEGCM
!
!-------------------------------------------------------------------------------

use        types_mod, only : r8, digits12, missing_r8, i4, PI,                      &
                             earth_radius, gravity, obstypelength

use time_manager_mod, only : time_type, set_calendar_type, set_time_missing,        &
                             set_time, get_time, print_time,                        &
                             set_date, get_date, print_date,                        &
                             operator(*),  operator(+), operator(-),                &
                             operator(>),  operator(<), operator(/),                &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_close_maxdist_init,                 &
                             get_close_obs_init, loc_get_close_obs => get_close_obs,&
                             set_location, get_location, query_location,            &
                             get_dist, vert_is_height, horiz_dist_only,             &
                             get_close_type, vert_is_undef, VERTISUNDEF,            &
                             VERTISPRESSURE, VERTISHEIGHT, vert_is_pressure

use    utilities_mod, only : file_exist, open_file, close_file, logfileunit,        &
                             error_handler, E_ERR, E_MSG, E_WARN, nmlfileunit,      &
                             do_output, find_namelist_in_file, check_namelist_read, &
                             do_nml_file, do_nml_term, nc_check, register_module,   &
                             file_to_text, find_textfile_dims

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,           &
                             KIND_V_WIND_COMPONENT,           &
                             KIND_TEMPERATURE,                &! neutral temperature obs
                             KIND_PRESSURE,                   &! neutral pressure obs
                             KIND_ELECTRON_DENSITY,           &! Ne obs
                             KIND_ATOMIC_OXYGEN_MIXING_RATIO, &! neutral composition obs
                             KIND_MOLEC_OXYGEN_MIXING_RATIO,  &! neutral composition obs
                             KIND_1D_PARAMETER,               &
                             KIND_GEOPOTENTIAL_HEIGHT,        & 
                             KIND_DENSITY_ION_OP,             &! Atomic oxygen ion obs
                             KIND_VERTICAL_TEC,               &! total electron content
                             get_raw_obs_kind_index

use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use mpi_utilities_mod,only : my_task_id

use typesizes
use netcdf

implicit none
private

!DART mandatory public interfaces
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

!TIEGCM specific routines
public :: get_restart_file_name, &
            read_TIEGCM_restart, &
          update_TIEGCM_restart

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   '$URL$'
character(len=32 ), parameter :: revision = '$Revision$'
character(len=128), parameter :: revdate  = '$Date$'

!-------------------------------------------------------------------------------
! namelist with default values
! output_state_vector = .true.  results in a "state-vector"   netCDF file
! output_state_vector = .false. results in a "prognostic-var" netCDF file

! IMPORTANT: Change output file names in tiegcm.nml to match these names
! i.e.  OUTPUT='tiegcm_restart_p.nc'
!       SECOUT='tiegcm_s.nc'
character(len=128) :: tiegcm_restart_file_name   = 'tiegcm_restart_p.nc'
character(len=128) :: tiegcm_secondary_file_name = 'tiegcm_s.nc'
character(len=128) :: tiegcm_namelist_file_name  = 'tiegcm.nml'
logical            :: output_state_vector = .false.
integer            :: debug = 0
logical            :: estimate_f10_7 = .false.

integer, parameter :: StateVarNmax = 80
integer, parameter :: Ncolumns = 2
character(len=NF90_MAX_NAME) :: state_variables(StateVarNmax * Ncolumns) = ' '
character(len=NF90_MAX_NAME) :: auxiliary_variables(StateVarNmax)

namelist /model_nml/ output_state_vector, tiegcm_restart_file_name, &
                     tiegcm_secondary_file_name, tiegcm_namelist_file_name, &
                     state_variables, auxiliary_variables, debug, &
                     estimate_f10_7

!-------------------------------------------------------------------------------
! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: storder
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens ! ntime, [nlev,] nlat, nlon
   integer  :: posdef
   integer  :: rank
   integer  :: varsize     ! prod(dimlens(1:rank))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   character(len=obstypelength) :: kind_string
   character(len=obstypelength) :: verticalvar
   integer  :: dart_kind
   integer  :: xtype
   integer  :: rangeRestricted
   real(r8) :: maxvalue
   real(r8) :: minvalue
   real(r8) :: missingR8
   logical  :: has_missing_value
   logical  :: prognostic  ! is this a prognostic variable (updateable)
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
end type progvartype

type(progvartype), dimension(StateVarNmax) :: progvar
integer :: nfields  ! number of prognostic variables

!-------------------------------------------------------------------------------
! define model parameters

integer                               :: nilev, nlev, nlon, nlat
real(r8),dimension(:),    allocatable :: lons, lats, levs, ilevs, plevs, pilevs
real(r8)                              :: TIEGCM_missing_value !! global attribute
real(r8)                              :: TIEGCM_reference_pressure
integer                               :: time_step_seconds
integer                               :: time_step_days
type(time_type)                       :: time_step

character(len=obstypelength) :: variable_table(StateVarNmax, Ncolumns)

! include_vTEC = .true.  vTEC must be calculated from other vars
! include_vTEC = .false. just ignore vTEC altogether

logical  :: include_vTEC = .true.
logical  :: include_vTEC_in_state = .false.

! IMPORTANT: 1 D model parameters (e.g., F107) are read in from "tiegcm.nml"
! When "estimate_f10_7 = .true.", "state_num_1d" should be greater than or
! equal to 1 so that 1 D model parameters will be included in the state vector
! (note "estimate_f10_7" option is still under
! development by Tomoko Matsuo as of June 24, 2011)
real(r8) :: f10_7

integer                               :: model_size
real(r8), allocatable                 :: ens_mean(:)

! FOR NOW OBS LOCATIONS ARE EXPECTED GIVEN IN HEIGHT [m],
! AND SO VERTICAL LOCALIZATION COORDINATE IS *always* HEIGHT
! (note that gravity adjusted geopotential height (ZG)
!  read in from "tiegcm_s.nc")
!integer                              :: vert_localization_coord = VERTISHEIGHT

real(r8), allocatable, dimension(:,:,:) :: ZG
integer :: ivarZG

logical                               :: first_pert_call = .true.
type(random_seq_type)                 :: random_seq


character(len = 129) :: string1, string2, string3
logical, save :: module_initialized = .false.

interface prog_var_to_vector
   module procedure var1d_to_vector
   module procedure var2d_to_vector
   module procedure var3d_to_vector
   module procedure var4d_to_vector
end interface

interface vector_to_prog_var
   module procedure vector_to_var1d
   module procedure vector_to_var2d
   module procedure vector_to_var3d
   module procedure vector_to_var4d
end interface

!===============================================================================
contains
!===============================================================================



subroutine static_init_model()
!-------------------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

integer :: iunit, io

character(len=128) :: restart_file
character(len=128) :: namelist_file
character(len=128) :: secondary_file

if (module_initialized) return ! only need to do this once

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the namelist entry for model_mod from input.nml
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Read in TIEGCM namelist input file (just for definition)
! Read in TIEGCM grid definition etc from TIEGCM restart file
! Read in TIEGCM auxiliary variables from TIEGCM 'secondary' file
restart_file   = adjustl(tiegcm_restart_file_name)
namelist_file  = adjustl(tiegcm_namelist_file_name)
secondary_file = adjustl(tiegcm_secondary_file_name)

call read_TIEGCM_namelist(trim(namelist_file))
call read_TIEGCM_definition(trim(restart_file))
call read_TIEGCM_secondary(trim(secondary_file))

! error-check and convert namelist input to variable_table
! fill variable table and the module 'progvar' database
call verify_state_variables( state_variables, restart_file, nfields, variable_table )

! Compute overall model size

model_size = progvar(nfields)%indexN

if (do_output()) then
   write(*,*) 'nlon  = ', nlon
   write(*,*) 'nlat  = ', nlat
   write(*,*) 'nlev  = ', nlev
   write(*,*) 'nilev = ', nilev
   write(*,*) 'model_size = ', model_size
   if (estimate_f10_7) &
   write(*,*) 'estimating f10.7 ... init value ',f10_7
endif

allocate (ens_mean(model_size))

! Might as well use the Gregorian Calendar
call set_calendar_type('Gregorian')

! The time_step in terms of a time type must also be initialized.
time_step = set_time(time_step_seconds, time_step_days)

! call test_module()

end subroutine static_init_model



subroutine init_conditions(x)
!-------------------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'no good way to specify initial conditions'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

x(:) = MISSING_R8 ! just to silence compiler messages

end subroutine init_conditions



subroutine adv_1step(x, time)
!-------------------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

write(string1,*) 'TIEGCM cannot be advanced as a subroutine from within DART.'
write(string2,*) 'check your input.nml setting of "async" and try again.'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

! just to silence compiler messages

x(:) = MISSING_R8
call print_time(time)

end subroutine adv_1step



function get_model_size()
!-------------------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!-------------------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'no good way to specify initial time'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time



subroutine model_interpolate(x, location, ikind, obs_val, istatus)
!-------------------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The ikind variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: ikind
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

integer  :: ivar
integer  :: i, vstatus, which_vert
integer  :: lat_below, lat_above, lon_below, lon_above
integer  :: zero_lon_index
real(r8) :: lon_fract, lat_fract
real(r8) :: lon, lat, height, lon_lat_lev(3)
real(r8) :: bot_lon, top_lon, delta_lon, bot_lat, top_lat, delta_lat
real(r8) :: val(2,2), a(2)

if ( .not. module_initialized ) call static_init_model

! Default for successful return
istatus = 0
vstatus = 0

! Get the position
! FOR NOW OBS VERTICAL LOCATION IS ALWAYS HEIGHT
lon_lat_lev = get_location(location)
lon = lon_lat_lev(1) ! degree
lat = lon_lat_lev(2) ! degree
if(vert_is_height(location)) then
   height = lon_lat_lev(3)
elseif ((vert_is_undef(location)) .or. (ikind == KIND_VERTICAL_TEC)) then
   ! vertical location is undefined - e.g., MIDAS_TEC
   height = missing_r8
else
   which_vert = nint(query_location(location))
   write(string1,*) 'vertical coordinate type:',which_vert,' cannot be handled'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
endif

! Get lon and lat grid specs
bot_lon   = lons(1)                         ! 180.
delta_lon = abs((lons(1)-lons(2)))          ! 5. or 2.5
zero_lon_index = int(bot_lon/delta_lon) + 1 ! 37 or 73
top_lon   = lons(nlon)                      ! 175. or 177.5
bot_lat   = lats(1)                         !
top_lat   = lats(nlat)                      !
delta_lat = abs((lats(1)-lats(2)))          !


! Compute bracketing lon indices:
! TIEGCM [-180 175]  DART [180, 185, ..., 355, 0, 5, ..., 175]
if(lon > top_lon .and. lon < bot_lon) then     ! at wraparound point [175 < lon < 180]
   lon_below = nlon
   lon_above = 1
   lon_fract = (lon - top_lon) / delta_lon
elseif (lon >= bot_lon) then                  ! [180 <= lon <= 360]
   lon_below = int((lon - bot_lon) / delta_lon) + 1
   lon_above = lon_below + 1
   lon_fract = (lon - lons(lon_below)) / delta_lon
else                                           ! [0 <= lon <= 175 ]
   lon_below = int((lon - 0.0_r8) / delta_lon) + zero_lon_index
   lon_above = lon_below + 1
   lon_fract = (lon - lons(lon_below)) / delta_lon
endif

! Compute neighboring lat rows: TIEGCM [-87.5, 87.5] DART [-90, 90]
! NEED TO BE VERY CAREFUL ABOUT POLES; WHAT'S BEING DONE IS NOT GREAT!
if(lat >= bot_lat .and. lat <= top_lat) then ! -87.5 <= lat <= 87.5
   lat_below = int((lat - bot_lat) / delta_lat) + 1
   lat_above = lat_below + 1
   lat_fract = (lat - lats(lat_below) ) / delta_lat
else if(lat < bot_lat) then ! South of bottom lat
   lat_below = 1
   lat_above = 1
   lat_fract = 1.0_r8
else                        ! North of top lat
   lat_below = nlat
   lat_above = nlat
   lat_fract = 1.0_r8
endif

! TJH FIXME ... this model_interpolate needs to support KIND_GEOPOTENTIAL_HEIGHT
! in order to be able to use obs_def_upper_atm_mod:get_expected_gnd_gps_vtec()

if ( ikind == KIND_GEOPOTENTIAL_HEIGHT ) then
   write(string1,*)'KIND_GEOPOTENTIAL_HEIGHT currently unsupported'
   call error_handler(E_ERR,'model_interpolate',string1,source, revision, revdate)
endif

call error_handler(E_ERR,'model_interpolate','TJH FIXME unfinished', source, revision, revdate)

! Now, need to find the values for the four corners
if ((vert_is_undef(location)) .or. (ikind == KIND_VERTICAL_TEC)) then !2D fields

  if (ikind == KIND_VERTICAL_TEC) then

    ivar = FindVar_by_type(ikind)
    val(1, 1) = x(get_index(ivar, lon_below, lat_below))
    val(1, 2) = x(get_index(ivar, lon_below, lat_above))
    val(2, 1) = x(get_index(ivar, lon_above, lat_below))
    val(2, 2) = x(get_index(ivar, lon_above, lat_above))
  endif

else !3D fields

      call get_val(val(1, 1), x, lon_below, lat_below, height, ikind, vstatus)
   if (vstatus /= 1) &
      call get_val(val(1, 2), x, lon_below, lat_above, height, ikind, vstatus)
   if (vstatus /= 1) &
      call get_val(val(2, 1), x, lon_above, lat_below, height, ikind, vstatus)
   if (vstatus /= 1) &
      call get_val(val(2, 2), x, lon_above, lat_above, height, ikind, vstatus)

endif

! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no

istatus = vstatus
if(istatus /= 1) then
   do i = 1, 2
      a(i) = lon_fract * val(2, i) + (1.0_r8 - lon_fract) * val(1, i)
   end do
   obs_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
else
   obs_val = missing_r8
endif

end subroutine model_interpolate



function get_model_time_step()
!-------------------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_kind)
!-------------------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_kind

integer  :: remainder
integer  :: myindx, lon_index, lat_index, lev_index
integer  :: ivar
real(r8) :: height

if ( .not. module_initialized ) call static_init_model

ivar   = Find_Variable_by_index(index_in,'get_state_meta_data')
myindx = index_in - progvar(ivar)%index1 + 1

if     (progvar(ivar)%rank == 0) then  ! scalars ... no location
   location = set_location(0.0_r8, 0.0_r8,  400000.0_r8, VERTISHEIGHT)

   write(*,*)'TJH ivar, index_in, myindx ', &
                  ivar, index_in, myindx

elseif (progvar(ivar)%rank == 1) then  ! time?
   write(string1,*)trim(progvar(ivar)%varname),'has unsupported shape (1D)'
   write(string2,*)'dimension ('//trim(progvar(ivar)%dimnames(1))//') ... unknown location' 
   call error_handler(E_ERR,'get_state_meta_data', string1, &
              source, revision, revdate, text2=string2)

elseif (progvar(ivar)%rank == 2) then  ! something, and time?
   write(string1,*)trim(progvar(ivar)%varname), &
                   'has unsupported shape (2D) ... unknown location'
   write(string2,*)'dimension 1 = ('//trim(progvar(ivar)%dimnames(1))//')'
   write(string3,*)'dimension 2 = ('//trim(progvar(ivar)%dimnames(2))//')'
   call error_handler(E_ERR,'get_state_meta_data', string1, &
              source, revision, revdate, text2=string2, text3=string3)

elseif (progvar(ivar)%rank == 3) then  ! [Fortran ordering] = lon, lat, time
   ! The time dimension is always length 1, so it really doesn't matter.

   lat_index = 1 + (myindx - 1) / nlon ! relies on integer arithmetic
   lon_index = myindx - (lat_index-1) * nlon

   write(*,*)'TJH ivar, index_in, myindx, lon_index, lat_index', &
                  ivar, index_in, myindx, lon_index, lat_index

   if (trim(progvar(ivar)%varname) == 'VTEC') then
      ! assign arbitrary height to allow for localization
      location = set_location(lons(lon_index), lats(lat_index), 300000.0_r8, VERTISHEIGHT)
   else
      location = set_location(lons(lon_index), lats(lat_index), 0.0_r8, VERTISUNDEF)
   endif

elseif (progvar(ivar)%rank == 4) then  ! [Fortran ordering] = lon, lat, lev, time

   lev_index = 1 + (myindx - 1) / (nlon * nlat)
   remainder = myindx - (lev_index-1) * nlon * nlat
   lat_index = 1 + (remainder - 1) / nlon
   lon_index = remainder - (lat_index-1) * nlon

!   write(*,*)'TJH ivar, index_in, myindx, loni, lati, levi', &
!                  ivar, index_in, myindx, lon_index, lat_index, lev_index

   height    = convert_lnpressure_to_height(ivar, lon_index, lat_index, lev_index) 
   location  = set_location(lons(lon_index), lats(lat_index), height, VERTISHEIGHT)

else
   write(string1,*)'Problem with DART variable ',trim(progvar(ivar)%varname)
   write(string2,*)'has unsupported number (',progvar(ivar)%rank,') of dimensions.'
   call error_handler(E_ERR,'get_state_meta_data', string1, &
                      source, revision, revdate, text2=string2)
endif

! If the type is wanted, return it
if(present(var_kind)) var_kind = progvar(ivar)%dart_kind

end subroutine get_state_meta_data




subroutine end_model()
!-------------------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if ( .not. module_initialized ) call static_init_model


end subroutine end_model



function nc_write_model_atts( ncid ) result (ierr)
!-------------------------------------------------------------------------------
! TJH 20 Dec 2013 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
!
! TJH 20 Dec 2013 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset


integer, intent(in)  :: ncid      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

integer :: myndims
integer :: ivar, VarID
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
character(len=NF90_MAX_NAME) :: varname

integer :: lonDimID, latDimID, levDimID, ilevDimID
integer :: lonVarID, latVarID, levVarID, ilevVarID

!-------------------------------------------------------------------------------
! variables for the namelist output
!-------------------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_tiegcm_namelist

!-------------------------------------------------------------------------------
! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.
!-------------------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncid refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_Inquire(ncid, nDimensions, nVariables, &
              nAttributes, unlimitedDimID), 'nc_write_model_atts','inquire')
call nc_check(nf90_Redef(ncid),'nc_write_model_atts','redef')

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncid, name='NMLlinelen', dimid = linelenDimID), &
       'nc_write_model_atts', 'inq_dimid NMLlinelen')
call nc_check(nf90_inq_dimid(ncid=ncid, name='copy', dimid=MemberDimID),&
       'nc_write_model_atts', 'inq_dimid copy')
call nc_check(nf90_inq_dimid(ncid=ncid, name='time', dimid=  TimeDimID),&
       'nc_write_model_atts', 'inq_dimid time')

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
                     ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncid, name='StateVariable', &
                        len=model_size, dimid = StateVarDimID),&
       'nc_write_model_atts', 'state def_dim')

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,str1    ),&
       'nc_write_model_atts', 'creation put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model_source'  ,source  ),&
       'nc_write_model_atts', 'source put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model_revision',revision),&
       'nc_write_model_atts', 'revision put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model_revdate' ,revdate ),&
       'nc_write_model_atts', 'revdate put')
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model','TIEGCM'         ),&
       'nc_write_model_atts', 'model put')

!-------------------------------------------------------------------------------
! Determine shape of namelist.
! long lines are truncated when read into textblock
!-------------------------------------------------------------------------------

call find_textfile_dims(tiegcm_namelist_file_name, nlines, linelen)
if (nlines > 0) then
  has_tiegcm_namelist = .true.
else
  has_tiegcm_namelist = .false.
endif

if (has_tiegcm_namelist) then
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncid, name='tiegcmNMLnlines', &
          len = nlines, dimid = nlinesDimID), &
          'nc_write_model_atts', 'def_dim tiegcmNMLnlines')

   call nc_check(nf90_def_var(ncid,name='tiegcm_nml', xtype=nf90_char, &
          dimids = (/ linelenDimID, nlinesDimID /), varid=nmlVarID), &
          'nc_write_model_atts', 'def_var tiegcm_namelist')

   call nc_check(nf90_put_att(ncid, nmlVarID, 'long_name', &
          'contents of '//trim(tiegcm_namelist_file_name)), &
          'nc_write_model_atts', 'put_att tiegcm_namelist')

endif

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncid,name='StateVariable', xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID), &
              'nc_write_model_atts', 'statevariable def_var')
   call nc_check(nf90_put_att(ncid, StateVarVarID, 'long_name', 'State Variable ID'), &
              'nc_write_model_atts', 'statevariable long_name')
   call nc_check(nf90_put_att(ncid, StateVarVarID, 'units',     'indexical'), &
              'nc_write_model_atts', 'statevariable units')
   call nc_check(nf90_put_att(ncid, StateVarVarID, 'valid_range', (/ 1, model_size /)), &
              'nc_write_model_atts', 'statevariable valid_range')

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncid, name='state', xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID), 'nc_write_model_atts', 'state def_var')
   call nc_check(nf90_put_att(ncid, StateVarID, 'long_name', 'model state or fcopy'), &
              'nc_write_model_atts', 'state long_name')

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncid), 'nc_write_model_atts', 'state enddef')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncid, StateVarVarID, (/ (i,i=1,model_size) /) ), &
              'nc_write_model_atts', 'state put_var')

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid= ncid ,  name='lon', &
             & len =  nlon,  dimid= lonDimID), 'nc_write_model_atts lon')
   call nc_check(nf90_def_dim(ncid= ncid ,  name='lat', &
             & len =  nlat,  dimid= latDimID), 'nc_write_model_atts lat')
   call nc_check(nf90_def_dim(ncid= ncid ,  name='lev', &
             & len =  nlev,  dimid= levDimID), 'nc_write_model_atts lev')
   call nc_check(nf90_def_dim(ncid= ncid ,  name='ilev', &
             & len = nilev,  dimid=ilevDimID), 'nc_write_model_atts ilev')

   !----------------------------------------------------------------------------
   ! Create the (empty) Variables and the Attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid, name='lon', &
             & xtype=nf90_double, dimids=lonDimID, varid=lonVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, lonVarID, &
             & 'long_name', 'geographic longitude (-west, +east)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, lonVarID, 'units', 'degrees_east'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, lonVarID, 'valid_range', &
             & (/ -180.0_r8, 180.0_r8 /)),'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid, name='lat', &
             & xtype=nf90_double, dimids=latDimID, varid=latVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, latVarID, &
             & 'long_name', 'geographic latitude (-south +north)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, latVarID, 'units', 'degrees_north'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, latVarID, 'valid_range', &
             & (/ -90.0_r8, 90.0_r8 /)),'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid, name='lev', &
             & xtype=nf90_double, dimids=levDimID, varid=levVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'long_name', 'midpoint levels'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'short_name', 'ln(p0/p)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'units', ''),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, levVarID, 'positive', 'up'),&
             'nc_write_model_atts')

   call nc_check(nf90_def_var(ncid, name='ilev', &
             & xtype=nf90_double, dimids=ilevDimID, varid=ilevVarID),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'long_name', 'interface levels'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'short_name', 'ln(p0/p)'),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'units', ''),&
             'nc_write_model_atts')
   call nc_check(nf90_put_att(ncid, ilevVarID, 'positive', 'up'),&
             'nc_write_model_atts')

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and their Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncid, MemberDimID, unlimitedDimID, myndims, mydimids)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncid, name=trim(varname), xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(varname)//' def_var' )

      call nc_check(nf90_put_att(ncid, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(varname)//' put_att long_name' )
      call nc_check(nf90_put_att(ncid, VarID, 'units',     trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(varname)//' put_att units' )
      call nc_check(nf90_put_att(ncid, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(varname)//' put_att dart_kind' )
   enddo

   call nc_check(nf90_enddef(ncid), 'nc_write_model_atts', 'prognostic enddef')

   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncid, lonVarID, lons),'nc_write_model_atts','put_var lons')
   call nc_check(nf90_put_var(ncid, latVarID, lats),'nc_write_model_atts','put_var lats')
   call nc_check(nf90_put_var(ncid, levVarID, levs),'nc_write_model_atts','put_var levs')
   call nc_check(nf90_put_var(ncid,ilevVarID,ilevs),'nc_write_model_atts','put_var ilevs')

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_tiegcm_namelist) then
   call file_to_text(tiegcm_namelist_file_name, textblock)
   call nc_check(nf90_put_var(ncid, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncid), 'nc_write_model_atts', 'sync')
if (do_output()) write (*,*) 'nc_write_model_atts: netCDF file ', ncid, ' is synched '

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncid, statevec, copyindex, timeindex ) result (ierr)
!-------------------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

integer,                intent(in) :: ncid      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames
integer :: i, ivar, VarID, ncNdims, dimlen, numdims, timedimcounter
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)     :: data_1d_array
real(r8), allocatable, dimension(:,:)   :: data_2d_array
real(r8), allocatable, dimension(:,:,:) :: data_3d_array

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncid refers to an open netCDF file,
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_Inquire(ncid, nDimensions, nVariables, &
                  nAttributes, unlimitedDimID), 'nc_write_model_vars', 'inquire')

call nc_check(nf90_inq_dimid(ncid, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy ')

call nc_check(nf90_inq_dimid(ncid, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time ')

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncid, 'state', StateVarID), &
          'nc_write_model_vars', 'state inq_varid' )
   call nc_check(NF90_put_var(ncid, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)), &
          'nc_write_model_vars', 'state put_var')

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields

      varname = trim(progvar(ivar)%varname)
      string2 = trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried.
      ! This requires that Time is the unlimited dimension (the last one in Fortran),
      ! and that 'copy' is the second-to-last. The variables declared in the DART
      ! diagnostic files are required to have the same shape as in the source
      ! restart file. If Time is present there, it must also be the 'last' one.

      ! FIXME ... somewhere I should ensure that IF time is present in the original
      ! prognostic variable from the model, it is the last/unlimited dimension.

      call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      timedimcounter = 0
      mystart(:) = 1
      mycount(:) = 1
      DimCheck : do i = 1,ncNdims

         write(string1,'(A,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimnames(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if (dimIDs(i) == CopyDimID) cycle DimCheck
         if (dimIDs(i) == TimeDimID) then
            timedimcounter = 1
            cycle DimCheck
         endif

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate)
         endif

         mycount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      ! If the original variable is shaped:
      !     XXXXXX(time, somedimension) -or- XXXXXX(somedimension)
      ! then it is a 1D variable in our context.
      ! If it is shaped
      !     XXXXXX(time, south_north, west_east) -or- XXXXXX(south_north, west_east)
      ! it really is 2D ...
      !
      ! this adjustment to numdims below is to remove the Time dimension

      numdims = progvar(ivar)%rank - timedimcounter

      if ( (do_output()) .and. debug > 99 ) then 
         write(*,*)'nc_write_model_vars '//trim(varname)//' numdims is ',numdims
         write(*,*)'nc_write_model_vars '//trim(varname)//' start   is ',mystart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count   is ',mycount(1:ncNdims)
         write(*,'(1x,A20)')'nc_write_model_vars ',dimnames(1:progvar(ivar)%rank)
      endif

      if ( numdims <= 0 ) then  ! handle true scalars

         if ( ncNdims /= 2 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 2 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         call nc_check(nf90_put_var(ncid, VarID, (/ statevec(progvar(ivar)%index1) /), &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))

      elseif ( numdims == 1 ) then ! 1D arrays

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ))
         call vector_to_prog_var(statevec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncid, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d_array)

      elseif ( numdims == 2 ) then ! 2D arrays, usually (lon,lat)

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2) ))
         call vector_to_prog_var(statevec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncid, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      elseif ( numdims == 3 ) then ! 3D arrays, usually (lon,lat,lev)

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_3d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2),  &
                                 progvar(ivar)%dimlens(3) ))
         call vector_to_prog_var(statevec, ivar, data_3d_array)
         call nc_check(nf90_put_var(ncid, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      else

         write(string1,*)'do not know how to handle tiegcm variables with more than 3 non-time dimensions'
         write(string2,*)trim(progvar(ivar)%varname),' has dimensions ', progvar(ivar)%dimnames(1:progvar(ivar)%rank)
         call error_handler(E_ERR,'nc_write_model_vars',string1, &
                    source,revision,revdate,text2=string2)

      endif

   enddo

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

if (do_output()) write (*,*) 'nc_write_model_vars: Finished filling variables '
call nc_check(nf90_sync(ncid), 'nc_write_model_vars', 'sync')
if (do_output()) write (*,*) 'nc_write_model_vars: netCDF file is synched '

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!-------------------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

integer                 :: i, variable_type
type(location_type)     :: temp_loc

if ( .not. module_initialized ) call static_init_model

! An interface is provided
interf_provided = .true.

! If first call initialize random sequence
! CAUTION: my_task_id is NOT emsemble member number
! For example, my_task_id will be in [0,N-1]
! if a single instance of the model using N MPI tasks.

if(first_pert_call) then
   call init_random_seq(random_seq,my_task_id())
   first_pert_call = .false.
endif

do i = 1, get_model_size()
   call get_state_meta_data(i, temp_loc, variable_type)
   if(variable_type == KIND_1D_PARAMETER) then
      pert_state(i) = random_gaussian(random_seq,state(i),20.0_r8)
   else
      pert_state(i) = state(i)
   endif
end do

end subroutine pert_model_state



subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                            num_close, close_ind, dist)
!-------------------------------------------------------------------------------
!
! Given a DART ob (referred to as "base") and a set of obs priors or
! state variables returns the subset of close ones to the "base" ob, their
! indices, and their distances to the "base" ob...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate.
! FOR NOW VERTICAL LOCALIZATION IS DONE ONLY IN HEIGHT (ZG)
! OBS VERTICAL LOCATION IS GIVEN IN HEIGHT (model_interpolate)
! STATE VERTICAL LOCATION IS GIVEN IN HEIGHT (get_state_meta_data)

! Note that both base_obs_loc and obs_loc are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling
! routine. The calling routine is always filter_assim and these arrays are local
! arrays within filter_assim. In other words, these modifications will only
! matter within filter_assim, but will not propagate backwards to filter.

type(get_close_type), intent(in)     :: gc
type(location_type),  intent(inout)  :: base_obs_loc, obs_loc(:)
integer,              intent(in)     :: base_obs_kind, obs_kind(:)
integer,              intent(out)    :: num_close, close_ind(:)
real(r8),             intent(out)    :: dist(:)

integer                              :: k, t_ind

call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                                           num_close, close_ind, dist)

! Localize
if (estimate_f10_7) then

   do k = 1, num_close

      t_ind  = close_ind(k)

      if ( (obs_kind(t_ind) == KIND_MOLEC_OXYGEN_MIXING_RATIO) &
         .or. (obs_kind(t_ind) == KIND_U_WIND_COMPONENT) &
         .or. (obs_kind(t_ind) == KIND_V_WIND_COMPONENT) &
         .or. (obs_kind(t_ind) == KIND_TEMPERATURE) ) then

         !set distance to a very large value so that it won't get updated
         dist(k) = 2.0_r8 * PI

      elseif (obs_kind(t_ind) == KIND_1D_PARAMETER) then
         if (estimate_f10_7) then
            dist(k) = dist(k)*0.25_r8
         else !not estimate_f10_7
            dist(k) = 2.0_r8 * PI
         endif

      endif

   enddo ! loop over k = 1, num_close

endif

end subroutine get_close_obs



subroutine ens_mean_for_model(filter_ens_mean)
!-------------------------------------------------------------------------------
!
! Not used in low-order models
! Stores provided ensemble mean within the module for later use

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model


!===============================================================================
! TIEGCM public routines 
!===============================================================================


subroutine read_TIEGCM_restart(file_name, statevec, model_time)
!-------------------------------------------------------------------------------
!
! Read TIEGCM restart file fields and pack it into a DART vector

character(len=*),       intent(in)  :: file_name
real(r8), dimension(:), intent(out) :: statevec
type(time_type),        intent(out) :: model_time

integer :: LonDimID, LatDimID, LevDimID, IlevDimID, TimeDimID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount

integer                                   :: ncid
integer                                   :: DimID, dimlen
integer                                   :: time_dimlen
integer                                   :: VarID, var_year_id
integer,  dimension(:), allocatable       :: yeartmp
integer,  dimension(:,:), allocatable     :: mtimetmp
integer,  parameter                       :: nmtime = 3
integer,  dimension(nmtime)               :: mtime  ! day, hour, minute
integer                                   :: year, utsec, doy
integer                                   :: nlevm1

real(r8), allocatable, dimension(:)       :: temp1D
real(r8), allocatable, dimension(:,:)     :: temp2D
real(r8), allocatable, dimension(:,:,:)   :: temp3D
real(r8), allocatable, dimension(:,:,:,:) :: temp4D

integer :: i, ivar, ncNdims

if ( .not. module_initialized ) call static_init_model

nlevm1 = nlev - 1

if( .not. file_exist(file_name)) then
   write(string1,*)trim(file_name)//' does not exist.'
   call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

if (do_output()) print *, 'read_TIEGCM_restart:reading restart:', file_name

! Open the netCDF file and then read all the static information.

call nc_check( nf90_open( file_name, NF90_NOWRITE, ncid ), &
                              'read_TIEGCM_restart', 'open')

!... check for matching dimensions
call nc_check( nf90_inq_dimid(ncid, 'lon', LonDimID), &
         'read_TIEGCM_restart', 'inq_dimid lon')
call nc_check( nf90_inquire_dimension(ncid, LonDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension lon')
if (dimlen .ne. nlon) then
  write(string1, *) trim(file_name), ' dim_lon = ',dimlen, ' DART expects ',nlon
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'lat', LatDimID), &
         'read_TIEGCM_restart', 'inq_dimid lat')
call nc_check( nf90_inquire_dimension(ncid, LatDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension lat')
if (dimlen .ne. nlat) then
  write(string1, *) trim(file_name), ' dim_lat = ',dimlen, ' DART expects ',nlat
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'lev', LevDimID), &
         'read_TIEGCM_restart', 'inq_dimid lev')
call nc_check( nf90_inquire_dimension(ncid, LevDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension lev')
if (dimlen .ne. nlev) then
  write(string1, *) trim(file_name), ' dim_lev = ',dimlen, ' DART expects ',nlev
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'ilev', IlevDimID), &
         'read_TIEGCM_restart', 'inq_dimid ilev')
call nc_check( nf90_inquire_dimension(ncid, IlevDimID, len=dimlen), &
         'read_TIEGCM_restart', 'inquire_dimension ilev')
if (dimlen .ne. nilev) then
  write(string1, *) trim(file_name), ' dim_ilev = ',dimlen, ' DART expects ',nilev
  call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'time', TimeDimID), &
         'read_TIEGCM_restart', 'inq_dimid time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'read_TIEGCM_restart', 'inquire_dimension time')

! Loop over all variables in DART state.
! Make sure the shape of the variable matches what we expect.
! Get the hyperslab with the most current (last) time.

ReadVariable: do ivar = 1,nfields

   string2 = trim(file_name)//' '//trim(progvar(ivar)%varname)

   if (trim(progvar(ivar)%varname) == 'f10_7') then
      statevec( progvar(ivar)%index1 ) = f10_7
      cycle ReadVariable
   endif

   if (trim(progvar(ivar)%varname) == 'ZG') then
      ! FIXME ... convert ZG to meters here and be done with it ...
      call prog_var_to_vector(ivar, ZG, statevec)
      deallocate(ZG)
      cycle ReadVariable
   endif

   if (trim(progvar(ivar)%varname) == 'VTEC') then
      allocate(temp2D(nlon,nlat))
      call create_vtec(ncid, time_dimlen, temp2D)
      call prog_var_to_vector(ivar, temp2D, statevec)
      deallocate(temp2D)
      cycle ReadVariable
   endif

   call nc_check(nf90_inq_varid(ncid, trim(progvar(ivar)%varname), VarID), &
          'read_TIEGCM_restart', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'read_TIEGCM_restart', 'inquire '//trim(string2))

   if (ncNdims /= progvar(ivar)%rank) then
      write(string1,*)trim(string2),'has ',ncNDims,'dimensions.'
      write(string3,*)'same thing in DART has ',progvar(ivar)%rank,' dimensions.'
      call error_handler(E_ERR, 'read_TIEGCM_restart', trim(string1), &
                         source, revision, revdate, text2=trim(string3))
   endif

   mystart(:) = 1
   mycount(:) = 1
   DimCheck : do i = 1,progvar(ivar)%rank

      write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
               'read_TIEGCM_restart', trim(string1))

      ! Check the shape of the variable
      if ( dimIDs(i) == TimeDimID ) then
         dimlen = 1
      elseif ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string3,*)' but it should be.'
         call error_handler(E_ERR, 'read_TIEGCM_restart', trim(string1), &
                         source, revision, revdate, text2=trim(string3))
      endif

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = time_dimlen
   where(dimIDs == TimeDimID) mycount = 1
   where(dimIDs ==  LonDimID) mycount = nlon
   where(dimIDs ==  LatDimID) mycount = nlat
   where(dimIDs ==  LevDimID) mycount = nlev
   where(dimIDs == iLevDimID) mycount = nilev

   if (do_output() .and. (debug > 99)) then
      write(*,*)'read_TIEGCM_restart:',trim(string2)
      write(*,*)'read_TIEGCM_restart: start ',mystart(1:ncNdims)
      write(*,*)'read_TIEGCM_restart: count ',mycount(1:ncNdims)
      write(*,*)
   endif

   if     (progvar(ivar)%rank == 1) then
      allocate(temp1D(mycount(1))) 
      call nc_check(nf90_get_var(ncid, VarID, values=temp1D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call prog_var_to_vector(ivar, temp1D, statevec)
      deallocate(temp1D) 

   elseif (progvar(ivar)%rank == 2) then
      allocate(temp2D(mycount(1), mycount(2)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp2D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call prog_var_to_vector(ivar, temp2D, statevec)
      deallocate(temp2D)

   elseif (progvar(ivar)%rank == 3) then
      allocate(temp3D(mycount(1), mycount(2), mycount(3)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp3D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call prog_var_to_vector(ivar, temp3D, statevec)
      deallocate(temp3D)

   elseif (progvar(ivar)%rank == 4) then
      allocate(temp4D(mycount(1), mycount(2), mycount(3), mycount(4)))
      call nc_check(nf90_get_var(ncid, VarID, values=temp4D, &
              start = mystart(1:ncNdims), count = mycount(1:ncNdims)), &
              'read_TIEGCM_restart', 'get_var '//trim(string2))
      call fill_top(progvar(ivar)%varname,temp4D,dimIDs,LevDimID,iLevDimID)
      call prog_var_to_vector(ivar, temp4D, statevec)
      deallocate(temp4D)
   else
      write(string1,*) trim(string2),' has unsupported number of dims (', &
                                     progvar(ivar)%rank,')'
      call error_handler(E_ERR, 'read_TIEGCM_restart', trim(string1), &
                         source, revision, revdate)
   endif

enddo ReadVariable

! Get the current time/date of the model and convert to a DART time_type
! Convert the year/doy/hour/minute to a dart time.
! Start by finding the dart time of the year and adding the rest to it.

call nc_check(nf90_inq_dimid(ncid, 'mtimedim', DimID), 'read_TIEGCM_restart', &
                  'inq_dimid mtimedim')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
        'read_TIEGCM_restart', 'inquire_dimension mtimedim')
if (dimlen /= nmtime) then
   write(string1, *) trim(file_name), ' mtimedim = ',dimlen, ' DART expects ', nmtime
   call error_handler(E_ERR,'read_TIEGCM_restart',string1,source,revision,revdate)
endif

!... get mtime
allocate(mtimetmp(dimlen, time_dimlen))
call nc_check(nf90_inq_varid(ncid, 'mtime', VarID), 'read_TIEGCM_restart', &
                  'inq_varid mtime')
call nc_check(nf90_get_var(ncid, VarID, values=mtimetmp), &
        'read_TIEGCM_restart', 'get_var mtime')
mtime = mtimetmp(:,time_dimlen)
deallocate(mtimetmp)

!... get year
allocate(yeartmp(time_dimlen))
call nc_check(nf90_inq_varid(ncid, 'year', var_year_id), 'read_TIEGCM_restart', &
                  'inq_varid year')
call nc_check(nf90_get_var(ncid, var_year_id, values=yeartmp), &
       'read_TIEGCM_restart', 'get_var year')
year = yeartmp(time_dimlen)
deallocate(yeartmp)

if (do_output()) print *, 'read_TIEGCM_restart: mtime (doy/hour/minute) and year:', mtime, year
doy   =  mtime(1)
utsec = (mtime(2)*60 + mtime(3))*60
model_time = set_time(utsec, doy-1) + set_date(year, 1, 1)  ! Jan 1 of whatever year.

if (do_output()) call print_time(model_time, str=' read_TIEGCM_restart: model_time ')
if (do_output()) call print_date(model_time, str=' read_TIEGCM_restart: model_date ')

call nc_check( nf90_close(ncid), 'update_TIEGCM_restart', 'close')

end subroutine read_TIEGCM_restart



subroutine update_TIEGCM_restart(file_name, var)
!-------------------------------------------------------------------------------
!
! Updates TIEGCM restart file fields

character(len=*), intent(in) :: file_name
real(r8),         intent(in) :: var

!TJH    integer                             :: ncerr
!TJH    integer                             :: ncid
!TJH    integer                             :: DimID, dimlen
!TJH    integer                             :: TimeDimID, time_dimlen
!TJH    integer                             :: VarID
!TJH    integer, parameter                  :: nmtime = 3
!TJH    integer, dimension(nmtime)          :: mtime  ! day, hour, minute
!TJH    integer                             :: utsec, doy !year
!TJH 
!TJH    real(r8), dimension(nlon,nlat,nlev) :: TN, TN_NM, UN, UN_NM, VN, VN_NM
!TJH    real(r8), dimension(nlon,nlat,nlev) :: O1, O1_NM, O2, O2_NM, OP
!TJH    real(r8), dimension(nlon,nlat,nilev):: NE
!TJH 
!TJH    integer                             :: nlevm1
!TJH    type(time_type)                     :: jan1, tbase
!TJH    integer                             :: year, month, day, hour, mins, sec
!TJH 
!TJH    if ( .not. module_initialized ) call static_init_model
!TJH 
!TJH    nlevm1 = nlev -1
!TJH 
!TJH    if( .not. file_exist(file_name)) then
!TJH       write(string1,*) trim(adjustl(file_name)),' not available.'
!TJH       call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
!TJH    endif
!TJH 
!TJH 
!TJH    if (do_output()) print *, 'update_TIEGCM_restart: opening restart'
!TJH    ncerr = nf90_open( file_name, NF90_WRITE, ncid )! open with read/write access
!TJH    call nc_check(ncerr, 'update_TIEGCM_restart','open')  ! will die if error
!TJH    if (do_output()) print *, 'update_TIEGCM_restart: opened with '//trim(nf90_strerror(ncerr))
!TJH 
!TJH 
!TJH    !... check for matching dimensions
!TJH    call nc_check( nf90_inq_dimid(ncid, 'lon', DimID), &
!TJH           'update_TIEGCM_restart', 'inq_dimid lon')
!TJH    call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
!TJH           'update_TIEGCM_restart', 'inquire_dimension lon')
!TJH    if (dimlen .ne. nlon) then
!TJH      write(string1, *) trim(file_name), ' dim_lon = ',dimlen, ' DART expects ',nlon
!TJH      call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
!TJH    endif
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_dimid(ncid, 'lat', DimID), &
!TJH           'update_TIEGCM_restart', 'inq_dimid lat')
!TJH    call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
!TJH           'update_TIEGCM_restart', 'inquire_dimension lat')
!TJH    if (dimlen .ne. nlat) then
!TJH      write(string1, *) trim(file_name), ' dim_lat = ',dimlen, ' DART expects ',nlat
!TJH      call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
!TJH    endif
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_dimid(ncid, 'lev', DimID), &
!TJH           'update_TIEGCM_restart', 'inq_dimid lev')
!TJH    call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
!TJH           'update_TIEGCM_restart', 'inquire_dimension lev')
!TJH    if (dimlen .ne. nlev) then
!TJH      write(string1, *) trim(file_name), ' dim_lev = ',dimlen, ' DART expects ',nlev
!TJH      call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
!TJH    endif
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_dimid(ncid, 'ilev', DimID), &
!TJH           'update_TIEGCM_restart', 'inq_dimid ilev')
!TJH    call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
!TJH           'update_TIEGCM_restart', 'inquire_dimension ilev')
!TJH    if (dimlen .ne. nilev) then
!TJH      write(string1, *) trim(file_name), ' dim_ilev = ',dimlen, ' DART expects ',nilev
!TJH      call error_handler(E_ERR,'update_TIEGCM_restart',string1,source,revision,revdate)
!TJH    endif
!TJH 
!TJH 
!TJH    call nc_check( nf90_inquire(ncid, unlimitedDimId = TimeDimID), &
!TJH           'update_TIEGCM_restart', 'inquire id of unlimited dimension time')
!TJH    call nc_check( nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
!TJH           'update_TIEGCM_restart', 'inquire_dimension time')
!TJH 
!TJH 
!TJH !... put variables into TIEGCM array
!TJH 
!TJH    TN(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,TYPE_local_TN+1)
!TJH    TN(:,:,  nlev)      = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH    TN_NM(:,:,1:nlevm1) = var%vars_3d(:,:,1:nlevm1,TYPE_local_TN_NM+1)
!TJH    TN_NM(:,:,   nlev ) = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH    O1                  = var%vars_3d(:,:,:,TYPE_local_O1+1)
!TJH    O1_NM               = var%vars_3d(:,:,:,TYPE_local_O1_NM+1)
!TJH    O2                  = var%vars_3d(:,:,:,TYPE_local_O2+1)
!TJH    O2_NM               = var%vars_3d(:,:,:,TYPE_local_O2_NM+1)
!TJH 
!TJH    if (.not. only_neutral_density) then
!TJH       UN(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,TYPE_local_UN+1)
!TJH       UN(:,:,   nlev)     = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH       UN_NM(:,:,1:nlevm1) = var%vars_3d(:,:,1:nlevm1,TYPE_local_UN_NM+1)
!TJH       UN_NM(:,:,nlev)     = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH       VN(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,TYPE_local_VN+1)
!TJH       VN(:,:,  nlev)      = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH       VN_NM(:,:,1:nlevm1) = var%vars_3d(:,:,1:nlevm1,TYPE_local_VN_NM+1)
!TJH       VN_NM(:,:,nlev)     = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH       OP(:,:,1:nlevm1)    = var%vars_3d(:,:,1:nlevm1,TYPE_local_OP+1)
!TJH       OP(:,:,nlev)        = TIEGCM_missing_value        !fill top slot with missing value
!TJH 
!TJH       NE                  = var%vars_3d(:,:,:,TYPE_local_NE+1)
!TJH    endif ! (.not. only_neutral_density)
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_varid(ncid, 'TN', VarID), &
!TJH           'update_TIEGCM_restart', 'inq_varid TN')
!TJH    call nc_check( nf90_put_var(ncid, VarID, values=TN, &
!TJH                   start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH           'update_TIEGCM_restart', 'put_var TN')
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_varid(ncid, 'TN_NM', VarID), &
!TJH           'update_TIEGCM_restart', 'inq_varid TN_NM')
!TJH    call nc_check( nf90_put_var(ncid, VarID, values=TN_NM, &
!TJH                   start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH           'update_TIEGCM_restart', 'put_var TN_NM')
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_varid(ncid, 'O1', VarID), &
!TJH           'update_TIEGCM_restart', 'inq_varid O1')
!TJH    call nc_check( nf90_put_var(ncid, VarID, values=O1, &
!TJH                   start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH           'update_TIEGCM_restart', 'put_var O1')
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_varid(ncid, 'O1_NM', VarID), &
!TJH           'update_TIEGCM_restart', 'inq_varid O1_NM')
!TJH    call nc_check( nf90_put_var(ncid, VarID, values=O1_NM, &
!TJH                   start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH           'update_TIEGCM_restart', 'put_var O1_NM')
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_varid(ncid, 'O2', VarID), &
!TJH           'update_TIEGCM_restart', 'inq_varid O2')
!TJH    call nc_check( nf90_put_var(ncid, VarID, values=O2, &
!TJH                   start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH           'update_TIEGCM_restart', 'put_var O2')
!TJH 
!TJH 
!TJH    call nc_check( nf90_inq_varid(ncid, 'O2_NM', VarID), &
!TJH           'update_TIEGCM_restart', 'inq_varid O2_NM')
!TJH    call nc_check( nf90_put_var(ncid, VarID, values=O2_NM, &
!TJH                   start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH           'update_TIEGCM_restart', 'put_var O2_NM')
!TJH 
!TJH    if (.not. only_neutral_density) then
!TJH       call nc_check( nf90_inq_varid(ncid, 'UN', VarID), &
!TJH              'update_TIEGCM_restart', 'inq_varid UN')
!TJH       call nc_check( nf90_put_var(ncid, VarID, values=UN, &
!TJH              start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH              'update_TIEGCM_restart', 'put_var UN')
!TJH 
!TJH 
!TJH       call nc_check( nf90_inq_varid(ncid, 'UN_NM', VarID), &
!TJH              'update_TIEGCM_restart', 'inq_varid UN_NM')
!TJH       call nc_check( nf90_put_var(ncid, VarID, values=UN_NM, &
!TJH              start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH              'update_TIEGCM_restart', 'put_var UN_NM')
!TJH 
!TJH       call nc_check( nf90_inq_varid(ncid, 'VN', VarID), &
!TJH              'update_TIEGCM_restart', 'inq_varid VN')
!TJH       call nc_check( nf90_put_var(ncid, VarID, values=VN, &
!TJH              start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH             'update_TIEGCM_restart', 'put_var VN')
!TJH 
!TJH       call nc_check( nf90_inq_varid(ncid, 'VN_NM', VarID), &
!TJH              'update_TIEGCM_restart', 'inq_varid VN_NM')
!TJH       call nc_check( nf90_put_var(ncid, VarID, values=VN_NM, &
!TJH              start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nlev,1/)), &
!TJH              'update_TIEGCM_restart', 'put_var VN_NM')
!TJH 
!TJH       call nc_check( nf90_inq_varid(ncid, 'OP', VarID), &
!TJH              'update_TIEGCM_restart', 'inq_varid OP')
!TJH       call nc_check( nf90_put_var(ncid, VarID, values=OP, &
!TJH              start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nilev,1/)), &
!TJH              'update_TIEGCM_restart', 'put_var OP')
!TJH 
!TJH       call nc_check( nf90_inq_varid(ncid, 'NE', VarID), &
!TJH              'update_TIEGCM_restart', 'inq_varid NE')
!TJH       call nc_check( nf90_put_var(ncid, VarID, values=NE, &
!TJH               start = (/1,1,1,time_dimlen/), count = (/nlon,nlat,nilev,1/)), &
!TJH               'update_TIEGCM_restart', 'put_var NE')
!TJH 
!TJH    endif ! (.not. only_neutral_density)
!TJH 
!TJH !... mtime and year
!TJH    call get_date(var%valid_time, year, month, day, hour, mins, sec )
!TJH    jan1  = set_date(year,1,1)
!TJH    tbase = var%valid_time - jan1    ! total time since the start of the year.
!TJH 
!TJH    call get_time(tbase, utsec, doy)
!TJH 
!TJH    mtime(1) = doy + 1  ! Have to add January 1 back in
!TJH    mtime(2) = hour
!TJH    mtime(3) = mins
!TJH 
!TJH    if (do_output()) print *, 'update_TIEGCM_restart: mtime (doy/hour/minute):', mtime
!TJH    call nc_check( nf90_inq_varid(ncid, 'mtime', VarID), &
!TJH           'update_TIEGCM_restart','inq_varid mtime')
!TJH    call nc_check( nf90_put_var( ncid, VarID, values=mtime, &
!TJH                   start = (/1,time_dimlen/), count = (/nmtime,1/)), &
!TJH           'update_TIEGCM_restart','get_var mtime')
!TJH 
!TJH    if (do_output()) print *, 'update_TIEGCM_restart: year:', year
!TJH    call nc_check( nf90_inq_varid(ncid, 'year', VarID), &
!TJH           'update_TIEGCM_restart','inq_varid year')
!TJH    call nc_check( nf90_put_var( ncid, VarID, values=year, &
!TJH                   start = (/time_dimlen/)) , &
!TJH           'update_TIEGCM_restart','get_var year')
!TJH 
!TJH    call nc_check( nf90_sync(ncid), 'update_TIEGCM_restart', 'sync')
!TJH 
!TJH    call nc_check( nf90_close(ncid), 'update_TIEGCM_restart', 'close')
!TJH 
end subroutine update_TIEGCM_restart


function get_restart_file_name()
character(len=128) :: get_restart_file_name

if ( .not. module_initialized ) call static_init_model
get_restart_file_name = adjustl(tiegcm_restart_file_name)

end function get_restart_file_name



!===============================================================================
! routines below here are private to the module
!===============================================================================



subroutine read_TIEGCM_namelist(file_name)
!-------------------------------------------------------------------------------
! Read in TIEGCM namelist input
!
! under certain situations, the value of f10.7 is a parameter to be estimated
! and needs to be added to state vector

character(len=*), intent(in) :: file_name
integer  :: iunit, io
integer  :: daysec = 86400

!-------------------------------------------------------------------------------
! 1/3/2011, the namelist definition taken from $TGCMROOT/tiegcm1.93/src/input.F
!           the following parameter values are from params.F
!           modify the namelist definition for future tiegcm updates

integer,parameter :: mxind_time = 500 ! max number of time-dependent solar index points
integer,parameter :: mxhvols = 100    ! max number of output history file
integer,parameter :: mxseries = 10    ! max number of time series for primary histories
integer,parameter :: mxfsech = 100    ! max number of fields on secondary histories

! Namelist user input variables:

character(len=80):: &
     &  label,           &! optional generic text label for this run
     &  tempdir,         &! temporary directory
     &  magvol,          &! file name or mss path to magnetic data file
     &  amievol           ! file or mss path of amie data file (optional)

! date and calday are no longer supported, and are replaced by start_day,
! start_year, and calendar_advance. Date and calday are retained here so
! error usage statements can be issued if user sets one of them.

integer :: &
     &  start_day,       &! starting day of year (integer 0->365)
     &  start_year,      &! starting year (4-digit integer yyyy)
     &  calendar_advance,&! if > 0, advance calendar day from start_day
     &  date(3),         &! old: model starting year, day ( 2 ints yyyy,dd)
     &  calday,          &! old: starting calendar day (0-mxday)
     &  mxday,           &! calendar day (0-mxday)
     &  step,            &! model time step (integer seconds)
     &  dispose,         &! dispose output files to mss if dispose==1 or 2
     &  eddy_dif,        &! 0/1 flag for DOY dependent eddy diffusion (difk, dift, xmue)
     &  dynamo,          &! 0/1 flag for dynamo
     &  tideann,         &! 0/1 flag for annual tide (deprecated as of May 2008)
     &  aurora,          &! 0/1 flag for aurora
     &  ntask_lat,       &! number of tasks in latitude  dimension
     &  ntask_lon         ! number of tasks in longitude dimension
real :: &
     &  tide(10),        &! semidiurnal tide amplitudes and phases
     &  tide2(2),        &! diurnal tide amplitude and phase
     &  tide3m3(2),      &! 2-day wave amplitude and phase
     &  f107,            &! 10.7 cm daily solar flux
     &  f107a,           &! 10.7 cm average (81-day) solar flux
     &  colfac            ! collision factor

! Input parameters that can be either constant or time-dependent:
real :: &
     &  power,           &! hemispheric power (gw) (hpower on histories)
     &  ctpoten,         &! cross-cap potential (volts)
     &  bximf,           &! BX component of IMF
     &  byimf,           &! BY component of IMF
     &  bzimf,           &! BZ component of IMF in nT
     &  swvel,           &! Solar wind velocity in km/s
     &  swden,           &! Solar wind density in #/cm3
     &  al,              &! AL lower magnetic auroral activity index in nT
     &  kp                ! Kp index
real,dimension(4,mxind_time) :: power_time,ctpoten_time,           &
     &  bximf_time,byimf_time,bzimf_time,swvel_time,swden_time,al_time,  &
     &  kp_time
integer :: &
     &  ntimes_ctpoten,ntimes_power,ntimes_bximf,ntimes_byimf,           &
     &  ntimes_bzimf,ntimes_swden,ntimes_swvel,ntimes_al,ntimes_kp
logical :: aluse    ! logical to use AL in Weimer 2001 model or not

! Parameters as read from namelist:
real :: rd_power,rd_ctpoten,rd_f107,rd_f107a,rd_bximf,rd_byimf,    &
     &  rd_bzimf,rd_swvel,rd_swden,rd_kp
!
! If indices_interp==1, time-dependent indices (power_time, ctpoten_time, etc)
! will be interpolated to model time, otherwise they will change only
! when the given values change. This has no effect on indices given as constants.

integer :: indices_interp=1

! Import data file names:

integer,parameter :: mxlen_filename=80
character(len=mxlen_filename) ::                                   &
!
! 4/2/08 btf: Introducing Weimer 2005 model (wei05sc.F).
!             Retain ability to call either the 2001 or 2005 weimer models
!             for now, to facilitate comparison runs, so potential_model
!             can be either WEIMER01 or WEIMER05.
!
     &  potential_model,   &! electric potential model used
                            ! Values can be 'HEELIS', 'WEIMER', or 'NONE'
                            ! If absent, the default value is set to 'HEELIS'
     &  weimer_ncfile,     &! path to netcdf weimer01 coefficients file
     &  wei05sc_ncfile,    &! path to netcdf data files for weimer05 model
     &  gpi_ncfile,        &! mss path or file path to netcdf gpi data file
     &  ncep_ncfile,       &! ncep data file (time-gcm only)
     &  see_ncfile,        &! mss path or file path to netcdf SEE flux data file
     &  imf_ncfile,        &! mss path or disk file path to netcdf IMF data file
     &  gswm_mi_di_ncfile, &! gswm migrating diurnal data file
     &  gswm_mi_sdi_ncfile,&! gswm migrating semi-diurnal data file
     &  gswm_nm_di_ncfile, &! gswm non-migrating diurnal data file
     &  gswm_nm_sdi_ncfile,&! gswm non-migrating semi-diurnal data file
     &  saber_ncfile,      &! SABER data (T,Z)
     &  tidi_ncfile         ! TIDI data (U,V)
!
!     integer,parameter :: ngpivars = 4
!     real :: gpi_vars(ngpivars) ! f107,f107a,power,ctpoten
!     character(len=16) ::
!    |  gpi_names(ngpivars)      ! names of gpi_vars

! Primary history user input (dimension parameters are in params.h):
character(len=80) :: &
     &  source,            &! file containing source history (optional)
        output(mxhvols)     ! output file(s) (required)
integer ::           &
     &  source_start(3),   &! source history model time
     &  start(3,mxseries), &! primary history model start time(s)
     &  stop(3,mxseries),  &! primary history model stop time(s)
     &  hist(3,mxseries),  &! primary history disk write frequency
     &  save(3,mxseries),  &! primary history file save frequency
     &  mxhist_prim,       &! max number of histories per primary file
     &  msreten,           &! retention period for history files
     &  noutput             ! number of output files given
!
! Secondary history user input (dimension parameters are in params.h):
character(len=80) ::   &
     &  secsource,           &! file containing source sec_history (for mhd)
     &  secout(mxhvols)       ! secondary history output file(s)
character(len=16) ::   &
     &  secflds(mxfsech)      ! secondary history output fields
integer ::             &
     &  secstart(3,mxseries),&! secondary history model start time(s)
     &  secstop(3,mxseries), &! secondary history model stop time(s)
     &  sechist(3,mxseries), &! secondary history disk write frequency
     &  secsave(3,mxseries), &! secondary history file save frequency
     &  mxhist_sech,         &! max number of histories per secondary file
     &  sech_nbyte            ! 4 or 8: write real or double values to secondary file
!
! Namelist for read:
namelist/tgcm_input/                                        &
     &  label,tempdir,magvol,amievol,date,calday,step,dispose,    &
     &  source,source_start,output,start,stop,hist,save,          &
     &  secout,secstart,secstop,sechist,secsave,secflds,          &
     &  potential_model,eddy_dif,dynamo,tide,tide2,tide3m3,       &
     &  f107,f107a,power,ctpoten,bximf,byimf,bzimf,swvel,swden,al,&
     &  kp,colfac,tideann,aurora,gpi_ncfile,gswm_mi_di_ncfile,    &
     &  gswm_mi_sdi_ncfile,gswm_nm_di_ncfile,gswm_nm_sdi_ncfile,  &
     &  mxhist_prim,mxhist_sech,msreten,ntask_lat,ntask_lon,      &
     &  start_day,start_year,calendar_advance,see_ncfile,         &
     &  ctpoten_time,power_time,bximf_time,byimf_time,bzimf_time, &
     &  kp_time,al_time,swden_time,swvel_time,indices_interp,     &
     &  imf_ncfile,saber_ncfile,tidi_ncfile,sech_nbyte

!-------------------------------------------------------------------------------

if( .not. file_exist(file_name)) then
   write(string1,*) trim(adjustl(file_name)),' not available.'
   call error_handler(E_ERR,'read_TIEGCM_namelist',string1,source,revision,revdate)
endif

if (do_output()) print *, 'read_TIEGCM_namelist: reading restart:', file_name

! Read the namelist entry tgcm_input from tiegcm.nml
! TJH It would be nice to read the namelist and skip all the ';' in column 1.
! Basically, we are just getting the value of f10.7 and saving it.

call find_namelist_in_file('tiegcm.nml', 'tgcm_input', iunit)
read(iunit, nml = tgcm_input, iostat = io)
call check_namelist_read(iunit, io, 'tgcm_input')

if (step >= daysec) then
    time_step_days    = int(step/daysec)
    time_step_seconds = mod(step,daysec)
else
    time_step_days    = 0
    time_step_seconds = step
endif

f10_7 = f107  ! save this in module storage

end subroutine read_TIEGCM_namelist



subroutine read_TIEGCM_definition(file_name)
!-------------------------------------------------------------------------------
! Read TIEGCM grid definition and Geopotential from a tiegcm restart file
! fills metadata storage variables:
! lons(:), nlon
! lats(:), nlat
! lev(:),  nlev
! ilev(:), nilev
! plevs(:)
! pilevs(:)

character(len=*), intent(in) :: file_name
integer  :: ncid, VarID, DimID, dimlen
integer  :: TimeDimID, time_dimlen, missing_value_len
real(r8) :: p0

integer, dimension(:,:), allocatable     :: mtimetmp
integer, parameter                       :: nmtime = 3
integer, dimension(nmtime)               :: mtime  ! day, hour, minute
integer                                  :: utsec, doy

if( .not. file_exist(file_name)) then
  write(string1,*) trim(adjustl(file_name)),' not available.'
  call error_handler(E_ERR,'read_TIEGCM_definition',string1,source,revision,revdate)
endif

if (do_output()) print *, 'read_TIEGCM_definition:reading restart:', file_name

call nc_check(nf90_open(file_name, NF90_NOWRITE, ncid), &
       'read_TIEGCM_definition','open '//trim(file_name))

! longitude

call nc_check(nf90_inq_dimid(ncid, 'lon', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid lon')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nlon), 'read_TIEGCM_definition', &
                  'inquire_dimension lon')
allocate(lons(nlon))
call nc_check(nf90_inq_varid(ncid, 'lon', VarID), 'read_TIEGCM_definition', &
                  'inq_varid lon')
call nc_check(nf90_get_var(ncid, VarID, values=lons), 'read_TIEGCM_definition', &
                  'get_var lon')

where (lons < 0.0_r8) lons = lons + 360.0_r8 ! DART [0, 360] TIEGCM [-180, 180]

! latitiude

call nc_check(nf90_inq_dimid(ncid, 'lat', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid lat')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nlat), 'read_TIEGCM_definition', &
                  'inquire_dimension lat')
allocate(lats(nlat))
call nc_check(nf90_inq_varid(ncid, 'lat', VarID), 'read_TIEGCM_definition', &
                  'inq_varid lat')
call nc_check(nf90_get_var(ncid, VarID, values=lats), 'read_TIEGCM_definition', &
                  'get_var lat')

! pressure

call nc_check(nf90_inq_varid(ncid, 'p0', VarID), 'read_TIEGCM_definition', &
                  'inq_varid p0')
call nc_check(nf90_get_var(ncid, VarID, values=p0), 'read_TIEGCM_definition', &
                  'get_var p0')

TIEGCM_reference_pressure = p0

call nc_check(nf90_inq_dimid(ncid, 'lev', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid lev')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nlev), 'read_TIEGCM_definition', &
                  'inquire_dimension lev')

allocate(levs(nlev), plevs(nlev))

call nc_check(nf90_inq_varid(ncid, 'lev', VarID), 'read_TIEGCM_definition', &
                  'inq_varid lev')
call nc_check(nf90_get_var(ncid, VarID, values=levs), 'read_TIEGCM_definition', &
                  'get_var lev')

plevs = p0 * exp(-levs) * 100.0_r8 ![Pa] = 100* [millibars] = 100* [hPa]

call nc_check(nf90_inq_dimid(ncid, 'ilev', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid ilev')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nilev), 'read_TIEGCM_definition', &
                  'inquire_dimension ilev')

allocate(ilevs(nilev), pilevs(nilev))

call nc_check(nf90_inq_varid(ncid, 'ilev', VarID), 'read_TIEGCM_definition', &
                  'inq_varid ilev')
call nc_check(nf90_get_var(ncid, VarID, values=ilevs), 'read_TIEGCM_definition', &
                  'get_var ilev')

pilevs = p0 * exp(-ilevs) * 100.0_r8 ! [Pa] = 100* [millibars] = 100* [hPa]

if (nlev .ne. nilev) then
  write(string1, *) ' nlev = ',nlev,' nilev = ',nilev, &
                      'are different; DART assumes them to be the same'
  call error_handler(E_ERR,'read_TIEGCM_definition',string1,source,revision,revdate)
endif

! TJH FIXME ... these appear to be unused from here to the end

call nc_check(nf90_inquire_attribute(ncid, nf90_global, 'missing_value', &
    len=missing_value_len), 'read_TIEGCM_definition', 'inquire attribute "missing_value"')
if (missing_value_len /= 1) then
  write(string1, *) ' global attribute missing_value length is ', missing_value_len, ' DART expects 1'
  call error_handler(E_ERR,'read_TIEGCM_definition',string1,source,revision,revdate)
endif

call nc_check(nf90_get_att(ncid, nf90_global, 'missing_value', TIEGCM_missing_value), &
       'read_TIEGCM_definition', 'get_att global attribute named missing_value')

!... get mtime
call nc_check(nf90_inq_dimid(ncid, 'mtimedim', DimID), 'read_TIEGCM_definition', &
                  'inq_dimid mtimedim')
call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), 'read_TIEGCM_definition', &
                  'inquire_dimension mtimedim')
if (dimlen .ne. nmtime) then
  write(string1, *) trim(file_name), ' mtimedim = ',dimlen, ' DART expects ', nmtime
  call error_handler(E_ERR,'read_TIEGCM_definition',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire(ncid, unlimitedDimId = TimeDimID), 'read_TIEGCM_definition', &
                  'inquire id of unlimited dimension time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), 'read_TIEGCM_definition', &
                  'inquire_dimension time')

allocate(mtimetmp(dimlen, time_dimlen))
call nc_check(nf90_inq_varid(ncid, 'mtime', VarID), 'read_TIEGCM_definition', &
                  'inq_varid mtime')
call nc_check(nf90_get_var(ncid, VarID, values=mtimetmp), 'read_TIEGCM_definition', &
                  'get_var mtime')
mtime = mtimetmp(:,time_dimlen)
deallocate(mtimetmp)

call nc_check(nf90_close(ncid),'read_TIEGCM_definition', 'close')

if (do_output()) print *, 'read_TIEGCM_definition: mtime (doy/hour/minute):', mtime

doy   =  mtime(1)
utsec = (mtime(2)*60 + mtime(3))*60

end subroutine read_TIEGCM_definition



subroutine read_TIEGCM_secondary(file_name)
!-------------------------------------------------------------------------------
!
! Read TIEGCM Geopotential (ZG) from a tiegcm secondary output file

character(len=*), intent(in):: file_name

integer :: ncid
integer :: VarID, DimID, dimlen
integer :: TimeDimID, time_dimlen

real(r8) :: spvalR8

if( .not. file_exist(file_name)) then
  write(string1,*) trim(adjustl(file_name)),' not available.'
  call error_handler(E_ERR,'read_TIEGCM_secondary',string1,source,revision,revdate)
endif

if (do_output()) print *, 'read_TIEGCM_secondary:reading restart:', file_name

call nc_check(nf90_open(file_name, NF90_NOWRITE, ncid), &
                        'read_TIEGCM_secondary', 'open')

!... check for matching dimensions
call nc_check( nf90_inq_dimid(ncid, 'lon', DimID), &
       'read_TIEGCM_secondary', 'inq_dimid lon')
call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
       'read_TIEGCM_secondary', 'inquire_dimension lon')
if (dimlen .ne. nlon) then
  write(string1, *) trim(file_name), ' dim_lon = ',dimlen, ' DART expects ',nlon
  call error_handler(E_ERR,'read_TIEGCM_secondary',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'lat', DimID), &
       'read_TIEGCM_secondary', 'inq_dimid lat')
call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
       'read_TIEGCM_secondary', 'inquire_dimension lat')
if (dimlen .ne. nlat) then
  write(string1, *) trim(file_name), ' dim_lat = ',dimlen, ' DART expects ',nlat
  call error_handler(E_ERR,'read_TIEGCM_secondary',string1,source,revision,revdate)
endif

call nc_check( nf90_inq_dimid(ncid, 'ilev', DimID), &
       'read_TIEGCM_secondary', 'inq_dimid ilev')
call nc_check( nf90_inquire_dimension(ncid, DimID, len=dimlen), &
       'read_TIEGCM_secondary', 'inquire_dimension ilev')
if (dimlen .ne. nilev) then
  write(string1, *) trim(file_name), ' dim_ilev = ',dimlen, ' DART expects ',nilev
  call error_handler(E_ERR,'read_TIEGCM_secondary',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire(ncid, unlimitedDimId = TimeDimID), &
       'read_TIEGCM_secondary', 'inquire id of unlimited dimension time')
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
       'read_TIEGCM_secondary', 'inquire_dimension time')

allocate(ZG(nlon,nlat,nilev)) ! comes from module storage

!... actually read the target variable
call nc_check(nf90_inq_varid(ncid, 'ZG', VarID),    &
       'read_TIEGCM_secondary', 'inq_varid ZG')
call nc_check(nf90_get_var(ncid, VarID, values=ZG,  &
                         start = (/ 1, 1, 1, time_dimlen /),    &
                         count = (/ nlon, nlat, nilev, 1 /)),    &
                         'read_TIEGCM_secondary', 'get_var ZG')

if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
   where(ZG == spvalR8) ZG = MISSING_R8
endif

! FIXME ... check units and convert them to meters if need be.
!if (nf90_get_att(ncid, VarID, 'units' , string1) == NF90_NOERR) then
!   if(trim(string1) == 'cm') then
!      where(ZG /= MISSING_R8) ZG = ZG/100.0_r8
!   endif
!endif

call nc_check(nf90_close(ncid),'read_TIEGCM_secondary', 'close')

end subroutine read_TIEGCM_secondary





subroutine verify_state_variables( state_variables, filename, ngood, table )
! This routine checks the user input against the variables available in the
! input netcdf file to see if it is possible to construct the DART state vector
! specified by the input.nml:model_nml:state_variables  variable.

character(len=*), dimension(:),   intent(in)  :: state_variables
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer  :: i, nrows, ncols
integer  :: ivar, index1, indexN, varsize
integer  :: ncid, VarID, dimlen
real(r8) :: spvalR8

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=obstypelength) :: dimname

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
              'verify_state_variables','open '//trim(filename))

nrows = size(table,1)
ncols = size(table,2)

! Convert the (input) 1D array "state_varibles" into a table with two colums.
! The number of rows in the table correspond to the number of variables in the
! DART state vector.

ngood = 0
MyLoop : do i = 1, nrows

   varname    = trim(state_variables(2*i -1))
   dartstr    = trim(state_variables(2*i   ))
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or.  table(i,2) == ' ' ) then
      string1 = 'model_nml:state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in netCDF file, as long as it is not vTEC.
   ! vTEC gets constructed. 

   if (varname == 'VTEC') then
      include_vTEC_in_state = .true.
   else
      write(string1,'(''No variable '',a,'' in '',a)') trim(varname), trim(filename)
      call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string1))
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''No obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 3) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
   endif

   ngood = ngood + 1

enddo MyLoop

if (ngood == (nrows-2)) then
   string1 = 'WARNING: There is a possibility you need to increase ''StateVarNmax'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG, 'verify_state_variables', string1, &
                      source, revision, revdate, text2=string2)
endif

! Augment the state vector with the geopotential height
ngood = ngood + 1
table(ngood,1) = 'ZG'
table(ngood,2) = 'KIND_GEOPOTENTIAL_HEIGHT'

! Do we need to augment the state vector with the parameter to estimate?
if ( estimate_f10_7 ) then
   ngood = ngood + 1
   table(ngood,1) = 'f10_7'
   table(ngood,2) = 'KIND_1D_PARAMETER'
endif

!-------------------------------------------------------------------------------
! Now that we know how many variables, etc., fill metadata structure 'progvar'
! read_TIEGCM_restart() uses this structure to figure out what to put in the
! DART state vector. 
!-------------------------------------------------------------------------------

index1  = 1;
indexN  = 0;
FillLoop : do ivar = 1, ngood

   varname                   = trim(table(ivar,1))
   progvar(ivar)%varname     = varname
   progvar(ivar)%long_name   = 'blank'
   progvar(ivar)%units       = 'blank'
   progvar(ivar)%xtype       = -1
   progvar(ivar)%storder     = 'blank'
   progvar(ivar)%dimlens     = -1
   progvar(ivar)%posdef      = -1
   progvar(ivar)%rank     = -1
   progvar(ivar)%varsize     = -1
   progvar(ivar)%index1      = -1
   progvar(ivar)%indexN      = -1
   progvar(ivar)%kind_string = trim(table(ivar,2))
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%verticalvar = 'undefined'
   progvar(ivar)%rangeRestricted = -1
   progvar(ivar)%maxvalue    = MISSING_R8
   progvar(ivar)%minvalue    = MISSING_R8
   progvar(ivar)%missingR8   = MISSING_R8
   progvar(ivar)%has_missing_value = .false.
   progvar(ivar)%prognostic  = .true.
   progvar(ivar)%dimnames    = 'blank'

   if (trim(varname) == 'ZG') then
      ivarZG = ivar
      ! must come from auxiliary file.  sort-of-cheating by hardwiring ...
      varsize = nlon*nlat*nilev*1
      progvar(ivar)%long_name         = 'Geometric height'
      progvar(ivar)%units             = 'cm'
      progvar(ivar)%xtype             = NF90_DOUBLE
      progvar(ivar)%missingR8         = 1.e+36
      progvar(ivar)%has_missing_value = .true.
      progvar(ivar)%prognostic        = .false.  ! do not update - localize far away
      progvar(ivar)%verticalvar       = 'ilev'
      progvar(ivar)%dimlens(1:4)      = (/nlon,  nlat,  nilev,      1 /)
      progvar(ivar)%dimnames(1:4)     = (/'lon ','lat ','ilev', 'time'/)
      progvar(ivar)%rank              = 4
      progvar(ivar)%varsize           = varsize
      progvar(ivar)%index1            = index1
      progvar(ivar)%indexN            = index1 + varsize - 1
      index1                          = index1 + varsize ! sets up for next variable
      cycle FillLoop
   endif

   if (trim(varname) == 'f10_7') then
      ! parameter to estimate, value comes from tiegcm.nml
      varsize = 1
      progvar(ivar)%long_name         = '10.7 cm daily solar flux'
      progvar(ivar)%units             = 'none'
      progvar(ivar)%xtype             = NF90_DOUBLE
      progvar(ivar)%missingR8         = MISSING_R8
      progvar(ivar)%dimlens(1)        = 0
      progvar(ivar)%dimnames(1)       = 'parameter'
      progvar(ivar)%rank              = 0
      progvar(ivar)%varsize           = varsize
   !  progvar(ivar)%rangeRestricted   = -1
   !  progvar(ivar)%maxvalue          = MISSING_R8
   !  progvar(ivar)%minvalue          = MISSING_R8
      progvar(ivar)%index1            = index1
      progvar(ivar)%indexN            = index1 + varsize - 1
      index1                          = index1 + varsize ! sets up for next variable

      cycle FillLoop
   endif

   if (trim(varname) == 'VTEC') then
      ! 2D variable not in netCDF file, but useful to be part of state.
      varsize = nlon*nlat*1
      progvar(ivar)%long_name         = 'Total Electron Content'
      progvar(ivar)%units             = 'TECU'
      progvar(ivar)%xtype             = NF90_DOUBLE
      progvar(ivar)%missingR8         = MISSING_R8
      progvar(ivar)%has_missing_value = .false.
   !  progvar(ivar)%rangeRestricted   = -1
   !  progvar(ivar)%maxvalue          = MISSING_R8
   !  progvar(ivar)%minvalue          = MISSING_R8
      progvar(ivar)%dimlens(1:3)      = (/nlon,  nlat,      1 /)
      progvar(ivar)%dimnames(1:3)     = (/'lon ','lat ','time'/)
      progvar(ivar)%rank              = 3
      progvar(ivar)%varsize           = varsize
      progvar(ivar)%index1            = index1
      progvar(ivar)%indexN            = index1 + varsize - 1
      index1                          = index1 + varsize ! sets up for next variable

      cycle FillLoop
   endif

   ! Now that the "special" cases are handled, just read the information
   ! from the input netCDF file and insert into our structure.

   string2 = trim(filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'verify_state_variables', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%rank, xtype=progvar(ivar)%xtype), &
            'verify_state_variables', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'verify_state_variables', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'verify_state_variables', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Saving any FillValue so I can use it when I read and write ...

   if (progvar(ivar)%xtype == NF90_DOUBLE) then
      if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
         progvar(ivar)%missingR8         = spvalR8
         progvar(ivar)%has_missing_value = .true.
      endif
   else
      ! FIXME ... do some error checking ... do we support variables other than 'double'
   endif

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%rank

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'verify_state_variables', string1)
      ! read_TIEGCM_restart only reads the latest time step, so no matter how many
      ! time steps are defined, the time dimension is only 1
      if (trim(dimname) == 'time') dimlen = 1

      ! record what kind of vertical coordinate system is used for this variable
      if (trim(dimname) ==  'lev') progvar(ivar)%verticalvar =  'lev'
      if (trim(dimname) == 'ilev') progvar(ivar)%verticalvar = 'ilev'

      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

enddo FillLoop

ReportLoop : do ivar = 1, ngood
   if ((debug > 8) .and. do_output()) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name         ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units             ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype             ',progvar(ivar)%xtype
      write(logfileunit,*) '  dimnames          ',progvar(ivar)%dimnames(1:progvar(ivar)%rank)
      write(logfileunit,*) '  dimlens           ',progvar(ivar)%dimlens( 1:progvar(ivar)%rank)
      write(logfileunit,*) '  numdims           ',progvar(ivar)%rank
      write(logfileunit,*) '  varsize           ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1            ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN            ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind         ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string       ',trim(progvar(ivar)%kind_string)
      write(logfileunit,*) '  xtype             ',progvar(ivar)%xtype
      write(logfileunit,*) '  rangeRestricted   ',progvar(ivar)%xtype
      write(logfileunit,*) '  maxvalue          ',progvar(ivar)%maxvalue
      write(logfileunit,*) '  minvalue          ',progvar(ivar)%minvalue
      write(logfileunit,*) '  missingR8         ',progvar(ivar)%missingR8
      write(logfileunit,*) '  has_missing_value ',progvar(ivar)%has_missing_value

      write(    *      ,*)
      write(    *      ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(    *      ,*) '  long_name         ',trim(progvar(ivar)%long_name)
      write(    *      ,*) '  units             ',trim(progvar(ivar)%units)
      write(    *      ,*) '  xtype             ',progvar(ivar)%xtype
      write(    *      ,*) '  dimnames          ',progvar(ivar)%dimnames(1:progvar(ivar)%rank)
      write(    *      ,*) '  dimlens           ',progvar(ivar)%dimlens( 1:progvar(ivar)%rank)
      write(    *      ,*) '  numdims           ',progvar(ivar)%rank
      write(    *      ,*) '  varsize           ',progvar(ivar)%varsize
      write(    *      ,*) '  index1            ',progvar(ivar)%index1
      write(    *      ,*) '  indexN            ',progvar(ivar)%indexN
      write(    *      ,*) '  dart_kind         ',progvar(ivar)%dart_kind
      write(    *      ,*) '  kind_string       ',trim(progvar(ivar)%kind_string)
      write(    *      ,*) '  xtype             ',progvar(ivar)%xtype
      write(    *      ,*) '  rangeRestricted   ',progvar(ivar)%xtype
      write(    *      ,*) '  maxvalue          ',progvar(ivar)%maxvalue
      write(    *      ,*) '  minvalue          ',progvar(ivar)%minvalue
      write(    *      ,*) '  missingR8         ',progvar(ivar)%missingR8
      write(    *      ,*) '  has_missing_value ',progvar(ivar)%has_missing_value
   endif

enddo ReportLoop

call nc_check(nf90_close(ncid),'verify_state_variables','close '//trim(filename))

end subroutine verify_state_variables




subroutine fill_top(varname,temp4D,dimIDs,LevDimID,iLevDimID)
! Some variables need to have the top level replaced by the level underneath it.
! At present, I do not know of a rule to identify which varibles.
! It is _not_ as simple as the variables on 'lev' get theirs replaced while
! the variables on 'ilev' do not. 
!
! perhaps we need to check the top level of everyone for all missing values
! and key on that ... TJH

character(len=*),             intent(in)    :: varname
real(r8), dimension(:,:,:,:), intent(inout) :: temp4D
integer, dimension(:),        intent(in)    :: dimIDs
integer,                      intent(in)    :: LevDimID, iLevDimID

integer :: i, nlevels, nlevelsm1
integer :: levdim

levdim = 0

! only certain variables get this treatment
if ((varname(1:2) == 'TN') .or. &
    (varname(1:2) == 'UN') .or. &
    (varname(1:2) == 'VN') .or. &
    (varname(1:2) == 'OP') ) then

   DIMLoop : do i = 1,size(dimIDs)
      if ((dimIDs(i) == LevDimID) .or. (dimIDs(i) == iLevDimID)) then
         levdim = i
         exit DIMLoop
      endif
   enddo DIMLoop

   if (levdim == 0) then
      ! TJH FIXME we have a problem
   elseif (levdim == 1) then
      nlevels   = size(temp4D,levdim)
      nlevelsm1 = nlevels - 1
      temp4D(nlevels,:,:,:) = temp4D(nlevelsm1,:,:,:)

   elseif (levdim == 2) then
      nlevels   = size(temp4D,levdim)
      nlevelsm1 = nlevels - 1
      temp4D(:,nlevels,:,:) = temp4D(:,nlevelsm1,:,:)

   elseif (levdim == 3) then
      nlevels   = size(temp4D,levdim)
      nlevelsm1 = nlevels - 1
      temp4D(:,:,nlevels,:) = temp4D(:,:,nlevelsm1,:)

   elseif (levdim == 4) then
      nlevels   = size(temp4D,levdim)
      nlevelsm1 = nlevels - 1
      temp4D(:,:,:,nlevels) = temp4D(:,:,:,nlevelsm1)
   endif

endif

end subroutine fill_top



subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)
!-------------------------------------------------------------------------------
! I am trying to preserve the original shape of the variable as much as possible.
!
! the netCDF declarations look like : variable(time,       level, lat, lon) becomes
!                                     variable(time, copy, level, lat, lon)
!
! Since 'time' or 'Time' is the unlimited dimension in both ... I can skip it
! in the DEFDIM loop.

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i, mydimid

ndims = 0

DEFDIM : do i = 1,progvar(ivar)%rank

   if ((trim(progvar(ivar)%dimnames(i)) == 'Time') .or. &
       (trim(progvar(ivar)%dimnames(i)) == 'time')) cycle DEFDIM

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
                           'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid

enddo DEFDIM

ndims = ndims + 1               ! The next-to-last dimension is 'copy'
dimids(ndims) = memberdimid
ndims = ndims + 1               ! The last dimension is unlimited == time
dimids(ndims) = unlimitedDimid

if ( (do_output()) .and. debug > 7 ) then
   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   write(logfileunit,*)' thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%rank)
   write(     *     ,*)' thus dimids ',dimids(1:ndims)
endif

return
end subroutine define_var_dims



subroutine create_vtec( ncid, last_time, vTEC)
! Create the vTEC from constituents in the netCDF file.

integer,                  intent(in)  :: ncid
integer,                  intent(in)  :: last_time
real(r8), dimension(:,:), intent(out) :: vTEC

real(r8), allocatable, dimension(:,:,:) :: NE, Z, TI, TE, OP
real(r8), allocatable, dimension(:,:,:) :: NEm_extended, ZGG, ZGG_extended
real(r8), allocatable, dimension(:,:)   :: GRAVITYtop, Tplasma, Hplasma
real(r8), allocatable, dimension(:)     :: delta_ZGG, NE_middle

real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
real(r8), PARAMETER :: omass      = 2.678e-26_r8 ! mass of atomic oxgen kg

real(r8) :: earth_radiusm
integer  :: VarID, nlev10, j, k

allocate( NE(nlon,nlat,nilev), NEm_extended(nlon,nlat,nilev+10), &
         ZGG(nlon,nlat,nilev), ZGG_extended(nlon,nlat,nilev+10), &
           Z(nlon,nlat,nilev))
allocate( TI(nlon,nlat,nlev), TE(nlon,nlat,nlev), OP(nlon,nlat,nlev) )
allocate( GRAVITYtop(nlon,nlat), Tplasma(nlon,nlat), Hplasma(nlon,nlat) )
allocate( delta_ZGG(nlev+9), NE_middle(nlev+9) )

!... NE (interfaces)
call nc_check(nf90_inq_varid(ncid, 'NE', VarID), 'create_vtec', 'inq_varid NE')
call nc_check(nf90_get_var(ncid, VarID, values=NE,     &
                   start = (/ 1, 1, 1, last_time /),      &
                   count = (/ nlon, nlat, nilev, 1 /)),      &
                   'create_vtec', 'get_var NE') 

!... Z (interfaces)
call nc_check(nf90_inq_varid(ncid, 'Z', VarID), 'create_vtec', 'inq_varid Z')
call nc_check(nf90_get_var(ncid, VarID, values=Z,      &
                   start = (/ 1, 1, 1, last_time /),      &
                   count = (/ nlon, nlat, nilev, 1 /)),      &
                   'create_vtec', 'get_var Z') 

!... TI (midpoints)
call nc_check(nf90_inq_varid(ncid, 'TI', VarID), 'create_vtec', 'inq_varid TI')
call nc_check(nf90_get_var(ncid, VarID, values=TI,     &
                   start = (/ 1, 1, 1, last_time /),      &
                   count = (/ nlon, nlat, nlev, 1 /)),       &                            
                   'create_vtec', 'get_var TI') 

!... TE (midpoints)
call nc_check(nf90_inq_varid(ncid, 'TE', VarID), 'create_vtec', 'inq_varid TE')
call nc_check(nf90_get_var(ncid, VarID, values=TE,     &
                   start = (/ 1, 1, 1, last_time /),      &
                   count = (/ nlon, nlat, nlev, 1 /)),       &                            
                   'create_vtec', 'get_var TE') 

!... OP (midpoints)
call nc_check(nf90_inq_varid(ncid, 'OP', VarID), 'create_vtec', 'inq_varid OP')
call nc_check(nf90_get_var(ncid, VarID, values=OP,     &
                   start = (/ 1, 1, 1, last_time /),      &
                   count = (/ nlon, nlat, nlev, 1 /)),       &
                   'create_vtec', 'get_var OP')

! Construct vTEC given the parts

Z             = Z * 1.0e-2_r8            ! Convert Z (geopotential height) in cm to m
earth_radiusm = earth_radius * 1000.0_r8 ! Convert earth_radius in km to m
NE            = NE * 1.0e+6_r8           ! Convert NE in #/cm^3 to #/m^3 


! Convert to geometric height from geopotential height 
ZGG  = (earth_radiusm * Z) / (earth_radiusm - Z)

! Gravity at the top layer
GRAVITYtop(:,:) = gravity * (earth_radiusm / (earth_radiusm + ZGG(:,:,nilev))) ** 2

! Plasma Temperature 
Tplasma(:,:) = (TI(:,:,nlev-1) + TE(:,:,nlev-1)) / 2.0_r8

! Compute plasma scale height
Hplasma = (2.0_r8 * k_constant / omass ) * Tplasma / GRAVITYtop   

! NE is extrapolated to 10 more layers                         
nlev10  = nlev + 10

ZGG_extended(:,:,1:nlev) = ZGG
NEm_extended(:,:,1:nlev) = NE 

do j = nlev, nlev10 
   NEm_extended(:,:,j) = 2.0_r8 * NEm_extended(:,:,j-1) * exp(-1.0_r8)
   ZGG_extended(:,:,j) = ZGG_extended(:,:,j-1) + Hplasma(:,:) / 2.0_r8 
enddo 

! finally calculate vTEC - one gridcell at a time.

do j = 1, nlat
do k = 1, nlon
   delta_ZGG(1:(nlev10-1)) = ZGG_extended(k,j,2:nlev10)-ZGG_extended(k,j,1:(nlev10-1))  
   NE_middle(1:(nlev10-1)) = NEm_extended(k,j,2:nlev10)+NEm_extended(k,j,1:(nlev10-1))/2.0_r8 
   vTEC(k,j) = sum(NE_middle * delta_ZGG) * 1.0e-16_r8 ! Convert to TECU (1.0e+16 #/m^2) 
enddo
enddo  
                     
deallocate( NE, NEm_extended, ZGG, ZGG_extended, Z)
deallocate( TI, TE, OP )
deallocate( GRAVITYtop, Tplasma, Hplasma )
deallocate( delta_ZGG, NE_middle )

end subroutine create_vtec



subroutine get_val(val, x, lon_index, lat_index, height, obs_kind, istatus)
!-------------------------------------------------------------------------------
!
real(r8),           intent(out) :: val
real(r8),           intent(in)  :: x(:)
integer,            intent(in)  :: lon_index, lat_index
real(r8),           intent(in)  :: height
integer,            intent(in)  :: obs_kind
integer,            intent(out) :: istatus

integer                         :: ivar
integer                         :: k, lev_top, lev_bottom
real(r8)                        :: zgrid, delta_z, zgrid_top, zgrid_bottom
real(r8)                        :: val_top, val_bottom, frac_lev

! No errors to start with
istatus = 0

! To find a layer height: what's the unit of height [m]
! pressure level ln(p0/p) -- interface [-7.0 7.0] and midlevel [-6.75 7.25]
! Ne and ZG are defined at interface
! T, U, V, O, O2 & OP are defined at midlevel
! T, U, V  at top midlevel pressure level are missing values in TIEGCM
! but filled in DART with the values at nlev -1

write(string1,*) 'TJH routine needs to be rewritten/checked'
call error_handler(E_ERR,'get_val',string1,source,revision,revdate)

val = MISSING_R8
ivar = -1

! FIXME There is some confusion about using the T-minus-1 variables
! in this construct. Both TN and TN_NM have the same dart_kind,
! so we use the first one ... but it is not guaranteed that TN
! must preceed TN_NM, for example. 

VARLOOP : do k = 1,nfields
   if (progvar(ivar)%dart_kind == obs_kind) then
      ivar = k
      exit VARLOOP
   endif
enddo VARLOOP

! TJH FIXME ... what do we do if no variable with the KIND of interest

if (ivar < 0 ) then ! Failed  
   istatus = 1
   val = 0.0_r8
   return
endif

! TJH FIXME ... convert ZG to m when we stick it in the dart state vector?
!               would remove all the divide by 100's

if (obs_kind == KIND_ELECTRON_DENSITY) then

   zgrid_bottom = x(get_index(ivarZG, lon_index, lat_index,    1)) / 100.0_r8  ![m] = /100 [cm]
   zgrid_top    = x(get_index(ivarZG, lon_index, lat_index, nlev)) / 100.0_r8
   if ((zgrid_bottom > height) .or. (zgrid_top < height)) then
      istatus = 1 !obs height is above or below the model boundary
      val = 0.0
      return
   endif
 
   h_loop_interface : do k = 2, nilev
 
      zgrid = x(get_index(ivarZG, lon_index, lat_index, k)) / 100.0_r8 ![m] = /100 [cm]
 
      if (height   <= zgrid) then
         lev_top    = k
         lev_bottom = lev_top -1
         delta_z    = zgrid - x(get_index(ivarZG,lon_index,lat_index,lev_bottom)) / 100.0_r8
         frac_lev   = (zgrid - height)/delta_z
         exit h_loop_interface
      endif
 
   enddo h_loop_interface

else

   !mid_level 1      [m] = /100 [cm]
   zgrid_bottom = 0.50_r8 / 100.0_r8 * &
                  (x(get_index(ivarZG, lon_index, lat_index, 1)) + &
                   x(get_index(ivarZG, lon_index, lat_index, 2)))
 
   !mid_level nlev-1
   zgrid_top    = 0.50_r8 / 100.0_r8 * &
                  (x(get_index(ivarZG, lon_index, lat_index, nlev-1)) + &
                   x(get_index(ivarZG, lon_index, lat_index, nlev)))
 
   if ((zgrid_bottom > height) .or. (zgrid_top < height)) then
     istatus = 1 !obs height is above or below the model boundary
     val = 0.0
     return
   endif
 
   h_loop_midpoint:do k = 2, nilev-1
 
     zgrid = 0.50_r8 / 100.0_r8 * &               ! [m] = ZGtiegcm/100 [cm]
             (x(get_index(ivarZG, lon_index, lat_index, k  )) + &
              x(get_index(ivarZG, lon_index, lat_index, k+1)))
 
     if (height  <= zgrid) then
       lev_top    = k
       lev_bottom = lev_top -1
       delta_z    = zgrid -  0.50_r8 / 100.0_r8 * &
                    (x(get_index(ivarZG, lon_index, lat_index, lev_bottom  )) + &
                     x(get_index(ivarZG, lon_index, lat_index, lev_bottom+1)))
         frac_lev = (zgrid - height)/delta_z
       exit h_loop_midpoint
     endif
 
   enddo h_loop_midpoint

endif

if (obs_kind == KIND_PRESSURE) then
   val_top    = plevs(lev_top)     !pressure at midpoint [Pa]
   val_bottom = plevs(lev_bottom)  !pressure at midpoint [Pa]
   val = exp(frac_lev * log(val_bottom) + (1.0 - frac_lev) * log(val_top))
else
   ! This is the part confounded by using the same dart_kind for similar variables.
   val_top    = x(get_index(ivar, lon_index, lat_index, lev_top))
   val_bottom = x(get_index(ivar, lon_index, lat_index, lev_bottom))
   val =     frac_lev *     val_bottom  + (1.0 - frac_lev) *     val_top
endif

end subroutine get_val



function convert_lnpressure_to_height(ivar, lonindex, latindex, levindex)
!-----------------------------------------------------
! TIEGCM's 'natural' vertical coordinate is pressure, DART needs it in height
! If it is geopotential height already ... this is easy.

integer, intent(in) :: ivar
integer, intent(in) :: lonindex
integer, intent(in) :: latindex
integer, intent(in) :: levindex
real(r8)            :: convert_lnpressure_to_height

integer  :: myindex
real(r8) :: height
real(r8) :: lnpressure

if (     trim(progvar(ivar)%verticalvar) ==  'lev') then ! "midpoint  levels"
    lnpressure =  levs(levindex)
elseif ( trim(progvar(ivar)%verticalvar) == 'ilev') then ! "interface levels"
    lnpressure = ilevs(levindex)
else
   write(string1,*)'unknown vertical coordinate system <', &
                                  trim(progvar(ivar)%verticalvar),'>'
   write(string2,*)'on variable ',trim(progvar(ivar)%varname)
   call error_handler(E_ERR,'convert_lnpressure_to_height', string1, &
                      source, revision, revdate, text2=string2)
endif

! Need the ensemble mean value of ZG for this location
! ZG exists on ilev coordinates ... "interface levels"

myindex = get_index(ivarZG, lonindex, latindex, levindex)

if ((progvar(ivar)%dart_kind == KIND_ELECTRON_DENSITY)     .or.  &
    (progvar(ivar)%dart_kind == KIND_GEOPOTENTIAL_HEIGHT)) then
    !NE defined at interface levels
    !ZG defined at interface levels

     height = ens_mean(myindex)

     ! TOMOKO ... FIXME ... you have ensemble mean height and the lnpressure at this point
     ! I just need to set a return value to satisfy the compiler
     convert_lnpressure_to_height = height * lnpressure  ! BOGUS algorithm

else
     !TN UN VN O1 OP defined at midpoints
     !TIEGCM: top midpoint slot contains missing values for TN UN VN OP
!    if (lev_index+2 > nlev) then
!       height = ens_mean(get_index(lon_index, lat_index, &
!                  nlev,TYPE_local_ZG)) / 100.0_r8
!    else
!       height = 0.50_r8 / 100.0_r8 * &
!              (ens_mean(get_index(lon_index, lat_index,lev_index,TYPE_local_ZG)) +              &
!               ens_mean(get_index(lon_index_1, lat_index+1, lev_index+2,TYPE_local_ZG)))
!    endif

     convert_lnpressure_to_height = height * lnpressure  ! BOGUS algorithm

endif

end function convert_lnpressure_to_height



function get_index(ivar, ilon, ilat, ilev)
!-------------------------------------------------------------------------------
! returns the statevector index for the associated variable at 
integer,           intent(in) :: ivar
integer, optional, intent(in) :: ilon, ilat, ilev
integer                       :: get_index

if     (progvar(ivar)%rank == 0) then ! scalars
    get_index = progvar(ivar)%index1

elseif (progvar(ivar)%rank == 3) then ! Fortran shape is lon.lat.time

    if (present(ilon) .and. present(ilat)) then
       
       get_index = progvar(ivar)%index1 + ilon + (ilat-1)*nlon + 1

       write(*,*)'get_index ivar, ilon, ilat, indx ',ivar, ilon, ilat, get_index

    else
       write(string1,*) trim(progvar(ivar)%varname)//' has both lon and lat dimensions.'
       write(string2,*) 'lon_index present is ',present(ilon)
       write(string3,*) 'lat_index present is ',present(ilat)
       call error_handler(E_ERR, 'get_index', string1, &
          source, revision, revdate, text2=string2, text3=string3)
    endif

elseif (progvar(ivar)%rank == 4) then ! Fortran shape is lon.lat.lev.time

    if (present(ilon) .and. present(ilat) .and. present(ilev)) then
       get_index = progvar(ivar)%index1
    else
       write(string1,*) trim(progvar(ivar)%varname)//' has shape(lon,lat,lev).'
       write(string2,*)'not requesting all three indices.'
       call error_handler(E_ERR, 'get_index', string1, &
          source, revision, revdate, text2=string2)
    endif

else
   write(string1,*) trim(progvar(ivar)%varname)//' has unsupported shape.'
   call error_handler(E_ERR, 'get_index', string1, source, revision, revdate )
endif

end function get_index



!===============================================================================
! All these routines are overloaded to 
! prog_var_to_vector() and
! vector_to_prog_var()
!===============================================================================

subroutine var1d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:)
real(r8), intent(out) ::   x(:)

integer :: icount, idim1

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:1D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim1 = 1,size(var,1)
   x(icount) = var(idim1)
   icount = icount + 1
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:1D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var1d_to_vector



subroutine vector_to_var1d(x, ifield, var)
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ifield
real(r8), intent(out) :: var(:)

integer :: icount, idim1

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'vector_to_prog_var:1D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim1 = 1,size(var,1)
   var(idim1) = x(icount)
   icount = icount + 1
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'vector_to_prog_var:1D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine vector_to_var1d

!-------------------------------------------------------------------------------

subroutine var2d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:,:)
real(r8), intent(out) :: x(:)

integer :: icount,idim1,idim2

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:2D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   x(icount) = var(idim1,idim2)
   icount = icount + 1 
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:2D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var2d_to_vector



subroutine vector_to_var2D(x, ifield, var)
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ifield
real(r8), intent(out) :: var(:,:)

integer :: icount,idim1,idim2

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'vector_to_prog_var:2D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   var(idim1,idim2) = x(icount)
   icount = icount + 1 
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'vector_to_prog_var:2D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine vector_to_var2D

!-------------------------------------------------------------------------------

subroutine var3d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:,:,:)
real(r8), intent(out) :: x(:)

integer :: icount,idim1,idim2,idim3

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:3D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim3 = 1,size(var,3)
do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   x(icount) = var(idim1,idim2,idim3)
   icount = icount + 1 
enddo
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:3D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var3d_to_vector


subroutine vector_to_var3d(x, ifield, var)
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ifield
real(r8), intent(out) :: var(:,:,:)

integer :: icount,idim1,idim2,idim3

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'vector_to_prog_var:3D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim3 = 1,size(var,3)
do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   var(idim1,idim2,idim3) = x(icount)
   icount = icount + 1 
enddo
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'vector_to_prog_var:3D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine vector_to_var3d

!-------------------------------------------------------------------------------

subroutine var4d_to_vector(ifield, var, x)
integer,  intent(in)  :: ifield
real(r8), intent(in)  :: var(:,:,:,:)
real(r8), intent(out) :: x(:)

integer :: icount,idim1,idim2,idim3,idim4

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'prog_var_to_vector:4D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim4 = 1,size(var,4)
do idim3 = 1,size(var,3)
do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   x(icount) = var(idim1,idim2,idim3,idim4)
   icount = icount + 1 
enddo
enddo
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'prog_var_to_vector:4D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine var4d_to_vector


subroutine vector_to_var4d(x, ifield, var)
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: ifield
real(r8), intent(out) :: var(:,:,:,:)

integer :: icount,idim1,idim2,idim3,idim4

if (size(var) /= progvar(ifield)%varsize) then
   write(string1, *) trim(progvar(ifield)%varname), &
              'is supposed to have size ',progvar(ifield)%varsize
   write(string2, *) 'it has size ',size(var)
   call error_handler(E_ERR, 'vector_to_prog_var:4D', string1, &
              source, revision, revdate, text2=string2)
endif

icount = progvar(ifield)%index1

do idim4 = 1,size(var,4)
do idim3 = 1,size(var,3)
do idim2 = 1,size(var,2)
do idim1 = 1,size(var,1)
   var(idim1,idim2,idim3,idim4) = x(icount)
   icount = icount + 1 
enddo
enddo
enddo
enddo

if ( (icount-1) /= progvar(ifield)%indexN ) then
   write(string1, *) trim(progvar(ifield)%varname), 'supposed to end at ', &
                          progvar(ifield)%indexN
   write(string2, *) 'it ended at ',icount-1
   call error_handler(E_ERR, 'vector_to_prog_var:4D', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine vector_to_var4d



function FindVar_by_type(ikind)
! Finds the first variable of the appropriate DART KIND

integer, intent(in) :: ikind
integer             :: FindVar_by_type

integer :: ivar

VARLOOP : do ivar = 1,nfields
   if (progvar(ivar)%dart_kind == ikind) then
      FindVar_by_type = ivar
      return
   endif
enddo VARLOOP

write(string1, *) 'unable to find a variable with a DART kind of ',ikind
call error_handler(E_ERR, 'FindVar_by_type', string1, source, revision, revdate )

end function FindVar_by_type



subroutine test_module()
! Take a stroll through the state vector and make sure I have all the right
! attributes for every element in the state vector

integer :: i, ivar, var_kind
type(location_type) :: location


do i = 1,model_size

   ivar = Find_Variable_by_index(i,'test_module')

   call get_state_meta_data(i, location, var_kind)

enddo


end subroutine test_module




function Find_Variable_by_index(myindx,msgstring)
! Given an index into the DART state vector, return the index of metadata
! variable 'progvar' responsible for this portion of the state vector
integer,          intent(in) :: myindx
character(len=*), intent(in) :: msgstring
integer                      :: Find_Variable_by_index

integer :: ivar

Find_Variable_by_index = -1

FindIndex : do ivar = 1,nfields
   if ((myindx >= progvar(ivar)%index1)  .and. &
       (myindx <= progvar(ivar)%indexN)) then
      Find_Variable_by_index = ivar
      exit FindIndex
   endif
enddo FindIndex

if (Find_Variable_by_index < 0) then
   write(string1,*)'index ',myindx,' is out of range of all variables.'
   write(string2,*)'model size is ',model_size
   call error_handler(E_ERR, 'Find_Variable_by_index'//trim(msgstring), string1, &
                      source, revision, revdate, text2=string2 )
endif

end function Find_Variable_by_index

!===============================================================================
! End of model_mod
!===============================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
