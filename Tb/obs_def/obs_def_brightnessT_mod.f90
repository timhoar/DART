! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Banner for checking maximum lengths (i.e. 32)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! BEGIN DART PREPROCESS KIND LIST
! AMSRE_BRIGHTNESS_T,             KIND_BRIGHTNESS_TEMPERATURE
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
!  case(AMSRE_BRIGHTNESS_T)
!     call get_brightness_temperature(state_time, ens_index, location, obs_key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!  case(AMSRE_BRIGHTNESS_T)
!     call read_amsre_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!  case(AMSRE_BRIGHTNESS_T)
!     call write_amsre_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!  case(AMSRE_BRIGHTNESS_T)
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
                             set_location, VERTISUNDEF, LocationDims
use time_manager_mod, only : time_type, get_date, set_date, print_date, print_time, &
                             get_time, set_time, operator(-), operator(/=)
use    utilities_mod, only : register_module, E_ERR, E_MSG, error_handler,     &
                             check_namelist_read, find_namelist_in_file,       &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             nc_check, file_exist, is_longitude_between,       &
                             ascii_file_format
use        model_mod, only : get_ncols_in_gridcell, &
                             get_colids_in_gridcell, &
                             get_clm_restart_filename, &
                             get_gridsize, get_grid_arrays

use typesizes
use netcdf

use radiative_transfer_mod, only : ss_model

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
integer            :: nlon, nlat, nlevgrnd, ncolumns, nlevtot, nlevsno, ens_size
type(time_type)    :: model_time
integer, PARAMETER :: LAKEUNIT = 3

character(len=129), allocatable, dimension(:) :: fname   ! instances of CLM restart files
integer,            allocatable, dimension(:) :: ncid    ! instances of CLM restart files
real(r8),           allocatable, dimension(:) :: lon, lat, levgrnd, area
integer,            allocatable, dimension(:) :: cols1d_ityplun
real(r8),           allocatable, dimension(:) :: cols1d_wtxy

real(r8), PARAMETER :: RAD2KM = 40030.0_r8/(2.0_r8 * PI) ! (mean radius of earth ~6371km)

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
real(r4),         parameter :: AMSRE_inc_angle = 55.0_r4 ! incidence angle (degrees)

integer, SAVE :: MAXamsrekey = 24*366  ! one year of hourly data - to start
integer, SAVE ::    amsrekey = 0       ! useful length of metadata arrays

!----------------------------------------------------------------------
! Properties required for a snow column
!----------------------------------------------------------------------

type snowprops
   private
   integer  :: nlayers         ! aux_ins(1)
   real(r4) :: t_grnd          ! aux_ins(2) ground temperature [K]
   real(r4) :: soilsat         ! aux_ins(3) soil saturation [fraction]
   real(r4) :: soilporos       ! aux_ins(4) soil porosity [fraction]
   real(r4) :: propconst       ! aux_ins(5) proportionality between grain size & correlation length.
   integer  :: nprops          ! [thickness, density, diameter, liqwater, temperature]
   real(r4), pointer, dimension(:) :: thickness      !  LAYER THICKNESS [M]
   real(r4), pointer, dimension(:) :: density        !  LAYER DENSITY [KG/M3]
   real(r4), pointer, dimension(:) :: grain_diameter !  LAYER GRAIN DIAMETER [M]
   real(r4), pointer, dimension(:) :: liquid_water   !  LAYER LIQUID WATER CONTENT [FRAC]
   real(r4), pointer, dimension(:) :: temperature    !  LAYER TEMPERATURE [K]
end type snowprops

type(snowprops) :: snowcolumn

!----------------------------------------------------------------------
! namelist items
!----------------------------------------------------------------------

character(len=256) :: casename = 'clm_dart'
logical            :: debug = .false.

namelist /obs_def_brightnessT_nml/ casename, debug

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


subroutine get_brightness_temperature(state_time, ens_index, location, key, obs_val, istatus)

! This is THE forward observation operator. Given the state and a location, return the value
! The parts of the state required for this forward operator are not required to
! be part of the DART state vector. They are currently directly harvested from the CLM
! restart file. As such, the posteriors are not informative.

type(time_type),     intent(in)  :: state_time ! valid time of DART state
integer,             intent(in)  :: ens_index  ! Ensemble member number
type(location_type), intent(in)  :: location   ! target location (i.e. obs)
integer,             intent(in)  :: key        ! pointer for obs custom metadata
real(r8),            intent(out) :: obs_val    ! model estimate of observation value
integer,             intent(out) :: istatus    ! status of the calculation

integer, parameter :: N_FREQ = 1  ! observations come in one frequency at a time
integer, parameter :: N_POL  = 2  ! code automatically computes both polarizations

! variables required by ss_snow() routine
real(r4), allocatable, dimension(:,:) :: y ! 2D array of snow properties
real(r4) :: aux_ins(5)     ! [nsnowlyrs, ground_T, soilsat, poros, proportionality]
integer  :: ctrl(4)        ! [n_lyrs, n_aux_ins, n_snow_ins, n_freq]
real(r4) :: freq( N_FREQ)  ! frequencies at which calculations are to be done
real(r4) :: tetad(N_FREQ)  ! incidence angle of satellite
real(r4) :: tb_ubc(N_POL,N_FREQ) ! upper boundary condition brightness temperature
real(r4) :: tb_out(N_POL,N_FREQ) ! calculated brightness temperature - output

! support variables 
character(len=256)                  :: filename
integer,  allocatable, dimension(:) :: columns_to_get
real(r4), allocatable, dimension(:) :: tb
real(r8), allocatable, dimension(:) :: weights
real(r8), dimension(LocationDims)   :: loc
real(r8)  :: loc_lon, loc_lat
character :: pol                    ! observation polarization
integer   :: ilonmin(1), ilatmin(1) ! need to be array-valued for minloc intrinsic
integer   :: ilon, ilat, icol, ncols, lccode

istatus  = 1
obs_val  = MISSING_R8

tetad(:) = AMSRE_inc_angle
freq(:)  = observation_metadata(key)%frequency
pol      = observation_metadata(key)%polarization
lccode   = observation_metadata(key)%landcovercode

! FIXME ... determine lon/lat indices
! Poor method ... will not work for single column case nor
! for irregular/unstructured grids might not work well over poles ...
loc      = get_location(location)
loc_lon  = loc(1)   ! longitude of observation (in degrees)
loc_lat  = loc(2)   ! latitude  of observation (in degrees)
ilonmin  = minloc( abs(loc_lon-LON) )
ilatmin  = minloc( abs(loc_lat-LAT) )
ilon     = ilonmin(1)
ilat     = ilatmin(1)

ncols    = get_ncols_in_gridcell(ilon,ilat)

! Early return if there are no CLM columns at this location.
! Forward operator returns 'failed', but life goes on.
if ((ncols == 0) .and. do_output() ) then
   write(string1, *) 'gridcell ilon/ilat (',ilon,ilat,') has no CLM columns.'
   write(string2, '(''lon,lat ('',f12.6,'','',f12.6,'')'')')
   call error_handler(E_MSG,'get_brightness_temperature',string1,text2=string2)
   return
endif

if (debug .and. do_output()) write(*,*)
if (debug .and. do_output()) write(*,*)
if (debug .and. do_output()) write(*,*)'TJH debug ... computing gridcell ',ilon,ilat

allocate( columns_to_get(ncols), tb(ncols), weights(ncols) )
columns_to_get(:) = -1
tb(:)             = 0.0_r4
weights(:)        = 0.0_r8
call get_colids_in_gridcell(ilon, ilat, columns_to_get)

! FIXME Presently skipping gridcells with lakes.
! get_column_snow() must also modified to use bulk snow formulation for lakes.
if ( any(cols1d_ityplun(columns_to_get) == LAKEUNIT))  then
   deallocate(columns_to_get, tb, weights)
   return
endif

! need to know which restart file to use to harvest information
call build_clm_instance_filename(ens_index, state_time, filename)

! Loop over all columns in the gridcell that has the right location.

SNOWCOLS : do icol = 1,ncols

   weights(icol) = cols1d_wtxy(columns_to_get(icol)) ! relative weight of column
   call get_column_snow(filename, columns_to_get(icol)) ! allocates snowcolumn

   if ( debug .and. do_output() ) then
      if (snowcolumn%nlayers < 1) then
         write(string1, *) 'column (',columns_to_get(icol),') has no snow'
         call error_handler(E_MSG,'get_brightness_temperature',string1)
      else
         write(*,*)'nprops   ',snowcolumn%nprops
         write(*,*)'nlayers  ',snowcolumn%nlayers
         write(*,*)'t_grnd   ',snowcolumn%t_grnd
         write(*,*)'soilsat  ',snowcolumn%soilsat
         write(*,*)'soilpor  ',snowcolumn%soilporos
         write(*,*)'proconst ',snowcolumn%propconst
         write(*,*)'thickness',snowcolumn%thickness
         write(*,*)'density  ',snowcolumn%density
         write(*,*)'diameter ',snowcolumn%grain_diameter
         write(*,*)'liqwater ',snowcolumn%liquid_water
         write(*,*)'temp     ',snowcolumn%temperature
      endif
   endif

   if ( snowcolumn%nlayers == 0 ) then
      ! If there is no snow, the ss_model will calculate the brightness
      ! temperature of the bare soil. To indicate this, aux_ins(1) must
      ! be 0 and ctrl(1) must be 1
      ctrl(1) = 1
   else
      ctrl(1) = snowcolumn%nlayers
   endif
   ctrl(2) = 0              ! not used as far as I can tell
   ctrl(3) = snowcolumn%nprops
   ctrl(4) = N_FREQ

   aux_ins(1) = real(snowcolumn%nlayers,r4)
   aux_ins(2) = snowcolumn%t_grnd
   aux_ins(3) = snowcolumn%soilsat
   aux_ins(4) = snowcolumn%soilporos
   aux_ins(5) = 0.5_r4                ! FIXME - hardwired

   ! TJH must convert liquid water to a fraction [0,1]
   ! TJH 1cm3 == 1g so 100cm * 100cm = 10,000g
   ! fully saturated then is 10kg/m^2

   allocate( y(ctrl(1), snowcolumn%nprops) )
   if ( aux_ins(1) > 0 ) then
      y(:,1)  = snowcolumn%thickness
      y(:,2)  = snowcolumn%density
      y(:,3)  = snowcolumn%grain_diameter / 1000000.0_r4 ! need meters (from microns)
      y(:,4)  = snowcolumn%liquid_water / 10.0_r4        ! convert to fraction
      y(:,5)  = snowcolumn%temperature
   else ! dummy values for bare ground
      y(:,1)  = 0.0_r4
      y(:,2)  = 0.0_r4
      y(:,3)  = 0.0_r4
      y(:,4)  = 0.0_r4
      y(:,5)  = 0.0_r4
   endif

   ! FIXME Ally ... if you have a better way to specify/determine,
   ! tb_ubc do it here. Call your atmospheric model to calculate it.
   tb_ubc(:,N_FREQ) = (/ 2.7_r4, 2.7_r4 /)

   ! If the landcovercode (lccode) indicates that you want to use
   ! a different radiative transfer model ... implement it here.
   ! this will involve changing the following 'if' statement.

   if (lccode > 0 ) then
      ! the tb_out array contains the calculated brightness temperature outputs
      ! at each polarization (rows) and frequency (columns).
      call ss_model(ctrl, freq, tetad, y, tb_ubc, aux_ins, tb_out)
   else
      ! call to alternative radiative transfer model goes here.
   endif

   write(*,*)'column ', columns_to_get(icol),' tb_out is ',tb_out
   write(*,*)

   if (pol == 'H') then
      tb(icol) = tb_out(1,1)   ! second dimension is only 1 frequency
   else
      tb(icol) = tb_out(2,1)   ! second dimension is only 1 frequency
   endif

   deallocate( y )
   call destroy_column_snow()

enddo SNOWCOLS

! FIXME ... account for heterogeneity somehow ...
! must aggregate all columns in the gridcell
! area-weight the average
obs_val = sum(tb * weights) / sum(weights)

if (debug .and. do_output()) then
   write(*,*)'tb      for all columns is ',tb
   write(*,*)'weights for all columns is ',weights
   write(*,*)'(weighted) obs value    is ',obs_val
endif

deallocate( columns_to_get, tb, weights )

istatus = 0

end subroutine get_brightness_temperature


!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of internal routines for obs_def_brightnessT_mod
!----------------------------------------------------------------------
!----------------------------------------------------------------------


subroutine initialize_module

! Called once to set values and allocate space.

integer :: iunit, io, myncid, varid
integer :: yyyymmdd, sssss
integer :: year, month, day, hour, minute, second, leftover

character(len=258) :: filename = 'no_clm_restart_file'

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
amsrekey = 0
observation_metadata(:) = missing_metadata

call get_clm_restart_filename( filename )

call nc_check(nf90_open(trim(filename), nf90_nowrite, myncid), &
          'obs_def_brightnessT_mod.initialize_routine','open '//trim(filename))

call GetDimensions(myncid,filename) ! sets nlon, nlat, ncolumns, etc.

allocate( lon(nlon), lat(nlat), levgrnd(nlevgrnd) )
call get_grid_arrays(lon, lat, levgrnd)

! Read in the column landunit type. Knowing if the column is a lake is useful.
allocate( cols1d_ityplun(ncolumns) )
call nc_check(nf90_inq_varid(myncid,'cols1d_ityplun', varid), &
              'obs_def_brightnessT_mod.initialize_routine', 'inq_varid cols1d_ityplun')
call nc_check(nf90_get_var(  myncid, varid, cols1d_ityplun ), &
              'obs_def_brightnessT_mod.initialize_routine', 'get_var   cols1d_ityplun')

! Read in the column weight relative to corresponding gridcell
allocate( cols1d_wtxy(ncolumns) )
call nc_check(nf90_inq_varid(myncid,'cols1d_wtxy', varid), &
              'obs_def_brightnessT_mod.initialize_routine', 'inq_varid cols1d_wtxy')
call nc_check(nf90_get_var(  myncid, varid, cols1d_wtxy ), &
              'obs_def_brightnessT_mod.initialize_routine', 'get_var   cols1d_wtxy')

! Read in the current model time
call nc_check(nf90_inq_varid(myncid, 'timemgr_rst_curr_ymd', varid), &
              'obs_def_brightnessT_mod.initialize_routine','inq_varid timemgr_rst_curr_ymd '//trim(filename))
call nc_check(nf90_get_var(  myncid, varid, yyyymmdd), &
              'obs_def_brightnessT_mod.initialize_routine','get_var yyyymmdd '//trim(filename))

call nc_check(nf90_inq_varid(myncid, 'timemgr_rst_curr_tod',  varid), &
              'obs_def_brightnessT_mod.initialize_routine','inq_varid timemgr_rst_curr_tod '//trim(filename))
call nc_check(nf90_get_var(  myncid, varid,    sssss), &
              'obs_def_brightnessT_mod.initialize_routine','get_var   sssss '//trim(filename))

year     = yyyymmdd/10000
leftover = yyyymmdd - year*10000
month    = leftover/100
day      = leftover - month*100

hour     = sssss/3600
leftover = sssss - hour*3600
minute   = leftover/60
second   = leftover - minute*60

model_time = set_date(year, month, day, hour, minute, second)
call get_time(model_time, second, day)

if (debug .and. do_output()) then
   write(*,*)'obs_def_brightnessT_mod.initialize_routine:nlon,nlat,ncolumns is ',nlon,nlat,ncolumns
   write(*,*)'obs_def_brightnessT_mod.initialize_routine: model time is ',yyyymmdd, sssss
   call print_date(model_time,'obs_def_brightnessT_mod.initialize_routine:model date is')
   call print_time(model_time,'obs_def_brightnessT_mod.initialize_routine:model time is')
endif

if (debug .and. do_output()) then
   call test_ss_model()
endif

end subroutine initialize_module


!======================================================================


subroutine GetDimensions(myncid, fname)

! Harvest information from the CLM restart file and model_mod.
! This routine sets the following module variables:
! nlon, nlat, nlevgrnd, ncolumn, nlevtot, nlevsno, nlevgrnd

integer,          intent(in) :: myncid
character(len=*), intent(in) :: fname

integer :: dimid

! Get the number of columns in the restart file
call nc_check(nf90_inq_dimid(myncid, 'column', dimid), &
              'obs_def_brightnessT.GetDimensions','inq_dimid column '//trim(fname))
call nc_check(nf90_inquire_dimension(myncid, dimid, len=ncolumns), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension column '//trim(fname))

! Get the number of total levels in the restart file
call nc_check(nf90_inq_dimid(myncid, 'levtot', dimid), &
              'obs_def_brightnessT.GetDimensions','inq_dimid levtot '//trim(fname))
call nc_check(nf90_inquire_dimension(myncid, dimid, len=nlevtot), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension levtot '//trim(fname))

! Get the number of snow levels in the restart file
call nc_check(nf90_inq_dimid(myncid, 'levsno', dimid), &
              'obs_def_brightnessT.GetDimensions','inq_dimid levsno '//trim(fname))
call nc_check(nf90_inquire_dimension(myncid, dimid, len=nlevsno), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension levsno '//trim(fname))

! Get the number of ground levels in the restart file
call nc_check(nf90_inq_dimid(myncid, 'levgrnd', dimid), &
              'obs_def_brightnessT.GetDimensions','inq_dimid levgrnd '//trim(fname))
call nc_check(nf90_inquire_dimension(myncid, dimid, len=nlevgrnd), &
              'obs_def_brightnessT.GetDimensions','inquire_dimension levgrnd '//trim(fname))

call get_gridsize(nlon, nlat, nlevgrnd)

end subroutine GetDimensions


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


subroutine build_clm_instance_filename(instance, state_time, filename)
! If the instance is 1, it could be a perfect model scenario
! or it could be the first instance of many. CLM has a different
! naming scheme for these.
!
! the casename is part of the module namelist

integer,          intent(in)  :: instance
type(time_type),  intent(in)  :: state_time
character(len=*), intent(out) :: filename

integer :: year, month, day, hour, minute, second

100 format (A,'.clm2_',I4.4,'.r.',I4.4,'-',I2.2,'-',I2.2,'-',I5.5,'.nc')
110 format (A,'.clm2'      ,'.r.',I4.4,'-',I2.2,'-',I2.2,'-',I5.5,'.nc')

call get_date(state_time, year, month, day, hour, minute, second)
second = second + minute*60 + hour*3600

write(filename,110) trim(casename),year,month,day,second

! Check if in a perfect model scenario
if( file_exist(filename) ) then
   ens_size = 1
   if (debug .and. do_output()) then
      write(string1,*)'Running in a perfect model configuration with ',trim(filename)
      call error_handler(E_MSG, 'obs_def_brightnessT:build_clm_instance_filename', string1)
   endif
   return
endif

! 'normal' situation
write(filename,100) trim(casename),instance,year,month,day,second
if( file_exist(filename) ) then
   return
else
   write(string1,*)'Unable to create viable CLM restart filename:'
   call error_handler(E_ERR, 'obs_def_brightnessT:build_clm_instance_filename', &
        string1, text2=trim(filename))
endif

end subroutine build_clm_instance_filename


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



subroutine get_column_snow(filename, snow_column )
! Read all the variables needed for the radiative transfer model as applied
! to a single CLM column.
!
! The treatment of snow-related variables is complicated.
! The SNLSNO variable defines the number of snow layers with valid values.
! HOWEVER, if the snow depth is < 0.01 m, the snow is not represented by a layer,
! so the SNLSNO(i) is zero even though there is a trace of snow.
! Even a trace amount of snow results in some sort of snow cover fraction.
!
! Lakes are treated differently.
! The SNLSNO(i) is always zero, even though there is snow.
! The snow over lakes is wholly contained in the bulk formulation variables
! as opposed to the snow layer variables.


! float WATSAT(levgrnd, lat, lon) ;
!       WATSAT:long_name = "saturated soil water content (porosity)" ;
!       WATSAT:units = "mm3/mm3" ;
!       WATSAT:_FillValue = 1.e+36f ;
!       WATSAT:missing_value = 1.e+36f ;

character(len=*), intent(in)  :: filename
integer,          intent(in)  :: snow_column

real(r8) :: t_grnd(1) ! ground temperature
integer  :: snlsno(1) ! number of snow layers
integer  :: ityplun   ! type of column unit (lake, etc.)

real(r8), allocatable, dimension(:) :: h2osoi_liq, h2osoi_ice, t_soisno
real(r8), allocatable, dimension(:) :: dzsno, zsno, zisno, snw_rds

integer               :: myncid, varid, ilayer, nlayers, ij
integer, dimension(2) :: ncstart, nccount

! Set some return values
snowcolumn%nprops    = 0
snowcolumn%nlayers   = 0
snowcolumn%t_grnd    = 0.0_r4
snowcolumn%soilsat   = 0.0_r4
snowcolumn%soilporos = 0.0_r4
snowcolumn%propconst = 0.0_r4

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, myncid), &
              'get_column_snow','open '//trim(filename))

! Get the (scalar) number of active snow layers for this column.
call nc_check(nf90_inq_varid(myncid,'SNLSNO', varid), 'get_column_snow', 'inq_varid SNLSNO')
call nc_check(nf90_get_var(  myncid, varid, snlsno, start=(/ snow_column /), count=(/ 1 /)), &
              'get_column_snow', 'get_var SNLSNO')

nlayers = abs(snlsno(1))

! Get the ground temperature for this column.
! double T_GRND(column); long_name = "ground temperature" ; units = "K" ;
call nc_check(nf90_inq_varid(myncid,'T_GRND', varid), 'get_column_snow', 'inq_varid T_GRND')
call nc_check(nf90_get_var(  myncid, varid, t_grnd, start=(/ snow_column /), count=(/ 1 /)), &
              'get_column_snow', 'get_var T_GRND')

ityplun = cols1d_ityplun(snow_column)   ! are we a lake

! FIXME ... lake columns use a bulk formula for snow
if (ityplun == LAKEUNIT ) return

! double H2OSOI_LIQ(column, levtot); long_name = "liquid water" ; units = "kg/m2" ;
! double H2OSOI_ICE(column, levtot); long_name = "ice lens"     ; units = "kg/m2" ;
! double T_SOISNO(  column, levtot); long_name = "soil-snow temperature" ; units = "K" ;

allocate(h2osoi_liq(nlevtot), h2osoi_ice(nlevtot), t_soisno(nlevtot))
ncstart = (/ 1, snow_column /)
nccount = (/ nlevtot,   1   /)

call nc_check(nf90_inq_varid(myncid,'T_SOISNO', varid), 'get_column_snow', 'inq_varid T_SOISNO')
call nc_check(nf90_get_var(  myncid, varid, t_soisno, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var T_SOISNO')

call nc_check(nf90_inq_varid(myncid,'H2OSOI_LIQ', varid), 'get_column_snow', 'inq_varid H2OSOI_LIQ')
call nc_check(nf90_get_var(  myncid, varid, h2osoi_liq, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var H2OSOI_LIQ')

call nc_check(nf90_inq_varid(myncid,'H2OSOI_ICE', varid), 'get_column_snow', 'inq_varid H2OSOI_ICE')
call nc_check(nf90_get_var(  myncid, varid, h2osoi_ice, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var H2OSOI_ICE')

! double   DZSNO(column, levsno); long_name = "snow layer thickness"        ; units = "m" ;
! double    ZSNO(column, levsno); long_name = "snow layer depth"            ; units = "m" ;
! double   ZISNO(column, levsno); long_name = "snow interface depth"        ; units = "m" ;
! double snw_rds(column, levsno); long_name = "snow layer effective radius" ; units = "um" ;

allocate(dzsno(nlevsno), zsno(nlevsno), zisno(nlevsno), snw_rds(nlevsno))
ncstart = (/ 1, snow_column /)
nccount = (/ nlevsno,   1   /)

call nc_check(nf90_inq_varid(myncid,'DZSNO', varid), 'get_column_snow', 'inq_varid DZSNO')
call nc_check(nf90_get_var(  myncid, varid, dzsno, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var DZSNO')

call nc_check(nf90_inq_varid(myncid,'ZSNO', varid), 'get_column_snow', 'inq_varid ZSNO')
call nc_check(nf90_get_var(  myncid, varid, zsno, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var ZSNO')

call nc_check(nf90_inq_varid(myncid,'ZISNO', varid), 'get_column_snow', 'inq_varid ZISNO')
call nc_check(nf90_get_var(  myncid, varid, zisno, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var ZISNO')

call nc_check(nf90_inq_varid(myncid,'snw_rds', varid), 'get_column_snow', 'inq_varid snw_rds')
call nc_check(nf90_get_var(  myncid, varid, snw_rds, start=ncstart, count=nccount), &
     'get_column_snow', 'get_var snw_rds')

! Print a summary so far
if (debug .and. do_output()) then
   write(*,*)'get_column_snow: raw CLM data for column ',snow_column
   write(*,*)'  # of snow layers, column ityp, ground temp :', snlsno, ityplun, t_grnd
   write(*,*)'  h2osoi_liq :', h2osoi_liq(1:nlevsno)
   write(*,*)'  h2osoi_ice :', h2osoi_ice(1:nlevsno)
   write(*,*)'  t_soisno   :',   t_soisno(1:nlevsno)
   write(*,*)'  dzsno      :',      dzsno(1:nlevsno)
   write(*,*)'  zsno       :',       zsno(1:nlevsno)
   write(*,*)'  zisno      :',      zisno(1:nlevsno)
   write(*,*)'  snw_rds    :',    snw_rds(1:nlevsno)
endif

allocate( snowcolumn%thickness(nlayers)     , &
          snowcolumn%density(nlayers)       , &
          snowcolumn%grain_diameter(nlayers), &
          snowcolumn%liquid_water(nlayers)  , &
          snowcolumn%temperature(nlayers)   )

! Fill the output array ... finally
! FIXME ... soilsat   is a hardwired value
! FIXME ... soiporos  is a hardwired value
! FIXME ... propconst is a hardwired value

snowcolumn%nprops    = 5
snowcolumn%nlayers   = nlayers
snowcolumn%t_grnd    = t_grnd(1) ! aux_ins(2) ground temperature [K]
snowcolumn%soilsat   = 0.3       ! aux_ins(3) soil saturation [fraction]
snowcolumn%soilporos = 0.4       ! aux_ins(4) soil porosity [fraction]
snowcolumn%propconst = 0.5       ! aux_ins(5) proportionality between grain size & correlation length.

ij = 0
do ilayer = (nlevsno-nlayers+1),nlevsno
   ij = ij + 1
   snowcolumn%thickness(ij)      = dzsno(ilayer)
   snowcolumn%density(ij)        = (h2osoi_liq(ilayer) + h2osoi_ice(ilayer)) / dzsno(ilayer)
   snowcolumn%grain_diameter(ij) = snw_rds(ilayer)
   snowcolumn%liquid_water(ij)   = h2osoi_liq(ilayer)
   snowcolumn%temperature(ij)    = t_soisno(ilayer)
   write(*,*)'   get_column_snow: filling layer ',ij,' with info from ilayer ',ilayer
enddo

deallocate(h2osoi_liq, h2osoi_ice, t_soisno)
deallocate(dzsno, zsno, zisno, snw_rds)

end subroutine get_column_snow



subroutine destroy_column_snow

if (associated(snowcolumn%thickness))      deallocate(snowcolumn%thickness)
if (associated(snowcolumn%density))        deallocate(snowcolumn%density)
if (associated(snowcolumn%grain_diameter)) deallocate(snowcolumn%grain_diameter)
if (associated(snowcolumn%liquid_water))   deallocate(snowcolumn%liquid_water)
if (associated(snowcolumn%temperature))    deallocate(snowcolumn%temperature)

snowcolumn%nprops    = 0
snowcolumn%nlayers   = 0
snowcolumn%t_grnd    = 0.0_r4
snowcolumn%soilsat   = 0.0_r4
snowcolumn%soilporos = 0.0_r4
snowcolumn%propconst = 0.0_r4

end subroutine destroy_column_snow



subroutine test_ss_model()

! test THE forward observation operator

integer, parameter :: N_FREQ = 1  ! observations come in one frequency at a time
integer, parameter :: N_POL  = 2  ! code automatically computes both polarizations
integer   :: nlayers              ! number of snow levels - 5 in this case
character :: pol                  ! observation polarization [H,V]

! variables required by ss_snow() routine
real(r4), allocatable, dimension(:,:) :: y ! 2D array
real(r4) :: aux_ins(5) ! properties: [nlyrs, ground_T, soilsat, poros, proportionality]
integer  :: ctrl(4)        ! N_LYRS, N_AUX_INS, N_SNOW_INS, N_FREQ
real(r4) :: freq( N_FREQ)  ! frequencies at which calculations are to be done
real(r4) :: tetad(N_FREQ)  ! incidence angle of satellite
real(r4) :: tb_ubc(N_POL,N_FREQ) ! UPPER BOUNDARY CONDITION BRIGHTNESS TEMPERATURE
real(r4) :: tb_out(N_POL,N_FREQ) ! brightness temperature

tetad(:) = 55.0   ! AMSR-E incidence angle
freq(:)  = 89.0   ! test frequency (GHz)
pol      = 'H'    ! test polarization
nlayers  = 5      ! 5 snow layers in test

allocate( y(nlayers,5) ) ! snow layers -x- 5 properties

tb_ubc(:,N_FREQ) = (/ 2.7, 2.7 /)  ! two polarizations

ctrl(1) = nlayers
ctrl(2) = 0         ! not used as far as I can tell
ctrl(3) = 5
ctrl(4) = N_FREQ

aux_ins(1) = real(nlayers,r4)
aux_ins(2) = 271.1123
aux_ins(3) = 0.3
aux_ins(4) = 0.4
aux_ins(5) = 0.5_r4

y(:,1) = (/   0.6213,   0.2071,   0.1033,   0.0497,   0.0200 /) ! snow thickness (meters)
y(:,2) = (/ 270.5943, 150.7856,  96.1940,  66.4903,  58.0377 /) ! density (kg/m3)
y(:,3) = (/  93.1651,  84.0811,  67.6715,  66.8529,  65.6391 /) / 1000000.0_r4  ! grain diameter (m)
y(:,4) = (/   0.0000,   0.0000,   0.0000,   0.0000,   0.0000 /) ! liquid water fraction
y(:,5) = (/ 266.1220, 256.7763, 247.9525, 240.4609, 235.8929 /) ! temperature (K)

! the tb_out array contains the calculated brightness temperature outputs
! at each polarization (rows) and frequency (columns).

call ss_model(ctrl, freq, tetad, y, tb_ubc, aux_ins, tb_out)

write(*,*)'TEST case: tb_out is ',tb_out

deallocate( y )

end subroutine test_ss_model

!======================================================================


end module obs_def_brightnessT_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
