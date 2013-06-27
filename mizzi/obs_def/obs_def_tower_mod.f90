! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS KIND LIST
!WATER_TABLE_DEPTH,              KIND_WATER_TABLE_DEPTH,        COMMON_CODE
!SOIL_TEMPERATURE,               KIND_SOIL_TEMPERATURE,         COMMON_CODE
!SOIL_MOISTURE,                  KIND_SOIL_MOISTURE,            COMMON_CODE
!LAYER_LIQUID_WATER,             KIND_LIQUID_WATER,             COMMON_CODE
!LAYER_ICE,                      KIND_ICE,                      COMMON_CODE
!SNOW_THICKNESS,                 KIND_SNOW_THICKNESS,           COMMON_CODE
!SNOW_WATER,                     KIND_SNOW_WATER,               COMMON_CODE
!MODIS_SNOWCOVER_FRAC,           KIND_SNOWCOVER_FRAC,           COMMON_CODE
!LEAF_CARBON,                    KIND_LEAF_CARBON,              COMMON_CODE
!LEAF_AREA_INDEX,                KIND_LEAF_AREA_INDEX,          COMMON_CODE
!TOWER_AIR_TEMPERATURE,          KIND_TEMPERATURE,              COMMON_CODE
!TOWER_SOIL_TEMPERATURE,         KIND_TEMPERATURE,              COMMON_CODE
!TOWER_U_WIND_COMPONENT,         KIND_U_WIND_COMPONENT,         COMMON_CODE
!TOWER_V_WIND_COMPONENT,         KIND_V_WIND_COMPONENT,         COMMON_CODE
!TOWER_GLOBAL_RADIATION,         KIND_RADIATION,                COMMON_CODE
!TOWER_NET_CARBON_FLUX,          KIND_NET_CARBON_FLUX,          COMMON_CODE
!TOWER_LATENT_HEAT_FLUX,         KIND_LATENT_HEAT_FLUX
!TOWER_SENSIBLE_HEAT_FLUX,       KIND_SENSIBLE_HEAT_FLUX
!TOWER_NETC_ECO_EXCHANGE,        KIND_NET_CARBON_PRODUCTION
! END DART PREPROCESS KIND LIST

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_tower_mod, only : get_scalar_from_3Dhistory
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(TOWER_LATENT_HEAT_FLUX)
!     call get_scalar_from_3Dhistory('EFLX_LH_TOT_R', state_time, ens_index, location, obs_time, obs_val, istatus)
!  case(TOWER_SENSIBLE_HEAT_FLUX)
!     call get_scalar_from_3Dhistory('FSH', state_time, ens_index, location, obs_time, obs_val, istatus)
!  case(TOWER_NETC_ECO_EXCHANGE)
!     call get_scalar_from_3Dhistory('NEP', state_time, ens_index, location, obs_time, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX,TOWER_SENSIBLE_HEAT_FLUX,TOWER_NETC_ECO_EXCHANGE)
!       continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX,TOWER_SENSIBLE_HEAT_FLUX,TOWER_NETC_ECO_EXCHANGE)
!       continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX,TOWER_SENSIBLE_HEAT_FLUX,TOWER_NETC_ECO_EXCHANGE)
!       continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_tower_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r4, r8, digits12, MISSING_R8, PI, deg2rad
use     location_mod, only : location_type, get_location
use time_manager_mod, only : time_type, get_date, set_date, print_date, print_time, &
                             get_time, set_time, operator(-), operator(/=)
use    utilities_mod, only : register_module, E_ERR, E_MSG, error_handler, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             nc_check, file_exist, is_longitude_between

use typesizes
use netcdf

implicit none
private

public :: get_scalar_from_3Dhistory

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical            :: module_initialized = .false.
character(len=129) :: string1, string2, string3
integer            :: nlon, nlat, ntime, ens_size
type(time_type)    :: initialization_time
real(r8)           :: edgeNorth,edgeEast,edgeSouth,edgeWest

character(len=129), allocatable, dimension(:) :: fname
integer,            allocatable, dimension(:) :: ncid
real(r8),           allocatable, dimension(:) :: lon, lat
real(digits12),     allocatable, dimension(:) :: rtime

! namelist items
character(len=256) :: casename = 'clm_tim'
logical            :: verbose = .false.
logical            :: debug = .false.

namelist /obs_def_tower_nml/ casename, verbose, debug

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module(model_time)
type(time_type), intent(in) :: model_time

! Called once to set values and allocate space, open all the CLM files
! that have the observations, etc.

integer :: iunit, io, i
integer :: dimid, varid
integer :: year, month, day, hour, minute, second, leftover
integer, allocatable, dimension(:) :: yyyymmdd,sssss
type(time_type) :: tower_time

! Prevent multiple calls from executing this code more than once.
if (module_initialized) then
   if (initialization_time /= model_time) then
      string1 = 'model time does not match initialization time'
      string2 = 'model time does not match initialization time'
      string3 = 'model time does not match initialization time'
      call error_handler(E_ERR, 'initialize_routine', string1, &
                     source, revision, revdate, text2=string2,text3=string3)
   endif
   return
else
   initialization_time = model_time
endif

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_tower_nml", iunit)
read(iunit, nml = obs_def_tower_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_tower_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_tower_nml)
if (do_nml_term()) write(     *     , nml=obs_def_tower_nml)

! Need to know what day we are trying to assimilate.
! The model stops at midnight, we want all the observations for THE PREVIOUS DAY.
! The CLM h0 files contain everything from 00:00 to 23:30 for the date in the filename.
! The data for [23:30 -> 00:00] get put in the file for the next day.

tower_time = model_time - set_time(0,1)
call get_date(tower_time, year, month, day, hour, minute, second)
second = second + minute*60 + hour*3600

! Figure out how many files (i.e. ensemble size) and construct their names.
! The CLM h0 files are constructed such that the midnight that starts the
! day is IN the file. The last time in the file is 23:30 ... 

100 format (A,'.clm2_',I4.4,'.h1.',I4.4,'-',I2.2,'-',I2.2,'-',I5.5,'.nc')

ens_size = 0
ENSEMBLESIZE : do i = 1,200

   write(string1,100) trim(casename),i,year,month,day,second

   if( file_exist(string1) ) then
      if(verbose .and. do_output()) write(*,*)'observation file "',trim(string1),'" exists.'
      ens_size = ens_size + 1
   else
      if(verbose .and. do_output()) write(*,*)'WARNING observation file "',trim(string1),'" does not exist.'
      exit ENSEMBLESIZE
   endif
enddo ENSEMBLESIZE

if (ens_size < 2) then

   write(string1,100) trim(casename),1,year,month,day,second
   write(string2,*)'cannot find files to use for observation operator.'
   write(string3,*)'trying files with names like "',trim(string1),'"'
   call error_handler(E_ERR, 'initialize_routine', string2, &
                  source, revision, revdate, text2=string3)

elseif (ens_size >= 200) then

   write(string2,*)'ensemble size (',ens_size,') is unnaturally large.'
   write(string3,*)'trying files with names like "',trim(string1),'"'
   call error_handler(E_ERR, 'initialize_routine', string2, &
                  source, revision, revdate, text2=string3)

else

   if (verbose .and. do_output()) write(*,*)'Ensemble size is believed to be ',ens_size

endif

allocate(fname(ens_size),ncid(ens_size))
ncid = 0

ENSEMBLE : do i = 1,ens_size
   write(fname(i),100) trim(casename),i,year,month,day,second
   call nc_check(nf90_open(trim(fname(i)), nf90_nowrite, ncid(i)), &
       'initialize_routine','open '//trim(fname(i)))
enddo ENSEMBLE

i = 1

! Harvest information from the first observation file.
! FIXME All other files will be opened to make sure they have the same dimensions.

call nc_check(nf90_inq_dimid(ncid(i), 'lon', dimid), &
       'initialize_routine','inq_dimid lon '//trim(fname(i)))
call nc_check(nf90_inquire_dimension(ncid(i), dimid, len=nlon), &
       'initialize_routine','inquire_dimension lon '//trim(fname(i)))

call nc_check(nf90_inq_dimid(ncid(i), 'lat', dimid), &
       'initialize_routine','inq_dimid lat '//trim(fname(i)))
call nc_check(nf90_inquire_dimension(ncid(i), dimid, len=nlat), &
       'initialize_routine','inquire_dimension lat '//trim(fname(i)))

call nc_check(nf90_inq_dimid(ncid(i), 'time', dimid), &
       'initialize_routine','inq_dimid time '//trim(fname(i)))
call nc_check(nf90_inquire_dimension(ncid(i), dimid, len=ntime), &
       'initialize_routine','inquire_dimension lat '//trim(fname(i)))

allocate(lon(nlon),lat(nlat))
allocate(rtime(ntime),yyyymmdd(ntime),sssss(ntime))

call nc_check(nf90_inq_varid(ncid(i), 'lon', varid), &
       'initialize_routine','inq_varid lon '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, lon), 'initialize_routine', 'get_var lon')

call nc_check(nf90_inq_varid(ncid(i), 'lat', varid), &
       'initialize_routine','inq_varid lat '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, lat), 'initialize_routine', 'get_var lat')

call nc_check(nf90_inq_varid(ncid(i), 'mcdate', varid), &
       'initialize_routine','inq_varid mcdate '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, yyyymmdd), 'initialize_routine', 'get_var yyyymmdd')

call nc_check(nf90_inq_varid(ncid(i), 'mcsec', varid), &
       'initialize_routine','inq_varid mcsec '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, sssss), 'initialize_routine', 'get_var sssss')

! Determine the geographic boundaries of the contents of the history file.

call nc_check(nf90_inq_varid(ncid(i), 'edgen', varid), &
       'initialize_routine','inq_varid edgen '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, edgeNorth), 'initialize_routine', 'get_var edgeNorth')

call nc_check(nf90_inq_varid(ncid(i), 'edgee', varid), &
       'initialize_routine','inq_varid edgee '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, edgeEast), 'initialize_routine', 'get_var edgeEast')

call nc_check(nf90_inq_varid(ncid(i), 'edges', varid), &
       'initialize_routine','inq_varid edges '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, edgeSouth), 'initialize_routine', 'get_var edgeSouth')

call nc_check(nf90_inq_varid(ncid(i), 'edgew', varid), &
       'initialize_routine','inq_varid edgew '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, edgeWest), 'initialize_routine', 'get_var edgeWest')

! call nc_check(nf90_inq_varid(ncid(i), 'time', varid), &
!        'initialize_routine','inq_varid time '//trim(fname(i)))
! call nc_check(nf90_get_var(ncid(i), varid, time), 'initialize_routine', 'get_var time')

! Convert time in file to a time compatible with the observation sequence file.
do i = 1,ntime

   year     = yyyymmdd(i)/10000
   leftover = yyyymmdd(i) - year*10000
   month    = leftover/100
   day      = leftover - month*100

   hour     = sssss(i)/3600
   leftover = sssss(i) - hour*3600
   minute   = leftover/60
   second   = leftover - minute*60

   tower_time = set_date(year, month, day, hour, minute, second)
   call get_time(tower_time, second, day)

   rtime(i) = real(day,digits12) + real(second,digits12)/86400.0_digits12

   if (debug .and. do_output()) then
      write(*,*)'timestep yyyymmdd sssss',i,yyyymmdd(i),sssss(i)
      call print_date(tower_time,'tower_mod date')
      call print_time(tower_time,'tower_mod time')
      write(*,*)'tower_mod time as a real ',rtime(i)
   endif

enddo

if (debug .and. do_output()) write(*,*)'obs_def_tower      lon',lon
if (debug .and. do_output()) write(*,*)'obs_def_tower      lat',lat

! FIXME
! check all other ensemble member history files to make sure metadata is the same.

deallocate(yyyymmdd, sssss)

end subroutine initialize_module



subroutine get_scalar_from_3Dhistory(varstring, state_time, ens_index, location, obs_time, obs_val, istatus )
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
type(time_type),     intent(in)  :: state_time
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
real(r8)               :: loc_lon, loc_lat
real(r4), dimension(1) :: hyperslab
real(r4)               :: spvalR4
real(r8)               :: scale_factor, add_offset
real(digits12)         :: otime
character(len=20)      :: strshort

if ( .not. module_initialized ) call initialize_module(state_time)

obs_val = MISSING_R8
istatus = 1

! if observation is outside region encompassed in the history file - fail
loc      = get_location(location) ! loc is in DEGREES
loc_lon  = loc(1)
loc_lat  = loc(2)

if ( .not. is_longitude_between(loc_lon, edgeWest, edgeEast, doradians=.FALSE.)) return
if ((loc_lat < edgeSouth) .or. (loc_lat > edgeNorth)) return

! Now that we know the observation operator is possible, continue ...

write(strshort,'(''ens_index '',i4,1x,A)')ens_index,trim(varstring)

if (ens_index > ens_size) then
   write(string1,*)'believed to have ',ens_size,'ensemble members for observation operator.'
   write(string2,*)'asking to use operator for ensemble member ',ens_index
   call error_handler(E_ERR, 'get_scalar_from_3Dhistory', &
              string1, source, revision, revdate, text2=string2)
endif

! bombproofing ... make sure the netcdf file is open.

write(*,*)'ncid(',ens_index,') is ',ncid(ens_index)
call nc_check(nf90_inquire(ncid(ens_index)), &
              'get_scalar_from_3Dhistory', 'inquire '//trim(strshort))

! bombproofing ... make sure the variable is the shape and size we expect

call nc_check(nf90_inq_varid(ncid(ens_index), trim(varstring), varid), &
        'get_scalar_from_3Dhistory', 'inq_varid '//trim(strshort))
call nc_check(nf90_inquire_variable(ncid(ens_index), varid, xtype=xtype, ndims=ndims, &
        dimids=dimids, natts=natts),'get_scalar_from_3Dhistory','inquire variable '//trim(strshort))

if (ndims /= 3) then
   write(string1,*)trim(varstring),' is supposed to have 3 dimensions, it has',ndims
   call error_handler(E_ERR, 'get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! If the variable is not a NF90_FLOAT, then the assumptions for processing
! the missing_value, _FillValue, etc., may not be correct.
if (xtype /= NF90_FLOAT) then
   write(string1,*)trim(varstring),' is supposed to be a 32 bit real. xtype = ',NF90_FLOAT,' it is ',xtype 
   call error_handler(E_ERR, 'get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 1 is longitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(1), len=dimlen), &
              'get_scalar_from_3Dhistory', 'inquire_dimension 1 '//trim(strshort))
if (dimlen /= nlon) then
   write(string1,*)'LON has length',nlon,trim(varstring),' has ',dimlen,'longitudes.'
   call error_handler(E_ERR, 'get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 2 is latitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(2), len=dimlen), &
              'get_scalar_from_3Dhistory', 'inquire_dimension 2 '//trim(strshort))
if (dimlen /= nlat) then
   write(string1,*)'LAT has length',nlat,trim(varstring),' has ',dimlen,'latitudes.'
   call error_handler(E_ERR, 'get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Dimension 3 is time
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(3), len=dimlen), &
              'get_scalar_from_3Dhistory', 'inquire_dimension 3'//trim(strshort))
if (dimlen /= ntime) then
   write(string1,*)'TIME has length',ntime,trim(varstring),' has ',dimlen,'times.'
   call error_handler(E_ERR, 'get_scalar_from_3Dhistory', &
              string1, source, revision, revdate)
endif

! Find the grid cell and timestep of interest 
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
   write(*,*)'get_scalar_from_3Dhistory:targetlon, lon, lon index is ', &
                                           loc_lon,lon(gridloni),gridloni
   write(*,*)'get_scalar_from_3Dhistory:targetlat, lat, lat index is ', &
                                           loc_lat,lat(gridlatj),gridlatj
   write(*,*)'get_scalar_from_3Dhistory:  targetT,   T,   T index is ', &
                                           otime,rtime(timei),timei
endif

if ( abs(otime - rtime(timei)) > 30*60 ) then
   if (debug .and. do_output()) then
      write(*,*)'get_scalar_from_3Dhistory: no close time ... skipping observation'
      call print_time(obs_time,'get_scalar_from_3Dhistory:observation time')
      call print_date(obs_time,'get_scalar_from_3Dhistory:observation date')
   endif
   istatus = 2
   return
endif

! Grab exactly the scalar we want.

ncstart = (/ gridloni, gridlatj, timei /)
nccount = (/        1,        1,     1 /)

call nc_check(nf90_get_var(ncid(ens_index), varid, hyperslab, start=ncstart, count=nccount), &
     'get_scalar_from_3Dhistory', 'get_var')

obs_val = hyperslab(1)

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


end module obs_def_tower_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
