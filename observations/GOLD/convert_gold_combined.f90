! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_gold_combined

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_gold_nmax - program that reads a GOLD netCDF TDISK and ON2
!                       observation files and writes a DART
!                       obs_seq file using the DART library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM (MADIS)
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM (MADIS)
!     modified Aug. 2020 George Bowden, UNSW Canberra
!
!     modified to use a common set of utilities, better netcdf error checks,
!     able to insert obs with any time correctly (not only monotonically
!     increasing times)    nancy collins,  ncar/image   11 march 2010
!     
!     keep original obs times, make source for all converters as similar
!     as possbile.   nancy collins,  ncar/image   26 march 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use      location_mod, only : VERTISHEIGHT, VERTISPRESSURE, VERTISUNDEF
use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read
use  time_manager_mod, only : time_type, set_calendar_type, set_date, operator(>=), &
                              increment_time, get_time, operator(-), GREGORIAN
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : GOLD_TEMPERATURE, GOLD_ON2COLUMN, GOLD_NEMAX
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              getvar_real_3d, &
                              create_3d_obs, getvar_int, getdimlen, set_missing_name

use netcdf

implicit none

!character(len=13),  parameter :: gold_netcdf_file = 'gold_input.nc'
!character(len=129), parameter :: gold_out_file    = 'obs_seq.gold'

! things which can/should be in the nc_to_obs_nml
character(len=128) :: gold_tdisk_netcdf_file = 'gold_tdisk_input.nc'
character(len=128) :: gold_on2_netcdf_file = 'gold_on2_input.nc'
character(len=128) :: gold_nemax_netcdf_file = 'gold_nemax_input.nc'
character(len=128) :: gold_out_file  = 'obs_seq.gold'
logical            :: use_tdisk = .true.
logical            :: use_on2 = .true.
logical            :: use_nemax = .true.
logical            :: debug = .true.

namelist /GOLD_combined_nc_to_obs_nml/  &
     gold_tdisk_netcdf_file, &
     gold_on2_netcdf_file,   &
     gold_nemax_netcdf_file, &
     gold_out_file,          &
     use_tdisk,              &
     use_on2,                &
     use_nemax,              &
     debug

logical, parameter :: use_input_qc              = .true.

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer :: ncid_tdisk, ncid_on2, ncid_nemax
integer :: nvars
integer :: i, oday, osec, nused, iunit, rcio
integer :: nobs_tdisk, nlat_tdisk, nlon_tdisk, nscan_tdisk, nmsl_tdisk
integer :: nobs_on2, nlat_on2, nlon_on2, nscan_on2, nmsl_on2
integer :: nobs_nemax, nlat_nemax, nlon_nemax, nscan_nemax, nmsl_nemax

integer :: k, l, m, n, j
integer :: varid
           
logical  :: file_exist, first_obs

real(r8) :: pres, height, qc, qerr, oerr
real(r8) :: tdisk_miss, on2_miss, nemax_miss

integer :: year_obs, month_obs, day_obs, hour_obs, minute_obs, second_obs

character, allocatable  :: tobs_tdisk(:,:,:,:)
real(r8), allocatable   :: lat_tdisk(:,:,:), lon_tdisk(:,:,:)
character, allocatable  :: tobs_on2(:,:,:,:)
real(r8), allocatable   :: lat_on2(:,:,:), lon_on2(:,:,:)
character, allocatable  :: tobs_nemax(:,:,:)
real(r8), allocatable   :: lat_nemax(:,:,:), lon_nemax(:,:,:)

real(r8), allocatable   :: tdisk(:,:,:), tdisk_unc_ran(:,:,:), &
     tdisk_unc_sys(:,:,:), tdisk_unc_mod(:,:,:)
real(r8), allocatable   :: on2(:,:,:), on2_unc_ran(:,:,:), &
     on2_unc_sys(:,:,:), on2_unc_mod(:,:,:)
real(r8), allocatable   :: nemax(:,:,:), nemax_unc_ran(:,:,:), &
     nemax_unc_sys(:,:,:), nemax_unc_mod(:,:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

character(len = 8), parameter :: varname_time_utc = 'time_utc'
character(len = 24)           :: tobs_string = '1970-01-01T00:00:00.000Z'


call initialize_utilities('convert_gold_combined')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'GOLD_combined_nc_to_obs_nml', iunit)
read(iunit, nml = GOLD_combined_nc_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'GOLD_combined_nc_to_obs_nml')

!print*,gold_netcdf_file

! put the reference date into DART format
call set_calendar_type(GREGORIAN)

first_obs = .true.

call nc_check( nf90_open(gold_tdisk_netcdf_file, nf90_nowrite, ncid_tdisk), &
     'convert_gold_combined', 'opening file '//trim(gold_tdisk_netcdf_file) )
call nc_check( nf90_open(gold_on2_netcdf_file, nf90_nowrite, ncid_on2), &
     'convert_gold_combined', 'opening file '//trim(gold_on2_netcdf_file) )
call nc_check( nf90_open(gold_nemax_netcdf_file, nf90_nowrite, ncid_nemax), &
     'convert_gold_combined', 'opening file '//trim(gold_nemax_netcdf_file) )

call getdimlen(ncid_tdisk, "nlats", nlat_tdisk)
call getdimlen(ncid_tdisk, "nlons", nlon_tdisk)
call getdimlen(ncid_tdisk, "nscans", nscan_tdisk)
call getdimlen(ncid_tdisk, "max_string_len ", nmsl_tdisk)
call set_missing_name("missing_value")

call getdimlen(ncid_on2, "nlats", nlat_on2)
call getdimlen(ncid_on2, "nlons", nlon_on2)
call getdimlen(ncid_on2, "nscans", nscan_on2)
call getdimlen(ncid_on2, "max_string_len ", nmsl_on2)
call set_missing_name("missing_value")

call getdimlen(ncid_nemax, "nlats", nlat_nemax)
call getdimlen(ncid_nemax, "nlons", nlon_nemax)
call getdimlen(ncid_nemax, "nscans", nscan_nemax)
call getdimlen(ncid_nemax, "max_string_len ", nmsl_nemax)
call set_missing_name("missing_value")

nobs_tdisk = nlat_tdisk * nlon_tdisk * nscan_tdisk
nobs_on2 = nlat_on2 * nlon_on2 * nscan_on2
nobs_nemax = nlat_nemax * nlon_nemax * nscan_nemax
nvars = 1
! Check meaning of nvars

! Allocate arrays for tdisk
allocate( lat_tdisk(nlon_tdisk,nlat_tdisk,nscan_tdisk))
allocate( lon_tdisk(nlon_tdisk,nlat_tdisk,nscan_tdisk))
allocate( tdisk(nlon_tdisk,nlat_tdisk,nscan_tdisk))
allocate( tdisk_unc_ran(nlon_tdisk,nlat_tdisk,nscan_tdisk))
allocate( tdisk_unc_sys(nlon_tdisk,nlat_tdisk,nscan_tdisk))
allocate( tdisk_unc_mod(nlon_tdisk,nlat_tdisk,nscan_tdisk))
!allocate( tobs(nmsl,nlon,nscan))
allocate( tobs_tdisk(24,nlon_tdisk,nlat_tdisk,nscan_tdisk))

! Allocate arrays for on2
allocate( lat_on2(nlon_on2,nlat_on2,nscan_on2))
allocate( lon_on2(nlon_on2,nlat_on2,nscan_on2))
allocate( on2(nlon_on2,nlat_on2,nscan_on2))
allocate( on2_unc_ran(nlon_on2,nlat_on2,nscan_on2))
allocate( on2_unc_sys(nlon_on2,nlat_on2,nscan_on2))
allocate( on2_unc_mod(nlon_on2,nlat_on2,nscan_on2))
!allocate( tobs(nmsl,nlon,nscan))
allocate( tobs_on2(24,nlon_on2,nlat_on2,nscan_on2))

! Allocate arrays for nemax
allocate( lat_nemax(nlon_nemax,nlat_nemax,nscan_nemax))
allocate( lon_nemax(nlon_nemax,nlat_nemax,nscan_nemax))
allocate( nemax(nlon_nemax,nlat_nemax,nscan_nemax))
allocate( nemax_unc_ran(nlon_nemax,nlat_nemax,nscan_nemax))
allocate( nemax_unc_sys(nlon_nemax,nlat_nemax,nscan_nemax))
allocate( nemax_unc_mod(nlon_nemax,nlat_nemax,nscan_nemax))
!allocate( tobs(nmsl,nlon,nscan))
allocate( tobs_nemax(24,nlon_nemax,nscan_nemax))

! read in the data arrays for tdisk
call getvar_real_3d(ncid_tdisk, "latitude",      lat_tdisk) ! latitudes
call getvar_real_3d(ncid_tdisk, "longitude",     lon_tdisk) ! longitudes
call getvar_real_3d(ncid_tdisk, "tdisk",         tdisk,         tdisk_miss) ! maximum electron concentration
call getvar_real_3d(ncid_tdisk, "tdisk_unc_ran", tdisk_unc_ran, tdisk_miss) ! maximum electron concentration
call getvar_real_3d(ncid_tdisk, "tdisk_unc_sys", tdisk_unc_sys, tdisk_miss) ! maximum electron concentration
call getvar_real_3d(ncid_tdisk, "tdisk_unc_mod", tdisk_unc_mod, tdisk_miss) ! maximum electron concentration

! read in the data arrays for on2
call getvar_real_3d(ncid_on2, "latitude",     lat_on2) ! latitudes
call getvar_real_3d(ncid_on2, "longitude",    lon_on2) ! longitudes
call getvar_real_3d(ncid_on2, "on2",          on2,         on2_miss) ! maximum electron concentration
call getvar_real_3d(ncid_on2, "on2_unc_ran", on2_unc_ran, on2_miss) ! maximum electron concentration
call getvar_real_3d(ncid_on2, "on2_unc_sys", on2_unc_sys, on2_miss) ! maximum electron concentration
call getvar_real_3d(ncid_on2, "on2_unc_mod", on2_unc_mod, on2_miss) ! maximum electron concentration

! read in the data arrays for nemax
call getvar_real_3d(ncid_nemax, "latitude",     lat_nemax) ! latitudes
call getvar_real_3d(ncid_nemax, "longitude",    lon_nemax) ! longitudes
call getvar_real_3d(ncid_nemax, "nmax",         nemax,         nemax_miss) ! maximum electron concentration
call getvar_real_3d(ncid_nemax, "nmax_unc_ran", nemax_unc_ran, nemax_miss) ! maximum electron concentration
call getvar_real_3d(ncid_nemax, "nmax_unc_sys", nemax_unc_sys, nemax_miss) ! maximum electron concentration
call getvar_real_3d(ncid_nemax, "nmax_unc_mod", nemax_unc_mod, nemax_miss) ! maximum electron concentration

! read the data for the time array tdisk
call nc_check( nf90_inq_varid(ncid_tdisk, varname_time_utc, varid), &
               'getvar_char', 'inquire var '// trim(varname_time_utc))
call nc_check( nf90_get_var(ncid_tdisk, varid, tobs_tdisk, &
               start = (/ 1, 1, 1, 1 /), count = (/ 24, nlon_tdisk, nlat_tdisk, nscan_tdisk /)), &
               'getvar_char', 'getting var '// trim(varname_time_utc))

! read the data for the time array on2
call nc_check( nf90_inq_varid(ncid_on2, varname_time_utc, varid), &
               'getvar_char', 'inquire var '// trim(varname_time_utc))
call nc_check( nf90_get_var(ncid_on2, varid, tobs_on2, &
               start = (/ 1, 1, 1, 1 /), count = (/ 24, nlon_on2, nlat_on2, nscan_on2 /)), &
               'getvar_char', 'getting var '// trim(varname_time_utc))

! read the data for the time array nemax
call nc_check( nf90_inq_varid(ncid_nemax, varname_time_utc, varid), &
               'getvar_char', 'inquire var '// trim(varname_time_utc))
call nc_check( nf90_get_var(ncid_nemax, varid, tobs_nemax, &
               start = (/ 1, 1, 1 /), count = (/ 24, nlon_nemax, nscan_nemax /)), &
               'getvar_char', 'getting var '// trim(varname_time_utc))

call nc_check( nf90_close(ncid_tdisk), &
     'convert_gold_combined', 'closing file '//trim(gold_tdisk_netcdf_file) )
call nc_check( nf90_close(ncid_on2), &
     'convert_gold_combined', 'closing file '//trim(gold_on2_netcdf_file) )
call nc_check( nf90_close(ncid_nemax), &
     'convert_gold_combined', 'closing file '//trim(gold_nemax_netcdf_file) )

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=gold_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
   call read_obs_seq(gold_out_file, 0, 0, nvars*(nobs_tdisk + nobs_on2 + nobs_nemax), obs_seq)

else

   ! create a new one
   call init_obs_sequence(obs_seq, num_copies, num_qc, nvars*(nobs_tdisk + nobs_on2 + nobs_nemax))
   do i = 1, num_copies
      call set_copy_meta_data(obs_seq, i, 'GOLD observation')
   end do
   do i = 1, num_qc
      call set_qc_meta_data(obs_seq, i, 'Data QC')
   end do

endif

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

year_obs = 1970
month_obs = 1
day_obs = 1
hour_obs = 0
minute_obs = 0
second_obs = 0

nused = 0
scanloop_tdisk: do k = 1, nscan_tdisk
   lonloop_tdisk: do l = 1, nlon_tdisk
      latloop_tdisk: do m = 1, nlat_tdisk

         if ('*' == tobs_tdisk(1,l,m,k)) cycle latloop_tdisk
         do j = 1, 24
            tobs_string(j:j) = tobs_tdisk(j,l,m,k)
         end do

         read(tobs_string(1:4),*)   year_obs
         read(tobs_string(6:7),*)   month_obs
         read(tobs_string(9:10),*)  day_obs
         read(tobs_string(12:13),*) hour_obs
         read(tobs_string(15:16),*) minute_obs
         read(tobs_string(18:23),*) second_obs

         time_obs = set_date(year_obs, month_obs, day_obs, hour_obs, minute_obs, second_obs)
         
         if (.not. (lon_tdisk(l,m,k) < 180.0_r8 .and. lon_tdisk(l,m,k) > -180.0_r8)) cycle latloop_tdisk
         if (.not. (lat_tdisk(l,m,k) <  90.0_r8 .and. lat_tdisk(l,m,k) >  -90.0_r8)) cycle latloop_tdisk

         if ( lon_tdisk(l,m,k) < 0.0_r8 )  lon_tdisk(l,m,k) = lon_tdisk(l,m,k) + 360.0_r8
         
         ! extract actual time of observation in file into oday, osec.
         call get_time(time_obs, osec, oday)
         
         ! Fixed height of 150 km
         height = 1.5E5
         
         !oerr = (tdisk_unc_ran(l,m,k)**2.0 + tdisk_unc_sys(l,m,k)**2.0)**0.5
         oerr = tdisk_unc_ran(l,m,k)

         !if (not(isnan(tdisk(l,m,k)))) then
         if ((not(isnan(tdisk(l,m,k)))) .and. (not(isnan(oerr)))) then
            call create_3d_obs(lat_tdisk(l,m,k), lon_tdisk(l,m,k), height, VERTISHEIGHT, tdisk(l,m,k), &
                 GOLD_TEMPERATURE, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
         endif
         
      end do latloop_tdisk
   end do lonloop_tdisk
end do scanloop_tdisk

scanloop_on2: do k = 1, nscan_on2
   lonloop_on2: do l = 1, nlon_on2
      latloop_on2: do m = 1, nlat_on2

         if ('*' == tobs_on2(1,l,m,k)) cycle latloop_on2
         do j = 1, 24
            tobs_string(j:j) = tobs_on2(j,l,m,k)
         end do

         read(tobs_string(1:4),*)   year_obs
         read(tobs_string(6:7),*)   month_obs
         read(tobs_string(9:10),*)  day_obs
         read(tobs_string(12:13),*) hour_obs
         read(tobs_string(15:16),*) minute_obs
         read(tobs_string(18:23),*) second_obs

         time_obs = set_date(year_obs, month_obs, day_obs, hour_obs, minute_obs, second_obs)
         
         if (.not. (lon_on2(l,m,k) < 180.0_r8 .and. lon_on2(l,m,k) > -180.0_r8)) cycle latloop_on2
         if (.not. (lat_on2(l,m,k) <  90.0_r8 .and. lat_on2(l,m,k) >  -90.0_r8)) cycle latloop_on2

         if ( lon_on2(l,m,k) < 0.0_r8 )  lon_on2(l,m,k) = lon_on2(l,m,k) + 360.0_r8
         
         ! extract actual time of observation in file into oday, osec.
         call get_time(time_obs, osec, oday)
         
         ! Not sure what to do about altitude coordinate, not included in GOLD file
         !pres = 5.0E-5
         ! Fixed height of 150 km
         height = 1.5E5
         
         oerr = (on2_unc_ran(l,m,k)**2.0 + on2_unc_sys(l,m,k)**2.0)**0.5

         if ((not(isnan(on2(l,m,k)))) .and. (not(isnan(oerr)))) then 
            !call create_3d_obs(lat_on2(l,m,k), lon_on2(l,m,k), pres, VERTISPRESSURE, on2(l,m,k), &
            !     GOLD_ON2COLUMN, oerr, oday, osec, qc, obs)
            call create_3d_obs(lat_on2(l,m,k), lon_on2(l,m,k), height, VERTISUNDEF , on2(l,m,k), &
                 GOLD_ON2COLUMN, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
         endif
         
      end do latloop_on2
   end do lonloop_on2
end do scanloop_on2

scanloop_nemax: do k = 1, nscan_nemax
   lonloop_nemax: do l = 1, nlon_nemax

      if ('*' == tobs_nemax(1,l,k)) cycle lonloop_nemax
      do j = 1, 24
         tobs_string(j:j) = tobs_nemax(j,l,k)
      end do

      read(tobs_string(1:4),*)   year_obs
      read(tobs_string(6:7),*)   month_obs
      read(tobs_string(9:10),*)  day_obs
      read(tobs_string(12:13),*) hour_obs
      read(tobs_string(15:16),*) minute_obs
      read(tobs_string(18:23),*) second_obs
      
      time_obs = set_date(year_obs, month_obs, day_obs, hour_obs, minute_obs, second_obs)
      
      latloop_nemax: do m = 1, nlat_nemax
         
         if (.not. (lon_nemax(l,m,k) < 180.0_r8 .and. lon_nemax(l,m,k) > -180.0_r8)) cycle latloop_nemax
         if (.not. (lat_nemax(l,m,k) <  90.0_r8 .and. lat_nemax(l,m,k) >  -90.0_r8)) cycle latloop_nemax

         if ( lon_nemax(l,m,k) < 0.0_r8 )  lon_nemax(l,m,k) = lon_nemax(l,m,k) + 360.0_r8
         
         ! extract actual time of observation in file into oday, osec.
         call get_time(time_obs, osec, oday)
         
         ! Not sure what to do about altitude coordinate, not included in GOLD file
         !pres = 5.0E-5
         ! Fixed height of 150 km
         height = 1.5E5
         
         oerr = (nemax_unc_ran(l,m,k)**2.0 + nemax_unc_sys(l,m,k)**2.0)**0.5

         if ((not(isnan(nemax(l,m,k)))) .and. (not(isnan(oerr)))) then 
            !call create_3d_obs(lat_nemax(l,m,k), lon_nemax(l,m,k), pres, VERTISPRESSURE, nemax(l,m,k), &
            !     GOLD_NEMAX, oerr, oday, osec, qc, obs)
            call create_3d_obs(lat_nemax(l,m,k), lon_nemax(l,m,k), height, VERTISUNDEF , nemax(l,m,k), &
                 GOLD_NEMAX, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
         endif
         
      end do latloop_nemax
   end do lonloop_nemax
end do scanloop_nemax

100 continue

print*,"Finished processing, about to write"

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, gold_out_file)

! end of main program
call finalize_utilities()


end program convert_gold_combined

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
