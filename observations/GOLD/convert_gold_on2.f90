! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_gold_on2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_gold_nmax - program that reads a GOLD netCDF O/N2
!                       observation file and writes a DART
!                       obs_seq file using the DART library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM (MADIS)
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM (MADIS)
!     modified Jun. 2020 George Bowden, UNSW Canberra
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
use      location_mod, only : VERTISPRESSURE, VERTISHEIGHT, VERTISUNDEF
use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read
use  time_manager_mod, only : time_type, set_calendar_type, set_date, operator(>=), &
                              increment_time, get_time, operator(-), GREGORIAN
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : GOLD_ON2COLUMN
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              getvar_real_3d, getvar_real_2d, &
                              create_3d_obs, getvar_int, getdimlen, set_missing_name

use netcdf

implicit none

!character(len=13),  parameter :: gold_netcdf_file = 'gold_input.nc'
!character(len=129), parameter :: gold_out_file    = 'obs_seq.gold'

! things which can/should be in the nc_to_obs_nml
character(len=128) :: gold_netcdf_file = 'gold_input.nc'
character(len=128) :: gold_out_file  = 'obs_seq.gold'
logical            :: debug = .true.

namelist /GOLD_on2_nc_to_obs_nml/  &
     gold_netcdf_file, &
     gold_out_file,    &
     debug

logical, parameter :: use_input_qc              = .true.

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer :: ncid, nvars
integer :: i, oday, osec, nused, iunit, rcio
integer :: nobs, nlat, nlon, nscan, nmsl
integer :: k, l, m, n, j
integer :: varid
           
logical  :: file_exist, first_obs

real(r8) :: pres, height, qc, qerr, oerr, on2_miss

integer :: year_obs, month_obs, day_obs, hour_obs, minute_obs, second_obs

character, allocatable  :: tobs(:,:,:,:)
real(r8), allocatable   :: lat(:,:), lon(:,:)
real(r8), allocatable   :: on2(:,:,:), on2_unc_ran(:,:,:), &
                           on2_unc_sys(:,:,:), on2_unc_mod(:,:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

character(len = 8), parameter :: varname_time_utc = 'time_utc'
character(len = 24)           :: tobs_string = '1970-01-01T00:00:00.000Z'


call initialize_utilities('convert_gold_on2')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'GOLD_on2_nc_to_obs_nml', iunit)
read(iunit, nml = GOLD_on2_nc_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, 'GOLD_on2_nc_to_obs_nml')

!print*,gold_netcdf_file

! put the reference date into DART format
call set_calendar_type(GREGORIAN)

first_obs = .true.

call nc_check( nf90_open(gold_netcdf_file, nf90_nowrite, ncid), &
               'convert_gold_on2', 'opening file '//trim(gold_netcdf_file) )

call getdimlen(ncid, "nlats", nlat)
call getdimlen(ncid, "nlons", nlon)
call getdimlen(ncid, "nscans", nscan)
call getdimlen(ncid, "max_string_len ", nmsl)
call set_missing_name("missing_value")

nobs = nlat * nlon * nscan
nvars = 1

allocate( lat(nlon,nlat))
allocate( lon(nlon,nlat))
allocate( on2(nlon,nlat,nscan))
allocate( on2_unc_ran(nlon,nlat,nscan))
allocate( on2_unc_sys(nlon,nlat,nscan))
allocate( on2_unc_mod(nlon,nlat,nscan))
!allocate( tobs(nmsl,nlon,nscan))
allocate( tobs(24,nlon,nlat,nscan))

! read in the data arrays
call getvar_real_2d(ncid, "latitude",     lat) ! latitudes
call getvar_real_2d(ncid, "longitude",    lon) ! longitudes
call getvar_real_3d(ncid, "on2",         on2,         on2_miss) ! Column O/N2
call getvar_real_3d(ncid, "on2_unc_ran", on2_unc_ran, on2_miss) ! Column O/N2 random uncertainty
call getvar_real_3d(ncid, "on2_unc_sys", on2_unc_sys, on2_miss) ! Column O/N2 systematic uncertainty
call getvar_real_3d(ncid, "on2_unc_mod", on2_unc_mod, on2_miss) ! Column O/N2 model uncertainty

! read the data for the time array
call nc_check( nf90_inq_varid(ncid, varname_time_utc, varid), &
               'getvar_char', 'inquire var '// trim(varname_time_utc))
!call nc_check( nf90_get_var(ncid, varid, tobs), &
!               'getvar_char', 'getting var '// trim(varname_time_utc))
call nc_check( nf90_get_var(ncid, varid, tobs, &
               start = (/ 1, 1, 1, 1 /), count = (/ 24, nlon, nlat, nscan /)), &
               'getvar_char', 'getting var '// trim(varname_time_utc))

call nc_check( nf90_close(ncid), &
     'convert_gold_on2', 'closing file '//trim(gold_netcdf_file) )

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=gold_out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
   call read_obs_seq(gold_out_file, 0, 0, nvars*nobs, obs_seq)

else

   ! create a new one
   call init_obs_sequence(obs_seq, num_copies, num_qc, nvars*nobs)
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
scanloop: do k = 1, nscan
   lonloop: do l = 1, nlon
      latloop: do m = 1, nlat

         if ('*' == tobs(1,l,m,k)) cycle latloop
         do j = 1, 24
            tobs_string(j:j) = tobs(j,l,m,k)
         end do

         read(tobs_string(1:4),*)   year_obs
         read(tobs_string(6:7),*)   month_obs
         read(tobs_string(9:10),*)  day_obs
         read(tobs_string(12:13),*) hour_obs
         read(tobs_string(15:16),*) minute_obs
         read(tobs_string(18:23),*) second_obs

         time_obs = set_date(year_obs, month_obs, day_obs, hour_obs, minute_obs, second_obs)
         
         !if (.not. (lon(l,m) < 180.0_r8 .and. lon(l,m) > -180.0_r8)) cycle latloop
         !if (.not. (lat(l,m) <  90.0_r8 .and. lat(l,m) >  -90.0_r8)) cycle latloop

         if ( lon(l,m) < 0.0_r8 )  lon(l,m) = lon(l,m) + 360.0_r8
         
         ! extract actual time of observation in file into oday, osec.
         call get_time(time_obs, osec, oday)
         
         ! Not sure what to do about altitude coordinate, not included in GOLD file
         !pres = 5.0E-5
         ! Fixed height of 150 km
         height = 1.5E5
         
         !oerr = (on2_unc_ran(l,m,k)**2.0 + on2_unc_sys(l,m,k)**2.0 + on2_unc_mod(l,m,k)**2.0)**0.5
         oerr = (on2_unc_ran(l,m,k)**2.0 + on2_unc_sys(l,m,k)**2.0)**0.5
         !oerr = on2_unc_ran(l,m,k)

         if ((not(isnan(on2(l,m,k)))) .and. (not(isnan(oerr)))) then 
            !call create_3d_obs(lat(l,m,k), lon(l,m,k), pres, VERTISPRESSURE, on2(l,m,k), &
            !     GOLD_ON2COLUMN, oerr, oday, osec, qc, obs)
            call create_3d_obs(lat(l,m), lon(l,m), height, VERTISUNDEF, on2(l,m,k), &
                 GOLD_ON2COLUMN, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
         endif
         
      end do latloop
   end do lonloop
end do scanloop

100 continue

print*,"Finished processing, about to write"

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, gold_out_file)

! end of main program
call finalize_utilities()


end program convert_gold_on2

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
