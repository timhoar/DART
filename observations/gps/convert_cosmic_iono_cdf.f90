! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> a version of the netcdf -> dart obs_seq converter for ionosphere profiles
!> the file type from the CDAAC data site is 'ionPrf'.
!>
!> i don't know where the obs errors are going to come from. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_cosmic_gps_cdf - program that reads a COSMIC GPS observation 
!                            profile and writes the data to a DART 
!                            observation sequence file. 
!
!     created June 2008 Ryan Torn, NCAR/MMM
!     modify  Fab. 2010  I-TE LEE, NCAR/HAO & NCU
!                       Modify to read the ionPrf data file of COSMIC
!                       and write the data to the DART format.  
!             Oct. 2010  I-TE LEE, NCAR/HAO & NCU
!                       Remove the code for lower atmospheric observations
!             Dec. 2010  I-TE LEE, NCAR/HAO & NCU
!                       Add new error valur for log scale testing
!             Jan. 2011  I-TE LEE, NCAR/HAO & NCU
!                       Modify for real GOX observations
!             Jan. 2011  I-TE LEE, NACR/HAO & NCU
!                       Add subroutine to calculate the observation varience 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!export MALLOC_CHECK_=0

program convert_cosmic_iono_cdf

use        types_mod, only : r8, metadatalength
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,&
                             increment_time, get_time, set_date, operator(-),  &
                             print_date
use    utilities_mod, only : initialize_utilities, find_namelist_in_file,    &
                             check_namelist_read, nmlfileunit, do_nml_file,   &
                             get_next_filename, error_handler, E_ERR, E_MSG, &
                             nc_check, find_textfile_dims, do_nml_term
use     location_mod, only : VERTISHEIGHT, set_location
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq,       &
                             static_init_obs_sequence, init_obs, destroy_obs, &
                             write_obs_seq, init_obs_sequence, get_num_obs,   &
                             insert_obs_in_seq, destroy_obs_sequence,         &
                             set_copy_meta_data, set_qc_meta_data, set_qc,    & 
                             set_obs_values, set_obs_def, insert_obs_in_seq
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, &
                             set_obs_def_key
use     obs_kind_mod, only : GPS_PROFILE 

use           netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
!
!Set the value equal to 0 for the OSSE testing, for real case it should be large then 1
!
integer, parameter ::   num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=metadatalength) :: meta_data
character (len=129) :: msgstring, next_infile
character (len=80)  :: name
character (len=19)  :: datestr
character (len=6)   :: subset
integer :: rcode, ncid, varid, nlevels, k, nfiles, num_new_obs,  &
           aday, asec, dday, dsec, oday, osec,                   &
           iyear, imonth, iday, ihour, imin, isec,               &
           glat, glon, zloc, obs_num, io, iunit, nobs, filenum, dummy
logical :: file_exist, first_obs, did_obs, from_list = .false.
real(r8) :: hght_miss,  elec_miss, err, oerr,               & 
            qc, lato, lono, hghto, eleco, wght, nx, ny,   & 
            nz, ds, htop, rfict, obsval, phs, obs_val(1), qc_val(1)
            

real(r8), allocatable :: lat(:), lon(:), hght(:), elec(:), & 
                         hghtp(:) 
!add F3/C ionprf observation error /add by ITL 2011.01.31
real(r8), dimension(15,37,25) :: f3coerr

type(obs_def_type)      :: obs_def
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, time_anal

!------------------------------------------------------------------------
!  Declare namelist parameters
!------------------------------------------------------------------------

integer, parameter :: nmaxlevels = 20000   !  max number of observation levels

logical  :: local_operator = .true.   ! see html file for more on non/local
logical  :: overwrite_time = .true.   !false.  ! careful - see note below
real(r8) :: obs_levels(nmaxlevels) = -1.0_r8
!real(r8) :: obs_window = 0.250      ! accept obs within +/- hours from anal time
real(r8) :: obs_window = 0.5
character(len=128) :: gpsro_netcdf_file     = 'cosmic_gps_input.nc'
character(len=128) :: gpsro_netcdf_filelist = 'cosmic_gps_input_list'
character(len=128) :: gpsro_out_file        = 'obs_seq.gpsro'

namelist /convert_cosmic_iono_nml/ obs_levels, local_operator, obs_window, &
                                  gpsro_netcdf_file,    &
                                  gpsro_netcdf_filelist, gpsro_out_file 


! 'overwrite_time' replaces the actual observation times with the
! analysis time for all obs.  this is intentionally not in the namelist
! because we think observations should preserve the original times in all
! cases.  but one of our users had a special request to overwrite the time
! with the synoptic time, so there is code in this converter to do that
! if this option is set to .true.  but we still don't encourage its use.

! initialize some values
obs_num = 1
qc = 0.0_r8

print*,'Enter the target assimilation time (yyyy-mm-dd_hh:mm:ss)'
read*,datestr

call set_calendar_type(GREGORIAN)
read(datestr(1:4),   fmt='(i4)') iyear
read(datestr(6:7),   fmt='(i2)') imonth
read(datestr(9:10),  fmt='(i2)') iday
read(datestr(12:13), fmt='(i2)') ihour
read(datestr(15:16), fmt='(i2)') imin
read(datestr(18:19), fmt='(i2)') isec
time_anal = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(time_anal, asec, aday)

!  read the necessary parameters from input.nml
call initialize_utilities()

!call find_namelist_in_file("work/input.nml", "convert_cosmic_iono_nml", iunit)
call find_namelist_in_file("input.nml", "convert_cosmic_iono_nml", iunit)
read(iunit, nml = convert_cosmic_iono_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_cosmic_iono_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=convert_cosmic_iono_nml)
if (do_nml_term()) write(     *     , nml=convert_cosmic_iono_nml)

! namelist checks for sanity

!  count observation levels, make sure observation levels increase from 0
nlevels = 0
do k = 1, nmaxlevels
  if ( obs_levels(k) == -1.0_r8 )  exit
  nlevels = k
end do
do k = 2, nlevels
  if ( obs_levels(k-1) >= obs_levels(k) ) then
    call error_handler(E_ERR, 'convert_cosmic_gps_cdf',       &
                       'Observation levels should increase',  &
                       source, revision, revdate)
  end if
end do

!  should error check the window some
if (obs_window <= 0.0_r8 .or. obs_window > 24.0_r8) then
    call error_handler(E_ERR, 'convert_cosmic_gps_cdf',       &
                       'Bad value for obs_window (hours)',    &
                       source, revision, revdate)
else
   ! convert to seconds
   obs_window = obs_window * 3600.0_r8
endif


! cannot have both a single filename and a list; the namelist must
! shut one off.
if (gpsro_netcdf_file /= '' .and. gpsro_netcdf_filelist /= '') then
  call error_handler(E_ERR, 'convert_cosmic_gps_cdf',                     &
                     'One of gpsro_netcdf_file or filelist must be NULL', &
                     source, revision, revdate)
endif
if (gpsro_netcdf_filelist /= '') from_list = .true.

! need to know a reasonable max number of obs that could be added here.
if (from_list) then
   call find_textfile_dims(gpsro_netcdf_filelist, nfiles, dummy)
   num_new_obs = nlevels * nfiles
else
   num_new_obs = nlevels
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
inquire(file=gpsro_out_file, exist=file_exist)
if ( file_exist ) then

print *, "found existing obs_seq file, appending to ", trim(gpsro_out_file)
   call read_obs_seq(gpsro_out_file, 0, 0, num_new_obs, obs_seq)

else

  print *, "no existing obs_seq file, creating ", trim(gpsro_out_file)
  print *, "max entries = ", num_new_obs
  call init_obs_sequence(obs_seq, num_copies, num_qc, num_new_obs)
  do k = 1, num_copies
    call set_copy_meta_data(obs_seq, k, 'observations')
  end do
  do k = 1, num_qc
    call set_qc_meta_data(obs_seq, k, 'Ionprf QC')
  end do

end if

!Reading ionprf observation error /add by ITL 2011.01.31
!open(16,FILE='f3coerr.dat',STATUS='old',FORM='FORMATTED')
!read(16,*) f3coerr
open(16,FILE='obserr.dat',STATUS='old',FORM='FORMATTED')
read(16,*) f3coerr
!!

allocate(hghtp(nlevels))  ;  
did_obs = .false.

! main loop that does either a single file or a list of files

filenum = 1
fileloop: do      ! until out of files

   ! get the single name, or the next name from a list
   if (from_list) then 
      next_infile = get_next_filename(gpsro_netcdf_filelist, filenum)
   else
      next_infile = gpsro_netcdf_file
      if (filenum > 1) next_infile = ''
   endif
   if (next_infile == '') exit fileloop
  
   !  open the occultation profile, check if it is within the window
   call nc_check( nf90_open(next_infile, nf90_nowrite, ncid), 'file open', next_infile)
   call nc_check( nf90_get_att(ncid,nf90_global,'year',  iyear) ,'get_att year')
   call nc_check( nf90_get_att(ncid,nf90_global,'month', imonth),'get_att month')
   call nc_check( nf90_get_att(ncid,nf90_global,'day',   iday)  ,'get_att day')
   call nc_check( nf90_get_att(ncid,nf90_global,'hour',  ihour) ,'get_att hour')
   call nc_check( nf90_get_att(ncid,nf90_global,'minute',imin)  ,'get_att minute')
   call nc_check( nf90_get_att(ncid,nf90_global,'second',isec)  ,'get_att second')
   
   time_obs = set_date(iyear, imonth, iday, ihour, imin, isec)
   call get_time(time_obs,  osec, oday)
   
   !time1-time2 is always positive no matter the relative magnitudes
   call get_time(time_anal-time_obs,dsec,dday)
   if ( real(dsec+dday*86400) > obs_window ) then
     call error_handler(E_MSG, 'convert_cosmic_gps_cdf: ',         &
                       'Input file '//trim(next_infile), &
                         source, revision, revdate)
     write(msgstring, '(A,F8.4,A)') 'Ignoed because obs time > ', &
                       obs_window / 3600.0, ' hours from analysis time'
     call error_handler(E_MSG, '', msgstring,        &
                        source, revision, revdate)
     call print_date(time_obs,  '    observation time: ')
     call print_date(time_anal, '  window center time: ')

     filenum = filenum + 1
     cycle fileloop
   end if
   print *,next_infile
 
   call nc_check( nf90_inq_dimid(ncid, "MSL_alt", varid), 'inq dimid MSL_alt')
   call nc_check( nf90_inquire_dimension(ncid, varid, name, nobs), 'inq dim MSL_alt')
   print *,nobs 
   allocate( lat(nobs))  ;  allocate( lon(nobs))
   allocate(hght(nobs))  ;  allocate(elec(nobs)) 
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "GEO_lat", varid) ,'inq varid GEO_lat')
   call nc_check( nf90_get_var(ncid, varid, lat)         ,'get var   GEO_lat')
   
   ! read the latitude array
   call nc_check( nf90_inq_varid(ncid, "GEO_lon", varid) ,'inq varid GEO_lon')
   call nc_check( nf90_get_var(ncid, varid, lon)         ,'get var   GEO_lon')
   
   ! read the altitude array
   call nc_check( nf90_inq_varid(ncid, "MSL_alt", varid) ,'inq varid MSL_alt')
   call nc_check( nf90_get_var(ncid, varid, hght)        ,'get_var   MSL_alt')
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', hght_miss) ,'get_att _FillValue MSL_alt')
   
   ! read the electorn density profiles (modify by ITL)
   call nc_check( nf90_inq_varid(ncid, "ELEC_dens", varid) ,'inq varid ELEC_dens')
   call nc_check( nf90_get_var(ncid, varid, elec)          ,'get var   ELEC_dens')
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', elec_miss) ,'get_att _FillValue Elec')

   ! check for the data quality
   if (maxval(elec) > 10000000) then
     call error_handler(E_MSG, 'Bad convert_cosmic_gps_cdf: ',         &
                       'Input file '//trim(next_infile), &
                         source, revision, revdate)
     deallocate( lat, lon, hght, elec )
     filenum = filenum + 1
     cycle fileloop
   end if

   call nc_check( nf90_close(ncid) , 'close file')
   
   obsloop: do k = 1, nlevels
   
     call interp_height_wght(hght, obs_levels(k), nobs, zloc, wght)
     if ( zloc < 1 )  cycle obsloop
     hghtp(nlevels-k+1) = obs_levels(k) * 1000.0_r8
   
   end do obsloop
   
   first_obs = .true.
   
   obsloop2: do k = 1, nlevels
   
     call interp_height_wght(hght, obs_levels(k), nobs, zloc, wght)
     if ( zloc < 1 )  cycle obsloop2
   
     lato  = wght * lat(zloc)  + (1.0_r8 - wght) * lat(zloc+1)
     lono  = wght * lon(zloc)  + (1.0_r8 - wght) * lon(zloc+1)
     if ( lono < 0.0_r8 )  lono = lono + 360.0_r8
     hghto = wght * hght(zloc) + (1.0_r8 - wght) * hght(zloc+1)
     hghto = hghto * 1000.0_r8 ![hghto(km)=hghto(m)x1000]
     eleco = wght * elec(zloc) + (1.0_r8 - wght) * elec(zloc+1)
   
        obsval = eleco
        !Real observation experiment
        !oerr   = ionprf_obserr_percent(lono,lato,hghto,ihour,imin,f3coerr)
        !oerr  = 0.01_r8 * ionprf_obserr_percent(lono,lato,hghto,ihour,imin,f3coerr) * obsval
        !oerr = 0.01_r8 * obsval
        !oerr=abs(oerr)
         
        !OSSE condition
        !oerr   = 1400.0 !This value is apply for linear scale only.
         oerr = 30000 - abs(hghto*0.105-30000)   
     
     subset = 'GPSREF'
     call set_obs_def_location(obs_def,set_location(lono,lato,hghto,3))
     call set_obs_def_kind(obs_def, GPS_PROFILE)
     if (overwrite_time) then   ! generally do not want to do this here.
        call set_obs_def_time(obs_def, set_time(asec, aday))
     else
        call set_obs_def_time(obs_def, set_time(osec, oday))
     endif
     call set_obs_def_error_variance(obs_def, oerr*oerr )
     call set_obs_def_key(obs_def, obs_num)
     call set_obs_def(obs, obs_def)
   
     obs_val(1) = obsval
     call set_obs_values(obs, obs_val)
     if (obsval < 1000.0) then
        qc=5.00000
     elseif (obsval > 100000000.0) then
        qc=5.00000
     elseif (err > 45) then
        qc=5.00000
     elseif (err <-25) then
        qc=5.00000
     else
        qc=0.00000
     endif
     qc_val(1)  = qc
     call set_qc(obs, qc_val)

     ! first one, insert with no prev.  otherwise, since all times will be the
     ! same for this column, insert with the prev obs as the starting point.
     ! (this code used to call append, but calling insert makes it work even if
     ! the input files are processed out of strict time order, which one message
     ! i got seemed to indicate was happening...)
     if (first_obs) then
        call insert_obs_in_seq(obs_seq, obs)
        first_obs = .false.
     else
        call insert_obs_in_seq(obs_seq, obs, prev_obs)
     endif
     obs_num = obs_num+1
     prev_obs = obs

     if (.not. did_obs) did_obs = .true.


   end do obsloop2

  ! clean up and loop if there is another input file
  deallocate( lat, lon, hght, elec )

  filenum = filenum + 1

end do fileloop

! done with main loop.  if we added any obs to the sequence, write it out.
if (did_obs) then
!print *, 'ready to write, nobs = ', get_num_obs(obs_seq)
   if (get_num_obs(obs_seq) > 0) &
      call write_obs_seq(obs_seq, gpsro_out_file)

   ! minor stab at cleanup, in the off chance this will someday get turned
   ! into a subroutine in a module.  probably not all that needs to be done,
   ! but a start.
!print *, 'calling destroy_obs'
   call destroy_obs(obs)
   call destroy_obs(prev_obs)
print *, 'skipping destroy_seq'
   ! get core dumps here, not sure why?
   !if (get_num_obs(obs_seq) > 0) call destroy_obs_sequence(obs_seq)
endif

! END OF MAIN ROUTINE

contains

! local subroutines/functions follow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   interp_height_wght - subroutine that finds the vertical levels
!                        closest to the desired height level and the
!                        weights to give to these levels to perform
!                        vertical interpolation.
!
!    hght  - height levels in column
!    level - height level to interpolate to
!    zgrid - index of lowest level for interpolation
!    wght  - weight to give to the lower level in interpolation
!    iz    - number of vertical levels
!
!     created June 2008, Ryan Torn NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interp_height_wght(hght, level, iz, zgrid, wght)

use        types_mod, only : r8

implicit none

integer,  intent(in)  :: iz
real(r8), intent(in)  :: hght(iz), level
integer,  intent(out) :: zgrid
real(r8), intent(out) :: wght

integer :: k, klev, kbot, ktop, kinc, kleva

if ( hght(1) > hght(iz) ) then
  kbot = iz  ;   ktop = 1   ;  kinc = -1   ;  kleva = 0
else
  kbot = 1   ;   ktop = iz  ;  kinc = 1    ;  kleva = 1
end if

if ( (hght(kbot) <= level) .AND. (hght(ktop) >= level) ) then

  do k = kbot, ktop, kinc  !  search for appropriate level
    if ( hght(k) > level ) then
      klev = k - kleva
      exit
    endif
  enddo

  ! compute the weights
  zgrid = klev
  wght  = (level-hght(klev+1)) / (hght(klev) - hght(klev+1))

else

  zgrid = -1
  wght  = 0.0_r8

endif

return
end subroutine interp_height_wght

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   ionprf_obserr_percent - function that computes the observation 
!                           error for a electron density observation.
!                           These numbers are taken from a Liu's and Yue's
!                           paper.
!
!    hghto - height of electron density observation
!    lono  - longitude of electron density observation
!    lato - latitude of electron density observation
!
!    created by I-TE LEE NCAR/HAO & NCU, 01/26/2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ionprf_obserr_percent(lone,late,hghte,houre,mine,oerr)

use   types_mod, only : r8

implicit none

!integer, parameter :: nobs_level = 22  !  maximum number of obs levels

real(r8), intent(in)  :: lone,late,hghte
integer , intent(in)  :: houre, mine 
real(r8)              :: mag_eq(73), slt, ionprf_obserr_percent, &
                         lonc, mq, nlat
real(r8)              :: err_top1, err_bottom1, err_obs1, &
                         err_top2, err_bottom2, err_obs2, err_obs
integer               :: altc, latc, ltc
real(r8), intent(in), dimension(15,37,25)   :: oerr


!Convert longitude to solar local time (simple method: longitude difference)
 slt= INT(NINT( (lone/15 + houre + mine/60) *10))
 if(slt>240) then 
     slt = slt - 240
 elseif(slt<0) then
     slt = slt + 240
 endif

!Calculate the magnetic equator latitude (refer to IGRF output in 2008)
data mag_eq/ 11.37_r8, 11.24_r8, 11.06_r8, 10.83_r8, 10.50_r8, 10.03_r8, &
              9.43_r8,  8.75_r8,  8.11_r8,  7.61_r8,  7.32_r8,  7.22_r8, &
              7.27_r8,  7.42_r8,  7.61_r8,  7.82_r8,  8.01_r8,  8.16_r8, &
              8.23_r8,  8.21_r8,  8.08_r8,  7.89_r8,  7.75_r8,  7.64_r8, & 
              7.64_r8,  7.72_r8,  7.81_r8,  7.85_r8,  7.77_r8,  7.55_r8, &
              7.18_r8,  6.64_r8,  5.93_r8,  5.07_r8,  4.12_r8,  3.16_r8, &
              2.26_r8,  1.47_r8,  0.80_r8,  0.23_r8, -0.28_r8, -0.78_r8, &
             -1.30_r8, -1.86_r8, -2.46_r8, -3.08_r8, -3.69_r8, -4.31_r8, &
             -4.97_r8, -5.71_r8, -6.56_r8, -7.52_r8, -8.58_r8, -9.69_r8, &
             -10.78_r8, -11.71_r8, -12.34_r8, -12.46_r8, -11.90_r8, -10.61_r8, &
             -8.63_r8, -6.08_r8, -3.17_r8, -0.14_r8,  2.73_r8,  5.23_r8, &
              7.28_r8,  8.85_r8,  9.98_r8, 10.75_r8, 11.20_r8, 11.38_r8, &
             11.38_r8/ 

 lonc=INT(lone/5+1)
 mq=mag_eq(lonc) + (mag_eq(lonc+1) - mag_eq(lonc)) / 5 * (lone+5-5*lonc)
 nlat=late-mq
 if(nlat > 90) then
    nlat=90
 elseif(nlat <-90) then
    nlat=-90
 endif

!Linear interpolation the error percentage
 latc = INT((nlat+90)/5+1)
 altc = INT((hghte-100000)/50000+1)
  ltc = INT(NINT(slt/10)+1)

! err_top    = oerr(altc+1,latc,ltc)+(oerr(altc+1,latc+1,ltc)-oerr(altc,latc,ltc))/5*(nlat+95-latc*5)
! err_bottom = oerr(altc,latc,ltc)+(oerr(altc,latc+1,ltc)-oerr(altc+1,latc,ltc))/5*(nlat+95-latc*5)
! err_obs    = err_bottom+ (err_top-err_bottom)/50000*(hghte-50000-altc*50000)

!!3-Dimensioonal Linear interpolation of the observation error percentage
err_top1    = oerr(altc+1,latc,ltc)+(oerr(altc+1,latc+1,ltc)-oerr(altc,latc,ltc))/5*(nlat+95-latc*5)
err_bottom1 = oerr(altc,latc,ltc)+(oerr(altc,latc+1,ltc)-oerr(altc+1,latc,ltc))/5*(nlat+95-latc*5)
err_obs1    = err_bottom1 + (err_top1-err_bottom1)/50000*(hghte-50000-altc*50000)

err_top2    = oerr(altc+1,latc,ltc+1)+(oerr(altc+1,latc+1,ltc+1)-oerr(altc,latc,ltc+1))/5*(nlat+95-latc*5)
err_bottom2 = oerr(altc,latc,ltc+1)+(oerr(altc,latc+1,ltc+1)-oerr(altc+1,latc,ltc+1))/5*(nlat+95-latc*5)
err_obs2    = err_bottom2 + (err_top2-err_bottom2)/50000*(hghte-50000-altc*50000)

err_obs     = 10.0_r8 + err_obs1 + (err_obs2-err_obs1)/10 * (slt+10-ltc*10)


ionprf_obserr_percent = err_obs

return
end function ionprf_obserr_percent

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
