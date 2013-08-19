! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program modis_subset_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   modis_subset_to_obs - 
!
!     created 3 May 2012   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), &
                              operator(>), operator(<), &
                              operator(==), operator(/=), operator(<=), operator(>=)

use      location_mod, only : VERTISSURFACE

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : MODIS_LEAF_AREA_INDEX, &
                              MODIS_FPAR

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=128) :: text_input_file = 'MOD15A2.fn_usbouldr.txt'
character(len=128) :: metadata_file   = 'MODIS_subset_metadata.txt'
character(len=128) :: obs_out_file    = 'obs_seq.out'
real(r8)           :: maxgoodqc       = 3.0_r8
logical            :: verbose         = .false.

namelist /modis_subset_to_obs_nml/ text_input_file, metadata_file, obs_out_file, maxgoodqc, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=600)      :: input_line, bigline
character(len=256)      :: string1, string2, string3
integer                 :: iline
logical                 :: first_obs
integer                 :: oday, osec, rcio, iunit
integer                 :: num_copies, num_qc, max_obs
real(r8)                :: oerr, qc, obval
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time
real(r8), parameter     :: umol_to_gC = (1.0_r8/1000000.0_r8) * 12.0_r8

type modistype
  character(len=20) :: HDFnamestring       = 'HDFname'
  character(len=20) :: Productstring       = 'Product'
  character(len=20) :: Datestring          = 'Date'
  character(len=20) :: Sitestring          = 'Site'
  character(len=20) :: ProcessDatestring   = 'ProcessDate'
  character(len=20) :: Bandstring          = 'Band'
  character(len=20) :: Fparstring          = 'Fpar_1km'
  character(len=20) :: Laistring           = 'Lai_1km'
  character(len=20) :: FparStdDevstring    = 'FparStdDev_1km'
  character(len=20) :: LaiStdDevstring     = 'LaiStdDev_1km'
  character(len=20) :: FparLai_QCstring    = 'FparLai_QC'
  character(len=20) :: FparExtra_QCstring  = 'FparExtra_QC'
  integer  :: HDFnameindex
  integer  :: Productindex
  integer  :: Dateindex
  integer  :: Siteindex
  integer  :: ProcessDateindex
  integer  :: Bandindex
  integer  :: Fparindex
  integer  :: FparStdDevindex
  integer  :: FparLai_QCindex
  integer  :: FparExtra_QCindex
  integer  :: Lai_index
  integer  :: LaiStdDev_index
  character(len=80) :: HDFname
  character(len=80) :: Product
  character(len=80) :: Date
  character(len=80) :: Site
  character(len=80) :: ProcessDate
  character(len=80) :: Band
  real(r8) :: Fpar
  real(r8) :: Lai
  real(r8) :: FparStdDev
  real(r8) :: FparLai_QC
  real(r8) :: FparExtra_QC
  real(r8) :: LaiStdDev
  
  type(time_type)   :: time_obs
  type(time_type)   :: Fpar_time
  type(time_type)   :: Lai_time
  type(time_type)   :: FparStdDev_time
  type(time_type)   :: FparLai_QC_time
  type(time_type)   :: FparExtra_QC_time
  type(time_type)   :: LaiStdDev_time

  character(len=40) :: Fpar_site
  character(len=40) :: Lai_site
  character(len=40) :: FparStdDev_site
  character(len=40) :: FparLai_QC_site
  character(len=40) :: FparExtra_QC_site
  character(len=40) :: LaiStdDev_site
end type modistype

type(modistype) :: modis

type metadata_type
   integer :: site_name_indx      ! Match(columns,'SITE_NAME')
   integer :: country_indx        ! Match(columns,'COUNTRY')
   integer :: latitude_indx       ! Match(columns,'LATITUDE')
   integer :: longitude_indx      ! Match(columns,'LONGITUDE')
   integer :: modis_site_id_indx  ! Match(columns,'MODIS_SITE_ID')
   integer :: fluxnet_id_indx     ! Match(columns,'FLUXNET_ID')
   integer :: fluxnet_key_id_indx ! Match(columns,'FLUXNET_KEY_ID')
   integer :: nstations = 0
   character(len=80), allocatable, dimension(:) :: site_name
   character(len=80), allocatable, dimension(:) :: country
   real(r8),          allocatable, dimension(:) :: latitude
   real(r8),          allocatable, dimension(:) :: longitude
   character(len=80), allocatable, dimension(:) :: modis_site_id
   character(len=80), allocatable, dimension(:) :: fluxnet_id
   character(len=80), allocatable, dimension(:) :: fluxnet_key_id
end type metadata_type

type(metadata_type) :: metadata

integer           :: siteindex
real(r8)          :: latitude
real(r8)          :: longitude
integer           :: landcover
character(len=80) :: IDstring
integer, dimension(49) :: values

type(time_type)   :: current_time

integer :: bob, bits567

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('modis_subset_to_obs')

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "modis_subset_to_obs_nml", iunit)
read(iunit, nml = modis_subset_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "modis_subset_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=modis_subset_to_obs_nml)
if (do_nml_term()) write(     *     , nml=modis_subset_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
prev_time = set_time(0, 0)

! Check to make sure the metadata filename exists.

call Get_MODIS_subset_metadata()

! We need to know the maximum number of observations in the input file.
! There are multiple lines for each observation, so file line count
! is more than enough. This program will only write out the actual number 
! created.
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(text_input_file)

max_obs    = count_file_lines(iunit)
num_copies = 1
num_qc     = 1
first_obs  = .true.

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(  obs_seq, 1, 'MODIS QC')

! The first line describes all the fields ... column headers, if you will

rewind(iunit)
call decode_header(iunit)

obsloop: do iline = 2,max_obs

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) bigline
   if (rcio < 0) exit obsloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline, trim(text_input_file)
      call error_handler(E_ERR,'main', string1, source, revision, revdate)
   endif

   input_line = adjustl(bigline)

   ! parse line into a temporary array ..

   ! parse the line into the modis structure (including the observation time)
   call stringparse(input_line, iline, IDstring, values)

   if (iline <= 2) then
      write(*,*)''
      write(*,*)'Check of the first observation: (column,string,value)'
      write(*,*)modis%HDFnameindex,     modis%HDFnamestring,     modis%HDFname
      write(*,*)modis%Productindex,     modis%Productstring,     modis%Product
      write(*,*)modis%Dateindex,        modis%Datestring,        modis%Date
      write(*,*)modis%Siteindex,        modis%Sitestring,        modis%Site
      write(*,*)modis%ProcessDateindex, modis%ProcessDatestring, modis%ProcessDate
      write(*,*)modis%Bandindex,        modis%Bandstring,        modis%Band

      call print_date(modis%time_obs, 'observation date is')
      call print_time(modis%time_obs, 'observation time is')
      current_time = modis%time_obs
   end if

   if (verbose) call print_date(modis%time_obs, 'obs time is')

   siteindex = get_metadata(modis%Site)

   ! extract the site from the filename and determine metadata 

   latitude  = metadata%latitude( siteindex)
   longitude = metadata%longitude(siteindex)
   landcover = -23

   if (verbose) print *, 'modis located at lat, lon =', latitude, longitude

   ! check the lat/lon values to see if they are ok
   if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

   if (( latitude > 90.0_r8 .or. latitude  <  -90.0_r8 ) .or. &
       (longitude <  0.0_r8 .or. longitude >  360.0_r8 )) then

      write (string2,*)'latitude  should be [-90, 90] but is ',latitude
      write (string3,*)'longitude should be [  0,360] but is ',longitude

      string1 ='modis location error in input.nml&modis_subset_to_obs_nml'
      call error_handler(E_ERR,'modis_subset_to_obs', string1, source, revision, &
                      revdate, text2=string2,text3=string3)

   endif

   ! make an obs derived type, and then add it to the sequence
   ! If the QC value is good, use the observation.
   ! Increasingly larger QC values are more questionable quality data.
   ! The observation errors are from page 183, Table 7.1(A) in 
   ! Chapter 7 of a book by A.D. Richardson et al. via Andy Fox.  

!  if (modis%FparExtra_QC <= maxgoodqc) then   ! Sensible Heat Flux [W m-2]

   select case ( IDstring )
   case ( 'Fpar_1km' )
     modis%Fpar      = check_modis_qc(values(25),'Fpar_1km')
     modis%Fpar_time = modis%time_obs

   case ( 'Lai_1km' )
     modis%Lai      = check_modis_qc(values(25),'Lai_1km')
     modis%Lai_time = modis%time_obs

   case ( 'FparExtra_QC' ) ! This is of no value to us at this time.
     modis%FparExtra_QC      = values(25)
     modis%FparExtra_QC_time = modis%time_obs

   case ( 'FparLai_QC' )
     modis%FparLai_QC      = check_modis_qc(values(25),'FparLai_QC')
     modis%FparLai_QC_time = modis%time_obs

   case ( 'FparStdDev_1km' )
     modis%FparStdDev      = check_modis_qc(values(25),'FparStdDev_1km')
     modis%FparStdDev_time = modis%time_obs

   case ( 'LaiStdDev_1km' )
     modis%LaiStdDev      = check_modis_qc(values(25),'LaiStdDev_1km')
     modis%LaiStdDev_time = modis%time_obs

   case default
   end select

   ! FIXME ... after observation is confirmed new time ... put that data into
   ! new array ... reinitialize all fields ... etc.

   if (modis%time_obs /= current_time) then
      ! there is a new time ... must be finished reading this obs 
   !  if ( consistent_time( modis%Fpar_time, modis%Lai_time,  &
   !          modis%FparExtra_QC_time, modis%FparLai_QC_time, &
   !          modis%FparStdDev_time, modis%LaiStdDev_time) )

      if ( .true. ) then

         ! QC flags are binary-coded ascii strings  10011101
         ! bits 5,6,7 (the last three) are decoded as follows:
         ! 000 ... Main(RT) method used, best result possible (no saturation)
         ! 001 ... Main(RT) method used with saturation, Good, very usable
         ! 010 ... Main(RT) method failed due to bad geometry, empirical algorithm used
         ! 011 ... Main(RT) method failed due to other problems
         ! 100 ... pixel not produced at all

         ! relying on integer arithmetic
         bob = 1000 * (modis%FparLai_QC / 1000)
         bits567 = modis%FparLai_QC - bob

         write(*,*)'modis%FparLai_QC and last 3',modis%FparLai_QC, bits567

         if (bits567 > 1) then
            modis%Fpar = MISSING_R8
            modis%Lai  = MISSING_R8
         endif

         if (modis%FparLai_QC <= maxgoodqc) then

            call get_time(modis%time_obs, osec, oday)

            oerr = modis%FparStdDev
            qc   = real(modis%FparLai_QC,r8)

            obval= real(modis%Fpar,r8)

            call create_3d_obs(latitude, longitude, 0.0_r8, VERTISSURFACE, obval, &
                               MODIS_FPAR, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, modis%time_obs, prev_obs, prev_time, first_obs)

         endif

      endif

      ! Call finalize_obs
      current_time = modis%time_obs
   endif


end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.
!
!       NOTE: assumes the code is using the threed_sphere locations module,
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error (in units of standard deviation)
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs)
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: okind, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation (will also
!                be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)

use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
use time_manager_mod, only : time_type, operator(>=)

type(obs_sequence_type), intent(inout) :: seq
type(obs_type),          intent(inout) :: obs, prev_obs
type(time_type),         intent(in)    :: obs_time
type(time_type),         intent(inout) :: prev_time
logical,                 intent(inout) :: first_obs

! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.

if(first_obs) then    ! for the first observation, no prev_obs
   call insert_obs_in_seq(seq, obs)
   first_obs = .false.
else
   if(obs_time >= prev_time) then  ! same time or later than previous obs
      call insert_obs_in_seq(seq, obs, prev_obs)
   else                            ! earlier, search from start of seq
      call insert_obs_in_seq(seq, obs)
   endif
endif

! update for next time
prev_obs  = obs
prev_time = obs_time

end subroutine add_obs_to_seq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   count_file_lines --
!           count the lines in a text file.
!           rewinds the unit after counting.
!
!     iunit - handle to the already-open text file
!
!     created May 2, 2012   Tim Hoar, NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function count_file_lines(iunit)

integer, intent(in) :: iunit
integer :: count_file_lines

integer :: i
character(len=128) :: oneline

integer, parameter :: tenmillion = 10000000
rewind(iunit)

count_file_lines = 0
countloop : do i = 1,tenmillion

   read(iunit,'(A)',iostat=rcio) oneline

   if (rcio < 0) exit countloop ! end of file
   if (rcio > 0) then
      write (string1,'('' read around line '',i8)')i
      call error_handler(E_ERR,'count_file_lines', string1, &
                         source, revision, revdate)
   endif
   count_file_lines = count_file_lines + 1

enddo countloop
rewind(iunit)

if (count_file_lines >= tenmillion) then
   write (string1,'('' suspiciously large number of lines '',i8)')count_file_lines
   call error_handler(E_MSG,'count_file_lines', string1, &
                         source, revision, revdate)
endif

end function count_file_lines




subroutine decode_header(iunit)
! Reads the first line of the header and parses the information.
! And by parse, I mean determine which columns are the columns
! of interest.

integer, intent(in) :: iunit

integer, parameter :: maxwordlength = 80
integer :: i,charcount,columncount,wordlength,maxlength
character(len=maxwordlength), dimension(:), allocatable :: columns
integer, dimension(10) :: qc = 0

! Read the line and strip off any leading whitespace.

read(iunit,'(A)',iostat=rcio) bigline
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(bigline(1:40)),'>'
  call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
endif

input_line = adjustl(bigline)

! Count how many commas are in the line - use this to determine how many columns

charcount = CountChar(input_line,',')
columncount = charcount + 1
allocate(columns(columncount))

columncount  = 0  ! track the number of columns
wordlength   = 0  ! number of characters in the column descriptor
charcount    = 0  ! the position of the (last) comma
do i = 1,len_trim(input_line)
   if (input_line(i:i) == ',') then
      columncount = columncount + 1
      if (wordlength > maxwordlength) then
         write(string1,*)'unexpected long word ... starts <',&
                           input_line((i-wordlength):(i-1)),'>'
         call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
      endif
      columns(columncount) = input_line((i-wordlength):(i-1)) 
      if (verbose) write(*,*)'word(',columncount,') is ',columns(columncount)
      wordlength = 0
      charcount = i
   else
      wordlength = wordlength + 1
   endif
enddo

! There is one more column after the last comma

columns(columncount+1) = input_line((charcount+1):len_trim(input_line))

! Finally get to the task at hand

modis%HDFnameindex     = Match(columns,modis%HDFnamestring)     ! used to be  1
modis%Productindex     = Match(columns,modis%Productstring)     ! used to be  2
modis%Dateindex        = Match(columns,modis%Datestring)        ! used to be  3
modis%Siteindex        = Match(columns,modis%Sitestring)        ! used to be  4
modis%ProcessDateindex = Match(columns,modis%ProcessDatestring) ! used to be  5
modis%Bandindex        = Match(columns,modis%Bandstring)        ! used to be  6

! Check to make sure we got all the indices we need

qc( 1) = CheckIndex( modis%HDFnameindex,     'HDFname' )
qc( 2) = CheckIndex( modis%Productindex,     'Product' )
qc( 3) = CheckIndex( modis%Dateindex,        'Date' )
qc( 4) = CheckIndex( modis%Siteindex,        'Site' )
qc( 9) = CheckIndex( modis%ProcessDateindex, 'ProcessDate' )
qc(10) = CheckIndex( modis%Bandindex,        'Band' )

if (any(qc /= 0) ) then
  write(string1,*)'Did not find all the required column indices.'
  call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
endif

! Summarize if desired

if (verbose) then
   write(*,*)'index is ', modis%HDFnameindex,    ' at one point it was  1'
   write(*,*)'index is ', modis%Productindex,    ' at one point it was  2'
   write(*,*)'index is ', modis%Dateindex,       ' at one point it was  3'
   write(*,*)'index is ', modis%Siteindex,       ' at one point it was  4'
   write(*,*)'index is ', modis%ProcessDateindex,' at one point it was  5'
   write(*,*)'index is ', modis%Bandindex,       ' at one point it was  6'
endif

deallocate(columns)

end subroutine decode_header



function CountChar(str1,solo)
! Count the number of instances of the single character in a character string.
! useful when parsing a comma-separated list, for example.
! Count the commas and add 1 to get the number of items in the list.

integer                      :: CountChar
character(len=*), intent(in) :: str1
character,        intent(in) :: solo

integer :: i

CountChar = 0
do i = 1,len_trim(str1)
   if (str1(i:i) == solo) CountChar = CountChar + 1
enddo

end function CountChar



function Match(sentence,word)
! Determine the first occurrence of the 'word' in a sentence.
! In this context, a sentence is a character array, the dimension
! of the array is the number of words in the sentence.
! This is a case-sensitive match. Trailing blanks are removed.

integer :: Match
character(len=*), dimension(:), intent(in) :: sentence
character(len=*),               intent(in) :: word

integer :: i

Match = 0
WordLoop : do i = 1,len(sentence)
   if (trim(sentence(i)) == trim(word)) then
      Match = i
      return
   endif
enddo WordLoop

end function Match



function CheckIndex( myindex, context )
! Routine to issue a warning if the index was not found.
! Returns an error code ... 0 means the index WAS found
! a negative number means the index was NOT found - an error condition.
! I want to check ALL the indexes before fatally ending.

integer                       :: CheckIndex
integer,          intent(in)  :: myindex
character(len=*), intent(in)  :: context

if (myindex == 0) then
   write(string1,*)'Did not find column header matching ',trim(context)
   call error_handler(E_MSG,'decode_header:CheckIndex',string1, source, revision, revdate)
   CheckIndex = -1 ! not a good thing
else
   CheckIndex = 0  ! Good to go
endif

end function CheckIndex



subroutine stringparse(str1, linenum, ObsString, quantities )
! just declare everything as reals and chunk it

character(len=*), intent(in) :: str1
integer         , intent(in) :: linenum
character(len=80), intent(out) :: ObsString
integer, dimension(49), intent(out) :: quantities

character(len=80), dimension(6) :: column_names
integer :: iyear, idoy
type(time_type) :: time0, time1, time2

quantities(:) = MISSING_R8

read(str1,*,iostat=rcio) column_names, quantities

if (rcio /= 0) then
   write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(str1(1:40)),'>'
   call error_handler(E_ERR,'stringparse',string1, source, revision, revdate)
endif

! Stuff what we want into the modis structure
!
! Convert to 'CLM-friendly' units AFTER we determine observation error variance.
! That happens in the main routine.
!
! NEE_or_fMDS    has units     [umolCO2 m-2 s-1] 
! H_f            has units     [W m-2]
! LE_f           has units     [W m-2]
!
! (CLM) NEE      has units     [gC m-2 s-1]

modis%HDFname     = trim(column_names(modis%HDFnameindex))
modis%Product     = trim(column_names(modis%Productindex))
modis%Date        = trim(column_names(modis%Dateindex))
modis%Site        = trim(column_names(modis%Siteindex))
modis%ProcessDate = trim(column_names(modis%ProcessDateindex))
modis%Band        = trim(column_names(modis%Bandindex))

read(modis%Date,'(1x,i4,i3)')iyear,idoy

time0   = set_date(iyear, 1, 1, 0, 0, 0)
time1   = set_time(0, idoy)
modis%time_obs = time0 + time1

write(*,*)'year is ',iyear
write(*,*)'idoy is ',idoy
call print_time(modis%time_obs)
call print_date(modis%time_obs)

IDstring = modis%Band

end subroutine stringparse



subroutine Get_MODIS_subset_metadata()
! Reads the list of station names and locations

character(len=80), dimension(7) :: columns

!SITE_NAME,COUNTRY,LATITUDE,LONGITUDE,MODIS_SITE_ID,FLUXNET_ID,FLUXNET_KEY_ID
! ,,,,,,
!Tomakomai_National_Forest,Japan,42.73697,141.51864,fn_jptomako,JP-Tom,jp.tomakomai.01

iunit = open_file(metadata_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(metadata_file)

metadata%nstations = count_file_lines(iunit) - 2

if (verbose) print *, 'there are ',metadata%nstations,' lines in the file.'

allocate(metadata%site_name(     metadata%nstations))
allocate(metadata%country(       metadata%nstations))
allocate(metadata%latitude(      metadata%nstations))
allocate(metadata%longitude(     metadata%nstations))
allocate(metadata%modis_site_id( metadata%nstations))
allocate(metadata%fluxnet_id(    metadata%nstations))
allocate(metadata%fluxnet_key_id(metadata%nstations))

rewind(iunit)
call decode_metadata(iunit)

! Skip the first two lines 
rewind(iunit)
read(iunit,'(A)',iostat=rcio) bigline
read(iunit,'(A)',iostat=rcio) bigline

readloop: do iline = 1,metadata%nstations

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) bigline
   if (rcio < 0) exit readloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline+2, trim(text_input_file)
      call error_handler(E_ERR,'Get_MODIS_subset_metadata', string1, source, revision, revdate)
   endif

   input_line = adjustl(bigline)

   read(input_line,*,iostat=rcio) columns
   if (rcio /= 0) then
      write(string1,*)'Cannot parse line',iline+2,'. Begins <',trim(input_line(1:40)),'>'
      call error_handler(E_ERR,'Get_MODIS_subset_metadata',string1, source, revision, revdate)
   endif

   metadata%site_name(iline)      = trim(columns(metadata%site_name_indx))
   metadata%country(iline)        = trim(columns(metadata%country_indx))
   metadata%modis_site_id(iline)  = trim(columns(metadata%modis_site_id_indx))
   metadata%fluxnet_id(iline)     = trim(columns(metadata%fluxnet_id_indx))
   metadata%fluxnet_key_id(iline) = trim(columns(metadata%fluxnet_key_id_indx))

   read(columns(metadata%latitude_indx ),*,iostat=rcio) metadata%latitude( iline)
   if (rcio/=0) metadata%latitude(iline) = MISSING_R8

   read(columns(metadata%longitude_indx),*,iostat=rcio) metadata%longitude(iline)
   if (rcio/=0) metadata%longitude(iline) = MISSING_R8

   if (verbose) then
      write(*,*)'line number ',iline
      write(*,*)'   site_name      is ', metadata%site_name(iline)
      write(*,*)'   country        is ', metadata%country(iline)
      write(*,*)'   latitude       is ', metadata%latitude( iline)
      write(*,*)'   longitude      is ', metadata%longitude(iline)
      write(*,*)'   modis_site_id  is ', metadata%modis_site_id(iline)
      write(*,*)'   fluxnet_id     is ', metadata%fluxnet_id(iline)
      write(*,*)'   fluxnet_key_id is ', metadata%fluxnet_key_id(iline)
      write(*,*)
   endif

enddo readloop

end subroutine Get_MODIS_subset_metadata



subroutine decode_metadata(iunit)
! Reads the first line of the header and parses the information.
! And by parse, I mean determine which columns are the columns
! of interest.

integer,              intent(in)  :: iunit

integer, parameter :: maxwordlength = 80
integer :: i,charcount,columncount,wordlength,maxlength
character(len=maxwordlength), dimension(:), allocatable :: columns
integer, dimension(10) :: qc = 0

! Read the line and strip off any leading whitespace.

read(iunit,'(A)',iostat=rcio) bigline
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(bigline(1:40)),'>'
  call error_handler(E_ERR,'decode_metadata',string1, source, revision, revdate)
endif

input_line = adjustl(bigline)

! Count how many commas are in the line - use this to determine how many columns

charcount = CountChar(input_line,',')
columncount = charcount + 1
allocate(columns(columncount))

columncount  = 0  ! track the number of columns
wordlength   = 0  ! number of characters in the column descriptor
charcount    = 0  ! the position of the (last) comma
do i = 1,len_trim(input_line)
   if (input_line(i:i) == ',') then
      columncount = columncount + 1
      if (wordlength > maxwordlength) then
         write(string1,*)'unexpected long word ... starts <',&
                           input_line((i-wordlength):(i-1)),'>'
         call error_handler(E_ERR,'decode_metadata',string1, source, revision, revdate)
      endif
      columns(columncount) = input_line((i-wordlength):(i-1)) 
      if (verbose) write(*,*)'word(',columncount,') is ',columns(columncount)
      wordlength = 0
      charcount = i
   else
      wordlength = wordlength + 1
   endif
enddo

! There is one more column after the last comma

columns(columncount+1) = input_line((charcount+1):len_trim(input_line))

metadata%site_name_indx      = Match(columns,'SITE_NAME')
metadata%country_indx        = Match(columns,'COUNTRY')
metadata%latitude_indx       = Match(columns,'LATITUDE')
metadata%longitude_indx      = Match(columns,'LONGITUDE')
metadata%modis_site_id_indx  = Match(columns,'MODIS_SITE_ID')
metadata%fluxnet_id_indx     = Match(columns,'FLUXNET_ID')
metadata%fluxnet_key_id_indx = Match(columns,'FLUXNET_KEY_ID')

! Check to make sure we got all the indices we need

qc( 1) = CheckIndex( metadata%site_name_indx,     'SITE_NAME' )
qc( 2) = CheckIndex( metadata%country_indx,       'COUNTRY' )
qc( 3) = CheckIndex( metadata%latitude_indx,      'LATITUDE' )
qc( 4) = CheckIndex( metadata%longitude_indx,     'LONGITUDE' )
qc( 5) = CheckIndex( metadata%modis_site_id_indx, 'MODIS_SITE_ID' )
qc( 6) = CheckIndex( metadata%fluxnet_id_indx,    'FLUXNET_ID' )
qc( 7) = CheckIndex( metadata%fluxnet_key_id_indx,'FLUXNET_KEY_ID' )

if (any(qc /= 0) ) then
  write(string1,*)'Did not find all the required column indices.'
  call error_handler(E_ERR,'decode_metadata',string1, source, revision, revdate)
endif

deallocate(columns)

end subroutine decode_metadata




function   get_metadata(observation_site)
character(len=*), intent(in) :: observation_site
integer :: get_metadata

integer :: istation
character(len=40) :: x,y

get_metadata = -1

StationLoop: do istation = 1,metadata%nstations

   x = trim(metadata%modis_site_id(istation))
   y = trim(observation_site)

!  write(*,*)'checking ',istation,' x is .',trim(x),'.  y is .',trim(y),'.'

   if (x == y) then
      get_metadata = istation
      exit StationLoop
   endif

enddo StationLoop

end function get_metadata


function check_modis_qc(value,valstring)
! Checks the appropriate QC values for the value in question
! and converts 'good' values to proper units.
   integer, intent(in) :: value
   character(len=*), intent(in) :: valstring
   real(r8) :: check_modis_qc

   select case ( valstring )
   case ('Fpar_1km') 
      if (value < 249) check_modis_qc = value * 0.01_r8

   case ('FparStdDev_1km') 
      if (value < 248) check_modis_qc = value * 0.01_r8

   case ('Lai_1km') 
      if (value < 249) check_modis_qc = value * 0.1_r8

   case ('LaiStdDev_1km') 
      if (value < 248) check_modis_qc = value * 0.1_r8
   case default
   end select

end function check_modis_qc



end program modis_subset_to_obs


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
