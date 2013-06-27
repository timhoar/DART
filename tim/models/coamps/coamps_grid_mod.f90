! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module coamps_grid_mod

!------------------------------
! MODULE:       coamps_grid_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the definition of a COAMPS grid data structure 
! and routines dealing with COAMPS grid information.
!------------------------------ 

  use coamps_intrinsic_mod, only : ij2ll,                &
                                   ll2ij,                &
                                   uvg2uv

  use coamps_util_mod,      only : C_REAL,               & 
                                   C_REAL4,              & 
                                   check_io_status,      &
                                   check_alloc_status,   &
                                   check_dealloc_status, &
                                   fix_for_platform4

  use location_mod,  only : get_location,  &
                            location_type, &
                            set_location

  use utilities_mod, only : do_output,     &
                            E_ERR,         &
                            E_WARN,        &
                            E_MSG,         &
                            error_handler, &
                            get_unit,      &
                            nc_check

  use coamps_flat_file_mod, only : read_flat_file, &
                                   generate_one_flat_file_name

  use types_mod,     only : r8

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! COAMPS grid data structure
  public :: coamps_grid
  public :: initialize_grid
  
  ! Location conversion functions
  public :: gridpt_to_latlon
  public :: latlon_to_gridpt
  public :: gridpt_to_location 
  public :: location_to_gridpt
  public :: grid_ij_to_vector_index

  ! Wind conversion functions
  public :: grid_wind_to_earth_wind

  ! Information about the grid size
  public :: get_grid_dims
  public :: get_grid_field_size
  public :: get_grid_num_levels
  public :: check_ij_within_grid

  ! Information about horizontal coordinates
  public :: get_grid_delta_x
  public :: get_grid_delta_y

  ! Information about vertical coordinates
  public :: get_grid_msigma
  public :: get_grid_wsigma
  public :: get_grid_dsigmaw
  public :: get_terrain_height_at_points

  ! Domain decomposition tools 
  public :: initialize_decomposition
  public :: decompose_domain
  public :: get_grid_iminf
  public :: get_grid_imaxf
  public :: get_grid_jminf
  public :: get_grid_jmaxf

  ! Debugging output
  public :: dump_grid_info
  public :: nc_dump_grid_info
  
  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  ! [none]
  !------------------------------
  ! END EXTERNAL INTERFACES
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS 
  !------------------------------

  ! Derived data type for the COAMPS grid information - this is all
  ! contained in the datahd_*_infofld files that the COAMPS analysis
  ! produces, so we can just read that rather than processing a
  ! namelist.  Most of these are pretty straighforward, but:
  !  anchor_i - the i-grid point where the grid is anchored to parent
  !  anchor_j - the j-grid point where the grid is anchored to parent
  !  iref     - the i-grid point where the grid is anchored to Earth
  !  jref     - the j-grid point where the grid is anchored to Earth
  !  dsigma   - the spacing between mass sigma levels
  !  msigma   - actual sigma values at mass levels
  !  terrain  - terrain height
  !  [ij]minf - the minimum grid point for process N
  !  [ij]maxf - the maximum grid point for process N
  ! NOTE: For now, "nests" is currently unused and the pts, anchor,
  ! and [ij]ref variables are scalars - when multiple grid nests are
  ! supported, these will be vectors.  Currently only a single grid
  ! is supported to make my life much, much, much easier.
  type :: coamps_grid
     real(kind=r8)                              :: pres_lvls
     real(kind=r8)                              :: sigm_lvls
     real(kind=r8)                              :: map_proj
     real(kind=r8)                              :: stdlat1
     real(kind=r8)                              :: stdlat2
     real(kind=r8)                              :: stdlon
     real(kind=r8)                              :: reflat
     real(kind=r8)                              :: reflon
     real(kind=r8)                              :: delta_x
     real(kind=r8)                              :: delta_y
     real(kind=r8)                              :: nests
     real(kind=r8)                              :: pts_x
     real(kind=r8)                              :: pts_y
     real(kind=r8)                              :: anchor_i
     real(kind=r8)                              :: anchor_j
     real(kind=r8)                              :: iref
     real(kind=r8)                              :: jref
     real(kind=r8)                              :: grid_constant
     real(kind=r8), dimension(:),   allocatable :: dsigma
     real(kind=r8), dimension(:),   allocatable :: dsigmaw
     real(kind=r8), dimension(:),   allocatable :: msigma
     real(kind=r8), dimension(:),   allocatable :: wsigma
     real(kind=r8), dimension(:,:), allocatable :: terrain

     ! These fields are used for dealing with domain decomposition
     integer, dimension(:), allocatable :: iminf
     integer, dimension(:), allocatable :: imaxf
     integer, dimension(:), allocatable :: jminf
     integer, dimension(:), allocatable :: jmaxf
  end type coamps_grid

  integer, parameter :: SINGLE_POINT = 1

  ! Model vertical coordinate indexes top-down
  integer, parameter :: TOP_LEVEL    = 1

  ! Location of entries in the datahd file
  integer, parameter :: DATAHD_NUM_P_LEVELS   = 1
  integer, parameter :: DATAHD_NUM_S_LEVELS   = 2
  integer, parameter :: DATAHD_MAP_PROJECTION = 3
  integer, parameter :: DATAHD_STANDARD_LAT_1 = 4
  integer, parameter :: DATAHD_STANDARD_LAT_2 = 5
  integer, parameter :: DATAHD_STANDARD_LON   = 6
  integer, parameter :: DATAHD_REFERENCE_LAT  = 7
  integer, parameter :: DATAHD_REFERENCE_LON  = 8
  integer, parameter :: DATAHD_DELTA_X        = 9
  integer, parameter :: DATAHD_DELTA_Y        = 10
  integer, parameter :: DATAHD_NUM_NESTS      = 11
  integer, parameter :: DATAHD_NUM_X_POINTS   = 30
  integer, parameter :: DATAHD_NUM_Y_POINTS   = 31
  integer, parameter :: DATAHD_ANCHOR_POINT_I = 32
  integer, parameter :: DATAHD_ANCHOR_POINT_J = 33
  integer, parameter :: DATAHD_REFERENCE_I    = 34
  integer, parameter :: DATAHD_REFERENCE_J    = 35
  integer, parameter :: DATAHD_DSIGMA_OFFSET  = 500
  integer, parameter :: DATAHD_MSIGMA_OFFSET  = 800

  ! POSSIBLE COAMPS PROJECTIONS
  integer, parameter :: PROJ_MERCATOR                = 1
  integer, parameter :: PROJ_LAMBERT_CONIC_CONFORMAL = 2
  integer, parameter :: PROJ_POLAR_STEREOGRAPHIC     = 3
  integer, parameter :: PROJ_NONE                    = 4

  !------------------------------
  ! END TYPES AND CONSTANTS 
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  logical :: module_initialized = .false.

  real(kind=r8), parameter :: pi         = 3.141592741012573_r8
  real(kind=r8), parameter :: DEG_TO_RAD = pi/180.0_r8

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------

  ! initialize_grid
  ! ---------------
  ! Populates a coamps_grid data structure using the COAMPS domain
  ! data (datahd) file and the COAMPS terrain height (terrht) file
  ! based on a given date-time group.
  ! N.B. "grid" parameter is marked as intent(inout) to ensure that
  !      grid will not be undefined on exit if we don't need to run
  !      this routine
  !  PARAMETERS
  !   IN  dtg               base date-time group for model run
  ! INOUT grid              coamps_grid structure to initialize
  subroutine initialize_grid(dtg, grid, dsnrff)
    character(len=*),  intent(in)    :: dtg
    type(coamps_grid), intent(inout) :: grid
    character(len=*),  intent(in)    :: dsnrff
    
    character(len=64) :: datahd_filename
    integer           :: datahd_unit
    character(len=64) :: terrht_filename
    integer           :: terrht_unit
    integer           :: terrht_record_len
    
    ! Terrain data is stored in a flat file, so it uses the COAMPS
    ! real number size.  Single dimension makes byteswapping easier
    real(kind=C_REAL), dimension(:), allocatable :: coamps_terrain

    integer :: field_size

    ! Error-checking info
    character(len=*), parameter :: routine = 'initialize_grid'
    integer :: io_status
    integer :: alloc_status
    integer :: dealloc_status

    if (module_initialized) return

    ! Initialize the domain data section of the coamps_grid structure
    call generate_datahd_filename(dtg, datahd_filename)
    datahd_unit = get_unit()
    print *,trim(dsnrff)//datahd_filename
    open(unit=datahd_unit, file=trim(dsnrff)//datahd_filename, &
         status='old', access='sequential', action='read',     &
         form='formatted', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Opening datahd file')
    call read_datahd_file(datahd_unit, grid)
    close(datahd_unit)
    
    ! Need to set up permanent and temporary storage for terrain data
    allocate( grid%terrain(int(grid%pts_x), int(grid%pts_y)), &
              stat = alloc_status )
    call check_alloc_status(alloc_status, routine, source,    &
                            revision, revdate, 'grid%terrain')

    call get_grid_field_size(grid, field_size)
    allocate (coamps_terrain(field_size), stat = alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'coamps_terrain (temporary)')

    ! Make sure when checking the record length to use what *COAMPS*
    ! wrote it as, not what we want to store it as
    call generate_terrht_filename(dtg, grid, terrht_filename)
    terrht_unit = get_unit()
    !inquire(iolength=terrht_record_len) coamps_terrain
    terrht_record_len=C_REAL4*field_size

    open(unit = terrht_unit, file =trim(dsnrff)//terrht_filename, &
         status = 'old', access = 'direct', action = 'read',      &
         form = 'unformatted',  recl = terrht_record_len,         &
         iostat = io_status) 
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Opening terrain height file')
!   call read_terrht_file(terrht_unit, coamps_terrain)
    call read_flat_file(terrht_unit, coamps_terrain)
    close(terrht_unit)

    ! COAMPS writes big-endian flat files no matter what platform is
    ! being used, so we may need to change the data we read.
    !call fix_for_platform(coamps_terrain, field_size, C_REAL)

    grid%terrain = reshape(real(coamps_terrain, kind=r8),        &
                           (/ int(grid%pts_x), int(grid%pts_y) /)) 

    deallocate(coamps_terrain, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'coamps_terrain')

    module_initialized = .true.
  end subroutine initialize_grid
  
  ! gridpt_to_latlon
  ! ----------------
  ! Given a COAMPS grid structure, convert an (i,j) representation
  ! of a point on that grid to a (lat, lon) representation.
  ! Note that this routine will happily accept and convert (i,j)
  ! points that are beyond the bounds of the domain.  For this 
  ! direction of conversion, that's not too big a problem.
  !  PARAMETERS
  !   IN  grid              coamps_grid to base conversion on
  !   IN  ir                point's x-direction index
  !   IN  jr                point's y-direction index
  !   OUT lat               point's latitude
  !   OUT lon               point's longitude
  subroutine gridpt_to_latlon(grid, ir, jr, lat, lon) 
    type(coamps_grid), intent(in)  :: grid
    real(kind=r8),     intent(in)  :: ir
    real(kind=r8),     intent(in)  :: jr
    real(kind=r8),     intent(out) :: lat
    real(kind=r8),     intent(out) :: lon

    ! Work-around for not being able to pass in scalars as arrays
    ! of length 1
    real(kind=r8), dimension(SINGLE_POINT) :: ir_array,  jr_array
    real(kind=r8), dimension(SINGLE_POINT) :: lat_array, lon_array

    ir_array(SINGLE_POINT) = ir
    jr_array(SINGLE_POINT) = jr

    call ij2ll(int(grid%map_proj), grid%reflat,    grid%reflon,  &
               int(grid%iref),     int(grid%jref), grid%stdlat1, &
               grid%stdlat2,       grid%stdlon,    grid%delta_x, &
               grid%delta_y,       ir_array,       jr_array,     &
               SINGLE_POINT,       lat_array,      lon_array)  

    lat = lat_array(SINGLE_POINT)
    lon = lon_array(SINGLE_POINT)
  end subroutine gridpt_to_latlon

  ! latlon_to_gridpt
  ! ----------------
  ! Given a COAMPS grid structure, convert a (lat,lon) 
  ! representation of a point on that grid to an (i, j) 
  ! representation.
  ! Note that this routine will happily accept and convert 
  ! (lat,lon) points that are beyond the bounds of the domain.
  ! This could happen in several ways, but especially if an obs
  ! gets included in the sequence that for some reason isn't in
  ! the domain - in this case, warn the user.
  !  PARAMETERS
  !   IN  grid              coamps_grid to base conversion on
  !   IN  lat               point's latitude
  !   IN  lon               point's longitude
  !   OUT ir                point's x-direction index
  !   OUT jr                point's y-direction index
  subroutine latlon_to_gridpt(grid, lat, lon, ir, jr) 
    type(coamps_grid), intent(in)  :: grid
    real(kind=r8),     intent(in)  :: lat
    real(kind=r8),     intent(in)  :: lon
    real(kind=r8),     intent(out) :: ir
    real(kind=r8),     intent(out) :: jr

!   logical :: is_within_grid

    ! Work-around for not being able to pass in scalars as arrays
    ! of length 1
    real(kind=r8), dimension(SINGLE_POINT) :: ir_array,  jr_array
    real(kind=r8), dimension(SINGLE_POINT) :: lat_array, lon_array

!   character(len=128) :: message
    
    ! If anything goes wrong in the conversion, pre-define values 
    ! that should set off alarm bells
    ir = -999.0_r8
    jr = -999.0_r8

    lat_array(SINGLE_POINT) = lat
    lon_array(SINGLE_POINT) = lon

    call ll2ij(int(grid%map_proj), grid%reflat,    grid%reflon,  &
               int(grid%iref),     int(grid%jref), grid%stdlat1, &
               grid%stdlat2,       grid%stdlon,    grid%delta_x, &
               grid%delta_y,       ir_array,       jr_array,     &
               SINGLE_POINT,       lat_array,      lon_array)  

    ir = ir_array(SINGLE_POINT)
    jr = jr_array(SINGLE_POINT)

!    call check_ij_within_grid(grid, ir, jr, is_within_grid)
!    if(.not.is_within_grid) then
!      write (message,'(2(A,1x,F7.3,1x),A)') 'LAT=',lat,'LON=',lon,'is outside domain.'
!      call error_handler(E_MSG, 'latlon_to_gridpt', trim(message),&
!                         source, revision, revdate)
!    end if
  end subroutine latlon_to_gridpt

  ! gridpt_to_location
  ! ------------------
  ! Given a COAMPS grid, an (i,j) horizontal location in that grid,
  ! and a vertical location and vertical coordinate type, returns
  ! a DART location structure for that point
  !  PARAMETERS
  !   IN  grid              coamps_grid for lat/lon conversion
  !   IN  ii                point's x-direction index
  !   IN  jj                point's y-direction index
  !   IN  vert_loc          point's vertical direction location
  !   IN  vert_coord        point's vertical direction coordinate
  !   OUT loc               point's location as a DART structure
  subroutine gridpt_to_location(grid, ii, jj, vert_loc, vert_coord, loc)
    type(coamps_grid),   intent(in)  :: grid
    integer,             intent(in)  :: ii
    integer,             intent(in)  :: jj
    real(kind=r8),       intent(in)  :: vert_loc
    integer,             intent(in)  :: vert_coord
    type(location_type), intent(out) :: loc
    
    real(kind=r8) :: real_ii, real_jj
    real(kind=r8) :: point_lat, point_lon 

    ! Get things squared away in 2D first
    real_ii = real(ii, kind=r8)
    real_jj = real(jj, kind=r8)
    call gridpt_to_latlon(grid, real_ii, real_jj, point_lat, &
                          point_lon)

    ! Do this with keyword arguments because I can't keep the order
    ! straight
    loc = set_location( lon        = point_lon, &
                        lat        = point_lat, &
                        vert_loc   = vert_loc,  &
                        which_vert = vert_coord )
  end subroutine gridpt_to_location

  ! location_to_gridpt
  ! ------------------
  ! Given a COAMPS grid and a DART location, return the location
  ! of a point as an (i, j) coordinate and a vertical coordinate.
  ! Note that DART's get_location function that's used to break
  ! this down does *not* return any information on what that
  ! vertical coordinate is, so no information is returned to the
  ! caller.
  !  PARAMETERS
  !   IN  grid              coamps_grid for lat/loon conversion
  !   IN  loc               DART location structure to unpack
  !   OUT ir                point's x-direction index
  !   OUT jr                point's y-direction index
  !   OUT vert_loc          point's vertical location
  subroutine location_to_gridpt(grid, loc, ir, jr, vert_loc)
    type(coamps_grid),   intent(in)  :: grid
    type(location_type), intent(in)  :: loc
    real(kind=r8),       intent(out) :: ir
    real(kind=r8),       intent(out) :: jr
    real(kind=r8),       intent(out) :: vert_loc

    ! The array returned by get_location is (lon, lat, vertical)
    real(kind=r8), dimension(3) :: loc_array

    integer, parameter :: DART_LOC_LON  = 1
    integer, parameter :: DART_LOC_LAT  = 2
    integer, parameter :: DART_LOC_VERT = 3
    
    loc_array = get_location(loc)

    ! Deal with the 2D component separately from the vertical
    ! Again, be explicit since I apparently like to mix these up
    call latlon_to_gridpt( grid = grid,                    &
                           lat  = loc_array(DART_LOC_LAT), &
                           lon  = loc_array(DART_LOC_LON), & 
                           ir   = ir,                      &
                           jr   = jr ) 
    vert_loc = loc_array(DART_LOC_VERT)
  end subroutine location_to_gridpt

  ! grid_ij_to_vector_index
  ! -----------------------
  ! Given a COAMPS grid, converts an horizontal (ii,jj) location 
  ! to an index in the vector representation
  !  PARAMETERS
  !   IN  grid              coamps_grid to use for the conversion
  !   IN  ii                zonal index component
  !   IN  jj                meridional index component
  !   OUT myindex           position that the point (ii,jj) has in
  !                         the field re-arranged in vector form
  subroutine grid_ij_to_vector_index(grid, ii, jj, myindex)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(in)  :: ii
    integer,           intent(in)  :: jj
    integer,           intent(out) :: myindex

    myindex = (jj - 1) * grid%pts_x + ii
  end subroutine grid_ij_to_vector_index

  ! get_grid_dims
  ! -------------
  ! Given a COAMPS grid structure, calculate the i (i.e. x) and 
  ! j (i.e. y) dimensions
  !  PARAMETERS
  !   IN  grid              coamps_grid to get limits from
  !   OUT grid_i            number of points on the i dimension
  !   OUT grid_j            number of points on the j dimension
  subroutine get_grid_dims(grid, grid_i, grid_j)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(out) :: grid_i
    integer,           intent(out) :: grid_j

    grid_i = int(grid%pts_x)
    grid_j = int(grid%pts_y)
  end subroutine get_grid_dims

  ! get_grid_field_size
  ! -------------------
  ! Given a coamps grid structure, calculate the field size (the
  ! total number of points in the i and j directions)
  !  PARAMETERS
  !   IN  grid              coamps_grid to get size from
  !   OUT field_size        number of grid points on a single level
  subroutine get_grid_field_size(grid, field_size)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(out) :: field_size

    field_size = int(grid%pts_x) * int(grid%pts_y)
  end subroutine get_grid_field_size

  ! get_grid_num_levels
  ! -------------------
  ! Returns the number of mass sigma levels on a given COAMPS grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to get levels from
  !   OUT num_sig_lvls      number of sigma levels on the grid
  subroutine get_grid_num_levels(grid, num_sig_lvls)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(out) :: num_sig_lvls
    
    num_sig_lvls = grid%sigm_lvls
  end subroutine get_grid_num_levels

  ! check_ij_within_grid
  ! --------------------
  ! Calculates if a given (i,j) index is within a given COAMPS grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to use for checking
  !   IN  ir                zonal index
  !   IN  jr                meridional index
  !   OUT is_within_grid    True if the (i,j) location listed is
  !                         inside the given grid
  subroutine check_ij_within_grid(grid, ir, jr, is_within_grid)
    type(coamps_grid), intent(in)  :: grid
    real(kind=r8),     intent(in)  :: ir
    real(kind=r8),     intent(in)  :: jr
    logical,           intent(out) :: is_within_grid

    is_within_grid = .true.
    if ( (ceiling(ir) < 1) .or. (floor(ir) > grid%pts_x) .or. &
         (ceiling(jr) < 1) .or. (floor(jr) > grid%pts_y) ) then
       is_within_grid = .false.
    end if
  end subroutine check_ij_within_grid

  ! get_grid_delta_x
  ! ----------------
  ! Returns the x-direction grid spacing (in meters) for the given grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to pull data from
  !   OUT grid_delta_x      x-direction grid spacing
  subroutine get_grid_delta_x(grid, grid_delta_x)
    type(coamps_grid), intent(in)  :: grid
    real(kind=r8),     intent(out) :: grid_delta_x

    grid_delta_x = grid%delta_x
  end subroutine get_grid_delta_x

  ! get_grid_delta_y
  ! ----------------
  ! Returns the y-direction grid spacing (in meters) for the given grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to pull data from
  !   OUT grid_delta_y      y-direction grid spacing
  subroutine get_grid_delta_y(grid, grid_delta_y)
    type(coamps_grid), intent(in)  :: grid
    real(kind=r8),     intent(out) :: grid_delta_y

    grid_delta_y = grid%delta_y
  end subroutine get_grid_delta_y

  ! get_grid_msigma
  ! ---------------
  ! Returns the mass sigma levels for a given COAMPS grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to pull msigma out of
  !   OUT msigma            mass sigma levels for the grid
  subroutine get_grid_msigma(grid, msigma)
    type(coamps_grid),           intent(in)  :: grid
    real(kind=r8), dimension(:), intent(out) :: msigma
    
    integer :: kka

    ! Don't attempt a copy if the data in the grid won't fit
    if ( size(grid%msigma) > size(msigma) ) then
       call error_handler(E_ERR, 'get_grid_msigma', 'msigma    &
                          &passed into function is too small', &
                          &source, revision, revdate)
    end if

    do kka=1,size(grid%msigma)
       msigma(kka) = grid%msigma(kka)
    end do
  end subroutine get_grid_msigma

  ! get_grid_wsigma
  ! ---------------
  ! Returns the w sigma levels for a given COAMPS grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to pull wsigma out of
  !   OUT wsigma            w sigma levels for the grid
  subroutine get_grid_wsigma(grid, wsigma)
    type(coamps_grid),           intent(in)  :: grid
    real(kind=r8), dimension(:), intent(out) :: wsigma
    
    integer :: kkw

    ! Bail if the data in the grid won't fit
    if ( size(grid%wsigma) > size(wsigma) ) then
       call error_handler(E_ERR, 'get_grid_wsigma', 'wsigma ' // &
                          'passed into function is too small',   &
                          source, revision, revdate)
    end if

    do kkw=1,size(grid%wsigma)
       wsigma(kkw) = grid%wsigma(kkw)
    end do
  end subroutine get_grid_wsigma

  ! get_grid_dsigmaw
  ! ----------------
  ! Returns the distance between w sigma levels for a given COAMPS
  ! grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to pull dsigmaw out of
  !   OUT dsigmaw           distance between w sigma levels
  subroutine get_grid_dsigmaw(grid, dsigmaw)
    type(coamps_grid),           intent(in)  :: grid
    real(kind=r8), dimension(:), intent(out) :: dsigmaw
    
    integer :: kkw

    ! Bail if the data in the grid won't fit
    if ( size(grid%dsigmaw) > size(dsigmaw) ) then
       call error_handler(E_ERR, 'get_grid_dsigmaw', 'dsigmaw ' // &
                          'passed into function is too small',     &
                          source, revision, revdate)
    end if

    do kkw=1,size(grid%dsigmaw)
       dsigmaw(kkw) = grid%dsigmaw(kkw)
    end do
  end subroutine get_grid_dsigmaw

  ! get_terrain_height_at_points
  ! ----------------------------
  ! Given a COAMPS grid and a vector representation of a set of 
  ! (i,j) points, returns the terrain height for those points
  !  PARAMETERS
  !   IN  grid              coamps_grid to pull terrht from
  !   IN  i_points          points's x-direction index
  !   IN  j_points          points's y-direction index
  !   OUT terrain_height    points's terrain height
  subroutine get_terrain_height_at_points(grid, i_points, j_points,&
                                          terrain_height)
    type(coamps_grid),           intent(in)  :: grid
    integer,       dimension(:), intent(in)  :: i_points
    integer,       dimension(:), intent(in)  :: j_points
    real(kind=r8), dimension(:), intent(out) :: terrain_height

    integer :: nn, num_points

    ! Sanity check
    if (size(i_points) .ne. size(j_points)) then
      call error_handler(E_ERR, 'get_terrain_height_at_points',   &
                         'Error in getting terrain height - ' //  &
                         'i/j point arrays are different sizes!', &
                         source, revision, revdate)
    end if

    num_points = size(i_points)
    do nn=1,num_points
       terrain_height(nn) = grid%terrain(i_points(nn), j_points(nn))
    end do
  end subroutine get_terrain_height_at_points

  ! initialize_decomposition
  ! ------------------------
  ! Given a COAMPS grid the total number of non-I/O processors for
  ! that grid, allocate space for each processor's [ij](min|max)f 
  ! extent
  ! N.B. The grid is not defined as intent(out) here since that 
  !      could mean that the other components of the coamps_grid 
  !      are undefined after this routine call.  This is avoided 
  !      by using intent(inout).
  !  PARAMETERS
  ! INOUT grid              coamps_grid to allocate space in
  !   IN  totalprocs        total number of non-I/O processors
  subroutine initialize_decomposition(grid, totalprocs)
    type(coamps_grid), intent(inout) :: grid
    integer,           intent(in)    :: totalprocs

    character(len=*), parameter :: routine = 'initialize_decomposition'
    integer :: alloc_status

    allocate( grid%iminf(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%iminf')
    allocate( grid%imaxf(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%imaxf')
    allocate( grid%jminf(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%jminf')
    allocate( grid%jmaxf(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%jmaxf')
  end subroutine initialize_decomposition

  ! decompose_domain
  ! ----------------
  ! Given the grid, the number of I/O processors, the number of
  ! processors in x and y, and the number of boundary points (nbnam
  ! in the COAMPS namelist), mimics the domain decomposition
  ! performed by COAMPS routine domdec.F (the math is taken from that
  ! module and variables renamed for use here).  The maximum and
  ! minimum [ij] points for each processor are given in arrays
  ! [ij](min|max)f within the coamps_grid structure.  This routine
  ! populates those arrays.  
  !  PARAMETERS
  ! INOUT grid              coamps_grid to decompose
  !   IN  num_io_procs      number of I/O processors (0 or 1)
  !   IN  dom_proc_x        number of processors in the x direction
  !   IN  dom_proc_y        number of processors in the y direction
  !   IN  nbound            number of halo boundary points
  subroutine decompose_domain(grid, num_io_procs, dom_proc_x, &
                              dom_proc_y, nbound) 
    type(coamps_grid), intent(inout) :: grid
    integer,           intent(in)    :: num_io_procs
    integer,           intent(in)    :: dom_proc_x
    integer,           intent(in)    :: dom_proc_y
    integer,           intent(in)    :: nbound

    ! Total number of processors we are iterating over, along with
    ! the current processor and some counters
    integer :: nproc_x, nproc_y
    integer :: totalprocs, curproc, x_count, y_count

    ! Number of points per processor (not including the boundary)
    integer :: xperproc, yperproc
    integer :: xextra, yextra

    ! The initial boundaries for the domain (not counting the
    ! halo points)
    integer, allocatable, dimension(:) :: begin_i
    integer, allocatable, dimension(:) :: begin_j
    integer, allocatable, dimension(:) :: end_i
    integer, allocatable, dimension(:) :: end_j

    ! Error handling
    character(len=*), parameter :: routine = 'decompose_domain'
    integer :: alloc_status, dealloc_status

    ! If we are only doing single processor I/O, it doesn't matter
    ! how many processors there are - we treat it as a single
    ! decomposition
    if (num_io_procs .gt. 0) then
       nproc_x = 1
       nproc_y = 1
    else
       nproc_x = dom_proc_x
       nproc_y = dom_proc_y
    endif
    totalprocs = nproc_x * nproc_y

    ! This sequence of allocations could potentially be offloaded to
    ! the initialize_decomposition subroutine but as the arrays in 
    ! question are local and temporary, I've kept them here.
    allocate( begin_i(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'begin_i')
    allocate( begin_j(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'begin_j')
    allocate( end_i(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'end_i')
    allocate( end_j(totalprocs), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'end_j')
    
    ! Get the processor-point counts
    xperproc = (grid%pts_x - 2) / nproc_x
    yperproc = (grid%pts_y - 2) / nproc_y
    xextra   = (grid%pts_x - 2) - (xperproc * nproc_x)
    yextra   = (grid%pts_y - 2) - (yperproc * nproc_y)

    do y_count = 1,nproc_y
       do x_count = 1,nproc_x

          ! 2D grid -> 1D number
          curproc = (y_count - 1) * nproc_x + x_count

          ! X points
          if (x_count .eq. 1) then
             begin_i(curproc) = 2
          else
             begin_i(curproc) = end_i(curproc-1) + 1
          endif

          end_i(curproc) = begin_i(curproc) + xperproc - 1
          if (x_count .gt. nproc_x-xextra) then
             end_i(curproc) = end_i(curproc) + 1
          endif
       
          ! Y Points
          if (y_count .eq. 1) then
             begin_j(curproc) = 2
          else
             begin_j(curproc) = end_j(curproc-nproc_x) + 1
          endif
          end_j(curproc) = begin_j(curproc) + yperproc - 1
          if (y_count .gt. nproc_y-yextra) then
             end_j(curproc) = end_j(curproc) + 1
          endif

          ! Take the boundaries into consideration
          grid%iminf(curproc) = begin_i(curproc) - nbound
          grid%imaxf(curproc) = end_i(curproc)   + nbound
          grid%jminf(curproc) = begin_j(curproc) - nbound
          grid%jmaxf(curproc) = end_j(curproc)   + nbound
       enddo
    enddo

    deallocate( begin_i, begin_j, end_i, end_j, stat=dealloc_status )
    call check_dealloc_status(dealloc_status, routine, source, &
                              revision, revdate, '[ij](min|max)f')
  end subroutine decompose_domain

  ! get_grid_iminf
  ! --------------
  ! Accessor method for the index of a subgrid's left extent
  !  PARAMETERS
  !   IN  grid              coamps_grid to get data from
  !   IN  subgrid_num       decomposed subgrid number
  !   OUT iminf             i-index in large domain for the
  !                         left edge of this subgrid 
  subroutine get_grid_iminf(grid, subgrid_num, iminf)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(in)  :: subgrid_num
    integer,           intent(out) :: iminf

    iminf = grid%iminf(subgrid_num)
  end subroutine get_grid_iminf

  ! get_grid_imaxf
  ! --------------
  ! Accessor method for the index of a subgrid's right extent
  !  PARAMETERS
  !   IN  grid              coamps_grid to get data from
  !   IN  subgrid_num       decomposed subgrid number
  !   OUT imaxf             i-index in large domain for the
  !                         right edge of this subgrid 
  subroutine get_grid_imaxf(grid, subgrid_num, imaxf)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(in)  :: subgrid_num
    integer,           intent(out) :: imaxf

    imaxf = grid%imaxf(subgrid_num)
  end subroutine get_grid_imaxf

  ! get_grid_jminf
  ! --------------
  ! Accessor method for the index of a subgrid's lower extent
  !  PARAMETERS
  !   IN  grid              coamps_grid to get data from
  !   IN  subgrid_num       decomposed subgrid number
  !   OUT jminf             j-index in large domain for the
  !                         lower edge of this subgrid
  subroutine get_grid_jminf(grid, subgrid_num, jminf)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(in)  :: subgrid_num
    integer,           intent(out) :: jminf

    jminf = grid%jminf(subgrid_num)
  end subroutine get_grid_jminf

  ! get_grid_jmaxf
  ! --------------
  ! Accessor method for the index of a subgrid's upper extent
  !  PARAMETERS
  !   IN  grid              coamps_grid to get data from
  !   IN  subgrid_num       decomposed subgrid number
  !   OUT jmaxf             j-index in large domain for the
  !                         upper edge of this subgrid
  subroutine get_grid_jmaxf(grid, subgrid_num, jmaxf)
    type(coamps_grid), intent(in)  :: grid
    integer,           intent(in)  :: subgrid_num
    integer,           intent(out) :: jmaxf

    jmaxf = grid%jmaxf(subgrid_num)
  end subroutine get_grid_jmaxf

  ! grid_wind_to_earth_wind
  ! ----------------
  ! Rotates the grid relative winds to a true wind
  !  PARAMETERS
  !   IN  grid              coamps_grid structure
  !   IN  dart_loc          dart_location
  !   INOUT v_wind          v wind component
  subroutine grid_wind_to_earth_wind(u_wind,v_wind,grid,dart_loc)
    real(kind=r8),       intent(inout) :: u_wind
    real(kind=r8),       intent(inout) :: v_wind
    type(coamps_grid),   intent(in)    :: grid
    type(location_type), intent(in)    :: dart_loc

    real(kind=r8), dimension(3) :: location
    real(kind=r8)               :: longitude
    real(kind=r8)               :: grid_rotation(1)

    real(kind=r8)               :: u_tmp(1,1)
    real(kind=r8)               :: v_tmp(1,1)

    u_tmp(1,1) = u_wind
    v_tmp(1,1) = v_wind
    location   = get_location(dart_loc)
    longitude  = location(1)
    grid_rotation(1)=calc_grid_rotation(grid,longitude)
    call uvg2uv(u_tmp,v_tmp,1,1,grid_rotation,u_tmp,v_tmp)
    u_wind     = u_tmp(1,1)
    v_wind     = v_tmp(1,1)
  end subroutine grid_wind_to_earth_wind

  ! dump_grid_info
  ! --------------
  ! Dumps a human-readable summary of the supplied grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to print out
  subroutine dump_grid_info(grid)
    type(coamps_grid), intent(in) :: grid

    if (do_output()) then
       write (*,*) "**** COAMPS GRID INFORMATION ****"
       write (*,*) "Number of pressure levels:", grid%pres_lvls
       write (*,*) "Number of sigma levels:   ", grid%sigm_lvls
       write (*,*) "Map projection type:      ", grid%map_proj
       write (*,*) "Latitude intersection 1:  ", grid%stdlat1
       write (*,*) "Latitude intersection 2:  ", grid%stdlat2
       write (*,*) "Longitude intersection:   ", grid%stdlon
       write (*,*) "Reference latitude:       ", grid%reflat
       write (*,*) "Reference longitude:      ", grid%reflon
       write (*,*) "X-direction spacing (m):  ", grid%delta_x
       write (*,*) "Y-direction spacing (m):  ", grid%delta_y
       write (*,*) "Number of nests:          ", grid%nests
       write (*,*) "Number of X grid points:  ", grid%pts_x
       write (*,*) "Number of Y grid points:  ", grid%pts_y
       write (*,*) "i point on parent grid:   ", grid%anchor_i
       write (*,*) "j point on parent grid:   ", grid%anchor_j
       write (*,*) "i point at ref. lat/lon:  ", grid%iref
       write (*,*) "j point at ref. lat/lon:  ", grid%jref
    end if
  end subroutine dump_grid_info

  ! nc_dump_grid_info
  ! --------------
  ! Dumps a netCDF summary of the supplied grid
  !  PARAMETERS
  !   IN  grid              coamps_grid to print out
  !   IN  ncFileID          netCDF file ID to print to
  subroutine nc_dump_grid_info(grid, ncFileID)

    use netcdf        ! for the routines
    use typeSizes     ! for NF90_GLOBAL attribute identifier

    type(coamps_grid), intent(in) :: grid
    integer,           intent(in) :: ncFileID

    character(len=*), parameter :: routine = 'nc_dump_grid_info'

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Number_of_pressure_levels", grid%pres_lvls), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Number_of_sigma_levels",    grid%sigm_lvls), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Map_projection_type",       grid%map_proj), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Latitude_intersection_1",   grid%stdlat1), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Latitude_intersection_2",   grid%stdlat2), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Longitude_intersection",    grid%stdlon), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Reference_latitude",        grid%reflat), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Reference_longitude",       grid%reflon), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "X_spacing_meters",   grid%delta_x), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Y_spacing_meters",   grid%delta_y), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Number_of_nests",           grid%nests), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Number_of_X_gridpoints",   grid%pts_x), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "Number_of_Y_gridpoints",   grid%pts_y), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "i_point_on_parent_grid",    grid%anchor_i), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "j_point_on_parent_grid",    grid%anchor_j), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "i_point_at_ref_lat_lon",    grid%iref), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
       "j_point_at_ref_lat_lon",    grid%jref), routine)
    
  end subroutine nc_dump_grid_info

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! generate_datahd_filename
  ! ------------------------
  ! Generates the COAMPS domain information file name for the first
  ! nest at a given date-time group.  This does *not* generate any
  ! path information - it only returns the file name.
  !  PARAMETERS
  !   IN  dtg               base date-time group for model run
  !   OUT datahd_filename   name of domain information file
  subroutine generate_datahd_filename(dtg, datahd_filename)
    character(len=10), intent(in)  :: dtg
    character(len=64), intent(out) :: datahd_filename

    ! The format of the datahd filename is fixed except for the date
    ! -time group: mimics a 2000x1 flat file with no level
    ! information or tau information
!   call generate_flat_file_name( var_name   = 'datahd',      &
!                                 level_type = 'sfc',         &
!                                 level1     = 0,             &
!                                 level2     = 0,             &
!                                 gridnum    = 1,             &
!                                 aoflag     = 'a',           &
!                                 xpts       = 2000,          &
!                                 ypts       = 1,             &
!                                 dtg        = dtg,           &
!                                 tau_dd     = 0,             &
!                                 tau_hh     = 0,             &
!                                 tau_mm     = 0,             &
!                                 tau_ss     = 0,             &
!                                 field_type = 'infofld',     &
!                                 file_name  = datahd_filename )

    call generate_one_flat_file_name( var_name   = 'datahd',      &
                                      level_type = 'sfc',         &
                                      level1     = 0,             &
                                      level2     = 0,             &
                                      gridnum    = 1,             &
                                      aoflag     = 'a',           &
                                      xpts       = 2000,          &
                                      ypts       = 1,             &
                                      dtg        = dtg,           &
                                      tau_hh     = 0,             &
                                      tau_mm     = 0,             &
                                      tau_ss     = 0,             &
                                      field_type = 'infofld',     &
                                     file_name  = datahd_filename )
  end subroutine generate_datahd_filename

  ! generate_terrht_filename
  ! ------------------------
  ! Generates the COAMPS terrain height flat file name for the first
  ! nest of a grid at a given date-time group.  This does *not* 
  ! generate any path information - it only returns the file name.
  !
  ! Assumes that the coamps_grid structure has already been populated
  ! with the domain information, since we need the number of x and y
  ! grid points to generate the file name
  !  PARAMETERS
  !   IN  dtg               base date-time group for model run 
  !   IN  grid              coamps_grid that we need terrain info for
  !   OUT terrht_filename   name of terrain height flat file
  subroutine generate_terrht_filename(dtg, grid, terrht_filename)
    character(len=10), intent(in)  :: dtg
    type(coamps_grid), intent(in)  :: grid
    character(len=64), intent(out) :: terrht_filename
    
    integer :: x_points
    integer :: y_points

    x_points = int(grid%pts_x)
    y_points = int(grid%pts_y)

    ! The format of the terrain height file is fixed except for
    ! the horizontal grid size and date-time group
!   call generate_flat_file_name( var_name   = 'terrht',      &
!                                 level_type = 'sfc',         &
!                                 level1     = 0,             &
!                                 level2     = 0,             &
!                                 gridnum    = 1,             &
!                                 aoflag     = 'a',           &
!                                 xpts       = x_points,      &
!                                 ypts       = y_points,      &
!                                 dtg        = dtg,           &
!                                 tau_dd     = 0,             &
!                                 tau_hh     = 0,             &
!                                 tau_mm     = 0,             &
!                                 tau_ss     = 0,             &
!                                 field_type = 'fcstfld',     &
!                                 file_name  = terrht_filename )

    call generate_one_flat_file_name( var_name   = 'terrht',      &
                                      level_type = 'sfc',         &
                                      level1     = 0,             &
                                      level2     = 0,             &
                                      gridnum    = 1,             &
                                      aoflag     = 'a',           &
                                      xpts       = x_points,      &
                                      ypts       = y_points,      &
                                      dtg        = dtg,           &
                                      tau_hh     = 0,             &
                                      tau_mm     = 0,             &
                                      tau_ss     = 0,             &
                                      field_type = 'fcstfld',     &
                                      file_name  = terrht_filename )
  end subroutine generate_terrht_filename

!  ! generate_flat_file_name
!  ! -----------------------
!  ! Given field, level, and grid information, generate the properly
!  ! formatted 64-character COAMPS flat file name.  Note that this
!  ! does *not* generate any path information - it only returns the
!  ! file name.
!  !  PARAMETERS
!  !   IN  var_name          the field the file contains
!  !   IN  level_type        vertical level type (height/pressure/etc)
!  !   IN  level1            lowest vertical level in the file
!  !   IN  level2            highest vertical level in the file
!  !                         (for files for a single level, level1 is 
!  !                          that level and level2 is left to 0)
!  !   IN  gridnum           nest number (only 1 supported for now)
!  !   IN  aoflag            field type: (a)tmosphere or (o)cean
!  !   IN  xpts              number of points in the x direction
!  !   IN  ypts              number of points in the y direction
!  !   IN  dtg               base date-time group
!  !   IN  tau_dd            forecast lead time - day component
!  !   IN  tau_hh            forecast lead time - hour component
!  !   IN  tau_mm            forecast lead time - minute component
!  !   IN  tau_ss            forecast lead time - second component
!  !   IN  field_type        type of field (e.g. fcstfld, infofld)
!  !   OUT file_name         COAMPS flat file name
!  subroutine generate_flat_file_name(var_name, level_type, level1, &
!                                     level2, gridnum, aoflag, xpts,&
!                                     ypts, dtg, tau_dd, tau_hh,    &
!                                     tau_mm, tau_ss, field_type,   &
!                                     file_name)
!    character(len=6), intent(in)   :: var_name
!    character(len=3), intent(in)   :: level_type
!    integer, intent(in)            :: level1
!    integer, intent(in)            :: level2
!    integer, intent(in)            :: gridnum
!    character(len=1), intent(in)   :: aoflag
!    integer, intent(in)            :: xpts
!    integer, intent(in)            :: ypts
!    character(len=10), intent(in)  :: dtg
!    integer, intent(in)            :: tau_dd
!    integer, intent(in)            :: tau_hh
!    integer, intent(in)            :: tau_mm
!    integer, intent(in)            :: tau_ss
!    character(len=7), intent(in)   :: field_type
!    character(len=64), intent(out) :: file_name
!
!    write(file_name, 100) var_name, level_type, level1, level2,&
!         & gridnum, aoflag, xpts, ypts, dtg, tau_dd, tau_hh, tau_mm,&
!         & tau_ss, field_type
!
!100 format(A6'_'A3'_'I6.6'_'I6.6'_'I1A1I4.4'x'I4.4'_'A10'_'I2.2I2.2&
!         I2.2I2.2'_'A7)
!  end subroutine generate_flat_file_name

  ! read_datahd_file
  ! ----------------
  ! Given the unit number of an *open* COAMPS domain information file,
  ! populate the portion of a COAMPS grid structure that corresponds
  ! to the data contained in the file 
  !  PARAMETERS
  !   IN  datahd_unit       unit number of open datahd file
  !   OUT grid              coamps_grid structure to be filled
  subroutine read_datahd_file(datahd_unit, grid)
    integer,           intent(in)  :: datahd_unit
    type(coamps_grid), intent(out) :: grid

    ! Domain information from file
    real(kind=r8), dimension(2000) :: domain_data
    
    ! Model sigma levels
    integer :: kka
    integer :: cur_lev 

    ! True if a the datahd file is opened & if memory's allocated
    logical :: is_opened, is_allocated

    ! Error checking
    character(len=*), parameter :: routine = 'read_datahd_file'
    integer :: io_status, alloc_status

    ! This is really an "assert" statement rather than error capture
    ! since the user should not have to worry about this
    inquire(unit=datahd_unit, opened=is_opened)
    if ( .not. is_opened ) then
      call error_handler(E_ERR,'read_datahd_file','Trying to read a &
                         &datahd file that is not open yet', source,&
                         &revision, revdate)
    end if
    
    ! Probably shouldn't do this as free-form I/O, but it works so
    ! stick with it for now - maybe add a real format statement later.
    read(unit=datahd_unit, fmt=*,iostat=io_status) domain_data
    call check_io_status(io_status, routine, source, revision, &
                         revdate, 'Reading datahd file')

    ! Populate the grid structure - do the scalar variables first
    ! This information can be found around line 875 in the COAMPS
    ! atmos/libsrc/aalib/coama.F file
    grid%pres_lvls = domain_data(DATAHD_NUM_P_LEVELS)
    grid%sigm_lvls = domain_data(DATAHD_NUM_S_LEVELS)
    grid%map_proj  = domain_data(DATAHD_MAP_PROJECTION)
    grid%stdlat1   = domain_data(DATAHD_STANDARD_LAT_1)
    grid%stdlat2   = domain_data(DATAHD_STANDARD_LAT_2)
    grid%stdlon    = domain_data(DATAHD_STANDARD_LON)
    grid%reflat    = domain_data(DATAHD_REFERENCE_LAT)
    grid%reflon    = domain_data(DATAHD_REFERENCE_LON)
    grid%delta_x   = domain_data(DATAHD_DELTA_X)
    grid%delta_y   = domain_data(DATAHD_DELTA_Y) 
    grid%nests     = domain_data(DATAHD_NUM_NESTS)

    ! Since COAMPS doesn't output the grid constant for the projection
    ! we need to caluclate it here.  Needed to rotate between grid wind
    ! and true wind.
    grid%grid_constant=calc_grid_constant(grid)

    ! Individual nest information - currently only 1 next so
    ! explicitly specify indices.  This would be a loop later if I
    ! allowed nested grids.
    grid%pts_x    = domain_data(DATAHD_NUM_X_POINTS)
    grid%pts_y    = domain_data(DATAHD_NUM_Y_POINTS)
    grid%anchor_i = domain_data(DATAHD_ANCHOR_POINT_I)
    grid%anchor_j = domain_data(DATAHD_ANCHOR_POINT_J)
    grid%iref     = domain_data(DATAHD_REFERENCE_I)
    grid%jref     = domain_data(DATAHD_REFERENCE_J)

    ! Again, this is an assert more than an error trap
    is_allocated = allocated(grid%dsigma)  .and. &
                   allocated(grid%dsigmaw) .and. &
                   allocated(grid%msigma)  .and. &
                   allocated(grid%wsigma)
    if (is_allocated) then
      call error_handler(E_ERR, routine, 'Grid sigma arrays are ' // &
                         'already allocated', source, revision, revdate)
    end if
              
    kka = grid%sigm_lvls
    allocate( grid%dsigma(kka),      stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%dsigma')
    allocate( grid%dsigmaw(kka + 1), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%dsigmaw')
    allocate( grid%msigma(kka),      stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%msigma')
    allocate( grid%wsigma(kka + 1),  stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'grid%wsigma')

    ! Now that we have the space, information about mass sigma 
    ! levels can just be read in from the file...
    do cur_lev=1,kka
       grid%dsigma(cur_lev) = domain_data(DATAHD_DSIGMA_OFFSET + cur_lev)
       grid%msigma(cur_lev) = domain_data(DATAHD_MSIGMA_OFFSET + cur_lev)
    end do

    ! ... but we have to calculate the w sigma level information
    ! (this is based on the COAMPS utility package source code)
    grid%dsigmaw(TOP_LEVEL) = grid%dsigma(TOP_LEVEL)   * 0.5_r8
    grid%dsigmaw(kka + 1)   = grid%dsigma(kka)         * 0.5_r8
    do cur_lev = 1, (kka - 1)
       grid%dsigmaw(cur_lev+1) = grid%msigma(cur_lev) - grid%msigma(cur_lev+1)
    end do
    grid%wsigma(kka + 1) = 0.0_r8
    do cur_lev=1,kka
       grid%wsigma(kka+1 - cur_lev) = grid%wsigma(kka+1 - cur_lev + 1) + &
                                      grid%dsigma(kka+1 - cur_lev)
    end do
  end subroutine read_datahd_file

  ! calc_grid_rotation
  ! ----------------
  ! Given a coamps grid, calculates the grid constant for the projection
  !  PARAMETERS
  !   IN  grid              coamps_grid structure
  !   OUT grid_constant     projection grid constant
  function calc_grid_constant(grid) result(grid_constant)
    type(coamps_grid), intent(in)  :: grid
    real(kind=r8)                  :: grid_constant

    select case(int(grid%map_proj))
      case(PROJ_MERCATOR)
        grid_constant=0.0_r8
      case(PROJ_LAMBERT_CONIC_CONFORMAL)
        if (grid%stdlat1 == grid%stdlat2) then
          grid_constant=sin(abs(grid%stdlat1)*DEG_TO_RAD)
        else
          grid_constant=(log(sin((90.0_r8-abs(grid%stdlat1))*DEG_TO_RAD))      &
                        -log(sin((90.0_r8-abs(grid%stdlat2))*DEG_TO_RAD)))     &
                       /(log(tan((90.0_r8-abs(grid%stdlat1))*DEG_TO_RAD*0.5_r8))  &
                        -log(tan((90.0_r8-abs(grid%stdlat2))*DEG_TO_RAD*0.5_r8)))
        endif
      case(PROJ_POLAR_STEREOGRAPHIC)
        grid_constant=1.0_r8
      case default
        grid_constant=1.0_r8
    end select
    return
  end function calc_grid_constant


  ! calc_grid_rotation
  ! ----------------
  ! Given a longitude and grid structure, calculates the 
  ! rotation angle between grid north and true north.
  !  PARAMETERS
  !   IN  grid              coamps_grid structure
  !   IN  longitude         longitude of point
  !   OUT grid_rotation     rotation angle.
  function calc_grid_rotation(grid,longitude) result(grid_rotation)
    type(coamps_grid),   intent(in)  :: grid
    real(kind=r8),       intent(in)  :: longitude

    real(kind=r8)                    :: grid_rotation

    select case(int(grid%map_proj))
      case(PROJ_LAMBERT_CONIC_CONFORMAL, PROJ_POLAR_STEREOGRAPHIC)

        grid_rotation=(grid%stdlon-longitude)*grid%grid_constant
        !case 1: standard longitude is east of Meridian, (between 0 and 90E)
        !and the grid longitude is west of Meridian (between 90W and Meridian)
        if (grid%stdlon >=   0.0_r8 .and. grid%stdlon <= 90.0_r8 .and.  &
            longitude   >= 270.0_r8 .and.   longitude < 360.0_r8) then
                grid_rotation=(grid%stdlon+360.0_r8-longitude)*grid%grid_constant
        endif
        !case 2: standard longitude is west of Meridian, (between 90W and Meridian)
        !and the grid longitude is east of Meridian (between 0 and 90E)
        if (grid%stdlon >= 270.0_r8 .and. grid%stdlon <  360.0_r8 .and. &
            longitude   >=   0.0_r8 .and. longitude   <=  90.0_r8) then
              grid_rotation=(grid%stdlon-longitude-360.0)*grid%grid_constant
        endif
      case default
        grid_rotation=0.0_r8
    end select
    return
  end function calc_grid_rotation

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_grid_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
