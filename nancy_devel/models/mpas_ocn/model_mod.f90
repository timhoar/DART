! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This is the interface between the model model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI, MISSING_I
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      & 
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_obs_init, get_close_obs_destroy,        &
                             loc_get_close_obs => get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper, nmlfileunit,       &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file, do_nml_file, do_nml_term

use     obs_kind_mod, only : paramname_length,        &
                             get_raw_obs_kind_index,  &
                             get_raw_obs_kind_name,   &
                             KIND_SURFACE_PRESSURE,   &
                             KIND_VERTICAL_VELOCITY,  &
                             KIND_POTENTIAL_TEMPERATURE, &
                             KIND_EDGE_NORMAL_SPEED,  &
                             KIND_TEMPERATURE,        &
                             KIND_U_WIND_COMPONENT,   &
                             KIND_V_WIND_COMPONENT,   &
                             KIND_PRESSURE,           &
                             KIND_DENSITY,            & 
                             KIND_VAPOR_MIXING_RATIO

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf
use get_geometry_mod
use get_reconstruct_mod


implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: get_model_analysis_filename,  &
          get_grid_definition_filename, &
          analysis_file_to_statevector, &
          statevector_to_analysis_file, &
          get_analysis_time,            &
          write_model_time,             &
          get_grid_dims

! version controlled file description for error handling, do not edit

character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

character(len=256) :: string1, string2
logical, save :: module_initialized = .false.

! Real (physical) constants as defined exactly in MPAS.
! redefined here for consistency with the model.
real(r8), parameter :: rgas = 287.0_r8
real(r8), parameter :: cp = 1003.0_r8
real(r8), parameter :: cv = 716.0_r8
real(r8), parameter :: p0 = 100000.0_r8
real(r8), parameter :: rcv = rgas/(cp-rgas)

! FIXME: one of the example ocean files had a global attr with 6371220.0
! instead of 1229.   ??
real(r8), parameter :: radius = 6371229.0 ! meters

! Storage for a random sequence for perturbing a single initial state

type(random_seq_type) :: random_seq

! Structure for computing distances to cell centers, and assorted arrays
! needed for the get_close code.
type(get_close_type)             :: cc_gc
type(location_type), allocatable :: cell_locs(:)
integer, allocatable             :: dummy(:), close_ind(:)

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.0001   ! tiny amounts
logical            :: output_state_vector = .true.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_analysis_filename = 'mpas_analysis.nc'
character(len=256) :: grid_definition_filename = 'mpas_analysis.nc'
logical            :: use_new_code = .true.

namelist /model_nml/             &
   model_analysis_filename,      &
   grid_definition_filename,     &
   output_state_vector,          &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_perturbation_amplitude, &
   calendar,                     &
   debug,                        &
   use_new_code

!------------------------------------------------------------------
! DART state vector are specified in the input.nml:mpas_vars_nml namelist.
!------------------------------------------------------------------

integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
character(len=NF90_MAX_NAME) :: mpas_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

namelist /mpas_vars_nml/ mpas_state_variables

integer :: nfields

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME), dimension(NF90_MAX_VAR_DIMS) :: dimname
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: xtype         ! netCDF variable type (NF90_double, etc.) 
   integer :: numdims       ! number of dims - excluding TIME
   integer :: numvertical   ! number of vertical levels in variable
   integer :: numcells      ! number of horizontal locations (cell centers)
   integer :: numedges      ! number of horizontal locations (edges for velocity components)
   logical :: ZonHalf       ! vertical coordinate has dimension nVertLevels
   integer :: varsize       ! prod(dimlens(1:numdims))
   integer :: index1        ! location in dart state vector of first occurrence
   integer :: indexN        ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=paramname_length) :: kind_string
   logical  :: clamping     ! does variable need to be range-restricted before 
   real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! Grid parameters - the values will be read from an mpas analysis file. 

integer :: nCells        = -1  ! Total number of cells making up the grid
integer :: nVertices     = -1  ! Unique points in grid that are corners of cells
integer :: nEdges        = -1  ! Straight lines between vertices making up cells
integer :: maxEdges      = -1  ! Largest number of edges a cell can have
integer :: nVertLevels   = -1  ! Vertical levels; count of vert cell centers
integer :: nVertLevelsP1 = -1  ! Vert levels plus 1; count of vert cell faces
integer :: vertexDegree  = -1  ! Max number of cells/edges that touch any vertex
integer :: nSoilLevels   = -1  ! Number of soil layers

! scalar grid positions

! FIXME: we read in a lot of metadata about the grids.  if space becomes an
! issue we could consider reading in only the x,y,z arrays for all the items
! plus the radius, and then compute the lat/lon for locations needed by 
! get_state_meta_data() on demand.  most of the computations we need to do
! are actually easier in xyz coords (no issue with poles).

! FIXME: it may be desirable to read in xCell(:), yCell(:), zCell(:)
! to keep from having to compute them on demand, especially since we
! have converted the radian lat/lon of the cell centers into degrees.
! we have to convert back, then take a few sin and cos to get xyz.
! time/space/accuracy tradeoff here.

real(r8), allocatable :: xVertex(:), yVertex(:), zVertex(:)
real(r8), allocatable :: xEdge(:), yEdge(:), zEdge(:)
real(r8), allocatable :: lonEdge(:) ! edge longitudes (degrees, original radians in file)
real(r8), allocatable :: latEdge(:) ! edge longitudes (degrees, original radians in file)
real(r8), allocatable :: lonCell(:) ! cell center longitudes (degrees, original radians in file)
real(r8), allocatable :: latCell(:) ! cell center latitudes  (degrees, original radians in file)
real(r8), allocatable :: zGridFace(:,:)   ! geometric height at cell faces   (nVertLevelsP1,nCells)
real(r8), allocatable :: zGridCenter(:,:) ! geometric height at cell centers (nVertLevels,  nCells)
real(r8), allocatable :: zGridEdge(:,:)   ! geometric height at edge centers (nVertLevels,  nEdges)

!real(r8), allocatable :: zEdgeFace(:,:)   ! geometric height at edges faces  (nVertLevelsP1,nEdges)
!real(r8), allocatable :: zEdgeCenter(:,:) ! geometric height at edges faces  (nVertLevels  ,nEdges)

integer,  allocatable :: cellsOnVertex(:,:) ! list of cell centers defining a triangle
integer,  allocatable :: verticesOnCell(:,:)

integer,  allocatable :: edgesOnCell(:,:) ! list of edges that bound each cell
integer,  allocatable :: cellsOnEdge(:,:) ! list of cells that bound each edge
integer,  allocatable :: nedgesOnCell(:) ! list of edges that bound each cell
real(r8), allocatable :: edgeNormalVectors(:,:)

! Boundary information might be needed ... regional configuration?
! Read if available.

integer,  allocatable :: boundaryEdge(:,:)
integer,  allocatable :: boundaryVertex(:,:)
integer,  allocatable :: maxLevelCell(:)

real(r8), allocatable :: ens_mean(:)         ! may be needed for forward ops

integer         :: model_size          ! the state vector length
type(time_type) :: model_timestep      ! smallest time to adv model

logical :: global_grid = .true.        ! true = the grid covers the sphere with no holes
logical :: all_levels_exist_everywhere = .true. ! true = cells defined at all levels
logical :: has_real_u = .false.        ! true = has original u on edges
logical :: has_uvreconstruct = .false. ! true = has reconstructed at centers

!------------------------------------------------------------------
! The model analysis manager namelist variables
!------------------------------------------------------------------

character(len= 64) :: ew_boundary_type, ns_boundary_type

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

INTERFACE prog_var_to_vector
      MODULE PROCEDURE prog_var_1d_to_vector
      MODULE PROCEDURE prog_var_2d_to_vector
      MODULE PROCEDURE prog_var_3d_to_vector
END INTERFACE

INTERFACE get_analysis_time
      MODULE PROCEDURE get_analysis_time_ncid
      MODULE PROCEDURE get_analysis_time_fname
END INTERFACE

INTERFACE get_index_range
      MODULE PROCEDURE get_index_range_int
      MODULE PROCEDURE get_index_range_string
END INTERFACE

!------------------------------------------------

! The regular grid used for triangle interpolation divides the sphere into
! a set of regularly spaced lon-lat boxes. The number of boxes in
! longitude and latitude are set by num_reg_x and num_reg_y. Making the
! number of regular boxes smaller decreases the computation required for
! doing each interpolation but increases the static storage requirements
! and the initialization computation (which seems to be pretty small).
integer, parameter :: num_reg_x = 90, num_reg_y = 90

! The max_reg_list_num controls the size of temporary storage used for
! initializing the regular grid. Two arrays
! of size num_reg_x*num_reg_y*max_reg_list_num are needed. The initialization
! fails and returns an error if max_reg_list_num is too small. A value of
! ??? is sufficient for ???
integer, parameter :: max_reg_list_num = 100

! The triangle interpolation keeps a list of how many and which triangles
! overlap each regular lon-lat box. The number is stored in
! array triangle_num. The allocatable array
! triangle_list lists the uniquen index 
! of each overlapping triangle. The entry in
! triangle_start for a given regular lon-lat box indicates
! where the list of triangles begins in the triangle_list.

integer :: triangle_start(num_reg_x, num_reg_y)
integer :: triangle_num  (num_reg_x, num_reg_y) = 0
integer, allocatable :: triangle_list(:)

contains

!==================================================================
! All the public REQUIRED interfaces come first - just by convention.
!==================================================================


!------------------------------------------------------------------

subroutine static_init_model()

! Called to do one time initialization of the model.
! 
! All the grid information comes from the initialization of
! the dart_model_mod module.

! Local variables - all the important ones have module scope


integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname,dimname
character(len=paramname_length)       :: kind_string
integer :: ncid, VarID, numdims, varsize, dimlen
integer :: iunit, io, ivar, i, index1, indexN, iloc, kloc
integer :: ss, dd
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
integer :: cel1, cel2

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Read the MPAS variable list to populate DART state vector
! Intentionally do not try to dump them to the nml unit because
! they include large character arrays which output pages of white space.
! The routine that reads and parses this namelist will output what
! values it found into the log.
call find_namelist_in_file('input.nml', 'mpas_vars_nml', iunit)
read(iunit, nml = mpas_vars_nml, iostat = io)
call check_namelist_read(iunit, io, 'mpas_vars_nml')

!---------------------------------------------------------------
! Set the time step ... causes mpas namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids 
! 3) read them from the analysis file

! read_grid_dims() fills in the following module global variables:
!  nCells, nVertices, nEdges, maxEdges, nVertLevels, nVertLevelsP1, vertexDegree, nSoilLevels
call read_grid_dims()

allocate(latCell(nCells), lonCell(nCells)) 
allocate(zGridFace(nVertLevelsP1, nCells))
allocate(zGridCenter(nVertLevels, nCells))

!allocate(zEdgeFace(  nVertLevelsP1, nEdges))
!allocate(zEdgeCenter(nVertLevels,   nEdges))
allocate(zGridEdge(nVertLevels,   nEdges))

allocate(cellsOnVertex(vertexDegree, nVertices))
allocate(nEdgesOnCell(nCells))
allocate(edgesOnCell(maxEdges, nCells))
allocate(cellsOnEdge(2, nEdges))
allocate(verticesOnCell(maxEdges, nCells))
allocate(edgeNormalVectors(3, nEdges))
allocate(latEdge(nEdges), lonEdge(nEdges)) 
allocate(xVertex(nVertices), yVertex(nVertices), zVertex(nVertices))
allocate(xEdge(nEdges), yEdge(nEdges), zEdge(nEdges))

! this reads in latCell, lonCell, zGridFace, cellsOnVertex
call get_grid()

! read in vert cell face locations and then compute vertical center locations
do kloc=1, nCells
 do iloc=1, nVertLevels
   zGridCenter(iloc,kloc) = (zGridFace(iloc,kloc) + zGridFace(iloc+1,kloc))*0.5_r8
 enddo
enddo

! FIXME: This code is supposed to check whether an edge has 2 neighbours or 1 neighbour and then
!        compute the height accordingly.  HOWEVER, the array cellsOnEdge does not change with 
!        depth, but it should as an edge may have 2 neighbour cells at the top but not at depth.
do kloc=1, nEdges
 do iloc=1, nVertLevels
   cel1 = cellsOnEdge(1,kloc)
   cel2 = cellsOnEdge(2,kloc)
   if (cel1>0 .and. cel2>0) then
      zGridEdge(iloc,kloc) = (zGridCenter(iloc,cel1) + zGridCenter(iloc,cel2))*0.5_r8
   else if (cel1>0) then
      zGridEdge(iloc,kloc) = zGridCenter(iloc,cel1)
   else if (cel2>0) then
      zGridEdge(iloc,kloc) = zGridCenter(iloc,cel2)
   else  !this is bad...
      write(string1,*)'Edge ',kloc,' at vertlevel ',iloc,' has no neighbouring cells!'
      call error_handler(E_ERR,'static_init_model', string1, source, revision, revdate)
   endif
 enddo
enddo
              
!---------------------------------------------------------------
! Compile the list of model variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the model analysis file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the model
! analysis file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.


call nc_check( nf90_open(trim(model_analysis_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(model_analysis_filename))

call verify_state_variables( mpas_state_variables, ncid, model_analysis_filename, &
                             nfields, variable_table )

TimeDimID = FindTimeDimension( ncid )

if (TimeDimID < 0 ) then
   write(string1,*)'unable to find a dimension named Time.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                    'static_init_model', 'inquire '//trim(model_analysis_filename))

if ( (TimeDimID > 0) .and. (unlimitedDimID > 0) .and. (TimeDimID /= unlimitedDimID)) then
   write(string1,*)'IF Time is not the unlimited dimension, I am lost.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%numdims     = 0
   progvar(ivar)%numvertical = 1
   progvar(ivar)%dimlens     = MISSING_I
   progvar(ivar)%numcells    = MISSING_I
   progvar(ivar)%numedges    = MISSING_I

   string2 = trim(model_analysis_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, xtype=progvar(ivar)%xtype, &
           dimids=dimIDs, ndims=numdims), 'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them. 
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Since we are not concerned with the TIME dimension, we need to skip it.
   ! When the variables are read, only a single timestep is ingested into
   ! the DART state vector.

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,numdims

      if (dimIDs(i) == TimeDimID) cycle DimensionLoop

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                                          'static_init_model', string1)

      progvar(ivar)%numdims    = progvar(ivar)%numdims + 1
      progvar(ivar)%dimlens(i) = dimlen
      progvar(ivar)%dimname(i) = trim(dimname)
      varsize = varsize * dimlen

      select case ( dimname(1:6) )
         case ('nCells')
            progvar(ivar)%numcells = dimlen
         case ('nEdges')
            progvar(ivar)%numedges = dimlen
         case ('nVertL')  ! nVertLevels, nVertLevelsP1, nVertLevelsP2
            progvar(ivar)%numvertical = dimlen
         case ('nSoilL')  ! nSoilLevels
            progvar(ivar)%numvertical = dimlen
      end select

   enddo DimensionLoop

   call set_variable_clamping(ivar)

   if (progvar(ivar)%numvertical == nVertLevels) then
      progvar(ivar)%ZonHalf = .TRUE.
   else
      progvar(ivar)%ZonHalf = .FALSE.
   endif

   if (varname == 'u') has_real_u = .true.
   if (varname == 'uReconstructZonal' .or. &
       varname == 'uReconstructMeridional') has_uvreconstruct = .true.

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   !if ( debug > 0 ) call dump_progvar(ivar)

enddo

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(model_analysis_filename))

model_size = progvar(nfields)%indexN

if ( debug > 0 .and. do_output()) then
  write(logfileunit,*)
  write(     *     ,*)
  write(logfileunit,'(" static_init_model: nCells, nVertices, nVertLevels =",3(1x,i6))') &
                                          nCells, nVertices, nVertLevels
  write(     *     ,'(" static_init_model: nCells, nVertices, nVertLevels =",3(1x,i6))') &
                                          nCells, nVertices, nVertLevels
  write(logfileunit, *)'static_init_model: model_size = ', model_size
  write(     *     , *)'static_init_model: model_size = ', model_size
  if ( global_grid ) then
     write(logfileunit, *)'static_init_model: grid is a global grid '
     write(     *     , *)'static_init_model: grid is a global grid '
  else
     write(logfileunit, *)'static_init_model: grid has boundaries '
     write(     *     , *)'static_init_model: grid has boundaries ' 
  endif
  if ( all_levels_exist_everywhere ) then
     write(logfileunit, *)'static_init_model: all cells have same number of vertical levels '
     write(     *     , *)'static_init_model: all cells have same number of vertical levels '
  else
     write(logfileunit, *)'static_init_model: cells have varying number of vertical levels ' 
     write(     *     , *)'static_init_model: cells have varying number of vertical levels '
  endif
endif

allocate( ens_mean(model_size) )

! Initialize the interpolation data structures
call init_interp()

end subroutine static_init_model


!------------------------------------------------------------------

subroutine get_state_meta_data(index_in, location, var_type)

! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with KIND_

! passed variables

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type
  
! Local variables

integer  :: nxp, nzp, iloc, vloc, nf, n
integer  :: myindx
real(r8) :: height

if ( .not. module_initialized ) call static_init_model

myindx = -1
nf     = -1

! Determine the right variable
FindIndex : do n = 1,nfields
    if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) THEN
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      exit FindIndex
    endif
enddo FindIndex

if( myindx == -1 ) then
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

! Now that we know the variable, find the cell or edge

if (     progvar(nf)%numcells /= MISSING_I) then
   nxp = progvar(nf)%numcells
elseif ( progvar(nf)%numedges /= MISSING_I) then
   nxp = progvar(nf)%numedges
else
     write(string1,*) 'ERROR, ',trim(progvar(nf)%varname),' is not defined on edges or cells'
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

nzp  = progvar(nf)%numvertical
iloc = 1 + (myindx-1) / nzp    ! cell index
vloc = myindx - (iloc-1)*nzp   ! vertical level index

! the zGrid array contains the location of the cell top and bottom faces, so it has one 
! more value than the number of cells in each column.  for locations of cell centers
! you have to take the midpoint of the top and bottom face of the cell.
if (progvar(nf)%numedges /= MISSING_I) then
   if ( progvar(nf)%ZonHalf ) then
      height = zGridEdge(vloc,iloc)
   else
      call error_handler(E_ERR, 'get_state_meta_data', 'no support for edges at face heights', &
                         source, revision, revdate)
   endif
else
   if ( progvar(nf)%ZonHalf ) then
      height = zGridCenter(vloc,iloc)
   else if (nzp <= 1) then
      height = zGridFace(1,iloc)
   else
      height = zGridFace(vloc,iloc)
   endif
endif



if (progvar(nf)%numedges /= MISSING_I) then
   if (nzp <= 1) then
      location = set_location(lonEdge(iloc),latEdge(iloc), height, VERTISSURFACE)
   else
      location = set_location(lonEdge(iloc),latEdge(iloc), height, VERTISHEIGHT)
   endif
else ! must be on cell centers
   if (nzp <= 1) then
      location = set_location(lonCell(iloc),latCell(iloc), height, VERTISSURFACE)
   else
      location = set_location(lonCell(iloc),latCell(iloc), height, VERTISHEIGHT)
   endif
endif

if (debug > 9) then

    write(*,'("INDEX_IN / myindx / IVAR / NX, NZ: ",2(i10,2x),3(i5,2x))') index_in, myindx, nf, nxp, nzp
    write(*,'("                       ILOC, KLOC: ",2(i5,2x))') iloc, vloc
    write(*,'("                      LON/LAT/HGT: ",3(f12.3,2x))') lonCell(iloc), latCell(iloc), height

endif

if (present(var_type)) then
   var_type = progvar(nf)%dart_kind
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------

subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

! given a state vector, a location, and a KIND_xxx, return the
! interpolated value at that location, and an error code.  0 is success,
! anything positive is an error.  (negative reserved for system use)
!
! for specific error codes, see the local_interpolate() comments

! passed variables

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: interp_val
integer,             intent(out) :: istatus

! local storage

integer  :: obs_kinds(3), ivar
real(r8) :: values(3)

! call the normal interpolate code.  if it fails because
! the kind doesn't exist directly in the state vector, try a few
! kinds we know how to convert.

interp_val = MISSING_R8
istatus    = 888888       ! must be positive (and integer)

obs_kinds(1) = obs_type

! debug only
if (.false.) then
   ivar = get_progvar_index_from_kind(obs_kinds(1))
   if ( ivar > 0 .and. debug > 7 ) call dump_progvar(ivar, x)
   ivar = get_progvar_index_from_kind(KIND_POTENTIAL_TEMPERATURE)
   if ( ivar > 0 .and. debug > 7 ) call dump_progvar(ivar, x)
   ivar = get_progvar_index_from_kind(KIND_DENSITY)
   if ( ivar > 0 .and. debug > 7 ) call dump_progvar(ivar, x)
   ivar = get_progvar_index_from_kind(KIND_VAPOR_MIXING_RATIO)
   if ( ivar > 0 .and. debug > 7 ) call dump_progvar(ivar, x)
   ivar = get_progvar_index_from_kind(KIND_EDGE_NORMAL_SPEED)
   if ( ivar > 0 .and. debug > 7 ) call dump_progvar(ivar, x)
endif


call local_interpolate(x, location, 1, obs_kinds, values, istatus)
if (istatus /= 88) then
   ! this is for debugging - when we're confident the code is
   ! returning consistent values and rc codes, both these tests can
   ! be removed for speed.  FIXME.
   if (istatus /= 0 .and. values(1) /= MISSING_R8) then
      write(string1,*) 'interp routine returned a bad status but not a MISSING_R8 value'
      write(string2,*) 'value = ', values(1), ' istatus = ', istatus
      call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate, &
                         text2=string2)
   endif
   if (istatus == 0 .and. values(1) == MISSING_R8) then
      write(string1,*) 'interp routine returned a good status but set value to MISSING_R8'
      call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
   endif
   interp_val = values(1)
if (debug > 5) print *, 'returning value ', interp_val
   return
endif

if (obs_type == KIND_TEMPERATURE) then
   obs_kinds(1) = KIND_POTENTIAL_TEMPERATURE
   obs_kinds(2) = KIND_DENSITY
   obs_kinds(3) = KIND_VAPOR_MIXING_RATIO

   call local_interpolate(x, location, 3, obs_kinds, values, istatus)
!print *, '1 local interpolate returns istatus = ', istatus
   if (istatus /= 0) then
      ! this is for debugging - when we're confident the code is
      ! returning consistent values and rc codes, both these tests can
      ! be removed for speed.  FIXME.
      if (istatus /= 0 .and. .not. any(values == MISSING_R8)) then
         write(string1,*) 'interp routine returned a bad status but good values'
         call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
      endif
      if (istatus == 0 .and. any(values == MISSING_R8)) then
         write(string1,*) 'interp routine returned a good status but bad values'
         call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
      endif
      interp_val = MISSING_R8
      return
   endif

   ! compute sensible temperature as return value
   interp_val = theta_to_tk(values(1), values(2), values(3))
   return
endif

!print *, '2 local interpolate returns istatus = ', istatus

! add other cases here for kinds we want to handle
! if (istatus == 88) then
!  could try different interpolating different kinds here.
!  if (obs_type == KIND_xxx) then
!    code goes here
!  endif
! endif

! shouldn't get here.
interp_val = MISSING_R8
istatus = 100

end subroutine model_interpolate


!------------------------------------------------------------------

subroutine local_interpolate(x, location, num_kinds, obs_kinds, interp_vals, istatus)

! THIS VERSION IS ONLY CALLED INTERNALLY, so we can mess with the arguments.
! For a given lat, lon, and height, interpolate the correct state values
! to that location for the filter from the model state vectors.  this version
! can take multiple kinds at a single location.  there is only a single return
! istatus code - 0 is all interpolations worked, positive means one or more
! failed.
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 88:  this kind is not in the state vector
!       ISTATUS = 11:  Could not find a triangle that contains this lat/lon
!       ISTATUS = 12:  Height vertical coordinate out of model range.
!       ISTATUS = 13:  Missing value in interpolation.
!       ISTATUS = 16:  Don't know how to do vertical velocity for now
!       ISTATUS = 17:  Unable to compute pressure values 
!       ISTATUS = 18:  altitude illegal
!       ISTATUS = 19:  could not compute u using RBF code
!       ISTATUS = 101: Internal error; reached end of subroutine without 
!                      finding an applicable case.
!

! Passed variables

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: num_kinds
integer,             intent(in)  :: obs_kinds(:)
real(r8),            intent(out) :: interp_vals(:)
integer,             intent(out) :: istatus

! Local storage

real(r8) :: loc_array(3), llon, llat, lheight, fract, v_interp
integer  :: i, j, ivar, ier, lower, upper, this_kind
integer  :: pt_base_offset, density_base_offset, qv_base_offset

real(r8) :: weights(3), lower_interp, upper_interp, ltemp, utemp
integer  :: tri_indices(3), base_offset(num_kinds)

if ( .not. module_initialized ) call static_init_model

! Assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value.  Make any error codes set here be in the 10s

interp_vals(:) = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! see if all observation kinds are in the state vector.  this sets an
! error code and returns without a fatal error if any answer is no.
! exceptions:  if the vertical location is specified in pressure
! we end up computing the sensible temperature as part of converting
! to height (meters), so if the kind to be computed is sensible temp
! (which is not in state vector; potential temp is), go ahead and
! say yes, we know how to compute it.  if the vert is any other units
! then fail and let the calling code call in for the components
! it needs to compute sensible temp itself.  also, if the obs is a
! wind (or ocean current) component, and the edge normal speeds are
! in the state vector, we can do it.

oktointerp: do i=1, num_kinds
   ivar = get_progvar_index_from_kind(obs_kinds(i))
   if (ivar <= 0) then
      ! exceptions 1 and 2:
      ! FIXME: the new code cannot yet compute sensible temperature in one go.  fail for that.
      ! remove the new code test once we add it.
      if ((obs_kinds(i) == KIND_TEMPERATURE) .and. vert_is_pressure(location) .and. &
          .not. use_new_code) cycle oktointerp
      if ((obs_kinds(i) == KIND_U_WIND_COMPONENT) .or. obs_kinds(i) == KIND_V_WIND_COMPONENT) then
         ivar = get_progvar_index_from_kind(KIND_EDGE_NORMAL_SPEED)
         if (ivar > 0) cycle oktointerp
      endif

      istatus = 88            ! this kind not in state vector
      return
   endif
enddo oktointerp

! Not prepared to do w interpolation at this time
do i=1, num_kinds
   if(obs_kinds(i) == KIND_VERTICAL_VELOCITY) then
      istatus = 16
      return
   endif
enddo

! Get the individual locations values
loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)

if (debug > 5) print *, 'requesting interpolation at ', llon, llat, lheight


!  if (debug > 2) print *, 'base offset now ', base_offset(1)

! FIXME: this needs to play well with multiple types in the kinds list.
! for now do only the first one.  if 'u' is part of the state vector
! and if they are asking for U/V components, use the new RBF code.
! if u isn't in the state vector, default to trying to interpolate
! in the reconstructed U/V components at the cell centers.

if (use_new_code) then
   kindloop: do i=1, num_kinds
      if ((obs_kinds(i) == KIND_U_WIND_COMPONENT .or. &
           obs_kinds(i) == KIND_V_WIND_COMPONENT) .and. has_real_u) then
         if (obs_kinds(i) == KIND_U_WIND_COMPONENT) then
            call compute_u_with_rbf(x, location, .TRUE., interp_vals(i), istatus)
         else
            call compute_u_with_rbf(x, location, .FALSE., interp_vals(i), istatus)
         endif

      else if (obs_kinds(i) /= KIND_TEMPERATURE) then
         ivar = get_progvar_index_from_kind(obs_kinds(i))
         call compute_scalar_with_barycentric(x, location, ivar, interp_vals(i), istatus)
         if (istatus /= 0) return
 
      else if (obs_kinds(i) == KIND_TEMPERATURE) then
         print *, 'need to add case in new code for sensible temperature'
         stop
         ! need to get potential temp, pressure, qv here, but can
         ! use same weights.  does this get pushed down into the
         ! scalar code?
      endif
   enddo kindloop
   return
endif

!! start of original code - should be unneeded if we decide to
!! use the new code.
if (debug > 4) print *, 'using original interpolation code'

! Find the start and end offsets for these fields in the state vector x(:)
do i=1, num_kinds
   if (obs_kinds(i) == KIND_TEMPERATURE) then
      base_offset(i) = 0  ! this won't be used in this case
   else
      call get_index_range(obs_kinds(i), base_offset(i))
   endif
enddo

! Find the indices of the three cell centers that surround this point in
! the horizontal along with the barycentric weights.
call get_cell_indices(llon, llat, tri_indices, weights, ier)
if (debug > 5) write(*, *) 'tri_inds ', tri_indices
if (debug > 5) write(*, *) 'weights ', weights

! If istatus is not zero couldn't find a triangle, fail
if(ier /= 0) then
   istatus = 11
   return
endif

! FIXME: cannot do both surface pressure and any other 3d field with
! this code (yet).  if num kinds > 1 and any are not surface pressure, this
! should fail.
! Surface pressure is a 2D field
if(obs_kinds(1) == KIND_SURFACE_PRESSURE) then
   lheight = 1.0_r8 
   ! the 1 below is max number of vert levels to examine
   call triangle_interp(x, base_offset(1), tri_indices, weights, &
      nint(lheight), 1, interp_vals(1), ier) 
   if(ier /= 0) then
      istatus = 13
   else
      istatus = 0
   endif
   return
endif

! the code below has vert undef returning an error, and vert surface
! returning level 1 (which i believe is correct - the vertical grid
! is terrain following).  FIXME:
! vert undef could return anything in the column in theory, right?

if (vert_is_undef(location)) then
   istatus = 12
   return
endif

if(vert_is_surface(location)) then
   do j=1, num_kinds
      lheight = 1.0_r8  ! first grid level is surface
      call triangle_interp(x, base_offset(j), tri_indices, weights, &
         nint(lheight), nVertLevels, interp_vals(j), ier)
      if(ier /= 0) then
         istatus = 13
         return
      endif
   enddo
   istatus = 0
   return
endif

! If vertical is on a model level don't need vertical interpolation either
! (FIXME: since the vertical value is a real, in the wrf model_mod they
! support non-integer levels and interpolate in the vertical just like
! any other vertical type.  we could easily add that here.)
if(vert_is_level(location)) then
   ! FIXME: if this is W, the top is nVertLevels+1
   if (lheight > nVertLevels) then
      istatus = 12
      return
   endif
   do j=1, num_kinds
      call triangle_interp(x, base_offset(j), tri_indices, weights, &
         nint(lheight), nVertLevels, interp_vals(j), ier)
      if(ier /= 0) then
         istatus = 13
         return
      endif
   enddo
   istatus = 0
   return
endif

! Vertical interpolation for pressure coordinates
  if(vert_is_pressure(location) ) then 
   ! Need to get base offsets for the potential temperature, density, and water 
   ! vapor mixing fields in the state vector
   call get_index_range(KIND_POTENTIAL_TEMPERATURE, pt_base_offset)
   call get_index_range(KIND_DENSITY, density_base_offset)
   call get_index_range(KIND_VAPOR_MIXING_RATIO, qv_base_offset)
!print *, 'bases: t/rho/v = ', pt_base_offset, density_base_offset, qv_base_offset
   call find_pressure_bounds(x, lheight, tri_indices, weights, nVertLevels, &
         pt_base_offset, density_base_offset, qv_base_offset, lower, upper, fract, &
         ltemp, utemp, ier)
   if(ier /= 0) then
      interp_vals = MISSING_R8
      istatus = 17
      return
   endif
   ! if any of the input kinds are sensible temperature, we had to compute
   ! that value already to convert the pressure to the model levels, so just
   ! interpolate in the vertical and we're done with that one.
   do j=1, num_kinds
      if (obs_kinds(j) == KIND_TEMPERATURE) then
         ! Already have both values, interpolate in the vertical here.
         if (ltemp == MISSING_R8 .or. utemp == MISSING_R8) then
            istatus = 12
            return
         endif
         interp_vals(j) = (1.0_r8 - fract) * ltemp + fract * utemp
      else
         ! Next interpolate the observed quantity to the horizontal point at both levels
         call triangle_interp(x, base_offset(j), tri_indices, weights, lower, nVertLevels, &
                              lower_interp, ier)
         if(ier /= 0) then
            istatus = 13
            return
         endif
         call triangle_interp(x, base_offset(j), tri_indices, weights, upper, nVertLevels, &
                              upper_interp, ier)
         if(ier /= 0) then
            istatus = 13
            return
         endif

         ! Got both values, interpolate and return
         if (lower_interp == MISSING_R8 .or. upper_interp == MISSING_R8) then
            istatus = 12
            return
         endif
         interp_vals(j) = (1.0_r8 - fract) * lower_interp + fract * upper_interp
      endif
   enddo
   istatus = 0
   return
endif


! in this section, unlike the others, we loop 3 times adding successive
! contributions to the return value.  if there's an error the 
! result value has been set to 0 so we have to reset it to the 
! missing value instead of only setting the istatus code.
if(vert_is_height(location)) then
   ! For height, can do simple vertical search for interpolation for now
   ! Get the lower and upper bounds and fraction for each column
   interp_vals(:) = 0.0_r8
   do i = 1, 3
      call find_height_bounds(lheight, nVertLevels, zGridCenter(:, tri_indices(i)), &
                              lower, upper, fract, ier)
      if(ier /= 0) then
         istatus = 12
         interp_vals(:) = MISSING_R8
         return
      endif
      do j=1, num_kinds
         call vert_interp(x, base_offset(j), tri_indices(i), nVertLevels, lower, fract, v_interp, ier)
         if(ier /= 0) then
            istatus = 13
            interp_vals(:) = MISSING_R8
            return
         endif
         interp_vals(j) = interp_vals(j) + weights(i) * v_interp
      enddo
   enddo
   istatus = 0
   return
endif

! shouldn't get here.
interp_vals(:) = MISSING_R8
istatus = 101

end subroutine local_interpolate


!------------------------------------------------------------------

function nc_write_model_atts( ncFileID ) result (ierr)

! TJH -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector.
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

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------------------------

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: nCellsDimID
integer :: nEdgesDimID, maxEdgesDimID
integer :: nVerticesDimID
integer :: VertexDegreeDimID
integer :: nSoilLevelsDimID
integer :: nVertLevelsDimID
integer :: nVertLevelsP1DimID


! for the prognostic variables
integer :: ivar, VarID, mpasFileID

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
integer :: myndims

character(len=128) :: filename

real(r8), allocatable, dimension(:)   :: data1d

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
           'nc_write_model_atts', 'creation put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
           'nc_write_model_atts', 'source put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
           'nc_write_model_atts', 'revision put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
           'nc_write_model_atts', 'revdate put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'MPAS_OCN' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

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

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   ! Leave define mode.
   call nc_check(nf90_enddef(ncFileID),'nc_write_model_atts','state enddef '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nCells', &
          len = nCells, dimid = nCellsDimID),'nc_write_model_atts', 'nCells def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nEdges', &
          len = nEdges, dimid = nEdgesDimID),'nc_write_model_atts', 'nEdges def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='maxEdges', &
          len = maxEdges, dimid = maxEdgesDimID),'nc_write_model_atts', 'maxEdges def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertices', &
          len = nVertices, dimid = nVerticesDimID),'nc_write_model_atts', &
               'nVertices def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='VertexDegree', &
          len = VertexDegree, dimid = VertexDegreeDimID),'nc_write_model_atts', &
               'VertexDegree def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertLevels', &
          len = nVertLevels, dimid = NVertLevelsDimID),'nc_write_model_atts', &
                                      'nVertLevels def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertLevelsP1', &
          len = nVertLevelsP1, dimid = NVertLevelsP1DimID),'nc_write_model_atts', &
                                      'nVertLevelsP1 def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nSoilLevels', &
          len = nSoilLevels, dimid = nSoilLevelsDimID),'nc_write_model_atts', &
               'nSoilLevels def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lonCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'lonCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center longitudes'), &
                 'nc_write_model_atts', 'lonCell long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lonCell units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lonCell valid_range '//trim(filename))

   ! Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='latCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'latCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center latitudes'), &
                 'nc_write_model_atts', 'latCell long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'latCell units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'latCell valid_range '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='xCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'xCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center x cartesian coordinates'), &
                 'nc_write_model_atts', 'xCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='yCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'yCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center y cartesian coordinates'), &
                 'nc_write_model_atts', 'yCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='zCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'zCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center z cartesian coordinates'), &
                 'nc_write_model_atts', 'zCell long_name '//trim(filename))

   ! Grid vertical information
   call nc_check(nf90_def_var(ncFileID,name='zgrid',xtype=nf90_double, &
                 dimids=(/ nVertLevelsP1DimID, nCellsDimID /) ,varid=VarID), &
                 'nc_write_model_atts', 'zgrid def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid zgrid'), &
                 'nc_write_model_atts', 'zgrid long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'zgrid units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'zgrid units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'zgrid cartesian_axis '//trim(filename))

   ! Vertex Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lonVertex', xtype=nf90_double, &
                 dimids=nVerticesDimID, varid=VarID),&
                 'nc_write_model_atts', 'lonVertex def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'vertex longitudes'), &
                 'nc_write_model_atts', 'lonVertex long_name '//trim(filename))

   ! Vertex Latitudes
   call nc_check(nf90_def_var(ncFileID,name='latVertex', xtype=nf90_double, &
                 dimids=nVerticesDimID, varid=VarID),&
                 'nc_write_model_atts', 'latVertex def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'vertex latitudes'), &
                 'nc_write_model_atts', 'latVertex long_name '//trim(filename))

   ! Grid relationship information
   call nc_check(nf90_def_var(ncFileID,name='nEdgesOnCell',xtype=nf90_int, &
                 dimids=nCellsDimID ,varid=VarID), &
                 'nc_write_model_atts', 'nEdgesOnCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid nEdgesOnCell'), &
                 'nc_write_model_atts', 'nEdgesOnCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='cellsOnVertex',xtype=nf90_int, &
                 dimids=(/ VertexDegreeDimID, nVerticesDimID /) ,varid=VarID), &
                 'nc_write_model_atts', 'cellsOnVertex def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid cellsOnVertex'), &
                 'nc_write_model_atts', 'cellsOnVertex long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='verticesOnCell',xtype=nf90_int, &
                 dimids=(/ maxEdgesDimID, nCellsDimID /) ,varid=VarID), &
                 'nc_write_model_atts', 'verticesOnCell def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid verticesOnCell'), &
                 'nc_write_model_atts', 'verticesOnCell long_name '//trim(filename))

   call nc_check(nf90_def_var(ncFileID,name='areaCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'areaCell def_var '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ncFileID, ivar, MemberDimID, unlimitedDimID, myndims, mydimids) 

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(string1)//' put_att units' )

   enddo

   !----------------------------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables that DART needs and has locally 
   !----------------------------------------------------------------------------

   call nc_check(NF90_inq_varid(ncFileID, 'lonCell', VarID), &
                 'nc_write_model_atts', 'lonCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, lonCell ), &
                'nc_write_model_atts', 'lonCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'latCell', VarID), &
                 'nc_write_model_atts', 'latCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, latCell ), &
                'nc_write_model_atts', 'latCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'zgrid', VarID), &
                 'nc_write_model_atts', 'zgrid inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, zGridFace ), &
                'nc_write_model_atts', 'zgrid put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'nEdgesOnCell', VarID), &
                 'nc_write_model_atts', 'nEdgesOnCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, nEdgesOnCell ), &
                'nc_write_model_atts', 'nEdgesOnCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'verticesOnCell', VarID), &
                 'nc_write_model_atts', 'verticesOnCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, verticesOnCell ), &
                'nc_write_model_atts', 'verticesOnCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'cellsOnVertex', VarID), &
                 'nc_write_model_atts', 'cellsOnVertex inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, cellsOnVertex ), &
                'nc_write_model_atts', 'cellsOnVertex put_var '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables needed for plotting only.
   ! DART has not read these in, so we have to read them from the input file
   ! and parrot them to the DART output file. 
   !----------------------------------------------------------------------------

   call nc_check(nf90_open(trim(grid_definition_filename), NF90_NOWRITE, mpasFileID), &
                 'nc_write_model_atts','open '//trim(grid_definition_filename))

   allocate(data1d(nCells))
   call nc_check(nf90_inq_varid(mpasFileID, 'xCell', VarID), &
                 'nc_write_model_atts',     'xCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'xCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'xCell', VarID), &
                 'nc_write_model_atts',     'xCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'xCell put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'yCell', VarID), &
                 'nc_write_model_atts',     'yCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'yCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'yCell', VarID), &
                 'nc_write_model_atts',     'yCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'yCell put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'zCell', VarID), &
                 'nc_write_model_atts',     'zCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'zCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'zCell', VarID), &
                 'nc_write_model_atts',     'zCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'zCell put_var '//trim(filename))

   call nc_check(nf90_inq_varid(mpasFileID, 'areaCell', VarID), &
                 'nc_write_model_atts',     'areaCell inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'areaCell get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'areaCell', VarID), &
                 'nc_write_model_atts',     'areaCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'areaCell put_var '//trim(filename))
   deallocate(data1d)

   allocate(data1d(nVertices))
   call nc_check(nf90_inq_varid(mpasFileID, 'lonVertex', VarID), &
                 'nc_write_model_atts',     'lonVertex inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'lonVertex get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'lonVertex', VarID), &
                 'nc_write_model_atts',     'lonVertex inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'lonVertex put_var '//trim(filename))
   
   call nc_check(nf90_inq_varid(mpasFileID, 'latVertex', VarID), &
                 'nc_write_model_atts',     'latVertex inq_varid ')
   call nc_check(nf90_get_var(mpasFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'latVertex get_var ')
   call nc_check(nf90_inq_varid(ncFileID,   'latVertex', VarID), &
                 'nc_write_model_atts',     'latVertex inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, data1d ), &
                 'nc_write_model_atts',     'latVertex put_var '//trim(filename))
   deallocate(data1d)

   call nc_check(nf90_close(mpasFileID),'nc_write_model_atts','close '//trim(grid_definition_filename))
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!------------------------------------------------------------------

function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)         

! TJH 29 Aug 2011 -- all errors are fatal, so the
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

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: state_vec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname 
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields  

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      mystart = 1   ! These are arrays, actually
      mycount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

     ! FIXME - wouldn't hurt to make sure each of these match something.
     !         could then eliminate the if ncndims /= xxx checks below.

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if ( debug > 1 ) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',mystart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',mycount(1:ncNdims)
      endif

      if (     progvar(ivar)%numdims == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ) )
         call vector_to_prog_var(state_vec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_1d_array)

      elseif ( progvar(ivar)%numdims == 2 ) then

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2) ))
         call vector_to_prog_var(state_vec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      elseif ( progvar(ivar)%numdims == 3) then

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_3d_array( progvar(ivar)%dimlens(1), &
                                 progvar(ivar)%dimlens(2), &
                                 progvar(ivar)%dimlens(3)))
         call vector_to_prog_var(state_vec, ivar, data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      else

         ! FIXME put an error message here
         write(string1,*)'no support (yet) for 4d fields'
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate)

      endif

   enddo


endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!------------------------------------------------------------------

function get_model_size()

! Returns the size of the model as an integer. 
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------

function get_model_time_step()

! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!------------------------------------------------------------------

subroutine ens_mean_for_model(filter_ens_mean)

! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles.

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model


!------------------------------------------------------------------

subroutine end_model()

! Does any shutdown and clean-up needed for model.

if (allocated(latCell))        deallocate(latCell)
if (allocated(lonCell))        deallocate(lonCell)
if (allocated(zGridFace))      deallocate(zGridFace)
if (allocated(zGridCenter))    deallocate(zGridCenter)
if (allocated(cellsOnVertex))  deallocate(cellsOnVertex)

end subroutine end_model


!------------------------------------------------------------------

subroutine pert_model_state(state, pert_state, interf_provided)

! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

real(r8)              :: pert_ampl
real(r8)              :: minv, maxv, temp
type(random_seq_type) :: random_seq
integer               :: i, j, s, e
integer, save         :: counter = 1


! generally you do not want to perturb a single state
! to begin an experiment - unless you make minor perturbations
! and then run the model free for long enough that differences
! develop which contain actual structure.
!
! the subsequent code is a pert routine which
! can be used to add minor perturbations which can be spun up.
!
! if all values in a field are identical (i.e. 0.0) this 
! routine will not change those values since it won't 
! make a new value outside the original min/max of that
! variable in the state vector.  to handle this case you can
! remove the min/max limit lines below.


! start of pert code

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! the first time through get the task id (0:N-1)
! and set a unique seed per task.  this won't
! be consistent between different numbers of mpi
! tasks, but at least it will reproduce with
! multiple runs with the same task count.
! best i can do since this routine doesn't have
! the ensemble member number as an argument
! (which i think it needs for consistent seeds).
!
! this only executes the first time since counter
! gets incremented after the first use and the value
! is saved between calls.
if (counter == 1) counter = counter + (my_task_id() * 1000)

call init_random_seq(random_seq, counter)
counter = counter + 1

do i=1, nfields
   ! starting and ending indices in the linear state vect
   ! for each different state kind.
   s = progvar(i)%index1
   e = progvar(i)%indexN
   ! original min/max data values of each type
   minv = minval(state(s:e))
   maxv = maxval(state(s:e))
   do j=s, e
      ! once you change pert_state, state is changed as well
      ! since they are the same storage as called from filter.
      ! you have to save it if you want to use it again.
      temp = state(j)  ! original value
      ! perturb each value individually
      ! make the perturbation amplitude N% of this value
      pert_ampl = model_perturbation_amplitude * temp
      pert_state(j) = random_gaussian(random_seq, state(j), pert_ampl)
      ! keep it from exceeding the original range
      pert_state(j) = max(minv, pert_state(j))
      pert_state(j) = min(maxv, pert_state(j))
   enddo
enddo

end subroutine pert_model_state


!------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, num_close, close_ind, dist)
!
! FIXME ... not tested.
!
! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

type(get_close_type),              intent(in)    :: gc
type(location_type),               intent(inout) :: base_obs_loc
integer,                           intent(in)    :: base_obs_kind
type(location_type), dimension(:), intent(inout) :: obs_loc
integer,             dimension(:), intent(in)    :: obs_kind
integer,                           intent(out)   :: num_close
integer,             dimension(:), intent(out)   :: close_ind
real(r8),            dimension(:), intent(out)   :: dist

integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc

! Initialize variables to missing status

num_close = 0
close_ind = -99
dist      = 1.0e9   !something big and positive (far away) in radians
istatus1  = 0
istatus2  = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_array = get_location(base_obs_loc)
base_which = nint(query_location(base_obs_loc))

! fixme ... 
if (.not. horiz_dist_only) then
!  if (base_which /= wrf%dom(1)%vert_coord) then
!     call vert_interpolate(ens_mean, base_obs_loc, base_obs_kind, istatus1)
!  elseif (base_array(3) == MISSING_R8) then
!     istatus1 = 1
!  endif
endif

if (istatus1 == 0) then

   ! Loop over potentially close subset of obs priors or state variables
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for obs_loc).
   call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind)

   do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = obs_loc(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
      ! This should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
      if (.not. horiz_dist_only) then
 !fixme       if (local_obs_which /= wrf%dom(1)%vert_coord) then
 !fixme           call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
            ! Store the "new" location into the original full local array
            obs_loc(t_ind) = local_obs_loc
 !fixme        endif
      endif

      ! Compute distance - set distance to a very large value if vert coordinate is missing
      ! or vert_interpolate returned error (istatus2=1)
      local_obs_array = get_location(local_obs_loc)
      if (( (.not. horiz_dist_only)             .and. &
            (local_obs_array(3) == MISSING_R8)) .or.  &
            (istatus2 == 1)                   ) then
            dist(k) = 1.0e9
      else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
      endif

   enddo
endif

end subroutine get_close_obs


!------------------------------------------------------------------

subroutine init_time(time)

! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

! this shuts up the compiler warnings about unused variables
time = set_time(0, 0)

write(string1,*) 'Cannot initialize MPAS time via subroutine call; start_from_restart cannot be F'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

end subroutine init_time


!------------------------------------------------------------------

subroutine init_conditions(x)

! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model
 
! this shuts up the compiler warnings about unused variables
x = 0.0_r8

write(string1,*) 'Cannot initialize MPAS state via subroutine call; start_from_restart cannot be F'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

end subroutine init_conditions


!------------------------------------------------------------------

subroutine adv_1step(x, time)

! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called IF the namelist parameter
! async is set to 0 in perfect_model_obs or filter -OR- if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If none of these options
! are used (the model will only be advanced as a separate 
! model-specific executable), this can be a NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*) 'Cannot advance MPAS with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

end subroutine adv_1step



!==================================================================
! The (model-specific) additional public interfaces come next
!  (these are not required by dart but are used by other programs)
!==================================================================


subroutine get_model_analysis_filename( filename )

! return the name of the analysis filename that was set
! in the model_nml namelist

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(model_analysis_filename)

end subroutine get_model_analysis_filename


!-------------------------------------------------------------------

subroutine get_grid_definition_filename( filename )

! return the name of the grid_definition filename that was set
! in the model_nml namelist

character(len=*), intent(out) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(grid_definition_filename)

end subroutine get_grid_definition_filename


!-------------------------------------------------------------------

subroutine analysis_file_to_statevector(filename, state_vector, model_time)

! Reads the current time and state variables from a mpas analysis
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: model_time

! temp space to hold data while we are reading it
integer  :: ndim1, ndim2, ndim3
integer  :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: ncid, TimeDimID, TimeDimLength
character(len=256) :: myerrorstring

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'analysis_file_to_statevector',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
             'analysis_file_to_statevector','open '//trim(filename))

model_time = get_analysis_time(ncid, filename)

if (do_output()) &
    call print_time(model_time,'time in restart file '//trim(filename))
if (do_output()) &
    call print_date(model_time,'date in restart file '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

TimeDimID = FindTimeDimension( ncid )

if ( TimeDimID > 0 ) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=TimeDimLength), &
            'analysis_file_to_statevector', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(filename)//' '//trim(varname)

   ! determine the shape of the netCDF variable

   call nc_check(nf90_inq_varid(ncid,   varname, VarID), &
            'analysis_file_to_statevector', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'analysis_file_to_statevector', 'inquire '//trim(myerrorstring))

   mystart = 1   ! These are arrays, actually.
   mycount = 1

   ! Only checking the shape of the variable - sans TIME
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'analysis_file_to_statevector', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'analysis_file_to_statevector',string1,source,revision,revdate)
      endif

      mycount(i) = dimlen

   enddo DimCheck

   where(dimIDs == TimeDimID) mystart = TimeDimLength  ! pick the latest time
   where(dimIDs == TimeDimID) mycount = 1              ! only use one time

   if ( debug > 1 ) then
      write(*,*)'analysis_file_to_statevector '//trim(varname)//' start = ',mystart(1:ncNdims)
      write(*,*)'analysis_file_to_statevector '//trim(varname)//' count = ',mycount(1:ncNdims)
   endif

   if (ncNdims == 1) then

      ! If the single dimension is TIME, we only need a scalar.
      ! Pretty sure this cannot happen ...
      ndim1 = mycount(1)
      allocate(data_1d_array(ndim1))
      call nc_check(nf90_get_var(ncid, VarID, data_1d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      write(string1, *) 'data min/max ', trim(varname), minval(data_1d_array), maxval(data_1d_array)
      call error_handler(E_MSG, '', string1, &
                        source,revision,revdate)

      call prog_var_to_vector(data_1d_array, state_vector, ivar)
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      allocate(data_2d_array(ndim1, ndim2))
      call nc_check(nf90_get_var(ncid, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      write(string1, *) 'data min/max ', trim(varname), minval(data_2d_array), maxval(data_2d_array)
      call error_handler(E_MSG, '', string1, &
                        source,revision,revdate)

      call prog_var_to_vector(data_2d_array, state_vector, ivar)
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ndim1 = mycount(1)
      ndim2 = mycount(2)
      ndim3 = mycount(3)
      allocate(data_3d_array(ndim1, ndim2, ndim3))
      call nc_check(nf90_get_var(ncid, VarID, data_3d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'analysis_file_to_statevector', 'get_var '//trim(varname))

      write(string1, *) 'data min/max ', trim(varname), minval(data_3d_array), maxval(data_3d_array)
      call error_handler(E_MSG, '', string1, &
                        source,revision,revdate)

      call prog_var_to_vector(data_3d_array, state_vector, ivar)
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'analysis_file_to_statevector', string1, &
                        source,revision,revdate)
   endif

enddo

call nc_check(nf90_close(ncid), &
             'analysis_file_to_statevector','close '//trim(filename))

end subroutine analysis_file_to_statevector


!-------------------------------------------------------------------

subroutine statevector_to_analysis_file(state_vector, filename, statetime)

! Writes the current time and state variables from a dart state
! vector (1d array) into a mpas netcdf analysis file.

real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statetime

! temp space to hold data while we are writing it
integer :: i, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array2
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array

! temp space to hold horizontal winds
real(r8), allocatable :: u(:,:), ucell_incr(:,:), vcell_incr(:,:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: varname
integer :: VarID, ncNdims, dimlen
integer :: zonal, meridional
integer :: ncFileID, TimeDimID, TimeDimLength
logical :: both
type(time_type) :: model_time

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for writing.'
   call error_handler(E_ERR,'statevector_to_analysis_file',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncFileID), &
             'statevector_to_analysis_file','open '//trim(filename))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the mpas analysis file, and state vector contents from a different
! time won't be consistent with the rest of the file.

model_time = get_analysis_time(ncFileID, filename)

if ( model_time /= statetime ) then
   call print_time( statetime,'DART current time',logfileunit)
   call print_time(model_time,'mpas current time',logfileunit)
   call print_time( statetime,'DART current time')
   call print_time(model_time,'mpas current time')
   write(string1,*)trim(filename),' current time must equal model time'
   call error_handler(E_ERR,'statevector_to_analysis_file',string1,source,revision,revdate)
endif

if (do_output()) &
    call print_time(statetime,'time of DART file '//trim(filename))
if (do_output()) &
    call print_date(statetime,'date of DART file '//trim(filename))

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

TimeDimID = FindTimeDimension( ncFileID )

if ( TimeDimID > 0 ) then
   call nc_check(nf90_inquire_dimension(ncFileID, TimeDimID, len=TimeDimLength), &
            'statevector_to_analysis_file', 'inquire timedimlength '//trim(filename))
else
   TimeDimLength = 0
endif

PROGVARLOOP : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   if ( varname == 'uReconstructZonal' .or. &
        varname == 'uReconstructMeridional' ) then
      call update_wind_components(ncFileID, state_vector)
      cycle PROGVARLOOP
   endif

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'statevector_to_analysis_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'statevector_to_analysis_file', 'inquire '//trim(string2))

   mystart = 1   ! These are arrays, actually.
   mycount = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'statevector_to_analysis_file', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'statevector_to_analysis_file', string1, &
                         source, revision, revdate, text2=string2)
      endif

      mycount(i) = dimlen

   enddo DimCheck


   where(dimIDs == TimeDimID) mystart = TimeDimLength
   where(dimIDs == TimeDimID) mycount = 1   ! only the latest one

   if ( debug > 1 ) then
      write(*,*)'statevector_to_analysis_file '//trim(varname)//' start is ',mystart(1:ncNdims)
      write(*,*)'statevector_to_analysis_file '//trim(varname)//' count is ',mycount(1:ncNdims)
   endif


   if (progvar(ivar)%numdims == 1) then
      allocate(data_1d_array(mycount(1)))
      call vector_to_prog_var(state_vector, ivar, data_1d_array)

      write(string1, *) 'data min/max ', trim(varname), minval(data_1d_array), maxval(data_1d_array)
      call error_handler(E_MSG, '', string1, source,revision,revdate)

      if ( progvar(ivar)%clamping ) then
         where ( data_1d_array < progvar(ivar)%range(1) ) data_1d_array = progvar(ivar)%range(1)
         where ( data_1d_array > progvar(ivar)%range(2) ) data_1d_array = progvar(ivar)%range(2)

         write(string1, *) 'after clamping min/max ', trim(varname), &
                            minval(data_1d_array), maxval(data_1d_array)
         call error_handler(E_MSG, '', string1, source,revision,revdate)

      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
            start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      allocate(data_2d_array(mycount(1), mycount(2)))
      call vector_to_prog_var(state_vector, ivar, data_2d_array)

      write(string1, *) 'data min/max ', trim(varname), minval(data_2d_array), maxval(data_2d_array)
      call error_handler(E_MSG, '', trim(string1), source,revision,revdate)

      if ( progvar(ivar)%clamping ) then
         where ( data_2d_array < progvar(ivar)%range(1) ) data_2d_array = progvar(ivar)%range(1)
         where ( data_2d_array > progvar(ivar)%range(2) ) data_2d_array = progvar(ivar)%range(2)

         write(string1, *) 'after clamping min/max ', trim(varname), &
                            minval(data_2d_array), maxval(data_2d_array)
         call error_handler(E_MSG, '', string1, source,revision,revdate)

      endif


      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   elseif (progvar(ivar)%numdims == 3) then

      allocate(data_3d_array(mycount(1), mycount(2), mycount(3)))
      call vector_to_prog_var(state_vector, ivar, data_3d_array)

      write(string1, *) 'data min/max ', trim(varname), minval(data_3d_array), maxval(data_3d_array)
      call error_handler(E_MSG, '', string1, source,revision,revdate)

      if ( progvar(ivar)%clamping ) then
         where ( data_3d_array < progvar(ivar)%range(1) ) data_3d_array = progvar(ivar)%range(1)
         where ( data_3d_array > progvar(ivar)%range(2) ) data_3d_array = progvar(ivar)%range(2)

         write(string1, *) 'after clamping min/max ', trim(varname), &
                            minval(data_3d_array), maxval(data_3d_array)
         call error_handler(E_MSG, '', string1, source,revision,revdate)

      endif

      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
        start=mystart(1:ncNdims), count=mycount(1:ncNdims)), &
            'statevector_to_analysis_file', 'put_var '//trim(varname))
      deallocate(data_3d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'statevector_to_analysis_file', string1, &
                        source,revision,revdate)
   endif

enddo PROGVARLOOP 

call nc_check(nf90_close(ncFileID), &
             'statevector_to_analysis_file','close '//trim(filename))

end subroutine statevector_to_analysis_file


!------------------------------------------------------------------

function get_analysis_time_ncid( ncid, filename )

! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename
type(time_type) :: get_analysis_time_ncid

! local variables
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, idims
integer           :: VarID, numdims

character(len=64) :: timestring

if ( .not. module_initialized ) call static_init_model

call nc_check( nf90_inq_varid(ncid, 'xtime', VarID), &
              'get_analysis_time', 'inquire xtime '//trim(filename))

call nc_check( nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'get_analysis_time', 'inquire TIME '//trim(filename))

if (numdims /= 2) then
   write(string1,*) 'xtime variable has unknown shape in ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

call nc_check( nf90_inquire_dimension(ncid, dimIDs(1), len=idims(1)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))
call nc_check( nf90_inquire_dimension(ncid, dimIDs(2), len=idims(2)), &
                 'get_analysis_time', 'inquire time dimension length '//trim(filename))

if (idims(2) /= 1) then
   write(string1,*) 'multiple timesteps (',idims(2),') in file ', trim(filename)
   write(string2,*) 'We are using the LAST one, presumably, the LATEST timestep.'
   call error_handler(E_MSG,'get_analysis_time',string1,source,revision,revdate,text2=string2)
endif

! Get the highest ranking time ... the last one, basically.

call nc_check( nf90_get_var(ncid, VarID, timestring, start = (/ 1, idims(2) /)), &
              'get_analysis_time', 'get_var xtime '//trim(filename))

get_analysis_time_ncid = string_to_time(timestring)

if (debug > 5) then
   call print_date(get_analysis_time_ncid, 'get_analysis_time:model date')
   call print_time(get_analysis_time_ncid, 'get_analysis_time:model time')
endif

end function get_analysis_time_ncid


!------------------------------------------------------------------

function get_analysis_time_fname(filename)

! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_analysis_time_fname

character(len=*), intent(in) :: filename

integer :: i

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif

! find the first number and use that as the start of the string conversion
i = scan(filename, "0123456789")
if (i <= 0) then
   write(string1,*) 'cannot find time string in name ', trim(filename)
   call error_handler(E_ERR,'get_analysis_time',string1,source,revision,revdate)
endif 

get_analysis_time_fname = string_to_time(filename(i:i+19))

end function get_analysis_time_fname


!------------------------------------------------------------------

subroutine write_model_time(time_filename, model_time, adv_to_time)
 character(len=*), intent(in)           :: time_filename
 type(time_type),  intent(in)           :: model_time
 type(time_type),  intent(in), optional :: adv_to_time

integer :: iunit
character(len=19) :: timestring
type(time_type)   :: deltatime

iunit = open_file(time_filename, action='write')

timestring = time_to_string(model_time)
write(iunit, '(A)') timestring

if (present(adv_to_time)) then
   timestring = time_to_string(adv_to_time)
   write(iunit, '(A)') timestring

   deltatime = adv_to_time - model_time
   timestring = time_to_string(deltatime, interval=.true.)
   write(iunit, '(A)') timestring
endif

call close_file(iunit)

end subroutine write_model_time


!------------------------------------------------------------------

subroutine get_grid_dims(Cells, Vertices, Edges, VertLevels, VertexDeg, SoilLevels)

! public routine for returning the counts of various things in the grid
!

integer, intent(out) :: Cells         ! Total number of cells making up the grid
integer, intent(out) :: Vertices      ! Unique points in grid which are corners of cells
integer, intent(out) :: Edges         ! Straight lines between vertices making up cells
integer, intent(out) :: VertLevels    ! Vertical levels; count of vert cell centers
integer, intent(out) :: VertexDeg     ! Max number of edges that touch any vertex
integer, intent(out) :: SoilLevels    ! Number of soil layers

if ( .not. module_initialized ) call static_init_model

Cells      = nCells
Vertices   = nVertices
Edges      = nEdges
VertLevels = nVertLevels
VertexDeg  = vertexDegree
SoilLevels = nSoilLevels

end subroutine get_grid_dims


!==================================================================
! The (model-specific) private interfaces come last
!==================================================================


!------------------------------------------------------------------

function time_to_string(t, interval)

! convert time type into a character string with the
! format of YYYY-MM-DD_hh:mm:ss

! passed variables
 character(len=19) :: time_to_string
 type(time_type), intent(in) :: t
 logical, intent(in), optional :: interval

! local variables

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: ndays, nsecs
logical :: dointerval

if (present(interval)) then
   dointerval = interval
else
   dointerval = .false.
endif

! for interval output, output the number of days, then hours, mins, secs
! for date output, use the calendar routine to get the year/month/day hour:min:sec
if (dointerval) then
   call get_time(t, nsecs, ndays)
   if (ndays > 99) then
      write(string1, *) 'interval number of days is ', ndays
      call error_handler(E_ERR,'time_to_string', 'interval days cannot be > 99', &
                         source, revision, revdate, text2=string1)
   endif
   ihour = nsecs / 3600
   nsecs = nsecs - (ihour * 3600)
   imin  = nsecs / 60
   nsecs = nsecs - (imin * 60)
   isec  = nsecs
   write(time_to_string, '(I2.2,3(A1,I2.2))') &
                        ndays, '_', ihour, ':', imin, ':', isec
else
   call get_date(t, iyear, imonth, iday, ihour, imin, isec)
   write(time_to_string, '(I4.4,5(A1,I2.2))') &
                        iyear, '-', imonth, '-', iday, '_', ihour, ':', imin, ':', isec
endif

end function time_to_string


!------------------------------------------------------------------

function string_to_time(s)

! parse a string to extract time.  the expected format of
! the string is YYYY-MM-DD_hh:mm:ss  (although the exact
! non-numeric separator chars are skipped and not validated.)

 type(time_type) :: string_to_time
 character(len=*), intent(in) :: s

integer :: iyear, imonth, iday, ihour, imin, isec

read(s,'(i4,5(1x,i2))') iyear, imonth, iday, ihour, imin, isec
string_to_time = set_date(iyear, imonth, iday, ihour, imin, isec)

end function string_to_time


!------------------------------------------------------------------

function set_model_time_step()

! the static_init_model ensures that the model namelists are read.

type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! these are from the namelist
!FIXME: sanity check these for valid ranges?
set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!------------------------------------------------------------------

subroutine read_grid_dims()

! Read the grid dimensions from the MPAS netcdf file.
!
! The file name comes from module storage ... namelist.

integer :: grid_id, dimid

if ( .not. module_initialized ) call static_init_model

! get the ball rolling ...

call nc_check( nf90_open(trim(grid_definition_filename), NF90_NOWRITE, grid_id), &
              'read_grid_dims', 'open '//trim(grid_definition_filename))

! nCells : get dimid for 'nCells' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nCells', dimid), &
              'read_grid_dims','inq_dimid nCells '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nCells), &
            'read_grid_dims','inquire_dimension nCells '//trim(grid_definition_filename))

! nVertices : get dimid for 'nVertices' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertices', dimid), &
              'read_grid_dims','inq_dimid nVertices '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertices), &
            'read_grid_dims','inquire_dimension nVertices '//trim(grid_definition_filename))

! nEdges : get dimid for 'nEdges' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nEdges', dimid), &
              'read_grid_dims','inq_dimid nEdges '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nEdges), &
            'read_grid_dims','inquire_dimension nEdges '//trim(grid_definition_filename))

! maxEdges : get dimid for 'maxEdges' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'maxEdges', dimid), &
              'read_grid_dims','inq_dimid maxEdges '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=maxEdges), &
            'read_grid_dims','inquire_dimension maxEdges '//trim(grid_definition_filename))

! nVertLevels : get dimid for 'nVertLevels' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertLevels', dimid), &
              'read_grid_dims','inq_dimid nVertLevels '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertLevels), &
            'read_grid_dims','inquire_dimension nVertLevels '//trim(grid_definition_filename))

! nVertLevelsP1 : get dimid for 'nVertLevelsP1' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertLevelsP1', dimid), &
              'read_grid_dims','inq_dimid nVertLevelsP1 '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertLevelsP1), &
            'read_grid_dims','inquire_dimension nVertLevelsP1 '//trim(grid_definition_filename))

! vertexDegree : get dimid for 'vertexDegree' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'vertexDegree', dimid), &
              'read_grid_dims','inq_dimid vertexDegree '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=vertexDegree), &
            'read_grid_dims','inquire_dimension vertexDegree '//trim(grid_definition_filename))

! nSoilLevels : get dimid for 'nSoilLevels' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nSoilLevels', dimid), &
              'read_grid_dims','inq_dimid nSoilLevels '//trim(grid_definition_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nSoilLevels), &
            'read_grid_dims','inquire_dimension nSoilLevels '//trim(grid_definition_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'read_grid_dims','close '//trim(grid_definition_filename) )

if (debug > 7) then
   write(*,*)
   write(*,*)'read_grid_dims: nCells        is ', nCells
   write(*,*)'read_grid_dims: nVertices     is ', nVertices
   write(*,*)'read_grid_dims: nEdges        is ', nEdges
   write(*,*)'read_grid_dims: maxEdges      is ', maxEdges
   write(*,*)'read_grid_dims: nVertLevels   is ', nVertLevels
   write(*,*)'read_grid_dims: nVertLevelsP1 is ', nVertLevelsP1
   write(*,*)'read_grid_dims: vertexDegree  is ', vertexDegree
   write(*,*)'read_grid_dims: nSoilLevels   is ', nSoilLevels
endif

end subroutine read_grid_dims


!------------------------------------------------------------------

subroutine get_grid()

! Read the actual grid values in from the MPAS netcdf file.
!
! The file name comes from module storage ... namelist.
! This reads in the following arrays:
!   latCell, lonCell, zGridFace, cellsOnVertex (all in module global storage)


integer  :: ncid, VarID

! Read the netcdf file data

call nc_check(nf90_open(trim(grid_definition_filename), nf90_nowrite, ncid), 'get_grid', 'open '//trim(grid_definition_filename))

! Read the variables

call nc_check(nf90_inq_varid(ncid, 'latCell', VarID), &
      'get_grid', 'inq_varid latCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, latCell), &
      'get_grid', 'get_var latCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lonCell', VarID), &
      'get_grid', 'inq_varid lonCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, lonCell), &
      'get_grid', 'get_var lonCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'zgrid', VarID), &
      'get_grid', 'inq_varid zgrid '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, zGridFace), &
      'get_grid', 'get_var zgrid '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'cellsOnVertex', VarID), &
      'get_grid', 'inq_varid cellsOnVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, cellsOnVertex), &
      'get_grid', 'get_var cellsOnVertex '//trim(grid_definition_filename))

! MPAS analysis files are in radians - at this point DART needs degrees.

latCell = latCell * rad2deg
lonCell = lonCell * rad2deg

! Read the variables

call nc_check(nf90_inq_varid(ncid, 'edgeNormalVectors', VarID), &
      'get_grid', 'inq_varid edgeNormalVectors '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, edgeNormalVectors), &
      'get_grid', 'get_var edgeNormalVectors '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'nEdgesOnCell', VarID), &
      'get_grid', 'inq_varid nEdgesOnCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, nEdgesOnCell), &
      'get_grid', 'get_var nEdgesOnCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'edgesOnCell', VarID), &
      'get_grid', 'inq_varid edgesOnCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, edgesOnCell), &
      'get_grid', 'get_var edgesOnCell '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'cellsOnEdge', VarID), &
      'get_grid', 'inq_varid cellsOnEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, cellsOnEdge), &
      'get_grid', 'get_var cellsOnEdge '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'latEdge', VarID), &
      'get_grid', 'inq_varid latEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, latEdge), &
      'get_grid', 'get_var latEdge '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'lonEdge', VarID), &
      'get_grid', 'inq_varid lonEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, lonEdge), &
      'get_grid', 'get_var lonEdge '//trim(grid_definition_filename))

latEdge = latEdge * rad2deg
lonEdge = lonEdge * rad2deg

call nc_check(nf90_inq_varid(ncid, 'xVertex', VarID), &
      'get_grid', 'inq_varid xVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, xVertex), &
      'get_grid', 'get_var xVertex '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'yVertex', VarID), &
      'get_grid', 'inq_varid yVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, yVertex), &
      'get_grid', 'get_var yVertex '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'zVertex', VarID), &
      'get_grid', 'inq_varid zVertex '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, zVertex), &
      'get_grid', 'get_var zVertex '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'xEdge', VarID), &
      'get_grid', 'inq_varid xEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, xEdge), &
      'get_grid', 'get_var xEdge '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'yEdge', VarID), &
      'get_grid', 'inq_varid yEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, yEdge), &
      'get_grid', 'get_var yEdge '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'zEdge', VarID), &
      'get_grid', 'inq_varid zEdge '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, zEdge), &
      'get_grid', 'get_var zEdge '//trim(grid_definition_filename))

call nc_check(nf90_inq_varid(ncid, 'verticesOnCell', VarID), &
      'get_grid', 'inq_varid verticesOnCell '//trim(grid_definition_filename))
call nc_check(nf90_get_var( ncid, VarID, verticesOnCell), &
      'get_grid', 'get_var verticesOnCell '//trim(grid_definition_filename))

! Get the boundary information if available. 
! Assuming the existence of this variable is sufficient to determine if
! the grid is defined everywhere or not.

if ( nf90_inq_varid(ncid, 'boundaryEdge', VarID) == NF90_NOERR ) then
   allocate(boundaryEdge(nVertLevels,nEdges))
   call nc_check(nf90_get_var( ncid, VarID, boundaryEdge), &
      'get_grid', 'get_var boundaryEdge '//trim(grid_definition_filename))
   global_grid = .false.
endif

if ( nf90_inq_varid(ncid, 'boundaryVertex', VarID) == NF90_NOERR ) then
   allocate(boundaryVertex(nVertLevels,nVertices))
   call nc_check(nf90_get_var( ncid, VarID, boundaryVertex), &
      'get_grid', 'get_var boundaryVertex '//trim(grid_definition_filename))
   global_grid = .false.
endif

if ( nf90_inq_varid(ncid, 'maxLevelCell', VarID) == NF90_NOERR ) then
   allocate(maxLevelCell(nCells))
   call nc_check(nf90_get_var( ncid, VarID, maxLevelCell), &
      'get_grid', 'get_var maxLevelCell '//trim(grid_definition_filename))
   all_levels_exist_everywhere = .false.
endif

call nc_check(nf90_close(ncid), 'get_grid','close '//trim(grid_definition_filename) )

! A little sanity check

if ( debug > 7 ) then

   write(*,*)
   write(*,*)'latCell           range ',minval(latCell),           maxval(latCell)
   write(*,*)'lonCell           range ',minval(lonCell),           maxval(lonCell)
   write(*,*)'zgrid             range ',minval(zGridFace),         maxval(zGridFace)
   write(*,*)'cellsOnVertex     range ',minval(cellsOnVertex),     maxval(cellsOnVertex)
   write(*,*)'edgeNormalVectors range ',minval(edgeNormalVectors), maxval(edgeNormalVectors)
   write(*,*)'nEdgesOnCell      range ',minval(nEdgesOnCell),      maxval(nEdgesOnCell)
   write(*,*)'EdgesOnCell       range ',minval(EdgesOnCell),       maxval(EdgesOnCell)
   write(*,*)'cellsOnEdge       range ',minval(cellsOnEdge),       maxval(cellsOnEdge)
   write(*,*)'latEdge           range ',minval(latEdge),           maxval(latEdge)
   write(*,*)'lonEdge           range ',minval(lonEdge),           maxval(lonEdge)
   write(*,*)'xVertex           range ',minval(xVertex),           maxval(xVertex)
   write(*,*)'yVertex           range ',minval(yVertex),           maxval(yVertex)
   write(*,*)'zVertex           range ',minval(zVertex),           maxval(zVertex)
   write(*,*)'xEdge             range ',minval(xEdge),             maxval(xEdge)
   write(*,*)'yEdge             range ',minval(yEdge),             maxval(yEdge)
   write(*,*)'zEdge             range ',minval(zEdge),             maxval(zEdge)
   write(*,*)'verticesOnCell    range ',minval(verticesOnCell),    maxval(verticesOnCell)
   if (allocated(boundaryEdge)) &
   write(*,*)'boundaryEdge      range ',minval(boundaryEdge),      maxval(boundaryEdge)
   if (allocated(boundaryVertex)) &
   write(*,*)'boundaryVertex    range ',minval(boundaryVertex),    maxval(boundaryVertex)
   if (allocated(maxLevelCell)) &
   write(*,*)'maxLevelCell      range ',minval(maxLevelCell),      maxval(maxLevelCell)

endif

end subroutine get_grid


!------------------------------------------------------------------


subroutine update_wind_components( ncid, state_vector )

! Special processing if the DART state uses the 'Reconstructed' winds (on the grid centers)
! as opposed to the components on the grid edge centers (components normal to and parallel 
! with the edge direction).  We can READ reconstructed winds directly from the analysis file,
! but we must UPDATE the grid edge arrays as well as the centers.
! 
! SYHA: We should update normal velocity even when one of the horizontal wind components
! exists because we no longer compute the full normal velocity from both u and v winds,
! but only update the original field with the analysis increments.
! Get the contents of the U, uReconstructZonal, uReconstructMeridional arrays.
! We need to read them to compute their increments.

! must first be allocated by calling code with the following sizes:
integer,  intent(in)  :: ncid                  ! netCDF handle for model_analysis_filename
real(r8), intent(in)  :: state_vector(:)

real(r8), allocatable :: u(:,:)                ! u(nVertLevels, nEdges) 
real(r8), allocatable :: ucell_incr(:,:)       ! uReconstructZonal(nVertLevels, nCells) 
real(r8), allocatable :: vcell_incr(:,:)       ! uReconstructMeridional(nVertLevels, nCells) 
real(r8), allocatable :: data_2d_array(:,:)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME) :: dimname 
integer :: VarID, numdims, dimlen, ntimes, i
integer :: zonal, meridional
logical :: both
logical :: already_updated = .false.

if ( .not. module_initialized ) call static_init_model

if ( already_updated ) return  ! only need to do this routine one time

allocate(         u(nVertLevels, nEdges))
allocate(ucell_incr(nVertLevels, nCells))
allocate(vcell_incr(nVertLevels, nCells))

!
! Read 'u' : the normal component of the wind defined on the grid edges.
! Read all the values all dimensions but the time dimension.
! Only read the last time (if more than 1 present)
!

call nc_check(nf90_inq_varid(ncid, 'u', VarID), &
              'update_wind_components', 'inq_varid u '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'update_wind_components', 'inquire u '//trim(model_analysis_filename))

do i=1, numdims
   write(string1,*)'inquire u length for dimension ',i
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                 'update_wind_components', trim(string1)//' '//trim(model_analysis_filename))
   if (trim(dimname) == 'Time') then
      mystart(i)       = dimlen
      mycount(numdims) = 1
   else
      mystart(i)       = 1
      mycount(i)       = dimlen
   endif
enddo

call nc_check( nf90_get_var(ncid, VarID, u, start=mystart(1:numdims), count=mycount(1:numdims)), &
              'update_wind_components', 'get_var u '//trim(model_analysis_filename))

!
! Read the original uReconstructZonal
! Read all the values all dimensions but the time dimension.
! Only read the last time (if more than 1 present)
!

call nc_check(nf90_inq_varid(ncid, 'uReconstructZonal', VarID), &
              'update_wind_components', 'inq_varid uReconstructZonal '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'update_wind_components', 'inquire uReconstructZonal '//trim(model_analysis_filename))

do i=1, numdims
   write(string1,*)'inquire uReconstructZonal length for dimension ',i
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                 'update_wind_components', trim(string1)//' '//trim(model_analysis_filename))
   if (trim(dimname) == 'Time') then
      mystart(i)       = dimlen
      mycount(numdims) = 1
   else
      mystart(i)       = 1
      mycount(i)       = dimlen
   endif
enddo

call nc_check( nf90_get_var(ncid, VarID, ucell_incr, start=mystart(1:numdims), count=mycount(1:numdims)), &
              'update_wind_components', 'get_var uReconstructZonal '//trim(model_analysis_filename))

!
! Read uReconstructMeridional
!

call nc_check(nf90_inq_varid(ncid, 'uReconstructMeridional', VarID), &
              'update_wind_components', 'inq_varid uReconstructMeridional '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'update_wind_components', 'inquire uReconstructMeridional '//trim(model_analysis_filename))

do i=1, numdims
   write(string1,*)'inquire uReconstructMeridional length for dimension ',i
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen, name=dimname), &
                 'update_wind_components', trim(string1)//' '//trim(model_analysis_filename))
   if (trim(dimname) == 'Time') then
      mystart(i)       = dimlen
      mycount(numdims) = 1
   else
      mystart(i)       = 1
      mycount(i)       = dimlen
   endif
enddo

call nc_check( nf90_get_var(ncid, VarID, vcell_incr, start=mystart(1:numdims), count=mycount(1:numdims)), &
              'update_wind_components', 'get_var uReconstructMeridional '//trim(model_analysis_filename))

if ( debug > 7 ) then
   write(*,*)
   write(*,*)'update_wind_components: org u          range ',minval(u),          maxval(u)
   write(*,*)'update_wind_components: org zonal      range ',minval(ucell_incr), maxval(ucell_incr)
   write(*,*)'update_wind_components: org meridional range ',minval(vcell_incr), maxval(vcell_incr)
endif

! The state vector has updated zonal and meridional wind components.
! (Implicit in just being IN this routine)

zonal      = get_index_from_varname('uReconstructZonal')
meridional = get_index_from_varname('uReconstructMeridional')

if (zonal > 0) then
   allocate(data_2d_array(progvar(zonal)%dimlens(1),progvar(zonal)%dimlens(2)))
   call vector_to_prog_var(state_vector, zonal, data_2d_array )
   ucell_incr = data_2d_array - ucell_incr
   deallocate(data_2d_array)
else
   ucell_incr(:,:) = 0.0_r8
endif

if (meridional > 0) then
   allocate(data_2d_array(progvar(meridional)%dimlens(1),progvar(meridional)%dimlens(2)))
   call vector_to_prog_var(state_vector, meridional, data_2d_array)
   vcell_incr = data_2d_array - vcell_incr
   deallocate(data_2d_array)
else
   vcell_incr(:,:) = 0.0_r8
endif

if ( debug > 7 ) then
   write(*,*)
   write(*,*)'update_wind_components: u increment    range ',minval(ucell_incr), maxval(ucell_incr)
   write(*,*)'update_wind_components: v increment    range ',minval(vcell_incr), maxval(vcell_incr)
endif

! Now that we have the U,V increments (at the cell centers) and 
! the prior 'normal' wind component ('u' - at the cell edges),
! convert the increments to increments at the edges.

allocate(data_2d_array(nVertLevels, nEdges))

call uv_cell_to_edges(ucell_incr, vcell_incr, data_2d_array)

! Update normal velocity 
u(:,:) = u(:,:) + data_2d_array(:,:)

if ( debug > 7 ) then
   write(*,*)
   write(*,*)'update_wind_components: u after update:',minval(u), maxval(u)
endif

! Finally update the normal wind component field.

call put_u(u)

already_updated = .true. ! Change flag so we only do this routine once.

deallocate(data_2d_array, ucell_incr, vcell_incr, u)

! TJH FIXME can remove handle_winds, winds_present

end subroutine update_wind_components


!------------------------------------------------------------------

subroutine put_u(u)

! Put the newly updated 'u' field back into the netcdf file.
!
! The file name comes from module storage ... namelist.

! must first be allocated by calling code with the following sizes:
real(r8), intent(in) :: u(:,:)       ! u(nVertLevels, nEdges) 

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount, numu
integer :: ncid, VarID, numdims, nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: ntimes, i

! Read the netcdf file data

if ( .not. module_initialized ) call static_init_model


call nc_check(nf90_open(trim(model_analysis_filename), nf90_write, ncid), &
              'put_u', 'open '//trim(model_analysis_filename))

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
              'put_u', 'inquire '//trim(model_analysis_filename))

call nc_check(nf90_inquire_dimension(ncid, unlimitedDimID, len=ntimes), &
              'put_u', 'inquire time dimension length '//trim(model_analysis_filename))

call nc_check(nf90_inq_varid(ncid, 'u', VarID), &
              'put_u', 'inq_varid u '//trim(model_analysis_filename))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
              'put_u', 'inquire u '//trim(model_analysis_filename))

do i=1, numdims
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=numu(i)), &
                 'put_u', 'inquire U dimension length '//trim(model_analysis_filename))
enddo

! for all but the time dimension, read all the values.   
! for time read only the last one (if more than 1 present)
mystart = 1
mystart(numdims) = ntimes
mycount = numu
mycount(numdims) = 1

call nc_check(nf90_put_var(ncid, VarID, u, start=mystart, count=mycount), &
              'put_u', 'get_var u '//trim(model_analysis_filename))


call nc_check(nf90_close(ncid), 'put_u','close '//trim(model_analysis_filename) )


! A little sanity check

if ( debug > 7 ) then

   write(*,*)
   write(*,*)'u       range ',minval(u),     maxval(u)

endif

end subroutine put_u


!------------------------------------------------------------------

subroutine vector_to_1d_prog_var(x, ivar, data_1d_array)

! convert the values from a 1d array, starting at an offset,
! into a 1d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array

integer :: idim1,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim1 = 1, size(data_1d_array, 1)
   data_1d_array(idim1) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------

subroutine vector_to_2d_prog_var(x, ivar, data_2d_array)

! convert the values from a 1d array, starting at an offset,
! into a 2d array.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: idim1,idim2,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim2 = 1,size(data_2d_array, 2)
   do idim1 = 1,size(data_2d_array, 1)
      data_2d_array(idim1,idim2) = x(ii)
      ii = ii + 1
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var


!------------------------------------------------------------------

subroutine vector_to_3d_prog_var(x, ivar, data_3d_array)

! convert the values from a 1d array, starting at an offset,
! into a 3d array.

real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: idim1,idim2,idim3,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim3 = 1,size(data_3d_array, 3)
   do idim2 = 1,size(data_3d_array, 2)
      do idim1 = 1,size(data_3d_array, 1)
         data_3d_array(idim1,idim2,idim3) = x(ii)
         ii = ii + 1
      enddo
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var


!------------------------------------------------------------------

subroutine prog_var_1d_to_vector(data_1d_array, x, ivar)

! convert the values from a 1d array into a 1d array
! starting at an offset.

real(r8), dimension(:),   intent(in)    :: data_1d_array
real(r8), dimension(:),   intent(inout) :: x
integer,                  intent(in)    :: ivar

integer :: idim1,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim1 = 1, size(data_1d_array, 1)
   x(ii) = data_1d_array(idim1)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_1d_to_vector', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_1d_to_vector


!------------------------------------------------------------------

subroutine prog_var_2d_to_vector(data_2d_array, x, ivar)

! convert the values from a 2d array into a 1d array
! starting at an offset.

real(r8), dimension(:,:), intent(in)    :: data_2d_array
real(r8), dimension(:),   intent(inout) :: x
integer,                  intent(in)    :: ivar

integer :: idim1,idim2,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim2 = 1,size(data_2d_array, 2)
   do idim1 = 1,size(data_2d_array, 1)
      x(ii) = data_2d_array(idim1,idim2)
      ii = ii + 1
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_2d_to_vector', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_2d_to_vector


!------------------------------------------------------------------

subroutine prog_var_3d_to_vector(data_3d_array, x, ivar)

! convert the values from a 2d array into a 1d array
! starting at an offset.

real(r8), dimension(:,:,:), intent(in)    :: data_3d_array
real(r8), dimension(:),     intent(inout) :: x
integer,                    intent(in)    :: ivar

integer :: idim1,idim2,idim3,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar(ivar)%index1

do idim3 = 1,size(data_3d_array, 3)
   do idim2 = 1,size(data_3d_array, 2)
      do idim1 = 1,size(data_3d_array, 1)
         x(ii) = data_3d_array(idim1,idim2,idim3) 
         ii = ii + 1
      enddo
   enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' read wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'prog_var_3d_to_vector', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine prog_var_3d_to_vector


!------------------------------------------------------------------

subroutine verify_state_variables( state_variables, ncid, filename, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i, j, VarID
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname, dimname
character(len=NF90_MAX_NAME) :: dartstr
integer :: dimlen, numdims
logical :: failure

if ( .not. module_initialized ) call static_init_model

failure = .FALSE. ! perhaps all with go well

nrows = size(table,1)
ncols = size(table,2)

ngood = 0
MyLoop : do i = 1, nrows

   varname    = trim(state_variables(2*i -1))
   dartstr    = trim(state_variables(2*i   ))
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' ) then
      string1 = 'mpas_vars_nml:model state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in model analysis variable list

   write(string1,'(''variable '',a,'' in '',a)') trim(varname), trim(filename)
   write(string2,'(''there is no '',a)') trim(string1)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string2))

   ! Make sure variable is defined by (Time,nCells) or (Time,nCells,vertical)
   ! unable to support Edges or Vertices at this time.

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims), &
                 'verify_state_variables', 'inquire '//trim(string1))

   DimensionLoop : do j = 1,numdims

      write(string2,'(''inquire dimension'',i2,'' of '',a)') j,trim(string1)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(j), len=dimlen, name=dimname), &
                                          'verify_state_variables', trim(string2))
      select case ( trim(dimname) )
         case ('Time')
            ! supported - do nothing
         case ('nCells')
            ! supported - do nothing
         case ('nEdges')
            ! supported - do nothing
         case ('nVertLevels')
            ! supported - do nothing
         case ('nVertLevelsP1')
            ! supported - do nothing
         case ('nSoilLevels')
            ! supported - do nothing
         case default
            write(string2,'(''unsupported dimension '',a,'' in '',a)') trim(dimname),trim(string1)
            call error_handler(E_MSG,'verify_state_variables',string2,source,revision,revdate)
            failure = .TRUE.
      end select

   enddo DimensionLoop

   if (failure) then
       string2 = 'unsupported dimension(s) are fatal'
       call error_handler(E_ERR,'verify_state_variables',string2,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector 

   if ( debug > 0 .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

! TJH FIXME need to add check so they cannot have both normal winds and reconstructed winds in
! DART state vector.

end subroutine verify_state_variables


!------------------------------------------------------------------

subroutine dump_progvar(ivar, x)

 integer,  intent(in)           :: ivar
 real(r8), intent(in), optional :: x(:)

! if present, x is a state vector.  dump the data min/max for this var.

!%! type progvartype
!%!    private
!%!    character(len=NF90_MAX_NAME) :: varname
!%!    character(len=NF90_MAX_NAME) :: long_name
!%!    character(len=NF90_MAX_NAME) :: units
!%!    integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
!%!    integer :: xtype         ! netCDF variable type (NF90_double, etc.) 
!%!    integer :: numdims       ! number of dims - excluding TIME
!%!    integer :: numvertical   ! number of vertical levels in variable
!%!    integer :: numcells      ! number of horizontal locations (typically cell centers)
!%!    integer :: numedges
!%!    logical :: ZonHalf       ! vertical coordinate has dimension nVertLevels
!%!    integer :: varsize       ! prod(dimlens(1:numdims))
!%!    integer :: index1        ! location in dart state vector of first occurrence
!%!    integer :: indexN        ! location in dart state vector of last  occurrence
!%!    integer :: dart_kind
!%!    character(len=paramname_length) :: kind_string
!%!    logical  :: clamping     ! does variable need to be range-restricted before 
!%!    real(r8) :: range(2)     ! being stuffed back into MPAS analysis file.
!%! end type progvartype

integer :: i

! take care of parallel runs where we only want a single copy of
! the output.
if (.not. do_output()) return

write(logfileunit,*)
write(     *     ,*)
write(logfileunit,*) 'variable number ',ivar,' is ',trim(progvar(ivar)%varname)
write(     *     ,*) 'variable number ',ivar,' is ',trim(progvar(ivar)%varname)
write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
write(logfileunit,*) '  numvertical ',progvar(ivar)%numvertical
write(     *     ,*) '  numvertical ',progvar(ivar)%numvertical
write(logfileunit,*) '  numcells    ',progvar(ivar)%numcells
write(     *     ,*) '  numcells    ',progvar(ivar)%numcells
write(logfileunit,*) '  numedges    ',progvar(ivar)%numedges
write(     *     ,*) '  numedges    ',progvar(ivar)%numedges
write(logfileunit,*) '  ZonHalf     ',progvar(ivar)%ZonHalf
write(     *     ,*) '  ZonHalf     ',progvar(ivar)%ZonHalf
write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
write(logfileunit,*) '  index1      ',progvar(ivar)%index1
write(     *     ,*) '  index1      ',progvar(ivar)%index1
write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string
write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
write(logfileunit,*) '  clamping    ',progvar(ivar)%clamping
write(     *     ,*) '  clamping    ',progvar(ivar)%clamping
write(logfileunit,*) '  clmp range  ',progvar(ivar)%range
write(     *     ,*) '  clmp range  ',progvar(ivar)%range
do i = 1,progvar(ivar)%numdims
   write(logfileunit,*) '  dimension/length/name ',i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimname(i))
   write(     *     ,*) '  dimension/length/name ',i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimname(i))
enddo

if (present(x)) then
   write(logfileunit, * )                 'min/max = ', &
              minval(x(progvar(ivar)%index1:progvar(ivar)%indexN)), &
              maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))
   write(    *      , * )                 'min/max = ', &
              minval(x(progvar(ivar)%index1:progvar(ivar)%indexN)), &
              maxval(x(progvar(ivar)%index1:progvar(ivar)%indexN))
endif

end subroutine dump_progvar


!------------------------------------------------------------------

function FindTimeDimension(ncid) result(timedimid)

! Find the Time Dimension ID in a netCDF file.
! If there is none - (spelled the obvious way) - the routine
! returns a negative number. You don't HAVE to have a TIME dimension.

integer                      :: timedimid
integer,          intent(in) :: ncid

integer :: nc_rc

TimeDimID = -1 ! same as the netCDF library routines. 
nc_rc = nf90_inq_dimid(ncid,'Time',dimid=TimeDimID)

end function FindTimeDimension


!------------------------------------------------------------------

subroutine winds_present(zonal,meridional,both)

integer, intent(out) :: zonal, meridional
logical, intent(out) :: both

! if neither of uReconstructZonal or uReconstructMeridional are in the
!   state vector, set both to .false. and we're done.
! if both are there, return the ivar indices for each
! if only one is there, it's an error.

zonal      = get_index_from_varname('uReconstructZonal')
meridional = get_index_from_varname('uReconstructMeridional')

if (zonal > 0 .and. meridional > 0) then
  both = .true.
  return
else if (zonal < 0 .and. meridional < 0) then
  both = .false. 
  return
endif

! only one present - error.
write(string1,*) 'both components for U winds must be in state vector'
call error_handler(E_ERR,'winds_present',string1,source,revision,revdate)

end subroutine winds_present


!------------------------------------------------------------------

subroutine handle_winds(u_prior , zonal_incr, meridional_incr)

 real(r8), intent(in) :: u_prior(:,:), zonal_incr(:,:), meridional_incr(:,:)
 
!  the current plan for winds is:
!  read the reconstructed zonal and meridional winds at cell centers
!  do the assimilation of U,V winds and update the centers
!  at write time, compute the increments at the centers, 
!  map the increments back to the edges rotating
!  to be normal and parallel to the edge directions.
!  We write out the updated U,V winds at cell centers to the analysis 
!  file although they are going to be reinitialized when the model
!  gets started (since they are diagnostic variables in the model).

! the 'cellsOnEdge' array has the ids of the two neighboring cell numbers
! the lonEdge/latEdge arrays have the ?midpoints of the edges.
! the ends of each edge are the lonVertex/latVertex arrays
! the location of the winds are at the lonCell/latCell arrays (cell centers)

real(r8), allocatable :: u(:,:), du(:,:)

allocate(u(nVertLevels, nEdges))
allocate(du(nVertLevels, nEdges))

if ( debug > 7 ) then
write(*,*)'u wind increment:',minval(zonal_incr), maxval(zonal_incr)
write(*,*)'v wind increment:',minval(meridional_incr), maxval(meridional_incr)
write(*,*)'u_prior: ',minval(u_prior), maxval(u_prior)
endif

call uv_cell_to_edges(zonal_incr, meridional_incr, du)

if ( debug > 7 ) then
write(*,*)
write(*,*)'du after uv_cell_to_edges:',minval(du), maxval(du)
endif

! Update normal velocity 
u(:,:) = u_prior(:,:) + du(:,:)

if ( debug > 7 ) then
write(*,*)
write(*,*)'u after update:',minval(u), maxval(u)
endif

call put_u(u)

deallocate(u, du) 

end subroutine handle_winds


!------------------------------------------------------------

subroutine set_variable_clamping(ivar)

! The model may behave poorly if some quantities are outside
! a physically realizable range.
!
! FIXME : add more DART types
! FIXME2 : need to be able to set just min or just max
!  and leave the other MISSING_R8

integer, intent(in) :: ivar

select case (trim(progvar(ivar)%kind_string))
   case ('KIND_VAPOR_MIXING_RATIO')
      progvar(ivar)%clamping = .true.
      progvar(ivar)%range    = (/ 0.0_r8, 1.0_r8 /)
   case default
      progvar(ivar)%clamping = .false.
      progvar(ivar)%range    = MISSING_R8
end select

end subroutine set_variable_clamping


!------------------------------------------------------------

subroutine define_var_dims(ncid,ivar, memberdimid, unlimiteddimid, ndims, dimids) 

! set the dimids array needed to augment the natural shape of the variable
! with the two additional dimids needed by the DART diagnostic output.
integer,               intent(in)  :: ncid
integer,               intent(in)  :: ivar
integer,               intent(in)  :: memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i,mydimid

ndims  = 0
dimids = 0

do i = 1,progvar(ivar)%numdims

   ! Each of these dimension names (originally from the MPAS analysis file)
   ! must exist in the DART diagnostic netcdf files. 

   call nc_check(nf90_inq_dimid(ncid, trim(progvar(ivar)%dimname(i)), mydimid), &
              'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimname(i)))

   ndims = ndims + 1

   dimids(ndims) = mydimid

enddo

ndims         = ndims + 1
dimids(ndims) = memberdimid
ndims         = ndims + 1
dimids(ndims) = unlimiteddimid

end subroutine define_var_dims


!------------------------------------------------------------

subroutine get_index_range_string(string,index1,indexN)

! Determine where a particular DART kind (string) exists in the 
! DART state vector.

character(len=*),  intent(in)  :: string
integer,           intent(out) :: index1
integer, optional, intent(out) :: indexN

integer :: i

index1 = 0
if (present(indexN)) indexN = 0

FieldLoop : do i=1,nfields
   if (progvar(i)%kind_string /= trim(string)) cycle FieldLoop
   index1 = progvar(i)%index1
   if (present(indexN)) indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

if (index1 == 0) then
   write(string1,*) 'Problem, cannot find indices for '//trim(string)
   call error_handler(E_ERR,'get_index_range_string',string1,source,revision,revdate)
endif
end subroutine get_index_range_string


!------------------------------------------------------------------

subroutine get_index_range_int(dartkind,index1,indexN)

! Determine where a particular DART kind (integer) exists in the 
! DART state vector.

integer,           intent(in)  :: dartkind
integer,           intent(out) :: index1
integer, optional, intent(out) :: indexN

integer :: i
character(len=paramname_length) :: string

index1 = 0
if (present(indexN)) indexN = 0

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   index1 = progvar(i)%index1
   if (present(indexN)) indexN = progvar(i)%indexN
   exit FieldLoop
enddo FieldLoop

string = get_raw_obs_kind_name(dartkind)

if (index1 == 0) then
   write(string1,*) 'Problem, cannot find indices for kind ',dartkind,trim(string)
   call error_handler(E_ERR,'get_index_range_int',string1,source,revision,revdate)
endif

end subroutine get_index_range_int


!------------------------------------------------------------------

function get_progvar_index_from_kind(dartkind)

! Determine what index a particular DART kind (integer) is in the
! progvar array.
integer :: get_progvar_index_from_kind
integer, intent(in) :: dartkind

integer :: i

FieldLoop : do i=1,nfields
   if (progvar(i)%dart_kind /= dartkind) cycle FieldLoop
   get_progvar_index_from_kind = i 
   return
enddo FieldLoop

get_progvar_index_from_kind = -1

end function get_progvar_index_from_kind


!------------------------------------------------------------------

function get_index_from_varname(varname)

! Determine what index corresponds to the given varname
! if name not in state vector, return -1 -- not an error.

integer :: get_index_from_varname
character(len=*), intent(in) :: varname

integer :: i

FieldLoop : do i=1,nfields
   if (trim(progvar(i)%varname) == trim(varname)) then
      get_index_from_varname = i
      return
   endif 
enddo FieldLoop

get_index_from_varname = -1
return

end function get_index_from_varname


!==================================================================
! The following (private) interfaces are used for triangle interpolation
!==================================================================


!------------------------------------------------------------------

subroutine find_pressure_bounds(x, p, tri_indices, weights, nbounds, &
   pt_base_offset, density_base_offset, qv_base_offset, lower, upper, fract, &
   ltemp, utemp, ier)

! Finds vertical interpolation indices and fraction for a quantity with 
! pressure vertical coordinate. Loops through the height levels and
! computes the corresponding pressure at the horizontal point.  nbounds is
! the number of vertical levels in the potential temperature, density,
! and water vapor grids.

real(r8),  intent(in)  :: x(:)
real(r8),  intent(in)  :: p
integer,   intent(in)  :: tri_indices(3)
real(r8),  intent(in)  :: weights(3)
integer,   intent(in)  :: nbounds
integer,   intent(in)  :: pt_base_offset, density_base_offset, qv_base_offset
integer,   intent(out) :: lower, upper
real(r8),  intent(out) :: fract
real(r8),  intent(out) :: ltemp, utemp
integer,   intent(out) :: ier

integer  :: i, gip_err
real(r8) :: pressure(nbounds)

! Default error return is 0
ier = 0

! Find the lowest pressure
call get_interp_pressure(x, pt_base_offset, density_base_offset, qv_base_offset, &
   tri_indices, weights, 1, nbounds, pressure(1), gip_err, ltemp)
if(gip_err /= 0) then
   ier = gip_err
   return
endif

! Get the highest pressure level
call get_interp_pressure(x, pt_base_offset, density_base_offset, qv_base_offset, &
   tri_indices, weights, nbounds, nbounds, pressure(nbounds), gip_err)
if(gip_err /= 0) then
   ier = gip_err
   return
endif

! Check for out of the column range
if(p > pressure(1) .or. p < pressure(nbounds)) then
   ier = 2
   return
endif

! Loop through the rest of the column from the bottom up
do i = 2, nbounds
   call get_interp_pressure(x, pt_base_offset, density_base_offset, qv_base_offset, &
      tri_indices, weights, i, nbounds, pressure(i), gip_err, utemp)
   if(gip_err /= 0) then
      ier = gip_err
      return
   endif

   ! Is pressure between i-1 and i level?
   if(p > pressure(i)) then
      lower = i - 1
      upper = i
      ! FIXME: should this be interpolated in log(p)??  yes.
      fract = (p - pressure(i-1)) / (pressure(i) - pressure(i-1))
      !fract = exp(log(p) - log(pressure(i-1))) / (log(pressure(i)) - log(pressure(i-1)))
      return
   endif
   
   ltemp = utemp
end do

! Shouldn't ever fall off end of loop
ier = 3

end subroutine find_pressure_bounds

!------------------------------------------------------------------

subroutine get_interp_pressure(x, pt_offset, density_offset, qv_offset, &
   tri_indices, weights, lev, nlevs, pressure, ier, temperature)

! Finds the value of pressure at a given point at model level lev

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: pt_offset, density_offset, qv_offset
integer,  intent(in)  :: tri_indices(3)
real(r8), intent(in)  :: weights(3)
integer,  intent(in)  :: lev, nlevs
real(r8), intent(out) :: pressure
integer,  intent(out) :: ier
real(r8), intent(out), optional :: temperature

integer  :: i, offset
real(r8) :: pt(3), density(3), qv(3), pt_int, density_int, qv_int, tk


! Get the values of potential temperature, density, and vapor at each corner
do i = 1, 3
   offset = (tri_indices(i) - 1) * nlevs + lev - 1
   pt(i) =      x(pt_offset + offset)
   density(i) = x(density_offset + offset)
   qv(i) =      x(qv_offset + offset)
   ! Error if any of the values are missing; probably will be all or nothing
   if(pt(i) == MISSING_R8 .or. density(i) == MISSING_R8 .or. qv(i) == MISSING_R8) then
      ier = 2
      return
   endif
end do

! Interpolate three state values in horizontal
pt_int =      sum(weights * pt)
density_int = sum(weights * density)
qv_int =      sum(weights * qv)

! Get pressure at the interpolated point
call compute_full_pressure(pt_int, density_int, qv_int, pressure, tk)

! if the caller asked for sensible temperature, return it
if (present(temperature)) temperature = tk

! Default is no error
ier = 0

end subroutine get_interp_pressure


!------------------------------------------------------------------

subroutine vert_interp(x, base_offset, cellid, nlevs, lower, fract, val, ier)

! Interpolates in vertical in column indexed by tri_index for a field
! with base_offset.  Vertical index is varying fastest here. Returns ier=0
! unless missing value is encounterd. 

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: base_offset
integer,  intent(in)  :: cellid
integer,  intent(in)  :: nlevs
integer,  intent(in)  :: lower
real(r8), intent(in)  :: fract
real(r8), intent(out) :: val
integer,  intent(out) :: ier

integer  :: offset
real(r8) :: lx, ux

! Default return is good
ier = 0

! Get the value at the lower and upper points
offset = base_offset + (cellid - 1) * nlevs + lower - 1
lx = x(offset)
ux = x(offset + 1)

! Check for missing value
if(lx == MISSING_R8 .or. ux == MISSING_R8) then
   ier = 2
   return
endif 

! Interpolate
val = (1.0_r8 - fract)*lx + fract*ux

end subroutine vert_interp


!------------------------------------------------------------------

subroutine find_height_bounds(height, nbounds, bounds, lower, upper, fract, ier)

! Finds position of a given height in an array of height grid points and returns
! the index of the lower and upper bounds and the fractional offset.  ier returns 0 
! unless there is an error. Could be replaced with a more efficient search if there 
! are many vertical levels.

real(r8), intent(in)  :: height
integer,  intent(in)  :: nbounds
real(r8), intent(in)  :: bounds(nbounds)
integer,  intent(out) :: lower, upper
real(r8), intent(out) :: fract
integer,  intent(out) :: ier

! Assume that the spacing on altitudes is arbitrary and do the simple thing
! which is a linear search. Probably not worth any fancier searching unless
! models get to be huge.

integer :: i

if(height < bounds(1) .or. height > bounds(nbounds)) then
   !JLA 
   !write(*, *) 'fail in find_height_bounds ', height, bounds(1), bounds(nbounds)
   ier = 2
   return
endif

do i = 2, nbounds
   if(height <= bounds(i)) then
      lower = i - 1
      upper = i
      fract = (height - bounds(lower)) / (bounds(upper) - bounds(lower))
      ier = 0
      return
   endif
end do

! Shouldn't ever fall off end of loop
ier = 3

end subroutine find_height_bounds

!------------------------------------------------------------------

subroutine find_vert_level(x, loc, oncenters, lower, upper, fract, ier)

! given a location and var types, return the two level numbers that 
! enclose the given vertical value plus the fraction between them.   

! FIXME:  this handles data at cell centers, at edges, but not
! data on faces.

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: loc
logical,             intent(in)  :: oncenters
integer,             intent(out) :: lower, upper
real(r8),            intent(out) :: fract
integer,             intent(out) :: ier

real(r8) :: lat, lon, vert, tmp(3)
integer  :: verttype, cellid, edgeid
integer  :: pt_base_offset, density_base_offset, qv_base_offset

! the plan is to take in: whether this var is on cell centers or edges,
! and the location so we can extract the vert value and which vert flag.
! compute and return the lower and upper index level numbers that
! enclose this vert, along with the fract between them.  ier is set
! in case of error (e.g. outside the grid, on dry land, etc).

! kinds we have to handle:
!vert_is_undef,    VERTISUNDEF
!vert_is_surface,  VERTISSURFACE
!vert_is_level,    VERTISLEVEL
!vert_is_pressure, VERTISPRESSURE
!vert_is_height,   VERTISHEIGHT

! unpack the location into local vars
tmp = get_location(loc)
lon  = tmp(1)
lat  = tmp(2)
vert = tmp(3)
verttype = nint(query_location(loc))

! these first 3 types need no cell/edge location information.

! no defined vertical location (e.g. vertically integrated vals)
if (vert_is_undef(loc)) then
   ier = 12
   return
endif

! vertical is defined to be on the surface (level 1 here)
if(vert_is_surface(loc)) then
   lower = 1
   upper = 2
   fract = 0.0_r8
   ier = 0
   return
endif

! model level numbers (supports fractional levels)
if(vert_is_level(loc)) then
   ! FIXME: if this is W, the top is nVertLevels+1
   if (vert > nVertLevels) then
      ier = 12
      return
   endif
   if (vert == nVertLevels) then
      lower = nint(vert) - 1   ! round down
      upper = nint(vert)
      fract = 1.0_r8
      ier = 0
      return
   endif
   lower = aint(vert)   ! round down
   upper = lower+1
   fract = vert - lower
!print *, '1 lower, upper = ', lower, upper
   ier = 0
   return
endif

! ok, now we need to know where we are in the grid for heights or pressures
! as the vertical coordinate.


! find the cell/edge that contains this point.
cellid = find_closest_cell_center(lat, lon)
!print *, 'found closest cell center to ', lon, lat, ' which is ', cellid
if (.not. inside_cell(cellid, lat, lon)) then
   ier = 13   
   return
endif
! FIXME: here seems like where we should handle data on cell face centers
if (.not. oncenters) then
   edgeid = find_closest_edge(cellid, lat, lon)
   !print *, 'found closest edge = ', edgeid
   !print *, 'size of edge height array = ', shape(zGridEdge)
   !print *, '(closest cellid id: ', cellid, ')'
   !print *, 'size of center height array = ', shape(zGridCenter)
endif

! Vertical interpolation for pressure coordinates
if(vert_is_pressure(loc) ) then 
   ! Need to get base offsets for the potential temperature, density, and water 
   ! vapor mixing fields in the state vector
   call get_index_range(KIND_POTENTIAL_TEMPERATURE, pt_base_offset)
   call get_index_range(KIND_DENSITY, density_base_offset)
   call get_index_range(KIND_VAPOR_MIXING_RATIO, qv_base_offset)
!print *, '2bases: t/rho/v = ', pt_base_offset, density_base_offset, qv_base_offset
   call find_pressure_bounds2(x, vert, cellid, nVertLevels, &
         pt_base_offset, density_base_offset, qv_base_offset,  &
         lower, upper, fract, ier)

!print *, '2 lower, upper, ier = ', lower, upper, ier
   return
endif

! grid is in height, so this needs to know which cell to index into
! for the column of heights and call the bounds routine.
if(vert_is_height(loc)) then
   ! For height, can do simple vertical search for interpolation for now
   ! Get the lower and upper bounds and fraction for each column
   if (oncenters) then
      call find_height_bounds(vert, nVertLevels, zGridCenter(:, cellid), &
                              lower, upper, fract, ier)
   else
      call find_height_bounds(vert, nVertLevels, zGridEdge(:, edgeid), &
                              lower, upper, fract, ier)
   endif
!print *, '3 lower, upper = ', lower, upper
   return
endif

! Shouldn't ever fall out of the 'if' before returning
ier = 3

end subroutine find_vert_level

!------------------------------------------------------------------

subroutine find_pressure_bounds2(x, p, cellid, nbounds, &
   pt_base_offset, density_base_offset, qv_base_offset, &
   lower, upper, fract, ier)

! Finds vertical interpolation indices and fraction for a quantity with 
! pressure vertical coordinate. Loops through the height levels and
! computes the corresponding pressure at the horizontal point.  nbounds is
! the number of vertical levels in the potential temperature, density,
! and water vapor grids.

real(r8),  intent(in)  :: x(:)
real(r8),  intent(in)  :: p
integer,   intent(in)  :: cellid
integer,   intent(in)  :: nbounds
integer,   intent(in)  :: pt_base_offset, density_base_offset, qv_base_offset
integer,   intent(out) :: lower, upper
real(r8),  intent(out) :: fract
integer,   intent(out) :: ier

integer  :: i, gip_err
real(r8) :: pressure(nbounds)

! Default error return is 0
ier = 0

! Find the lowest pressure
call get_interp_pressure2(x, pt_base_offset, density_base_offset, qv_base_offset, &
   cellid, 1, nbounds, pressure(1), gip_err)
!print *, 'find p bounds2, pr(1) = ', pressure(1), gip_err
if(gip_err /= 0) then
   ier = gip_err
   return
endif

! Get the highest pressure level
call get_interp_pressure2(x, pt_base_offset, density_base_offset, qv_base_offset, &
   cellid, nbounds, nbounds, pressure(nbounds), gip_err)
!print *, 'find p bounds2, pr(n) = ', pressure(nbounds), gip_err
if(gip_err /= 0) then
   ier = gip_err
   return
endif

! Check for out of the column range
if(p > pressure(1) .or. p < pressure(nbounds)) then
!print *, 'find p bounds2, p, pr(1), pr(n) = ', p, pressure(1), pressure(nbounds)
   ier = 2
   return
endif

! Loop through the rest of the column from the bottom up
do i = 2, nbounds
   call get_interp_pressure2(x, pt_base_offset, density_base_offset, qv_base_offset, &
      cellid, i, nbounds, pressure(i), gip_err)
!print *, 'find p bounds i, pr(i) = ', i, pressure(i), gip_err
   if(gip_err /= 0) then
      ier = gip_err
      return
   endif

   ! Is pressure between i-1 and i level?
   if(p > pressure(i)) then
      lower = i - 1
      upper = i
      ! FIXME: should this be interpolated in log(p)??  yes.
      if (pressure(i) == pressure(i-1)) then
         fract = 0.0_r8
      else
         fract = (p - pressure(i-1)) / (pressure(i) - pressure(i-1))
         !fract = exp(log(p) - log(pressure(i-1))) / (log(pressure(i)) - log(pressure(i-1)))
      endif
!print *, "looping pr: i, p, pr(i), lower, upper, fract = ", i, p, pressure(i), lower, upper, fract
      return
   endif

end do

!print *, 'fell off end, find pressure bounds 2'
! Shouldn't ever fall off end of loop
ier = 3

end subroutine find_pressure_bounds2

!------------------------------------------------------------------

subroutine get_interp_pressure2(x, pt_offset, density_offset, qv_offset, &
   cellid, lev, nlevs, pressure, ier)

! Finds the value of pressure at a given point at model level lev

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: pt_offset, density_offset, qv_offset
integer,  intent(in)  :: cellid
integer,  intent(in)  :: lev, nlevs
real(r8), intent(out) :: pressure
integer,  intent(out) :: ier

integer  :: i, offset
real(r8) :: pt, density, qv, tk


! Get the values of potential temperature, density, and vapor 
offset = (cellid - 1) * nlevs + lev - 1
pt =      x(pt_offset + offset)
density = x(density_offset + offset)
qv =      x(qv_offset + offset)
! Error if any of the values are missing; probably will be all or nothing
if(pt == MISSING_R8 .or. density == MISSING_R8 .or. qv == MISSING_R8) then
   ier = 2
   return
endif
!print *, 'offset: base pt, dens, qv = ', offset, pt_offset, density_offset, qv_offset
!print *, 'vals: pt, dens, qv = ', pt, density, qv

! Get pressure at the cell center
call compute_full_pressure(pt, density, qv, pressure, tk)

! Default is no error
ier = 0

end subroutine get_interp_pressure2


!------------------------------------------------------------------

subroutine get_cell_indices(lon, lat, indices, weights, istatus)

! Finds the indices of the three cell centers that form the vertices of
! the triangle that contains the given (lon, lat) point. Also returns
! the barycentric weights for the three corners for this point. Returns istatus
! 0 if successful, otherwise nonzero.

real(r8),            intent(in) :: lon, lat
integer,            intent(out) :: indices(3)
real(r8),           intent(out) :: weights(3)
integer,            intent(out) :: istatus

! Local storage
integer  :: num_inds, start_ind
integer  :: x_ind, y_ind

! Succesful return has istatus of 0
istatus = 0

! Figure out which of the regular grid boxes this is in
call get_reg_box_indices(lon, lat, x_ind, y_ind)
num_inds =  triangle_num  (x_ind, y_ind)
start_ind = triangle_start(x_ind, y_ind)

! If there are no triangles overlapping, can't do interpolation
if(num_inds == 0) then
   istatus = 1
   return
endif

! Search the list of triangles to see if (lon, lat) is in one
call get_triangle(lon, lat, num_inds, start_ind, indices, weights, istatus)
!print *, 'got triangle.  indices = ', indices
   

if(istatus /= 0) istatus = 2

end subroutine get_cell_indices


!------------------------------------------------------------

subroutine get_triangle(lon, lat, num_inds, start_ind, indices, weights, istatus)

! Given a latitude longitude point, and the starting address in the list for the 
! triangles that might contain this lon lat in the list, finds the
! triangle that contains the point and returns the indices and barycentric
! weights. Returns istatus of 0 if a value is found and istatus of 1 if
! there is no value found.

real(r8), intent(in)  :: lon, lat
integer,  intent(in)  :: num_inds, start_ind
integer,  intent(out) :: indices(3)
real(r8), intent(out) :: weights(3)
integer,  intent(out) :: istatus

integer :: i, j, ind
real(r8) :: clons(3), clats(3)
real(r8) :: lonmax, lonmin, lonfix

! Assume successful until proven otherwise
istatus = 0

! Loop through the candidate triangles
do i = start_ind, start_ind + num_inds - 1
   ! Get the index of the triangle
   ind = triangle_list(i)
   ! Get corner lons and lats
   call get_triangle_corners(ind, clons, clats)

   ! Need to deal with longitude wraparound before doing triangle computations.
   ! Begin by finding the range of the longitudes (including the target point).
   lonmax = max(lon, maxval(clons))
   lonmin = min(lon, minval(clons))
   ! lonfix is used to store the target longitude
   lonfix = lon

   ! If the range is more than 180.0, assume that points wrapped around 0 longitude
   ! Move the points to a 180 to 540 degree representation
   if((lonmax - lonmin) > 180.0_r8) then
      do j = 1, 3
         if(clons(j) < 180.0_r8) clons(j) = clons(j) + 360.0_r8
         if(lonfix < 180.0_r8) lonfix = lonfix + 360.0_r8
      end do
   endif

   ! Get the barycentric weights 
   call get_barycentric_weights(lonfix, lat, clons, clats, weights)
   ! Is point in this triangle? Yes, if weights are in range [0 1]
   ! If so, return the indices of corners. 
   if(maxval(weights) <= 1.0_r8 .and. minval(weights) >= 0.0_r8) then
      ! Get the indices for this triangle
      indices(:) = cellsOnVertex(:, ind)
      return
   endif
end do

! Falling off the end means failure for now (could weakly extrapolate)
weights = -1.0_r8
indices = -1
istatus = 1

end subroutine get_triangle


!------------------------------------------------------------

subroutine triangle_interp(x, base_offset, tri_indices, weights, &
   level, nlevels, interp_val, ier) 

! Given state, offset for start of horizontal slice, the indices of the
! triangle vertices in that slice, and the barycentric weights, computes
! the interpolated value. Returns ier=0 unless a missing value is found.

real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: base_offset
integer,  intent(in)  :: tri_indices(3)
real(r8), intent(in)  :: weights(3)
integer,  intent(in)  :: level, nlevels
real(r8), intent(out) :: interp_val
integer,  intent(out) :: ier

integer  :: i, offset
real(r8) :: corner_val(3)

! Find the three corner values
do i = 1, 3
   offset = base_offset + (tri_indices(i) -1) * nlevels + level - 1
   corner_val(i) = x(offset)
   if(corner_val(i) == MISSING_R8) then
      ier = 1
      interp_val = MISSING_R8
      return
   endif
end do

interp_val = sum(weights * corner_val)
ier = 0

end subroutine triangle_interp


!------------------------------------------------------------

subroutine get_3d_weights(p, v1, v2, v3, weights)

! Given a point p (x,y,z) inside a triangle, and the (x,y,z)
! coordinates of the triangle corner points (v1, v2, v3), 
! find the weights for a barycentric interpolation.  this
! computation only needs two of the three coordinates, so figure
! out which quadrant of the sphere the triangle is in and pick
! the 2 axes which are the least planar:
!  (x,y) near the poles,
!  (y,z) near 0 and 180 longitudes near the equator,
!  (x,z) near 90 and 270 longitude near the equator.

real(r8), intent(in)  :: p(3)
real(r8), intent(in)  :: v1(3), v2(3), v3(3)
real(r8), intent(out) :: weights(3)

real(r8) :: lat, lon
real(r8) :: cxs(3), cys(3)

! FIXME: this costs - if we already have the lat/lon, just
! pass them down?
call xyz_to_latlon(p(1), p(2), p(3), lat, lon)

! above or below 45 in latitude, where -90 < lat < 90:
if (lat >= 45.0_r8 .or. lat <= -45.0_r8) then
   cxs(1) = v1(1)
   cxs(2) = v2(1)
   cxs(3) = v3(1)
   cys(1) = v1(2)
   cys(2) = v2(2)
   cys(3) = v3(2)
   call get_barycentric_weights(p(1), p(2), cxs, cys, weights)
   return
endif

! nearest 0 or 180 in longitude, where 0 < lon < 360:
if ( lon <= 45.0_r8 .or. lon >= 315.0_r8 .or. &
    (lon >= 135.0_r8 .and. lon <= 225.0_r8)) then
   cxs(1) = v1(2)
   cxs(2) = v2(2)
   cxs(3) = v3(2)
   cys(1) = v1(3)
   cys(2) = v2(3)
   cys(3) = v3(3)
   call get_barycentric_weights(p(2), p(3), cxs, cys, weights)
   return
endif

! last option, nearest 90 or 270 in lon:
cxs(1) = v1(1)
cxs(2) = v2(1)
cxs(3) = v3(1)
cys(1) = v1(3)
cys(2) = v2(3)
cys(3) = v3(3)
call get_barycentric_weights(p(1), p(3), cxs, cys, weights)

end subroutine get_3d_weights


!------------------------------------------------------------

subroutine get_barycentric_weights(x, y, cxs, cys, weights)

! Computes the barycentric weights for a 2d interpolation point 
! (x,y) in a 2d triangle with the given (cxs,cys) corners.

real(r8), intent(in)  :: x, y, cxs(3), cys(3)
real(r8), intent(out) :: weights(3)

real(r8) :: denom

! Get denominator
denom = (cys(2) - cys(3)) * (cxs(1) - cxs(3)) + &
   (cxs(3) - cxs(2)) * (cys(1) - cys(3))

weights(1) = ((cys(2) - cys(3)) * (x - cxs(3)) + &
   (cxs(3) - cxs(2)) * (y - cys(3))) / denom

weights(2) = ((cys(3) - cys(1)) * (x - cxs(3)) + &
   (cxs(1) - cxs(3)) * (y - cys(3))) / denom

weights(3) = 1.0_r8 - weights(1) - weights(2)

end subroutine get_barycentric_weights


!------------------------------------------------------------

subroutine get_reg_box_indices(lon, lat, x_ind, y_ind)

! Given a longitude and latitude in degrees returns the index of the regular
! lon-lat box that contains the point.

real(r8), intent(in)  :: lon, lat
integer,  intent(out) :: x_ind, y_ind

call get_reg_lon_box(lon, x_ind)
call get_reg_lat_box(lat, y_ind)

end subroutine get_reg_box_indices


!------------------------------------------------------------

subroutine get_reg_lon_box(lon, x_ind)

! Determine which regular longitude box a longitude is in.

real(r8), intent(in)  :: lon
integer,  intent(out) :: x_ind

x_ind = int(num_reg_x * lon / 360.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lon == 360.0_r8) x_ind = num_reg_x

end subroutine get_reg_lon_box


!------------------------------------------------------------

subroutine get_reg_lat_box(lat, y_ind)

! Determine which regular latitude box a latitude is in.

real(r8), intent(in)  :: lat
integer,  intent(out) :: y_ind

y_ind = int(num_reg_y * (lat + 90.0_r8) / 180.0_r8) + 1

! Watch out for exactly at top; assume all lats and lons in legal range
if(lat == 90.0_r8)  y_ind = num_reg_y

end subroutine get_reg_lat_box


!------------------------------------------------------------

subroutine init_interp()

! Initializes data structures needed for MPAS interpolation.
! This should be called at static_init_model time to avoid 
! having all this temporary storage in the middle of a run.

! Build the data structure for interpolation for a triangle grid.
! Need a temporary data structure to build this.
! This array keeps a list of the indices of cell center triangles
! that potentially overlap each regular boxes.

! Current version assumes that triangles are on lat/lon grid. This
! leads to interpolation errors that increase for near the poles 
! and for coarse grids. At some point need to at least treat 
! triangles as being locally inscribed in the sphere. Actual
! exact spherical geometry is almost certainly not needed to
! forward operator interpolation.

integer :: reg_list(num_reg_x, num_reg_y, max_reg_list_num)

real(r8) :: c_lons(3), c_lats(3)
integer  :: i, j, k, ier
integer  :: reg_lon_ind(2), reg_lat_ind(2), total, ind


! Loop through each of the triangles
do i = 1, nVertices

   ! Set up array of lons and lats for the corners of this triangle
   call get_triangle_corners(i, c_lons, c_lats)

   ! Get list of regular boxes that cover this triangle.
   call reg_box_overlap(c_lons, c_lats, reg_lon_ind, reg_lat_ind, ier)

   ! Update the temporary data structures of triangles that overlap regular quads
   ! If this triangle had pole, don't add it
   if(ier == 0) &
      call update_reg_list(triangle_num, reg_list, reg_lon_ind, reg_lat_ind, i)

enddo

!if (do_output()) write(*,*)'to determine (minimum) max_reg_list_num values for new grids ...'
!if (do_output()) write(*,*)'triangle_num is ',maxval(triangle_num)

! Invert the temporary data structure. The total number of entries will be 
! the sum of the number of triangle cells for each regular cell. 
total = sum(triangle_num)

! Allocate storage for the final structures in module storage
allocate(triangle_list(total))

! Fill up the long list by traversing the temporary structure. Need indices 
! to keep track of where to put the next entry.
ind = 1

! Loop through each regular grid box
do i = 1, num_reg_x
   do j = 1, num_reg_y

      ! The list for this regular box starts at the current indices.
      triangle_start(i, j) = ind

      ! Copy all the close triangles for regular box(i, j)
      do k = 1, triangle_num(i, j)
         triangle_list(ind) = reg_list(i, j, k)
         ind = ind + 1
      enddo
   enddo
enddo

! Confirm that the indices come out okay as debug
if(ind /= total + 1) then
   string1 = 'Storage indices did not balance: : contact DART developers'
   call error_handler(E_ERR, 'init_interp', string1, source, revision, revdate)
endif

end subroutine init_interp


!------------------------------------------------------------

subroutine get_triangle_corners(ind, lon_corners, lat_corners)

integer,  intent(in)  :: ind
real(r8), intent(out) :: lon_corners(3), lat_corners(3)

integer :: i, cell

! Grabs the corner lons and lats for a given triangle.
do i = 1, 3
   ! Loop through the three cells adjacent to the vertex at the triangle center.
   cell = cellsOnVertex(i, ind)
   ! Get the lats and lons of the centers
   lon_corners(i) = lonCell(cell)
   lat_corners(i) = latCell(cell)
end do

end subroutine get_triangle_corners


!------------------------------------------------------------

subroutine reg_box_overlap(x_corners, y_corners, reg_lon_ind, reg_lat_ind, ier)

! Find a set of regular lat lon boxes that covers all of the area possibley covered by 
! a triangle whose corners are given by the dimension three x_corners 
! and y_corners arrays.  The two dimensional arrays reg_lon_ind and reg_lat_ind
! return the first and last indices of the regular boxes in latitude and
! longitude respectively. These indices may wraparound for reg_lon_ind.  
! A special computation is needed for a triangle that contains one of the poles.
! If the longitude boxes overlap 0
! degrees, the indices returned are adjusted by adding the total number of
! boxes to the second index (e.g. the indices might be 88 and 93 for a case
! with 90 longitude boxes).
! For this version, any triangle that has a vertex at the pole or contains a pole
! returns an error.

real(r8), intent(in)  :: x_corners(3), y_corners(3)
integer,  intent(out) :: reg_lon_ind(2), reg_lat_ind(2)
integer,  intent(out) :: ier

real(r8) :: lat_min, lat_max, lon_min, lon_max
integer  :: i

! Default is success
ier = 0

! Finding the range of latitudes is cake, caveat a pole triangle
lat_min = minval(y_corners)
lat_max = maxval(y_corners)

! For now, will not allow interpolation into triangles that have corners at pole
if (lat_min <= -89.999 .or. lat_max >= 89.999) then
   ier = 1
   return
endif

! Lons are much trickier. Need to make sure to wraparound the
! right way. 
! All longitudes for non-pole rows have to be within 180 degrees
! of one another while a triangle containing the pole will have more
! than a 180 degree span. Need to confirm that there are no exceptions
! to this.

lon_min = minval(x_corners)
lon_max = maxval(x_corners)

if((lon_max - lon_min) > 180.0_r8) then
   ! If the max longitude value is more than 180 
   ! degrees larger than the min, then there must be wraparound or a pole.
   ! Then, find the smallest value > 180 and the largest < 180 to get range.
   lon_min = 360.0_r8
   lon_max = 0.0_r8
   do i=1, 3
      if(x_corners(i) > 180.0_r8 .and. x_corners(i) < lon_min) lon_min = x_corners(i)
      if(x_corners(i) < 180.0_r8 .and. x_corners(i) > lon_max) lon_max = x_corners(i)
   enddo
   ! See if this is a triangle containing the pole.
   ! This happens if the difference after wraparound is also greater than 180 degrees.
   if((360.0_r8 - lon_min) + (lon_max) > 180.0_r8) then
      ! Set the min and max lons and lats for a pole triangle
      lon_min = 0.0_r8
      lon_max = 360.0_r8
      ! North or south pole?
      if(lat_min > 0.0_r8) then
         lat_max = 90.0_r8
      else 
         lat_min = -90.0_r8
      endif
      ! For now will fail on pole overlap
      ier = 1
      return
   endif
endif

! Get the indices for the extreme longitudes
call get_reg_lon_box(lon_min, reg_lon_ind(1))
call get_reg_lon_box(lon_max, reg_lon_ind(2))

! Figure out the indices of the regular boxes for min and max lats
call get_reg_lat_box(lat_min, reg_lat_ind(1))
call get_reg_lat_box(lat_max, reg_lat_ind(2))

! Watch for wraparound again; make sure that second index is greater than first
if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + num_reg_x

end subroutine reg_box_overlap


!------------------------------------------------------------

subroutine update_reg_list(reg_list_num, reg_list, reg_lon_ind, reg_lat_ind, triangle_index)

! Updates the data structure listing dipole quads that are in a given regular box

integer, intent(inout) :: reg_list_num(:, :), reg_list(:, :, :)
integer, intent(inout) :: reg_lon_ind(2), reg_lat_ind(2)
integer, intent(in)    :: triangle_index

integer :: ind_x, index_x, ind_y

! Loop through indices for each possible regular rectangle
! Have to watch for wraparound in longitude
if(reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + num_reg_x

do ind_x = reg_lon_ind(1), reg_lon_ind(2)
   ! Inside loop, need to go back to wraparound indices to find right box
   index_x = ind_x
   if(index_x > num_reg_x) index_x = index_x - num_reg_x

   do ind_y = reg_lat_ind(1), reg_lat_ind(2)
      ! Make sure the list storage isn't full
      if(reg_list_num(index_x, ind_y) >= max_reg_list_num) then
         write(string1,*) 'max_reg_list_num (',max_reg_list_num,') is too small ... increase'
         call error_handler(E_ERR, 'update_reg_list', string1, source, revision, revdate)
      endif

      ! Increment the count
      reg_list_num(index_x, ind_y) = reg_list_num(index_x, ind_y) + 1
      ! Store this quad in the list for this regular box
      reg_list(index_x, ind_y, reg_list_num(index_x, ind_y)) = triangle_index
   enddo
enddo

end subroutine update_reg_list


!------------------------------------------------------------
! new code below here.  nsc 10jan2012
!------------------------------------------------------------

subroutine compute_scalar_with_barycentric(x, loc, ival, dval, ier)
real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: loc
integer,             intent(in)  :: ival
real(r8),            intent(out) :: dval
integer,             intent(out) :: ier

! compute the values at the correct vertical level for each
! of the 3 cell centers defining a triangle that encloses the
! the interpolation point, then interpolate once in the horizontal 
! using barycentric weights to get the value at the interpolation point.

integer, parameter :: listsize = 30 
integer  :: nedges, edgelist(listsize), i, j, neighborcells(maxEdges), edgeid, nvert
real(r8) :: xdata(listsize), ydata(listsize), zdata(listsize)
real(r8) :: t1(3), t2(3), t3(3), r(3), fdata(3), weights(3)
integer  :: vertindex, index1, progindex, cellid, verts(listsize), closest_vert
real(r8) :: lat, lon, vert, tmp(3), fract, lowval(3), uppval(3), p(3)
integer  :: verttype, lower, upper, c(3), vindex, v, vp1
logical  :: inside, foundit


! unpack the location into local vars
tmp = get_location(loc)
lon  = tmp(1)
lat  = tmp(2)
vert = tmp(3)
verttype = nint(query_location(loc))


cellid = find_closest_cell_center(lat, lon)
if (cellid < 1) then
   dval = MISSING_R8
   ier = 11
   return
endif
  
c(1) = cellid

if (on_boundary(cellid)) then
   dval = MISSING_R8
   ier = 11
   return
endif

if (.not. inside_cell0(lat, lon, cellid)) then
   dval = MISSING_R8
   ier = 11
   return
endif

! closest vertex to given point.
closest_vert = closest_vertex_ll(cellid, lat, lon)

! collect the neighboring cell ids and vertex numbers
! this 2-step process avoids us having to read in the
! cellsOnCells() array which i think we only need here.
! if it comes up in more places, we can give up the space
! and read it in and then this is a direct lookup.
! also note which index is the closest vert and later on
! we can start the triangle search there.
vindex = 1
nedges = nEdgesOnCell(cellid)
do i=1, nedges
   edgeid = edgesOnCell(i, cellid)
   if (cellsOnEdge(1, edgeid) /= cellid) then
      neighborcells(i) = cellsOnEdge(1, edgeid)
   else
      neighborcells(i) = cellsOnEdge(2, edgeid)
   endif
   verts(i) = verticesOnCell(i, cellid) 
   if (verts(i) == closest_vert) vindex = i
   call latlon_to_xyz(latCell(neighborcells(i)), lonCell(neighborcells(i)), &
      xdata(i), ydata(i), zdata(i))
enddo


! get the cartesian coordinates in the cell plane for the closest center
call latlon_to_xyz(latCell(cellid), lonCell(cellid), t1(1), t1(2), t1(3))

! and the observation point
call latlon_to_xyz(lat, lon, r(1), r(2), r(3))

! find the cell-center-tri that encloses the obs point
! figure out which way vertices go around cell?
foundit = .false.
findtri: do i=vindex, vindex+nedges
   v = mod(i-1, nedges) + 1
   vp1 = mod(i, nedges) + 1
   t2(1) = xdata(v)
   t2(2) = ydata(v)
   t2(3) = zdata(v)
   t3(1) = xdata(vp1)
   t3(2) = ydata(vp1)
   t3(3) = zdata(vp1)
   call inside_triangle(t1, t2, t3, r, inside, p)
   if (inside) then
      ! p is the xyz of the intersection point in this plane
      ! t2 and t3 are corners, v and vp1 are vert indices
      ! which are same indices for cell centers
      c(2) = neighborcells(v)
      c(3) = neighborcells(vp1)
      foundit = .true.
      exit findtri  
   endif
enddo findtri
if (.not. foundit) then
   dval = MISSING_R8
   ier = 11
   return
endif

! need vert index for the vertical level
call find_vert_level(x, loc, .true., lower, upper, fract, ier)
if (ier /= 0) then
   return
endif

! get the starting index in the state vector
index1 = progvar(ival)%index1
nvert = progvar(ival)%numvertical

! go around triangle and interpolate in the vertical
! t1, t2, t3 are the xyz of the cell centers
! c(3) are the cell ids
do i = 1, 3
   lowval(i) = x(index1 + (c(i)-1) * nvert + lower-1)
   uppval(i) = x(index1 + (c(i)-1) * nvert + upper-1)
   if (vert_is_pressure(loc)) then
      !fdata(i) = exp(log(lowval(i))*(1.0_r8 - fract) + log(uppval(i))*fract)
      fdata(i) = lowval(i)*(1.0_r8 - fract) + uppval(i)*fract
   else
      fdata(i) = lowval(i)*(1.0_r8 - fract) + uppval(i)*fract
   endif
enddo

! now have vertically interpolated values at cell centers.
! get weights and compute value at requested point.
call get_3d_weights(r, t1, t2, t3, weights)

dval = sum(weights * fdata)

ier = 0

end subroutine compute_scalar_with_barycentric

!------------------------------------------------------------

subroutine compute_u_with_rbf(x, loc, zonal, uval, ier)
real(r8),            intent(in)  :: x(:)
type(location_type), intent(in) :: loc
logical,             intent(in)  :: zonal
real(r8),            intent(out) :: uval
integer,             intent(out) :: ier


integer, parameter :: listsize = 30  ! max edges is 10, times 3 cells
logical, parameter :: on_a_sphere = .true.
integer  :: nedges, edgelist(listsize), i, j, nvert
real(r8) :: xdata(listsize), ydata(listsize), zdata(listsize)
real(r8) :: edgenormals(3, listsize)
real(r8) :: veldata(listsize)
real(r8) :: xreconstruct, yreconstruct, zreconstruct
real(r8) :: ureconstructx, ureconstructy, ureconstructz
real(r8) :: ureconstructzonal, ureconstructmeridional
real(r8) :: datatangentplane(3,2)
real(r8) :: coeffs_reconstruct(3,listsize)
integer  :: vertindex, index1, progindex, cellid
real(r8) :: lat, lon, vert, tmp(3), fract, lowval, uppval
integer  :: verttype, lower, upper


! FIXME: make this cache the last value and if the location is
! the same as before and it's asking for V now instead of U,
! skip the expensive computation.

! unpack the location into local vars
tmp = get_location(loc)
lon = tmp(1)
lat = tmp(2)
vert = tmp(3)
verttype = nint(query_location(loc))

call find_surrounding_edges(lat, lon, nedges, edgelist, cellid)
if (nedges <= 0) then
   ! we are on a boundary, no interpolation
   uval = MISSING_R8
   ier = 18
   return
endif
!print *, 'obs lon, lat: ', lon, lat
!print *, 'rbf: surrounding edge list, nedges = ', nedges
!print *, edgelist(1:nedges)
!print *, 'closest cell = ', cellid

! the rbf code needs (their names == our names):
! nData == nedges
! xyz data == xyzEdge
! normalDirectionData == edgeNormalVectors
! velocitydata = U field

! need vert index for the vertical level
call find_vert_level(x, loc, .false., lower, upper, fract, ier)
if (ier /= 0) return

!print *, 'find_vert_level returns l, u, f = ', lower, upper, fract

progindex = get_index_from_varname('u')
if (progindex < 0) then
   ! cannot compute u if it isn't in the state vector
   uval = MISSING_R8
   ier = 18
   return
endif
index1 = progvar(progindex)%index1
nvert = progvar(progindex)%numvertical

do i = 1, nedges
   if (edgelist(i) > size(xEdge)) then
      print *, 'edgelist has index larger than edge count', i, edgelist(i), size(xEdge)
      stop
   endif
   xdata(i) = xEdge(edgelist(i))
   ydata(i) = yEdge(edgelist(i))
   zdata(i) = zEdge(edgelist(i))

   do j=1, 3
      edgenormals(j, i) = edgeNormalVectors(j, edgelist(i))
   enddo

!print *, 'index1, edgelist(i), nVertLevels, lower = ', index1, edgelist(i), nvert, lower
!print *, 'index1, edgelist(i), nVertLevels, upper = ', index1, edgelist(i), nvert, upper
   lowval = x(index1 + (edgelist(i)-1) * nvert + lower-1)
   uppval = x(index1 + (edgelist(i)-1) * nvert + upper-1)
!print *, 'lowval, uppval = ', lowval, uppval
   if (vert_is_pressure(loc)) then
      veldata(i) = lowval*(1.0_r8 - fract) + uppval*fract
      !veldata(i) = exp(log(lowval)*(1.0_r8 - fract) + log(uppval)*fract)
   else
      veldata(i) = lowval*(1.0_r8 - fract) + uppval*fract
   endif
!print *, 'veldata at right vert height for edge: ', i, edgelist(i), veldata(i)
enddo



! get the cartesian coordinates in the cell plane for the reconstruction point
!call latlon_to_xyz_on_plane(lat, lon, cellid, &
!              xreconstruct,yreconstruct,zreconstruct)
! FIXME: on plane is more expensive than just the intersection with the
! sphere.  if it doesn't matter, use this one.
!print *, 'xyz on plane: ', xreconstruct,yreconstruct,zreconstruct
call latlon_to_xyz(lat, lon, xreconstruct,yreconstruct,zreconstruct)
!print *, 'xyz only: ', xreconstruct,yreconstruct,zreconstruct

!! FIXME: DEBUG, remove
!if (lon == 270.00) then
!   print *, xreconstruct, yreconstruct, zreconstruct
!   do i=1, nedges
!      print *, xdata(i), ydata(i), zdata(i)
!   enddo
!endif

! call a simple subroutine to define vectors in the tangent plane
call get_geometry(nedges, xdata, ydata, zdata, &
              xreconstruct, yreconstruct, zreconstruct, edgenormals, &
              on_a_sphere, datatangentplane)

! calculate coeffs_reconstruct
call get_reconstruct_init(nedges, xdata, ydata, zdata, &
              xreconstruct, yreconstruct, zreconstruct, edgenormals, &
              datatangentplane, coeffs_reconstruct)

! do the reconstruction
call get_reconstruct(nedges, lat*deg2rad, lon*deg2rad, &
              coeffs_reconstruct, on_a_sphere, veldata, &
              ureconstructx, ureconstructy, ureconstructz, &
              ureconstructzonal, ureconstructmeridional)

!print *, 'U,V vals from reconstruction: ', ureconstructzonal, ureconstructmeridional
!print *, 'XYZ vals from reconstruction: ', ureconstructx, ureconstructy, ureconstructz

! FIXME: it would be nice to return both and not have to call this
! code twice.  crap.
if (zonal) then
   uval = ureconstructzonal
else
   uval = ureconstructmeridional
endif

ier = 0

end subroutine compute_u_with_rbf

!------------------------------------------------------------

subroutine find_surrounding_edges(lat, lon, nedges, edge_list, cellid)
real(r8), intent(in)  :: lat, lon
integer,  intent(out) :: nedges, edge_list(:)
integer,  intent(out) :: cellid

! given an arbitrary lat/lon location, find the edges of the
! cells that share the nearest vertex.

integer :: vertexid

! find the cell id that has a center point closest
! to the given point.
cellid = find_closest_cell_center(lat, lon)
if (cellid < 1) then
   nedges = 0
   edge_list(:) = -1
   return
endif
   
if (.not. inside_cell(cellid, lat, lon)) then
   nedges = 0
   edge_list(:) = -1
   return
endif

! inside this cell, find the vertex id that the point
! is closest to.
vertexid = closest_vertex_ll(cellid, lat, lon)
if (vertexid <= 0) then
   ! call error handler?  unexpected
   nedges = 0
   edge_list(:) = -1
   return
endif

! fill in the number of unique edges and fills the
! edge list with the edge ids.  the code that detects
! boundary edges for the ocean or regional atmosphere
! is incorporated here.  nedges can come back 0 in that case.
call make_edge_list(vertexid, nedges, edge_list)

end subroutine find_surrounding_edges

!------------------------------------------------------------

subroutine init_closest_center()

! use nCells, latCell, lonCell to initialize a GC structure
! to be used later in find_closest_cell_center().

! set up a GC in the locations mod

integer :: i

allocate(cell_locs(nCells), dummy(nCells), close_ind(nCells))
dummy = 0

do i=1, nCells
   cell_locs(i) = set_location(lonCell(i), latCell(i), 0.0_r8, VERTISSURFACE)
enddo

! FIXME: should be smaller; now slightly less then 1/2 the sphere.
call get_close_maxdist_init(cc_gc, PI* 0.95_r8)  
call get_close_obs_init(cc_gc, nCells, cell_locs)

end subroutine init_closest_center

!------------------------------------------------------------

function find_closest_cell_center(lat, lon)

! Determine the cell index for the closest center to the given point
! 2D calculation only.

real(r8), intent(in)  :: lat, lon
integer               :: find_closest_cell_center

type(location_type) :: pointloc
integer :: i, closest_cell, num_close
real(r8) :: closest_dist, dist
logical, save :: search_initialized = .false.

real(r8) :: l(3)

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_center()
   search_initialized = .true.
endif

pointloc = set_location(lon, lat, 0.0_r8, VERTISSURFACE)

! set up a GC in the locations mod
! call get_close()
! find the closest distance
! return the cell index number

call loc_get_close_obs(cc_gc, pointloc, 0, cell_locs, dummy, num_close, close_ind)
!print *, 'num close = ', num_close

closest_cell = -1
closest_dist = 1.0e9_r8  ! something large in radians
do i=1, num_close
l = get_location(cell_locs(close_ind(i)))
   dist = get_dist(pointloc, cell_locs(close_ind(i)), no_vert = .true.)
!print *, i, close_ind(i), l(1)*deg2rad, l(2)*deg2rad, dist
   if (dist < closest_dist) then
      closest_dist = dist
      closest_cell = close_ind(i)
   endif
enddo
!print *, 'closest ind, dist = ', closest_cell, dist
 
! decide what to do if we don't find anything.
if (closest_cell < 0) then
   if (debug > 5) print *, 'cannot find nearest cell to lon, lat: ', lon, lat
   find_closest_cell_center = -1
   return
endif

! this is the cell index for the closest center
find_closest_cell_center = closest_cell

end function find_closest_cell_center

!------------------------------------------------------------

subroutine finalize_closest_center()

! get rid of storage associated with GC for cell centers.

call get_close_obs_destroy(cc_gc)

end subroutine finalize_closest_center

!------------------------------------------------------------

function find_closest_edge(cellid, lat, lon)

! given a cellid and a lat/lon, find which edge is closest
! to the location.  2D calculation only.

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: lat, lon
integer               :: find_closest_edge

type(location_type) :: pointloc, edgeloc
integer :: i, closest_edge, edgeid
real(r8) :: closest_dist, dist

pointloc = set_location(lon, lat, 0.0_r8, VERTISSURFACE)

closest_edge = -1
closest_dist = 9.0e9  ! something large in radians
nedges = nEdgesOnCell(cellid)
do i=1, nedges
   edgeid = edgesOnCell(i, cellid)
   edgeloc = set_location(lonEdge(edgeid), latEdge(edgeid), 0.0_r8, VERTISSURFACE)
   dist = get_dist(pointloc, edgeloc, no_vert = .true.)
   if (dist < closest_dist) then
      closest_dist = dist
      closest_edge = edgeid
   endif
enddo
 
! decide what to do if we don't find anything.
if (closest_edge < 0) then
   if (debug > 5) print *, 'cannot find nearest edge to lon, lat: ', lon, lat
   find_closest_edge = -1
   return
endif

! this is the edge index for the closest edge center
find_closest_edge = edgeid

end function find_closest_edge

!------------------------------------------------------------

function on_boundary(cellid)

! use the surface (level 1) to determine if any edges (or vertices?)
! are on the boundary, and return true if so.   if the global flag
! is set, skip all code and return false immediately.

integer,  intent(in)  :: cellid
logical               :: on_boundary

! do this completely with topology of the grid.  if any of
! the cell edges are marked as boundary edges, return no.
! otherwise return yes.

integer :: nedges, i, edgeid, vertical

if (global_grid) then
   on_boundary = .false.
   return
endif

! how many edges (same # for verts) to check
nedges = nEdgesOnCell(cellid)

! go around the edges and check the boundary array.
! if any are boundaries, return true.  else, false.

do i=1, nedges
   edgeid = edgesOnCell(i, cellid)

   vertical = 1

   ! FIXME: this is an int array.  is it 0=false,1=true?
   if (boundaryEdge(edgeid, vertical) > 0) then
      on_boundary = .true.
      return
   endif

enddo

on_boundary = .false.

end function on_boundary

!------------------------------------------------------------

function inside_cell(cellid, lat, lon)

! this function no longer really determines if we are inside
! the cell or not.  what it does do is determine if the nearest
! cell is on the grid boundary in any way and says no if it is
! a boundary.  if we have a flag saying this a global grid, we
! can avoid doing any work and immediately return true.  for a
! global atmosphere this is always so; for a regional atmosphere
! and for the ocean (which does not have cells on land) this is
! necessary test.

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: lat, lon
logical               :: inside_cell

! do this completely with topology of the grid.  if any of
! the cell edges are marked as boundary edges, return no.
! otherwise return yes.

integer :: nedges, i, edgeid, vert
real(r8) :: v1(3), v2(3), p(3), vec1(3), vec2(3), r(3), m

! if we're on a global grid, skip all this code
if (global_grid) then
   inside_cell = .true.
   return
endif

nedges = nEdgesOnCell(cellid)

! go around the edges and check the boundary array.
! if any are true, return false.  even if we are inside
! this cell, we aren't going to be able to interpolate it
! so shorten the code path.

! FIXME: at some point we can be more selective and try to
! interpolate iff the edges of the three cells which are
! going to contribute edges to the RBF exist, even if some
! of the other cell edges are on the boundary.  so this
! decision means we won't be interpolating some obs that in
! theory we have enough information to interpolate.  but it
! is conservative for now - we certainly won't try to interpolate
! outside the existing grid.

do i=1, nedges
   edgeid = edgesOnCell(i, cellid)

   ! FIXME: this is an int array.  is it 0=false,1=true?
   ! BOTHER - we need the vert for this and we don't have it
   ! and in fact can't compute it if the interpolation point
   ! has pressure or height as its vertical coordinate.
   vert = 1

   if (boundaryEdge(edgeid, vert) > 0) then
      inside_cell = .false.
      return
   endif

enddo

inside_cell = .true.

end function inside_cell

!------------------------------------------------------------

function inside_cell0(lat, lon, cellid)

! CURRENTLY UNUSED.

! Determine if the given lat/lon is inside the specified cell:
! 1. convert the point to x,y,z on plane defined by cell.
! 2. take the cross product of v1 to the point and v1 to v2
! 3. if the sign is always positive, you're inside.  if it's
!    ever negative, you're outside and you can return.


real(r8), intent(in)  :: lat, lon
integer,  intent(in)  :: cellid
logical               :: inside_cell0

! convert lat/lon to cartesian coordinates and then determine
! from the xyz vertices if this point is inside the cell.

integer :: nverts, i, vertexid
real(r8) :: v1(3), v2(3), p(3), vec1(3), vec2(3), r(3), m

! for the global atmosphere it should always be true. 
! for regional atmosphere and for the ocean, it can be false.
if (global_grid) then
   inside_cell0 = .true.
   return
endif

! cartesian location of point on surface of sphere
call latlon_to_xyz(lat, lon, p(1), p(2), p(3))
!call latlon_to_xyz_on_plane(lat, lon, cellid, p(1), p(2), p(3))

! nedges and nverts is same
nverts = nEdgesOnCell(cellid)

! go around the edges and take the cross product with
! the point.  if all the signs are the same it's inside.
! (or something like this.)
do i=1, nverts
   vertexid = verticesOnCell(i, cellid)
   v1(1) = xVertex(vertexid)
   v1(2) = yVertex(vertexid)
   v1(3) = zVertex(vertexid)
   if (i /= nverts) then
      vertexid = verticesOnCell(i+1, cellid)
   else
      vertexid = verticesOnCell(1, cellid)
   endif
   v2(1) = xVertex(vertexid)
   v2(2) = yVertex(vertexid)
   v2(3) = zVertex(vertexid)

   ! compute the vectors we need here - vertex to point,
   ! v1 to v2, etc.
   vec1 = p - v1
   vec2 = v2 - v1

   call vector_cross_product(vec1, vec2, r)

   call vector_magnitude(r, m)

   if (m < 0.0_r8) then
      inside_cell0 = .false.
      return
   endif
    
enddo

inside_cell0 = .true.

! see also:
! http://tog.acm.org/resources/GraphicsGems/gems/RayPolygon.c

end function inside_cell0

!------------------------------------------------------------

function closest_vertex_ll(cellid, lat, lon)

! Return the vertex id of the closest one to the given point
! this version uses lat/lon.  see closest_vertex_xyz for the
! cartesian version.

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: lat, lon
integer               :: closest_vertex_ll

integer :: nverts, i, closest, vertexid
real(r8) :: distsq, closest_dist, x, y, z, px, py, pz

! use the same radius as MPAS for computing this
call latlon_to_xyz(lat, lon, px, py, pz)

closest_vertex_ll = closest_vertex_xyz(cellid, px, py, pz)
if (closest_vertex_ll < 0) then
   if (debug > 5) print *, 'cannot find nearest vertex to lon, lat: ', lon, lat
endif

end function closest_vertex_ll

!------------------------------------------------------------

function closest_vertex_xyz(cellid, px, py, pz)

! Return the vertex id of the closest one to the given point
! see closest_vertex_ll for the lat/lon version (which calls this)

integer,  intent(in)  :: cellid
real(r8), intent(in)  :: px, py, pz
integer               :: closest_vertex_xyz

integer :: nverts, i, closest, vertexid
real(r8) :: distsq, closest_dist, x, y, z, dx, dy, dz

! nedges and nverts is same in a closed figure
!print *, 'cellid = ', cellid
nverts = nEdgesOnCell(cellid)

closest_dist = 1.0e38   ! something really big; these are meters not radians
closest_vertex_xyz = -1

!print *, 'close: nverts ', nverts
!print *, 'close: point ', px, py, pz
do i=1, nverts
   vertexid = verticesOnCell(i, cellid)
!print *, 'close: v ', i, xVertex(vertexid), yVertex(vertexid), zVertex(vertexid)
   dx = xVertex(vertexid) - px
   dy = yVertex(vertexid) - py
   dz = zVertex(vertexid) - pz
   distsq = (dx * dx) + (dy * dy) + (dz * dz)
!print *, 'close: d ', dx, dy, dx, distsq
   if (distsq < closest_dist) then
      closest_dist = distsq
      closest_vertex_xyz = vertexid
   endif
enddo

end function closest_vertex_xyz

!------------------------------------------------------------

subroutine make_edge_list(vertexid, nedges, edge_list)

! FIXME: will need a vertical level number input arg here to detect
! the boundary edges/vertices correctly.

! given a vertexid, look up the N cells which share it as
! a vertex, and then for those cells look up the edge ids.
! return a list with all edges listed exactly once (shared edges 
! are detected and not replicated in the output list).  
! the edge_list output should be at least 10x the Ncells to 
! guarentee it will be large enough if all cells are disjoint.

integer, intent(in)  :: vertexid
integer, intent(out) :: nedges, edge_list(:)

integer :: edgecount, c, e, listlen, l, nextedge
integer :: ncells, cellid_list(3)
logical :: found

! use the cellsOnVertex() array to find the three cells
! which share this vertex.  note that if we wanted to change
! the number of edges - for example only return the 3 immediate
! edges to this vertex, or if we wanted to expand the region
! to include the ~12 second-nearest-cell-neighbors, you can
! change this code here.
ncells = 3
cellid_list(1) = cellsOnVertex(1, vertexid)
cellid_list(2) = cellsOnVertex(2, vertexid)
cellid_list(3) = cellsOnVertex(3, vertexid)

! use nEdgesOnCell(nCells) and edgesOnCell(nCells, 10) to
! find the edge numbers.  add them to the list, skipping if
! the edge is already there.  increment nedges each time a
! new edge is added.  check arrays for enough length before
! starting to work.


! FIXME: the ocean files have:
!  integer boundaryEdge(nVertLevels, nEdges)
!  integer boundaryVertex(nVertLevels, nVertices)
! as a first pass, if ANY of the edges or vertices are on
! the boundary, punt and return 0 as the edge count.  later
! once this is working, decide if a single boundary vertex or
! edge is ok if it's the exterior of the edges we are including
! and if it has good data values.

listlen = 0
do c=1, ncells
   edgecount = nEdgesOnCell(cellid_list(c))
   do e=1, edgecount
      nextedge = edgesOnCell(e, cellid_list(c))
      ! FIXME: 
      ! if (boundaryEdge(nextedge, vert)) then
      !    nedges = 0
      !    edge_list(:) = -1
      !    return
      ! endif
      found = .false.
      addloop: do l=1, listlen
         if (edge_list(l) == nextedge) then
            found = .true.
            exit addloop
         endif
      enddo addloop
      if ( .not. found) then
         listlen = listlen + 1
         edge_list(listlen) = nextedge
      endif
   enddo
enddo

nedges = listlen

end subroutine make_edge_list

!------------------------------------------------------------

subroutine latlon_to_xyz(lat, lon, x, y, z)

! Given a lat, lon in degrees, return the cartesian x,y,z coordinate 
! on the surface of a specified radius relative to the origin 
! at the center of the earth.  (this radius matches the one
! used at MPAS grid generation time and must agree in order
! to be consistent with the cartisian coordinate arrays in
! the MPAS data files.)

real(r8), intent(in)  :: lat, lon
real(r8), intent(out) :: x, y, z

real(r8) :: rlat, rlon

rlat = lat * deg2rad
rlon = lon * deg2rad

x = radius * cos(rlon) * cos(rlat)
y = radius * sin(rlon) * cos(rlat)
z = radius * sin(rlat)

end subroutine latlon_to_xyz

!------------------------------------------------------------

subroutine xyz_to_latlon(x, y, z, lat, lon)

! Given a cartesian x, y, z coordinate relative to the origin
! at the center of the earth, using a fixed radius specified
! by MPAS (in the grid generation step), return the corresponding
! lat, lon location in degrees.

real(r8), intent(in)  :: x, y, z
real(r8), intent(out) :: lat, lon

real(r8) :: rlat, rlon

! right now this is only needed for debugging messages.
! the arc versions of routines are expensive.

rlat = PI/2.0_r8 - acos(z/radius)
rlon = atan2(y,x)
if (rlon < 0) rlon = rlon + PI*2

lat = rlat * rad2deg
lon = rlon * rad2deg

end subroutine xyz_to_latlon

!------------------------------------------------------------

subroutine inside_triangle(t0, t1, t2, s, inside, intp)

! given 3 corners of a triangle and an xyz point, compute
! whether the ray from the origin to s intersects the triangle
! in the plane defined by the vertices.  sets t/f flag for inside
! and returns intersection point in xyz if true.

! Uses the parametric form description from:
!  http://en.wikipedia.org/wiki/Line-plane_intersection
! in this case, point a is the origin (0,0,0) point b is 's',
! and points 0,1,2 are t0,t1,t2.  the vector [t,u,v] is 'v'.

real(r8), intent(in)  :: t0(3), t1(3), t2(3)
real(r8), intent(in)  :: s(3)
logical,  intent(out) :: inside
real(r8), intent(out) :: intp(3)

real(r8) :: m(3,3), v(3) ! intermediates to compute intersection
real(r8) :: mi(3,3)      ! invert of m
integer  :: i


m(1,1) = -s(1)
m(1,2) = t1(1) - t0(1)
m(1,3) = t2(1) - t0(1)

m(2,1) = -s(2)
m(2,2) = t1(2) - t0(2)
m(2,3) = t2(2) - t0(2)

m(3,1) = -s(3)
m(3,2) = t1(3) - t0(3)
m(3,3) = t2(3) - t0(3)

call invert3(m, mi)

m = matmul(m, mi)

v = matmul(mi, -t0) 

intp = 0.0_r8
inside = .false.

! first be sure the triangle intersects the line
! between [0,1].  compute the intersection point.
! then test that v(2) and v(3) are both between 
! [0,1], and that v(2)+v(3) <= 1.0 
! if all true, intersection pt is inside tri.
if (v(1) >= 0.0_r8 .and. v(1) <= 1.0_r8) then
   intp(1) = s(1) * v(1)
   intp(2) = s(2) * v(1)
   intp(3) = s(3) * v(1)

   if ((v(2) >= 0.0_r8 .and. v(2) <= 1.0_r8) .and. &
       (v(3) >= 0.0_r8 .and. v(3) <= 1.0_r8) .and. &
       (v(2)+v(3) <= 1.0_r8)) inside = .true.

endif

end subroutine inside_triangle

!------------------------------------------------------------

subroutine latlon_to_xyz_on_plane(lat, lon, cellid, x, y, z)

! Given a lat, lon in degrees, and the id of a cell in the 
! MPAS grid, return the cartesian x,y,z coordinate of that
! location ON THE PLANE defined by the vertices of that cell.
! This will be different from the x,y,z of the surface of the
! sphere.  Uses the parametric form description from
! http://en.wikipedia.org/wiki/Line-plane_intersection

real(r8), intent(in)  :: lat, lon
integer,  intent(in)  :: cellid
real(r8), intent(out) :: x, y, z

integer  :: nverts, i, vertexid
real(r8) :: s(3)         ! location of point on surface
real(r8) :: p(3,3)       ! first 3 vertices of cell, xyz
real(r8) :: intp(3)      ! intersection point with plane
logical  :: inside       ! if true, intersection is inside tri


call latlon_to_xyz(lat, lon, s(1), s(2), s(3))

! get the first 3 vertices to define plane
! intersect with sx,sy,sz to get answer

! nedges and nverts is same
nverts = nEdgesOnCell(cellid)
if (nverts < 3) then
   print *, 'nverts is < 3', nverts
   stop
endif

! use first 3 verts to define plane
do i=1, 3
   vertexid = verticesOnCell(i, cellid)

   p(1,i) = xVertex(vertexid)
   p(2,i) = yVertex(vertexid)
   p(3,i) = zVertex(vertexid)
enddo

! in this case we don't care about inside this tri, just where
! it intersects the plane.
call inside_triangle(p(:,1), p(:,2), p(:,3), s, inside, intp)

x = intp(1)
y = intp(2)
z = intp(3)

end subroutine latlon_to_xyz_on_plane

!------------------------------------------------------------

subroutine vector_magnitude(a, r)

! Given a cartesian vector, compute the magnitude

real(r8), intent(in)  :: a(3)
real(r8), intent(out) :: r

r = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

end subroutine vector_magnitude

!------------------------------------------------------------

subroutine vector_cross_product(a, b, r)

! Given 2 cartesian vectors, compute the cross product of a x b

real(r8), intent(in)  :: a(3), b(3)
real(r8), intent(out) :: r(3)

r(1) = a(2)*b(3) - a(3)*b(2)
r(2) = a(3)*b(1) - a(1)*b(3)
r(3) = a(1)*b(2) - a(2)*b(1)

end subroutine vector_cross_product

!------------------------------------------------------------

subroutine determinant3(a, r)

! Given a 3x3 matrix, compute the determinant

real(r8), intent(in)  :: a(3,3)
real(r8), intent(out) :: r

r = a(1,1)*(a(2,2)*a(3,3) - (a(3,2)*a(2,3))) + &
    a(2,1)*(a(3,2)*a(1,3) - (a(3,3)*a(1,2))) + &
    a(3,1)*(a(1,2)*a(2,3) - (a(2,2)*a(1,3)))

end subroutine determinant3

!------------------------------------------------------------

subroutine invert3(a, r)

! Given a 3x3 matrix, compute the inverse

real(r8), intent(in)  :: a(3,3)
real(r8), intent(out) :: r(3,3)

real(r8) :: det, b(3,3)

call determinant3(a, det)
if (det == 0.0_r8) then
   print *, 'matrix cannot be inverted'
   r = 0.0_r8
   return
endif

b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
b(2,1) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)

b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
b(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)

b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
b(2,3) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)

r = b / det

end subroutine invert3

!------------------------------------------------------------

!==================================================================
! The following (private) routines were borrowed from the MPAS code
!==================================================================

!------------------------------------------------------------------

subroutine uv_cell_to_edges(zonal_wind, meridional_wind, du)

! Project u, v wind increments at cell centers onto the edges.
! FIXME:
!        we can hard-code R3 here since it comes from the (3d) x/y/z cartesian coordinate.
!        We define nEdgesOnCell in get_grid_dims, and read edgesOnCell in get_grid.
!        We read edgeNormalVectors in get_grid to use this subroutine.
!        Here "U" is the prognostic variable in MPAS, and we update it with the wind
!        increments at cell centers.

real(r8), intent(in) :: zonal_wind(:,:)             ! u wind updated from filter
real(r8), intent(in) :: meridional_wind(:,:)        ! v wind updated from filter
real(r8), intent(out):: du(:,:)                     ! normal velocity increment on the edges

! Local variables
integer, parameter :: R3 = 3
real(r8) :: east(R3,nCells), north(R3,nCells)
real(r8) :: lonCell_rad(nCells), latCell_rad(nCells)
integer  :: iCell, iEdge, jEdge, k

if ( .not. module_initialized ) call static_init_model

! Initialization
du(:,:) = 0.0_r8

! Back to radians (locally)
lonCell_rad = lonCell*deg2rad
latCell_rad = latCell*deg2rad

! Compute unit vectors in east and north directions for each cell:
do iCell = 1, nCells

    east(1,iCell) = -sin(lonCell_rad(iCell))
    east(2,iCell) =  cos(lonCell_rad(iCell))
    east(3,iCell) =  0.0_r8
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell_rad(iCell))*sin(latCell_rad(iCell))
    north(2,iCell) = -sin(lonCell_rad(iCell))*sin(latCell_rad(iCell))
    north(3,iCell) =  cos(latCell_rad(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

enddo

if ( debug > 7 ) then
write(*,*)
write(*,*)'u,v incr before uv_cell_to_edges:',minval(zonal_wind), maxval(zonal_wind), &
                                              minval(meridional_wind), maxval(meridional_wind)
endif

! Project analysis increments from the cell centers to the edges
do iCell = 1, nCells
   do jEdge = 1, nEdgesOnCell(iCell)
      iEdge = edgesOnCell(jEdge, iCell)
      do k = 1, nVertLevels
         du(k,iEdge) = du(k,iEdge) + 0.5_r8 * zonal_wind(k,iCell)   &
                     * (edgeNormalVectors(1,iEdge) * east(1,iCell)  &
                     +  edgeNormalVectors(2,iEdge) * east(2,iCell)  &
                     +  edgeNormalVectors(3,iEdge) * east(3,iCell)) &
                     + 0.5_r8 * meridional_wind(k,iCell)            &
                     * (edgeNormalVectors(1,iEdge) * north(1,iCell) &
                     +  edgeNormalVectors(2,iEdge) * north(2,iCell) &
                     +  edgeNormalVectors(3,iEdge) * north(3,iCell)) 
      enddo
   enddo
enddo

if ( debug > 7 ) then
write(*,*)
write(*,*)'du inside uv_cell_to_edges:',minval(du), maxval(du)
endif

end subroutine uv_cell_to_edges


!------------------------------------------------------------------

subroutine r3_normalize(ax, ay, az)

!normalizes the vector (ax, ay, az)

real(r8), intent(inout) :: ax, ay, az
real(r8) :: mi

 mi = 1.0_r8 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

end subroutine r3_normalize


!------------------------------------------------------------------

function theta_to_tk (theta, rho, qv)

! Compute sensible temperature [K] from potential temperature [K].
! matches computation done in MPAS model

real(r8), intent(in)  :: theta    ! potential temperature [K]
real(r8), intent(in)  :: rho      ! dry density
real(r8), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
real(r8)  :: theta_to_tk          ! sensible temperature [K]

! Local variables
real(r8) :: theta_m               ! potential temperature modified by qv
real(r8) :: exner                 ! exner function

theta_m = (1.0_r8 + 1.61_r8 * qv)*theta
exner = ( (rgas/p0) * (rho*theta_m) )**rcv

! Temperature [K]
theta_to_tk = theta * exner

end function theta_to_tk


!------------------------------------------------------------------

subroutine compute_full_pressure(theta, rho, qv, pressure, tk)

! Compute full pressure from the equation of state.
! matches computation done in MPAS model

real(r8), intent(in)  :: theta    ! potential temperature [K]
real(r8), intent(in)  :: rho      ! dry density
real(r8), intent(in)  :: qv       ! water vapor mixing ratio [kg/kg]
real(r8), intent(out) :: pressure ! full pressure [Pa]
real(r8), intent(out) :: tk       ! return sensible temperature to caller

tk = theta_to_tk(theta, rho, qv)
pressure = rho * rgas * tk * (1.0_r8 + 1.61_r8 * qv)

!print *, 'compute_full_pressure: ', pressure, theta, rho, qv, tk
end subroutine compute_full_pressure



!===================================================================
! End of model_mod
!===================================================================
end module model_mod
