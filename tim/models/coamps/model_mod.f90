! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module model_mod

!------------------------------
! MODULE:       model_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module for the interface between DART and the U. S. Navy's COAMPS
! mesoscale model.  COAMPS was developed by the Naval Research Laboratory,
! Monterey, California.  COAMPS is a registered trademark of the Naval
! Research Laboratory.
!------------------------------ 

  use coamps_grid_mod,     only : coamps_grid,                  &
                                  dump_grid_info,               &
                                  nc_dump_grid_info,            &
                                  get_grid_dims,                &
                                  get_grid_dsigmaw,             &
                                  get_grid_field_size,          &
                                  get_grid_msigma,              &
                                  get_grid_num_levels,          &
                                  get_grid_wsigma,              &
                                  get_terrain_height_at_points, &
                                  gridpt_to_latlon,             &
                                  location_to_gridpt
                              
  use coamps_interp_mod,   only : interpolate,                  &
                                  set_interp_diag
  use coamps_restart_mod,  only : dump_restart_vars,            &
                                  get_num_vars,                 &
                                  get_pert_magnitude_by_index,  &
                                  get_pert_type_by_index,       &
                                  get_restart_grid,             &
                                  get_var_type_by_index,        &
                                  get_vert_coord_by_index,      &
                                  initialize_state_info,      &
                                  PERT_TYPE_INDIVID,            &  
                                  PERT_TYPE_NOPERTS,            &
                                  PERT_TYPE_UNIFORM
  use coamps_util_mod,     only : check_alloc_status,           &
                                  check_dealloc_status

  use location_mod,        only : get_close_type,                &
                                  get_dist,                      &
                                  get_location,                  &
                                  location_type,                 &
                                  loc_get_close_maxdist_init  => &
                                      get_close_maxdist_init,    &
                                  loc_get_close_obs           => &
                                      get_close_obs,             &
                                  loc_get_close_obs_init      => &
                                      get_close_obs_init,        &
                                  set_location,                  &
                                  vert_is_pressure,              &
                                  horiz_dist_only,               &
                                  query_location,                &
                                  VERTISLEVEL
  use netcdf
  use typeSizes
  use obs_kind_mod
  use random_seq_mod,      only : init_random_seq,               &
                                  random_gaussian,               &
                                  random_seq_type
  use time_manager_mod,    only : set_time,                      &
                                  time_type,                     &
                                  print_time,                    &
                                  print_date,                    &
                                  set_time_missing
  use types_mod,           only : MISSING_R8,                    &
                                  r8
  use utilities_mod,       only : check_namelist_read,           &
                                  do_output,                     &
                                  E_ERR,                         &
                                  E_MSG,                         &
                                  error_handler,                 &
                                  find_namelist_in_file,         &
                                  logfileunit,                   &
                                  get_unit,                      &
                                  nc_check,                      &
                                  register_module,               &
                                  do_nml_file,                   &
                                  do_nml_term,                   &
                                  nmlfileunit,                   &
                                  find_textfile_dims,            &
                                  file_to_text

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------
  
  ! Initialization/finalization
  public :: static_init_model 
  public :: end_model      

  ! NetCDF
  public :: nc_write_model_atts 
  public :: nc_write_model_vars 

  ! Ensemble generation
  public :: pert_model_state 

  ! Forward operator
  public :: model_interpolate
  public :: ens_mean_for_model

  ! Localization
  public :: get_close_maxdist_init
  public :: get_close_obs_init
  public :: get_close_obs

  ! Information about model setup
  public :: get_model_size 
  public :: get_state_meta_data 
  public :: get_model_time_step 

  ! Null interfaces
  public :: init_conditions
  public :: init_time      
  public :: adv_1step 

  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACE
  !------------------------------
  !  [none]
  !------------------------------
  ! END EXTERNAL INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS
  !------------------------------

  type state_vec_iterator
     private
     integer  :: nitems
     integer  :: cur_item
     character(len=128),dimension(10) :: something
     real(r8) :: bob
     logical  :: more
  end type state_vec_iterator

  INTERFACE get_var 
      MODULE PROCEDURE get_var_3d
      MODULE PROCEDURE get_var_2d
      MODULE PROCEDURE get_var_1d
  END INTERFACE


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

  character(len=128) :: msgstring   ! general purpose string for printing

  ! Main model_mod namelist - not too much here as we read most of
  ! the data we need in from the COAMPS files themselves
  character(len=10) :: cdtg = '1999083100' ! Date-time group
  integer           :: y_bound_skip = 3    ! How many x and y boundary
  integer           :: x_bound_skip = 3    ! points to skip when
                                           ! perturbing the model
                                           ! state
  logical           :: need_mean = .true.  ! Do we need the ensemble
                                           ! mean for for forward
                                           ! operator computation?
  character(len=80) :: dsnrff = './'       ! Path to data files
  logical           :: output_interpolation = .false. 
  logical           :: output_state_vector = .false.
  integer           :: debug = 0    ! turn up for more and more debug messages

  namelist /model_nml/ cdtg, y_bound_skip, x_bound_skip, need_mean, dsnrff, &
                       output_interpolation, output_state_vector, debug

  ! Locations of state variables
  type(location_type), dimension(:), allocatable :: all_locs
  integer,             dimension(:), allocatable :: all_kinds

  ! Grid information structure
  type(coamps_grid) :: restart_grid

  ! Ensemble mean
  real(r8), dimension(:), allocatable :: ensemble_mean

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------

  ! static_init_model
  ! -----------------
  ! One-time initialization of the model.  For COAMPS, this:
  !  1. Reads in the model_mod namelist
  !  2. Initializes the pressure levels for the state vector
  !  3. Generate the location data for each member of the state
  !  PARAMETERS
  !   [none]
  subroutine static_init_model()

    integer :: nml_unit
    
    character(len=*), parameter :: routine = 'static_init_model'
    integer :: io_status, alloc_status

    integer :: num_vars, max_i, max_j

    integer :: state_ii, ii,jj, var_index

    real(kind=r8) :: lon, lat
    real(kind=r8) :: vert_loc
    integer       :: vert_type


    call register_module(source, revision, revdate)

    ! Read in the namelist information
    call find_namelist_in_file('input.nml', 'model_nml', nml_unit)
    read (nml_unit,nml=model_nml,iostat=io_status)
    call check_namelist_read(nml_unit, io_status, 'model_nml')
    if (do_nml_file()) write (nmlfileunit,model_nml)
    if (do_nml_term()) write (     *     ,model_nml)

    call initialize_state_info(cdtg, 'state.vars', dsnrff)

    call set_interp_diag(output_interpolation)

    call get_restart_grid(restart_grid)
    if (do_output()) call dump_grid_info(restart_grid)
    
    allocate( all_locs(get_model_size()), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'all_locs')
    allocate( all_kinds(get_model_size()), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'all_kinds')
    if (need_mean) then
       allocate( ensemble_mean(get_model_size()), stat=alloc_status )
       call check_alloc_status(alloc_status, routine, source, revision, &
                               revdate, 'ensemble_mean')
    end if

    ! Calculate loop limits
    call get_num_vars(num_vars)
    call get_grid_dims(restart_grid, max_i, max_j)

    ! Location is easy provided the i/j coordinate and
    ! the sigma level, and the type is given by the restart file
    ! information structure
    state_ii = 1
    do var_index = 1, num_vars
       do jj = 1, max_j
          do ii =1, max_i
             call gridpt_to_latlon(grid = restart_grid,    &
                                   ir   = real(ii, kind=r8), &
                                   jr   = real(jj, kind=r8), &
                                   lat  = lat, lon=lon)

             ! Store the variable type for easy access from get_state_meta_data.
             call get_var_type_by_index(var_index, all_kinds(state_ii))

             ! Vertical coordinate information is based on the 
             ! state vector definition since it's a field-wide constant
             call get_vert_coord_by_index(var_index, vert_type, vert_loc)
             all_locs(state_ii) = set_location(lon, lat, vert_loc, vert_type)

             state_ii = state_ii + 1
          enddo
       enddo
    enddo

  end subroutine static_init_model

  ! end_model
  ! ---------
  ! Clean up the workspace once the program is ending
  !  PARAMETERS
  !   [none]
  subroutine end_model()

    character(len=*), parameter :: routine = 'end_model'
    integer                     :: dealloc_status

    deallocate(all_locs, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'all_locs')
    deallocate(all_kinds, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'all_locs')
    if (need_mean) then
       deallocate(ensemble_mean, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source,       &
                                  revision, revdate, 'all_locs')
    end if
  end subroutine end_model



  ! nc_write_model_atts
  ! -------------------
  ! Write model-specific global attributes to a NetCDF file.
  !  PARAMETERS
  !   IN  ncFileID          numeric ID of an *open* NetCDF file
  !   OUT ierr              0 if writing was successful
  function nc_write_model_atts( ncFileID ) result (ierr)
    integer, intent(in) :: ncFileID      ! netCDF file 
    integer             :: ierr          

    integer :: grid_i_size
    integer :: grid_j_size
 
    ! NetCDF variables                                               
    integer :: n_dims, n_vars, n_atts  ! NetCDF counts               
    integer :: time_dimid              ! Time and state dimensions   
    integer :: mem_dimid               ! Ensemble member dimension   
    integer :: ulim_dimid              ! unlimited dimension         

    ! variables for the namelist output

    character(len=129), allocatable, dimension(:) :: textblock
    integer :: LineLenDimID, nlinesDimID, nmlVarID
    integer :: nlines, linelen
    logical :: has_coamps_namelist

    ! Date and time
    integer, dimension(8)        :: dt_values 
    character(len=NF90_MAX_NAME) :: dt_string 

    call get_grid_dims(restart_grid, grid_i_size, grid_j_size)

    ! Default to no errors                                           
    ierr = 0 

    ! Ensure that we're dealing with an open & current  NetCDF file
    call nc_check(nf90_inquire(ncFileID, n_dims, n_vars, n_atts, ulim_dimid), routine)
    call nc_check(nf90_sync(ncFileID), routine)

    ! Go into define mode
    call nc_check(nf90_redef(ncFileID), routine)

    ! Find the ensemble member and time dimensions
    call nc_check(nf90_inq_dimid(ncFileID, 'copy',  mem_dimid), routine)
    call nc_check(nf90_inq_dimid(ncFileID, 'time', time_dimid), routine)

    ! Get the 'creation date' using the intrinsic F90 function
    call DATE_AND_TIME(values=dt_values)
    write (dt_string, '(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') & 
                          dt_values(1), dt_values(2), dt_values(3), &
                          dt_values(5), dt_values(6), dt_values(7) 
 
    ! Global attributes                                   
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'creation_date',  dt_string), routine) 
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'model_source',   source),    routine) 
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'model_revision', revision),  routine)  
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'model_revdate',  revdate),   routine) 
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'model',          'COAMPS'),  routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'grid_x',      grid_i_size), routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,     &
                  'grid_y',      grid_j_size), routine)

    call nc_dump_grid_info(restart_grid, ncFileID)

    !---------------------------------------------------------------------------
    ! Determine shape of COAMPS namelist - if it is called 'namelist'
    !---------------------------------------------------------------------------

    call find_textfile_dims('namelist', nlines, linelen)
    if (nlines > 0) then
      has_coamps_namelist = .true.
    else
      has_coamps_namelist = .false.
    endif

    if (debug > 0) print *, 'COAMPS namelist: nlines, linelen = ', nlines, linelen

    if (has_coamps_namelist) then
       allocate(textblock(nlines))
       textblock = ''

       call nc_check(nf90_def_dim(ncFileID, 'nlines', nlines, nlinesDimID), &
                           routine, 'def_dim nlines ')
    
       call nc_check(nf90_def_var(ncFileID, 'COAMPSnamelist', nf90_char, &
                          (/ linelenDimID, nlinesDimID /), nmlVarID),    &
                                 routine, 'def_var coamps_namelist')
       call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                     'contents of COAMPS namelist'), routine, 'put_att coamps_namelist')
    endif

    !-------------------------------------------------------------------------------
    ! Here is the extensible part. The simplest scenario is to output the state vector,
    ! parsing the state vector into model-specific parts is complicated, and you need
    ! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
    ! complicated part.
    !-------------------------------------------------------------------------------

    if ( output_state_vector ) then
       call nc_write_statearray_atts( ncFileID )
    else
       call nc_write_prognostic_atts( ncFileID )
    endif

    ! Need to get out of define mode to fill the variables 
    call nc_check(nf90_enddef(ncFileID), routine) 
 
    !---------------------------------------------------------------------------
    ! Fill the variables we can
    !---------------------------------------------------------------------------

    if (has_coamps_namelist) then
       call file_to_text('namelist', textblock)
       call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock), &
                          routine, 'put_var nmlVarID')
       deallocate(textblock)
    endif

    ! Sync the NetCDF file on disk, but don't close 
    call nc_check(nf90_sync(ncFileID), routine) 
    if (do_output()) then
       write (*,*) 'state_vector attributes written and file synced...' 
    end if

  end subroutine nc_write_state_vector_atts



  ! nc_write_prognostic_atts
  ! -------------------
  ! Write prognostic-variable -specific attributes to a NetCDF file.
  !  PARAMETERS
  !   IN  ncFileID          numeric ID of an *open* NetCDF file

  subroutine nc_write_prognostic_atts( ncFileID )
    integer, intent(in)         :: ncFileID      ! netCDF file 
    integer                     :: ierr          

    integer                     :: model_size 
    integer                     :: grid_i_size
    integer                     :: grid_j_size
    type(location_type)         :: loc        
    real(kind=r8), dimension(3) :: loc3d      ! lon/lat/hgt
    integer                     :: vartype    ! state variable type
    integer                     :: ii         
 
    ! Sync the NetCDF file, but don't close it                                                       
    call nc_check(nf90_sync(ncFileID), routine) 
    if (do_output()) write (*,*) 'Model attributes written and file synced...' 
    

  end subroutine nc_write_prognostic_atts



  ! nc_write_model_vars
  ! -------------------
  ! Writes the model variables (for now just the whole state vector)
  ! to the NetCDF file
  !  PARAMETERS
  !   IN  ncFileID          numeric *open* NetCDF file ID
  !   IN  statevec          large DART state vector
  !   IN  copyindex         which 'copy' to write - this is used for
  !                         indexing ensemble members
  !   IN  timeindex         which time to write (as an index)
  !   OUT ierr              0 if successful
  function nc_write_model_vars(ncFileID, statevec, copyindex, &
                               timeindex ) result (ierr)          
    integer,                intent(in) :: ncFileID
    real(r8), dimension(:), intent(in) :: statevec
    integer,                intent(in) :: copyindex
    integer,                intent(in) :: timeindex
    integer                            :: ierr   

    type(state_vec_iterator)     :: iterator
    integer, dimension(3)        :: var_dims
    integer                      :: ndims 
    character(len=NF90_MAX_NAME) :: var_name
    character(len=NF90_MAX_NAME), dimension(10) :: variable_list
    integer                      :: alloc_status
    integer                      :: varid 

    real(r8), allocatable, dimension(:,:,:) :: var3d
    real(r8), allocatable, dimension(:,:)   :: var2d
    real(r8), allocatable, dimension(:)     :: var1d

    ! Error handling
    character(len=*), parameter :: routine = 'nc_write_model_vars'

    ! Default to no errors - since any errors encountered here will
    ! just throw us out of the program
    ierr = 0

    if ( output_state_vector ) then

       ! Get the numerical index of the state variable, then dump the
       ! state vector to the corresponding ensemble member and time
       call nc_check(nf90_inq_varid(ncFileID, 'state', varid), &
                           routine, 'inq_varid state') 
       call nc_check(nf90_put_var(ncFileID, varid, statevec,   &
                                 start=(/ 1,copyindex,timeindex /)), &
                           routine, 'put_var statevec')

    else

       iterator = get_iterator(variable_list)    
       output_vars:  do while(has_next(iterator))

          call get_var_name(var_name,        iterator)
          call get_var_dims(var_dims, ndims, iterator)

          call nc_check(nf90_inq_varid(ncFileID, trim(var_name), varid), &
                          routine, 'inq_varid '//trim(var_name)) 

          if      ( ndims == 1 ) then

             allocate(var1d(var_dims(1)), stat=alloc_status)
             call check_alloc_status(alloc_status, routine, source, revision,   &
                            revdate, '1D allocate '//trim(var_name))
             call get_var(statevec, var1d, iterator)
             call nc_check(nf90_put_var(ncFileID, varid, var1d,                 &
                                      start=(/ 1,copyindex,timeindex /)),       &
                                      routine, 'put_var '//trim(var_name))

             deallocate(var1d, stat=alloc_status)
             call check_dealloc_status(alloc_status, routine, source, revision, &
                              revdate, '1D deallocate '//trim(var_name))

          else if ( ndims == 2 ) then

             allocate(var2d(var_dims(1),var_dims(2)), stat=alloc_status)
             call check_alloc_status(alloc_status, routine, source, revision,   &
                            revdate, '2D allocate '//trim(var_name))
             call get_var(statevec, var2d, iterator)
             call nc_check(nf90_put_var(ncFileID, varid, var2d,                 &
                                    start=(/ 1,1,copyindex,timeindex /)),       &
                                    routine, 'put_var '//trim(var_name))

             deallocate(var2d, stat=alloc_status)
             call check_dealloc_status(alloc_status, routine, source, revision, &
                              revdate, '2D deallocate '//trim(var_name))

          else if ( ndims == 3 ) then

             allocate(var3d(var_dims(1),var_dims(2),var_dims(3)), stat=alloc_status)
             call check_alloc_status(alloc_status, routine, source, revision,   &
                            revdate, '3D allocate '//trim(var_name))
             call get_var(statevec, var3d, iterator)
             call nc_check(nf90_put_var(ncFileID, varid, var3d,                 &
                                    start=(/ 1,1,1,copyindex,timeindex /)),     &
                                    routine, 'put_var '//trim(var_name))

             deallocate(var3d, stat=alloc_status)
             call check_dealloc_status(alloc_status, routine, source, revision, &
                              revdate, '3D deallocate '//trim(var_name))

          else

             write(msgstring,*)'not built to handle ndims == ',ndims
             call error_handler(E_ERR, routine, msgstring, source, revision, revdate) 

          endif

          call get_next(iterator)

       end do output_vars
    endif
 
    ! Flush the buffer to disk                                       
    call nc_check(nf90_sync(ncFileId), routine) 

  end function nc_write_model_vars

  ! pert_model_state
  ! ----------------
  ! Perturb the model state, field by field.  This can be done 3
  ! different ways:
  !  1. No perturbation
  !  2. Uniform perturbation - each element of the field has the
  !                            same additive perturbation
  !  3. Individual perturbation - each element of the field has a 
  !                               different additive perturbation
  ! The perturbation magnitude and option are supplied out of the
  ! dynamic restart vector definition - this allows us to supply a
  ! variance appropriate for each type of variable at each level.
  !  PARAMETERS
  !   IN  state             DART state vector
  !   OUT pert_state        state vector after perturbations
  !   OUT interf_provided   true if this routine did the perturbation
  subroutine pert_model_state(state, pert_state, interf_provided)
    ! Avoid having to link *everything* with the MPI utilities
    use mpi_utilities_mod, only : my_task_id

    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided

    type(random_seq_type), save :: random_sequence
    logical,               save :: initialized_random_sequence = .false. 
    integer,               save :: rand_seq_seed

    integer, parameter          :: DT_MILLISECOND = 8
    integer, dimension(8), save :: datetime

    real(r8) :: pert
    integer  :: ii,jj,kk
    integer  :: ijarea,gridii,gridjj
    integer  :: fullii
    integer  :: num_vars

    character(len=*), parameter :: routine = 'pert_model_state'
    integer                     :: alloc_status, dealloc_status
    
    ! Base the perturbation on the state.vars file
    real(kind=r8) :: pertsize
    integer       :: perttype

    ! Do the perturbations based on a 2-D grid: this'll make it
    ! easier to skip boundaries
    real(r8), dimension(:,:), allocatable :: temp_state
    real(r8), dimension(:,:), allocatable :: temp_pert

    if (.not. initialized_random_sequence) then

        ! Seed the random number generator based on the number
        ! SSmmmP, where:
        !  SS  = current time's seconds value
        !  mmm = current time's milliseconds value
        !  P   = current process ID
        ! This should give us unique perturbations between different
        ! processors and different runs of the software.
        call DATE_AND_TIME(values=datetime)
        rand_seq_seed = 1E1 * datetime(DT_MILLISECOND) + & 
                        1E0 * my_task_id()
        call init_random_seq(r=random_sequence, seed=rand_seq_seed)
    end if

    call error_handler(E_MSG, routine, 'Perturbing model state', &
                       source, revision, revdate) 

    call get_grid_field_size(restart_grid, ijarea)
    call get_grid_dims(restart_grid, gridii, gridjj)
    call get_num_vars(num_vars)

    allocate( temp_state(gridii,gridjj), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'temp_state')
    allocate( temp_pert(gridii,gridjj) , stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'temp_pert')

    ! Each field in the restart vector may be treated differently
    do kk=1,num_vars
       call get_pert_type_by_index(kk,perttype)
       call get_pert_magnitude_by_index(kk, pertsize)
 
       ! Pull out just this field from the state vector      
       fullii = 1 + (kk - 1) * ijarea
       temp_state = reshape(state(fullii:(kk * ijarea)), &
                           (/ gridii, gridjj /))

       ! Apply the perturbation - may not perturb the boundaries
       if (perttype .eq. PERT_TYPE_NOPERTS) then

          ! Perturbed state is the same as the state
          do jj=1,gridjj
             do ii=1,gridii
                temp_pert(ii,jj) = temp_state(ii,jj)
             end do
          end do

       elseif (perttype .eq. PERT_TYPE_UNIFORM) then

          ! Single perturbation value for this field
          pert = random_gaussian(random_sequence, 0.0_r8, pertsize)

          do    jj=(1 + y_bound_skip), (gridjj - y_bound_skip)
             do ii=(1 + x_bound_skip), (gridii - x_bound_skip)
                temp_pert(ii,jj) = temp_state(ii,jj) + pert
             enddo
          end do

       elseif (perttype .eq. PERT_TYPE_INDIVID) then

          ! Each point in the field gets its own perturbation
          do    jj=(1 + y_bound_skip), (gridjj - y_bound_skip)
             do ii=(1 + x_bound_skip), (gridii - x_bound_skip)
                pert = random_gaussian(random_sequence, 0.0_r8, pertsize)
                temp_pert(ii,jj) = temp_state(ii,jj) + pert
             enddo
          end do

       end if

       pert_state(fullii:(kk*ijarea)) = reshape(temp_pert, (/ ijarea /))
     enddo

    interf_provided = .true.

    deallocate(temp_state, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'temp_state')
    deallocate(temp_pert,  stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'temp_pert')
  end subroutine pert_model_state

  ! model_interpolate
  ! -----------------
  ! Given the DART state vector, a location, and a raw variable type
  ! interpolates the state to that location
  ! This implementation currently only supports pressure coordinates.
  !  PARAMETERS
  !   IN  x                 full state vector
  !   IN  location          DART location structure for interpolation
  !                         target location
  !   IN  obs_type          raw variable type
  !   OUT obs_val           interpolated value
  !   OUT interp_status     status of interpolation (0 is success)
  subroutine model_interpolate(x, location, obs_type, obs_val, &
                               interp_status)
    real(r8), dimension(:), intent(in)  :: x
    type(location_type),    intent(in)  :: location
    integer,                intent(in)  :: obs_type
    real(r8),               intent(out) :: obs_val
    integer,                intent(out) :: interp_status

    logical :: successful_interpolation

    ! Just hand this off to the COAMPS interpolation module - it will
    ! handle getting all the proper variables out - note that we are
    ! currently not supporting sea level pressure
    call interpolate(x, restart_grid, location, obs_type, obs_val, &
                     successful_interpolation)
    if (successful_interpolation) then
      interp_status = 0
    else
      obs_val = MISSING_R8
      interp_status = 1
    end if
  end subroutine model_interpolate

  ! ens_mean_for_model
  ! ------------------
  ! Allow the ensemble mean to be passed in and stored if we need it
  ! (can be handy for forward operators)
  !  PARAMETERS
  !   IN  ens_mean          ensemble mean state vector
  subroutine ens_mean_for_model(ens_mean)
    real(r8), intent(in), dimension(:) :: ens_mean

    if (need_mean) then
       ensemble_mean = ens_mean
    end if
  end subroutine ens_mean_for_model

  ! get_close_maxdist_init
  ! ----------------------
  ! Set the maximum distance for the processor for finding nearby 
  ! points. Wrapper for location module's get_close_maxdist_init 
  ! subroutine.
  !  PARAMETERS
  ! INOUT gc                get_close_type structure to initialize
  !   IN  maxdist           the maximum distance to process  
  subroutine get_close_maxdist_init (gc, maxdist)
    type(get_close_type), intent(inout) :: gc
    real(r8), intent(in)                :: maxdist

    call loc_get_close_maxdist_init(gc, maxdist)
  end subroutine get_close_maxdist_init

  ! get_close_obs_init
  ! ------------------
  ! Initializes part of get_close accelerator that depends on the
  ! particular observation(s).  Wrapper for location module's
  ! get_close_obs_init subroutine.
  !  PARAMETERS
  ! INOUT  gc               get_close_type accelerator to initialize
  !   IN   num              number of observations in the set
  !   IN   obs              set of observation locations
  subroutine get_close_obs_init(gc, num, obs)
    type(get_close_type), intent(inout) :: gc
    integer, intent(in)                 :: num
    type(location_type), intent(in)     :: obs(num)
    
    call loc_get_close_obs_init(gc, num, obs)
  end subroutine get_close_obs_init

  ! get_close_obs
  ! -------------
  ! Gets the number of close observations.  Wrapper for location
  ! module's get_close_obs subroutine.
  !  PARAMETERS
  !   IN  gc                get_close_type accelerator
  !   IN  base_obs_loc      location of the base observation
  !   IN  obs_loc           location of all the observations
  !   IN  base_obs_kind     raw type of the base observation
  !   IN  obs_kind          raw type of all the observations
  !   OUT num_close         how many observations are close to the
  !                         base observation
  !   OUT close_ind         which of the observations are close to
  !                         the base observation
  !   OUT dist              OPTIONAL distance from the observations
  !                         to the base observation
  subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc,&
       & obs_kind, num_close, close_ind, dist)

    type(get_close_type), intent(in)    :: gc
    type(location_type), intent(inout)  :: base_obs_loc, obs_loc(:)
    integer, intent(in)                 :: base_obs_kind, obs_kind(:)
    integer, intent(out)                :: num_close, close_ind(:)
    real(r8), optional, intent(out)     :: dist(:)

    integer                          :: t_ind, istatus1, istatus2, k
    integer                          :: base_which, local_obs_which
	real(r8), dimension(3)           :: base_array, local_obs_array
	real(r8)                         :: obs_val
    type(location_type)              :: local_obs_loc

    ! Initialize variables to missing status
    num_close = 0 ; close_ind = -99 ; dist = 1.0e9
    istatus1  = 0 ; istatus2  = 0

    ! Convert base_obs vertical coordinate to requested vertical coordinate if necessary
    base_array = get_location(base_obs_loc)
    base_which = nint(query_location(base_obs_loc))

    if (.not. horiz_dist_only) then
      if (base_which /= VERTISLEVEL) then
        call model_interpolate(ensemble_mean, base_obs_loc, KIND_VERTLEVEL, obs_val, istatus1)
        base_array(3)=obs_val ; base_which=VERTISLEVEL
        base_obs_loc = set_location(base_array(1),base_array(2),base_array(3),base_which)
      elseif (base_array(3) == MISSING_R8) then
        istatus1 = 1
      endif
    endif

    if (istatus1 == 0) then
      ! Get all the potentially close obs but no dist (optional argument dist(:) is not present)
      ! This way, we are decreasing the number of distance computations that will follow.
      ! This is a horizontal-distance operation and we don't need to have the relevant vertical
      ! coordinate information yet (for obs_loc).
      call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                             num_close, close_ind)

      ! Loop over potentially close subset of obs priors or state variables
      do k = 1, num_close

        t_ind           = close_ind(k)
        local_obs_loc   = obs_loc(t_ind)
        local_obs_which = nint(query_location(local_obs_loc))
        local_obs_array = get_location(local_obs_loc)

        ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
        ! This should only be necessary for obs priors, as state location information already
        ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
        if (.not. horiz_dist_only) then
          if (local_obs_which /= VERTISLEVEL) then
            call model_interpolate(ensemble_mean, obs_loc(t_ind), KIND_VERTLEVEL, obs_val, istatus2)
            local_obs_array(3)=obs_val ; local_obs_which=VERTISLEVEL

            ! Store the "new" location into the original full local array
            local_obs_loc = set_location(local_obs_array(1),local_obs_array(2), &
                                         local_obs_array(3),local_obs_which)
            obs_loc(t_ind) = local_obs_loc
          endif
        endif

        ! Compute distance - set distance to a very large value if vert coordinate is missing
        ! or vert_interpolate returned error (istatus2=1)
        if (((.not. horiz_dist_only).and.(local_obs_array(3) == missing_r8)).or.(istatus2 == 1)) then
          dist(k) = 1.0e9
        else
          dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
        endif
      end do
    end if

  end subroutine get_close_obs

  ! get_model_size
  ! --------------
  ! Returns the size of the DART state vector
  !  PARAMETERS
  !   OUT get_model_size    length of the DART state vector
  function get_model_size()
    integer :: get_model_size

    integer :: level_gridpoints, total_levels

    call get_grid_field_size(restart_grid, level_gridpoints)
    call get_num_vars(total_levels)
    get_model_size = level_gridpoints * (total_levels)
  end function get_model_size

  ! get_state_meta_data
  ! -------------------
  ! Get the location and variable kind for a particular index in
  ! the state vector
  !  PARAMETERS
  !   IN  index_in          position in state vector to query
  !   OUT location          DART location_type for that index
  !   OUT var_kind          OPTIONAL numeric variable type
  subroutine get_state_meta_data(index_in, location, var_kind)
    integer,             intent(in)            :: index_in
    type(location_type), intent(out)           :: location
    integer,             intent(out), optional :: var_kind

    ! All this has been pre-calculated
    location = all_locs(index_in)
    if (present(var_kind)) then
       var_kind = all_kinds(index_in)
    end if
  end subroutine get_state_meta_data

  ! get_model_time_step
  ! -------------------
  ! Returns the smallest increment in time that the model is capable
  ! of advancing the state - just call it a minute for now
  !  PARAMETERS
  !   OUT get_model_time_step  model time step as a DART time_type
  function get_model_time_step()
    type(time_type) :: get_model_time_step

    get_model_time_step = set_time(60,0)
  end function get_model_time_step

  ! init_conditions
  ! ---------------
  ! NULL INTERFACE
  ! (sets up initial conditions for the model, but we're using
  ! already existing restart file data)
  !  PARAMETERS
  !   OUT x                 state vector initial condition 
  subroutine init_conditions(x)
    real(r8), intent(out) :: x(:)

    call error_handler(E_ERR, 'init_conditions', 'inoperable', &
                       source, revision, revdate) 
    x = MISSING_R8

  end subroutine init_conditions

  ! init_time
  ! ---------
  ! NULL INTERFACE
  ! (sets up initial time information for the model, but we're using
  ! already existing restart file data)
  !  PARAMETERS
  !   OUT time              time initial condition
  subroutine init_time(time)
    type(time_type), intent(out) :: time

    time = set_time_missing()

  end subroutine init_time

  ! adv_1step
  ! ---------
  ! NULL INTERFACE
  ! (advances the model with function calls, but we're doing it
  ! asynchronously)
  !  PARAMETERS
  ! INOUT x                 (in)  state vector analysis
  !                         (out) state vector forecast 
  ! INOUT time              (in)  analysis time
  !                         (out) forecast time 
  subroutine adv_1step(x, time)
    real(r8),        intent(inout) :: x(:)
    type(time_type), intent(in)    :: time

    if (do_output()) then
       call print_time(time,'NULL interface adv_1step (no advance) DART time is')
       call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
       call print_date(time,'NULL interface adv_1step (no advance) DART date is')
       call print_date(time,'NULL interface adv_1step (no advance) DART date is',logfileunit)
    endif
    call error_handler(E_ERR, 'adv_1step', 'cannot advance synchronously', &
                       source, revision, revdate) 
    x = MISSING_R8  ! Just to satisfy compiler/complainer

  end subroutine adv_1step

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! nc_write_statearray_atts
  ! -------------------
  ! Write model-specific global attributes to a NetCDF file.
  !  PARAMETERS
  !   IN  ncFileID          numeric ID of an *open* NetCDF file
  subroutine nc_write_statearray_atts( ncFileID )
    integer, intent(in)         :: ncFileID      ! netCDF file 

    integer                     :: model_size 
    type(location_type)         :: loc        
    real(kind=r8), dimension(3) :: loc3d      ! lon/lat/hgt
    integer                     :: ii         
 
    ! NetCDF variables                                               
    integer :: state_varid             ! State vars     
    integer :: state_coord_varid       ! State vector coordinate var 
    integer :: type_varid              ! State var type variable
    integer :: state_dimid             ! state KIND dimension
    integer ::   mem_dimid             ! Ensemble member/copy dimension   
    integer ::  time_dimid             ! Time/Unlimited dimension   
    integer :: lat_dimid,lat_varid     ! Latitude coordinate
    integer :: lon_dimid,lon_varid     ! Longitude coordinate
    integer :: vert_dimid,vert_varid   ! Vertical coordinate

    ! Error handling
    character(len=*), parameter :: routine = 'nc_write_statearray_atts'

    integer, allocatable        :: vartype(:)    ! state variable type
    real(kind=r8), allocatable  :: lon_var(:)    ! 
    real(kind=r8), allocatable  :: lat_var(:)    ! 
    real(kind=r8), allocatable  :: vrt_var(:)    ! 

    ! Find the ensemble member and time dimensions
    call nc_check(nf90_inq_dimid(ncFileID, 'copy',  mem_dimid), routine)
    call nc_check(nf90_inq_dimid(ncFileID, 'time', time_dimid), routine)

    ! State dimension                                     
    model_size = get_model_size() 
    call nc_check(nf90_def_dim(ncFileID, 'StateVariable', model_size, state_dimid), &
                  routine)

    ! Location dimensions
    call nc_check(nf90_def_dim(ncFileID, 'lat',    model_size,  lat_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID, 'lon',    model_size,  lon_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID, 'height', model_size, vert_dimid), routine)

    ! Location coordinates - latitude
    call nc_check(nf90_def_var(ncFileID, 'lat', nf90_double, lat_dimid, lat_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'long_name',       'Latitude'), routine)
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'units',      'degrees_north'), routine)
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'valid_range', &
                                         (/ -90.0_r8, 90.0_r8 /)), routine)

    ! Location coordinates - longitude
    call nc_check(nf90_def_var(ncFileID, 'lon', nf90_double, lon_dimid, lon_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'long_name',      'Longitude'), routine)
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'units',       'degrees_east'), routine) 
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'valid_range', &
                                         (/ 0.0_r8, 360.0_r8 /)), routine)

    ! Location coordinates - height
    call nc_check(nf90_def_var(ncFileID, 'height', nf90_double, vert_dimid, vert_varid), routine)
    call nc_check(nf90_put_att(ncFileID, vert_varid, 'long_name', 'Sigma Height'), routine)
    call nc_check(nf90_put_att(ncFileID, vert_varid, 'units', 'meters'), routine)
 
    ! Location coordinates - state    
    call nc_check(nf90_def_var(ncFileID, 'StateVariable', nf90_int, &
                                         state_dimid, state_coord_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, state_coord_varid, 'long_name',     &
                                         'State Variable Index'), routine) 
    call nc_check(nf90_put_att(ncFileID, state_coord_varid, 'units', 'index'), routine) 
    call nc_check(nf90_put_att(ncFileID, state_coord_varid, 'valid_range', &
                                         (/1, model_size /)), routine) 

    ! State variable              
    call nc_check(nf90_def_var(ncFileID, 'state', nf90_double,           &
                               (/ state_dimid, mem_dimid, time_dimid /), &
                                 state_varid), routine)  
    call nc_check(nf90_put_att(ncFileID, state_varid, 'long_name', 'model state'), routine) 

    ! Variable types
    call nc_check(nf90_def_var(ncFileID, 'type', nf90_int, state_dimid, type_varid), routine)
    call nc_check(nf90_put_att(ncFileID, type_varid, 'long_name', 'Variable Type'), routine)

    !-------------------------------------------------------------------
    ! Need to get out of define mode to fill the variables 
    !-------------------------------------------------------------------
    call nc_check(nf90_enddef(ncFileID), routine) 
 
    ! Fill the state variable coordinate                             
    call nc_check(nf90_put_var(ncFileID, state_coord_varid, &
                               (/ (ii,ii=1, model_size) /)), routine)  
  
    allocate(lon_var(model_size))
    allocate(lat_var(model_size))
    allocate(vrt_var(model_size))
    allocate(vartype(model_size))

    do ii=1, model_size
       call get_state_meta_data(ii,loc,vartype(ii))
       loc3d = get_location(loc)
       lon_var(ii) = loc3d(1)
       lat_var(ii) = loc3d(2)
       vrt_var(ii) = loc3d(3)
    end do 

    call nc_check(nf90_put_var(ncFileID,  lon_varid, lon_var), routine)
    call nc_check(nf90_put_var(ncFileID,  lat_varid, lat_var), routine)
    call nc_check(nf90_put_var(ncFileID, vert_varid, vrt_var), routine)
    call nc_check(nf90_put_var(ncFileID, type_varid, vartype), routine)

    deallocate(lon_var)
    deallocate(lat_var)
    deallocate(vrt_var)
    deallocate(vartype)
 
    ! Sync the NetCDF file and buffer, but don't close it                                                       
    call nc_check(nf90_sync(ncFileID), routine) 
    if (do_output()) write (*,*) 'Model attributes written and file synced...' 

  end subroutine nc_write_statearray_atts

  ! nc_write_prognostic_atts
  ! -------------------
  ! Write model-specific global attributes to a NetCDF file.
  !  PARAMETERS
  !   IN  ncFileID          numeric ID of an *open* NetCDF file
  !   OUT ierr              0 if writing was successful
  !
  ! If the dimension of the domain is specified to by nx, ny, nx ...
  ! then var_stagger is related to var_dims in the following way: 
  !
  ! T-stagger has dimensions (nx,ny,nz); 
  ! W-stagger has dimensions (nx,ny,nz+1); 
  ! U-stagger has dimensions (nx-1,ny,nz); and 
  ! V-stagger has dimensions (nx,ny-1,nz) 
  ! (COAMPS is on a C-grid with p-u-p ordering in the horizontal).
 
  subroutine nc_write_prognostic_atts( ncFileID )
    integer, intent(in)         :: ncFileID      ! netCDF file 

    integer                     :: model_size, nx, ny, nz
    type(location_type)         :: loc        
    real(kind=r8), dimension(3) :: loc3d      ! lon/lat/hgt
    integer                     :: ii 
 
    ! NetCDF variables                                               
    integer :: time_dimid              ! Time/Unlimited dimension
    integer ::  mem_dimid              ! Ensemble member/copy dimension   
    integer ::  lat_dimid,  lat_varid  ! Latitude  coordinate
    integer ::  lon_dimid,  lon_varid  ! Longitude coordinate
    integer ::  lev_dimid,  lev_varid  ! Vertical  coordinate
    integer :: slat_dimid, slat_varid  ! Staggered Latitude  coordinate
    integer :: slon_dimid, slon_varid  ! Staggered Longitude coordinate
    integer :: slev_dimid, slev_varid  ! Staggered Vertical  coordinate

    ! Error handling
    character(len=*), parameter :: routine = 'nc_write_prognostic_atts'

    type(state_vec_iterator)     :: iterator
    integer                      :: ndims, varid
    integer, dimension(3)        :: var_dims, dimids
    integer, dimension(6,2)      :: localdimids
    character(len=1)             :: var_stagger ! (U|V|W|T)
    character(len=NF90_MAX_NAME) :: var_name
    character(len=NF90_MAX_NAME), dimension(10) :: variable_list

    real(kind=r8), allocatable ::  lon_var(:)
    real(kind=r8), allocatable ::  lat_var(:)
    real(kind=r8), allocatable ::  lev_var(:)
    real(kind=r8), allocatable :: slon_var(:)
    real(kind=r8), allocatable :: slat_var(:)
    real(kind=r8), allocatable :: slev_var(:)

    call get_grid_dims(restart_grid, nx, ny)
    call get_grid_num_levels(restart_grid, nz) ! mass sigma levels

    ! Find the ensemble member and time dimensions
    call nc_check(nf90_inq_dimid(ncFileID, 'copy',  mem_dimid), routine)
    call nc_check(nf90_inq_dimid(ncFileID, 'time', time_dimid), routine)
 
    ! Location dimensions
    call nc_check(nf90_def_dim(ncFileID,  'lat', ny,    lat_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID,  'lon', nx,    lon_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID,  'lev', nz,    lev_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID, 'slat', ny-1, slat_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID, 'slon', nx-1, slon_dimid), routine)
    call nc_check(nf90_def_dim(ncFileID, 'slev', nz+1, slev_dimid), routine)

    ! create local copy of dimension info to speed lookup

    localdimids(1,:) = (/ ny  ,  lat_dimid /)
    localdimids(2,:) = (/ nx  ,  lon_dimid /)
    localdimids(3,:) = (/ nz  ,  lev_dimid /)
    localdimids(4,:) = (/ ny-1, slat_dimid /)
    localdimids(5,:) = (/ nx-1, slon_dimid /)
    localdimids(6,:) = (/ nz+1, slev_dimid /)

    ! Location coordinates - latitude
    call nc_check(nf90_def_var(ncFileID, 'lat', nf90_double, lat_dimid, lat_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'long_name',       'Latitude'), routine)
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'units',      'degrees_north'), routine)
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'valid_range',  &
                               (/ -90.0_r8, 90.0_r8 /)), routine)

    ! Location coordinates - longitude
    call nc_check(nf90_def_var(ncFileID, 'lon', nf90_double, lon_dimid, lon_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'long_name', 'Longitude'), routine)
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'units', 'degrees_east'), routine) 
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'valid_range', &
                                (/ 0.0_r8, 360.0_r8 /)), routine)

    ! Location coordinates - height
    call nc_check(nf90_def_var(ncFileID, 'height', nf90_double, lev_dimid, lev_varid), routine)
    call nc_check(nf90_put_att(ncFileID, lev_varid, 'long_name', 'Sigma Height'), routine)
    call nc_check(nf90_put_att(ncFileID, lev_varid, 'units', 'meters'), routine)

    ! Location coordinates - staggered latitude
    call nc_check(nf90_def_var(ncFileID, 'slat', nf90_double, slat_dimid, slat_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, slat_varid, 'long_name', 'Staggered Latitude'), routine)
    call nc_check(nf90_put_att(ncFileID, slat_varid, 'units',      'degrees_north'), routine)
    call nc_check(nf90_put_att(ncFileID, slat_varid, 'valid_range',  &
                               (/ -90.0_r8, 90.0_r8 /)), routine)

    ! Location coordinates - staggered longitude
    call nc_check(nf90_def_var(ncFileID, 'slon', nf90_double, slon_dimid, slon_varid), routine) 
    call nc_check(nf90_put_att(ncFileID, slon_varid, 'long_name', 'Staggered Longitude'), routine)
    call nc_check(nf90_put_att(ncFileID, slon_varid, 'units', 'degrees_east'), routine) 
    call nc_check(nf90_put_att(ncFileID, slon_varid, 'valid_range', &
                                (/ 0.0_r8, 360.0_r8 /)), routine)

    ! Location coordinates - staggered height
    call nc_check(nf90_def_var(ncFileID, 'sheight', nf90_double, slev_dimid, slev_varid), routine)
    call nc_check(nf90_put_att(ncFileID, slev_varid, 'long_name', 'Staggered Sigma Height'), routine)
    call nc_check(nf90_put_att(ncFileID, slev_varid, 'units', 'meters'), routine)

    ! ------------------------------------------------------------------
    ! Loop over all the variables in the state vector
    ! ------------------------------------------------------------------

    iterator = get_iterator(variable_list)          ! FIXME
    define_atts:  do while(has_next(iterator))

      call get_var_dims(var_dims, ndims, iterator)  ! FIXME
      call get_var_name(var_name,        iterator)  ! FIXME
      call get_var_stagger(var_stagger,  iterator)  ! FIXME

      dimids = get_dimids( ndims, var_dims, localdimids ) 

      call nc_check(nf90_def_var(ncFileID, trim(var_name), nf90_double, &
                (/ dimids(1:ndims), mem_dimid, time_dimid /), varid),  &
                                   routine, 'def '//trim(var_name) ) 

      call nc_check(nf90_put_att(ncFileID, varid, 'grid', var_stagger//' type'), &
                                         routine, 'stagger for '//trim(var_name) )

      call get_next(iterator)
    end do define_atts

    ! ------------------------------------------------------------------
    ! Need to get out of define mode to fill the coordinate variables 
    ! FIXME: This is completely wrong for the prognostic variables.
    ! I don't know how to stride through get_state_meta_data to fill the
    ! appropriate arrays.
    ! ------------------------------------------------------------------
    call nc_check(nf90_enddef(ncFileID), routine) 
 
    allocate( lon_var(nx  ))
    allocate( lat_var(ny  ))
    allocate( lev_var(nz  ))
    allocate(slon_var(nx-1))
    allocate(slat_var(ny-1))
    allocate(slev_var(nz+1))

    do ii=1, model_size
       call get_state_meta_data(ii,loc)
       loc3d = get_location(loc)
        lon_var(ii) = loc3d(1)
        lat_var(ii) = loc3d(2)
        lev_var(ii) = loc3d(3)
       slon_var(ii) = loc3d(1)
       slat_var(ii) = loc3d(2)
       slev_var(ii) = loc3d(3)
    end do 

    call nc_check(nf90_put_var(ncFileID,  lon_varid,  lon_var), routine)
    call nc_check(nf90_put_var(ncFileID,  lat_varid,  lat_var), routine)
    call nc_check(nf90_put_var(ncFileID,  lev_varid,  lev_var), routine)
    call nc_check(nf90_put_var(ncFileID, slon_varid, slon_var), routine)
    call nc_check(nf90_put_var(ncFileID, slat_varid, slat_var), routine)
    call nc_check(nf90_put_var(ncFileID, slev_varid, slev_var), routine)

    deallocate( lon_var,  lat_var,  lev_var)
    deallocate(slon_var, slat_var, slev_var)
 
    ! Sync the NetCDF file, but don't close it                                                       
    call nc_check(nf90_sync(ncFileID), routine) 
    if (do_output()) write (*,*) 'Model attributes written and file synced...' 

  end subroutine nc_write_prognostic_atts

  ! get_iterator
  ! -------------------
  ! FIXME -  no idea what you do here
  !  PARAMETERS
  !   OUT var_name          string containing variable name
  !   IN  iterator          iterator for variable list
  function get_iterator(variable_list)
    character(len=*), dimension(:), intent(in) :: variable_list
    type(state_vec_iterator)                   :: get_iterator

    get_iterator%something(1) = variable_list(1) ! FIXME just to satisfy compiler/complainer

  end function get_iterator

  ! has_next
  ! -------------------
  ! FIXME - something like this, I suppose
  !  PARAMETERS
  !   OUT has_next          logical ... there is more to do
  !   IN  iterator          iterator for variable list
  function has_next(iterator)
    type(state_vec_iterator), intent(in) :: iterator
    logical                              :: has_next

    has_next = iterator%more ! FIXME just to satisfy compiler/complainer

  end function has_next

  ! get_next
  ! -------------------
  ! FIXME - completely in the dark here ...
  !  PARAMETERS
  !   INOUT  iterator          iterator for variable list
  subroutine get_next(iterator)
    type(state_vec_iterator), intent(inout) :: iterator

    if ( iterator%nitems < iterator%cur_item ) then
       iterator%cur_item = iterator%cur_item + 1
       iterator%more     = .true.
    else
       iterator%more     = .false.
    endif

  end subroutine get_next


  ! get_var_name
  ! -------------------
  ! Given an iterator to a variable list returns the current variable name
  !  PARAMETERS
  !   OUT var_name          string containing variable name
  !   IN  iterator          iterator for variable list
  subroutine get_var_name(var_name, iterator)
    character(len=NF90_MAX_NAME), intent(out) :: var_name
    type(state_vec_iterator),     intent(in)  :: iterator

    var_name = iterator%something( iterator%cur_item )
  end subroutine get_var_name

  ! get_var_dims
  ! -------------------
  ! Given an iterator to a variable list returns the an array containing
  ! the dimensions and the number of non-singleton dimensions
  !  PARAMETERS
  !   OUT var_dims(3)       3D array containing spatial dimensions   
  !   OUT ndims             number of non-singleton dimensions
  !   IN  iterator          iterator for variable list
  subroutine get_var_dims(var_dims, ndims, iterator)
    integer, dimension(3),    intent(out) :: var_dims
    integer,                  intent(out) :: ndims
    type(state_vec_iterator), intent(in)  :: iterator

    write(*,*) iterator%something(1) ! FIXME just to satisfy compiler/complainer
    ndims    = 3                     ! FIXME just to satisfy compiler/complainer
    var_dims = (/ 20, 30, 40 /)      ! FIXME just to satisfy compiler/complainer
  end subroutine

  ! get_var_stagger
  ! -------------------
  ! Given an iterator to a variable list returns the stagering of that
  ! variable
  !  PARAMETERS
  !   OUT var_stagger       Staggering of field (U|V|W|T)
  !   IN  iterator          iterator for variable list
  subroutine get_var_stagger(var_stagger,iterator)
    character(len=1),         intent(out) :: var_stagger
    type(state_vec_iterator), intent(in)  :: iterator

    write(*,*) iterator%something(1) ! FIXME just to satisfy compiler/complainer
    var_stagger = 'T'                ! FIXME just to satisfy compiler/complainer
  end subroutine get_var_stagger

  ! get_var_3d
  ! -------------------
  ! Given an iterator to a variable list and the big state vector
  ! returns the variable
  !  PARAMETERS
  !   OUT var_out           The output variable
  !   IN  statevec          The big state vector
  !   IN  iterator          iterator for variable list
  subroutine get_var_3d(statevec,var_out,iterator)
    real(r8), dimension(:,:,:), intent(out) :: var_out
    real(r8), dimension(:),     intent(in)  :: statevec
    type(state_vec_iterator),   intent(in)  :: iterator

    write(*,*) iterator%something  ! FIXME just to satisfy compiler/complainer
    write(*,*) statevec(1)         ! FIXME just to satisfy compiler/complainer
    var_out = 0.0_r8               ! FIXME just to satisfy the compiler/complainer
  end subroutine get_var_3d

  ! get_var_2d
  ! -------------------
  ! Given an iterator to a variable list and the big state vector
  ! returns the variable
  !  PARAMETERS
  !   OUT var_out           The output variable
  !   IN  statevec          The big state vector
  !   IN  iterator          iterator for variable list
  subroutine get_var_2d(statevec,var_out,iterator)
    real(r8), dimension(:,:),   intent(out) :: var_out
    real(r8), dimension(:),     intent(in)  :: statevec
    type(state_vec_iterator),   intent(in)  :: iterator

    write(*,*) iterator%something  ! FIXME just to satisfy compiler/complainer
    write(*,*) statevec(1)         ! FIXME just to satisfy compiler/complainer
    var_out = 0.0_r8               ! FIXME just to satisfy the compiler/complainer
  end subroutine get_var_2d

  ! get_var_1d
  ! -------------------
  ! Given an iterator to a variable list and the big state vector
  ! returns the variable
  !  PARAMETERS
  !   OUT var_out           The output variable
  !   IN  statevec          The big state vector
  !   IN  iterator          iterator for variable list
  subroutine get_var_1d(statevec,var_out,iterator)
    real(r8), dimension(:),     intent(out) :: var_out
    real(r8), dimension(:),     intent(in)  :: statevec
    type(state_vec_iterator),   intent(in)  :: iterator

    write(*,*) iterator%something  ! FIXME just to satisfy compiler/complainer
    write(*,*) statevec(1)         ! FIXME just to satisfy compiler/complainer
    var_out = 0.0_r8               ! FIXME just to satisfy the compiler/complainer
  end subroutine get_var_1d



  ! get_dimids
  ! -------------------
  ! Given a list of variable dimensions, try to figure out
  ! the matching netCDF dimension IDs (the table of netCDF
  ! dimension IDs and sizes is stored in a local variable
  ! to avoid incessant netCDF queries).
  !  PARAMETERS
  !   IN  ndims           the rank of the variable
  !   IN  var_dims        the array of the dimension of the variable
  !   IN  localdimds      the table of netCDF dim IDs and sizes
  !   OUT mydims          the array of netCDF dimension IDs 
  function get_dimids( ndims, var_dims, localdimids ) result(mydims)
    integer,                 intent(in) :: ndims
    integer, dimension(:),   intent(in) :: var_dims
    integer, dimension(:,:), intent(in) :: localdimids
    integer, dimension(SIZE(var_dims))  :: mydims

    integer :: i,j

    mydims = -999  ! set array to a bad value

    do i=1,ndims
       dimloop: do j=1,SIZE(localdimids,1)
          if ( var_dims(i) == localdimids(j,1) ) then
                 mydims(i) =  localdimids(j,2)
             exit dimloop
          endif
       enddo dimloop
    enddo

    if ( any(mydims(1:ndims) < 0 )) then
       write(msgstring,*) 'ERROR : cannot find dimids to match ',mydims(1:ndims)
       call error_handler(E_ERR, 'get_dimids', msgstring, &
                       source, revision, revdate) 
    endif
  end function get_dimids

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
