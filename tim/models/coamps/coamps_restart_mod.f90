! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module coamps_restart_mod

!------------------------------
! MODULE:       coamps_restart_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module containing the data structure and routines for dealing with
! COAMPS restart files and the dynamic DART state vector definition
!------------------------------ 

  use coamps_grid_mod, only : coamps_grid,         &
                              get_grid_dims,       &
                              get_grid_field_size, &
                              get_grid_num_levels, &
                              initialize_grid
  use coamps_util_mod, only : check_alloc_status,  &
                              uppercase,lowercase, & 
                              check_io_status 

  use location_mod,    only : VERTISLEVEL,   &
                              VERTISSURFACE
  use obs_kind_mod
  use types_mod,       only : r8
  use utilities_mod,   only : do_output,     &
                              E_ERR,         &
                              E_MSG,         &
                              E_WARN,        &
                              error_handler, &
                              get_unit

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  ! Initialization
  public :: initialize_state_info

  ! Size information
  public :: get_file_type
  public :: get_num_lvls
  public :: get_num_flds
  public :: get_num_vars

  ! Grid information
  public :: get_restart_grid

  ! Search functions
  public :: get_restart_index_by_properties

  ! Accessors by index
  public :: get_full_record_num_by_index
  public :: get_var_name_by_index
  public :: get_sigma_by_index
  public :: get_var_type_by_index
  public :: get_mean_flag_by_index
  public :: get_posdef_flag_by_index
  public :: get_vert_coord_by_index
  public :: get_pert_magnitude_by_index
  public :: get_pert_type_by_index
  public :: get_kind_by_index
  public :: get_update_flag_by_index
  public :: get_var_info_by_abs_index

  ! Numeric constants
  public :: PERT_TYPE_NOPERTS
  public :: PERT_TYPE_UNIFORM
  public :: PERT_TYPE_INDIVID
  public :: VAR_NAME_LEN

  ! Diagnostics
  public :: dump_restart_vars
  
  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACE
  !------------------------------

  ! This is the external routine that grabs where fields are in the
  ! restart file - it's generated via shell scripts and is therefore
  ! included as a separate procedure.  Pass in the constants so we
  ! don't need to repeat them inside the subroutine.
  interface
    subroutine get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,&
                             SINGLEIO, MULTIIO, var_name,           &
                             var_dim_type, var_record_num)
      integer,               intent(in)  :: DIM_TYPE_2D
      integer,               intent(in)  :: DIM_TYPE_3D
      integer,               intent(in)  :: DIM_TYPE_3DW
      integer,               intent(in)  :: SINGLEIO
      integer,               intent(in)  :: MULTIIO
      character(len=*),      intent(in)  :: var_name
      integer,               intent(out) :: var_dim_type
      integer, dimension(2), intent(out) :: var_record_num
    end subroutine get_name_info
  end interface

  !------------------------------
  ! END EXTERNAL INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS
  !------------------------------

  ! Dimension information - pretty straightforward except that 3Dw
  ! includes the surface, so there are kka+1 levels there
  integer, parameter :: DIM_TYPE_2D  = 1
  integer, parameter :: DIM_TYPE_3D  = 2
  integer, parameter :: DIM_TYPE_3DW = 3

  ! Length of some string variables
  integer, parameter :: VAR_NAME_LEN  = 10
  integer, parameter :: VAR_TYPE_LEN  = 32
  integer, parameter :: PERT_TYPE_LEN = 7

  ! Perturbation types - don't perturb at all, perturb whole sigma
  ! level field uniformly, or perturb each grid point individually.
  integer, parameter :: PERT_TYPE_NOPERTS = 0
  integer, parameter :: PERT_TYPE_UNIFORM = 1
  integer, parameter :: PERT_TYPE_INDIVID = 2
  character(len=PERT_TYPE_LEN), parameter :: NOPERTS_NAME = 'NOPERTS' 
  character(len=PERT_TYPE_LEN), parameter :: UNIFORM_NAME = 'UNIFORM'
  character(len=PERT_TYPE_LEN), parameter :: INDIVID_NAME = 'INDIVID'

  ! If the variable is defined on a mass level or a w level
  character, parameter :: MASS_LEVEL = 'M'
  character, parameter :: W_LEVEL    = 'W'

  ! Indices for the variable record array to take into account the
  ! difference between the numbering in the single I/O case and
  ! multiple process I/O case
  integer, parameter :: SINGLEIO = 1
  integer, parameter :: MULTIIO  = 2

  ! Whether or not to write the field back to the COAMPS restart file
  logical, parameter :: FLAG_UPDATE_FIELD  = .true.
  logical, parameter :: FLAG_FREEZE_FIELD  = .false.

  ! Define the strings used to signify whether we'll do the updates
  ! or not
  integer, parameter :: FIELD_UPDATE_LEN = 6
  character(len=FIELD_UPDATE_LEN), parameter :: FIELD_UPDATE='UPDATE'
  character(len=FIELD_UPDATE_LEN), parameter :: FIELD_FREEZE='FREEZE'

  ! Define the strings used to signify whether we'll check for
  ! positive definiteness or not and their corresponding flags
  integer, parameter :: POS_DEF_LEN = 8
  character(len=POS_DEF_LEN), parameter :: FORCE_POSITIVE='ISPOSDEF'
  character(len=POS_DEF_LEN), parameter :: ALLOW_NEGATIVE='NOPOSDEF'
  logical, parameter :: FLAG_FORCE_POSITIVE = .true.
  logical, parameter :: FLAG_ALLOW_NEGATIVE = .false.

  ! Define the strings used to signify whether flat or restart
  ! files are used for i/o
  integer, parameter :: FILE_TYPE_LEN  = 12
  character(len=FILE_TYPE_LEN), parameter :: FLAT_FILE_STRING='FLAT_FILE'
  character(len=FILE_TYPE_LEN), parameter :: RESTART_FILE_STRING='RESTART_FILE'
  logical, parameter :: FLAG_FLAT = .true.
  logical, parameter :: FLAG_RESTART = .false.


  ! Define an entry for a single field in the DART state vector.
  ! This encompasses several pieces of data
  !  Where to find the field in the COAMPS restart file:
  !   var_name, dim_type, var_record, sigma_record
  ! How to perturb the variable if necessary:
  !   pert_mag, pert_type
  ! Whether or not we actually assimilate the field (e.g. mean
  ! fields are needed in here to do interpolation but should
  ! not have their assimilated alterations written back to the
  ! COAMPS restart file):  
  !   update_field
  ! Help find specific variables for interpolation:
  !   mean_field, mass_level
  ! Forbid negative values and set anything less than zero to
  ! zero (e.g. for mixing ratios):
  !   positive_definite
  type :: restart_var
     character(len=VAR_NAME_LEN) :: var_name          ! variable name
     integer                     :: dim_type          ! Restart: 2D or 3D. Flat: Not used.
     integer,dimension(2)        :: var_record        ! Restart: Start/Stop location in file. Flat: Not used.
     integer                     :: sigma_record      ! vertical level index
     real(kind=r8)               :: pert_mag          ! magnitude of initial perturbation
     integer                     :: pert_type         ! perturbation type
     logical                     :: update_field      ! should this field be updated?
     integer                     :: var_type          ! variable KIND       
     logical                     :: mean_field        ! is this a mean field?
     logical                     :: mass_level        ! is this a mass level?
     logical                     :: positive_definite ! is the field positive definate?
     integer                     :: num_lvls          ! number of levels in the field.
!     character(len=1)            :: staggering        ! staggering (U|V|W|T)          
  end type restart_var

  ! Derived data type that contains the information used to read and
  ! store the field information for the state vector
  type :: coamps_restart
     private
     character(len=80)                            :: dsnrff
     character(len=10)                            :: cdtgm1
     character(len=10)                            :: cdtg
     type(coamps_grid)                            :: grid
     type(restart_var), dimension(:), allocatable :: restart_vars
     integer                                      :: num_flds
     logical                                      :: flat_file_io
  end type coamps_restart

  ! Iterator Types
  type :: state_vec_list
     private
     type(state_vec_node), pointer :: head
     type(state_vec_node), pointer :: tail
     integer                       :: list_length
  end type state_vec_list

  type :: state_vec_iterator
     private
     type(state_vec_node), pointer :: current
     integer                       :: cur_index
  end type state_vec_iterator

  type :: state_vec_node
     private
     type(restart_var)             :: value
     type(state_vec_node), pointer :: next
  end type state_vec_node
   ! End Iteratory Types

  !------------------------------
  ! END TYPES AND CONSTANTS
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------

  type(coamps_restart), target :: restart_info
  type(state_vec_list)         :: state_list
  type(state_vec_list)         :: var_list

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------

  ! initialize_state_info
  ! -----------------------
  ! Populate the restart_info structure by reading in the grid data
  ! and reading the state vector definition file
  !  PARAMETERS
  !   IN  dtg                  base date-time group for this run
  !   IN  restart_var_filename name of state vector definition file
  subroutine initialize_state_info(dtg,state_var_filename, dsnrff)
    character(len=10), intent(in) :: dtg
    character(len=*), intent(in)  :: state_var_filename
    character(len=*), intent(in)  :: dsnrff

    integer           :: state_var_unit

    character(len=*), parameter :: routine = 'initialize_state_info'
    integer           :: io_status

    restart_info%dsnrff   = dsnrff
    restart_info%cdtg     = dtg
    restart_info%num_flds = 0

    call initialize_grid(restart_info%cdtg, restart_info%grid, restart_info%dsnrff)

    state_var_unit = get_unit()
    open(unit=state_var_unit, file=state_var_filename, &
         status='old', access='sequential', action='read', &
         form='formatted', iostat=io_status)
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Opening file ' // state_var_filename)

    call read_state_var_file(state_var_unit)

    close(state_var_unit)

  end subroutine initialize_state_info

  ! get_file_type
  ! ------------
  ! Returns logical specifying restart or flat files
  !  PARAMETERS
  !  OUT file_type  Flat files or restart files?
  subroutine get_file_type(file_type) 
    logical, intent(out) :: file_type
    file_type = restart_info%flat_file_io
  end subroutine get_file_type

  ! get_num_lvls
  ! ------------
  ! Returns how many levels there are for this field
  !  PARAMETERS
  !   OUT num_lvls          number of lvls for this field
  subroutine get_num_lvls(nlvls,var_index)
    integer, intent(out) :: nlvls
    integer, intent(in)  :: var_index
    nlvls = restart_info%restart_vars(var_index)%num_lvls
  end subroutine get_num_lvls

  ! get_num_flds
  ! ------------
  ! Returns how many fields there are in the state vector
  !  PARAMETERS
  !   OUT num_vars          number of fields in state vector
  subroutine get_num_flds(num_flds)
    integer, intent(out) :: num_flds
    num_flds = restart_info%num_flds
  end subroutine get_num_flds

  ! get_num_vars
  ! ------------
  ! Returns how many levels there are in the state vector
  !  PARAMETERS
  !   OUT num_vars          number of levels in state vector
  subroutine get_num_vars(num_vars)
    integer, intent(out) :: num_vars
    num_vars = size(restart_info%restart_vars)
  end subroutine get_num_vars

  ! get_restart_grid
  ! ----------------
  ! Returns the grid structure associated with the state vector
  ! definition
  !  PARAMETERS
  !   OUT restart_grid      coamps_grid this state vector holds
  subroutine get_restart_grid(restart_grid)
    type(coamps_grid), intent(out) :: restart_grid
    restart_grid = restart_info%grid
  end subroutine get_restart_grid

  ! ----------------------------
  ! get_restart_index_by_properties
  ! -------------------------------
  ! Given the type of a variable, whether it's a mean field, and the
  ! type of level it is defined on plus the corresponding sigma index
  ! returns the position of that variable in the restart_vars array.
  !  PARAMETERS
  !   IN  var_type          integer form of variable type
  !   IN  is_mean           true if variable is a mean field
  !   IN  on_mass_level     true if variable is on a mass level
  !   IN  sig_index         sigma level index for variable
  !   OUT restart_index     position in the restart_vars array
  subroutine get_restart_index_by_properties(var_type, is_mean, &
                                             on_mass_level,     &
                                             sig_index,         &
                                             restart_index)
    integer, intent(in)  :: var_type
    logical, intent(in)  :: is_mean
    logical, intent(in)  :: on_mass_level
    integer, intent(in)  :: sig_index
    integer, intent(out) :: restart_index

    integer :: cur_index, total_vars
    integer :: cur_var_type
    integer :: cur_sigma
    logical :: cur_mean_flag
    logical :: cur_mass_flag

    character(len=128) :: message

    restart_index = -1

    ! Assume there is only one - if there are multiples, return the
    ! first one found
    call get_num_vars(total_vars)
    do cur_index = 1, total_vars
       cur_var_type  = restart_info%restart_vars(cur_index)%var_type
       cur_sigma     = restart_info%restart_vars(cur_index)%sigma_record
       cur_mean_flag = restart_info%restart_vars(cur_index)%mean_field
       cur_mass_flag = restart_info%restart_vars(cur_index)%mass_level

       if (cur_var_type .eq. var_type) then
          if (cur_sigma .eq. sig_index) then
             if (cur_mean_flag .eqv. is_mean) then
                if (cur_mass_flag .eqv. on_mass_level) then
                    restart_index = cur_index
                    exit
                end if
             end if
          end if
       end if
    end do

    ! The proper error handling in this case is simply to allow the
    ! negative restart index to stick around - let the caller handle
    ! dealing with this.  Since some variables might not be available
    ! at all levels, quitting would be a bad idea.
    if (restart_index .le. 0) then 
      write (message,*) "Could not find variable type ", var_type, &
                        "on (mass?",on_mass_level,") level",       &
                        sig_index
       call error_handler(E_WARN, "get_restart_index_by_properties",&
                         message, source, revision, revdate) 
    end if
  end subroutine get_restart_index_by_properties

  ! get_full_record_num_by_index
  ! ----------------------------
  ! Call get_full_record_num using a variable at the specified index
  ! in the state vector definition
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   IN  use_singleio          true if restart file is written using
  !                         a single I/O processor
  !   OUT recnum            variable's record number within the
  !                         restart file
  subroutine get_full_record_num_by_index(var_index, use_singleio, &
                                          recnum)
    integer, intent(in)  :: var_index
    logical, intent(in)  :: use_singleio
    integer, intent(out) :: recnum

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_full_record_num_by_index', &
                          'Variable index out of range', source, &
                          revision, revdate)
    else
       call get_full_record_num(restart_info%restart_vars(var_index),&
                                use_singleio, recnum)
    end if
  end subroutine get_full_record_num_by_index

  ! get_var_name_by_index
  ! ------------------
  ! Returns the variable name of a field at the specified index in the
  ! state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT var_name          variable's name
  subroutine get_var_name_by_index(var_index, var_name)
    integer, intent(in)                      :: var_index
    character(len=VAR_NAME_LEN), intent(out) :: var_name

    var_name = restart_info%restart_vars(var_index)%var_name
  end subroutine get_var_name_by_index

  ! get_sigma_by_index
  ! ------------------
  ! Returns the sigma level of a field at the specified index in the
  ! state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT sigma_index       variable's sigma level
  subroutine get_sigma_by_index(var_index, sigma_index)
    integer, intent(in)  :: var_index
    integer, intent(out) :: sigma_index

    sigma_index = restart_info%restart_vars(var_index)%sigma_record
  end subroutine get_sigma_by_index

  ! get_var_type_by_index
  ! ---------------------
  ! Return the variable type for the variable at the given index in
  ! the state vector definition
  !  PARAMETERS
  !   IN  restart_index     index of variable to interrogate
  !   OUT var_type          type of given variable
  subroutine get_var_type_by_index(restart_index, var_type)
    integer, intent(in)  :: restart_index
    integer, intent(out) :: var_type

    var_type = restart_info%restart_vars(restart_index)%var_type
  end subroutine get_var_type_by_index

  ! get_mean_flag_by_index
  ! ----------------------
  ! Returns if the field at the specified index in the state vector
  ! is a mean field
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT mean_field        true if the field is a mean field
  subroutine get_mean_flag_by_index(var_index, mean_field)
    integer, intent(in)  :: var_index
    logical, intent(out) :: mean_field

    mean_field = restart_info%restart_vars(var_index)%mean_field
  end subroutine get_mean_flag_by_index

  ! get_posdef_flag_by_index
  ! ----------------------------
  ! Call get_posdef_flag usign a variable at the specified index in
  ! the state vector definition
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   OUT update            true if this variable gets updated
  subroutine get_posdef_flag_by_index(var_index, posdef)
    integer, intent(in)  :: var_index
    logical, intent(out) :: posdef

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_posdef_flag_by_index',      &
                          'Variable index out of range', source,  &
                          revision, revdate)
    else
       call get_posdef_flag(restart_info%restart_vars(var_index), &
                            posdef)
    end if
  end subroutine get_posdef_flag_by_index

  ! get_vert_coord_by_index
  ! -----------------------
  ! Call get_vert_coord using a variable at the specified index in 
  ! the state vector definition
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   OUT vert_type         variable's vertical coordinate type
  !   OUT vert_level        variable's vertical coordinate
  subroutine get_vert_coord_by_index(var_index,vert_type,vert_level)
    integer, intent(in)        :: var_index
    integer, intent(out)       :: vert_type
    real(kind=r8), intent(out) :: vert_level

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_vert_coord_by_index',      &
                          'Variable index out of range', source, &
                          revision, revdate)
    else
       call get_vert_coord(restart_info%restart_vars(var_index), & 
                           restart_info%grid, vert_type, vert_level)
    end if
  end subroutine get_vert_coord_by_index

  ! get_pert_magnitude_by_index
  ! ---------------------------
  ! Returns the perturbation magnitude of a field at the specified
  ! index in the state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT pert_mag          perturbation magnitude
  subroutine get_pert_magnitude_by_index(var_index, pert_mag)
    integer, intent(in)        :: var_index
    real(kind=r8), intent(out) :: pert_mag

    pert_mag = restart_info%restart_vars(var_index)%pert_mag
  end subroutine get_pert_magnitude_by_index

  ! get_pert_type_by_index
  ! ----------------------
  ! Returns the perturbation type (no perturbation, entire field 
  ! perturbed by same value, or each point perturbed with a random
  ! value) of a field at the specified index in the state vector.
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT pert_type         perturbation type (defined as integer
  !                         constants in this module)
  subroutine get_pert_type_by_index(var_index, pert_type)
    integer, intent(in)  :: var_index
    integer, intent(out) :: pert_type
    
    pert_type = restart_info%restart_vars(var_index)%pert_type
  end subroutine get_pert_type_by_index

  ! get_kind_by_index
  ! -----------------
  ! Returns the variable kind of a field at the specified index in
  ! the state vector
  !  PARAMETERS
  !   IN  var_index         variable index in state vector to query
  !   OUT var_kind          numeric variable kind for that variable
  subroutine get_kind_by_index(var_index, var_kind)
    integer, intent(in)  :: var_index
    integer, intent(out) :: var_kind
    
    var_kind = restart_info%restart_vars(var_index)%var_type
  end subroutine get_kind_by_index

  ! get_update_flag_by_index
  ! ------------------------
  ! Call get_update_flag using a variable at the specified index in
  ! the state vector definition.
  !  PARAMETERS
  !   IN  var_index         index of variable to interrogate
  !   OUT update            true if this variable gets updated
  subroutine get_update_flag_by_index(var_index, update)
    integer, intent(in)  :: var_index
    logical, intent(out) :: update

    if (var_index > size(restart_info%restart_vars)) then
       call error_handler(E_ERR, 'get_update_flag_by_index',      &
                          'Variable index out of range', source,  &
                          revision, revdate)
    else
       call get_update_flag(restart_info%restart_vars(var_index), &
                            update)
    end if
  end subroutine get_update_flag_by_index

  ! get_var_info_by_abs_index
  ! -------------------------
  ! Maps an entry in the long state vector to a particular field,
  ! then returns that field's name and sigma level index.
  !  PARAMETERS
  !   IN  abs_index         Arbitrary index of an element in the
  !                         state vector
  !   OUT var_name          Name of the field containing the element
  !   OUT sigma_index       Sigma index of the field containing the 
  !                         element
  subroutine get_var_info_by_abs_index(abs_index, var_name, level)
    integer, intent(in)                      :: abs_index
    character(len=VAR_NAME_LEN), intent(out) :: var_name
    integer, intent(out)                     :: level

    integer :: invar
    integer :: gridsize

    call get_grid_field_size(restart_info%grid, gridsize)

    invar = int(abs_index / gridsize) + 1

    var_name = restart_info%restart_vars(invar)%var_name
    level = restart_info%restart_vars(invar)%sigma_record
  end subroutine get_var_info_by_abs_index

  ! dump_restart_vars
  ! -----------------
  ! Write out all the components of the state vector
  !  PARAMETERS
  !   [none]
  subroutine dump_restart_vars()
    integer :: ii

    if (do_output()) then
       write (*,*) "*** Restart File Contents ***"
       write (*,*) "-----------------------------"
    end if
    
    do ii=1,size(restart_info%restart_vars)
       call dump_restart_var(restart_info%restart_vars(ii))
    end do
  end subroutine dump_restart_vars

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! read_state_var_file
  ! ---------------------
  ! Reads in the state vector definition file - the first line 
  ! is the number of variables to follow, then each line after that
  ! contains information on each variables within the state vector
  !  PARAMETERS
  !   IN  statevec_def_unit    Unit number for the state vector
  !                            definition file
  subroutine read_state_var_file(statevec_def_unit)
    integer, intent(in) :: statevec_def_unit

    integer                         :: num_flds
    integer                         :: num_recs

    type(restart_var)               :: cur_var
    type(state_vec_iterator)        :: iterator

    character(len=*), parameter     :: routine = 'read_state_var_file'
    logical                         :: is_opened
    integer                         :: io_status
    integer                         :: alloc_status

    ! Things we read in from the file
    integer                         :: num_lvls  ! num of levels
    real(kind=r8)                   :: pert_mag  ! perturbation 
    logical                         :: mean_flag ! Mean or not
                                                 ! magnitude
    character                       :: lvl_type  ! mass or w level
    character(len=VAR_NAME_LEN)     :: name      ! Variable name
    character(len=PERT_TYPE_LEN)    :: type_name ! text pert type
    character(len=VAR_TYPE_LEN)     :: var_type  ! variable type 
    character(len=FILE_TYPE_LEN)    :: file_type ! file type 
    character(len=FIELD_UPDATE_LEN) :: write_fld ! write into restart
    character(len=POS_DEF_LEN)      :: pos_def   ! Positive definite?
    
    ! Things that we convert what we read in the file to
    integer                         :: pert_type    ! Numeric pert type
    logical                         :: lvl_flag     ! true if mass level
    logical                         :: write_flag   ! Write/Skip flag
    logical                         :: pos_flag     ! positive definite?
    integer                         :: var_type_num ! numeric var type
    logical                         :: flat_file_io ! numeric var type

    ! Counter indicies
    integer :: lvl_beg, lvl_end
    integer :: cur_fld
    integer :: cur_lvl
    integer :: cur_index = 0

    ! variables to deal with special case mean fields
    integer, parameter    :: NUM_DEFAULT_VARS=3
    character(len=VAR_NAME_LEN), dimension(NUM_DEFAULT_VARS) :: default_vars
    data default_vars /'THBM','EXBM','EXBW'/

    ! Assert
    inquire(unit=statevec_def_unit, opened=is_opened)
    if (.not. is_opened) then
       call error_handler(E_ERR, 'read_state_var_file',          &
                          'State vector definition file not open', &
                          source, revision, revdate)
    end if

    read(unit=statevec_def_unit,fmt='(A)',iostat=io_status) file_type
    call check_io_status(io_status, routine, source, revision, revdate, &
                         'Read file type in state.vars')
    call file_type_to_flag(file_type, flat_file_io)
    
    num_flds = 0
    call initialize_list(state_list)
    call initialize_list(var_list)

    do while(.true.) 
    ! Format is:
    !  1. Variable name - from the model
    !  2. Sigma Level or total number of levels
    !  3. Perturbation type (UNIFORM, INDIVID, NOPERTS)
    !  4. Perturbation magnitude
    !  5. Mass/W level (should be either 'M' or 'W')
    !  6. Variable type (from the raw observation kinds)
    !  7. Update field in restart vector or no?
    !  8. Whether field is a mean or not (TRUE/FALSE)
    !  9. Positive definiteness (ISPOSDEF, NOPOSDEF)
       read(unit=statevec_def_unit, fmt=*, iostat=io_status, end=300) &
            name, num_lvls, type_name, pert_mag, lvl_type, var_type,  &
            write_fld, mean_flag, pos_def

       num_flds = num_flds + 1
       call pert_name_to_type(type_name, pert_type)
       call level_type_to_flag(lvl_type, lvl_flag)
       call update_field_to_flag(write_fld, write_flag)
       call pos_def_to_flag(pos_def, pos_flag)
       var_type_num = get_raw_obs_kind_index(var_type)

       if(flat_file_io) then
         lvl_beg=1 ; lvl_end=num_lvls
       else
         lvl_beg=num_lvls ; lvl_end=num_lvls ; num_lvls = 1
       end if

       do cur_lvl=lvl_beg,lvl_end 
         call populate_info(name, cur_lvl, num_lvls, pert_type, pert_mag,  &
                            var_type_num, write_flag, mean_flag,           &
                            lvl_flag, pos_flag, flat_file_io, cur_var)
         call add_to_list(cur_var, state_list)
         if(cur_lvl == lvl_beg) call add_to_list(cur_var, var_list) 
       end do
    end do
300 continue

    ! 
    restart_info%num_flds     = num_flds
    restart_info%flat_file_io = flat_file_io

    if(flat_file_io) then
      do cur_fld = 1,NUM_DEFAULT_VARS
        call set_mean_state_var(default_vars(cur_fld))
      end do
    end if
    
    num_recs = state_list%list_length
    allocate(restart_info%restart_vars(num_recs), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,    &
                            revdate, 'restart_info%restart_vars')

    iterator = get_iterator(state_list)
    do while( has_next(iterator) )
       cur_index = cur_index + 1
       restart_info%restart_vars(cur_index) = get_next(iterator)
    end do

    call print_list(var_list)

  end subroutine read_state_var_file

  ! level_type_to_flag
  ! ------------------
  ! Convert character representation of [M]ass or [W] velocity level
  ! type to an internal logical type 
  !  PARAMETERS
  !   IN  level_type_char      level type as character (M/W)
  !   OUT is_mass_level        True if levels is [M]ass
  subroutine level_type_to_flag(level_type_char, is_mass_level)
    character, intent(in)  :: level_type_char
    logical, intent(out)   :: is_mass_level

    select case (adjustl(level_type_char))
       case (adjustl(MASS_LEVEL))
          is_mass_level = .true.
       case (adjustl(W_LEVEL))
          is_mass_level = .false.
       case default
          call error_handler(E_ERR, 'level_type_to_flag', 'Level ' //    &
                             'type' // level_type_char // ' not found!', &
                             source, revision, revdate)
       end select
  end subroutine level_type_to_flag

  ! update_field_to_flag
  ! --------------------
  ! Convert character representation of whether the field should be
  ! written back into the COAMPS restart file (UPDATE) or skipped
  ! (FREEZE) to the internal logical representation
  !  PARAMETERS
  !   IN  update_string     UPDATE or FREEZE
  !   OUT update_flag       consistent with module-defined flag
  subroutine update_field_to_flag(update_string, update_flag)
    character(len=FIELD_UPDATE_LEN), intent(in) :: update_string
    logical, intent(out)                        :: update_flag

    select case (adjustl(update_string))
    case (adjustl(FIELD_UPDATE))
       update_flag = FLAG_UPDATE_FIELD
    case (adjustl(FIELD_FREEZE))
       update_flag = FLAG_FREEZE_FIELD
    case default
       call error_handler(E_ERR, 'update_field_to_flag',          &
                          'Update type' // update_string //       &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine update_field_to_flag

  ! file_type_to_flag
  ! ---------------
  ! Converts a character string describing whether the entire field
  ! is read from a flat file (FLAT) or a restart file (RESTART) to 
  ! the interna llogical representation 
  !  PARAMETERS
  !   IN  file_type_string     FLAT or RESTART
  !   OUT flat_flag       consistent with module-defined flag
  subroutine file_type_to_flag(file_type_string, file_type_flag)
    character(len=FILE_TYPE_LEN), intent(in) :: file_type_string
    logical, intent(out)                     :: file_type_flag
    select case (adjustl(file_type_string))
    case (adjustl(FLAT_FILE_STRING))
       file_type_flag = FLAG_FLAT
    case (adjustl(RESTART_FILE_STRING))
       file_type_flag = FLAG_RESTART
    case default
       call error_handler(E_ERR, 'file_type_to_flag', &
                          'FILE_TYPE ' // file_type_string // &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine file_type_to_flag

  ! pos_def_to_flag
  ! ---------------
  ! Converts a character string describing whether the entire field
  ! is positive definite (ISPOSDEF) or not (NOPOSDEF) to the internal
  ! logical representation 
  !  PARAMETERS
  !   IN  posdef_string     ISPOSDEF or NOPOSDEF
  !   OUT posdef_flag       consistent with module-defined flag
  subroutine pos_def_to_flag(posdef_string, posdef_flag)
    character(len=POS_DEF_LEN), intent(in) :: posdef_string
    logical, intent(out)                   :: posdef_flag

    select case (adjustl(posdef_string))
    case (adjustl(ALLOW_NEGATIVE))
       posdef_flag = FLAG_ALLOW_NEGATIVE
    case (adjustl(FORCE_POSITIVE))
       posdef_flag = FLAG_FORCE_POSITIVE
    case default
       call error_handler(E_ERR, 'pos_def_to_flag', &
                          'Positive definite ' // posdef_string // &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine pos_def_to_flag

  ! pert_name_to_type
  ! -----------------
  ! Converts a character string describing a perturbation type to the
  ! internal numeric representation
  !  PARAMETERS
  !   IN  pert_name         NOPERTS, INDIVID, UNIFORM
  !   OUT pert_type         consistent with module-defined values
  subroutine pert_name_to_type(pert_name, pert_type)
    character(len=PERT_TYPE_LEN), intent(in) :: pert_name
    integer, intent(out)                     :: pert_type

    select case (pert_name)
    case (NOPERTS_NAME)
       pert_type = PERT_TYPE_NOPERTS
    case (UNIFORM_NAME)
       pert_type = PERT_TYPE_UNIFORM
    case (INDIVID_NAME)
       pert_type = PERT_TYPE_INDIVID
    case default
       call error_handler(E_ERR, 'pert_name_to_type',             &
                          'Perturbation type ' // pert_name //    &
                          ' not found!', source, revision, revdate)
    end select
  end subroutine pert_name_to_type

  ! populate_info
  ! -------------
  ! Populates a restart_var structure with the information read in
  ! from the state vector definition file.  The dimension type and 
  ! variable record in that dimension are looked up given the name
  ! specified.
  !  PARAMETERS
  !   IN  var_name          the variable name
  !   IN  sig_index         sigma level index for this variable
  !   IN  num_lvls          Number of contiguous levels
  !   IN  pert_type         perturbation type (after conversion)
  !   IN  pert_magnitude    perturbation magnitude
  !   IN  var_type          the variable type (see the list of raw
  !                         types in obs_kind_mod definition)
  !   IN  update_flag       True if field is to be written back to
  !                         COAMPS restart file
  !   IN  mean_flag         True if field is a mean field
  !   IN  level_flag        True if field is a mass level
  !   IN  pos_flag          True if field is restricted to be non-negative
  !   IN  file_flag         True if input_fields are field from flat files
  !   OUT cur_var           restart_var structure containing the 
  !                         supplied information as well as the 
  !                         dimension type/location in restart file
  subroutine populate_info(var_name, sig_index, num_lvls, pert_type,   &
                           pert_magnitude, var_type, update_flag,      &
                           mean_flag, level_flag, pos_flag, file_flag, &
                           cur_var) 
    character(len=VAR_NAME_LEN), intent(in)   :: var_name
    integer, intent(in)                       :: sig_index
    integer, intent(in)                       :: num_lvls
    integer, intent(in)                       :: pert_type
    real(kind=r8), intent(in)                 :: pert_magnitude
    integer, intent(in)                       :: var_type
    logical, intent(in)                       :: update_flag
    logical, intent(in)                       :: mean_flag
    logical, intent(in)                       :: level_flag
    logical, intent(in)                       :: pos_flag
    logical, intent(in)                       :: file_flag
    type(restart_var), intent(inout)          :: cur_var

    integer               :: cur_dim_type
    integer, dimension(2) :: cur_var_record

    ! These values we already have
    cur_var%var_name = var_name
    cur_var%sigma_record = sig_index
    cur_var%pert_mag = pert_magnitude
    cur_var%pert_type = pert_type
    cur_var%var_type = var_type
    cur_var%update_field = update_flag
    cur_var%mean_field = mean_flag
    cur_var%mass_level = level_flag
    cur_var%positive_definite = pos_flag

    if(sig_index .eq. 1) then
      cur_var%num_lvls = num_lvls
    else
      cur_var%num_lvls = 0
    end if

    ! Match the name up with the corresponding dimension/record
    ! information
    if(file_flag) then
      cur_var%dim_type = 0
      cur_var%var_record = 0
    else
      call get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,    &
                         SINGLEIO, MULTIIO, var_name, cur_dim_type, &
                         cur_var_record)
      cur_var%dim_type = cur_dim_type
      cur_var%var_record = cur_var_record
    end if
  end subroutine populate_info

  ! get_full_record_num
  ! -------------------
  ! Queries the given variable and returns the record number of a
  ! particular sigma level in the large COAMPS restart file
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   IN  use_singleio      true if restart file is written using
  !                         a single I/O processor
  !   OUT recnum            variable's record number within the
  !                         restart file
  subroutine get_full_record_num(cur_var, use_singleio, recnum)
    type(restart_var), intent(in) :: cur_var
    logical, intent(in)           :: use_singleio
    integer, intent(out)          :: recnum

    ! How many total variables of each kind there are - use this for
    ! calculating offsets between dimensions
    integer :: N1D, N3D

    ! Temporary shorthand so we don't need to keep referencing the
    ! structure in the select statement
    integer :: kka, k, n

    ! How many of each kind of variable we have depends on if we are
    ! doing single or multiple process I/O
    if (use_singleio) then
       N1D = 101
       N3D = 64
    else
       N1D = 112
       N3D = 56
    end if

    kka = restart_info%grid%sigm_lvls
    k   = cur_var%sigma_record
    
    if (use_singleio) then
       n = cur_var%var_record(SINGLEIO)
    else
       n = cur_var%var_record(MULTIIO)
    end if

    ! Calculate offsets - do this by taking into account the number
    ! of variables that occur before the target variable and all
    ! their associated sigma levels, then all the sigma levels of
    ! THIS variable that have occured before the target level.
    select case (cur_var%dim_type)
    case (DIM_TYPE_2D)
       recnum = n
    case (DIM_TYPE_3D)
       recnum = (N1D) + ((kka)*(n-1) + 1) + (k - 1)
    case (DIM_TYPE_3DW)
       recnum = (N1D + N3D*kka) + ((kka + 1)*(n - 1) + 1) + (k - 1)
    case default
       print *,cur_var%dim_type,DIM_TYPE_2D,DIM_TYPE_3D,DIM_TYPE_3DW
       call error_handler(E_ERR,'get_full_record_num',            &
                          'Unrecognized dimension type!', source, &
                          revision, revdate)
    end select
  end subroutine get_full_record_num

  ! get_posdef_flag
  ! ---------------
  ! Queries the given variable and returns whether the variable is
  ! constrained to never be negative (e.g. mixing ratios)
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   OUT posdef            true if field shouldn't be negative
  subroutine get_posdef_flag(cur_var, posdef)
    type(restart_var), intent(in) :: cur_var
    logical, intent(out)          :: posdef
    posdef = cur_var%positive_definite
  end subroutine get_posdef_flag

  ! get_vert_coord
  ! --------------
  ! Queries the given variable and returns whether the variable's
  ! vertical coordinate type and vertical coordinate location.
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   IN  grid              coamps_grid to use for coordinate calcs
  !   OUT vert_type         vertical coordinate type - currently
  !                         either surface or model level
  !   OUT vert_level        0 at surface or the sigma level
  subroutine get_vert_coord(cur_var, grid, vert_type, vert_level)
    type(restart_var), intent(in) :: cur_var
    type(coamps_grid), intent(in) :: grid
    integer, intent(out)          :: vert_type
    real(kind=r8), intent(out)    :: vert_level

    integer :: sigma_index

    sigma_index = cur_var%sigma_record

    ! Handle three cases - surface variables, on vertical mass levels
    ! and on vertical w levels
    if (sigma_index .eq. 0) then
       vert_type = VERTISSURFACE
       vert_level = 0.0
    else
       vert_type = VERTISLEVEL
       if (cur_var%mass_level) then
          vert_level = grid%msigma(sigma_index)
       else
          vert_level = grid%wsigma(sigma_index)
       end if
    end if
  end subroutine get_vert_coord

  ! get_update_flag
  ! ---------------
  ! Queries the given variable and returns whether the variable is
  ! written to the COAMPS restart file or not.
  !  PARAMETERS
  !   IN  cur_var           restart_var to pull data from
  !   OUT update            true if variable is written back to the
  !                         restart file
  subroutine get_update_flag(cur_var, update)
    type(restart_var), intent(in) :: cur_var
    logical, intent(out)          :: update
    update = cur_var%update_field
  end subroutine get_update_flag

  ! dump_restart_var
  ! ----------------
  ! Writes out the contents of a restart variable
  !  PARAMETERS
  !   IN  var               restart_var to dump contents of
  subroutine dump_restart_var(var)
    type(restart_var), intent(in) :: var

    character(len=32) :: kind_name

    kind_name = get_raw_obs_kind_name(var%var_type)

    if (do_output()) then
       write (*, '(A9 T15A10)')  "VAR NAME:",var%var_name
       write (*, '(A9 T15I8.8)') "NUM LVLS:",var%num_lvls
!       write (*, '(A9 T15I8.8)') "DIM TYPE:",var%dim_type
!       write (*, '(A9 T15I8.8)') "VAR RCRD:",var%var_record
       write (*, '(A9 T15I8.8)') "SGM RCRD:",var%sigma_record
       write (*, '(A9 T15F8.6)') "PTRB PCT:",var%pert_mag
       write (*, '(A9 T15I8.8)') "PTB TYPE:",var%pert_type
       write (*, '(A9 T15L1)')   "UPDATE??:",var%update_field
       write (*, '(A9 T15A32)')  "VAR TYPE:",kind_name
       write (*, '(A9 T15I8.8)') "VAR NMBR:",var%var_type
       write (*, '(A9 T15L1)')   "MASSLVL?:",var%mass_level
       write (*, '(A9 T15L1)')   "MEANFLD?:",var%mean_field
       write (*, *) "-----------------------------------------"
    end if
  end subroutine dump_restart_var

  ! set_mean_state_var
  ! ---------------------
  ! Sets state variable info for mean fields that are not in flat
  ! files and are defined internally (e.g. mean fields)
  ! 
  !  PARAMETERS
  !   IN     var_name    Variable name
  !   INOUT  var_index   Begining index location of variable
  !                      in state var structure.
  subroutine set_mean_state_var(var_name)
    character(len=VAR_NAME_LEN), intent(in)    :: var_name 
    type(restart_var) :: cur_var

    integer                         :: num_lvls  ! num of levels
    character                       :: lvl_type  ! mass or w level
    character(len=VAR_TYPE_LEN)     :: var_type  ! variable type 
    logical                         :: lvl_flag  ! true if mass level
    integer                         :: var_type_num !numeric var type

    ! Variables that have fixed default values
    real(kind=r8)                   :: pert_mag=0.0     ! perturbation magnitude
    logical                         :: mean_flag=.true. ! Mean or not
    integer                         :: pert_type=PERT_TYPE_NOPERTS  ! Numeric pert type
    logical                         :: write_flag=FLAG_FREEZE_FIELD ! Write/Skip flag
    logical                         :: pos_flag=FLAG_ALLOW_NEGATIVE ! is positive definite?
    logical                         :: file_type=FLAG_FLAT          ! is flat_file?

    ! Counter variable
    integer                         :: cur_lvl

    call get_grid_num_levels(restart_info%grid,num_lvls)
    select case (trim(uppercase(var_name)))
      case ('EXBM')
        lvl_type = 'M' 
        var_type = 'KIND_EXNER_FUNCTION'
      case ('THBM')
        lvl_type = 'M' 
        var_type = 'KIND_POTENTIAL_TEMPERATURE'
      case ('EXBW')
        lvl_type = 'W' 
        var_type = 'KIND_EXNER_FUNCTION'
        num_lvls = num_lvls + 1
    end select

    call level_type_to_flag(lvl_type, lvl_flag)
    var_type_num = get_raw_obs_kind_index(var_type)

    do cur_lvl = 1, num_lvls
      call populate_info(var_name, cur_lvl, num_lvls, pert_type,   &
                         pert_mag, var_type_num, write_flag,       &
                         mean_flag, lvl_flag, pos_flag, file_type, &
                         cur_var)
      call add_to_list(cur_var, state_list)
    end do
  end subroutine set_mean_state_var

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LINKED LIST ROUTINES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! initialize_state_vec_list
    ! ---------------------------
    ! Constructor for state_vec linked list
    !  PARAMETERS
    ! INOUT list    List to initialize
    subroutine initialize_list(list)
      type(state_vec_list), intent(inout) :: list
      nullify(list%head)
      nullify(list%tail)
    end subroutine initialize_list

    ! add_state_vec_to_list
    ! -------------------
    ! Adds the given state_vec to a given linked list
    !  PARAMETERS
    !   IN  new_entry       Integer to add to the list
    ! INOUT list            List to add the entry to
    subroutine add_to_list(new_entry, list)
      type(restart_var),    intent(in)    :: new_entry
      type(state_vec_list), intent(inout) :: list
      ! Add a new node
      if (.not. associated(list%head)) then
          allocate(list%head) ; list%tail => list%head
          list%list_length = 1
      else
          allocate(list%tail%next) ; list%tail => list%tail%next
          list%list_length = list%list_length + 1
      end if
      ! Assign the contents
      nullify(list%tail%next) 
      list%tail%value = new_entry
    end subroutine add_to_list

    ! get_state_vec_list_iterator
    ! -----------------------------
    ! Returns an iterator for the given linked list
    !  PARAMETERS
    !   IN  list            list to iterate over
    function get_iterator(list) result(iterator_out)
        type(state_vec_list),     intent(in)  :: list
        type(state_vec_iterator)              :: iterator_out
        iterator_out%current   => list%head
        iterator_out%cur_index =  1 
    end function get_iterator

    ! has_next
    ! -----------------------------
    ! Returns true if the iterator can return another value
    !  PARAMETERS
    !   IN  iterator        the iterator to check the state of
    function has_next(iterator)  result(has_next_out)
        type(state_vec_iterator), intent(in)  :: iterator
        logical                               :: has_next_out
        has_next_out = associated(iterator%current)
    end function has_next

    ! iterator_get_next
    ! ------------------------------
    ! Returns the current value at the iterator location, ***and advances the
    ! iterator to the next entry***
    !  PARAMETERS
    !   IN  iterator        List iterator to get value of and advance
    function get_next(iterator) result(get_next_out)
        type(state_vec_iterator), intent(inout) :: iterator
        type(restart_var)                       :: get_next_out
        get_next_out       =  iterator%current%value
        iterator%current   => iterator%current%next  
        iterator%cur_index =  iterator%cur_index + 1
    end function get_next

    ! print_state_vec_list
    ! ------------------
    ! Traverses the linked list and prints the nest ID of each element
    !  PARAMETERS
    !   IN  list            linked list to iterate over
    subroutine print_list(list)
        type(state_vec_list), intent(in) :: list
        type(state_vec_iterator)         :: iterator
        iterator = get_iterator(list)
        do while (has_next(iterator))
            call dump_restart_var(get_next(iterator))
        end do
    end subroutine print_list

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------
end module coamps_restart_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
