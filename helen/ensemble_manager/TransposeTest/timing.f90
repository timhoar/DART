!----------------------------------------------------
! Timing of transpose
!
!> @author
!> Helen 
!
!> DESCRIPTION
!> Aim to look at the scaling of the transpose 
!> all_vars_to_all_copies and all_copies_all_vars
!> my_task_id() 0 timings output to standard out
!> All my_task_id()s timings output to standard error
!> can use cpu time or MPI_WTIME (namelist flag)
!-----------------------------------------------------      
program timing
      use mpi_utilities_mod,    only : initialize_mpi_utilities,finalize_mpi_utilities,                &
                                       my_task_id, task_sync, task_count
      use utilities_mod,        only: find_namelist_in_file, check_namelist_read,                      &
                                      set_nml_output
      use ensemble_manager_mod, only: init_ensemble_manager, ensemble_type, all_vars_to_all_copies,    &
                                      all_copies_to_all_vars, end_ensemble_manager, get_my_num_vars,   &
                                      get_my_num_copies!, my_pe

       use mpi

implicit none

type(ensemble_type) :: ens_handle, obs_ens_handle, forward_op_ens_handle !> ensemble type

!namelist with default values
integer :: timing_ens_size = 80 !> ensemble size 
integer :: ens_size_extras = 6 !> number of extra columns to add on to ens_size
integer :: n_repeats = 1 !> number of repetitions of loop for a given array size
integer :: min_array_length = 20 !> minimum number of vars
integer :: inc_array_length = 10 !> step size for number of vars loop
integer :: max_array_length = 50 !> max number of vars
logical :: record_transpose = .false.

!timing variables
integer :: array_length !> length of array to transpose
integer :: i, k, j !> for loops
real:: cpu_start, cpu_finish !> for cpu_time
double precision:: start, finish !> for MPI_Wtime

! utilites
integer :: iunit !> for reading in namelist from flle
integer :: async = 0 !> what does this do for mpi_finalize
integer :: io !> for iostat
integer :: stderr = 0 !> for output to standard error
integer :: ierr !> for MPI_BCAST
character*120 :: task_str, file0, file1, file2
! namelist : debug flag,  ensemble size (80), number of extras (4 or 6 normally), array length( min, inc, max), number of repitions
! use_cpu_time (true/false) 
namelist /timing_nml/  &
 timing_ens_size, ens_size_extras, &
 min_array_length, inc_array_length, max_array_length, n_repeats, record_transpose

! initialize mpi
 call initialize_mpi_utilities('timing')

 call find_namelist_in_file("input.nml", "timing_nml", iunit)
 read(iunit, nml = timing_nml, iostat = io)

 call check_namelist_read(iunit, io, "timing_nml")

! print run info
 if (my_task_id() == 0) print*, 'ens_size', timing_ens_size, '+ extras', ens_size_extras

! print output format to screen
 if (my_task_id() == 0) then
   print*, 'Timing output: routine  array_length  repetition#: time'
   write(stderr, *) 'Timing output: my_task_id() routine array_length repetition# : time'
 end if
 
     numvar : do array_length = min_array_length, max_array_length, inc_array_length

        numRep : do i=1, n_repeats

           call task_sync()

           ! setup problem: initialize the ensemble manager storage - Need enough copies for ensemble plus extras
           call init_ensemble_manager(ens_handle, timing_ens_size + ens_size_extras, array_length, 1)

           ! open output files: my_pe needs to be public in ensemble manager
           if (i == 1 .and. record_transpose .eqv. .true. ) then
             !write(task_str, '(i10)') my_pe
             write(task_str, '(i10)') ens_handle%my_pe
             file0 = TRIM('inital' // TRIM(ADJUSTL(task_str)) // '.trans')
             file1 = TRIM('intermediate_output' // TRIM(ADJUSTL(task_str)) // '.trans')
             file2 = TRIM('final_output' // TRIM(ADJUSTL(task_str)) // '.trans')

             open(15, file=file0, status ='new')
             open(20, file=file1, status ='new') ! error if you already have results files
             open(30, file=file2, status ='new')
             !print *, 'my_pe', my_pe, 'my_task_id', my_task_id()

           endif

           ! set up the ensemble - var complete
           !  Note the use of my_pe not my_task_id()
           do j = 1, ens_handle%my_num_copies
             do k = 1, ens_handle%num_vars
!               ens_handle%vars(k,j) = j*100000 + (my_pe + 1)*1000000 + k
               ens_handle%vars(k,j) = j*100000 + (ens_handle%my_pe + 1)*1000000 + k
             enddo
           enddo

           if (i==1) then
             if( my_task_id()==0) then
               print *, 'task', my_task_id(), 'num copies = ',  get_my_num_copies(ens_handle), 'my num vars = ', get_my_num_vars(ens_handle)
               print *, 'num_vars               = ', array_length
               print *, 'ens_size               = ', timing_ens_size + ens_size_extras
               print *, 'average message_length = ', array_length / task_count()
             endif
           endif

           if ( record_transpose .eqv. .true. ) then
             do j = 1, ens_handle%num_vars
               write(15, *) ens_handle%vars(j,:)
              enddo
           endif

           call task_sync()

           ! call transpose all_vars_to_all_copies
           start = MPI_WTIME()
           call all_vars_to_all_copies(ens_handle)
           finish = MPI_WTIME()
           if (my_task_id() == 0 ) print *, 'all_vars_to_all_copies ', array_length, i,' : ', finish - start
           write(stderr, *)my_task_id(),' vars_to_copies',  array_length, i,' : ', finish - start

           ! write transpose results to file
           if ( record_transpose .eqv. .true. ) then
             do k = 1, ens_handle%num_copies
               write(20, *) ens_handle%copies(k,:)
              enddo
           endif

           call task_sync()

           ! call transpose all_copies_to_all_vars
           start = MPI_WTIME()
           call all_copies_to_all_vars(ens_handle)
           finish = MPI_WTIME()
           if (my_task_id() == 0 ) print *, 'all_copies_to_all_vars ', array_length, i,' : ', finish - start
           write(stderr, *)my_task_id(),' copies_to_vars',  array_length, i,' : ', finish - start

           call task_sync()

           ! write transpose results to file
           if ( record_transpose .eqv. .true. ) then
             do j = 1, ens_handle%num_vars
               write(30, *) ens_handle%vars(j, :)
             enddo
           endif

           ! destroy storage
           call end_ensemble_manager(ens_handle)
 
        end do numRep

     end do numvar

  close(15)
  close(20)
  close(30)

! finialize mpi
 call finalize_mpi_utilities(async=async)

end program timing
