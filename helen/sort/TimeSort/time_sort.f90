program time_sort
! written April 17th 2013 by Helen Kershaw
! Aim: to time the the sort routine used when performing the in code task geometery 
! in init_ensemble_manager
!  The task_task_list file is created using taskGeomYellowstone.m
use ensemble_manager_mod, only : sort_task_list

implicit none

integer ::  i, num_pes
integer, allocatable :: list(:), idx(:)
character*32 :: filename
real :: start, finish

! read master_task_list from file
call getarg(1, filename)
open(20, file = filename, status='old')

 read(20, *) num_pes
 allocate(list(num_pes))
 allocate(idx(num_pes))
 do i = 1, num_pes
   read(20, *) list(i)
 enddo

close(20)

! time sort of master_task_list

call cpu_time(start)
 call sort_task_list(list, idx, num_pes)
call cpu_time(finish)

print*, 'Time for sort_task_list = ',  finish - start, 'seconds'

deallocate(list, idx)

end program time_sort


