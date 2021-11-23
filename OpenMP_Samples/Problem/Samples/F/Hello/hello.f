       program main
       use omp_lib
       
       integer i

!$omp parallel do
       do i=1, 10
         print *, "Hello parallel world!  i:", i
       enddo
!$omp end parallel do

       stop
       end

