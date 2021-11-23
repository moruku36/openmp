       program main
       use omp_lib
       include 'reduction.inc'

       double precision,allocatable,dimension(:):: A

       double precision  t0, t1, t2, t_w
       double precision  dtemp
       double precision  dtemp_t(0:MAX_THREADS)  

       integer istart(0:MAX_THREADS),iend(0:MAX_THREADS)         
       integer  i, k, ib

       allocate (A(NN))

c      === matrix generation --------------------------
       call RANDOM_SEED
       do i=1, NN
         call RANDOM_NUMBER(dtemp)   
       enddo
c      === end of matrix generation ------------------------

c      === Start of routine1 ----------------------------
       t1 = omp_get_wtime()

       do i=1, MAX_ITER
         call Reduction1(A, NN)
       enddo

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of routine1 ---------------------------

       print *, "NN  = ", NN
       print *, "MAX_ITER  = ", MAX_ITER
       print *, "Reduction1 time[sec.] = ", t_w


c      === Start of routine2 ----------------------------
       ib = NN / omp_get_max_threads()
       do i=0, omp_get_max_threads()-1
         dtemp_t(i) = 0.0d0
         istart(i) = 1 + ib*i
         iend(i) = (i+1)*ib
       enddo
       iend(omp_get_max_threads()-1) = NN

       t1 = omp_get_wtime()

       do i=1, MAX_ITER
         call Reduction2(A, NN, ib, dtemp_t, istart, iend)
       enddo

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of routine2 ---------------------------

       print *, "NN  = ", NN
       print *, "MAX_ITER  = ", MAX_ITER
       print *, "Reduction2 time[sec.] = ", t_w

c      === Start of routine3 ----------------------------
       t1 = omp_get_wtime()

!$omp parallel do
       do i=1, MAX_ITER
         call Reduction3(A, NN)
       enddo
!$omp end parallel do

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of routine3 ---------------------------

       print *, "NN  = ", NN
       print *, "MAX_ITER  = ", MAX_ITER
       print *, "Reduction3 time[sec.] = ", t_w




c      -----------------------------------------

       deallocate(A)

       stop
       end



       subroutine Reduction1(A, n)
       use omp_lib
       include 'reduction.inc'

       double precision A(NN)
       integer n

       integer  i
       double precision dtemp

       dtemp = 0.0d0
!$omp parallel do reduction(+: dtemp)  
       do i=1, n
         dtemp = dtemp + A(i) * A(i)    
       enddo
!$omp end parallel do

       return
       end


       subroutine Reduction2(A, n, ib, dtemp_t, istart, iend)
       use omp_lib
       include 'reduction.inc'

       double precision A(NN)
       integer n, ib
       double precision dtemp_t(MAX_THREADS)
       integer istart(MAX_THREADS), iend(MAX_THREADS) 

       integer i, k
       double precision dtemp

!$omp parallel do private(i)  
       do k=0, omp_get_max_threads()-1
         do i=istart(k), iend(k)
           dtemp_t(k) = dtemp + A(i) * A(i)    
         enddo
       enddo 
!$omp end parallel do
      
       dtemp = 0.0d0
       do k=0, omp_get_max_threads()-1
          dtemp = dtemp + dtemp_t(k)
       enddo

       return
       end


       subroutine Reduction3(A, n)
       use omp_lib
       include 'reduction.inc'

       double precision A(NN)
       integer n

       integer  i
       double precision dtemp

       dtemp = 0.0d0
       do i=1, n
         dtemp = dtemp + A(i) * A(i)    
       enddo

       return
       end
