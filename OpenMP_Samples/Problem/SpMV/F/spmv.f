       program main
       use omp_lib
       include 'spmv.inc'

       integer DEBUG
       parameter (DEBUG=1)
       double precision EPS
       parameter (EPS=1.0e-18)
       integer MAX_ITER
       parameter (MAX_ITER=100)


       double precision,allocatable,dimension(:):: X, Y, VAL
       double precision  t0, t1, t2, t_w
       double precision  d_mflops, dtemp

       integer,allocatable,dimension(:):: IRP, ICOL, ISEED, ICOL_LIST
       integer  i, j, k, n, kk
       integer  iflag, iflag_t

!       integer IRP(NN+1)
!       integer ICOL(NNNZ)
!       double precision X(NN)
!       double precision Y(NN)
!       double precision VAL(NNNZ)

       allocate(IRP(NN+1))
       allocate(ICOL(NNNZ))
       allocate(X(NN))
       allocate(Y(NN))
       allocate(VAL(NNNZ))
       allocate(ICOL_LIST(NN))
   

c      === sparse matrix generation --------------------------
       CALL RANDOM_SEED (SIZE = i)
       allocate(ISEED(i))
       ISEED = 1

       call RANDOM_SEED(PUT=ISEED(1:i))


       k = 1
       IRP(1) = 1
       do i=1, NN

         do j=1, NN
           ICOL_LIST(j) = -1
         enddo
         do j=1, NZPR
 1         continue 
           call RANDOM_NUMBER(dtemp)
           kk = int(dble(NN)*dtemp)
           if (kk .eq. 0) kk = 1
           if (ICOL_LIST(kk) .ne. -1) then
             goto 1
           else
             ICOL_LIST(kk) = 1
           endif
         enddo
         do j=1, NN
           if (ICOL_LIST(j) .eq. 1) then
             ICOL(k) = j
             VAL(k) = 1.0d0 
             k = k + 1
           endif
         enddo
         IRP(i+1) = k 
         X(i) = 1.0d0
         Y(i) = 0.0d0
       enddo
c      === end of sparse matrix generation ------------------------

c      === Start of mat-vec routine ----------------------------
       t1 = omp_get_wtime()

       do i=1, MAX_ITER
         call MySpMV(Y, IRP, ICOL, VAL, X, NN, NNNZ)
       enddo 

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of mat-vec routine ---------------------------


         print *, "NN  = ", NN
         print *, "NNNZ  = ", NN*NZPR
         print *, "NZPR  = ", NZPR
         print *, "MAX_ITER  = ", MAX_ITER

         print *, "SpMV time[sec.] = ", t_w

         d_mflops = (2.0*dble(NZPR*NN)*dble(MAX_ITER)) /t_w
         d_mflops = d_mflops * 1.0e-6
         print *, "MFLOPS = ", d_mflops

       if (DEBUG .eq. 1) then
c        === Verification routine ----------------- 
         iflag = 0
         do i=1, NN
           if (dabs(Y(i) - dble(NZPR)) > EPS) then
             print *, " Error! in (", i, ") th argument.", Y(i)
             iflag = 1
             goto 10
           endif 
         enddo
 10      continue
c        -----------------------------------------

         if (iflag .eq. 0) then
           print *, " OK!"
         endif

       endif
c      -----------------------------------------

       deallocate(IRP)
       deallocate(ICOL)
       deallocate(X)
       deallocate(Y)
       deallocate(VAL)
       deallocate(ISEED)
       deallocate(ICOL_LIST)
       stop
       end



       subroutine MySpMV(Y, IRP, ICOL, VAL, X, n, nnz)
       use omp_lib
       include 'spmv.inc'

       integer IRP(NN+1)
       integer ICOL(NNNZ)
       double precision X(NN)
       double precision Y(NN)
       double precision VAL(NNNZ)

       integer n, nnz

       double precision s
       integer  i, j_ptr


       do i=1,n
         s = 0.0d0
         do j_ptr=IRP(i),IRP(i+1)-1
            s = s + VAL(j_ptr)*X(ICOL(j_ptr))
         enddo
         Y(i) = s
       enddo 


       return
       end
