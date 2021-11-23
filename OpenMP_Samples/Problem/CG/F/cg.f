       program main
       use omp_lib
       include 'cg.inc'

       integer DEBUG
       parameter (DEBUG=1)

       double precision EPS
       parameter (EPS=1.0e-8)

       double precision EPS_CG
       parameter (EPS_CG=1.0e-8)

       integer MAX_ITER
       parameter (MAX_ITER=100)


       double precision,allocatable,dimension(:):: X, B, VAL
       double precision,allocatable,dimension(:):: R, P, AX, AP, ADIAG
       double precision  t0, t1, t2, t_w
       double precision  d_mflops, dtemp
       double precision  alpha, beta, rdot, rdot_old, pAp
       double precision  xnorm, bnorm 

       integer,allocatable,dimension(:):: IRP, ICOL, ISEED
       integer  i, j, k, n
       integer  iflag, iflag_t

!       integer IRP(NN+1)
!       integer ICOL(NNNZ)
!       double precision X(NN)
!       double precision Y(NN)
!       double precision VAL(NNNZ)

       allocate(IRP(NN+1))
       allocate(ICOL(NNNZ))
       allocate(X(NN))
       allocate(B(NN))
       allocate(VAL(NNNZ))

       allocate(R(NN))
       allocate(P(NN))
       allocate(AX(NN))
       allocate(AP(NN))
       allocate(ADIAG(NN))
   
c      === sparse matrix generation --------------------------
       IRP(1) = 1
       VAL(1) = 4.0d0
       ICOL(1) = 1
       VAL(2) = -1.0d0
       ICOL(2) = 2
       IRP(2) = 3 
      
       k = 3    
       do i=2, NN-1
         ICOL(k) = i - 1
         VAL(k) = -1.0d0
         ICOL(k+1) = i 
         VAL(k+1) = 4.0d0
         ICOL(k+2) = i + 1 
         VAL(k+2) = -1.0d0
         
         k = k + 3
         IRP(i+1) = k 
       enddo
       VAL(k) = -1.0d0
       ICOL(k) = NN - 1
       VAL(k+1) = 4.0d0
       ICOL(k+1) = NN 
       IRP(NN+1) = k + 2

c      --- set diag elements of A
       do i=1, NN
         ADIAG(i) = 4.0d0
       enddo       
c      === end of sparse matrix generation ------------------------

c      === make right hand vector b -------------------------------
       do i=1, NN
         X(i) = 1.0d0
       enddo
       call MySpMV(B, IRP, ICOL, VAL, X, NN, NNNZ)    
c      === end of right hand vector b -----------------------------

c      === Start of CG ----------------------------
       t1 = omp_get_wtime()

c      === diagonal preconditioning
       do i=1, NN
         do j=IRP(i), IRP(i+1)-1
           VAL(j) = VAL(j) / ADIAG(i)
         enddo
         B(i) = B(i) / ADIAG(i)
       enddo 

c      === make x_0
       call RANDOM_SEED (SIZE = i)
       allocate(ISEED(i))
       ISEED = 1
       call RANDOM_SEED(PUT=ISEED(1:i))
       do i=1, NN
         call RANDOM_NUMBER(dtemp)
         X(i) = dtemp 
       enddo

c      === make p_0 = r_0 = b - A x_0  
       call MySpMV(AX, IRP, ICOL, VAL, X, NN, NNNZ)       
       do i=1, NN
         P(i) = B(i) - AX(i)
         R(i) = P(i)
       enddo


c      === CG main loop
       do i=1, MAX_ITER

c        --- alpha_k = (r_k, r_k) / (p_k, A p_k)
         call MySpMV(AP, IRP, ICOL, VAL, P, NN, NNNZ)
         pAp = 0.0d0 
         do j=1, NN
           pAp = pAp + P(j) * AP(j)  
         enddo 

         if (i .eq. 1) then
           rdot = 0.0d0
           do j=1, NN
             rdot = rdot + R(j) * R(j)
           enddo 
         endif

         alpha = rdot / pAp

c        --- x_{k+1} = x_{k} + alpha_k p_k
         do j=1, NN
           X(j) = X(j) + alpha * P(j)
         enddo 

c        --- Convergence check
         call MySpMV(AX, IRP, ICOL, VAL, X, NN, NNNZ)       
         do j=1, NN
           AX(j) = B(j) - AX(j)
         enddo          
         xnorm = 0.0d0
         do j=1, NN
           xnorm = xnorm + AX(j) * AX(j)
         enddo 
         xnorm = dsqrt(xnorm)

         print *, " iter=",i, "xnorm=",xnorm

         if (xnorm .lt. EPS_CG) then
           print *, "CG iteration is converged in residual ", EPS_CG
           goto 10
         endif


c        --- r_{k+1} = r_{k} - alpha_k A p_k
         do j=1, NN
           R(j) = R(j) - alpha * AP(j)
         enddo 

c        --- beta_k = (r_{k+1}, r_{k+1}) / (r_k, r_k)
         rdot_old = rdot
         rdot = 0.0d0
         do j=1, NN
           rdot = rdot + R(j) * R(j)
         enddo 
         beta = rdot / rdot_old

c        --- p_{k+1} = r_{k+1} + beta_k p_k
         do j=1, NN
           P(j) = R(j) + beta * P(j)
         enddo   

       enddo 
c      === End of CG main loop
       print *, "CG Iteration is not converged within " 
     &       ,MAX_ITER, " times." 

  10   continue

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of CG ---------------------------

         print *, "NN  = ", NN
         print *, "NNNZ  = ", NNNZ
         print *, "NZPR  = ", NZPR
         print *, "MAX_ITER  = ", MAX_ITER

         print *, "CG time[sec.] = ", t_w


       if (DEBUG .eq. 1) then
c        === Verification routine ----------------- 
         iflag = 0
         do i=1, NN
           if (dabs(X(i) - 1.0d0) > EPS) then
             print *, " Error! in (", i, ") th argument.", X(i)
             iflag = 1
             goto 100
           endif 
         enddo
 100     continue
c        -----------------------------------------

         if (iflag .eq. 0) then
           print *, " OK!"
         endif

       endif
c      -----------------------------------------

       deallocate(IRP)
       deallocate(ICOL)
       deallocate(X)
       deallocate(B)
       deallocate(VAL)
       deallocate(ISEED)

       deallocate(R)
       deallocate(P)
       deallocate(AX)
       deallocate(AP)
       deallocate(ADIAG)

       stop
       end



       subroutine MySpMV(Y, IRP, ICOL, VAL, X, n, nnz)
       use omp_lib
       include 'cg.inc'

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
