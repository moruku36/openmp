!      -------------------------------
!      Poisson equation with k=1, f= sin(i*j)
!      ------------------------------


       program main
       use omp_lib
       include 'poisson.inc'

       integer DEBUG
       parameter (DEBUG=1)

       double precision EPS
       parameter (EPS=1.0e-4)

       integer MAX_ITER
       parameter (MAX_ITER=10000)

       integer MAX_HEAT
       parameter (MAX_HEAT=100.0)


       double precision,allocatable,dimension(:,:):: U, U_old
       double precision  t0, t1, t2, t_w
       double precision  h, h_pow, dmax
     
       integer i, j
       integer ii, jj
    
       allocate(U(0:M+1,0:M+1))
       allocate(U_old(0:M+1,0:M+1))


       do j=0, M+1
         do i=0, M+1
           U(i, j) = 0.0d0
           U_old(i, j) = 0.0d0
         enddo  
       enddo
   
c      === Set given temperature ----------------------------------
       do j=0, M+1
         U(0, j) = MAX_HEAT
       enddo 
       do i=1, M
         U(i, 0) = MAX_HEAT - dble(i)/dble(M)*MAX_HEAT
         U(i, M+1) = dble(i)/dble(M)*MAX_HEAT
       enddo
       do j=0, M+1
         U(M+1, j) = 0.0d0
       enddo 
c      === End of setting given temperature ------------------------


c      === Start of Solving Poisson equation ----------------------------
       t1 = omp_get_wtime()

       h = 1.0d0 / (dble(M+1))
       h_pow = h*h

c      === Main loop
       do i=1, MAX_ITER

c        --- perform explicit method ( Gauss-Seidel Method )
         call MyPoisson(U, h_pow, M) 
c        --- compute maximum differences
         call CalcErr(U, U_old, M, dmax)

         if (mod(i,100) .eq. 0) then
            print *, " iter=",i, "dmax=",dmax
         endif

         if (dmax .lt. EPS) then
           print *, i, " Iteration is converged in residual ", EPS
           goto 10
         endif

c        --- copy back to U_old
         do jj=1, M
           do ii=1, M
             U_old(ii, jj) = U(ii, jj) 
           enddo  
         enddo        

       enddo 
c      === End of Main loop

       print *, "Iteration is not converged within " 
     &       ,MAX_ITER, " times." 

  10   continue

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of Solving Poisson equation ---------------------------

         print *, "M  = ", M
         print *, "MAX_ITER  = ", MAX_ITER

         print *, "time[sec.] = ", t_w

c        === output results  
!         do ii=1, M
!           do jj=1, M
!              print *, ii, jj, U(ii,jj)
!           enddo
!           print *, ""
!         enddo


       deallocate(U)
       deallocate(U_old)

       stop
       end


       subroutine MyPoisson(U_rhs, h_pow, M)
       use omp_lib

       double precision U_rhs(0:M+1, 0:M+1)
       double precision h_pow
       integer M 

       integer i, j  

!      === u_{i,j} = 1/4 (h^2 f_{i,j} + u_{i,j-1} +u_{i-1,j}+u_{i+1,j}+u_{i,j+1})
       do j=1, M
         do i=1, M
           U_rhs(i, j) = 0.25d0 * 
     &      ( h_pow * dsin(dble(i)*dble(j)) +
     &       U_rhs(i,j-1) + U_rhs(i-1,j) + U_rhs(i+1,j) + U_rhs(i,j+1) )
 
         enddo
       enddo

       return
       end


       subroutine CalcErr(U_rhs, U_lhs, M, dmax)
       use omp_lib

       double precision U_rhs(0:M+1, 0:M+1)
       double precision U_lhs(0:M+1, 0:M+1)
       double precision dmax
       integer M 

       integer i,j
       double precision dtemp 
      

       dmax = 0.0
       do j=1, M
         do i=1, M
           dtemp = dabs(U_rhs(i, j) - U_lhs(i, j))
           if (dtemp .gt. dmax ) dmax = dtemp
         enddo
       enddo 
 
       return
       end
