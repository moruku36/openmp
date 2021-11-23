       program main
       use omp_lib
       include 'dem.inc'

       integer MAX_ITER
!       === for test
!       parameter (MAX_ITER=300)
       parameter (MAX_ITER=30)

       double precision,allocatable,dimension(:):: X
       double precision,allocatable,dimension(:):: V
       double precision,allocatable,dimension(:):: V_old
       double precision,allocatable,dimension(:):: F

       double precision  t0, t1, t2, t_w
       double precision  t3, t4, t_s
       double precision  d_mesh_dist, d_max_velocity 
       double precision  dtemp

       integer,allocatable,dimension(:)::   ISEED
       integer,allocatable,dimension(:,:):: MC
       integer,allocatable,dimension(:)::   MC_ind
       integer,allocatable,dimension(:)::   Ind_MC
       integer,allocatable,dimension(:)::   MC_compact_ind

       integer  i, j, k, n
       integer  imax_num_mesh, icompact_num_mesh
       integer  itotal_num_collision
       integer  i_color
 
       character(len=20) filename
       character(len=10) cbuf1
       character(len=4)  cbuf2
       character(len=20) cbuf3

       allocate(X(NN))
       allocate(V(NN))
       allocate(V_old(NN))
       allocate(F(NN))
       allocate(Ind_MC(NN))
      
       dtemp = 1.0d0/DPD
       imax_num_mesh = int(dtemp) + 1
       d_mesh_dist = DPD
       d_max_velocity = d_mesh_dist / DTS

       print *, "Number of Particles =", NN
       print *, "Particle Diameter =", DPD
       print *, "Particle Half Diameter =", DPD_HALF

       print *, "Time Step =", DTS
       print *, "Particle weight = ", PARTICLE_WEIGHT
 
       print *, "imax_num_mesh =", imax_num_mesh 
       print *, "d_mesh_dist =", d_mesh_dist
       print *, "d_max_velocity =", d_max_velocity

      
       if (imax_num_mesh .lt. NN ) then
         print *, " !!The number of particles needs to reduce."
         stop
       endif


       allocate(MC(1:imax_num_mesh, MAX_PARTICLES_IN_MESH))
       allocate(MC_ind(1:imax_num_mesh))
       allocate(MC_compact_ind(1:imax_num_mesh))

       do i=1, imax_num_mesh
         MC_ind(i) = 0
       enddo

c      === make inital parameters
       call RANDOM_SEED (SIZE = i)
       allocate(ISEED(i))
       ISEED = 1
       call RANDOM_SEED(PUT=ISEED(1:i))
       
c      == location of particles
       do i=1, NN 
         dtemp = dble(i)/dble(NN)
         j = 1 + dtemp * imax_num_mesh
         X(i) = dtemp

         MC_ind(j) = MC_ind(j) + 1
         MC(j,MC_ind(j)) = i
         Ind_MC(i) = j 
       enddo 

c      == set inital velocity
       do i=1, NN
         call RANDOM_NUMBER(dtemp)
         dtemp = (0.5d0 - dtemp) * 2.0d0 
         V(i) = d_max_velocity * dtemp * 0.9d0
         V_old(i) = 0.0d0
       enddo


!      === Two particles
!       X(1) = 0.0d0
!       X(2) = 1.0d0
!       MC_ind(1) = 1
!       MC_ind(imax_num_mesh) = 1
!       MC(1,1) = 1
!       MC(imax_num_mesh,1) = 2
!       Ind_MC(1) = 1
!       Ind_MC(2) = imax_num_mesh
!       V(1) = d_max_velocity * 0.9d0
!       V(2) = -d_max_velocity * 0.9d0
!       V_old(1) = V(1)
!       V_old(2) = V(2)   
      

 
c      === Start of DEM step ----------------------------
       t1 = omp_get_wtime()

c      === Main loop of time step
       do i=1, MAX_ITER

c        === output current position for files 
!         cbuf1 = "./GNUPLOT/"
!         call Int2Char(i, cbuf2)
!         j = len_trim(cbuf2)
!         cbuf3 = cbuf2(1:j) // ".dat"
!         filename = cbuf1 // cbuf3
!         open(20,file=filename,status='replace')
!         do  j=1, NN
!           write (20, "(F6.4,X,F6.4)") X(j), 1.0d0  
!         enddo
!         close(20)

         t3 = omp_get_wtime()
 
c        === Clear Force
         itotal_num_collision = 0
         do j=1, NN
           F(j) = 0.0d0
         enddo

c        === Move particles
         call MoveParticles(X, V, V_old, 
     &          MC, MC_ind, imax_num_mesh, 
     &          Ind_MC, d_mesh_dist)
c        === End of move particles

c       === multicolor loop
        do i_color=0, 1

c          === Remove meshes that have no particle
           icompact_num_mesh = 1

           do j=1+i_color, imax_num_mesh, 2
             if (MC_ind(j) .ne. 0) then
!              --- there is particle
               MC_compact_ind(icompact_num_mesh) = j
               icompact_num_mesh = icompact_num_mesh + 1
             endif
           enddo

           icompact_num_mesh = icompact_num_mesh - 1
c          === End of removing meshes that have no particle

c          === Calc collisions and forces
           call CollisionAndCalcForce(X, V, V_old, F, 
     &        MC, MC_ind, MC_compact_ind,
     &        imax_num_mesh, icompact_num_mesh,
     &        itotal_num_collision)
c          === End of calc collisions and forces

!           print *, i, X(1), X(2)
!           if (i .eq. 11) stop
  
         enddo
!        == End of multicolor loop

!        == Update Velocity
!$omp parallel do
         do j=1, NN
           if (F(j) .ne. 0.0d0) then
             V(j) = V(j) + F(j)*DTS
           endif
         enddo
!$omp end parallel do


         t4 = omp_get_wtime()
         t_s = t4 - t3

         print *, i, "Time Step =", dble(i)*DTS
         print *, "   Execution Time [sec.] =", t_s, 
     &         ", Number of Collisions =", itotal_num_collision

       enddo
c      === End of Main loop of time step

       t2 = omp_get_wtime()
       t_w =  t2 - t1 
c      === End of DEM step ---------------------------

       print *, "Total execution time [sec.] = ", t_w



       deallocate(X)
       deallocate(V)
       deallocate(V_old)
       deallocate(F)
       deallocate(MC)
       deallocate(MC_ind)
       deallocate(MC_compact_ind)
       deallocate(Ind_MC)
       deallocate(ISEED)

       stop
       end

      subroutine Int2Char (num, char)
      integer num
      character*4 char

      integer i

      char = '    '
      write( char, '( I4 )' ) num

      char = adjustl(char)

      return
      end



       subroutine MoveParticles(X, V, V_old,
     &              MC, MC_ind, imax_num_mesh,
     &              Ind_MC, d_mesh_dist) 
       use omp_lib
       include 'dem.inc'

       double precision X(NN)
       double precision V(NN)
       double precision V_old(NN)
       integer MC(imax_num_mesh, MAX_PARTICLES_IN_MESH)
       integer MC_ind(imax_num_mesh)
       integer imax_num_mesh
       integer Ind_MC(NN)
       double precision d_mesh_dist

       integer i, j, k
       integer kk
       double precision dtemp

!      === Loop for particle NO.

!$omp parallel do private(j,k,kk,dtemp) 
       do i=1, NN

!        == remove particle data
         j = Ind_MC(i)
         do k=1, MC_ind(j)
           if (MC(j, k) .eq. i) then
!            == find the particle
             do kk=k, MC_ind(j)-1
               MC(j, kk) = MC(j, kk+1)
             enddo
             MC_ind(j) = MC_ind(j) - 1
             goto 10
           endif
         enddo
 10      continue

!        == move the particle 
         dtemp = V(i)*DTS 
         if (dabs(dtemp) >= d_mesh_dist) then
           print *, "Particle", i, "exceeds max mesh dist.", d_mesh_dist
           stop 
         endif
         X(i) = X(i) + dtemp

!        == check collison with wall
         if (X(i) <= DPD_HALF) then
!          === collision with left wall
           X(i) = DPD_HALF
           V(i) = -1.0d0 * V(i)
           V_old(i) = 0.0d0
         endif
         if (X(i) >= (1.0d0-DPD_HALF)) then
!          === collision with right wall
           X(i) = 1.0d0-DPD_HALF
           V(i) = -1.0d0 * V(i)
           V_old(i) = 0.0d0
         endif

!        == set new infoamation
         j = 1 + X(i) / d_mesh_dist
         Ind_MC(i) = j
         MC_ind(j) = MC_ind(j) + 1
         if (MC_ind(j) .gt. MAX_PARTICLES_IN_MESH) then
           print *, "Particle ",i, 
     &       ": Overflow of particles is happened in ",j,"th mesh."     
           print *, "Maximum particles per mesh =",MAX_PARTICLES_IN_MESH
           stop
         endif
 
         MC(j, MC_ind(j)) = i
          
       enddo 
!      === End of loop for particle NO. 
!$omp end parallel do

       return
       end


       subroutine CollisionAndCalcForce(X, V, V_old, F, 
     &      MC, MC_ind, MC_compact_ind,
     &      imax_num_mesh, icompact_num_mesh,
     &      itotal_num_collision  )
       use omp_lib
       include 'dem.inc'

       double precision X(NN)
       double precision V(NN)
       double precision V_old(NN)
       double precision F(NN)
       integer MC(imax_num_mesh, MAX_PARTICLES_IN_MESH)
       integer MC_ind(imax_num_mesh)
       integer MC_compact_ind(imax_num_mesh)
       integer imax_num_mesh
       integer icompact_num_mesh
       integer itotal_num_collision 

       integer i, j, k
       integer i_target
       integer i_particle, k_particle



!$omp parallel do private(i_target,j,i_particle,k,k_particle) 
!$omp&           reduction(+:itotal_num_collision)

!      === Loop for mesh NO. 
       do i=1, icompact_num_mesh

!        == get NO. of the target mesh  
         i_target = MC_compact_ind(i)

!        == Particles in i-th mesh
         do j=1, MC_ind(i_target)

!          --- tartget particle NO.
           i_particle = MC(i_target, j)

c          === for particles in i-th mesh
           do k=j+1, MC_ind(i_target)
             k_particle = MC(i_target, k)
             call CollisionDetect
     &           (i_particle, k_particle, X, V, V_old, F, 
     &            itotal_num_collision)
           enddo

c          === for particles in (i+1)-th mesh
           if ((i_target+1) <= imax_num_mesh) then
             do k=1, MC_ind(i_target+1)
               k_particle = MC(i_target+1, k)
               call CollisionDetect
     &           (i_particle, k_particle, X, V, V_old, F, 
     &            itotal_num_collision)
             enddo
           endif

         enddo

       enddo
!$omp end parallel do

!      === End of loop for mesh NO. 

       return
       end


       subroutine CollisionDetect
     &         (i_particle, k_particle, X, V, V_old, F, 
     &          itotal_num_collision)
       use omp_lib
       include 'dem.inc'

       integer i_particle, k_particle
       double precision X(NN)
       double precision V(NN)
       double precision V_old(NN)
       double precision F(NN)
       integer itotal_num_collision

       double precision dtemp

       dtemp = dabs(X(i_particle) - X(k_particle))
       if (dtemp <= DPD ) then
!        == make collision
         itotal_num_collision = itotal_num_collision + 1

!         print *, i_particle, k_particle, V(i_particle),  V(k_particle)

!        == Force calculation
         dtemp = PARTICLE_WEIGHT * 
     &           ( V(i_particle) - V(k_particle) ) / DTS   


         F(i_particle) = F(i_particle) - dtemp 
         F(k_particle) = F(k_particle) + dtemp 

!         print *, i_particle, k_particle, F(i_particle),  F(k_particle)

         V_old(i_particle) = 0.0d0
         V_old(k_particle) = 0.0d0
         

       endif

       return
       end 
