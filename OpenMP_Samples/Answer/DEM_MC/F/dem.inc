!      == Number of Particles
       integer  NN
!         -- for test
!         parameter (NN=20)
        parameter (NN=10000000)

!      == Particle Diameter 
       double precision DPD
!        -- for test
!        parameter (DPD=1.0d-2)
       parameter (DPD=1.0d-7)

       double precision DPD_HALF
       parameter (DPD_HALF=DPD/2.0d0)

!      == Particle Weight
       double precision PARTICLE_WEIGHT
!       --- for test
       parameter (PARTICLE_WEIGHT=1.0d-2)

!      == Time Step
       double precision DTS
!        -- for test 
        parameter (DTS=1.0d-3)

!      === Maximum number of particles in collision detection mesh
       integer MAX_PARTICLES_IN_MESH
       parameter (MAX_PARTICLES_IN_MESH=20)