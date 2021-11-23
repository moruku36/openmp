#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>


//!      == Number of Particles
//!   -- for test
//#define   N   2
//#define   N   20
#define   N    10000000

//!      == Particle Diameter
//!        -- for test
//#define   DPD  1.0e-2
#define   DPD  1.0e-7

#define   DPD_HALF   DPD/2.0

//!      == Particle Weight
//        -- for test
#define	  PARTICLE_WEIGHT  1.0e-2

//!      == Time Step
//!        -- for test
#define   DTS  1.0e-3

//!      === Maximum number of particles in collision detection mesh
#define   MAX_PARTICLES_IN_MESH  20

//!       === for test
//#define   MAX_ITER  300
#define   MAX_ITER  30


/* Please define the matrices in here */
static double X[N];
static double V[N];
static double V_old[N];
static double F[N];

static unsigned long *MC;
static unsigned long *MC_ind;
static unsigned long Ind_MC[N];
static unsigned long *MC_compact_ind;


void MoveParticles(double X[N], double V[N], double V_old[N],
		   unsigned long *MC, unsigned long *MC_ind, unsigned long imax_num_mesh,
		   unsigned long Ind_MC[N], double d_mesh_dist);
unsigned long CollisionAndCalcForce(double X[N], 
         double V[N], double V_old[N], double F[N],
	 unsigned long *MC, unsigned long *MC_ind, unsigned long *MC_compact_ind,
	 unsigned long imax_num_mesh, unsigned long icompact_num_mesh);

unsigned long CollisionDetect(unsigned long i_particle, unsigned long k_particle, 
       double X[N], double V[N], double V_old[N], double F[N]);


int main(int argc, char* argv[]) {

     double t0, t1, t2, t_w;
     double t3, t4, t_s;
     double d_mesh_dist, d_max_velocity;
     double dtemp;
     double dc_inv;

     unsigned long i, j, k, n;
     unsigned long imax_num_mesh, icompact_num_mesh;
     unsigned long itotal_num_collision;
     unsigned long itemp;

     char  cbuf[5];
     char  *filename;
     FILE* fp;


     dtemp = 1.0e0/DPD;
     imax_num_mesh = (unsigned long)dtemp + 1;
     d_mesh_dist = DPD;
     d_max_velocity = d_mesh_dist / DTS;

     printf("Number of Particles = %d \n", N);
     printf("Particle Diameter = %lf \n", DPD);
     printf("Particle Half Diameter = %lf \n", DPD_HALF);

     printf("Time Step = %lf \n", DTS);
     printf("Particle weight = %lf \n", PARTICLE_WEIGHT);
     printf("imax_num_mesh = %d \n", imax_num_mesh);
     printf("d_mesh_dist = %lf \n", d_mesh_dist);
     printf("d_max_velocity = %lf \n", d_max_velocity);

     if (imax_num_mesh < N ) {
       printf(" !!The number of particles needs to reduce. \n");
       exit(1);
     }

     MC_compact_ind = (unsigned long *)calloc(imax_num_mesh,sizeof(unsigned long));
     MC = (unsigned long *)calloc(imax_num_mesh*MAX_PARTICLES_IN_MESH,sizeof(unsigned long));
     MC_ind = (unsigned long *)calloc(imax_num_mesh,sizeof(unsigned long));

     for (i=0; i<imax_num_mesh; i++) {
       MC_ind[i] = 0;
       MC_compact_ind[i] = 0;
     }


     //=== make inital parameters
     srand(1);
     dc_inv = 1.0/(double)RAND_MAX;

     //== location of particles
     for (i=0; i<N; i++) {
       dtemp = (double)i/(double)N;
       j = (unsigned long)(dtemp * (double)imax_num_mesh);
       X[i] = dtemp;

       MC_ind[j] ++;
       MC[MAX_PARTICLES_IN_MESH*j + MC_ind[j]] = i;

       Ind_MC[i] = j;
     
     }

     // == set inital velocity
     for (i=0; i<N; i++) {
       dtemp = rand() * dc_inv;
       dtemp = (0.5 - dtemp) * 2.0;
       V[i] = d_max_velocity * dtemp * 0.9;
       V_old[i] = 0.0;
     }
     
     //=== Two particles
     /*
     X[0] = 0.0;
     X[1] = 1.0;
     MC_ind[0] = 1;
     MC_ind[imax_num_mesh-1] = 1;
     MC[MAX_PARTICLES_IN_MESH*0 +0] = 0;
     MC[MAX_PARTICLES_IN_MESH*(imax_num_mesh-1) +0] = 1;
     Ind_MC[0] = 0;
     Ind_MC[1] = imax_num_mesh-1;
     V[0] = d_max_velocity * 0.9;
     V[1] = -d_max_velocity * 0.9;
     V_old[0] = V[0];
     V_old[1] = V[1];
     */

     /* Start of DEM step ----------------------------*/
     t1 = omp_get_wtime();

     //=== Main loop of time step
     for (i=1; i<=MAX_ITER; i++) {


       //=== output current position for files
       /*       
       filename = (char *)calloc(25, sizeof(char)); 
       strncat(filename,"./GNUPLOT/", 10);
       sprintf(cbuf, "%lu", i);
       strncat(filename,cbuf, 4);
       strncat(filename,".dat",4);

       if ((fp = fopen(filename, "w")) == NULL) {
          printf("fopen err! \n");
          exit(1);
       }

       for (j=0; j<N; j++) {
          fprintf(fp, "%lf  1.000 \n",X[j]);
       }
  
       fclose(fp);
       free(filename);
       */      

     t3 = omp_get_wtime();

     //=== Clear Force
     itotal_num_collision = 0;
     for (j=0; j<N; j++) {
       F[j] = 0.0;
     }

     //=== Move particles
     MoveParticles(X, V, V_old, 
                   MC, MC_ind, imax_num_mesh,
		   Ind_MC, d_mesh_dist);
     //=== End of move particles


     //=== Remove meshes that have no particle
     itemp = 1;
     for (j=0; j<imax_num_mesh; j++) {

       if (MC_ind[j] != 0) {
	 //--- there is particle
	 MC_compact_ind[itemp-1] = j;
         itemp ++;
       }
     }

     icompact_num_mesh = itemp - 1;
     //=== End of removing meshes that have no particle

     //=== Calc collisions and forces
     itotal_num_collision = 
         CollisionAndCalcForce(X, V, V_old, F,
	   MC, MC_ind, MC_compact_ind,
	   imax_num_mesh, icompact_num_mesh);
     //=== End of calc collisions and forces


     //== Update Velocity
#pragma omp parallel for
     for (j=0; j<N; j++) {
       if (F[j] != 0.0) {
         V[j] = V[j]+F[j]*DTS;
       }
     }



     t4 = omp_get_wtime();
     t_s = t4 - t3;

     printf("%lu Time Step = %lf \n", i, (double)i*DTS);
     printf("   Execution Time [sec.] = %e, Number of Collisions = %lu \n",
   	      t_s, itotal_num_collision);

     } //=== End of main loop of time step

     t2 = omp_get_wtime();
     t_w =  t2 - t1; 
     /* End of Main Loop --------------------------- */

     printf("Total execution time = %e [sec.] \n",t_w);

     
     free(MC_compact_ind);
     free(MC);
     free(MC_ind);
     

     return 0;
}

void MoveParticles(double X[N], double V[N], double V_old[N],
      unsigned long *MC, unsigned long *MC_ind, unsigned long imax_num_mesh,
      unsigned long Ind_MC[N], double d_mesh_dist)
{

     unsigned long i, j, k;
     unsigned long kk;
     double dtemp;


     // === Loop for particle NO.
#pragma omp parallel for private(j,k,kk,dtemp)
     for (i=0; i<N; i++) {


       //== remove particle data
       j = Ind_MC[i];

       for (k=0; k<MC_ind[j]; k++) {
         if ( MC[MAX_PARTICLES_IN_MESH*j + k] == i) {
	   //== find the particle
	   for (kk=k; kk<MC_ind[j]-1; kk++) {
             MC[MAX_PARTICLES_IN_MESH*j + kk] = MC[MAX_PARTICLES_IN_MESH*j + kk+1];
           }
	   MC_ind[j] --;

           goto outloop;
         }
       }

   outloop:  


      // == move the particle
      dtemp = V[i]*DTS;
      if (fabs(dtemp) >= d_mesh_dist) {
         printf("Particle %lu exceeds max mesh dist. %lf", i, d_mesh_dist);
         exit(1);
      }

      X[i] = X[i] + dtemp;

      //== check collison with wall
      if (X[i] <= DPD_HALF) {
	 //=== collision with left wall
         X[i] = DPD_HALF;
	 V[i] = -1.0 * V[i];
         V_old[i] = 0.0;
       }
       if (X[i] >= (1.0-DPD_HALF)) {
         //=== collision with right wall
         X[i] = 1.0-DPD_HALF;
         V[i] = -1.0 * V[i];
         V_old[i] = 0.0;
       }
       // == set new infoamation

       j = X[i] / d_mesh_dist;
       Ind_MC[i] = j;

       MC_ind[j] ++;


       if (MC_ind[j] > MAX_PARTICLES_IN_MESH) {
         printf("Particle %lu :Overflow of particles　is happened in %lu th mesh. %lu \n",i, j, MC_ind[j]);
	 printf("Maximum particles per mesh = %lu \n", MAX_PARTICLES_IN_MESH);
         exit(1);
       }
       MC[MAX_PARTICLES_IN_MESH*j + MC_ind[j]-1] = i;

     } //! === End of loop for particle NO.

}

unsigned long CollisionAndCalcForce(double X[N], 
         double V[N], double V_old[N], double F[N],
	 unsigned long *MC, unsigned long *MC_ind, unsigned long *MC_compact_ind,
	 unsigned long imax_num_mesh, unsigned long icompact_num_mesh)
{

    unsigned long i, j, k;
    unsigned long i_target;
    unsigned long i_particle, k_particle;
    unsigned long inum_collision;



    inum_collision = 0;

#pragma omp parallel for private(i_target,j,i_particle,k_particle,k) reduction(+:inum_collision)
    //=== Loop for mesh NO.
    for(i=0; i<icompact_num_mesh; i++) {

      //== get NO. of the target mesh
      i_target = MC_compact_ind[i];
        

      //== Particles in i-th mesh
      for (j=0; j < MC_ind[i_target]; j++) {

	//--- tartget particle NO.
	i_particle = MC[MAX_PARTICLES_IN_MESH*i_target + j];


	//=== for particles in i-th mesh
        for (k=j+1; k< MC_ind[i_target]; k++) {
	   k_particle = MC[MAX_PARTICLES_IN_MESH*i_target + k];

           inum_collision += 
               CollisionDetect(i_particle, k_particle, X, V, V_old, F);
	}

	//=== for particles in (i+1)-th mesh
	if ((i_target+1) < imax_num_mesh) {
	  for (k=0; k< MC_ind[i_target+1]; k++) {
 
	    k_particle = MC[MAX_PARTICLES_IN_MESH*(i_target+1) + k];
            inum_collision +=
              CollisionDetect(i_particle, k_particle, X, V, V_old, F);
	  }
	}

      } 
 
    }//=== End of loop for mesh NO.



    return inum_collision;

}    


unsigned long CollisionDetect(unsigned long i_particle, unsigned long k_particle, 
       double X[N], double V[N], double V_old[N], double F[N])
{

    double dtemp;
    unsigned long inum_collision;


    inum_collision = 0;
 
    dtemp = fabs(X[i_particle] - X[k_particle]);
    if (dtemp <= DPD ) {

      //== make collision
      inum_collision ++;

      //printf("%lu %lu %lf %lf",i_particle,k_particle,V[i_particle],V[k_particle]);

#pragma omp critical

      //== Force calculation
      dtemp = PARTICLE_WEIGHT * ( V[i_particle] - V[k_particle] ) / DTS;
      F[i_particle] = F[i_particle] - dtemp;
      F[k_particle] = F[k_particle] + dtemp;

      //printf("%lu %lu %lf %lf",i_particle,k_particle,F[i_particle],F[k_particle]);
      V_old[i_particle] = 0.0;
      V_old[k_particle] = 0.0;

#pragma omp end critical

      //printf("%lu %lu %lf %lf",i_particle,k_particle,V[i_particle],V[k_particle]);
    }

    return inum_collision;

} 
