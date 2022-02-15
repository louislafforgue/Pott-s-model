// compile with: 
// path to the file: /mnt/c/Users/louis/Documents/"master innsbruck"/"[lecture] 3rd semester"/statistical_physics/potts_model
//   g++ Ising_Metropolis.cpp lattice.cpp nrutil.cpp -o test -O2 -w
/* We onside a 2D lattice subject to the hamiltonian: 
H=-J Sum_{<ij>} d(s_i,s_j)

*/ 
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include "potts.h"
#include "lattice.h"
#include "nrutil.h"

int main (){
  // set parameters
  L = 100;
  N = L*L;
  int q = 6; // number of spin degree of freedom
  double h=0.0; 
  double T = 0.1;
  int nmcstep = 10000;
  
  //printf("critical temperature = %.4f \n", 2./log(1.+sqrt(2.)));
  //double ma = pow((1.- pow(sinh(2./T),-4.)), 0.125);
  
  // array of length N
  s = dvector(0,N-1);
  s_initial_values = dvector(0,N-1);
  s0 = dvector(0,N-1);
  s1 = dvector(0,N-1);
  s2 = dvector(0,N-1);
  s3 = dvector(0,N-1);
  s4 = dvector(0,N-1);
  s5 = dvector(0,N-1);
  
  // random configuration
  for (int i =0; i<N; i++){
	int xa = rand()%q; // generate a random spin number between -q and q
  double x= xa*1.0;

  //comment for random configuration
  //x = 0.0; 

  //printf("%1.f \n", xa);
  s_initial_values[i]=x; // set the initial value to each lattice site
  //printf("%i \n", s_initial_values[i]);
  }

  //create the lattice
  class Lattice *lattice;
  lattice = new Lattice(N, q, h, s, s_initial_values);
  
  //save the initial configuration
  FILE * out;
  out = fopen("output1.dat","w");
  
  for (int k =0; k<N; k++) {
    s0[k]=s[k];
    //printf("%.1f \n",s0[k]);
    fprintf(out,"%1.f %.1f \n",k,s0[k]);
  }
  fclose(out);

  // open output file for the energy and magnetization
  out = fopen("output3.dat","w");
  double m =  lattice->get_magnetization(s);
  //printf("%.4f \n",m);
  fprintf(out,"%d %.1f %.4f \n",0,lattice->get_energy(s), m);
 

  
  // Monte Carlo - Metropolis rule - equilibration
  int l=0;
  int j=40;
  int i_mc=0;
  while (i_mc<3*j){
    i_mc++;
    lattice->Metropolis (s, T, nmcstep);
    m =  lattice->get_magnetization(s);
    fprintf(out,"%d %.1f %.4f  \n",
                        (i_mc+1)*nmcstep,lattice->get_energy(s), m);
    //printf("%i \n",i_mc); 
    //copy some configurations for the debbugging
    if (i_mc==1) {
	  for (int k =0; k<N; k++) {s1[k]=s[k];}
	} else if (i_mc==5) {
	  for (int k =0; k<N; k++) {s2[k]=s[k];}
	} else if (i_mc==30) {
	  for (int k =0; k<N; k++) {s3[k]=s[k];}
	} else if (i_mc==60) {
	  for (int k =0; k<N; k++) {s4[k]=s[k];}
	} else if (i_mc==100) {
	  for (int k =0; k<N; k++) {s5[k]=s[k];}
	}
	// if we reach the theoritical value we stop the simulation
/*
	if (m>ma && l==0){
      printf("relaxation time = %d \n", (i+1)*nmcstep);
      j=i;
      l=1;
	}
  */
  }
  // printf("%i",i_mc); 
  fclose(out);

// copy in output 4 the different configurations 
  out = fopen("output4.dat","w");
  for (int i =0; i<N; i++) {
    fprintf(out,"%.1f %.1f %.1f %.1f %.1f %.1f\n",
                        s0[i],s1[i],s2[i],s3[i],s4[i],s5[i]);
  }
  fclose(out);
  
  // Magnetization autocorrelation function
  //~ nmcstep=1000;
  //~ int nsave = 500000;
  //~ mvector = dvector(0,nsave-1);
  //~ for (int k=0; k<nsave; k++){
    //~ lattice->Metropolis (s, T, nmcstep);
    //~ mvector[k] =  lattice->get_magnetization(s);
  //~ }
  //~ out = fopen("output5.dat","w");
  //~ for (int k =0; k<nsave-10; k++){
	//~ fprintf(out,"%d %f\n",k*nmcstep, autocorrelation (mvector,nsave,k));
  //~ }
  //~ fclose(out);
  
  
  return 0;
}

double mean (double x[],int NN){
  // define the mean function
  double xm=0.0;
  for (int i =0; i<NN; i++){
    xm += x[i];
  }
  xm = xm/NN;
  return xm;
}

double variance (double x[],int NN){
  // define the variance function
  double xm=0.0;
  for (int i =0; i<NN; i++){
    xm += x[i]*x[i];
  }
  xm = xm/NN - mean(x,NN)*mean(x,NN);
  return xm;
}

double autocorrelation (double x[], int NN, int tau){
  // define the autocorrelation function
  // input:
  //   x = 1D array of data
  //   tau = time lapse
  // output:
  //   c_n = autocorrelation of xs at time lapse n
  double mu1 = mean (x,NN);
  double sigma2 = variance (x,NN);
  double c = 0.0;
  for (int i =0; i<NN-tau; i++){
    c+=(x[i]-mu1)*(x[i+tau]-mu1);
  }
  c = c/(NN-tau);
  c = c/sigma2;
  return c;	
}
