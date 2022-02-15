#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include "lattice.h"
#include "nrutil.h"

Lattice::Lattice(int N, int Q, double H, double s[], double s_initial_values[]){
  // call constructur make a lattice with periodic boundary,
  // spins initially setted equal to the array "s_initial_values"
  L = sqrt(N);
  q=Q; //number of different spins
  h=H; 
  for (int i =0; i<N; i++){

    s[i] = s_initial_values[i];
    //printf("%.1f \n", s[i]);
  }
}

int Lattice::get_spin_index(int i, int j){
  // i (j) = row (column) index from 0 to L-1
  // Output: the spin index corresponding to (i,j)
  int k = i*L + j;
  return k;
}

double Lattice::get_spin_value(double s[], int i, int j){
  // i (j) = row (column) index from 0 to L-1
  // Output: the spin value corresponding to (i,j)
  int k = i*L + j;
  return s[k];
}

int Lattice::get_left(int i){
  // return the index of the site at the left of site i
  int j=-1;
  if (i%L == 0){
    j = i+L-1;
  } else {
    j = i-1;
  }
  return j;
}

int Lattice::get_right(int i){
  // return the index of the site at the right of site i
  int j=-1;
  if ((i+1)%L == 0){
    j = i+1-L;
  } else {
    j = i+1;
  }
  return j;
}

int Lattice::get_above(int i){
  // return the index of the site above of site i
  int j=-1;
  if (i<L){
    j = i + (L-1)*L;
  } else {
    j = i-L;
  }
  return j;
}

int Lattice::get_below(int i){
  // return the index of the site below of site i
  int j=-1;
  if (i>=L*(L-1)){
    j = (i+L) % (L*L);
  } else {
    j = i+L;
  }
  return j;
}

double Lattice::get_energy(double s[]){
  // return the energy E of a given configuration
  int N = L*L;
  double E = 0.0;
  for (int i =0; i<N; i++){
    E -= int(s[i]==s[get_right(i)]) + int(s[i]==s[get_below(i)]);
  } 
  return E;
}

double Lattice::get_magnetization(double s[]){
  // return the magnetization m of a given configuration
  int N = L*L;
  double m =0.0;
  for (int i =0; i<N; i++){
    m += s[i];
  }
  m = fabs(m)/N;
  return m;
}

void Lattice::Metropolis (double s[], double T, int nmcstep){

  // perform nmcstep Monte Carlo time steps 
  // using the Metropolis Algorithm 
  // T = temperature
  int N = L*L;
  double beta = 1./T;
  double * precalc;
  //how many value can energy difference be? 
  precalc = dvector(-q,q);
  for (int i=-q; i<=q; i++){
    precalc[i] = exp(-beta*i);
  }
  for (int t=0; t<nmcstep; t++){
    //select a random spin
    int i = rand() % N;
    double E_in=0;
    //printf("%.1f %.1f %.1f %.1f \n", get_right(i), get_left(i), get_below(i),get_above(i)); 
    if(s[i]==s[get_right(i)]){E_in-=1;}
    if(s[i]==s[get_left(i)]){E_in-=1;}
    if(s[i]==s[get_below(i)]){E_in-=1;}
    if(s[i]==s[get_above(i)]){E_in-=1;}
    E_in-=s[i]*h;
    //select a random orientation
    double E_fin=0; 
    double xa = rand()%q;
    if(xa==s[get_right(i)]){E_fin-=1;}
    if(xa==s[get_left(i)]){E_fin-=1;}
    if(xa==s[get_below(i)]){E_fin-=1;}
    if(xa==s[get_above(i)]){E_fin-=1;}
    E_fin-=h*xa;

    //calculate the energy difference
    int j = int((E_fin-E_in));
    //printf("%.1f %.1f %.1f \n", xa, E_fin, E_in); 
    //if smaller than 0, we accept the move
    if (j<=0) {
      s[i] = xa;
	}
  //if not we accept the move with a certain probability 
  else {
	  if (precalc[j] > rand()/((double) RAND_MAX)) {
        s[i] = xa;
      }
	}
  }
}
