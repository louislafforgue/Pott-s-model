#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <complex>
#include "lattice.h"
#include "nrutil.h"

template <class T> class complex;

#define PI 3.1415926


Lattice::Lattice(int N, double s[], double s_initial_values[]){
  // call constructur make a lattice with periodic boundary,
  // spins initially setted equal to the array "s_initial_values"
  L = sqrt(N);
  for (int i =0; i<N; i++){
    s[i] = s_initial_values[i];
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

double Lattice::get_energy(double s[], double J, double h){
  // return the energy E of a given configuration
  int N = L*L;
  double E_int = 0.0;
  double E_ext = 0.0;
  for (int i =0; i<N; i++){
    if(s[get_right(i)] == s[i]){
	  E_int ++;
	}
	if(s[get_below(i)] == s[i]){
	  E_int ++;
	}
	E_ext += s[i];
  }
  E_int *= -J;
  E_ext *= -h;
  E_int += E_ext;
  return E_int;
}

std::complex<double> Lattice::get_magnetization(double s[], int q){
  // return the magnetization m of a given configuration
  int N = L*L;
  std::complex<double> m =(0.0,0.0);
  std::complex<double> np[q];
  for (int p = 0; p<q; p++){
	np[p] = (0.0,0.0);
  }
  for (int i=0; i<N; i++){
	np[int(s[i])] += std::complex<double>(1.0,0.0);
  }
  for (int p = 0; p<q; p++){
	np[p] /= std::complex<double> (N*1.0, 0.0);
  }
  
  for (int p =0; p<q; p++){
	m += exp(std::complex<double>(0.0,(2*PI*(p+1)/q))) * np[p];

  
  }
  return m;
}

void Lattice::Metropolis (double s[], double T, int q, double J, double h, int nmcstep){
  // perform nmcstep Monte Carlo time steps 
  // using the Metropolis Algorithm 
  // T = temperature
  int N = L*L;
  double beta = 1./T;
  double * precalc1;
  double * precalc2;
  precalc1 = dvector(-4,4);
  precalc2 = dvector(-(q-1),(q-1));
  for (int i=-4; i<=4; i++){
    precalc1[i] = exp(J*beta*i);
  }
  for (int i=-(q-1); i<q; i++){
    precalc2[i] = exp(h*beta*i);
  }
  for (int t=0; t<nmcstep; t++){
    int i = rand() % N;
    int c = rand() % q;
    while (s[i] == c){
	  c = rand() % q;
    }    
    
    double E_int_in = 0;
    double E_int_fin = 0;
    
    if(s[get_below(i)] == s[i]){
	  E_int_in ++;
    }
    if(s[get_below(i)] == c){
	  E_int_fin ++;
    }
    if(s[get_above(i)] == s[i]){
	  E_int_in ++;
    }
    if(s[get_above(i)] == c){
	  E_int_fin ++;
    }
    if(s[get_right(i)] == s[i]){
	  E_int_in ++;
    }
    if(s[get_right(i)] == c){
	  E_int_fin ++;
    }
    if(s[get_left(i)] == s[i]){
	  E_int_in ++;
    }
    if(s[get_left(i)] == c){
	  E_int_fin ++;
    }
    
    
    int k = int(E_int_fin-E_int_in);
    int l = int(c-s[i]);
	if ((precalc1[k]*precalc2[l]) > rand()/((double) RAND_MAX)) {
      s[i] = c;
    }
  }
}
