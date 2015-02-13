/* DMFT with IPT solver on real axis based on Dr. Haule's lecture materials.   
                                                   07/16/2014 Tsung-Han Lee */
#include "mpt.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <omp.h>

int main(int argc, char* argv[])
{
  int Nw=8192*2; //number of tau and omega points should be power of 2
  double U=2.2; // coulomb potential
  double T=0.02; // temperature
  double t=0.5;//hoping constant

  int Nit=400; // number of DMFT iteration
  double mix=1.0; // mixing constant for old and new solution

  mpt Mpt( Nw, U, T, t);
  Mpt.init_gw();
  //for(int i=0; i<Nw; i++) cout <<imag(Mpt.omega[i])<<"\t"<<real(Mpt.G0w[i]) <<"\t"<< imag(Mpt.G0w[i]) <<endl;
  Mpt.printG0w("Gw.dat");
  Mpt.InvFourier(Mpt.G0w,Mpt.G0t);
  Mpt.printG0t("Gt.dat");
  //for(int i=0; i<Nw; i+=1) cout<< Mpt.tau[i] << "\t" << real(Mpt.G0t[i])<<"\t"<<imag(Mpt.G0t[i]) <<endl;
  Mpt.Fourier(Mpt.G0t,Mpt.G0w);
  Mpt.printG0w("Gw2.dat");
  //for(int i=0; i<Nw; i++) cout <<imag(Mpt.omega[i])<<"\t"<< real(Mpt.G0w[i]) <<"\t"<< imag(Mpt.G0w[i]) <<endl;
  //double old_diff=1e10;
  //cout << old_diff << endl;

  /*for( int it=0; it< Nit; it++)   
  { 
    Pt.get_sig();
    double diff = Pt.get_diff(); //calculate difference
    //cout << diff <<endl;

    // mixing solutions
    if (diff<1e-4) 
      break;
    if (diff>old_diff)
        mix = mix/1.;//2;
    old_diff = diff;
    Pt.mix_sig(mix);
    printf("it= %i \t diff= %f \t mix= %f \n", it, diff, mix);

    Pt.get_G_G0(); // calculate Green's functions
  }*/

  //Mpt.printdata("ipt.dat"); 

  return 0;
}
