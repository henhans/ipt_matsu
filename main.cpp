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
  int Nw=1024*16; //number of tau and omega points should be power of 2
  double U=2.9; // coulomb potential
  double T=0.01; // temperature
  double t=0.5;//hoping constant
  int in=1; // initialize option 0 for metal, 1 for insulator, 2 from file
  double mix=1; // mixing between previous and current iterations

  int Nit=1000; // number of DMFT iteration

  mpt Mpt( Nw, U, T, t);
  Mpt.init_gw(in);
  Mpt.printG0w("Gw_init.dat");
  //Mpt.InvFourier(Mpt.G0w,Mpt.G0t);
  //Mpt.printG0t("Gt_init.dat");
  //Mpt.Fourier(Mpt.G0t,Mpt.G0w);
  //Mpt.printG0w("Gw2.dat");
  double old_diff=1e10;

  for( int it=0; it< Nit; it++)   
  { 
    Mpt.im_solver();
    double diff = Mpt.get_diff(); //calculate difference
    //cout << diff <<endl;

    // mixing solutions
    if (diff<1e-7) 
      break;
    Mpt.update_G(mix);
    printf("it= %i \t diff= %f \t G(iw_0)= %f \t S(iw_0)= %f \n", it, diff, imag(Mpt.Gw[0]), imag(Mpt.Sw[0]));

  }

  Mpt.printGw("Gw_final.dat"); 
  Mpt.printSw("Sw_final.dat");
  Mpt.InvFourier(Mpt.Gw,Mpt.Gt);
  Mpt.printGt("Gt_final.dat");
  Mpt.printSt("St_final.dat");

  return 0;
}
