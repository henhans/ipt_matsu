/* class for perturbation calculations 07/16/2014 Tsung-Han Lee*/

#ifndef MPT_
#define MPT_
#include <complex>

using namespace std;

class mpt
{
  public:
  mpt(int Nw,double U,double T, double t);

  void init_gw(int in); // initiate green function in t space
  void im_solver(); // impurity solver
  void update_G(double mix); //update Green's function
  double get_diff(); // calculate self energy difference betwen consecutive iterations
  //void get_G_G0(); // calculate Green's function G and G0
  void printG0w(const char* filename); // print out Gw
  void printG0t(const char* filename); // print out Gt
  void printGw(const char* filename); // print out Gw
  void printSw(const char* filename); // print out Sw
  void printSt(const char* filename); // print out St
  void Fourier(complex<double> Gt[],complex<double> Gw[]);// fourier transofrm
  void InvFourier(complex<double> Gw[],complex<double> Gt[]);// inverse fourier transofrm
  

  complex<double>* omega;
  double* tau;
  //complex<double>* sigma;
  //complex<double>* sigma_new;
  complex<double>* G0t;
  complex<double>* G0w;
  complex<double>* Gw;
  complex<double>* Sw;
  complex<double>* Sw_old;
  complex<double>* St;
  complex<double>* exp1;
  complex<double>* exp2;
  double* corr_wtot;//correction term
  double* A_n;//attenuation factor
  //complex<double>* G;

  //void init_funcs();
  
  private:
  int _Nw; //number of omega point
  double _U; // coulomb potential
  double _T; // temperature
  double _Beta;// inverse temperature
  double _t;// hoping constant
  double _domega;// the difference between omega grids
  double _dtau;// the difference between tau grids

  void FFT(double data[], unsigned long N, int isign);//Fast fourier routine fomr numerical recipie
  double cal_corr_wtot(double tau); // calculate high frequency correction
};

#endif
