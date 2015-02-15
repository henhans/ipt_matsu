#include "mpt.h"
#include "routines.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

#ifdef _OMP
#include <omp.h>
#endif

using namespace std;

mpt::mpt( int Nw, double U, double T, double t )
{
  _Nw = Nw;
  _U = U;
  _T = T;
  _Beta=1/T;
  _domega=2*pi/_Beta;
  _dtau=_Beta/_Nw;
  _t=t;

  omega = new complex<double> [Nw];
  tau = new double [Nw];
  //sigma = new complex<double> [Nw];
  //sigma_new = new complex<double> [Nw];
  G0t = new complex<double> [Nw];
  G0w = new complex<double> [Nw];
  Gw = new complex<double> [Nw];
  Sw  = new complex<double> [Nw];
  Sw_old = new complex<double> [Nw];
  St  = new complex<double> [Nw];
  exp1= new complex<double> [Nw];
  exp2= new complex<double> [Nw];
  A_n= new double [Nw];
  corr_wtot= new double [Nw];

  for (int i=0; i<_Nw; i++)
  {
     //sigma[i] =0;
     G0t[i] = 0.0;
     G0w[i] = 0.0;
     Gw[i] = 0.0;
     Sw[i]  = 0.0;
     Sw_old[i] = 0.0;
     St[i]  = 0.0;
     //omega[i] = double((2*i+1)*pi/_Beta);
     tau[i] = double( (i+0.5)*_Beta/_Nw );
     exp1[i] = exp( complex<double> (0.0,-(i+0.0)*pi/Nw) );
     exp2[i] = exp( complex<double>( 0.0,-(i+0.5)*pi/_Nw) );
     //cout<< omega[i]<<"\t"<<exphalf[i]<<endl;
  }
  // initial the omega (0,1,...,N/2-1,-N/2,-N/2+1,....-1)
  for (int i=0; i<_Nw/2; i++){
     omega[i] = complex<double>(0.0,(2*i+1)*pi/_Beta);
     omega[_Nw/2+i]=complex<double>(0.0,2*(_Nw/2+i-_Nw)+1.)*pi/_Beta;
     A_n[i] = sqr( sin((i+0.5)*pi/_Nw ) / ( (i+0.5)*pi/_Nw) );
     A_n[_Nw/2+i] = sqr( sin((_Nw/2+i-_Nw + 0.5)*pi/_Nw ) / ( (_Nw/2+i-_Nw+0.5)*pi/_Nw) );
  }
  // 
  for (int i=0; i<_Nw; i++){
    corr_wtot[i] = -0.5 * sign(tau[i]) + 2*_T*cal_corr_wtot(tau[i]);
    //cout<< imag(omega[i])<<"\t"<<corr_wtot[i]<<endl;
  }
  //omega[_Nw/2]=complex<double>(1.0,0.0);
  //for (int i=0; i<+Nw; i++) cout<<imag(omega[i])<<endl;
  //cout << Init.omega[0]<< endl;
}

double mpt::cal_corr_wtot(double tau) {
  double sum = 0.0;

  for (int i=0; i <_Nw/2; i++) {
    sum += sin( imag(omega[i])*tau ) / imag(omega[i]);
  }
 
  return sum;
}

void mpt::init_gw(int in) {
  double d=4;
  //double V=0.7;
  double n=0.5;
  double mu=_U/2;

  if(in==0) {
     for (int i=0; i<_Nw; i++){
        G0w[i]= /*V*V * 2/(d*d)*/ _t*_t * ( omega[i] - 1.0j*sign(imag(omega[i])) *
                  sqrt( - omega[i]*omega[i] + d*d) );//1/(omega[i]);
    }
  }
  else if(in==1) {
     for (int i=0; i<_Nw; i++) {
       //G0w[i]=1/omega[i];
       // initial for insulator in atomic limit
       G0w[i]=( 1.0 - n ) / (omega[i] + mu ) + n / (omega[i] + mu - _U );
       //cout << imag(omega[i]) <<"\t"<<imag(G0w[i])<< endl;
     }
  }
  else if(in==2) {
     FILE *f;
     f = fopen("Gw_final.dat", "r");
     double temp1, temp2, temp3;
     for(int i=0; i<_Nw; i++) {
       fscanf(f,"%le %le %le",&temp1,&temp2,&temp3);
       G0w[i]=complex<double> (temp2 , temp3);
       //cout << temp1<<endl;
     }
     cout<<omega[0]<<"\t"<<G0w[0]<<endl;
     fclose(f);
  }

}

void mpt::im_solver()
{
  // calculate the green's function from the input bath
  for (int i=0; i<_Nw; i++) {
     G0w[i]=1.0/(omega[i]-_t*_t*G0w[i]);
     Sw_old[i]= Sw[i];
  }

  // fourier transform to tau space
  InvFourier(G0w,G0t);
  //cout << "Gw to Gt done" << endl;
  // calculate Sig(tau) using second order perturbation
  double xjump=0.25; //the jump in selfenergy
  for (int i=0; i<_Nw; i++) {
     St[i]=complex<double> ( real(G0t[i]*G0t[i]*G0t[_Nw-1-i]), 0.0)*xjump;// -G(tau)^2*G(-tau)=G(tau)^2*G(Beta-tau)
     // take the real part of G0t to enforce the symmetry of fouriere can increase the accuracy
     // this can make the calculation close to critical point as well as close to Uc2!
     //St[i]=G0t[i]*G0t[i]*G0t[i]*xjump;
     //cout<< G0t[i] <<"\t" << G0t[_Nw-1-i] << endl;
  }
  Fourier(St,Sw);
  //cout << "St to Sw done" << endl;
  // calculate G(w) with selfenergy
  for (int i=0; i<_Nw; i++) {
     Sw[i]=_U*_U*Sw[i]/xjump;
     Sw[i]=complex<double>(0.0, imag(Sw[i]));
     // Here we also enforce the symmetry of the half filling (at halfilling the real part is zero)
     // This can also makes the low temperature and Uc2 calculation more stable.
     //Gw[i] = complex<double>(0.0,imag(1.0/( (1/G0w[i]) - Sw[i]) ) );
     Gw[i] = 1.0/( (1/G0w[i]) - Sw[i]) ;
  }
}


double mpt::get_diff()
{
  double diff=0.0;//diff[_Nw];
  for (int i=0; i<_Nw; i++)
  {
    diff += abs(Sw_old[i] -Sw[i] );
  }
  
  return diff;//TrapezIntegralMP( _Nw , diff , omega ); 
}

void mpt::update_G(double mix)
{
  for (int i=0; i<_Nw; i++)
  {
    G0w[i] = (1-mix)*G0w[i]+mix*Gw[i] ;
  }

}

void mpt::printG0w(const char* filename)
{
  FILE *f;
  f = fopen(filename, "w");  
  for(int i=0; i<_Nw; i++) 
     fprintf(f,"%.15le %.15le %.15le\n", imag(omega[i]), real(G0w[i]), imag(G0w[i]));
  fclose(f);

}

void mpt::printGw(const char* filename)
{
  FILE *f;
  f = fopen(filename, "w");
  for(int i=0; i<_Nw; i++)
     fprintf(f,"%.15le %.15le %.15le\n", imag(omega[i]), real(Gw[i]), imag(Gw[i]));
  fclose(f);

}

void mpt::printG0t(const char* filename)
{
  FILE *f;
  f = fopen(filename, "w");
  for(int i=0; i<_Nw; i++)
     fprintf(f,"%.15le %.15le %.15le\n", tau[i], real(G0t[i]), imag(G0t[i]));
  fclose(f);

}

void mpt::printSw(const char* filename)
{
  FILE *f;
  f = fopen(filename, "w");
  for(int i=0; i<_Nw; i++)
     fprintf(f,"%.15le %.15le %.15le\n", imag(omega[i]), real(Sw[i]), imag(Sw[i]));
  fclose(f);

}

void mpt::printSt(const char* filename)
{
  FILE *f;
  f = fopen(filename, "w");
  for(int i=0; i<_Nw; i++)
     fprintf(f,"%.15le %.15le %.15le\n", tau[i], real(St[i]), imag(St[i]));
  fclose(f);

}

// try fftw3 package
void mpt::InvFourier(complex<double> Gw[],complex<double> Gt[])
{
  fftw_complex in[_Nw], out[_Nw]; /* double [2] */
  fftw_plan p;

  for (int i=0; i<_Nw/2; i++)
  {
    in[i][0] = real(Gw[i]*exp1[i] );
    in[i][1] = imag(Gw[i]*exp1[i] );
  }
  for (int i=_Nw/2; i<_Nw; i++)
  {
    in[i][0] = -1*real(Gw[i]*exp1[i] );
    in[i][1] = -1*imag(Gw[i]*exp1[i] );
  }

  //cout<<"Inverse Fourier Transofrm"<<endl;
  p = fftw_plan_dft_1d(_Nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  for (int i=0; i<_Nw; i++) {
    Gt[i] = _T*(exp2[i] * complex<double>( out[i][0], out[i][1] )) +  complex<double>(corr_wtot[i], 0.0)  ;//.
    //cout<< real(G_t[i])<<endl;
  }

  fftw_destroy_plan(p);
  //for (int i=_Nw/2; i<_Nw; i++) {
  //  Gt[i] = _T*(exp2[i-_Nw] * complex<double>( data[2*i+1], data[2*i+2] )) ;//.
    //cout<< real(G_t[i])<<endl;
  //}

}

//--------------T to F: Prepare input data-----------

void mpt::Fourier(complex<double> Gt[],complex<double> Gw[])
{
  fftw_complex in[_Nw], out[_Nw]; /* double [2] */
  fftw_plan p;

  for (int i=0; i<_Nw; i++)
  {
    in[i][0] = real( Gt[i] *conj(exp1[i]) );
    in[i][1] = imag( Gt[i] *conj(exp1[i]) );
  }

  //cout<<"Fourier Transofrm"<<endl;
  p = fftw_plan_dft_1d(_Nw, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  for (int i=0; i<_Nw/2; i++) {
    Gw[i] = A_n[i]*_Beta/(_Nw) * ( (complex<double>( out[i][0], out[i][1] ))*conj(exp2[i]) );
//                                 - exp( complex<double>( 0.0, imag(omega[i])*tau[0] ) )* Gt[0]
//                                  - exp( complex<double>( 0.0, imag(omega[i])*tau[_Nw-1] )) * Gt[_Nw-1] ) ;
    //cout<< G_f[i]<<endl;
  }
  for (int i=_Nw/2; i<_Nw; i++) {
    Gw[i] = A_n[i]*(-1)*_Beta/(_Nw) * ( (complex<double>( out[i][0], out[i][1] ))*conj(exp2[i]) );
//                                  - exp( complex<double>( 0.0, imag(omega[i])*tau[0] ) )* Gt[0]
//                                  - exp( complex<double>( 0.0, imag(omega[i])*tau[_Nw-1] )) * Gt[_Nw-1] ) ;
    //cout<< G_f[i]<<endl;
  }

  fftw_destroy_plan(p);
}


//--------------F to T: Prepare input data-----------

/*void mpt::InvFourier(complex<double> Gw[],complex<double> Gt[])
{ 
  int isign=-1;
  double* data = new double [2*_Nw+1];
  
  for (int i=0; i<_Nw/2; i++)
  {
    data[2*i+1] = real(Gw[i]*exp1[i] );
    data[2*i+2] = imag(Gw[i]*exp1[i] );
  }
  for (int i=_Nw/2; i<_Nw; i++)
  {
    data[2*i+1] = -1*real(Gw[i]*exp1[i] );
    data[2*i+2] = -1*imag(Gw[i]*exp1[i] );
  } 

  //cout<<"Inverse Fourier Transofrm"<<endl;
  FFT(data, _Nw, isign);

  for (int i=0; i<_Nw; i++) {
    Gt[i] = _T*(exp2[i] * complex<double>( data[2*i+1], data[2*i+2] )) +  complex<double>(corr_wtot[i], 0.0)  ;//.
    //cout<< real(G_t[i])<<endl;
  }
  //for (int i=_Nw/2; i<_Nw; i++) {
  //  Gt[i] = _T*(exp2[i-_Nw] * complex<double>( data[2*i+1], data[2*i+2] )) ;//.
    //cout<< real(G_t[i])<<endl;
  //}
      
  delete [] data;
}

//--------------T to F: Prepare input data-----------

void mpt::Fourier(complex<double> Gt[],complex<double> Gw[])
{
  int isign=1;
  double* data = new double [2*_Nw+1];
  
  for (int i=0; i<_Nw; i++) 
  {
    data[2*i+1] = real( Gt[i] *conj(exp1[i]) );
    data[2*i+2] = imag( Gt[i] *conj(exp1[i]) );
  } 

  //cout<<"Fourier Transofrm"<<endl;
  FFT(data, _Nw, isign);

  for (int i=0; i<_Nw/2; i++) {
    Gw[i] = (A_n[i])*_Beta/(_Nw) * ((complex<double>( data[2*i+1], data[2*i+2] ))*conj(exp2[i])) ;
    //cout<< G_f[i]<<endl;
  }
  for (int i=_Nw/2; i<_Nw; i++) {
    Gw[i] = A_n[i]*(-1)*_Beta/(_Nw) * ((complex<double>( data[2*i+1], data[2*i+2] ))*conj(exp2[i])) ;
    //cout<< G_f[i]<<endl;
  }

  delete [] data;

}*/


//----------Fast Fourier-----------------------------
/*#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void mpt::FFT(double data[], unsigned long N, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n = N << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP*/
