#include <complex>
#include <vector>

using namespace std;

//======================== CONSTANTS ===============================//

const double pi = 3.14159265358979323846;
const double e = 2.71;
const complex<double> ii = complex<double>(0.0,1.0);

//======================= ROUTINES ==================================//

double sign(double x);
int int_sign(double x);

double sqr(double x);
int pow(int base, int exp);
complex<double> sqr(complex<double> x);
/*double abs(double x);*/
/*double abs(complex<double> x);*/

//--- integral ---//

double TrapezIntegral(int N, double Y[], double X[]);
complex<double> TrapezIntegral(int N, complex<double> Y[], double X[]);
double TrapezIntegralMP(int N, double Y[], double X[]);
complex<double> TrapezIntegralMP(int N, complex<double> Y[], double X[]);
//double TrapezIntegral(std::vector< double > Y, std::vector<double> X);
complex<double> TrapezIntegral(std::vector< complex<double> > Y, std::vector<double> X);
double EllipticIntegralFirstKind(double x);
double SI(double x);
complex<double> EllipticIntegralFirstKind(complex<double> x);
double interpl(int N, double* Y, double* X, double x);

complex<double> Hilbert(int N, complex<double> z, double w[], double d[]);
void KramarsKronig(int N, double w[], double imf[], double ref[]);

//======================== IO =======================================//

void PrintFunc(const char* FileName, int N, int M, double** Y, double* X);
void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X);
void PrintFunc(const char* FileName, int N, complex<double>* Y);
void PrintFunc(const char* FileName, int N, double* Y);
void PrintFunc(const char* FileName, int N, double* Y, double* X);
void PrintFunc(const char* FileName, std::vector< complex<double> > Y, std::vector<double> X);
void PrintFunc3D(const char* FileName, int N, complex<double>** Y, double* X);
void ReadFunc(const char* FileName, int &N, int &M, double** &X);

//===================vectors and matrices=============================//

void MultiplyByMatrix(int N, double* v, double** m);
void CreateRotationMatrix(double** m, int N, double angle, int* plane);
void RotateVector(int N, double* v, double angle, int* plane);
void RotateVector2D(double* v, double angle);
void InvertMatrix(int N, double** A, double** invA, double &det);

//==================== Init DOSes and Fermi ========================//

double FermiFunction(double x, double T);
double SemiCircleDOS(double x );
