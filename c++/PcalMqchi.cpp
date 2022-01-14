#define _USE_MATH_DEFINES

#include <sys/time.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#define BETA 0.36
#define KK 3.3
#define SIGMA0 0.08
#define NUTH 24
#define GAMMA 0.85
#define NA 0.9999088

double Ph(double h);
double Cnu(double M, double nu);
double Na(double M, double nu);
double Panu(double a, double M, double nu);
double Pnu(double nu);
double LL(double a1, double a2, double chi, double q);
  
double Minnu(double nu);
double nuinM(double M);


int main()
{
  // ---------------- start stop watch -----------------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // ---------------------------------------------------

  double nu = NUTH + 1e-2;

  cout << Panu(1e-10,Minnu(nu),nu) << ' '
       << Panu(1e-8,Minnu(nu),nu) << ' '
       << Panu(1e-6,Minnu(nu),nu) << ' '
       << Panu(1e-4,Minnu(nu),nu) << ' '
       << Panu(1e-2,Minnu(nu),nu) << endl;

  // ---------------- return elapsed time --------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // ---------------------------------------------------
}

double PlncalMqchi(double calM, double q, double chi) {
  
}


double Ph(double h) {
  return 563*h*h*exp(-12*h+2.5*pow(h,1.5)+8-3.2*pow(1500+pow(h,16),1./8));
}

double Cnu(double M, double nu) {
  return (1.56e-3)*sqrt(1-GAMMA*GAMMA)*pow(M,-1./3)*pow(2./5*nu/8,-2);
}

double Na(double M, double nu) {
  return NA;
}

double Panu(double a, double M, double nu) {
  return Ph(a/Cnu(M,nu))/Na(M,nu)/Cnu(M,nu);
}

double Pnu(double nu) {
  return sqrt(2./M_PI)*exp(-nu*nu/2)/erfc(NUTH/sqrt(2));
}

double LL(double a1, double a2, double chi, double q) {
  return min(a1,q*a2+(1+q)+chi) + min(a1,q*a2-(1+q)*chi);
}


double Minnu(double nu) {
  return KK*pow(nu*SIGMA0-NUTH*SIGMA0,BETA);
}

double nuinM(double M) {
  return pow(M,1./BETA)/SIGMA0 + NUTH;
}
