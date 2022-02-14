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

#define AMIN 1e-8
#define AMAX 1
#define LNASTEP 1840

#define CALMMIN 1e-3
#define CALMMAX 10
//#define LNCALMSTEP 9 //920

#define CHIMIN 1e-7
#define CHIMAX 1e-2
#define LNCHISTEP 115 //1150

#define QMIN 0.01
#define QMAX 1
#define QSTEP 99


double PcalMqchi(double calM, double q, double chi);
double PcalMqchiint(double a1, double a2, double M1, double M2, double nu1, double nu2, double q, double chi);
double Ph(double h);
double Cnu(double M, double nu);
double Na(double M, double nu);
double Panu(double a, double M, double nu);
double Pnu(double nu);
double LL(double a1, double a2, double chi, double q);
  
double Minnu(double nu);
double nuinM(double M);
double M1incalMq(double calM, double q);
double M2incalMq(double calM, double q);


int main(int argc, char *argv[])
{
  if (argc != 3) {
    cout << "iM, lncalMstep を正しく指定してください" << endl;
    return 1;
  }

  int iM = atoi(argv[1]);
  int lncalMstep = atoi(argv[2]);

  if (iM < 0 || iM > lncalMstep) {
    cout << "iM, lncalMstep を正しく指定してください" << endl;
    return 1;
  }
  
  // ---------------- start stop watch -----------------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // ---------------------------------------------------

  double dlncalM = (log(CALMMAX)-log(CALMMIN))/lncalMstep;
  double dlnchi = (log(CHIMAX)-log(CHIMIN))/LNCHISTEP;
  double dq = (QMAX-QMIN)/QSTEP;

  double calM = CALMMIN*exp(iM*dlncalM);

  cout << "calM = " << calM << endl;

#ifdef _OPENMP
  cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << endl;
#endif

  //double data[LNCALMSTEP+1][QSTEP+1][LNCHISTEP+1];
  double data[QSTEP+1][LNCHISTEP+1];
  
  int done = 0;
  //cout << "\r" << setw(3) << 100*done/(QSTEP+1) << "%" << flush;
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  //for (int iq=0; iq<=QSTEP; iq++) {
  int iq = 50;
  double q = QMIN+iq*dq;
    
  //for (int ic=0; ic<=LNCHISTEP; ic++) {
  int ic = 50;
  double chi = CHIMIN*exp(ic*dlnchi);
      
  data[iq][ic] = PcalMqchi(calM,q,chi);
  //}
    
#ifdef _OPENMP
#pragma omp critical
#endif
  //{
  //done++;
  //cout << "\r" << setw(3) << 100*done/(QSTEP+1) << "%" << flush;
  //}
  //}
  
  //cout << endl;

  cout << calM << ' ' << q << ' ' << chi << ' ' << data[iq][ic] << endl;
  
  /*
  string str = "PcalMqchi.dat";
  ofstream ofs(str,std::ios::app);

  for (int iq=0; iq<=QSTEP; iq++) {
    double q = QMIN+iq*dq;
    
    for (int ic=0; ic<=LNCHISTEP; ic++) {
      double chi = CHIMIN*exp(ic*dlnchi);
      
      ofs << calM << ' ' << q << ' ' << chi << ' ' << data[iq][ic] << endl;
    }
  }
  */

    

  // ---------------- return elapsed time --------------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // ---------------------------------------------------
}


double PcalMqchi(double calM, double q, double chi) {
  double PcalMqchi = 0;
  double dlna = (log(AMAX)-log(AMIN))/LNASTEP;

  double M1 = M1incalMq(calM,q);
  double M2 = M2incalMq(calM,q);
  double nu1 = nuinM(M1);
  double nu2 = nuinM(M2);
    
  double a1, a2;

  for (int i1 = 0; i1 <= LNASTEP; i1++) {
    a1 = AMIN*exp(dlna*i1);
      
    for (int i2 = 0; i2 <= LNASTEP; i2++) {
      a2 = AMIN*exp(dlna*i2);

      if (LL(a1,a2,chi,q) > 0) {
	if ((i1 == 0 || i1 == LNASTEP) && (i2 == 0 || i2 == LNASTEP)) {
	  PcalMqchi += PcalMqchiint(a1,a2,M1,M2,nu1,nu2,q,chi)/4.*dlna*dlna;
	} else if (i1 == 0 || i1 == LNASTEP || i2 == 0 || i2 == LNASTEP) {
	  PcalMqchi += PcalMqchiint(a1,a2,M1,M2,nu1,nu2,q,chi)/2*dlna*dlna;
	} else {
	  PcalMqchi += PcalMqchiint(a1,a2,M1,M2,nu1,nu2,q,chi)*dlna*dlna;
	}
      }
    }
  }

  double coeff = (1+q)/4/q/q/BETA/BETA/SIGMA0/SIGMA0*pow(pow(1+q,2./5)*calM*calM/pow(q,1./5)/KK/KK,1./BETA);
  return chi*coeff*PcalMqchi;
}

double PcalMqchiint(double a1, double a2, double M1, double M2, double nu1, double nu2, double q, double chi) {
  return LL(a1,a2,chi,q)*Panu(a1,M1,nu1)*Pnu(nu1)*Panu(a2,M2,nu2)*Pnu(nu2);
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
  return pow(M/KK,1./BETA)/SIGMA0 + NUTH;
}

double M1incalMq(double calM, double q) {
  return pow(q,-3./5)*pow(1+q,1./5)*calM;
}

double M2incalMq(double calM, double q) {
  return pow(q,2./5)*pow(1+q,1./5)*calM;
}
