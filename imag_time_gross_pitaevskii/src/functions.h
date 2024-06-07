
//this includes all the preconditoning plans for the funtions


#include <armadillo>
#include <complex>
#include <fftw3.h>

using namespace std;
using namespace arma;

void linearterm(cx_mat & ,cx_mat &);

cx_mat nonlinearterm(cx_mat , double, cx_mat ,cx_mat,cx_mat &,cx_mat &,fftw_plan,fftw_plan);


void read(cx_mat &,cx_mat &, double& , int&, fftw_plan);
void energy(double, int,cx_mat , cx_mat, double &);

void intensity(double, int, cx_mat );

void spectrum( cx_mat, rowvec &, rowvec & ,int, int );



void record_parameters();

void setup_ETDRK(cx_mat,cx_mat &,cx_mat &,cx_mat &,cx_mat &, cx_mat &, cx_mat &,cx_mat &, cx_mat &, cx_mat &, cx_mat &, double);