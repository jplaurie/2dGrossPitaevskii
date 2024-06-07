
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
void flux( cx_mat,cx_mat, rowvec &, rowvec & ,int, int);
//void forcing(cx_mat,mat);

void forcing(mat & , double, double);
void forcing2(mat & , double, double);


//cx_mat Filter(cx_mat,cx_mat, fftw_plan, fftw_plan);

//cx_mat Filter_Physical(cx_mat,cx_mat, fftw_plan, fftw_plan);
//void VortexNum(double, cx_mat);
void record_parameters();

void setup_ETDRK(cx_mat,cx_mat &,cx_mat &,cx_mat &,cx_mat &, cx_mat &, cx_mat &,cx_mat &, cx_mat &, cx_mat &, cx_mat &, double);