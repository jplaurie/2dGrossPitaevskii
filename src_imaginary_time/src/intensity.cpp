
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include "functions.h"
#include <fftw3.h>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;

void intensity(double runtime, int out_count, cx_mat psi){
    
ostringstream out_I;
out_I << "./data/psi." << setw(5) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
string filename = out_I.str();
ofstream fout_I(filename.c_str());
fout_I << scientific;
fout_I.precision(12);
    
for(int i=0; i < Nx; i++){
	for(int j=0; j < Ny; j++){
		fout_I << i*dx << " " << j*dy << " " << real(psi(i,j)) << " " << imag(psi(i,j)) << endl;
	}
	fout_I << endl;
}

fout_I.close();

return;
    
}