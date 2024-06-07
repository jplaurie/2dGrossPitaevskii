
//this creates the linear operator in fourier-space

//L = 0.5 * Delta PSI  + 0.5* mu * PSI  - i* U * Nabla PSI. 


#include <armadillo>
#include <fftw3.h>
#include "const.h"
#include "functions.h"
#include <cmath>
#include <complex>
using namespace std;
using namespace arma;


void linearterm(cx_mat & L,cx_mat & G){
    
    /*
     subroutine to define the linear operator of the laplacian, chemical potential and comoving velocity
    */

    double Ux = 0.0;
    double Uy = 0.0;
    if(FLAG_DIPOLE == 1){
        Ux = -2.0*pi/(2.0*pi*dipole_length);
    }

    double kx = 0.0;
    double ky = 0.0;

for(int j= 0; j < Ny/2; j++){
	for(int i = 0; i < Nx/2; i++){

        kx = 2.0 * pi * double(i) / Lx;
        ky = 2.0 * pi * double(j) / Ly;
		L(i,j) = complex<double>(-0.5* (pow(kx,2.0) +pow(ky,2.0)) ,0.0);
        L(i,j) += complex<double>( 0.5*chem_pot,0.0);
        L(i,j) += complex<double>(Ux*kx + Uy*ky ,0.0);

        kx = 2.0 * pi * double(-i-1) / Lx;
        ky = 2.0 * pi * double(j) / Ly;
		L(Nx-i-1,j) = complex<double>(-0.5* (pow(kx,2.0) +pow(ky,2.0)),0.0);
        L(Nx-i-1,j) += complex<double>( 0.5*chem_pot,0.0);
        L(Nx-i-1,j) += complex<double>(Ux*kx + Uy*ky ,0.0);

        kx = 2.0 * pi * double(i) / Lx;
        ky = 2.0 * pi * double(-j-1) / Ly;
		L(i,Ny-j-1) = complex<double>(-0.5* (pow(kx,2.0) +pow(ky,2.0)),0.0);
        L(i,Ny-j-1) += complex<double>(0.5*chem_pot,0.0);
        L(i,Ny-j-1) += complex<double>(Ux*kx + Uy*ky ,0.0);

        kx = 2.0 * pi * double(-i-1) / Lx;
        ky = 2.0 * pi * double(-j-1) / Ly;
		L(Nx-i-1,Ny-j-1) = complex<double>(-0.5* (pow(kx,2.0) +pow(ky,2.0)),0.0);
        L(Nx-i-1,Ny-j-1) += complex<double>( 0.5*chem_pot,0.0);
        L(Nx-i-1,Ny-j-1) += complex<double>(Ux*kx + Uy*ky ,0.0);
	}
}




return;
}






void setup_ETDRK(cx_mat L_E,cx_mat & Q1,cx_mat & Q2,cx_mat & Q3,cx_mat & Q4,cx_mat & Q5, cx_mat & F1, cx_mat & F2, cx_mat & F3, cx_mat & E1, cx_mat & E2, double h){

/*
Routine to set up commonly used operators for ETDRK4 timestepping routine
*/

cx_mat LR(Nx,Ny);
int M = 32;									//number of complex contour points


/*
This is the code for the original ETDRK4 from Cox and Matthews and also ETDRK-B (same)
*/

E1 = exp(h*L_E);
E2 = exp(0.5*h*L_E);


if(FLAG_ETDRK_ORDER == 2){
    
    for(int i  = 0; i < M; i++){

        LR = h*L_E + exp(complex<double>(0.0,1.0)*2.0*pi*(double(i)-0.5)/ double(M));
        Q1 += (exp(LR) - 1.0) / LR;
        F1 += (exp(LR) - 1.0 - LR) / pow(LR,2.0);

    }

    Q1 *= h / double(M);
    F1 *= h / double(M);

}
else if(FLAG_ETDRK_ORDER == 3){
    
    for(int i  = 0; i < M; i++){

        LR = h*L_E + exp(complex<double>(0.0,1.0)*2.0*pi*(double(i)-0.5)/ double(M));
    
        Q1 += (exp(0.5*LR) - 1.0) / LR;
        Q2 += (exp(LR) - 1.0) / LR; 
        F1 += (-4.0 - LR + exp(LR)%(4.0-3.0*LR+pow(LR,2.0))) / pow(LR,3.0);
        F2 += (2.0 + LR + exp(LR)%(-2.0 + LR)) / pow(LR,3.0);
        F3 += (-4.0 - 3.0*LR - pow(LR,2.0) + exp(LR)%(4.0-LR)) / pow(LR,3.0);
    }

    Q1 *= h / double(M);
    Q2 *= h / double(M);
    F1 *= h / double(M);
    F2 *= h / double(M);
    F3 *= h / double(M);
}
else if(FLAG_ETDRK_ORDER == 4){
for(int i  = 0; i < M; i++){

    LR = h*L_E + exp(complex<double>(0.0,1.0)*2.0*pi*(double(i)-0.5)/ double(M));
    
    Q1 += (exp(0.5*LR) - 1.0) / LR;
    Q2 += (exp(0.5*LR)%(LR-4.0) + LR + 4.0) / pow(LR,2.0); 
    Q3 += 2.0*(2.0*exp(0.5*LR) - LR - 2.0) / pow(LR,2.0); 
    Q4 += (exp(LR)%(LR-2.0) + LR + 2.0) / pow(LR,2.0); 
    Q5 += 2.0*(exp(LR) - LR - 1.0) / pow(LR,2.0); 
    F1 += (-4.0 - LR + exp(LR)%(4.0-3.*LR+pow(LR,2.0))) / pow(LR,3.0);
    F2 += (2.0 + LR + exp(LR)%(-2.0 + LR)) / pow(LR,3.0);
    F3 += (-4.0 - 3.0*LR - pow(LR,2.0) + exp(LR)%(4.0-LR)) / pow(LR,3.0);
}

Q1 *= h / double(M);
Q2 *= h / double(M);
Q3 *= h / double(M);
Q4 *= h / double(M);
Q5 *= h / double(M);
F1 *= h / double(M);
F2 *= h / double(M);
F3 *= h / double(M);
}
return;
}