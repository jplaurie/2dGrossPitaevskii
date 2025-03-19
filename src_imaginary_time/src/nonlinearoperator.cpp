#include <armadillo>
#include <fftw3.h>
#include "const.h"
#include "functions.h"
#include <complex>
#include <cmath>

using namespace std;
using namespace arma;

static cx_mat Npsi_temp_M(Mx,My,fill::zeros), psi_hat_N(Nx,Ny,fill::zeros);

//function conditioners for embedding data in larger arrays anti-aliasing
void embed_N_M( cx_mat, cx_mat &);
void embed_M_N( cx_mat, cx_mat &);
void dealias(cx_mat &);

cx_mat nonlinearterm( cx_mat A , double da, cx_mat L,cx_mat G, cx_mat & psi_M,cx_mat & psi_hat_M, fftw_plan FFT, fftw_plan IFFT){
    
/*
	Subroutine to compute the nonlinear cubic term in the 2DGP equation.  Use a 3/2 dealias rule but applied twice to preserve momentum conservation

*/

psi_hat_N.zeros();
psi_hat_N = A;


if(FLAG_TIMESTEP_ETDRK == false){
	psi_hat_N %= exp(da*L);
}
    
Npsi_temp_M.zeros();
psi_hat_M.zeros();
psi_M.zeros();
   
embed_N_M(psi_hat_N,psi_hat_M);             //embeds data in larger array
  
fftw_execute(IFFT);                  //inverse fft size M

Npsi_temp_M = psi_M;	//temp stores psi_M

psi_M %= psi_M;
 	
fftw_execute(FFT);

psi_hat_M /= double(Mx*My);            //normalize

dealias(psi_hat_M);

fftw_execute(IFFT);

psi_M %= -0.5*g*conj(Npsi_temp_M);
 
fftw_execute(FFT);                  //forward fft size M

psi_hat_M /= double(Mx*My);            //normalize
   
embed_M_N(psi_hat_M,psi_hat_N);             //puts data back into smaller array size N

//psi_hat_N /= (complex<double>(0.0,1.0));

if(FLAG_TIMESTEP_ETDRK == false){
	psi_hat_N %= exp(-da*L);
}

return psi_hat_N;                   //returns nonlinear term
}




//================================================================================================================
void embed_N_M( cx_mat A, cx_mat & B){
    
    B.zeros();
    
    for(int j =0; j < Ny/2 ; j++){
        for(int i =0; i < Nx/2 ; i++){
       
            B(i,j) = A(i,j);
            B(Mx-i-1,j)= A(Nx-i-1,j);
            B(i,My-j-1) = A(i,Ny-j-1);
            B(Mx-i-1,My-j-1)= A(Nx-i-1,Ny-j-1);
        }
    }
     
    return;
}


void embed_M_N( cx_mat B, cx_mat & A){
    
A.zeros();
    
for(int j =0; j < Ny/2 ; j++){
    for(int i =0; i < Nx/2 ; i++){
        A(i,j) = B(i,j);
        A(Nx-i-1,j)= B(Mx-i-1,j);
        A(i,Ny-j-1) = B(i,My-j-1);
        A(Nx-i-1,Ny-j-1)= B(Mx-i-1,My-j-1);
    }
}
    
return;
}
    
 
void dealias(cx_mat & A){
     
for(int j = Ny/2; j < My - (Ny/2) ; j++){
    for(int i = 0; i < Mx  ; i++){
        A(i,j) = complex<double>(0.0,0.0);
    }
}

for(int j = 0; j < My ; j++){
    for(int i = Nx/2; i < Mx - (Nx/2) ; i++){
        A(i,j) = complex<double>(0.0,0.0);
    }
}

return;
}
       
    
    
