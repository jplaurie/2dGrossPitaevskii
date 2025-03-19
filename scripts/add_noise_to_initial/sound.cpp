
//standalone file that creates an initial condition


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fstream>
#include <fftw3.h>

using namespace std;
using namespace arma;

const int Nx = 512;
const int Ny = 512;
const double pi = 3.14159265358979323846;
const double Lx = 2.0*pi;
const double dx= Lx/ (double) Nx;
const double Ly = 2.0*pi;
const double dy= Ly / (double) Ny;
const double mu = 400.0;


const double T= 2.e-4;
const double kf = 32.0; // center of gaussian noise
const double sigma = 8.0;

const int FLAG_SOUND_TYPE = 1;//  0==no sound, 1==RJ, 2==expoential 

void define_k2(mat & );


int main(){
	
    cx_mat psi(Nx,Ny,fill::zeros), sound_hat(Nx,Ny,fill::zeros), sound(Nx,Ny,fill::zeros);
    mat k2(Nx,Ny,fill::zeros);
   

    fftw_plan IFFT;


     IFFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) sound_hat.memptr(), (fftw_complex*) sound.memptr(), FFTW_BACKWARD, FFTW_PATIENT);


     define_k2(k2);


 fstream filein("./psi.initial");
    filein.precision(12);
    

    double ignore = 0.0;
    double theta = 0.0;
    double real_part =0.0;
    double imag_part = 0.0;


 for(int i=0; i< Nx; i++){
        for(int j =0; j < Ny; j++){
   filein >> ignore >> ignore >> real_part >> imag_part;
   psi(i,j) = complex<double>(real_part, imag_part);
    }
}



unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
mt19937_64 generator(seed1);
uniform_real_distribution<double> distribution(0.0,2.0*pi) ;

if(FLAG_SOUND_TYPE==1){
     for(int i=0; i< Nx; i++){
     	for(int j =0; j < Ny; j++){

     		theta = distribution(generator);
     		sound_hat(i,j) =  sqrt(T / (k2(i,j) + mu))*exp(complex<double>(0.0,theta));
     		sound_hat(0,0) = complex<double>(0.0,0.0);
     	}
     }
}

if(FLAG_SOUND_TYPE==2){
	for(int i=0; i< Nx; i++){
		for(int j =0; j < Ny; j++){
			theta = distribution(generator);
			sound_hat(i,j) =  (T/sqrt(2.0*pi))*exp( -0.5*pow( (sqrt(k2(i,j)) - kf)/sigma,2.0) )*exp(complex<double>(0.0,theta));
			sound_hat(0,0) = complex<double>(0.0,0.0);
		}
         }
}

    
fftw_execute(IFFT);

  for(int i=0; i< Nx; i++){
        for(int j =0; j < Ny; j++){
            psi(i,j) *= (1.0 + sound(i,j));

}
}


  
    ofstream fout("./psi_w_sound.initial");				//open up file
	fout.precision(12);
    for(int i=0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            fout << i*dx << " " << j*dy << " " << std::real(psi(i,j)) << " " << std::imag(psi(i,j)) << endl;
       
        }
        fout << endl;
    }

    fout.close();

    fftw_destroy_plan(IFFT);

    return 0;

}



void define_k2(mat & k2){



  for(int j =0; j < Ny/2 ; j++){
        for(int i =0; i < Nx/2 ; i++){

        	k2(i,j) =   pow( 2.0*pi*double(i) / Lx, 2.0) + pow( 2.0*pi*double(j) / Ly, 2.0)    ; 
        	k2(i,Ny-j-1) =  pow( 2.0*pi*double(i) / Lx, 2.0) + pow( 2.0*pi*double(-j-1) / Ly, 2.0)    ; 
        	k2(Nx-i-1,j) =  pow(2.0*pi*double(-i-1) / Lx , 2.0) + pow( 2.0*pi*double(j) / Ly, 2.0)     ; 
        	k2(Nx-i-1,Ny-j-1) =  pow( 2.0*pi*double(-i-1) / Lx, 2.0) + pow(2.0*pi*double(-j-1) / Ly , 2.0)    ; 

 }
}
return;
}
