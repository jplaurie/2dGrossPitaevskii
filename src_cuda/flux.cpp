
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include "functions.h"
#include <cufftw.h>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;

static int ib = 0;
static rowvec wave_flux(NBIN),energy_flux(NBIN);
static double k_bin,k,fwave;


void flux(cx_mat A, cx_mat B, rowvec & wave_flux_avg, rowvec & energy_flux_avg, int out_count, int avg_count){

    
    k=0.0;
    wave_flux.zeros();
    energy_flux.zeros();
    k_bin = min(2.0*pi / Lx, 2.0*pi/ Ly);
    fwave=0.0;
   

    for(int j =0; j< Ny/2; j++){
        for(int i =0; i< Nx/2; i++){
        
        
        
            k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            fwave = -2.0 * ( real(A(i,j))*real(B(i,j)) + imag(A(i,j))*imag(B(i,j)) );
  
            wave_flux(ib) += fwave;
            energy_flux(ib) += k*k*fwave;
      

            k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            fwave = -2.0 * ( real(A(i,Ny-j-1))*real(B(i,Ny-j-1)) + imag(A(i,Ny-j-1))*imag(B(i,Ny-j-1)));
           
            wave_flux(ib) += fwave;
            energy_flux(ib) += k*k*fwave;

            k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            fwave = -2.0 * ( real(A(Nx-i-1,j))*real(B(Nx-i-1,j)) + imag(A(Nx-i-1,j))*imag(B(Nx-i-1,j)));

            wave_flux(ib) += fwave;
            energy_flux(ib) += k*k*fwave;


            k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            fwave = -2.0 * ( real(A(Nx-i-1,Ny-j-1))*real(B(Nx-i-1,Ny-j-1)) + imag(A(Nx-i-1,Ny-j-1))*imag(B(Nx-i-1,Ny-j-1)));

            wave_flux(ib) += fwave;
            energy_flux(ib) += k*k*fwave;

        }
    }

    wave_flux_avg += wave_flux;
    energy_flux_avg += energy_flux;

    
    ostringstream out_flux;
    out_flux << "./output/flux." << setw(5) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
    string filename = out_flux.str();
    ofstream fout_flux(filename.c_str());
    fout_flux << scientific;
    fout_flux.precision(12);
    
    for(int i=0; i < NBIN; i++){
            
            fout_flux << (i*k_bin) << " " << wave_flux(i) << " " << energy_flux(i) << " " << wave_flux_avg(i) / double(avg_count) << " " << energy_flux_avg(i) / double(avg_count) <<  endl;
    }
    
    fout_flux.close();
   
    
    return;
}
