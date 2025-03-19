
// This compute the linear and nonlinear energies, and total wave action


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

static double Lin_E, Pot_E,Non_E, Int,kx,ky,k, psi_hat2, psi2, Momx,Momy, Ux, Uy;

void energy(double runtime, int out_count, cx_mat psi_hat, cx_mat psi, double & E_out){
    
   /*
    Subroutine that computes the components of the total energy and also energy and waveaction dissipation rates


   */

    k=0.0;
    kx=0.0;
    ky=0.0;

    if(FLAG_DIPOLE == true){
        Ux = - 2.0*pi/(2.0*pi*dipole_length);
        Uy = 0.0;
    }
    
    Lin_E = 0.0;
    Pot_E = 0.0;
    Non_E = 0.0;

    Momx=0.0;
    Momy=0.0;

    Int = 0.0;
   
    for(int j = 0; j < Ny/2; j++){
        for(int i = 0; i < Nx/2; i++){


        kx = 2.0 * pi * double(i) / Lx;
        ky = 2.0 * pi * double(j) / Ly;
        k = pow( pow( kx ,2.0) +pow(ky,2.0)  ,0.5);

        psi_hat2 = pow( abs(psi_hat(i,j)),2.0);
        psi2 = pow(abs(psi(i,j)),2.0);
        Lin_E += Lx * Ly * 0.5 * k*k  * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow(psi2,2.0);

        Momx += Lx * Ly * kx * psi_hat2; 
        Momx += Lx * Ly * ky * psi_hat2; 

        Int += dx * dy * psi2;
        
   
        kx = 2.0 * pi * double(-i-1) / Lx;
        ky = 2.0 * pi * double(j) / Ly;
        k = pow( pow( kx ,2.0) +pow(ky,2.0)  ,0.5);
      
        
        psi_hat2 = pow( abs(psi_hat(Nx-i-1,j)),2.0);
        psi2 = pow(abs(psi(Nx-i-1,j)),2.0);
        
        Lin_E += Lx * Ly * 0.5 * k*k * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow(psi2,2.0);
        Momx += Lx * Ly * kx * psi_hat2; 
        Momx += Lx * Ly * ky * psi_hat2; 
        Int += dx * dy * psi2;


        kx = 2.0 * pi * double(i) / Lx;
        ky = 2.0 * pi * double(-j-1) / Ly;
        k = pow( pow( kx ,2.0) +pow(ky,2.0)  ,0.5);

        psi_hat2 = pow( abs(psi_hat(i,Ny-j-1)),2.0);
        psi2 = pow(abs(psi(i,Ny-j-1)),2.0);
       
        Lin_E += Lx * Ly * 0.5 * k*k * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow(psi2,2.0);
        Momx += Lx * Ly * kx * psi_hat2; 
        Momx += Lx * Ly * ky * psi_hat2; 
        Int += dx * dy * psi2;
        

        kx = 2.0 * pi * double(-i-1) / Lx;
        ky = 2.0 * pi * double(-j-1) / Ly;
        k = pow( pow( kx ,2.0) +pow(ky,2.0)  ,0.5);
        
        psi_hat2 = pow( abs(psi_hat(Nx-i-1,Ny-j-1)),2.0);
        psi2 = pow(abs(psi(Nx-i-1,Ny-j-1)),2.0);

        Lin_E += Lx * Ly * 0.5 * k*k  * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow( psi2,2.0);
        Momx += Lx * Ly * kx * psi_hat2; 
        Momx += Lx * Ly * ky * psi_hat2; 
        Int += dx * dy * psi2;




        E_out = Lin_E + Pot_E + Non_E- Ux*Momx - Uy*Momy;


        }
    }
    



    ostringstream out_E;
    out_E << "./output/energy." << setw(5) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
    string filenameE = out_E.str();
    ofstream fout_E(filenameE.c_str());
    fout_E << scientific;
    fout_E.precision(12);
    
        fout_E << runtime << " " << Lin_E << " " << Pot_E << " " << Non_E << " " <<  Momx << " " << Momy << " " <<  Lin_E + Pot_E + Non_E - Ux*Momx - Uy*Momy << endl;
    
    fout_E.close();
    
    ostringstream out_W;
    out_W << "./output/waveaction." << setw(5) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
    string filenameW = out_W.str();
    ofstream fout_W(filenameW.c_str());
    fout_W << scientific;
    fout_W.precision(12);

    fout_W << runtime << " " << Int << " " << pow( abs(psi_hat(0,0)),2.0) << endl;
    
    fout_W.close();
    
    
   ostringstream out_L;
    out_L << "./output/coherencelength." << setw(5) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
    string filenameL = out_L.str();
    ofstream fout_L(filenameL.c_str());
    fout_L << scientific;
    fout_L.precision(12);

    fout_L << runtime << " " << 1.0/sqrt(g*pow( abs(psi_hat(0,0)),2.0)) << endl;
    
    fout_L.close();

    
    
    return;
}