
// This compute the linear and nonlinear energies, and total wave action


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

static double Lin_E, Pot_E,Non_E, Int,k,wave_diss_alpha,wave_diss_nu,energy_diss_alpha,energy_diss_nu, psi_hat2, psi2;

void energy(double runtime, int out_count, cx_mat psi_hat, cx_mat psi, double & E_out){
    
   /*
    Subroutine that computes the components of the total energy and also energy and waveaction dissipation rates


   */

    k=0.0;
    wave_diss_alpha = 0.0;
    wave_diss_nu = 0.0;
    energy_diss_alpha = 0.0;
    energy_diss_nu = 0.0;
    
    Lin_E = 0.0;
    Pot_E = 0.0;
    Non_E = 0.0;
    Int = 0.0;
   
    for(int j = 0; j < Ny/2; j++){
        for(int i = 0; i < Nx/2; i++){

        k = pow(  pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)  ,0.5);
        psi_hat2 = pow( abs(psi_hat(i,j)),2.0);
        psi2 = pow(abs(psi(i,j)),2.0);
        Lin_E += Lx * Ly * 0.5 * k*k  * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow(psi2,2.0);
        Int += dx * dy * psi2;
        
        wave_diss_alpha += 2.0*alpha*pow(k*k,alphapower) * psi_hat2 ;
        wave_diss_nu += 2.0*nu*pow(k*k,nupower) * psi_hat2 ;
        energy_diss_alpha += 2.0*alpha*pow(k*k,alphapower+1.0) * psi_hat2;
        energy_diss_nu += 2.0*nu*pow(k*k,nupower+1.0) * psi_hat2;

        k = pow(  pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)    ,0.5);
        
        psi_hat2 = pow( abs(psi_hat(Nx-i-1,j)),2.0);
        psi2 = pow(abs(psi(Nx-i-1,j)),2.0);
        
        Lin_E += Lx * Ly * 0.5 * k*k * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow(psi2,2.0);
        Int += dx * dy * psi2;

        wave_diss_alpha += 2.0*alpha*pow(k*k,alphapower) * psi_hat2 ;
        wave_diss_nu += 2.0*nu*pow(k*k,nupower) * psi_hat2 ;
        energy_diss_alpha += 2.0*alpha*pow(k*k,alphapower+1.0) * psi_hat2;
        energy_diss_nu += 2.0*nu*pow(k*k,nupower+1.0) * psi_hat2;

        k = pow(  pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)    ,0.5);
        psi_hat2 = pow( abs(psi_hat(i,Ny-j-1)),2.0);
        psi2 = pow(abs(psi(i,Ny-j-1)),2.0);
       
        Lin_E += Lx * Ly * 0.5 * k*k * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow(psi2,2.0)  ;
        Int += dx * dy * psi2;
        
        wave_diss_alpha += 2.0*alpha*pow(k*k,alphapower) * psi_hat2 ;
        wave_diss_nu += 2.0*nu*pow(k*k,nupower) * psi_hat2 ;
        energy_diss_alpha += 2.0*alpha*pow(k*k,alphapower+1.0) * psi_hat2;
        energy_diss_nu += 2.0*nu*pow(k*k,nupower+1.0) * psi_hat2;

        k = pow(  pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)    ,0.5);
        psi_hat2 = pow( abs(psi_hat(Nx-i-1,Ny-j-1)),2.0);
        psi2 = pow(abs(psi(Nx-i-1,Ny-j-1)),2.0);

        Lin_E += Lx * Ly * 0.5 * k*k  * psi_hat2;
        Pot_E -= dx * dy * 0.5 * chem_pot * psi2;
        Non_E += dx * dy * 0.25 * g * pow( psi2,2.0)   ;
        Int += dx * dy * psi2;

        wave_diss_alpha += 2.0*alpha*pow(k*k,alphapower) * psi_hat2 ;
        wave_diss_nu += 2.0*nu*pow(k*k,nupower) * psi_hat2 ;
        energy_diss_alpha += 2.0*alpha*pow(k*k,alphapower+1.0) * psi_hat2;
        energy_diss_nu += 2.0*nu*pow(k*k,nupower+1.0) * psi_hat2;


        E_out = Lin_E + Pot_E + Non_E;


        }
    }
    


    ostringstream out_D;
    out_D << "./output/dissipation." << setw(5) << setfill('0') << out_count << ends;         //creates file name for outputting data at time slice
    string filenameD = out_D.str();
    ofstream fout_D(filenameD.c_str());
    fout_D << scientific;
    fout_D.precision(12);

    fout_D << runtime << " " << wave_diss_alpha << " " << wave_diss_nu << " " << energy_diss_alpha << " " << energy_diss_nu << endl;
    
    fout_D.close();

    ostringstream out_E;
    out_E << "./output/energy." << setw(5) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
    string filenameE = out_E.str();
    ofstream fout_E(filenameE.c_str());
    fout_E << scientific;
    fout_E.precision(12);
    
        fout_E << runtime << " " << Lin_E << " " << Pot_E << " " << Non_E << " " << Lin_E + Pot_E + Non_E << endl;
    
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
