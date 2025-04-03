
// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include "const.h"
#include "nonlinear.h"

using namespace std;
using namespace arma;

static double k_bin = min(2.0*pi / Lx, 2.0*pi/ Ly);
static int ib;
static double k, k2, flux, eflux, psi2, psi_hat2;
static cx_mat nonlinear_hat(Nx,Ny,fill::zeros), temp_hat(Nx,Ny,fill::zeros), psi_dot_hat(Nx,Ny,fill::zeros),rho_dot_hat(Nx,Ny,fill::zeros),rho_hat(Nx,Ny,fill::zeros),rho_dot_M(Nx,Ny,fill::zeros),temp_M(Mx,My,fill::zeros);
static cx_mat lap(Nx,Ny, fill::zeros), gradx(Nx,Ny, fill::zeros),grady(Nx,Ny, fill::zeros), gradx_psi_hat(Nx,Ny,fill::zeros), grady_psi_hat(Nx,Ny,fill::zeros), gradx_psi_dot_hat(Nx,Ny,fill::zeros), grady_psi_dot_hat(Nx,Ny,fill::zeros);
//static complex<double> phi_hat;

void generateLaplacian(cx_mat & );
void generateGradx(cx_mat & );
void generateGrady(cx_mat & );


void computeEnergy(cx_mat psi_hat, cx_mat psi, double & linear_energy_out, double & potential_energy_out, double & nonlinear_energy_out){
    
   /*
    Subroutine that computes the components of the total energy and also energy and waveaction dissipation rate
   */

    k=0.0;
    
    linear_energy_out = 0.0;
    potential_energy_out = 0.0;
    nonlinear_energy_out = 0.0;
    
   
    for(int j = 0; j < Ny/2; j++){
        for(int i = 0; i < Nx/2; i++){

            k = pow(  pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)  ,0.5);
        
            psi_hat2 = pow( abs(psi_hat(i,j)),2.0);
            psi2 = pow(abs(psi(i,j)),2.0);
            
            linear_energy_out -= Lx * Ly * c * k*k  * psi_hat2;
            potential_energy_out += dx * dy * mu * psi2;
            nonlinear_energy_out += dx * dy * 0.5 * g * pow(psi2,2.0);
        


            k = pow(  pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)    ,0.5);
        
            psi_hat2 = pow( abs(psi_hat(Nx-i-1,j)),2.0);
            psi2 = pow(abs(psi(Nx-i-1,j)),2.0);
        
            linear_energy_out -= Lx * Ly * c * k*k * psi_hat2;
            potential_energy_out += dx * dy * mu * psi2;
            nonlinear_energy_out += dx * dy * 0.5 * g * pow(psi2,2.0);


            k = pow(  pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)    ,0.5);
            psi_hat2 = pow( abs(psi_hat(i,Ny-j-1)),2.0);
            psi2 = pow(abs(psi(i,Ny-j-1)),2.0);
       
            linear_energy_out -= Lx * Ly * c * k*k * psi_hat2;
            potential_energy_out += dx * dy * mu * psi2;
            nonlinear_energy_out += dx * dy * 0.5 * g * pow(psi2,2.0)  ;
       


        k = pow(  pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)    ,0.5);
        psi_hat2 = pow( abs(psi_hat(Nx-i-1,Ny-j-1)),2.0);
        psi2 = pow(abs(psi(Nx-i-1,Ny-j-1)),2.0);

        linear_energy_out -= Lx * Ly * c * k*k  * psi_hat2;
        potential_energy_out += dx * dy * mu * psi2;
        nonlinear_energy_out += dx * dy * 0.5 * g * pow( psi2,2.0)   ;
    
     
        }
    }
    return;
}

void computeWaveaction(cx_mat psi_hat, double & waveaction_out){
   /*
    Subroutine that computes the total waveaction
   */
        waveaction_out = 0.0;

        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++){
                waveaction_out +=  Lx * Ly * pow(abs(psi_hat(i,j)),2.0);
            }
        }
return;
}

void computeDissipationRate( cx_mat psi_hat, double & waveaction_dissipation_rate_alpha, double & waveaction_dissipation_rate_nu, double & energy_dissipation_rate_alpha, double & energy_dissipation_rate_nu ){


    k2=0.0;
    waveaction_dissipation_rate_alpha = 0.0;
    waveaction_dissipation_rate_nu = 0.0;
    energy_dissipation_rate_alpha = 0.0;
    energy_dissipation_rate_nu = 0.0;


    for(int j = 0; j < Ny/2; j++){
        for(int i = 0; i < Nx/2; i++){

        

            k2 =   pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0) ;
            psi_hat2 = pow( abs(psi_hat(i,j)),2.0);
        
            if(FLAG_HYPER_DISSIPATION == true){
                energy_dissipation_rate_nu -= 2.0* c *nu*pow(k2,nupower+1.0) * psi_hat2;
                waveaction_dissipation_rate_nu += 2.0*nu*pow(k2,nupower) * psi_hat2 ;
            }
            if(FLAG_HYPO_DISSIPATION == true && (i != 0 || j != 0)){
                energy_dissipation_rate_alpha -= 2.0* c *alpha*pow(k2,alphapower+1.0) * psi_hat2;          
                waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2 ;
            }

            k2 =   pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)  ;
            psi_hat2 = pow( abs(psi_hat(Nx-i-1,j)),2.0);
          
            if(FLAG_HYPER_DISSIPATION == true){
                energy_dissipation_rate_alpha -= 2.0* c *alpha*pow(k2,alphapower+1.0) * psi_hat2;          
                waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2 ;
            }    
            if(FLAG_HYPO_DISSIPATION == true){
                energy_dissipation_rate_alpha -= 2.0* c *alpha*pow(k2,alphapower+1.0) * psi_hat2;
                waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2 ;
            }

            k2 = pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)  ;
            psi_hat2 = pow( abs(psi_hat(i,Ny-j-1)),2.0);
         
            if(FLAG_HYPER_DISSIPATION == true){
                energy_dissipation_rate_nu -= 2.0* c *nu*pow(k2,nupower+1.0) * psi_hat2;
                waveaction_dissipation_rate_nu += 2.0*nu*pow(k2,nupower) * psi_hat2 ;
            }
            if(FLAG_HYPO_DISSIPATION == true){
                energy_dissipation_rate_alpha -= 2.0* c *alpha*pow(k2,alphapower+1.0) * psi_hat2;
                waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2 ;
            }

            k2 =  pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)  ;
            psi_hat2 = pow( abs(psi_hat(Nx-i-1,Ny-j-1)),2.0);
 
            if(FLAG_HYPER_DISSIPATION == true){
                energy_dissipation_rate_nu -= 2.0* c *nu*pow(k2,nupower+1.0) * psi_hat2;
                waveaction_dissipation_rate_nu += 2.0*nu*pow(k2,nupower) * psi_hat2 ;
            }
            if(FLAG_HYPO_DISSIPATION == true){
                energy_dissipation_rate_alpha -= 2.0* c *alpha*pow(k2,alphapower+1.0) * psi_hat2;
                waveaction_dissipation_rate_alpha += 2.0*alpha*pow(k2,alphapower) * psi_hat2 ;
            }
        }   
    }
    return;
}

void computeSpectrum( cx_mat psi_hat, rowvec & wave_spectrum){

  /*
    Subroutine that computes wave action spectrum
   */
  
    k = 0.0;
    
    wave_spectrum.zeros();

    
    
    for(int j =0; j< Ny/2; j++){
	    for(int i =0; i< Nx/2; i++){
        
		    k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
		    ib = int( k / k_bin);
		    wave_spectrum(ib) += pow( abs(psi_hat(i,j)) , 2.0);
		
		    k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
		    ib = int(k/k_bin);
			wave_spectrum(ib) += pow( abs(psi_hat(Nx-i-1,j)) , 2.0);
		   
		    k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
		    ib = int(k/k_bin);
		    wave_spectrum(ib) += pow( abs(psi_hat(i,Ny-j-1)) , 2.0);
	   
		    k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
		    ib = int(k/k_bin);
			wave_spectrum(ib) += pow( abs(psi_hat(Nx-i-1,Ny-j-1)) , 2.0);   
	    }
    }
return;

}

void computeFlux(cx_mat psi_hat, cx_mat L,cx_mat G,cx_mat & psi_M ,cx_mat & psi_hat_M, rowvec & wave_flux ,rowvec & energy_flux,fftw_plan FFT, fftw_plan IFFT ){


    nonlinear_hat = nonlinear(psi_hat,0.0, L ,G, psi_M, psi_hat_M, FFT, IFFT); 
    generateLaplacian(lap);
    generateGradx(gradx);
    generateGrady(grady);

    gradx_psi_hat = gradx % psi_hat;
    grady_psi_hat = grady % psi_hat;
    
    psi_dot_hat = ( -complex<double>(0.0,1.0) * c * lap % psi_hat - complex<double>(0.0,1.0) * mu * psi_hat + nonlinear_hat);

    gradx_psi_dot_hat = gradx % psi_dot_hat;
    grady_psi_dot_hat = grady % psi_dot_hat;
    
  

    // compute rho_hat
    embed_N_M(psi_hat,psi_hat_M);
    fftw_execute(IFFT); 
    temp_M = psi_M; // this is for rho_dot_hat
    psi_M %= conj(psi_M);
    fftw_execute(FFT); 
    psi_hat_M /= double(Mx*My);
    embed_M_N(psi_hat_M,rho_hat);
    


    // compute rho_dot_hat
    embed_N_M(psi_dot_hat , psi_hat_M);
    fftw_execute(IFFT); 
    psi_M %= 2.0*conj(temp_M);
    fftw_execute(FFT); 
    psi_hat_M /= double(Mx*My);
    embed_M_N(psi_hat_M,rho_dot_hat);


    k = 0.0;
    wave_flux.zeros();
    energy_flux.zeros();
    
    flux = 0.0;
    eflux = 0.0;

    for(int j =0; j< Ny/2; j++){
        for(int i =0; i< Nx/2; i++){
        
        
            k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            

            //flux = -2.0 * ( real(psi_hat(i,j))*real(nonlinear_hat(i,j)) + imag(psi_hat(i,j))*imag(nonlinear_hat(i,j)) );
            flux = -2.0 * real( psi_hat(i,j) * conj(psi_dot_hat(i,j))); //- this works to as above - this is form the definition of zhu et al
            eflux = 2.0 * c * real( gradx_psi_hat(i,j) * conj(gradx_psi_dot_hat(i,j)) +  grady_psi_hat(i,j) * conj(grady_psi_dot_hat(i,j))) + (mu * flux) - g * real( rho_hat(i,j) * conj(rho_dot_hat(i,j)) ); 

            for(int l = ib; l < NBIN; l++){
                wave_flux(l) += flux;
                energy_flux(l) += eflux;
            }

            k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            //flux = -2.0 * ( real(psi_hat(i,Ny-j-1))*real(nonlinear_hat(i,Ny-j-1)) + imag(psi_hat(i,Ny-j-1))*imag(nonlinear_hat(i,Ny-j-1)));
            flux = -2.0 * real( psi_hat(i,Ny-j-1) * conj(psi_dot_hat(i,Ny-j-1)));
            eflux = 2.0 * c * real( gradx_psi_hat(i,Ny-j-1) * conj(gradx_psi_dot_hat(i,Ny-j-1)) +  grady_psi_hat(i,Ny-j-1) * conj(grady_psi_dot_hat(i,Ny-j-1))) + (mu * flux) - g * real(rho_hat(i,Ny-j-1) * conj(rho_dot_hat(i,Ny-j-1)) ); 
            


            for(int l = ib; l < NBIN; l++){
                wave_flux(l) += flux;
                energy_flux(l) += eflux;
            }

            k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            //flux = -2.0 * ( real(psi_hat(Nx-i-1,j))*real(nonlinear_hat(Nx-i-1,j)) + imag(psi_hat(Nx-i-1,j))*imag(nonlinear_hat(Nx-i-1,j)));
            flux = -2.0 * real( psi_hat(Nx-i-1,j) * conj(psi_dot_hat(Nx-i-1,j)));
            eflux = 2.0 * c * real( gradx_psi_hat(Nx-i-1,j) * conj(gradx_psi_dot_hat(Nx-i-1,j)) +  grady_psi_hat(Nx-i-1,j) * conj(grady_psi_dot_hat(Nx-i-1,j))) + (mu * flux)  - g  * real(rho_hat(Nx-i-1,j) * conj(rho_dot_hat(Nx-i-1,j)) ); 
            

            for(int l = ib; l < NBIN; l++){
                wave_flux(l) += flux;
                energy_flux(l) += eflux;
            }


            k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            ib = int(k/k_bin);
            
            //flux = -2.0 * ( real(psi_hat(Nx-i-1,Ny-j-1))*real(nonlinear_hat(Nx-i-1,Ny-j-1)) + imag(psi_hat(Nx-i-1,Ny-j-1))*imag(nonlinear_hat(Nx-i-1,Ny-j-1)));
            flux = -2.0 * real( psi_hat(Nx-i-1,Ny-j-1) * conj(psi_dot_hat(Nx-i-1,Ny-j-1)));
            eflux = 2.0 * c * real( gradx_psi_hat(Nx-i-1,Ny-j-1) * conj(gradx_psi_dot_hat(Nx-i-1,Ny-j-1)) +  grady_psi_hat(Nx-i-1,Ny-j-1) * conj(grady_psi_dot_hat(Nx-i-1,Ny-j-1))) + (mu * flux) - g * real(rho_hat(Nx-i-1,Ny-j-1) * conj(rho_dot_hat(Nx-i-1,Ny-j-1)) ); 
            

            for(int l = ib; l < NBIN; l++){
                wave_flux(l) += flux;
                energy_flux(l) += eflux;
            }

        }
    }
    return;
}


void generateLaplacian(cx_mat & A){

    for(int j= 0; j < Ny/2; j++){
	    for(int i = 0; i < Nx/2; i++){

	        A(i,j) = complex<double>(- (pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)),0.0);
		    A(Nx-i-1,j) = complex<double>(-  (pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(j) / Ly),2.0)) ,0.0);
		    A(i,Ny-j-1) = complex<double>(-  (pow((2.0 * pi * double(i) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)) ,0.0);
		    A(Nx-i-1,Ny-j-1) = complex<double>(- (pow((2.0 * pi * double(-i-1) / Lx),2.0) +pow((2.0 * pi * double(-j-1) / Ly),2.0)),0.0);
	    }
    }

    return;
}

void generateGradx(cx_mat & A){

    for(int j= 0; j < Ny; j++){
	    for(int i = 0; i < Nx/2; i++){

	        A(i,j) = complex<double>(0.0, 2.0 * pi * double(i) / Lx);
		    A(Nx-i-1,j) = complex<double>( 0.0,  2.0 * pi * double(-i-1) / Lx);
		   	   
        }
    }

    return;
}

void generateGrady(cx_mat & A){

    for(int j= 0; j < Ny/2; j++){
	    for(int i = 0; i < Nx; i++){

	        A(i,j) = complex<double>(0.0, 2.0 * pi * double(j) / Ly);
            A(i,Ny-j-1) = complex<double>(0.0, 2.0 * pi * double(-j-1) / Ly);
		   	   
        }
    }

    return;
}






