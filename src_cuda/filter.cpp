

// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "Const.h"
#include "Functions.h"
#include <cufftw.h>
#include <fstream>

using namespace std;
using namespace arma;



static cx_mat G(Nx,Ny);
static blitz::Array<int,2> find_vortex(Nx,Ny);
static blitz::Array<complex<double>,2> psi_temp(Nx,Ny);
static blitz::Array<complex<double>,2> psi_hat_temp(Nx,Ny);
static double waveaction,dist,ddx,ddy;
static int ii,jj,vortex_near;


blitz::Array<complex<double>,2> Filter(blitz::Array<complex<double>,2> psi_hat, blitz::Array<complex<double>,2> psi, fftw_plan IFFTN, fftw_plan FFTN){
    


    G=0.0;
    find_vortex = 0;
    psi_hat_temp = complex<double>(0.,0.);
    psi_temp = complex<double>(0.,0.);
    waveaction =0.0;
    dist= 0.0;
    vortex_near=0;
    psi_hat_temp = psi_hat;
    
    //===========defines gaussian kernel =======================================
    
    blitz::firstIndex I;
    blitz::secondIndex J;
	
	G(Range(0,Nx/2 -1), Range(0, Ny/2-1)) =  exp( - pow(    (pow((2.0 * pi * I / Lx),2.0) +pow((2.0 * pi * J / Ly),2.0)) /(sigma*sigma)  ,2.0));
	G(Range(Nx/2 ,Nx-1), Range(0, Ny/2-1)) = exp( -pow((pow((2.0 * pi * (I-(Nx/2) ) / Lx),2.0) +pow((2.0 * pi * J / Ly),2.0))/(sigma*sigma),2.0)) ;
    
    G(Range(0,Nx/2 -1), Range(Ny/2, Ny-1)) =  exp(-pow((pow((2.0 * pi * I / Lx),2.0) +pow((2.0 * pi * (J-(Ny/2)) / Ly),2.0))/(sigma*sigma),2.0)) ;
    G(Range(Nx/2 ,Nx-1), Range(Ny/2, Ny-1)) = exp( -pow((pow((2.0 * pi * (I-(Nx/2) ) / Lx),2.0) +pow((2.0 * pi * (J-(Ny/2)) / Ly),2.0))/(sigma*sigma),2.0));
    
    
    
    psi_hat *= G;
    
    
    
    //=========================================================================
    
    fftw_execute(IFFTN);
    
    psi_temp = psi;    // this is now the smoothed field
    
    psi_hat = psi_hat_temp;
    
    fftw_execute(IFFTN);
    
    waveaction=0.0;
    
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            
            waveaction += abs(psi(i,j));
        }
    }
	waveaction /= double(Nx*Ny);

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            
            if ( (abs(real(psi(i,j))) / waveaction) < eps  &&  (abs(imag(psi(i,j)))/ waveaction) < eps){

            	find_vortex(i,j) = 1;
            	
            	for(int ii = 0; ii < Nx; ii++){
        			for(int jj = 0; jj < Ny; jj++){

        				dist = 0.0;
        				if(i-ii < -Nx/2){
        					dist = pow(abs(i-ii+Nx)*dx,2.0);
        				}
        				else if(i-ii >= Nx/2){
        					dist = pow(abs(i-ii-Nx)*dx,2.0);
        				}

        				else{
        					dist = pow(abs(i-ii)*dx,2.0);
        				}

        				if(j-jj < -Ny/2){
        					dist += pow(abs(j-jj+Ny)*dy,2.0);
        				}
        				else if	(j-jj >= Ny/2){ 
        					dist += pow(abs(j-jj-Ny)*dy,2.0);
        				}
        				else{
        				 dist += pow(abs(j-jj)*dy,2.0);
        				}


            			if(dist < vortex_ring_dist){

            				find_vortex(ii,jj) = 1;
            			}
            		}

            	}
            }
        }
    }




    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){

        	if(find_vortex(i,j)==1){
				psi_temp(i,j) = psi(i,j);
        	}
                
            
        }
    }
    
    
    psi = psi_temp;
    
    fftw_execute(FFTN);
    
    psi_hat /= (double) Nx*Ny;
    

    
    return psi_hat;
    
}


blitz::Array<complex<double>,2> Filter_Physical(blitz::Array<complex<double>,2> psi_hat,blitz::Array<complex<double>,2> psi, fftw_plan IFFTN, fftw_plan FFTN){
    
    find_vortex = 0;
 
    psi_temp = complex<double>(0.,0.);
    waveaction =0.0;
    ii=0;jj=0;
    //=========================================================================
    
    fftw_execute(IFFTN);
    
   
    
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            
            waveaction += abs(psi(i,j));
        }
    }

	waveaction /= double(Nx*Ny);

    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            
            if ( (abs(real(psi(i,j))) / waveaction) < eps  &&  (abs(imag(psi(i,j)))/ waveaction) < eps){

            	find_vortex(i,j) = 1;
            	
            	}
        }
    }
         


    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){

       // 	if(find_vortex(i,j) == 1){
        //		psi_temp(i,j) = psi(i,j);
        //	}
        
            vortex_near = 0;

        			for(int i2 = -3*(int)sigma; i2 <= 3*(int)sigma; i2++){
        				for(int j2 = -3*(int)sigma; j2 <= 3*(int)sigma; j2++){

                            ii = i2 + i;
                            jj = j2 + j;
        					
        					if(ii < 0) ii += Nx;
        					else if(ii > Nx - 1) ii -= Nx;
        					if(jj < 0) jj += Ny;
        					else if(jj > Ny - 1) jj -= Ny;
        					
                            if(find_vortex(ii,jj) == 1){
                                vortex_near = 1;
                            }

        					ddx =  i2;
                            ddy = j2;

        					psi_temp(i,j) += (1./(2.*pi*sigma*sigma)) * exp(-(0.5/(sigma*sigma))*(pow(ddx,2.0) + pow(ddy,2.0)))*psi(ii,jj);
        				}
        			}
                    if(vortex_near == 1){
                        psi_temp(i,j) = psi(i,j);
                    }
       	

        			

    	}
    }
        
    
    psi = psi_temp;
    
    fftw_execute(FFTN);
    
    psi_hat /= (double) Nx*Ny;
    
    return psi_hat;
    
}


void VortexNum(double runtime, blitz::Array<complex<double>,2> psi){
    
    int countvor=0;


    waveaction=0.0;
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            
            waveaction += abs(psi(i,j));
        }
    }
	waveaction /= (double) Nx*Ny;
    
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            
            if ( (abs(real(psi(i,j)))/ waveaction) < eps   &&  (abs(imag(psi(i,j)))/ waveaction) < eps){
                
                countvor++;

            }
        }
    }
        
            
            ofstream fout_vor("./Output/Vortexnumber.dat", ios::out | ios::app);
            fout_vor.precision(12);
            
            fout_vor << runtime << " " << countvor << endl;
            
            fout_vor.close();
            
    
            
            return;
    
}

