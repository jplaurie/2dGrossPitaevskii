#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>

#include <fftw3.h>
#include <fstream>
#include "const.h"
using namespace std;
using namespace arma;

static double k,kx,ky,waveaction_inj, energy_inj;
static int forcing_modes;


void initialiseForcing(mat & famp){
    
if(FLAG_FORCING_TYPE == "annulus"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;

	for(int j = 0; j < Ny/2; j++){
		for(int i = 0; i < Nx/2; i++){
        
		if((i==0) && (j == 0)) continue;                

		//define k        
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            
		if( abs(k-kf) < dk){
			famp(i,j) = force_amplitude;
			waveaction_inj += pow(famp(i,j),2.0);
			energy_inj += pow(famp(i,j),2.0)*(-c * k*k);
			forcing_modes += 1;
		}

		//define k
		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
               	 
		if( abs(k-kf) < dk){
			famp(Nx-i-1,j) = force_amplitude;
			waveaction_inj += pow(famp(Nx-i-1,j),2.0);
			energy_inj += pow(famp(Nx-i-1,j),2.0)*(-c * k*k);
			forcing_modes += 1;
		}

		//define k
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
              
		if( abs(k-kf) < dk){
			famp(i,Ny-j-1) = force_amplitude;
			waveaction_inj += pow(famp(i,Ny-j-1),2.0);
			energy_inj += pow(famp(i,Ny-j-1),2.0)*(-c * k*k);
			forcing_modes += 1;
		}

		//define k
		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
              	
		if( abs(k-kf) < dk){
			famp(Nx-i-1,Ny-j-1) = force_amplitude;
			waveaction_inj += pow(famp(Nx-i-1,Ny-j-1),2.0);
			energy_inj += pow(famp(Nx-i-1,Ny-j-1),2.0)*( - c * k*k);
			forcing_modes += 1;
		}

	

		}
	}
}
else if(FLAG_FORCING_TYPE == "gaussian"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;

	for(int j = 0; j < Ny/2; j++){
		for(int i = 0; i < Nx/2; i++){
	
		//define k
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);

		//set forcing amplitude 
		famp(i,j) = force_amplitude * exp(-0.5 * pow( (k-kf)/sigmaf,2.0) );
		waveaction_inj += pow(famp(i,j),2.0);
		energy_inj += pow(famp(i,j),2.0)*(-c * k*k);
		forcing_modes += 1;
        
		//define k   
		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
           
		famp(Nx-i-1,j) = force_amplitude * exp(-0.5 * pow( (k-kf)/sigmaf,2.0) );
		waveaction_inj += pow(famp(Nx-i-1,j),2.0);
		energy_inj += pow(famp(Nx-i-1,j),2.0)*( - c * k*k);
		forcing_modes += 1;

		//define k    
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            
		famp(i,Ny-j-1) = force_amplitude * exp(-0.5 * pow( (k-kf)/sigmaf,2.0) );
		waveaction_inj += pow(famp(i,Ny-j-1),2.0);
		energy_inj += pow(famp(i,Ny-j-1),2.0)*(- c * k*k);
		forcing_modes += 1;

		//define k
		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
          
		famp(Nx-i-1,Ny-j-1) = force_amplitude * exp(-0.5 * pow( (k-kf)/sigmaf,2.0) );
		waveaction_inj += pow(famp(Nx-i-1,Ny-j-1),2.0);
		energy_inj += pow(famp(Nx-i-1,Ny-j-1),2.0)*(- c * k*k);
		forcing_modes += 1;
          		
		}
	}

}
else if(FLAG_FORCING_TYPE == "exponential"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;
	
    
	for(int j = 0; j < Ny/2; j++){
		for(int i = 0; i < Nx/2; i++){
        
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            
		famp(i,j) = pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		waveaction_inj += pow(famp(i,j),2.0);
		energy_inj += pow(famp(i,j),2.0)*( - c *k*k);
		forcing_modes += 1;
            
		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
           
		famp(Nx-i-1,j) = pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		waveaction_inj += pow(famp(Nx-i-1,j),2.0);
		energy_inj += pow(famp(Nx-i-1,j),2.0)*( - c * k*k);
		forcing_modes += 1;
            
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            
		famp(i,Ny-j-1) = pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		waveaction_inj += pow(famp(i,Ny-j-1),2.0);
		energy_inj += pow(famp(i,Ny-j-1),2.0)*(-c * k*k);
		forcing_modes += 1;

		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
          
		famp(Nx-i-1,Ny-j-1) = pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		waveaction_inj += pow(famp(Nx-i-1,Ny-j-1),2.0);
		energy_inj += pow(famp(Nx-i-1,Ny-j-1),2.0)*(-c * k*k);
		forcing_modes += 1;
          		
		}       		
	}
}
else if(FLAG_FORCING_TYPE == "exp-log"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;
	
    
	for(int j = 0; j < Ny/2; j++){
		for(int i = 0; i < Nx/2; i++){
        
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
            
		famp(i,j) = force_amplitude * exp(-0.5*pow(log(k/kf)/sigmaf,2.0));
		waveaction_inj += pow(famp(i,j),2.0);
		energy_inj += pow(famp(i,j),2.0)*( - c *k*k);
		forcing_modes += 1;
            
		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);
           
		famp(Nx-i-1,j) =  force_amplitude * exp(-0.5*pow(log(k/kf)/sigmaf,2.0));
		waveaction_inj += pow(famp(Nx-i-1,j),2.0);
		energy_inj += pow(famp(Nx-i-1,j),2.0)*( - c * k*k);
		forcing_modes += 1;
            
		k = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
            
		famp(i,Ny-j-1) = force_amplitude *  exp(-0.5*pow(log(k/kf)/sigmaf,2.0));
		waveaction_inj += pow(famp(i,Ny-j-1),2.0);
		energy_inj += pow(famp(i,Ny-j-1),2.0)*(-c * k*k);
		forcing_modes += 1;

		k = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);
          
		famp(Nx-i-1,Ny-j-1) =  force_amplitude * exp(-0.5*pow(log(k/kf)/sigmaf,2.0));
		waveaction_inj += pow(famp(Nx-i-1,Ny-j-1),2.0);
		energy_inj += pow(famp(Nx-i-1,Ny-j-1),2.0)*(-c * k*k);
		forcing_modes += 1;
          		
		}       		
	}
}
else{
    cout << "FLAG_FORCING_TYPE not defined" << endl;
    exit(1);
}

    
if(FLAG_FORCING_RESCALE == true){
	famp *= sqrt(2.0*alpha/waveaction_inj);
	energy_inj *= (2.0*alpha/waveaction_inj);
	waveaction_inj = 2.0*alpha;
}
    
cout << "number of forcing modes = " << forcing_modes << endl;
cout << "energy injection rate = " << energy_inj << endl;
cout << "waveaction injection rate = " << waveaction_inj << endl;
    	
ofstream fout_force_data("./output/force_values.txt");
fout_force_data.precision(12);

fout_force_data << "number of forcing modes = " << forcing_modes << endl;
fout_force_data << "energy injection rate = " << energy_inj << endl;
fout_force_data << "waveaction injection rate = " << waveaction_inj << endl;
fout_force_data << "forcing mode kf = " << kf << endl;
fout_force_data << "forcing width dk = " << dk << endl;

fout_force_data.close();
    
ofstream fout_Amp("./output/force_spectrum.dat");
fout_Amp << scientific;
fout_Amp.precision(12);
 
for(int i =Nx/2; i< Nx; i++){
       
	kx = 2.0*pi * double(i-Nx)/ Lx;
      
	for(int j =Ny/2; j< Ny; j++){
             
		ky = 2.0*pi * double(j-Ny)/ Ly;
		fout_Amp << kx << " " << ky << " " << famp(i,j) << endl;
	}
	for(int j =0; j< Ny/2; j++){
            
		ky = 2.0*pi * double(j)/ Ly;
		fout_Amp << kx << " " << ky << " " << famp(i,j) << endl;
	}
	fout_Amp << endl;
}
for(int i =0; i< Nx/2; i++){
    		
	kx = 2.0*pi * double(i)/ Lx;
        
	for(int j =Ny/2; j< Ny; j++){
                  
		ky = 2.0*pi * double(j-Ny)/ Ly;
		fout_Amp << kx << " " << ky << " " << famp(i,j) << endl;     
	}
	for(int j =0; j< Ny/2; j++){
             
		ky = 2.0*pi * double(j)/ Ly;
		fout_Amp << kx << " " << ky << " " << famp(i,j) << endl;          
	}
	fout_Amp << endl;
}

fout_Amp.close();


return;
}


