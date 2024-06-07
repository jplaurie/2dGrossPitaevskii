
//this code reads in the initial data from file ./Initial.dat

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

static double ignore_in, real_part, imag_part;

void read(cx_mat & psi,cx_mat & psi_hat, double& runtime, int& out_count, fftw_plan FFTN){
    
/*
    reads in initial data located in file with number given by curframe.dat
*/

ignore_in =0.0;
real_part=0.0;
imag_part=0.0;
    
fstream filein("./data/curframe.dat");
filein.precision(12);
filein >> runtime >> out_count;
    
if(out_count < 0){

	out_count=0;

	ofstream fout_read("./data/psi.00000");
	fout_read.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
			fout_read << i*dx << " " << j*dy << " " << real(psi(i,j)) << " " << imag(psi(i,j)) << endl;
		}
		fout_read << endl;
	}
	fout_read.close();

	cout << "out_count < 0" << endl; 
	cout << "Simulation starting...from zero state" << endl;
	cout << "time = " << runtime << endl;
	cout << "out_count = 0" << endl; 

}
else{

	ostringstream in_data;
	in_data << "./data/psi." << setw(5) << setfill('0') << out_count << ends;
	string filename = in_data.str();
	ifstream filein2(filename.c_str());
	filein2.precision(12);
    
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			filein2 >> ignore_in >> ignore_in >> real_part >> imag_part;
			psi(i,j)= complex<double>(real_part,imag_part);
		}    
	}
	cout << "Loading data from file psi." << setw(5) << setfill('0') << out_count << endl;
	cout << "runtime of simulation = " << runtime << endl;
	cout << "Loading successful" << endl;
}
    
fftw_execute(FFTN);
psi_hat /= double(Nx*Ny);

return;
}


void record_parameters(){


ofstream fout_data("./parameters.txt");

fout_data << "Nx = " << Nx << " Ny = " << Ny << endl;
fout_data << "Lx = " << Lx << " Ly = " << Ly << endl; 
fout_data << "g = " << g << endl;
fout_data << "chemical potential = " << chem_pot << endl;
fout_data <<"================ Dissipation =====================" << endl;
fout_data << "FLAG_GINZBURG_LANDAU_DISS = " << FLAG_GINZBURG_LANDAU_DISS << endl;
fout_data << "Ginzburg dissipation gamma = " << gamma_0 << endl;
fout_data << "nupower = " << nupower << " nu = " << nu << endl;
fout_data << "alphapower = " << alphapower << " alpha = " << alpha << endl;
fout_data <<"================ Timestepping =========================" << endl;
fout_data << "FLAG_TIMESTEP_ETDRK = " << FLAG_TIMESTEP_ETDRK << endl;
fout_data << "FLAG_ETDRK_ORDER = " << FLAG_ETDRK_ORDER << endl;
fout_data << "dt = " << dt << endl;
fout_data << "nsteps = " << nsteps << endl;
fout_data << "outstep = " << outstep << endl;
fout_data <<"================ Forcing =========================" << endl;
fout_data << "FLAG_FORCE_TYPE = " << FLAG_FORCE_TYPE << endl;
fout_data << "FLAG_FORCE_EXP = " << FLAG_FORCE_EXP << endl;
fout_data << "FLAG_FORCE_AMP_RESCALE = " << FLAG_FORCE_AMP_RESCALE << endl;
fout_data << "kf = " << kf << " dk = " << dk << " Amp = " << Amp << endl;
fout_data << "force_power = " << force_power << endl;
fout_data <<"================ Sound Filter ====================" << endl;
fout_data << "FLAG_SOUND_FILTER= " << FLAG_SOUND_FILTER << endl;
fout_data << "sigma = " << sigma << endl;

fout_data.close();

cout << "FLAG_TIMESTEP_ETDRK = " << FLAG_TIMESTEP_ETDRK << endl;
cout << "FLAG_ETDRK_ORDER = " << FLAG_ETDRK_ORDER << endl;
cout << "FLAG_FORCE_TYPE = " << FLAG_FORCE_TYPE << endl;
cout << "FLAG_FORCE_EXP = " << FLAG_FORCE_EXP << endl;
cout << "FLAG_FORCE_AMP_RESCALE = " << FLAG_FORCE_AMP_RESCALE << endl;
cout << "FLAG_GINZBURG_LANDAU_DISS = " << FLAG_GINZBURG_LANDAU_DISS << endl;  
cout << "FLAG_SOUND_FILTER = " << FLAG_SOUND_FILTER << endl;

cout << "================ Parameters =====================" << endl;
cout << "Nx = " << Nx << " Ny = " << Ny << endl;
cout << "Lx = " << Lx << " Ly = " << Ly << endl; 
cout << "g = " << g << endl;
cout << "chemical potential = " << chem_pot << endl;
cout <<"================ Dissipation =====================" << endl;
  
cout << "nupower = " << nupower << " nu = " << nu << endl;
cout << "alphapower = " << alphapower << " alpha = " << alpha << endl;
if(FLAG_GINZBURG_LANDAU_DISS == true){
cout << "Ginzburg dissipation gamma = " << gamma_0 << endl; 
} 
cout <<"================ Timestepping =========================" << endl;
cout << "dt = " << dt << endl;
cout << "nsteps = " << nsteps << endl;
cout << "outstep = " << outstep << endl;

if(FLAG_FORCE_TYPE != 0){
cout <<"================ Forcing =========================" << endl;
cout << "kf = " << kf << " dk = " << dk << " Amp = " << Amp << endl;
cout << "force_power = " << force_power << endl;
}
if(FLAG_SOUND_FILTER == true){
cout <<"================ Sound Filter ====================" << endl;
cout << "sigma = " << sigma << endl;
}
return;
}
