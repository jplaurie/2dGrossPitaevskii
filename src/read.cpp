
//this code reads in the initial data from file ./Initial.dat

#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip> 
#include "const.h"

using namespace std;
using namespace arma;

static double ignore_in, real_part, imag_part;

void readData(cx_mat & psi,cx_mat & psi_hat, double & run_time_start, int & file_number_start, fftw_plan FFTN){
    
/*
    reads in initial data located in file with number given by curframe.dat
*/

ignore_in =0.0;
real_part=0.0;
imag_part=0.0;
    
fstream filein("./data/curframe.dat");
filein.precision(12);
filein >> run_time_start >> file_number_start;
    
if(file_number_start < 0){

	file_number_start = 0;

	ofstream fout_read("./data/psi.000000");
	fout_read.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
			fout_read << i*dx << " " << j*dy << " " << real(psi(i,j)) << " " << imag(psi(i,j)) << endl;
		}
		fout_read << endl;
	}
	fout_read.close();

	cout << "file_number < 0" << endl; 
	cout << "Simulation starting...from zero state" << endl;
	cout << "time = " << run_time_start << endl;
	cout << "file_number = 0" << endl; 

}
else{

	ostringstream in_data;
	in_data << "./data/psi." << setw(6) << setfill('0') << file_number_start << ends;
	string filename = in_data.str();
	ifstream filein2(filename.c_str());
	filein2.precision(12);
    
	for(int i = 0; i < Nx; i++){
		for(int j = 0; j < Ny; j++){
			filein2 >> ignore_in >> ignore_in >> real_part >> imag_part;
			psi(i,j)= complex<double>(real_part,imag_part);
		}    
	}
	cout << "Loading data from file psi." << setw(6) << setfill('0') << file_number_start << endl;
	cout << "runtime of simulation = " << run_time_start << endl;
	cout << "Loading successful" << endl;
}
    
fftw_execute(FFTN);
psi_hat /= double(Nx*Ny);

return;
}


void recordParameters(){


ofstream fout_data("./parameters.txt");

fout_data << "Nx = " << Nx << " Ny = " << Ny << endl;
fout_data << "Lx = " << Lx << " Ly = " << Ly << endl; 
fout_data << "c = " << c << " g = " << g << " mu = " << mu << endl;
fout_data <<"================ Dissipation =====================" << endl;
fout_data << "nupower = " << nupower << " nu = " << nu << endl;
fout_data << "alphapower = " << alphapower << " alpha = " << alpha << endl;
fout_data << "FLAG_GINZBURG_LANDAU_DISS = " << FLAG_GINZBURG_LANDAU_DISS << endl;
fout_data << "Ginzburg dissipation gamma = " << gamma_0 << endl;
fout_data <<"================ Timestepping =========================" << endl;
fout_data << "FLAG_TIMESTEP_METHOD = " << FLAG_TIMESTEP_METHOD << endl;
fout_data << "dt = " << dt << endl;
fout_data << "total_steps = " << total_steps << endl;
fout_data << "output_time = " << output_time << endl;
fout_data <<"================ Forcing =========================" << endl;
fout_data << "FLAG_FORCING_TYPE = " << FLAG_FORCING_TYPE << endl;
if(FLAG_FORCING_TYPE == "annulus"){
	fout_data << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	fout_data << "force_amplitude = " << force_amplitude << " kf = " << kf << " dk = " << dk << endl;
}
else if(FLAG_FORCING_TYPE == "gaussian"){
	fout_data << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	fout_data << "force_amplitude = " << force_amplitude << " kf = " << kf << " sigmaf = " << sigmaf << endl;
}
if(FLAG_FORCING_TYPE == "exponential"){
	fout_data << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	fout_data << "force_amplitude = " << force_amplitude << " kf = " << kf << " force_power = " << force_power << endl;
}
if(FLAG_SOUND_FILTER == true){
	fout_data <<"================ Sound Filter ====================" << endl;
	fout_data << "FLAG_SOUND_FILTER= " << FLAG_SOUND_FILTER << endl;
	fout_data << "sigma = " << sigma << endl;
}
fout_data.close();

cout << "Nx = " << Nx << " Ny = " << Ny << endl;
cout << "Lx = " << Lx << " Ly = " << Ly << endl; 
cout << "c = " << c << " g = " << g << " mu = " << mu << endl;
cout <<"================ Dissipation =====================" << endl;
cout << "nupower = " << nupower << " nu = " << nu << endl;
cout << "alphapower = " << alphapower << " alpha = " << alpha << endl;
cout << "FLAG_GINZBURG_LANDAU_DISS = " << FLAG_GINZBURG_LANDAU_DISS << endl;
cout << "Ginzburg dissipation gamma = " << gamma_0 << endl;
cout <<"================ Timestepping =========================" << endl;
cout << "FLAG_TIMESTEP_METHOD = " << FLAG_TIMESTEP_METHOD << endl;
cout << "dt = " << dt << endl;
cout << "total_steps = " << total_steps << endl;
cout << "output_time = " << output_time << endl;
cout <<"================ Forcing =========================" << endl;
cout << "FLAG_FORCING_TYPE = " << FLAG_FORCING_TYPE << endl;
if(FLAG_FORCING_TYPE == "annulus"){
	cout << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	cout << "force_amplitude = " << force_amplitude << " kf = " << kf << " dk = " << dk << endl;
}
else if(FLAG_FORCING_TYPE == "gaussian"){
	cout << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	cout << "force_amplitude = " << force_amplitude << " kf = " << kf << " sigmaf = " << sigmaf << endl;
}
if(FLAG_FORCING_TYPE == "exponential"){
	cout << "FLAG_FORCING_RESCALE = " << FLAG_FORCING_RESCALE << endl;
	cout << "force_amplitude = " << force_amplitude << " kf = " << kf << " force_power = " << force_power << endl;
}
if(FLAG_SOUND_FILTER == true){
	cout <<"================ Sound Filter ====================" << endl;
	cout << "FLAG_SOUND_FILTER= " << FLAG_SOUND_FILTER << endl;
	cout << "sigma = " << sigma << endl;
}
return;
}
