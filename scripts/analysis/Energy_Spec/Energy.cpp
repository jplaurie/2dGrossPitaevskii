
// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <blitz/array.h>
#include <complex>
#include <fftw3.h>
#include <fstream>

using namespace std;
using namespace blitz;

const int Nx=512;
const int Ny=512;
const int Mx = 2*Nx;                     
const int My = 2*Ny;                      
const double pi = 3.14159265358979323846;       //pi
const double Lx = 32.0*pi;                   //length of the box
const double dx = Lx/ (double) Nx;                      //grid size
const double Ly = 32.0*pi;                   //length of the box
const double dy = Ly/ (double) Ny;  
const double g=1.0;
const double chem_pot=0.0;


const int startfile = 700;
const int endfile = 900;
const int NBIN = 2.0*max(Nx,Ny);
const double k_bin = min(2.0*pi / Lx, 2.0*pi/ Ly);
const double dxm = Lx/ (double) Mx;
const double dym = Ly/ (double) My;
static blitz::Array<int,1> specnum(NBIN);


void embed_N_M(blitz::Array<complex<double>,2> , blitz::Array<complex<double>,2>);
void embed_M_N(blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>);
void embed_N_M_real(blitz::Array<complex<double>,2> , blitz::Array<complex<double>,2>);
void embed_M_N_real(blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>);
void define_k(blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>);
void compute_spectrum(blitz::Array<complex<double>,2>, blitz::Array<complex<double>,2> , blitz::Array<double,1>);


int main(){
    

blitz::Array<double,2> V(Nx,Ny),V_M(Mx,My),Vx(Nx,Ny), Vy(Nx,Ny), Vx_inc(Nx,Ny), Vy_inc(Nx,Ny),Vx_com(Nx,Ny), Vy_com(Nx,Ny), rhox(Nx,Ny), rhoy(Nx,Ny);
blitz::Array<complex<double>,2> psi(Nx,Ny),psi_hat(Nx,Ny),temp(Nx,Ny),temp_hat(Nx,Ny),kx(Mx,My),ky(Mx,My),k2(Mx,My), temp_M(Mx,My), temp_hat_M(Mx,My), V_hat_M(Mx,My/2 +1), psi_x(Mx,My), psi_y(Mx,My),Div_hat_M(Mx,My/2+1), specx(Nx,Ny/2+1),specy(Nx,Ny/2+1),kxr(Mx,My/2+1),kyr(Mx,My/2+1),k2r(Mx,My/2+1),V_hat(Nx,Ny/2+1),real_hat_M(Mx,My/2+1);

blitz::Array<double,1> kin_spec(NBIN), inc_spec(NBIN), com_spec(NBIN),qua_spec(NBIN),spectrum(NBIN);



    double E_kin_tot, E_inc_tot,E_com_tot, E_qua_tot,E_int_tot, E_kin_spec_tot,E_inc_spec_tot,E_com_spec_tot,E_qua_spec_tot,E_pot_tot;

    E_com_tot=0.0;E_inc_tot=0.0;E_kin_tot=0.0;E_qua_tot=0.0;E_int_tot=0.0,E_pot_tot=0.0;;
    E_kin_spec_tot=0.0;E_inc_spec_tot=0.0;E_com_spec_tot=0.0;E_qua_spec_tot=0.0;
    double real_in,imag_in,ignore;
   


    fftw_plan FFT;
    fftw_plan IFFT;
    fftw_plan IFFTM;
    fftw_plan FFTVM;
    fftw_plan FFTV;
    fftw_plan IFFTV;

    FFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) temp.data(), (fftw_complex*) temp_hat.data(), FFTW_FORWARD, FFTW_PATIENT);
    IFFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) temp_hat.data(), (fftw_complex*) temp.data(), FFTW_BACKWARD, FFTW_PATIENT);
    IFFTM = fftw_plan_dft_2d(Mx,My, (fftw_complex*) temp_hat_M.data(), (fftw_complex*) temp_M.data(), FFTW_BACKWARD, FFTW_PATIENT);
    FFTVM = fftw_plan_dft_r2c_2d(Mx,My, (double*) V_M.data(), (fftw_complex*) V_hat_M.data(), FFTW_PATIENT);
    FFTV = fftw_plan_dft_r2c_2d(Nx,Ny, (double*) V.data(), (fftw_complex*) V_hat.data(), FFTW_PATIENT);
    IFFTV = fftw_plan_dft_c2r_2d(Nx,Ny, (fftw_complex*) V_hat.data(), (double*) V.data(), FFTW_PATIENT);

    kin_spec = 0.0;
    inc_spec =0.0;
    com_spec= 0.0;
    qua_spec = 0.0;
    Vx=0.0;
    Vy=0.0;
    rhox=0.0;
    rhoy=0.0;
    Vx_inc=0.0;
    Vy_inc=0.0;
    Vx_com=0.0;
    Vy_com=0.0;
    kx=complex<double>(0.0,0.0);
    ky=complex<double>(0.0,0.0);
    k2=complex<double>(0.0,0.0);
    psi=complex<double>(0.0,0.0);
    psi_hat=complex<double>(0.0,0.0);
    temp=complex<double>(0.0,0.0);
    temp_M=complex<double>(0.0,0.0);
    temp_hat=complex<double>(0.0,0.0);
    temp_hat_M=complex<double>(0.0,0.0);
    psi_x=complex<double>(0.0,0.0);
    psi_y=complex<double>(0.0,0.0);
    real_hat_M = complex<double>(0.0,0.0);
    V_M=0.0;
    V_hat_M=complex<double>(0.0,0.0);
    Div_hat_M=complex<double>(0.0,0.0);
    temp = complex<double>(0.0,0.0);

define_k(kx, ky,k2,kxr,kyr,k2r);        //defines the arrays to compute x and y derivative and laplacian


for(int k=startfile; k <= endfile; k++){        //increment over file number

E_com_tot=0.0;E_inc_tot=0.0;E_kin_tot=0.0;E_qua_tot=0.0;E_int_tot=0.0,E_pot_tot=0.0;;
    E_kin_spec_tot=0.0;E_inc_spec_tot=0.0;E_com_spec_tot=0.0;E_qua_spec_tot=0.0;
//=====================================================================================================================
//================================== loads data from file - wave function psi ======================================
//=====================================================================================================================
    cout << "file = " << k << endl;

	ostringstream laurie;
	laurie << "../../gamma_0=1em2/data/psi."<<  setw(5) << setfill('0') << k << ".dat" << ends;
	
	string filename = laurie.str();
	ifstream data_file(filename.c_str());
	data_file.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
 
       data_file >> ignore >> ignore >> ignore >> real_in >> imag_in ;
            psi(i,j) = complex<double>(real_in,imag_in);
		//	psi(i,j) = complex<double>(cos(1.*i*dx+2.*j*dy),sin(3*i*dx+5*j*dy));

         //   cout << i << " " << j << " " << psi(i,j) << endl;
      }
    }


//=====================================================================================================================
//=====================================================================================================================
//============================computes the kinetic energy ========================

temp = psi;         //transfers to array temp.  Keeps copy of data in psi AT ALL TIMES
   	
fftw_execute(FFT);  //transforms original data to Fourier space

temp_hat /= (double) Nx * Ny; // normalizes

psi_hat = temp_hat;      //transfers FFT data to psi_hat.  Keeps copy of wavenumbers of original data at all times.

//======================computes x-derivative of psi ========================================
embed_N_M(temp_hat,temp_hat_M);         //embeds data to larger array

temp_hat_M *= kx;                // computes x derivative in fourier space

fftw_execute(IFFTM);            // inverse fourier transform

psi_x = temp_M;                     //this is now d_x\psi
//=============================================================================================

//======================computes y-derivative of psi ========================================

embed_N_M(psi_hat,temp_hat_M);      

temp_hat_M *= ky;            //computes y-derivative



fftw_execute(IFFTM);  //inverse fourier transform to physical space

psi_y = temp_M;          //this stores d_y \psi


//=============================================================================================

//==================  computes psi in larger mesh =================================================
embed_N_M(psi_hat,temp_hat_M);      //embeds data to larger array

fftw_execute(IFFTM);            // inverse fourier transform to larger array temp_M
//==============================================================================================

//=============================compute the x-velocity =======================

V_M= real(complex<double>(0.0,-0.5) * (conj(temp_M)*psi_x - temp_M*conj(psi_x)) / abs(temp_M));   //computes x-velocity in physical space

fftw_execute(FFTVM);            //fourier transform velocity

V_hat_M /= (double) Mx * My;     //normalize

embed_M_N_real(V_hat_M,V_hat);

specx = V_hat;
//================================================================
//=============record the velocity divergence =========================
Div_hat_M = V_hat_M * kxr;        //this is the first part of the velocity divergence: dx Vx

//===================record the x-velocity in smaller array ====================

//embed_N_M_real(V_hat_M,V_hat);


fftw_execute(IFFTV);     //inverse fourier transform 

Vx = V;     //record

embed_N_M(psi_hat,temp_hat_M);

fftw_execute(IFFTM);

V_M = real(complex<double>(0.0,-0.5) * (conj(temp_M)*psi_y - temp_M*conj(psi_y)) / abs(temp_M)); //zip(2.0*imag(conj(temp_M)*psi_y)/ abs(temp_M),0.0,complex<double>());  //computes y velocity

fftw_execute(FFTVM);

V_hat_M /= (double) Mx * My;     //normalize

embed_M_N_real(V_hat_M,V_hat);

specy = V_hat;

compute_spectrum(specx,specy,spectrum);  // kinetic energy

kin_spec = spectrum;


Div_hat_M += V_hat_M * kyr;        //divergence of the velocity nabla,V

//embed_M_N_real(V_hat_M,V_hat);     //embed x velocity to smaller array

fftw_execute(IFFTV);         //inverse fourier transform

Vy = V;         //save 

real_hat_M = kxr*Div_hat_M/k2r;

embed_M_N_real(real_hat_M,V_hat);

specx -= V_hat;

real_hat_M = kyr*Div_hat_M/k2r;

embed_M_N_real(real_hat_M,V_hat);

specy -= V_hat;

compute_spectrum(specx,specy,spectrum); // incompressible

inc_spec = spectrum;

real_hat_M = kxr*Div_hat_M/k2r;
embed_M_N_real(real_hat_M,specx);
real_hat_M = kyr*Div_hat_M/k2r;
embed_M_N_real(real_hat_M,specy);
//specx = kxr*Div_hat_M/k2r;
//specy = kyr*Div_hat_M/k2r; 

compute_spectrum(specx,specy,spectrum); // compressible

com_spec = spectrum;

V_hat_M = kyr*Div_hat_M/k2r;

embed_M_N_real(V_hat_M,V_hat);     //embed vorticity to smaller array

fftw_execute(IFFTV);     //inverse fourier transform

Vy_com = V;    //save

V_hat_M = kxr*Div_hat_M/k2r;
embed_M_N_real(V_hat_M,V_hat);     //embed vorticity to smaller array

fftw_execute(IFFTV);     //inverse fourier transform

Vx_com = V;    //save

Vx_inc = Vx - Vx_com;
Vy_inc = Vy - Vy_com;

//==========compute the quantum energy===========================================
//===============================================================================
V=0.0;
temp = psi;

fftw_execute(FFT);

temp_hat /= (double) Nx * Ny; // normalizes

embed_N_M(temp_hat,temp_hat_M);         //embeds data to larger array


fftw_execute(IFFTM);            // inverse fourier transform


V_M = abs(temp_M); // = rho

fftw_execute(FFTVM);

V_hat_M /= (double) Mx * My;

//embed_N_M_real(V_hat,V_hat_M);

real_hat_M = V_hat_M * kxr;
embed_M_N_real(real_hat_M,specx);
real_hat_M = V_hat_M * kyr;
embed_M_N_real(real_hat_M,specy);

compute_spectrum(specx,specy,spectrum);


qua_spec = spectrum;

V_hat=specy;

fftw_execute(IFFTV);
rhoy=V;


V_hat=specx;
fftw_execute(IFFTV);
rhox=V;


for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){

        E_kin_tot += 0.5*dx*dy*(pow(Vx(i,j),2.0)+pow(Vy(i,j),2.0));
        E_inc_tot += 0.5*dx*dy*(pow(Vx_inc(i,j),2.0)+pow(Vy_inc(i,j),2.0));
        E_com_tot += 0.5*dx*dy*(pow(Vx_com(i,j),2.0)+pow(Vy_com(i,j),2.0));
        E_qua_tot += 0.5*dx*dy*(pow(rhox(i,j),2.0)+pow(rhoy(i,j),2.0));
        E_pot_tot -= dx*dy*chem_pot*pow(abs(psi(i,j)),2.0);
        E_int_tot += dx*dy*g*0.5*pow(abs(psi(i,j)),4.0);
    }


}

for(int i=0; i < NBIN; i++){

	E_kin_spec_tot += 0.5*Lx*Ly*kin_spec(i);
	E_inc_spec_tot += 0.5*Lx*Ly*inc_spec(i);
	E_com_spec_tot += 0.5*Lx*Ly*com_spec(i);
	E_qua_spec_tot += 0.5*Lx*Ly*qua_spec(i);
}

//cout << "kinetic = " << E_kin_tot << " " << E_kin_spec_tot << endl; 
//cout << "incompr = " << E_inc_tot << " " << E_inc_spec_tot << endl; 
//cout << "compres = " << E_com_tot << " " << E_com_spec_tot << endl; 
//cout << "quantum = " << E_qua_tot << " " << E_qua_spec_tot << endl; 
//cout << "linear  = " << E_kin_tot+E_qua_spec_tot << endl;
//cout << "interaction  = " << E_int_tot << endl;

cout << "kinetic = " << E_kin_tot << " inc = " << E_inc_tot << " comp = " << E_com_tot << endl; 
cout << "quantum = " << E_qua_tot << endl; 
cout << "potential  = " << E_pot_tot << endl;
cout << "linear  = " << E_kin_tot+E_qua_spec_tot << " nonlinear  = " << E_int_tot << endl;
cout << "total = " << E_kin_tot + E_qua_tot + E_pot_tot +E_int_tot << endl;
cout << endl;
ostringstream out_energy;
     out_energy << "../../gamma_0=1em2/Output/kin_inc_com_total." << setw(5) << setfill('0') << k << ends;                        
string filenameenergy = out_energy.str();
ofstream fout_energy(filenameenergy.c_str());
fout_energy.precision(12);
           

fout_energy << E_kin_tot << " " << E_inc_tot << " " << E_com_tot << " " << E_qua_tot << " " << E_kin_tot + E_qua_tot << " " << E_pot_tot << " " << E_int_tot << " " << E_kin_tot + E_qua_tot + E_pot_tot + E_int_tot << endl;
        


ostringstream out_spec;
     out_spec << "../../gamma_0=1em2/Output/kin_inc_com_spectra." << setw(5) << setfill('0') << k << ends;                        
string filenamespec = out_spec.str();
ofstream fout_spec(filenamespec.c_str());
fout_spec.precision(12);

        for(int i=0; i < NBIN; i++){
           

                fout_spec << i*k_bin << " " << 0.5*kin_spec(i) << " " << 0.5*inc_spec(i) << " " << 0.5*com_spec(i) << " " << 0.5*qua_spec(i) << endl;
            }
         
    




}


fftw_destroy_plan(FFT);fftw_destroy_plan(IFFT);fftw_destroy_plan(FFTV);fftw_destroy_plan(FFTVM);fftw_destroy_plan(IFFTV);fftw_destroy_plan(IFFTM);
    
    return 0;
}


void embed_N_M( blitz::Array<complex<double>,2>A, blitz::Array<complex<double>,2>B){
    
    B=complex<double>(0.0,0.0);

    for(int i =0; i < Nx/2 ; i++){
        for(int j =0; j < Ny/2 ; j++){
            B(i,j) = A(i,j);
            B(Mx-i-1,j)= A(Nx-i-1,j);
            B(i,My-j-1) = A(i,Ny-j-1);
            B(Mx-i-1,My-j-1)= A(Nx-i-1,Ny-j-1);
        }
    }
    
    
    
    return;
}


void embed_N_M_real( blitz::Array<complex<double>,2>A, blitz::Array<complex<double>,2>B){
    
    B=complex<double>(0.0,0.0);

    for(int i =0; i < Nx/2 ; i++){
        for(int j =0; j < Ny/2 + 1 ; j++){
            B(i,j) = A(i,j);
            B(Mx-i-1,j)= A(Nx-i-1,j);
    
        }
    }
    
    
    
    return;
}



void embed_M_N( blitz::Array<complex<double>,2>B, blitz::Array<complex<double>,2>A){
    
    A=complex<double>(0.0,0.0);
    
    for(int i =0; i < Nx/2 ; i++){
      for(int j =0; j < Ny/2 ; j++){
        A(i,j) = B(i,j);
        A(Nx-i-1,j)= B(Mx-i-1,j);
        A(i,Ny-j-1) = B(i,My-j-1);
        A(Nx-i-1,Ny-j-1)= B(Mx-i-1,My-j-1);
      }
    }


 


    
    return;
}

void embed_M_N_real( blitz::Array<complex<double>,2>B, blitz::Array<complex<double>,2>A){
    
    A=complex<double>(0.0,0.0);
    
    for(int i =0; i < Nx/2 ; i++){
      for(int j =0; j < Ny/2 +1 ; j++){
        A(i,j) = B(i,j);
        A(Nx-i-1,j)= B(Mx-i-1,j);
      }
    }
    /*probably some symmetry */
  
    A(0,0) = complex<double>(real(A(0,0)),0.0);
    A(Nx/2,0) = complex<double>(real(A(Nx/2,0)),0.0);
    A(Nx/2,Ny/2) = complex<double>(real(A(Nx/2,Ny/2)),0.0);
    A(0,Ny/2) = complex<double>(real(A(0,Ny/2)),0.0);


       for(int i = 1; i< Nx/2; i++){
     A(Nx-i,0) = conj(A(i,0));                          // makes sure that w_hat(-kx, 0 ) = w_hat(kx,0)
     A(Nx-i,Ny/2) = conj(A(i,Ny/2));      // makes sure that w_hat(-kx,-Ny/2) = w_hat(kx, -Ny/2) 

}



    return;
}

    
void define_k(blitz::Array<complex<double>,2> kx,blitz::Array<complex<double>,2> ky,blitz::Array<complex<double>,2> k2,blitz::Array<complex<double>,2> kxr,blitz::Array<complex<double>,2> kyr,blitz::Array<complex<double>,2> k2r){


    blitz::firstIndex I;
    blitz::secondIndex J;

    kx(Range(0,Mx/2 -1), Range(0, My/2-1)) =  zip(0.0, 2.0 * pi * I / Lx, complex<double>());
    kx(Range(Mx/2 ,Mx-1), Range(0, My/2-1)) =  zip(0.0,2.0 * pi * (I-(Mx/2) ) / Lx, complex<double>());
    kx(Range(0,Mx/2 -1), Range(My/2, My-1)) =  zip(0.0,2.0 * pi * I / Lx, complex<double>());
    kx(Range(Mx/2 ,Mx-1), Range(My/2, My-1)) =  zip(0.0, 2.0 * pi * (I-(Mx/2) ) / Lx, complex<double>());
    
    ky(Range(0,Mx/2 -1), Range(0, My/2-1)) =   zip(0.0, 2.0 * pi * J / Ly, complex<double>());
    ky(Range(Mx/2 ,Mx-1), Range(0, My/2-1)) =  zip(0.0,2.0 * pi * J / Ly, complex<double>()); 
    ky(Range(0,Mx/2 -1), Range(My/2, My-1)) =   zip(0.0,2.0 * pi * (J-(My/2)) / Ly, complex<double>());
    ky(Range(Mx/2 ,Mx-1), Range(My/2, My-1)) =  zip(0.0,2.0 * pi * (J-(My/2)) / Ly, complex<double>());


    k2(Range(0,Mx/2 -1), Range(0, My/2-1)) =  zip(-pow(2.0 * pi * I / Lx ,2.0) - pow(2.0 * pi * J / Ly,2.0) , 0.0, complex<double>());
    k2(Range(Mx/2 ,Mx-1), Range(0, My/2-1)) =  zip(-pow(2.0 * pi * (I-(Mx/2) ) / Lx ,2.0) - pow(2.0 * pi * J / Ly,2.0) ,0.0, complex<double>());
    k2(Range(0,Mx/2 -1), Range(My/2, My-1)) =  zip(-pow(2.0 * pi * I / Lx ,2.0) - pow(2.0 * pi * (J-(My/2)) / Ly,2.0) ,0.0, complex<double>());
    k2(Range(Mx/2 ,Mx-1), Range(My/2, My-1)) =  zip(-pow(2.0 * pi * (I-(Mx/2) ) / Lx ,2.0) - pow(2.0 * pi * (J-(My/2)) / Ly,2.0),0.0 , complex<double>());
    k2(0,0) = complex<double>(1.0,0.0);


    kxr(Range(0,Mx/2 -1), Range(0, My/2)) =  zip(0.0, 2.0 * pi * I / Lx, complex<double>());
    kxr(Range(Mx/2 ,Mx-1), Range(0, My/2)) =  zip(0.0,2.0 * pi * (I-(Mx/2) ) / Lx, complex<double>());
 
    
    kyr(Range(0,Mx/2 -1), Range(0, My/2)) =   zip(0.0, 2.0 * pi * J / Ly, complex<double>());
    kyr(Range(Mx/2 ,Mx-1), Range(0, My/2)) =  zip(0.0,2.0 * pi * J / Ly, complex<double>()); 
    

    k2r(Range(0,Mx/2 -1), Range(0, My/2)) =  zip(-pow(2.0 * pi * I / Lx ,2.0) - pow(2.0 * pi * J / Ly,2.0) , 0.0, complex<double>());
    k2r(Range(Mx/2 ,Mx-1), Range(0, My/2)) =  zip(-pow(2.0 * pi * (I-(Mx/2) ) / Lx ,2.0) - pow(2.0 * pi * J / Ly,2.0) ,0.0, complex<double>());
    
    k2r(0,0) = complex<double>(1.0,0.0);










    return;
}

/*
void symmetry(blitz::Array<complex<double>,2> A){


    A(0,0,Range::all()) = complex<double>(0.0,0.0); //forces the zeroth wave number to be zero
    imag(A(Nx/2,0,Range::all())) = 0.0;
    imag(A(Nx/2,Ny_f-1,Range::all())) = 0.0;
    imag(A(0,Ny_f-1,Range::all())) = 0.0;

    for(int i = 1; i< Nx/2; i++){
        A(Nx-i,0,Range::all()) = conj(A(i,0,Range::all()));                             // makes sure that w_hat(-kx, 0 ) = w_h    at(kx,0)
        A(Nx-i,Ny_f-1,Range::all()) = conj(A(i,Ny_f-1,Range::all()));
    }
 



    return;
}*/

void compute_spectrum(blitz::Array<complex<double>,2> A,blitz::Array<complex<double>,2> B,blitz::Array<double,1> spectrum){

int ib =0;
double k =0.0;
double factor =0.0;
spectrum =0.0;
specnum =0;
 for(int i =0; i< Nx/2; i++){
        for(int j =0; j< Ny/2+1; j++){
                

            if(j==0 || j==Ny/2){
                factor = 1.0;
            }
            else{factor = 2.0;}
        
            k = pow( pow(   2.0*pi * i/ Lx,2.0)  + pow(   2.0*pi * j/ Ly,2.0) , 0.5);
            ib = (int) (k/k_bin);   
            spectrum(ib) +=  factor*(pow( abs(A(i,j)) , 2.0) + pow( abs(B(i,j)) , 2.0));
            specnum(ib)++;
            
            
            k = pow( pow(   2.0*pi * (-i-1)/ Lx,2.0)  + pow(   2.0*pi * j/ Ly,2.0) , 0.5);
            ib = (int) (k/k_bin);
            spectrum(ib) += factor*(pow( abs(A(Nx-i-1,j)) , 2.0)+ pow( abs(B(Nx-i-1,j)) , 2.0));
            specnum(ib)++;
           
        /*    k = pow( pow(   2.0*pi * i/ Lx,2.0)  + pow(   2.0*pi * (-j-1)/ Ly,2.0) , 0.5);
            ib = (int) (k/k_bin);
            spectrum(ib) += factor*(pow( abs(A(i,My-j-1)) , 2.0)  + pow( abs(B(i,My-j-1)) , 2.0) );
            specnum(ib)++;
            
            k = pow( pow(   2.0*pi * (-i-1)/ Lx,2.0)  + pow(   2.0*pi * (-j-1)/ Ly,2.0) , 0.5);
            ib = (int) (k/k_bin);
            spectrum(ib) += factor*( pow( abs(A(Mx-i-1,My-j-1)) , 2.0) + pow( abs(B(Mx-i-1,My-j-1)) , 2.0));
            specnum(ib)++;
          */  
        }
    }
  /*  
    for(int i=0; i < NBIN; i++){
        if(specnum(i) != 0){
            spectrum(i) /= (double) specnum(i);
       
        }
        else{
           spectrum(i)=0.0;
       
        }
      
    }
*/





}
