
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
const int startfile = 900;
const int endfile = 900;



void embed_N_M(blitz::Array<complex<double>,2> , blitz::Array<complex<double>,2>);
void embed_M_N(blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>);
void embed_N_M_real(blitz::Array<complex<double>,2> , blitz::Array<complex<double>,2>);
void embed_M_N_real(blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>);
void define_k(blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>,blitz::Array<complex<double>,2>);
int main(){
    

blitz::Array<double,2> Vx(Nx,Ny), Vy(Nx,Ny), W(Nx,Ny),V_M(Mx,My),V(Nx,Ny);
blitz::Array<complex<double>,2> psi(Nx,Ny),psi_hat(Nx,Ny),kx(Mx,My),ky(Mx,My),kxr(Mx,My/2+1),kyr(Mx,My/2+1), psi_M(Mx,My), psi_hat_M(Mx,My) , V_hat_M(Mx,My/2+1), psi_x(Mx,My), psi_y(Mx,My),W_hat_M(Mx,My/2+1), V_hat(Nx,Ny/2+1);

   

    double real_in,imag_in,ignore;
   


    fftw_plan FFT;
    fftw_plan IFFT;
    fftw_plan IFFTM;
    fftw_plan FFTVM;
    fftw_plan FFTV;
    fftw_plan IFFTVM;
    fftw_plan IFFTV;
  

    FFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) psi.data(), (fftw_complex*) psi_hat.data(), FFTW_FORWARD, FFTW_PATIENT);
    IFFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) psi_hat.data(), (fftw_complex*) psi.data(), FFTW_BACKWARD, FFTW_PATIENT);
    IFFTM = fftw_plan_dft_2d(Mx,My, (fftw_complex*) psi_hat_M.data(), (fftw_complex*) psi_M.data(), FFTW_BACKWARD, FFTW_PATIENT);
    FFTVM = fftw_plan_dft_r2c_2d(Mx,My, (double*) V_M.data(), (fftw_complex*) V_hat_M.data(), FFTW_PATIENT);
    FFTV = fftw_plan_dft_r2c_2d(Nx,Ny, (double*) V.data(), (fftw_complex*) V_hat.data(), FFTW_PATIENT); 
    IFFTVM = fftw_plan_dft_c2r_2d(Mx,My,  (fftw_complex*) V_hat_M.data(),(double*) V_M.data(), FFTW_PATIENT);
    IFFTV = fftw_plan_dft_c2r_2d(Nx,Ny,  (fftw_complex*) V_hat.data(),(double*) V.data(), FFTW_PATIENT);


    Vx=0.0;
    Vy=0.0;
    V=0.0;
    V_M=0.0;


    W=0.0;
    kx=complex<double>(0.0,0.0);
    ky=complex<double>(0.0,0.0);
    psi=complex<double>(0.0,0.0);
    psi_M=complex<double>(0.0,0.0);
    psi_hat=complex<double>(0.0,0.0);
    psi_hat_M=complex<double>(0.0,0.0);
    psi_x=complex<double>(0.0,0.0);
    psi_y=complex<double>(0.0,0.0);

    V_hat=complex<double>(0.0,0.0);
    V_hat_M=complex<double>(0.0,0.0);
    W_hat_M=complex<double>(0.0,0.0);

   
define_k(kx, ky,kxr,kyr);


for(int k=startfile; k <= endfile; k++){

    cout << "file = " << k << endl;
	ostringstream laurie;
	laurie << "../../gamma_0=1em3/data/psi."<<  setw(5) << setfill('0') << k << ".dat" << ends;
	
	string filename = laurie.str();
	ifstream data_file(filename.c_str());
	data_file.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
 
        data_file >> ignore >> ignore >> ignore >> real_in >> imag_in ;
            psi(i,j) = complex<double>(real_in,imag_in);
           // cout << i << " " << j << " " << psi(i,j) << endl;
      }
    }

   	
fftw_execute(FFT);  //transforms original data to Fourier space

psi_hat /= (double) Nx *Ny; // normalizes

embed_N_M(psi_hat,psi_hat_M);

psi_hat_M *= kx;                // computes x derivative in fourier space

fftw_execute(IFFTM);            // inverse fourier transform

 psi_x = psi_M;                     //record x-derivative

embed_N_M(psi_hat,psi_hat_M);

psi_hat_M *= ky;            //computes y-derivative
    
fftw_execute(IFFTM);  //inverse fourier transform to physical space

psi_y = psi_M;          //records to file

embed_N_M(psi_hat,psi_hat_M);      //embeds data to larger array

fftw_execute(IFFTM);            // inverse fourier transform to larger array

V_M=  real(complex<double>(0,-0.5)* (conj(psi_M)*psi_x - psi_M*conj(psi_x))/ pow(abs(psi_M),2.0) );  //computes x-velocity in physical space

fftw_execute(FFTVM);            //fourier transform velocity

V_hat_M /= (double) Mx * My;     //normalize

W_hat_M = -V_hat_M * kyr;        //computes y-derivative of velocity and saves to w_hat

embed_M_N_real(V_hat_M,V_hat);     //embeds velocity to smaller array


fftw_execute(IFFTV);     //inverse fourier transform 

Vx=V;     //record

V_M=  real(complex<double>(0,-0.5)* (conj(psi_M)*psi_y - psi_M*conj(psi_y))/ pow(abs(psi_M),2.0) );//computes y velocity

fftw_execute(FFTVM);

V_hat_M /= (double) Mx * My;     //normalize

W_hat_M += V_hat_M * kxr;        //this gives the vorticity

embed_M_N_real(V_hat_M,V_hat);     //embed x velocity to smaller array

fftw_execute(IFFTV);         //inverse fourier transform

Vy=V;         //save 

embed_M_N(W_hat_M,V_hat);     //embed vorticity to smaller array

fftw_execute(IFFTV);     //inverse fourier transform

W=V;    //save


ostringstream out_V;
     out_V << "../../gamma_0=1em3/Output/Velocity." << setw(5) << setfill('0') << k << ends;                        
string filenameV = out_V.str();
ofstream fout_V(filenameV.c_str());
fout_V.precision(12);

        for(int i=0; i < Nx; i++){
            for(int j=0; j < Ny; j++){

                fout_V << i*dx << " " << j*dy << " " << Vx(i,j) << " " << Vy(i,j) << endl;
            }
            fout_V << endl;
        }
     fout_V.close();


ostringstream out_W;
     out_W << "../../gamma_0=1em3/Output/w." << setw(5) << setfill('0') << k << ends;                        
string filenameW = out_W.str();
ofstream fout_W(filenameW.c_str());
fout_W.precision(12);

        for(int i=0; i < Nx; i++){
            for(int j=0; j < Ny; j++){

                fout_W << i*dx << " " << j*dy << " " << W(i,j) << endl;
            }
            fout_W << endl;
        }
     fout_W.close();

}


fftw_destroy_plan(FFT);fftw_destroy_plan(IFFT);fftw_destroy_plan(FFTV);fftw_destroy_plan(IFFTM);fftw_destroy_plan(FFTVM);fftw_destroy_plan(IFFTVM);
    
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
    
void define_k(blitz::Array<complex<double>,2> kx,blitz::Array<complex<double>,2> ky,blitz::Array<complex<double>,2> kxr,blitz::Array<complex<double>,2> kyr){


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


    kxr(Range(0,Mx/2 -1), Range(0, My/2)) =  zip(0.0, 2.0 * pi * I / Lx, complex<double>());
    kxr(Range(Mx/2 ,Mx-1), Range(0, My/2)) =  zip(0.0,2.0 * pi * (I-(Mx/2) ) / Lx, complex<double>());
   
    
    kyr(Range(0,Mx/2 -1), Range(0, My/2)) =   zip(0.0, 2.0 * pi * J / Ly, complex<double>());
    kyr(Range(Mx/2 ,Mx-1), Range(0, My/2)) =  zip(0.0,2.0 * pi * J / Ly, complex<double>()); 



    return;
}

