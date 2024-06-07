
// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace arma;

const int Nx=512;
const int Ny=512;
const int Mx = 3*Nx/2;                     
const int My = 3*Ny/2;                      
const double pi = 3.14159265358979323846;       //pi
const double Lx = 2.0*pi;                   //length of the box
const double dx = Lx/ double(Nx);                      //grid size
const double Ly = 2.0*pi;                   //length of the box
const double dy = Ly/ double(Ny);  
const int startfile = 1;
const int endfile = 310;
const int eta = 9;


void embed_N_M(cx_mat & , cx_mat &);
void embed_M_N(cx_mat &,cx_mat &);
void embed_N_M_real(cx_mat & , cx_mat &);
void embed_M_N_real(cx_mat &,cx_mat &);
void define_k(cx_mat &,cx_mat &);
void dealias_M(cx_mat &);
void symmetry(cx_mat & );
//complex<double> reconstruct(cx_mat, rowvec);
//rowvec find_vortex_position(rowvec,cx_mat, cx_mat, cx_mat,double);
bool is_vortex(int , int , mat );

int main(){
    

mat velx(Nx,Ny,fill::zeros), vely(Nx,Ny,fill::zeros), pseudo_vorticity(Nx,Ny,fill::zeros), real(Nx,Ny,fill::zeros), real_M(Mx,My,fill::zeros), circulation(Nx,Ny,fill::zeros);
cx_mat psi(Nx,Ny,fill::zeros),psi_hat(Nx,Ny,fill::zeros);

cx_mat temp(Nx,Ny,fill::zeros), temp_M(Mx,My,fill::zeros), temp_hat(Nx,Ny,fill::zeros), temp_hat_M(Mx,My,fill::zeros);
cx_mat real_hat(Nx/2+1,Ny,fill::zeros), real_hat_M(Mx/2+1,My, fill::zeros);


cx_mat ikx(Mx,My,fill::zeros),iky(Mx,My,fill::zeros),psi_M(Mx,My,fill::zeros), psi_hat_M(Mx,My,fill::zeros),psi_x(Mx,My,fill::zeros),psi_y(Mx,My,fill::zeros)   ;
 

    double real_in,imag_in,ignore;

    mat positive_vortex(Nx*Ny,2,fill::zeros),negative_vortex(Nx*Ny,2,fill::zeros);
rowvec circ(Nx*Ny,fill::zeros);

int pos_counter=0;
int neg_counter=0;
bool test;




    fftw_plan FFT;
    fftw_plan IFFT;
    fftw_plan FFTM;
    fftw_plan IFFTM;
    
    fftw_plan FFT_real;
    fftw_plan IFFT_real;
    fftw_plan FFTM_real;
    fftw_plan IFFTM_real;

    FFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) temp.memptr(), (fftw_complex*) temp_hat.memptr(), FFTW_FORWARD, FFTW_PATIENT);
    IFFT = fftw_plan_dft_2d(Nx,Ny, (fftw_complex*) temp_hat.memptr(), (fftw_complex*) temp.memptr(), FFTW_BACKWARD, FFTW_PATIENT);
    
    FFTM = fftw_plan_dft_2d(Mx,My, (fftw_complex*) temp_M.memptr(), (fftw_complex*) temp_hat_M.memptr(), FFTW_FORWARD, FFTW_PATIENT);
    IFFTM = fftw_plan_dft_2d(Mx,My, (fftw_complex*) temp_hat_M.memptr(), (fftw_complex*) temp_M.memptr(), FFTW_BACKWARD, FFTW_PATIENT);
    


    FFTM_real = fftw_plan_dft_r2c_2d(Mx,My, (double*) real_M.memptr(), (fftw_complex*) real_hat_M.memptr(), FFTW_PATIENT);
    FFT_real = fftw_plan_dft_r2c_2d(Nx,Ny, (double*) real.memptr(), (fftw_complex*) real_hat.memptr(), FFTW_PATIENT); 
    IFFTM_real = fftw_plan_dft_c2r_2d(Mx,My,  (fftw_complex*) real_hat_M.memptr(),(double*) real_M.memptr(), FFTW_PATIENT);
    IFFT_real = fftw_plan_dft_c2r_2d(Nx,Ny,  (fftw_complex*) real_hat.memptr(),(double*) real.memptr(), FFTW_PATIENT);

   
define_k(ikx, iky);


for(int k=startfile; k <= endfile; k++){

    cout << "file = " << k << endl;
	ostringstream laurie;
	laurie << "../data/psi."<<  setw(5) << setfill('0') << k << ends;
	
	string filename = laurie.str();
	ifstream data_file(filename.c_str());
	data_file.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
 
        data_file >> ignore >> ignore >> real_in >> imag_in ;
            psi(i,j) = complex<double>(real_in,imag_in);
    
      }
    }
pos_counter=0;
neg_counter=0;


temp = psi;  	
fftw_execute(FFT);  //transforms original data to Fourier space
temp_hat /= double(Nx * Ny); // normalizes
embed_N_M(temp_hat,temp_hat_M);
temp_hat_M %= ikx;                // computes x derivative in fourier space
fftw_execute(IFFTM);            // inverse fourier transform
psi_x = temp_M;                     //record x-derivative

embed_N_M(temp_hat,temp_hat_M);
temp_hat_M %= iky;            //computes y-derivative  
fftw_execute(IFFTM);  //inverse fourier transform to physical space
psi_y = temp_M;          //records to file

embed_N_M(temp_hat,temp_hat_M);      //embeds data to larger array
fftw_execute(IFFTM);            // inverse fourier transform to larger array
real_M =  arma::imag(conj(temp_M)%psi_x)/pow(abs(temp_M),2.0) ;//real(complex<double>(0,-0.5)* (conj(temp_M)%psi_x - temp_M%conj(psi_x)));   //computes x-velocity in physical space
fftw_execute(FFTM_real);            //fourier transform velocity
real_hat_M /= double(Mx * My);     //normalize
embed_M_N_real(real_hat_M,real_hat);     //embeds velocity to smaller array
fftw_execute(IFFT_real);     //inverse fourier transform 
velx=real;     //record

//embed_N_M(temp_hat,temp_hat_M);      //embeds data to larger array
//fftw_execute(IFFTM);     
real_M=  arma::imag(conj(temp_M)%psi_y)/pow(abs(temp_M),2.0) ;//real(complex<double>(0,-0.5)* (conj(psi_M)%psi_y - psi_M%conj(psi_y))) ;//computes y velocity
fftw_execute(FFTM_real);            //fourier transform velocity
real_hat_M /= double(Mx * My);     //normalize
embed_M_N_real(real_hat_M,real_hat);     //embeds velocity to smaller array
fftw_execute(IFFT_real);     //inverse fourier transform 
vely=real;     //record


real_M = arma::imag( conj(psi_x)%psi_y);

fftw_execute(FFTM_real);            //fourier transform velocity
dealias_M(real_hat_M);
real_hat_M /= double(Mx * My);     //normalize
embed_M_N_real(real_hat_M,real_hat);     //embeds velocity to smaller array
fftw_execute(IFFT_real);     //inverse fourier transfo
pseudo_vorticity=real;


pseudo_vorticity /= pow(abs(psi),2.0);

ostringstream out_V;
 out_V << "./vel_vort." << setw(5) << setfill('0') << k << ends;                        
string filenameV = out_V.str();
ofstream fout_V(filenameV.c_str());
fout_V.precision(12);

        for(int i=0; i < Nx; i++){
            for(int j=0; j < Ny; j++){

                fout_V << i*dx << " " << j*dy << " " << velx(i,j) << " " << vely(i,j) << " " << pseudo_vorticity(i,j) << endl;
            }
            fout_V << endl;
        }
     fout_V.close();




double max_vor = pseudo_vorticity.max();
double min_vor = pseudo_vorticity.min();

for(int i =0; i< Nx; i++){
    for(int j=0; j< Ny; j++){

        test = is_vortex(i,j,pseudo_vorticity);

        if(test == true){
            if(pseudo_vorticity(i,j)>0){
                positive_vortex(pos_counter,0)= i*dx;
                positive_vortex(pos_counter,1)= j*dy;
                pos_counter++;
            }
            else{
                negative_vortex(neg_counter,0)= i*dx;
                negative_vortex(neg_counter,1)= j*dy;
                neg_counter++;
            }
        }

    }
}




ostringstream out_positive_vortexpos;
out_positive_vortexpos << "./pos_vortices." << setw(5) << setfill('0') << k << ends;                        
string filenamep = out_positive_vortexpos.str();
ofstream fout_positive_vortexpos(filenamep.c_str());
fout_positive_vortexpos.precision(12);


 for(int i=0; i < pos_counter; i++){  
    fout_positive_vortexpos << k << " " << positive_vortex(i,0) << " " << positive_vortex(i,1) << " " << 1.0 << endl;
}
            
fout_positive_vortexpos.close();

ostringstream out_negative_vortexpos;
out_negative_vortexpos << "./neg_vortices." << setw(5) << setfill('0') << k << ends;                        
string filenamen = out_negative_vortexpos.str();
ofstream fout_negative_vortexpos(filenamen.c_str());
fout_negative_vortexpos.precision(12);


 for(int i=0; i < neg_counter; i++){  
    fout_negative_vortexpos << k << " " << negative_vortex(i,0) << " " << negative_vortex(i,1) << " " << -1.0 << endl;
}
            
fout_negative_vortexpos.close();



}




fftw_destroy_plan(FFT);fftw_destroy_plan(IFFT);fftw_destroy_plan(FFTM);fftw_destroy_plan(IFFTM);
fftw_destroy_plan(FFT_real);fftw_destroy_plan(IFFT_real);fftw_destroy_plan(FFTM_real);fftw_destroy_plan(IFFTM_real);



    return 0;
}



void embed_N_M( cx_mat & A, cx_mat & B){
    
    B.zeros();
    
    
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


void embed_M_N( cx_mat & B, cx_mat & A){
    
    A.zeros();
    
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

void embed_N_M_real(  cx_mat & A, cx_mat & B){
    
    B.zeros();

    for(int i =0; i < Nx/2+1 ; i++){
        for(int j =0; j < Ny/2  ; j++){
            B(i,j) = A(i,j);
            B(i,My-j-1)= A(i,Ny-j-1);
    
        }
    }
    
    
    
    return;
}




void embed_M_N_real( cx_mat & B , cx_mat & A){
    
    A.zeros();
    
    for(int i =0; i < Nx/2+1 ; i++){
      for(int j =0; j < Ny/2 ; j++){
      		A(i,j) = B(i,j);
        	A(i,Ny-j-1)= B(i,My-j-1);
      }
    }
    symmetry(A);
    return;
}

void symmetry(cx_mat & A){
  

  	A(0,0) = complex<double>(real(A(0,0)),0.0);
    A(0,Ny/2) = complex<double>(real(A(0,Ny/2)),0.0);
    A(Nx/2,Ny/2) = complex<double>(real(A(Nx/2,Ny/2)),0.0);
    A(Nx/2,0) = complex<double>(real(A(Nx/2,0)),0.0);


      

    for(int j = 1; j< Ny/2; j++){
     A(0,Ny-j) = conj(A(0,j));                          // makes sure that w_hat(-kx, 0 ) = w_hat(kx,0)
     A(Nx/2,Ny-j) = conj(A(Nx/2,j)); 
     }     // makes sure that w_hat(-kx,-Ny/2) = w_hat(kx, -Ny/2) 
return;

}


void dealias_M( cx_mat & A){

for(int i =Nx/2 +1; i < Mx/2 +1 ; i++){
      for(int j =Ny/2; j < My/2  ; j++){
            A(i,j) = complex<double>(0.0,0.0);
            A(0,Mx-j-1)= complex<double>(0.0,0.0);
      }
}
    return;
}


    
void define_k(cx_mat & ikx,cx_mat & iky){


  for(int j =0; j < My/2 ; j++){
        for(int i =0; i < Mx/2 ; i++){
        	ikx(i,j) = complex<double>(0.0, 2.0*pi*double(i) / Lx);
        	ikx(i,My-j-1) = complex<double>(0.0, 2.0*pi*double(i) / Lx);
        	ikx(Mx-i-1,j) = complex<double>(0.0, 2.0*pi*double(-i-1) / Lx);
        	ikx(Mx-i-1,My-j-1) = complex<double>(0.0, 2.0*pi*double(-i-1) / Lx);

        	iky(i,j) = complex<double>(0.0, 2.0*pi*double(j) / Ly);
        	iky(i,My-j-1) = complex<double>(0.0, 2.0*pi*double(-j-1) / Ly);
        	iky(Mx-i-1,j) = complex<double>(0.0, 2.0*pi*double(j) / Ly);
        	iky(Mx-i-1,My-j-1) = complex<double>(0.0, 2.0*pi*double(-j-1) / Ly);


}
}



    return;
}



bool is_vortex(int i, int j, mat pseudo_vorticity){


static double max_vort;
static double min_vort;
static double test_vort;
static int ii, jj;
test_vort = pseudo_vorticity(i,j);
max_vort = pseudo_vorticity.max();
min_vort = pseudo_vorticity.min();
if(test_vort > 0.8*max_vort){

    for(int i2 = i -eta; i2 < i+eta+1; i2++){
        if(i2<0) ii=i2+Nx;
        else if(i2>Nx-1) ii = i2-Nx;
        else ii=i2;
            for(int j2 = j -eta; j2 < j+eta+1; j2++){
                if(j2<0) jj=j2+Ny;
                else if(j2>Ny-1) jj = j2-Ny;
                else jj=j2;
                 if(i2==i && j2==j){
                    continue;
                 }
                 else if(pseudo_vorticity(ii,jj) > test_vort){
                    return false;
                }
            }

        }
    return true;
}
else if (test_vort < 0.8*min_vort){

    for(int i2 = i -eta; i2 < i+eta+1; i2++){
        if(i2<0) ii=i2+Nx;
        else if(i2>Nx-1) ii = i2-Nx;
        else ii=i2;
            for(int j2 = j -eta; j2 < j+eta+1; j2++){
                if(j2<0) jj=j2+Ny;
                else if(j2>Ny-1) jj = j2-Ny;
                else jj=j2;
                 if(i2==i && j2==j){
                    continue;
                 }
                 else if(pseudo_vorticity(ii,jj) < test_vort){
                    return false;
                }
            }

        }
return true;
}
else{
    return false;
}
}

