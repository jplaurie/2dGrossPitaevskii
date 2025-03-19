
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

const int Nx=256;
const int Ny=256;
const int Mx = 3*Nx/2;                     
const int My = 3*Ny/2;                      
const double pi = 3.14159265358979323846;       //pi
const double Lx = 2.0*pi;                   //length of the box
const double dx = Lx/ double(Nx);                      //grid size
const double Ly = 2.0*pi;                   //length of the box
const double dy = Ly/ double(Ny);  
const int startfile = 0;
const int endfile = 0;
const int eta = 1;

const double epss = 0.01;
const double pos_eps = 1.e-6;



void embed_N_M(cx_mat & , cx_mat &);
void embed_M_N(cx_mat &,cx_mat &);
void embed_N_M_real(cx_mat & , cx_mat &);
void embed_M_N_real(cx_mat &,cx_mat &);
void define_k(cx_mat &,cx_mat &);
double LineInt(mat &,mat &, int, int);
void dealias_M(cx_mat &);
void symmetry(cx_mat & );
complex<double> reconstruct(cx_mat, rowvec);
rowvec find_vortex_position(rowvec,cx_mat, cx_mat, cx_mat,double);

int main(){
    

mat velx(Nx,Ny,fill::zeros), vely(Nx,Ny,fill::zeros), pseudo_vorticity(Nx,Ny,fill::zeros), real(Nx,Ny,fill::zeros), real_M(Mx,My,fill::zeros);
cx_mat psi(Nx,Ny,fill::zeros),psi_hat(Nx,Ny,fill::zeros);

cx_mat temp(Nx,Ny,fill::zeros), temp_M(Mx,My,fill::zeros), temp_hat(Nx,Ny,fill::zeros), temp_hat_M(Mx,My,fill::zeros);
cx_mat real_hat(Nx/2+1,Ny,fill::zeros), real_hat_M(Mx/2+1,My, fill::zeros);


cx_mat ikx(Mx,My,fill::zeros),iky(Mx,My,fill::zeros),psi_M(Mx,My,fill::zeros), psi_hat_M(Mx,My,fill::zeros),psi_x(Mx,My,fill::zeros),psi_y(Mx,My,fill::zeros)   ;
 

    double real_in,imag_in,ignore;
   


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
	laurie << "../src/data_dipole/psi."<<  setw(5) << setfill('0') << k << ends;
	
	string filename = laurie.str();
	ifstream data_file(filename.c_str());
	data_file.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
 
        data_file >> ignore >> ignore >> real_in >> imag_in ;
            psi(i,j) = complex<double>(real_in,imag_in);
    
      }
    }

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
real_M =  arma::imag(conj(temp_M)%psi_x);//real(complex<double>(0,-0.5)* (conj(temp_M)%psi_x - temp_M%conj(psi_x)));   //computes x-velocity in physical space
fftw_execute(FFTM_real);            //fourier transform velocity
real_hat_M /= double(Mx * My);     //normalize
embed_M_N_real(real_hat_M,real_hat);     //embeds velocity to smaller array
fftw_execute(IFFT_real);     //inverse fourier transform 
velx=real;     //record

//embed_N_M(temp_hat,temp_hat_M);      //embeds data to larger array
//fftw_execute(IFFTM);     
real_M=  arma::imag(conj(temp_M)%psi_y);//real(complex<double>(0,-0.5)* (conj(psi_M)%psi_y - psi_M%conj(psi_y))) ;//computes y velocity
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










double max_psi = abs(psi).max();

cx_mat psi_hat(Nx,Ny,fill::zeros),psi_x_hat(Nx,Ny,fill::zeros),psi_y_hat(Nx,Ny,fill::zeros);
mat vortex_position(Nx*Ny,2,fill::zeros);
rowvec circ(Nx*Ny,fill::zeros);
rowvec pos(2,fill::zeros);
double counter = 0.0;


temp = psi;
fftw_execute(FFT); 
temp_hat /= double(Nx * Ny); 
psi_hat=temp_hat;

embed_N_M(temp_hat,temp_hat_M);
temp_hat_M %= ikx;                // computes x derivative in fourier space
embed_M_N(temp_hat_M,psi_x_hat);
temp_hat_M %= iky;                // computes x derivative in fourier space
embed_M_N(temp_hat_M,psi_y_hat);


for(int i=0; i< Nx; i++){
	for(int j=0; j< Ny; j++){

		if( std::abs(psi(i,j)) < 0.1*max_psi){
			pos(0) = i*dx;
			pos(1) = j*dy;
			
			pos =  find_vortex_position(pos,psi_hat, psi_x_hat, psi_y_hat,max_psi);
		
			vortex_position(counter,0) = pos(0);
			vortex_position(counter,1) = pos(1);
			circ(counter) = copysign(1.0,pseudo_vorticity(i,j));
			counter++;
         }

	}
}





ostringstream out_vortexpos;
     out_vortexpos << "./vortex_location." << setw(5) << setfill('0') << k << ends;                        
string filenamep = out_vortexpos.str();
ofstream fout_vortexpos(filenamep.c_str());
fout_vortexpos.precision(12);

      

            for(int m=0; m < counter; m++){  
                if(circ(m)!= 0){
                     fout_vortexpos << vortex_position(m,0) << " " << vortex_position(m,1) << " " << circ(m) << endl;
                }
            }
           
        
     fout_vortexpos.close();




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

double LineInt(mat & Vx, mat & Vy, int i, int j){

	double l1,l2,l3,l4;
	int jj,ii;
	l1=0.0;l2=0.0;l3=0.0;l4=0.0;

	jj = j + eta;
	if(jj < 0) jj += Ny;
	else if(jj > Ny-1) jj -= Ny;

	for(int i2 = i; i2 < i+eta+1; i2++){

        ii=i2;
		if(ii < 0) ii += Nx;
		else if(ii > Nx-1) ii -= Nx;
		
		l1 -= Vx(ii,jj);
	}

	jj = j;//- eta;
	if(jj < 0) jj += Ny;
	else if(jj > Ny-1) jj -= Ny;

	for(int i2 = i; i2 < i+eta+1; i2++){
        ii=i2;
		if(ii < 0) ii +=Nx;
		else if(ii > Nx-1) ii -= Nx;
		
		l3 += Vx(ii,jj);
	}

	ii = i + eta;
	if(ii < 0) ii += Nx;
	else if(ii > Nx-1) ii -= Nx;

	for(int j2 = j; j2 < j+eta+1; j2++){

        jj=j2;
		if(jj < 0) jj +=Ny;
		else if(jj > Ny-1) jj -= Ny;
		
		l2 += Vy(ii,jj);
	}

	ii = i ;//- eta;
	if(ii < 0) ii += Nx;
	else if(ii > Nx-1) ii -= Nx;

	for(int j2 = j; j2 < j+eta+1; j2++){

        jj=j2;
		if(jj < 0) jj +=Ny;
		else if(jj > Ny-1) jj -= Ny;
		
		l4 -= Vy(ii,jj);
	}


	return l1+l2+l3+l4;

}




rowvec find_vortex_position(rowvec pos, cx_mat psi_hat, cx_mat psi_x_hat, cx_mat psi_y_hat, double max_psi){

	mat jacobian(2,2,fill::zeros),inv_jacobian(2,2,fill::zeros);
	rowvec newpos(2,fill::zeros), psi_vec(2,fill::zeros);
	complex<double> temp_value=complex<double>(0.0,0.0);
 	double diff = 1.e6;
 	newpos = pos;
	
	while( diff > 1.e-6*max_psi){
	
	pos = newpos;

	temp_value = reconstruct(psi_x_hat,pos);


	jacobian(0,0) = std::real(temp_value);
	jacobian(0,1) = std::imag(temp_value);

	temp_value = reconstruct(psi_y_hat,pos);

	jacobian(1,0) = std::real(temp_value);
	jacobian(1,1) = std::imag(temp_value);

	inv_jacobian = inv(jacobian);

	temp_value = reconstruct(psi_hat,pos);
	
	cout << "psi =  " << temp_value << endl;

	newpos(0) = pos(0) - inv_jacobian(0,0)*std::real(temp_value) - inv_jacobian(0,1)*std::imag(temp_value);
	newpos(1) = pos(1) - inv_jacobian(1,0)*std::real(temp_value) - inv_jacobian(1,1)*std::imag(temp_value);


	//diff = sum( pow(newpos - pos, 2.0));
	diff = abs(temp_value);
	cout << pos << " " << newpos << "  " << diff <<   endl;



}

return newpos;
}



complex<double> reconstruct(cx_mat A_hat, rowvec x){

complex<double> value=complex<double>(0.0,0.0);
double kx, ky;

for(int i=0; i< Nx/2; i++){
	for(int j=0; j< Ny/2; j++){
		kx = 2.0*pi*double(i)/Lx;
		ky = 2.0*pi*double(j)/Ly;
		value += A_hat(i,j)*exp(complex<double>(0.0,1.0)*(kx*x(0)+ky*x(1)));

		kx = 2.0*pi*double(-i-1)/Lx;
		ky = 2.0*pi*double(j)/Ly;
		value += A_hat(Nx-i-1,j)*exp(complex<double>(0.0,1.0)*(kx*x(0)+ky*x(1)));

		kx = 2.0*pi*double(i)/Lx;
		ky = 2.0*pi*double(-j-1)/Ly;
		value += A_hat(i,Ny-j-1)*exp(complex<double>(0.0,1.0)*(kx*x(0)+ky*x(1)));

		kx = 2.0*pi*double(-i-1)/Lx;
		ky = 2.0*pi*double(-j-1)/Ly;
		value += A_hat(Nx-i-1,Ny-j-1)*exp(complex<double>(0.0,1.0)*(kx*x(0)+ky*x(1)));
		}
}

return value;

}


