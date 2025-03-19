
//standalone file that creates an initial condition


#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fstream>
#include <fftw3.h>

using namespace std;
using namespace arma;

const int Nx = 1024;
const int Ny = 1024;
const int N_vortices = 2;
const double pi = 3.14159265358979323846;
const double Lx = 2.0*pi;
const double dx= Lx/ (double) Nx;
const double Ly = 2.0*pi;
const double dy= Ly / (double) Ny;

const double g = 1000.0;
const double m_eps = 1.e-8;
const double dipole_distance = pi/8.0;  // original distance pi/4.0;
const double dk = 2.0;
const double kf = 32.0;
const double Amp = 0.015;
const bool FLAG_INITIAL_TYPE = 0;
const bool FLAG_READ_IN_FROM_FILE = 0;


void define_k(mat & );
double Heaviside(double x);

int main(){
	
    cx_mat psi(Nx,Ny,fill::zeros), psi_hat(Nx,Ny,fill::zeros);
    mat k(Nx,Ny,fill::zeros);
    double theta = 0.0;
    double kk= 0.0;
    double phase = 0.0;
fftw_plan IFFT;
IFFT = fftw_plan_dft_2d(Ny,Nx, (fftw_complex*) psi_hat.memptr(), (fftw_complex*) psi.memptr(), FFTW_BACKWARD, FFTW_PATIENT);


     define_k(k);


unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
mt19937_64 generator(seed1);
uniform_real_distribution<double> distribution(0.0,2.0*pi);
psi_hat.zeros();


if(FLAG_INITIAL_TYPE==0){
for(int j =0; j < Ny/2; j++){
        for(int i =0; i < Nx/2; i++){

                if((i==0) && (j == 0)) continue;
                phase= distribution(generator);

                kk = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);

                if( abs(kk-kf) < dk){
                        psi_hat(i,j) = Amp*complex<double>(cos(phase),sin(phase));
                }
                   phase= distribution(generator);

                kk = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(j)/ Ly,2.0) , 0.5);

                if( abs(kk-kf) < dk){
                        psi_hat(Nx-i-1,j)= Amp*complex<double>(cos(phase),sin(phase));
                }
                phase= distribution(generator);
                kk = pow( pow(   2.0*pi * double(i)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);

                if( abs(kk-kf) < dk){
                        psi_hat(i,Ny-j-1) = Amp*complex<double>(cos(phase),sin(phase));
                }
                   phase= distribution(generator);
                kk = pow( pow(   2.0*pi * double(-i-1)/ Lx,2.0)  + pow(   2.0*pi * double(-j-1)/ Ly,2.0) , 0.5);

                if( abs(kk-kf) < dk){
                        psi_hat(Nx-i-1,Ny-j-1) = Amp*complex<double>(cos(phase),sin(phase));
                }
        }
}

fftw_execute(IFFT);

}
else{







 
   //pade approximate at 4th order by caliari and zuccher 2016
    double a1 = 0.34010790700196714760;
    double b1 = (2304.0*pow(a1,3.0) + 656.*pow(a1,2.0) - 421.*a1 -28.) / (7680.*pow(a1,2.0) - 1680.*a1  -330.);
    double a2 = a1*(b1 - 0.25);
    double b2 = ((737280.*pow(a1,3.0) + 209920.*pow(a1,2.0) - 134720.*a1 - 8960.)*b1 - 364544.*pow(a1,3.0) + 70144.*pow(a1,2.0)+18256.*a1+393.0)/(2457600.*pow(a1,2.0)-537600.*a1-105600.);
    double a3 = (a1*(192.*b2-48.*b1+16.*a1+5.))/(192.);
    double b3 = ((61440.*a1-9600.)*b2 + (-30720.*pow(a1,2.0)+640.*a1 + 560.)*b1 + 8448.*pow(a1,2.0) - 1056.*a1 -21.) / (368640.*a1 - 92160.);
    double a4 = (4608.*a1*b3 - 1152.*a1*b2 + b1*(384.*pow(a1,2.0)+120.*a1) - 128.*pow(a1,2.0) - 7.*a1)/4608.;
    double r = 0.0;
    double pos_x = 0.0;
    double pos_y = 0.0;
    double rho = 0.0;
    double dtheta = 0.0;
    

    double gamma = 0.0;
    double xdiff = 0.0;
    double ydiff = 0.0;
    int m;
  


rowvec position_x(N_vortices,fill::zeros), position_y(N_vortices,fill::zeros), circulation(N_vortices,fill::zeros);


psi.fill(complex<double>(1.0,0.0));

if(FLAG_READ_IN_FROM_FILE ==1){
	
    fstream filein("./vortex_xy.initial");
    filein.precision(12);
    

 for(int k = 0; k < N_vortices; k++){
   filein >> position_x(k) >> position_y(k) >> circulation(k);
    }

}
else{
       //defines a dipole state
        position_x(0) = pi;
        position_y(0) = pi - 0.5*dipole_distance;
        circulation(0) = 1.0;

        position_x(1) = pi;
        position_y(1) = pi + 0.5*dipole_distance;
        circulation(1) = -1.0;


}

    for(int k = 0; k < N_vortices; k++){


            pos_x = position_x(k);
            pos_y = position_y(k);
            gamma = copysign(1.0,circulation(k));

       

        theta = 0.0;

        for(int j=0; j < Ny; j++){
        	for(int i=0; i < Nx; i++){
          
                m=1;
                dtheta = 1.e0;
                xdiff = dx*double(i) - pos_x;
                ydiff =  dy*double(j) - pos_y;

                theta = -gamma*atan( tanh(0.5*ydiff)*tan(0.5*(xdiff - pi))) + gamma*pi*Heaviside(xdiff);
               
               while(dtheta > m_eps){
            
                    dtheta = -gamma*atan( tanh(0.5*(ydiff + 2.0*pi*double(m)))*tan(0.5*(xdiff - pi)))  - gamma*atan( tanh(0.5*(ydiff - 2.0*pi*double(m)))*tan(0.5*(xdiff - pi)));
                    
                    theta +=  2.0*gamma*pi*Heaviside(xdiff);
                  
                    theta += dtheta;
                    m++;
                }
                
                theta -= 0.5*gamma*double(j)*dy*xdiff/pi;

                //theta += i*dx * Heaviside(abs(ydiff) - pi) + j*dy * Heaviside(abs(xdiff) - pi);

                r = pow(pow(dx*i - pos_x ,2.0) + pow(dy*j - pos_y ,2.0),0.5) ;        
                //rho = g*r*r*(a1+a2*g*r*r) / (1.0 + b1*g*r*r + b2*g*g*r*r*r*r);
                r *= (sqrt(g));
                rho = (a1*pow(r,2.0) + a2*pow(r,4.0) + a3*pow(r,6.0) + a4*pow(r,8.0))/(1+b1*pow(r,2.0)+b2*pow(r,4.0)+ b3*pow(r,6.0)+a4*pow(r,8.0));
                
                psi(i,j) *= rho*exp(complex<double>(0.0,theta));
                        
            }
        }
    }
}
    ofstream fout("./psi.initial");				//open up file
	fout.precision(12);
    for(int i=0; i < Nx; i++){
        for(int j=0; j < Ny; j++){
            fout << i*dx << " " << j*dy << " " << real(psi(i,j)) << " " << imag(psi(i,j)) << endl;
       
        }
        fout << endl;
    }

    fout.close();

fftw_destroy_plan(IFFT);
    return 0;

}


double Heaviside(double x){
    return copysign(0.5,x)+0.5;
}


void define_k(mat & k){



  for(int j =0; j < Ny/2 ; j++){
        for(int i =0; i < Nx/2 ; i++){

        	k(i,j) = pow(  pow( 2.0*pi*double(i) / Lx, 2.0) + pow( 2.0*pi*double(j) / Ly, 2.0)    ,0.5) ; 
        	k(i,Ny-j-1) = pow(  pow( 2.0*pi*double(i) / Lx, 2.0) + pow( 2.0*pi*double(-j-1) / Ly, 2.0)    ,0.5) ; 
        	k(Nx-i-1,j) = pow(  pow(2.0*pi*double(-i-1) / Lx , 2.0) + pow( 2.0*pi*double(j) / Ly, 2.0)    ,0.5) ; 
        	k(Nx-i-1,Ny-j-1) = pow(  pow( 2.0*pi*double(-i-1) / Lx, 2.0) + pow(2.0*pi*double(-j-1) / Ly , 2.0)    ,0.5) ; 

 }
}
}
