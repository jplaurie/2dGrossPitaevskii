
// This compute the linear and nonlinear energies, and total wave action


#include <iostream>
#include <cmath>
#include <blitz/array.h>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <random>

using namespace std;
using namespace blitz;

const int Nx=512;
const int Ny=512;
                   
const double pi = 3.14159265358979323846;       //pi
const double Lx = 32.0*pi;                   //length of the box
const double dx = Lx/ (double) Nx;                      //grid size
const double Ly = 32.0*pi;                   //length of the box
const double dy = Ly/ (double) Ny;  

const int N_avg = 10000;

const int startfile = 900;
const int endfile = 900;


const int rtot = 128;

static blitz::Array<double,1> cir_r(rtot), cir2_r(rtot);

static blitz::Array<double,2> w(Nx,Ny);
int main(){
 

 int i2,j2;   
w=0.0;
cir_r = 0.0;
cir2_r = 0.0;

double circ, circ2,dr;
circ=0.0;
circ2=0.0;

double real_in,ignore;
int i_center, j_center;


random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_int_distribution<int> distribution(0, Nx-1);

for(int k=startfile; k <= endfile; k++){        //increment over file number


    cout << "file = " << k << endl;

	ostringstream laurie;
	laurie << "../../gamma_0=1em3/Output/w."<<  setw(5) << setfill('0') << k << ends;
	
	string filename = laurie.str();
	ifstream data_file(filename.c_str());
	data_file.precision(12);

	for(int i=0; i < Nx; i++){
		for(int j=0; j < Ny; j++){
 
       data_file >> ignore >> ignore >> real_in ;
            w(i,j) = real_in;
            //cout << w(i,j) << endl;
      }
    }


    for(int number=1; number < N_avg; number++ ){
    
        cout << "number = " << number << endl;
       i_center = distribution(mt);                 // defines the center of the circle
        j_center = distribution(mt);
       // cout << i_center << " " << j_center << endl;
     


    for( int r=1; r < rtot; r++){       // iterates over the radius of ring (in terms of points)



        circ = 0.0;         
        circ2 = 0.0;

        for(int i=-r; i < r+1; i++){
            for(int j=-r; j < r+1; j++){


                i2 = i_center + i;      //determine the point of interest
                j2 = j_center + j;
                dr = pow(i*i + j*j,0.5);    // distance from centre
               
                if( dr  < r ){              //this is required to make it a circle
                   
                    if(i2 > Nx-1) i2 -= Nx;       //moves point if it goes through the periodic boundary
                    if(i2 < 0) i2 += Nx; 

                    if(j2 > Ny-1) j2 -= Ny;
                    if(j2 < 0) j2 += Ny;

             
                        circ += w(i2,j2);
                     //   circ2 += pow(w(i2,j2),2.0);

                    }
                }
            }

       //    circ /= (double) pi*pow(r,2.0);
         //   circ2 /= (double) pi*pow(r,2.0);    //normalizes by the area of the circle.
            
          //  cout << rr*dx << " " << circ << " " << circ2 << endl; 
            cir_r(r) += circ;      //adds to array
            cir2_r(r) += pow(circ,2.0);
           
        }
        

    }

        cir_r /= (double) N_avg;
        cir2_r /= (double) N_avg;



        ostringstream out_spec;
 out_spec << "./polar." << setw(5) << setfill('0') << k << ends;
string filenamespec = out_spec.str();
 ofstream fout_spec(filenamespec.c_str());
fout_spec.precision(12);

        for(int r=1; r < rtot; r++){
               fout_spec << r*dx<< " " << cir_r(r) << " " << cir2_r(r) << endl;
            }

 
}
    
    return 0;
}

