
const int Nx = 512;                      //number of spatial grid points
const int Mx = 3*Nx/2;                      //defined the size of the anti-aliasing array M=3N/2
const int Ny = 512;                      //number of spatial grid points
const int My = 3*Ny/2;                      //defined the size of the anti-aliasing array M=3N/2

const double g = 400.0;
const double chem_pot = g;


const int num_threads = 1;



const int NBIN = fmax(Nx,Ny);
//===============set timestepping parameters============================

const bool FLAG_TIMESTEP_ETDRK = true; //if true routine is ETDRK4 else RK2
const int FLAG_ETDRK_ORDER = 4;
const double dt = 1.e-5;                //time step
const int nsteps = 200000;                    //number of total time steps
const int outstep = 100;                    //outputs data at these time steps
//======================parameters for domain ===============================
const double pi = 3.14159265358979323846;       //pi
const double Lx = 2.0*pi;                   //length of the box
const double dx = Lx/ (double) Nx;                      //grid size
const double Ly = 2.0*pi;                   //length of the box
const double dy = Ly/ (double) Ny;                      //grid size


const bool FLAG_DIPOLE = true;
const double dipole_length = pi/8.0; // original distance pi/4.0;
