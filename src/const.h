
/*
Programme to solve the 2D Schrödinger / Gross-Pitaevskii equation

(i-Gamma) PSI_t =  c * (PSI_xx + PSI_yy)    +  g * PSI |PSI|^2  +   mu * PSI    + dissipation + forcing

Author: Jason Laurie
Date: 27/10/2023

*/


//define pi
const double pi = 3.14159265358979323846;

//number of grid points
const int Nx = 512;
const int Ny = 512;

//length of periodic box
const double Lx = 2.0 * pi;
const double Ly = 2.0 * pi;

/* equation parameters */

//linear parameter
const double c = -1.0;


//nonlinear parameter
const double g = 1000.0;

//chemical potential
const double mu = 0.0;


//grid spacing
const double dx = Lx / double(Nx);
const double dy = Ly / double(Ny);

//define the size of the anti-aliasing array for 3/2 - rule 
const int Mx = 3*Nx/2;
const int My = 3*Ny/2;

//openmp threads (not in use)
const int num_threads = 1;


const int NBIN = fmax(Nx,Ny);



//=============== timestepping ============================

/*
set timestepping method
4th order exponential time differencing Runge-Kutta -B = "ETDRK4"
2nd order Runge-Kutta with integrating factors = "RK2"     */
const std::string FLAG_TIMESTEP_METHOD = "ETDRK4";

//time step
const double dt = 1.e-5;

//total time steps
const int total_steps = 999999999;       //number of total time steps

//time at which you would like state and diagnostics to be recorded
const double output_time = 0.001;

//set this flag to be true to record diagnostics (code will run faster if this == "false")
const bool FLAG_OUTPUT_DIAGONOSTICS = "true";


//=========== dissipation =============================
//turn hyperviscosity on/off
const bool FLAG_HYPER_DISSIPATION = true;

//hyperviscosity coefficient
const double nu = 1.e-36;//1.e-32;      

//positive power of the laplacian       
const double nupower = 8.0;

//turns on a hyperviscosity cutoff
const bool FLAG_HYPER_CUTOFF = false;   

//value of cutoff the hyperviscosity (wavenumber)
const double nucutoff = 1024.0;             



//turn hypoviscosity on/off
const bool FLAG_HYPO_DISSIPATION = true;

//friction coefficient
const double alpha = 5.e3;

//negative power of the laplacian (alphapower < 0 or friction = 0)                              
const int alphapower = -2.0;

//turns on a hypoviscosity cutoff
const bool FLAG_HYPO_CUTOFF = false;    

//value of cutoff the hypoviscosity (wavenumber)
const double alphacutoff = 4.0;          

//adds ginzburg landau dissipation term
const bool FLAG_GINZBURG_LANDAU_DISS = false;

//strength of the ginzburg landau dissipation
const double gamma_0 = 0.e-3;


//=================== forcing  ===============================
/*
 "off" = no forcing
"annulus" = additive forcing in annulus in fourier space:      kf-dk <= kf <= kf+dk with constant force_amplitude
"gaussian" = additive forcing using Gaussian distribution:    force_amplitude * exp (-0.5 * [(k-kf)/sigmaf)]^2 ) 
"exponential" = additive forcing using Gaussian distribution:    force_amplitude * pow(k/kf,force_power)*exp(-pow(k/kf,force_power))
"exp-log = additive forcing using exponetial of log^2":           force_amplitude *  exp(-0.5*pow(log(k/kf)/sigmaf,2,0))
*/
const std::string FLAG_FORCING_TYPE = "exp-log"; //Flag to turn on additive forcing

//FLAG for having forcing amplitude rescaled to give unit energy density
const bool FLAG_FORCING_RESCALE = false; // Flag for having forcing amplitude rescaled to give unit energy density

//forcing amplitude 
const double force_amplitude = 1.e-1;//0.001;   
//forcing wavenumber
const double kf = 32;
//width of forcing annulus (in terms of wavenumber)                                     
const double dk = 2.0;
//forcing standard devation for "gaussian/exp-log" profile
const double sigmaf = 0.05;
// forcing power for "exponential" forcing
const double force_power = 4.0;

//const double mean_waveaction = 1.0; //paramter for rescaling forcing wih large-scale dissipation




//===================parameters for filter routine ===============================
const bool FLAG_SOUND_FILTER = false;					// FLAG_filter = 0 then no noise filter is applied
const int filterstep = 0;
const double sigma= 3.0;					// variance of Gaussian filter
const double eps_fil = 0.5;						// define vortex below this intensity (real and imag parts separately)
const double vortex_ring_dist = dx* 3.0;	//cookie cutter distance

