
const int Nx = 512;                      //number of spatial grid points
const int Mx = 3*Nx/2;                      //defined the size of the anti-aliasing array M=3N/2
const int Ny = 512;                      //number of spatial grid points
const int My = 3*Ny/2;                      //defined the size of the anti-aliasing array M=3N/2

const double g = 400.0; //self-similar g = 2000
const double chem_pot = g;


const int num_threads = 1;


const int NBIN = fmax(Nx,Ny);
//===============set timestepping parameters============================

const bool FLAG_TIMESTEP_ETDRK = true; //if true routine is ETDRK4 else RK2
const int FLAG_ETDRK_ORDER = 4;
const double dt = 1.e-5;                //time step
const int nsteps = 99999999;                    //number of total time steps
const int outstep = 100000;                    //outputs data at these time steps
//======================parameters for domain ===============================
const double pi = 3.14159265358979323846;       //pi
const double Lx = 2.0*pi;                   //length of the box
const double dx = Lx/ (double) Nx;                      //grid size
const double Ly = 2.0*pi;                   //length of the box
const double dy = Ly/ (double) Ny;                      //grid size


//===========parameters from equation =============================
const bool FLAG_GINZBURG_LANDAU_DISS = false;
const double gamma_0 = 0.e-3;
const double nu = 0.0;//e-3;					//hyperviscosity coefficient
const double nupower = 1.0;					//power of the laplacian 	
const double alpha = 0.0;					//friction coefficient
const double alphapower = 0.0;				//negative power of the laplacian


//===================parameters for forcing routine ===============================
//const bool FLAG_FORCE_ON = true; //Flag to turn on forcing
const int FLAG_FORCE_TYPE = 0;	// 0 = no forcing, 1 = additive, 2 = fixed amplitude, 3= fixed from inital condition
const bool FLAG_FORCE_AMP_RESCALE = false; // Flag for having forcing amplitude rescaled to give unit energy density
const bool FLAG_FORCE_EXP = false;  // turn to true to have exponential forcing spectrum **default is annulus***
const double kf = 64;							//forcing wavenumber
const double dk = 2;//10.0;							//width of forcing annulus (in terms of number of modes)
const double Amp = 5.e-3;						//amplitude of forcing
const double force_power = 4.0; // for exp force only f_k = (k/k_f)^force_power * exp( -(k / k_f)^force_power )
const double mean_waveaction = 10.0;

//===================parameters for filter routine ===============================
const bool FLAG_SOUND_FILTER = false;					// FLAG_filter = 0 then no noise filter is applied
const int filterstep = 0;
const double sigma= 3.0;					// variance of Gaussian filter
const double eps_fil = 0.5;						// define vortex below this intensity (real and imag parts separately)
const double vortex_ring_dist = dx* 3.0;	//cookie cutter distance
