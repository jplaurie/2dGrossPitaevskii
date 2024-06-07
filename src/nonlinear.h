using namespace arma;

void embed_N_M( cx_mat, cx_mat &);
void embed_M_N( cx_mat, cx_mat &);
void dealias(cx_mat &);
cx_mat nonlinear( cx_mat , double , cx_mat ,cx_mat , cx_mat & ,cx_mat & , fftw_plan, fftw_plan );