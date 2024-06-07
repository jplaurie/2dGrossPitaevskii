void computeEnergy(cx_mat , cx_mat , double & , double & , double & );
void computeWaveaction(cx_mat , double & );
void computeDissipationRate( cx_mat , double & , double & , double & , double & );
void computeSpectrum( cx_mat , rowvec & );
void computeFlux(cx_mat , cx_mat ,cx_mat,cx_mat & ,cx_mat &, rowvec &  ,rowvec & ,fftw_plan , fftw_plan );