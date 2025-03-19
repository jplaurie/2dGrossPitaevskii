import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 0.5*1e-3
file_start = 150 # 150
file_end = 200   # 250
wf=512.**2.0

total_wave_action = 168.9616
total_energy = 1.50605e7
Q = 64.24 # 32.12 /2.0 # wave action flux   16.06 #32.12/2.0
energy_flux = 2.0 * 8.41969e6


N_ensemble=10
L=2.0*np.pi

kmin_low = 0
kmax_low = 10

kmin_high = 100
kmax_high = 400


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)





filename = './data/spectrum.00125'
spec_avg = np.loadtxt(filename)
#spec_avg = 0.0*spec_avg





fig, axs = plt.subplots(1, 1,figsize=(8,5))

# f = t^a N()
# eta = omega* t^b
a = -2./3.  
b= -1./3.


for i in range(file_start,file_end+1,25):
	spec_avg *= 0.0
	for j in range(0,N_ensemble,1):
		filename = './data/spectrum.%.5d' % (i+j);
		#filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec

	spec_avg /= N_ensemble
	spec=spec_avg

	tau_a = power(i*dt,a)
	tau_b = power(i*dt,b)
	eta =  spec[1:,0]**2 * tau_b ; # w = k^2
	f = tau_a * spec[1:,1]/(2.0*spec[1:,0]);
	
	axs.plot(eta, f,linewidth=5,label=r'$t =$ %.4f' % round(dt*i,4))	


axs.set_xlabel(r'$\eta$')
axs.set_ylabel(r'$f^{RJ}(\eta)$')

axs.set_yscale('log')
axs.set_xscale('log')

axs.legend(fontsize=12, loc="best", ncol=1)
axs.set_ylim(1.e-16,1.e-6)	
plt.tight_layout()
plt.savefig('selfsimilar_powerlaw_scaling.pdf')
