import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 2.*1e-3

file_start = 150
file_end = 250


#wf=512**2.0

#total_wave_action = 168.9616
#total_energy = 1.50605e7
#wave_action_flux = 16.06 #32.12/2.0
#energy_flux = 2.0 * 8.41969e6


#N_ensemble=10
L=2.0*np.pi


#fit regions for T/mu
kmin_low=1
kmax_low = 6

#fit regions for T/\omega
kmin_high=200
kmax_high=400


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)

def scaling_fit(x,a1,b1):
    return a1*(x**b1)
	
def scaling_constant_fit(x,a2):
    return a2



filename = './data/spectrum.00125'
spec_avg = np.loadtxt(filename)
spec_avg = 0.0*spec_avg



# T and mu vs times

fig1, axs1 = plt.subplots(1, 2,figsize=(16,5))
fig2, axs2 = plt.subplots(1, 2,figsize=(16,5))

T_mu_data=np.zeros((file_end+1,5))

for i in range(file_start,file_end+1,1):
	spec_avg *= 0.0
	for j in range(0,N_ensemble,1):
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec
	
	spec_avg /= N_ensemble
	spec=spec_avg

	w = spec[kmin_low:kmax_low,0]**2; # w = k^2
	n = spec[kmin_low:kmax_low,1]/(2.0*spec[kmin_low:kmax_low,0]);
	popt1, pcov1 = curve_fit(scaling_constant_fit, w, n)
	popt11, pcov11 = curve_fit(scaling_fit, w, n)

	
	w = spec[kmin_high:kmax_high,0]**2;
	n = spec[kmin_high:kmax_high,1]*w/(2.0*spec[kmin_high:kmax_high,0]);
	popt2, pcov2 = curve_fit(scaling_constant_fit, w, n)
	popt21, pcov21 = curve_fit(scaling_fit, w, n)
	
	T  = popt2[0]
	mu = T/popt1[0]
	scaling_low=popt11[1]
	scaling_high = popt21[1]

	T_mu_data[i,0] = i*dt
	T_mu_data[i,1] = T
	T_mu_data[i,2] = mu
	T_mu_data[i,3] = scaling_low
	T_mu_data[i,4] = scaling_high-1.0

mu_0 = T_mu_data[file_start,2]
T_0 = T_mu_data[file_start,1]
t0 = T_mu_data[file_start,0]
total_energy = total_energy/(L*L*wf)

total_energy = wf*  0.6760918994815828

print(mu_0)
print(T_0)

#mu_0 = mu_0*np.exp(wave_action_flux*t0/T_0)

def scaling_Ts(x,a4):
    return wf/ (np.exp((wave_action_flux*(x)/a4)*np.exp(wave_action_flux*(x)/(2.*a4)))-1.0) 

popTT, pcovTT = curve_fit(scaling_Ts, T_mu_data[file_start:file_end,0], T_mu_data[file_start:file_end,2])

print('Ts = ' + str(popTT[0]))

axs1[0].plot( T_mu_data[:,0],T_mu_data[:,1], linestyle='-',linewidth=5)#,label=r'Numerics')	
axs1[1].plot( T_mu_data[:,0],T_mu_data[:,2], linestyle='-',linewidth=5)#,label=r'Numerics')	
axs2[0].plot( T_mu_data[:,0],T_mu_data[:,3], linestyle='-',linewidth=5)#,label=r'$x_{low}$')	
axs2[1].plot( T_mu_data[:,0],T_mu_data[:,4], linestyle='-',linewidth=5)#,label=r'$x_{high}$')

#axs1[0].plot( T_mu_data[:,0], (total_energy /wf ) * (1.0 - (mu_0/wf)*np.exp(-wave_action_flux* (T_mu_data[:,0]-t0) / T_0 ) *(np.log(mu_0/wf)-(wave_action_flux*  (T_mu_data[:,0]-t0) / T_0) )), linestyle='--',color='black',linewidth=5,label=r'$T(t)=\frac{E}{\omega_f}\left[1 - \frac{\mu_0e^{-\frac{Qt}{T_0}}}{\omega_f}\left(\ln\left(\frac{\mu_0}{\omega_f}\right) - \frac{Qt}{T_0}\right) \right]$')

axs1[1].plot( T_mu_data[:,0],mu_0*np.exp(-wave_action_flux*(T_mu_data[:,0]-t0)/(T_0)), linestyle='--',color='black',linewidth=5)#,label=r'$\mu(t) = \mu_0\exp\left( - \frac{Qt}{T_0} \right)$')
Ts=6.75
#axs1[1].plot( T_mu_data[1:,0],  wf/ (np.exp((wave_action_flux*(T_mu_data[1:,0])/Ts)*np.exp(wave_action_flux*(T_mu_data[1:,0])/(2.*Ts)))-1.0), linestyle='dotted',color='red',linewidth=5,label=r'$\mu(t)=\frac{\omega_f}{e^\frac{Qte^{Qt/2T^*}}{T^*}-1}$')
axs1[0].set_xlim(file_start*dt,file_end*dt)	
#axs[0].set_ylim(1.e-8,10)	
axs1[1].set_xlim(file_start*dt,file_end*dt)		
axs2[0].set_xlim(file_start*dt,file_end*dt)	
axs2[1].set_xlim(file_start*dt,file_end*dt)	

axs1[0].set_ylim(0,2)	
#axs1[1].set_ylim(0,20000)	

axs2[0].set_ylim(-.2,.2)	

axs2[1].set_ylim(-1.2,-0.8)	


	
axs1[0].set_xlabel(r'$t$')
axs1[0].set_ylabel(r'$T(t)$')
axs1[1].set_xlabel(r'$t$')
axs1[1].set_ylabel(r'$\mu(t)$')

axs2[0].set_xlabel(r'$t$')
axs2[0].set_ylabel(r'$x_{low}$')
axs2[1].set_xlabel(r'$t$')
axs2[1].set_ylabel(r'$x_{high}$')




#axs1[0].legend(fontsize=12, loc="best", ncol=1)
#axs1[1].legend(fontsize=12, loc="best", ncol=1)
#axs2[0].legend(fontsize=12, loc="best", ncol=1)
#axs2[1].legend(fontsize=12, loc="best", ncol=1)

#axs[0].set_yscale('log')
#axs[1].set_xscale('log')
axs1[1].set_yscale('log')

fig1.tight_layout()

fig1.savefig('rescaled_selfsimilar_T_mu_evolution.pdf')
fig2.savefig('rescaled_selfsimilar_scalings.pdf')












