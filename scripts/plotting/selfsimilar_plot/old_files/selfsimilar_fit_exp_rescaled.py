import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 2.*1e-3
file_start = 150
file_end = 250
wf=512**2.0

total_wave_action = 168.9616
total_energy = 1.50605e7
wave_action_flux = 16.06 #32.12/2.0
energy_flux = 2.0 * 8.41969e6


N_ensemble=10
L=2.0*np.pi

kmin_low=1
kmax_low = 6

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
		#filename = './data/spectrum.%.5d' % (i+j);
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



















#linear low/high k fits for Rayleigh-Jeans distribution


fig, axs = plt.subplots(1, 2,figsize=(15,5))


for i in range(file_start,file_end+1,25):
	spec_avg *= 0.0
	#for j in range(0,N_avg,1):
	for j in range(0,N_ensemble,1):
		#filename = './data/spectrum.%.5d' % (i+j);
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec
	
	spec_avg /= N_ensemble
	spec=spec_avg

	w = spec[kmin_low:kmax_low,0]**2;
	n = spec[kmin_low:kmax_low,1]/(2.0*spec[kmin_low:kmax_low,0]);
	popt1, pcov1 = curve_fit(scaling_constant_fit, w, n)
	
	w = spec[kmin_high:kmax_high,0]**2;
	n = spec[kmin_high:kmax_high,1]*w /(2.0*spec[kmin_high:kmax_high,0]);
	popt2, pcov2 = curve_fit(scaling_constant_fit, w, n)
	
	T  = popt2[0]
	mu = T/popt1[0]
	
	def scaling_flux_fit(x,a3,a4):
	    return a3/(a4+x)
	
	w = spec[kmin_low:kmax_high,0]**2;
	n = spec[kmin_low:kmax_high,1]/(2.0*spec[kmin_low:kmax_high,0]);
	popt3, pcov3 = curve_fit(scaling_flux_fit, w, n)
	
	T_RJ=popt3[0]
	mu_RJ=popt3[1]
	
	print('t = ' + str(round(i*dt,3)) + ' T = ' + str(T) + ' mu = ' + str(mu) + ' T_RJ = ' + str(T_RJ) + ' mu_RJ = ' + str(mu_RJ))
	
	x00 = spec[1:,0]**2.0
	y00 = spec[1:,1]/(2.0*spec[1:,0])  # n_\omega
	y01 = (spec[1:,1]/(2.0*spec[1:,0]) -  T*( mu + spec[1:,0]**2.0)**(-1.0)) / (T*( mu + spec[1:,0]**2.0)**(-1.0))
	y02 = (spec[1:,1]/(2.0*spec[1:,0]) -  T*( mu + spec[1:,0]**2.0)**(-1.0)) 
	y10 = spec[1:,1]/(2.0*spec[1:,0])
	y11 = (spec[1:,1]/(2.0*spec[1:,0]) -  T_RJ*( mu_RJ + spec[1:,0]**2.0)**(-1.0)) / (T_RJ*( mu_RJ + spec[1:,0]**2.0)**(-1.0)) 
	
	y12 = max(abs(y01[:kmax_high])) # phi_max
	x12 = dt*i # time

	if(i==file_start):
		axs[1].plot( spec[:,0]**2.0, 0.0*spec[:,0]**2.0,color='black', linestyle='-',linewidth=1)


	axs[0].plot( x00,y00,linewidth=5,label=r'%.2f' % round(dt*i,3))	
	axs[1].plot( x00, y01,linewidth=5,label=r'%.2f' % round(dt*i,3))	

	#axs[1][0].plot( x00, y02,linewidth=5,label=r'%.2f' % round(dt*i,3))	
	#axs[1][1].scatter(x12,y12,s=20,marker='o',c='black')



	#axs[1][0].plot( x00,y10,linewidth=5,label=r'%.3f' % round(dt*i,3))	
	#axs[1][1].plot( x00,y11,linewidth=5,label=r'%.3f' % round(dt*i,3))	
  
  # Rayleigh jeans
	if(i==250):
		axs[0].plot( spec[kmin_low:kmax_high,0]**2.0, T*( mu + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3,label=r'$n^{RJ}_\omega = \frac{T(t)}{\omega+\mu(t)}$')
	####	axs[0][1].plot( spec[kmin_low:kmax_high,0]**2.0, 0.0*spec[kmin_low:kmax_high,0]**2.0,color='black', linestyle='-',linewidth=1)
		#axs[1].plot( spec[kmin_low:kmax_high,0]**2.0, (T*spec[kmin_low:kmax_high,0]**2.0)/(mu+ spec[kmin_low:kmax_high,0]**2.0),color='black', linestyle=':',linewidth=3,label=r'$n_\omega = \frac{T}{\mu+\omega}$')
#		axs[1][0].plot( spec[kmin_low:kmax_high,0]**2.0, T_RJ*( mu_RJ + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3,label=r'$n^{\rm full}_{RJ} = \frac{T}{\mu+\omega}$')
		
		#shade background
		axs[0].fill_between(spec[kmin_low:kmax_low,0]**2.0, 0, 10, facecolor='gray', alpha=0.3)
		axs[0].fill_between(spec[kmin_high:kmax_high,0]**2.0, 0, 10, facecolor='gray', alpha=0.3)
		axs[1].fill_between(spec[kmin_low:kmax_low,0]**2.0, -1, 1, facecolor='gray', alpha=0.3)
		axs[1].fill_between(spec[kmin_high:kmax_high,0]**2.0, -1, 1, facecolor='gray', alpha=0.3)
			
	else:
		axs[0].plot( spec[kmin_low:kmax_high,0]**2.0, T*( mu + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3)
		#axs[1].plot( spec[kmin_low:kmax_high,0]**2.0, (T*spec[kmin_low:kmax_high,0]**2.0)/(mu+ spec[kmin_low:kmax_high,0]**2.0),color='black', linestyle=':',linewidth=3)
	#	axs[1][0].plot( spec[kmin_low:kmax_high,0]**2.0, T_RJ*( mu_RJ + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3)
'''		
# guide lines
axs[0].plot( spec[2:60,0]**2.0,5.e-6*spec[2:60,0]/spec[2:60,0],color='grey', linestyle='-',linewidth=5,label=r'$n_\omega \propto \omega^{0}$')	
axs[1].plot( spec[30:300,0]**2.0,4.*spec[30:300,0]/spec[30:300,0],color='grey', linestyle='-',linewidth=5,label=r'$n_\omega \propto \omega^{-1}$')

# Rayleigh jeans
axs[0].plot( spec[kmin_low:kmax_high,0]**2.0, popt3[0]*(popt3[1] + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=5,label=r'$n_\omega^{RJ},\ t=175$')
axs[1].plot( spec[kmin_low:kmax_high,0]**2.0, (popt3[0]*spec[kmin_low:kmax_high,0]**2.0)/(popt3[1]+ spec[kmin_low:kmax_high,0]**2.0),color='black', linestyle=':',linewidth=5,label=r'$n_\omega^{RJ},\ t=175$')

# linear best fits
axs[0].plot( spec[kmin_low:kmax_low,0]**2.0,popt1[0]*(spec[kmin_low:kmax_low,0]/spec[kmin_low:kmax_low,0])**popt1[1],color='black', linestyle='--',linewidth=5,label=r'$n_\omega \propto \omega^{-0.032379}$')		
axs[1].plot( spec[kmin_high:kmax_high,0]**2.0,popt2[0]*(spec[kmin_high:kmax_high,0]/spec[kmin_high:kmax_high,0])**popt2[1],color='black', linestyle='--',linewidth=5,label=r'$n_\omega \propto \omega^{-0.015361-1}$')	
'''	
axs[0].set_xlim(1,1000**2)	
axs[0].set_ylim(1.e-8,10)	
axs[1].set_xlim(1,1000**2)	

#axs[0][2].set_xlim(1,1000**2)	

axs[1].set_ylim(-1,1)	
#axs[1][0].set_xlim(1,1000**2)	
#axs[1][0].set_ylim(1.e-8,10)
#axs[1][1].set_xlim(1,1000**2)	
#axs[1][1].set_ylim(-1,1)


#axs[1][2].set_ylim(1e-6,0.1)	

	
axs[0].set_xlabel(r'$\omega$')
axs[0].set_ylabel(r'$n_\omega$')
axs[1].set_xlabel(r'$\omega$')
axs[1].set_ylabel(r'$(n_\omega - n^{RJ}_\omega)/n^{RJ}_{\omega}$')
#axs[1][0].set_xlabel(r'$\omega$')
#axs[1][0].set_ylabel(r'$n_\omega - n^{\rm lim}_{RJ}$')

#axs[1][1].set_xlabel(r'$t$')
#axs[1][1].set_ylabel(r'$\max\{|\phi_0|\}$')

#axs[1][1].set_xlabel(r'$\omega$')
#axs[1][1].set_ylabel(r'$(n_\omega - n^{\rm full}_{RJ})/n^{\rm full}_{RJ}$')



axs[0].legend(fontsize=12, loc="best", ncol=3)
#axs[1][0].legend(fontsize=12, loc="best", ncol=3)
axs[1].legend(fontsize=12, loc="best", ncol=3)
#axs[1][1].legend(fontsize=12, loc="best", ncol=3)
axs[0].set_yscale('log')
axs[0].set_xscale('log')
#axs[1][0].set_yscale('log')
#axs[1][0].set_xscale('log')
#axs[0][1].set_yscale('log')
axs[1].set_xscale('log')
#axs[1][1].set_yscale('log')
#axs[1][1].set_yscale('log')


#axs[0][2].set_xscale('log')
#axs[1][2].set_yscale('log')

plt.tight_layout()
plt.savefig('rescaled_selfsimilar_RJ_fit.pdf')



#spectrum fits at low/high k


fig, axs = plt.subplots(1, 2,figsize=(15,6))


for i in range(file_start,file_end+1,25):
	spec_avg *= 0.0
	#for j in range(0,N_avg,1):
	#	filename = './data/spectrum.%.5d' % (i+j);
	for j in range(0,N_ensemble,1):
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec
	
	spec_avg /= N_ensemble
	spec=spec_avg

	w = spec[kmin_low:kmax_low,0]**2;
	n = spec[kmin_low:kmax_low,1]/(2.0*spec[kmin_low:kmax_low,0])
	popt1, pcov1 = curve_fit(scaling_fit, w, n)
	
	w = spec[kmin_high:kmax_high,0]**2;
	n = spec[kmin_high:kmax_high,1]*w /(2.0*spec[kmin_high:kmax_high,0])
	popt2, pcov2 = curve_fit(scaling_fit, w, n)
	
	x = -(popt2[1]-1.0)
	a= x/(2.0*(x-1.0))
	b= 1.0/(2.0*(x-1.0))
	
	
	
	print('t = ' + str(round(i*dt,3)) + ' x_low = ' + str(popt1[1]) + ' x = ' + str(popt2[1]-1) + ' a = ' + str(a) + ' b = ' + str(b))
	
	axs[0].plot( spec[1:,0]**2.0,spec[1:,1]/(2.0*spec[1:,0]),linewidth=5,label=r'%.2f' % round(dt*i,3))	
	axs[1].plot( spec[1:,0]**2.0,spec[1:,1]*spec[1:,0]/(2.0),linewidth=5,label=r'%.2f' % round(dt*i,3))	
  
  #best fit lines
	axs[0].plot( spec[kmin_low:kmax_low,0]**2.0, popt1[0]*(spec[kmin_low:kmax_low,0]**(2.0))**popt1[1],color='black', linestyle=':',linewidth=3)
	axs[1].plot( spec[kmin_high:kmax_high,0]**2.0, popt2[0]*(spec[kmin_high:kmax_high,0]**(2.0))**popt2[1],color='black', linestyle=':',linewidth=3)
  

axs[0].set_xlim(1,500**2)	
axs[0].set_ylim(1.e-8,10)	
axs[1].set_xlim(10**2,500**2)	
axs[1].set_ylim(1.e-4,100)	

axs[0].set_xlabel(r'$\omega$')
axs[0].set_ylabel(r'$n_\omega$')
axs[1].set_xlabel(r'$\omega$')
axs[1].set_ylabel(r'$\omega n_\omega$')

axs[0].legend(fontsize=12, loc="best", ncol=3)
axs[1].legend(fontsize=12, loc="best", ncol=3)
axs[0].set_yscale('log')
axs[0].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_xscale('log')

plt.tight_layout()
plt.savefig('rescaled_selfsimilar_best_fit.pdf')




fig, axs = plt.subplots(1, 1,figsize=(8,5))

a= wave_action_flux / (T_mu_data[file_start,1])

for i in range(file_start,file_end+1,25):
	spec_avg *= 0.0
	for j in range(0,N_ensemble,1):
		#filename = './data/spectrum.%.5d' % (i+j);
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec

	spec_avg /= N_ensemble
	spec=spec_avg

	scale = np.exp(a*i*dt)
	
	eta = scale * spec[1:,0]**2 ; # w = k^2
	f = (1./scale) * spec[1:,1]/(2.0*spec[1:,0]);
	
	axs.plot(eta, f,linewidth=5,label=r'$t =$ %.2f' % round(dt*i,3))	



	
#axs.set_xlabel(r'$\eta = \omega\exp(at),\ a = $%.3f' % a)
#axs.set_ylabel(r'$f(\eta) =\exp(-at) n_\omega,\ a = $%.3f' % a)

axs.set_xlabel(r'$\eta$')
axs.set_ylabel(r'$f^{RJ}(\eta)$')

axs.set_yscale('log')
axs.set_xscale('log')

axs.legend(fontsize=12, loc="best", ncol=1)
axs.set_ylim(1.e-12,1.e-4)	
plt.tight_layout()
plt.savefig('rescaled_selfsimilar_exp_scaling.pdf')
