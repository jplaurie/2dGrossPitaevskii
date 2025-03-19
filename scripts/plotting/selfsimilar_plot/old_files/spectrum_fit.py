import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 1.e-1
N=250
N_ensemble=1
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




#linear low/high k fits for Rayleigh-Jeans distribution


fig, axs = plt.subplots(2, 2,figsize=(16,10))

for i in range(150,N+1,25):
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
	n = spec[kmin_low:kmax_low,1]/spec[kmin_low:kmax_low,0];
	popt1, pcov1 = curve_fit(scaling_constant_fit, w, n)
	
	w = spec[kmin_high:kmax_high,0]**2;
	n = spec[kmin_high:kmax_high,1]*spec[kmin_high:kmax_high,0];
	popt2, pcov2 = curve_fit(scaling_constant_fit, w, n)
	
	T  = popt2[0]
	mu = T/popt1[0]
	
	def scaling_flux_fit(x,a3,a4):
	    return a3/(a4+x)
	
	w = spec[kmin_low:kmax_high,0]**2;
	n = spec[kmin_low:kmax_high,1]/spec[kmin_low:kmax_high,0];
	popt3, pcov3 = curve_fit(scaling_flux_fit, w, n)
	
	T_RJ=popt3[0]
	mu_RJ=popt3[1]
	
	print('t = ' + str(round(i*dt,1)) + ' T = ' + str(T) + ' mu = ' + str(mu) + ' T_RJ = ' + str(T_RJ) + ' mu_RJ = ' + str(mu_RJ))
	
	
	axs[0][0].plot( spec[1:,0]**2.0,spec[1:,1]/spec[1:,0],linewidth=5,label=r'%.2f' % round(dt*i,1))	
	axs[0][1].plot( spec[1:,0]**2.0, (spec[1:,1]/spec[1:,0] -  T*( mu + spec[1:,0]**2.0)**(-1.0)) / (T*( mu + spec[1:,0]**2.0)**(-1.0))  ,linewidth=5,label=r'%.2f' % round(dt*i,1))	
	axs[1][0].plot( spec[1:,0]**2.0,spec[1:,1]/spec[1:,0],linewidth=5,label=r'%.2f' % round(dt*i,1))	
	axs[1][1].plot( spec[1:,0]**2.0,  (spec[1:,1]/spec[1:,0] -  T_RJ*( mu_RJ + spec[1:,0]**2.0)**(-1.0)) / (T_RJ*( mu_RJ + spec[1:,0]**2.0)**(-1.0)) ,linewidth=5,label=r'%.2f' % round(dt*i,1))	
  # Rayleigh jeans
	if(i==250):
		axs[0][0].plot( spec[kmin_low:kmax_high,0]**2.0, T*( mu + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3,label=r'$n^{\rm lim}_{RJ} = \frac{T}{\mu+\omega}$')
		#axs[1].plot( spec[kmin_low:kmax_high,0]**2.0, (T*spec[kmin_low:kmax_high,0]**2.0)/(mu+ spec[kmin_low:kmax_high,0]**2.0),color='black', linestyle=':',linewidth=3,label=r'$n_\omega = \frac{T}{\mu+\omega}$')
		axs[1][0].plot( spec[kmin_low:kmax_high,0]**2.0, T_RJ*( mu_RJ + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3,label=r'$n^{\rm full}_{RJ} = \frac{T}{\mu+\omega}$')
		
		#shade backgorund
		axs[0,0].fill_between(spec[kmin_low:kmax_low,0]**2.0, 0, 10, facecolor='gray', alpha=0.3)
		axs[0,0].fill_between(spec[kmin_high:kmax_high,0]**2.0, 0, 10, facecolor='gray', alpha=0.3)
		axs[0,1].fill_between(spec[kmin_low:kmax_low,0]**2.0, -1, 1, facecolor='gray', alpha=0.3)
		axs[0,1].fill_between(spec[kmin_high:kmax_high,0]**2.0, -1, 1, facecolor='gray', alpha=0.3)
			
	else:
		axs[0][0].plot( spec[kmin_low:kmax_high,0]**2.0, T*( mu + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3)
		#axs[1].plot( spec[kmin_low:kmax_high,0]**2.0, (T*spec[kmin_low:kmax_high,0]**2.0)/(mu+ spec[kmin_low:kmax_high,0]**2.0),color='black', linestyle=':',linewidth=3)
		axs[1][0].plot( spec[kmin_low:kmax_high,0]**2.0, T_RJ*( mu_RJ + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3)
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
axs[0][0].set_xlim(1,1000**2)	
axs[0][0].set_ylim(1.e-8,10)	
axs[0][1].set_xlim(1,1000**2)	
axs[0][1].set_ylim(-1,1)	
axs[1][0].set_xlim(1,1000**2)	
axs[1][0].set_ylim(1.e-8,10)
axs[1][1].set_xlim(1,1000**2)	
axs[1][1].set_ylim(-1,1)


	
axs[0][0].set_xlabel(r'$\omega$')
axs[0][0].set_ylabel(r'$n_\omega$')
axs[0][1].set_xlabel(r'$\omega$')
axs[0][1].set_ylabel(r'$(n_\omega - n^{\rm lim}_{RJ})/n^{\rm lim}_{RJ}$')
axs[1][0].set_xlabel(r'$\omega$')
axs[1][0].set_ylabel(r'$n_\omega$')
axs[1][1].set_xlabel(r'$\omega$')
axs[1][1].set_ylabel(r'$(n_\omega - n^{\rm full}_{RJ})/n^{\rm full}_{RJ}$')



axs[0][0].legend(fontsize=12, loc="best", ncol=3)
axs[1][0].legend(fontsize=12, loc="best", ncol=3)
axs[0][1].legend(fontsize=12, loc="best", ncol=3)
axs[1][1].legend(fontsize=12, loc="best", ncol=3)
axs[0][0].set_yscale('log')
axs[0][0].set_xscale('log')
axs[1][0].set_yscale('log')
axs[1][0].set_xscale('log')
#axs[0][1].set_yscale('log')
axs[0][1].set_xscale('log')
#axs[1][1].set_yscale('log')
axs[1][1].set_xscale('log')

plt.tight_layout()
plt.savefig('selfsimilar_RJ_fit.pdf')



#spectrum fits at low/high k


fig, axs = plt.subplots(1, 2,figsize=(15,5))


for i in range(150,N+1,25):
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
	n = spec[kmin_low:kmax_low,1]/spec[kmin_low:kmax_low,0];
	popt1, pcov1 = curve_fit(scaling_fit, w, n)
	
	w = spec[kmin_high:kmax_high,0]**2;
	n = spec[kmin_high:kmax_high,1]*spec[kmin_high:kmax_high,0];
	popt2, pcov2 = curve_fit(scaling_fit, w, n)
	
	x = -(popt2[1]-1.0)
	a= x/(2.0*(x-1.0))
	b= 1.0/(2.0*(x-1.0))
	
	
	
	print('t = ' + str(round(i*dt,1)) + ' x_low = ' + str(popt1[1]) + ' x = ' + str(popt2[1]-1) + ' a = ' + str(a) + ' b = ' + str(b))
	
	axs[0].plot( spec[1:,0]**2.0,spec[1:,1]/spec[1:,0],linewidth=5,label=r'%.2f' % round(dt*i,1))	
	axs[1].plot( spec[1:,0]**2.0,spec[1:,1]*spec[1:,0],linewidth=5,label=r'%.2f' % round(dt*i,1))	
  
  #best fit lines
	axs[0].plot( spec[kmin_low:kmax_low,0]**2.0, popt1[0]*spec[kmin_low:kmax_low,0]**(2.0*popt1[1]),color='black', linestyle=':',linewidth=3)
	axs[1].plot( spec[kmin_high:kmax_high,0]**2.0, popt2[0]*spec[kmin_high:kmax_high,0]**(2.0*popt2[1]),color='black', linestyle=':',linewidth=3)
  

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
plt.savefig('selfsimilar_best_fit.pdf')


#estimating t*

M=250
zeroMode=np.zeros((M-1,2))

counter=0
for i in range(1,M,1):
	spec_avg = 0.0*spec_avg
	#filename = './data/spectrum.%.5d' % (i+j);
	for j in range(0,N_ensemble,1):
		#filename = './data/spectrum.%.5d' % (i+j);
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec
	
	spec_avg /= N_ensemble
	spec=spec_avg
	zeroMode[counter,1] = np.sum(spec[0:4,1])/5.0
	zeroMode[counter,0] = i*dt
	counter= counter + 1

#x=1.01
#a= 0.5*x/(x-1)
def t_fit(t,a1,t1):
    return t1*t**a1


#popt4, pcov4 = curve_fit(t_fit, zeroMode[:,0], zeroMode[:,1])
#print('c = ' + str(popt4[0]) + ' t* = ' + str(popt4[1]))
fig, axs = plt.subplots(1, 1,figsize=(8,5))
axs.plot( zeroMode[:,0],zeroMode[:,1],linewidth=5)
#axs.plot( zeroMode[:,0],1.e-40*(30.0-zeroMode[:,0])**20.0,linewidth=5)
#axs.plot( zeroMode[:,0],popt4[1]*(zeroMode[:,0])**popt4[0],linewidth=5)

axs.set_xlim(0.1,100)	
#axs.set_ylim(1.e-8,10)	
axs.set_xlabel(r'$t$')
axs.set_ylabel(r'$n_\omega(0)$')
axs.set_yscale('log')
axs.set_xscale('log')
plt.tight_layout()
plt.savefig('selfsimilar_t_fit.pdf')

