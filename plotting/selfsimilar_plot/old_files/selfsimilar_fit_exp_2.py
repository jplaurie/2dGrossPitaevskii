import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 0.5*1e-3
file_start = 150 # 150
file_end = 250   # 250
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




#def scaling_fit(x,a1,b1):
 #   return a1*(x**b1)
	
#def scaling_constant_fit(x,a2):
 #   return a2



filename = './data/spectrum.00125'
spec_avg = np.loadtxt(filename)
#spec_avg = 0.0*spec_avg


def penalty_fnc(c1,c2):
	pen = (np.arange(kmax_high-1)+1)**2.0
	pen[0:kmax_low] = c1*np.ones(kmax_low )
	pen[kmax_low:kmin_high] = 1.e-12
	pen[kmin_high:] = c2*pen[kmin_high:]
	return pen

def scaling_rayleigh_jeans_fit(x,T,mu):
    return ((T/(mu+x))*penalty_fnc(c1,c2))

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

	w = spec[1:kmax_high,0]**2; # w = k^2
	n = spec[1:kmax_high,1]/(2.0*spec[1:kmax_high,0]);
	
	c1 = 1.0/n[0]
	c2 = 1.0/(n[kmin_high]*kmin_high**2.0)

	nn = n * penalty_fnc(c1,c2)

	parameters_fit, covariance_fit = curve_fit(scaling_rayleigh_jeans_fit, w, n * penalty_fnc(c1,c2))

	T  = parameters_fit[0]
	mu = parameters_fit[1]

	T_mu_data[i,0] = i*dt
	T_mu_data[i,1] = T
	T_mu_data[i,2] = mu


mu0 = T_mu_data[file_start,2]
T0 = T_mu_data[file_start,1]
time0 = T_mu_data[file_start,0]


muf = T_mu_data[file_end,2]
Tf = T_mu_data[file_end,1]
timef = T_mu_data[file_end,0]


total_energy = total_energy/(L*L*wf)

total_energy = wf*  0.6760918994815828

print("muf = " + str(muf))
print("Tf = " + str(Tf))


axs1[0].plot( T_mu_data[:,0],T_mu_data[:,1], linestyle='-',linewidth=5)
axs1[1].plot( T_mu_data[:,0],T_mu_data[:,2], linestyle='-',linewidth=5)

# wf = inf formula
#cf = muf*np.exp(Q*timef/Tf) 
#axs1[1].plot( T_mu_data[:,0],cf*np.exp(-Q*(T_mu_data[:,0])/Tf), linestyle=':',color='black',linewidth=5)#,label=r'$\mu(t) = \mu_0\exp\left( - \frac{Qt}{T_0} \right)$')



axs1[0].plot( T_mu_data[:,0], Tf * (T_mu_data[:,0]/T_mu_data[:,0]), linestyle='--',color='black',linewidth=5)

# finite correction formula
cf = (wf + mu0)/(mu0*np.exp(Q*time0/Tf))
axs1[1].plot( T_mu_data[:,0], wf / (cf*np.exp(Q*T_mu_data[:,0]/Tf)-1.0), linestyle='--',color='black',linewidth=5)#,label=r'$\mu(t) = \mu_0\exp\left( - \frac{Qt}{T_0} \right)$')


#axs1[1].plot( T_mu_data[1:,0],  wf/ (np.exp((wave_action_flux*(T_mu_data[1:,0])/Ts)*np.exp(wave_action_flux*(T_mu_data[1:,0])/(2.*Ts)))-1.0), linestyle='dotted',color='red',linewidth=5,label=r'$\mu(t)=\frac{\omega_f}{e^\frac{Qte^{Qt/2T^*}}{T^*}-1}$')
axs1[0].set_xlim(file_start*dt,file_end*dt)	
#axs[0].set_ylim(1.e-8,10)	
axs1[1].set_xlim(file_start*dt,file_end*dt)		
#axs2[0].set_xlim(file_start*dt,file_end*dt)	
#axs2[1].set_xlim(file_start*dt,file_end*dt)	

axs1[0].set_ylim(0,2)	
axs1[1].set_ylim(0,30000)	

#axs2[0].set_ylim(-.2,.2)	

#axs2[1].set_ylim(-1.2,-0.8)	


	
axs1[0].set_xlabel(r'$t$')
axs1[0].set_ylabel(r'$T(t)$')
axs1[1].set_xlabel(r'$t$')
axs1[1].set_ylabel(r'$\mu(t)$')



# inset axes....
###axins0 = axs1[0].inset_axes([0.4, 0.4, 0.57, 0.55])
###axins0.plot( T_mu_data[:,0],T_mu_data[:,1], linestyle='-',linewidth=5)
###axins0.plot( T_mu_data[:,0], Tf *T_mu_data[:,0]/T_mu_data[:,0], linestyle='--',color='black',linewidth=5)
###axins0.set_xlim(0.0625, 0.125)
###axins.set_ylim(10, 100000)
#axins.set_xticklabels('')
#axins.set_yticklabels('')
###axins0.set_xlabel(r'$t$')
###axins0.set_ylabel(r'$T(t)$')
###axins0.set_yscale('log')



# inset axes....
axins = axs1[1].inset_axes([0.4, 0.4, 0.57, 0.55])
axins.plot( T_mu_data[:,0],T_mu_data[:,2], linestyle='-',linewidth=5)
axins.plot( T_mu_data[:,0], wf / (cf*np.exp(Q*T_mu_data[:,0]/Tf)-1.0), linestyle='--',color='black',linewidth=5)
axins.set_xlim(file_start*dt+0.001, file_end*dt)
###axins.set_ylim(10, 100000)
#axins.set_xticklabels('')
#axins.set_yticklabels('')
axins.set_xlabel(r'$t$')
axins.set_ylabel(r'$\mu(t)$')
axins.set_yscale('log')

#axs2[0].set_xlabel(r'$t$')
#axs2[0].set_ylabel(r'$x_{low}$')
#axs2[1].set_xlabel(r'$t$')
#axs2[1].set_ylabel(r'$x_{high}$')




#axs1[0].legend(fontsize=12, loc="best", ncol=1)
#axs1[1].legend(fontsize=12, loc="best", ncol=1)
#axs2[0].legend(fontsize=12, loc="best", ncol=1)
#axs2[1].legend(fontsize=12, loc="best", ncol=1)

#axs[0].set_yscale('log')
#axs[1].set_xscale('log')
#axs1[1].set_yscale('log')

fig1.tight_layout()

fig1.savefig('rescaled_selfsimilar_T_mu_evolution_v2.pdf')
#fig2.savefig('rescaled_selfsimilar_scalings.pdf')














#linear low/high k fits for Rayleigh-Jeans distribution


fig, axs = plt.subplots(1, 2,figsize=(15,5))


for i in range(file_start,file_end+1,25):
	spec_avg *= 0.0
	for j in range(0,N_ensemble,1):
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spec = np.loadtxt(filename)
		spec_avg += spec
	
	spec_avg /= N_ensemble
	spec=spec_avg

	
	w = spec[1:kmax_high,0]**2; # w = k^2
	n = spec[1:kmax_high,1]/(2.0*spec[1:kmax_high,0]);

	c1 = 1./n[0]
	c2 = 1.0/(n[kmin_high]*kmin_high**2.0)

	parameters_fit, covariance_fit = curve_fit(scaling_rayleigh_jeans_fit, w, (n*penalty_fnc(c1,c2)))


	T  = parameters_fit[0]
	mu = parameters_fit[1]
	

	x00 = spec[1:,0]**2.0 # w
	y00 = spec[1:,1]/(2.0*spec[1:,0]) # n
	y01 = (spec[1:,1]/(2.0*spec[1:,0]) -  T*( mu + spec[1:,0]**2.0)**(-1.0)) / (T*( mu + spec[1:,0]**2.0)**(-1.0))
	y02 = (spec[1:,1]/(2.0*spec[1:,0]) -  T*( mu + spec[1:,0]**2.0)**(-1.0)) 

	if(i==file_start):
		axs[1].plot( spec[:,0]**2.0, 0.0*spec[:,0]**2.0,color='black', linestyle='-',linewidth=1)


	axs[0].plot( x00,y00,linewidth=5,label=r'%.4f' % round(dt*i,4))	# n_w
	axs[1].plot( x00, y01,linewidth=5,label=r'%.4f' % round(dt*i,4))	# relative error



  
  # Rayleigh jeans
	if(i==250):
		axs[0].plot( spec[kmin_low:kmax_high,0]**2.0, T*( mu + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3,label=r'$n^{RJ}_\omega = \frac{T(t)}{\omega+\mu(t)}$')
			
	else:
		axs[0].plot( spec[kmin_low:kmax_high,0]**2.0, T*( mu + spec[kmin_low:kmax_high,0]**2.0)**(-1.0),color='black', linestyle=':',linewidth=3)
		
axs[0].set_xlim(1,1000**2)	
axs[0].set_ylim(1.e-7,10)	
axs[1].set_xlim(1,1000**2)	



axs[1].set_ylim(-1,1)	

axs[0].set_xlabel(r'$\omega$')
axs[0].set_ylabel(r'$n_\omega$')
axs[1].set_xlabel(r'$\omega$')
axs[1].set_ylabel(r'$(n_\omega - n^{RJ}_\omega)/n^{RJ}_{\omega}$')


axs[0].legend(fontsize=12, loc="best", ncol=4)
axs[1].legend(fontsize=12, loc="best", ncol=3)
axs[0].set_yscale('log')
axs[0].set_xscale('log')

axs[1].set_xscale('log')


plt.tight_layout()
plt.savefig('rescaled_selfsimilar_RJ_fit_v2.pdf')






fig, axs = plt.subplots(1, 1,figsize=(8,5))

a = Q/ Tf
print("a = Q/Tf = " + str(a))
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
	
	axs.plot(eta, f,linewidth=5,label=r'$t =$ %.4f' % round(dt*i,4))	


axs.set_xlabel(r'$\eta$')
axs.set_ylabel(r'$f^{RJ}(\eta)$')

axs.set_yscale('log')
axs.set_xscale('log')

axs.legend(fontsize=12, loc="best", ncol=1)
axs.set_ylim(1.e-16,1.e-6)	
plt.tight_layout()
plt.savefig('rescaled_selfsimilar_exp_scaling_v2.pdf')
