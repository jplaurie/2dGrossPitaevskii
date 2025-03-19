import numpy as np
import matplotlib.pyplot as plt


dt = 1.e-3

L=2.0*np.pi
mode_index = 100
start=1
end=559

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

time = np.zeros(end-start+1)
mode = np.zeros((end-start+1,7))



index = 0
for file in range(start, end+1):
	
	filename = '/data2/2d_gp_selfsimilar/decay_at_256_long/output/spectrum.%.5d' % (file);
	spec = np.loadtxt(filename)

	time[index] = dt*file
	mode[index,0] = spec[4,1]
	mode[index,1] = spec[16,1]
	mode[index,2] = spec[32,1]
	mode[index,3] = spec[64,1]
	mode[index,4] = spec[128,1]
	mode[index,5] = spec[512,1]
	mode[index,6] = spec[1024,1]
	index = index+1

fig, axs = plt.subplots()

# w = 0.5*k**2
axs.loglog(time,mode[:,0] ,linewidth=5, label=r'$\omega=8$')	
axs.loglog(time,mode[:,1] ,linewidth=5, label=r'$\omega=128$')	
axs.loglog(time,mode[:,2] ,linewidth=5, label=r'$\omega=512$')	
axs.loglog(time,mode[:,3] ,linewidth=5, label=r'$\omega=2048$')	
axs.loglog(time,mode[:,4] ,linewidth=5, label=r'$\omega=8192$')	
axs.loglog(time,mode[:,5] ,linewidth=5, label=r'$\omega=131072$')	
#axs.loglog(time,mode[:,6] ,linewidth=5, label=r'$\omega=524288$')	
axs.loglog(time,0.01*time**(-2./3.) ,linewidth=5, color='black', ls='--', label=r'$\propto t^{-2/3}$')	
#axs[0].set_xlim(0,0.15)	
#axs[0].set_ylim(1e-6,1)	


axs.set_xlabel(r'$t$')
axs.set_ylabel(r'$n_\omega$')
axs.legend(loc=0,ncol=2)

plt.tight_layout()
plt.savefig('mode_evolution.pdf')


fig, axs = plt.subplots()


axs.loglog(time,mode[:,0] ,linewidth=5, label=r'$k=1$')	
axs.loglog(time,mode[:,1] ,linewidth=5, label=r'$k=5$')	
axs.loglog(time,mode[:,2] ,linewidth=5, label=r'$k=10$')	
axs.loglog(time,mode[:,3] ,linewidth=5, label=r'$k=50$')	
axs.loglog(time,mode[:,4] ,linewidth=5, label=r'$k=100$')	
axs.loglog(time,mode[:,5] ,linewidth=5, label=r'$k=500$')	
axs.loglog(time,mode[:,6] ,linewidth=5, label=r'$k=1000$')	
axs.loglog(time,50*time**(-2./3.) ,linewidth=5, color='black', ls='--', label=r'$\propto t^{-2/3}$')	
#axs[0].set_xlim(0,0.15)	
#axs[0].set_ylim(:,:)	


axs.set_xlabel(r'$t$')
axs.set_ylabel(r'$n_k$')
axs.legend(loc=0,ncol=2)

plt.tight_layout()
plt.savefig('mode_evolution_time_rescaled.pdf')
