import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 0.5*1.e-3
N=250
N_ensemble=1
L=2.0*np.pi



plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)




fig, axs = plt.subplots(1, 3,figsize=(25,5))



		#filename = './data/spectrum.%.5d' % (i+j);
filename = '/data2/2d_gp_selfsimilar/ensemble_1/output/energy.tot';
energy = np.loadtxt(filename)
filename = '/data2/2d_gp_selfsimilar/ensemble_1/output/waveaction.tot';
waveaction = np.loadtxt(filename)


	
axs[0].plot( 0.5*energy[:,0],energy[:,4]/L**2,linewidth=5)	
axs[1].plot( 0.5*waveaction[:,0],waveaction[:,1]/L**2,linewidth=5,label=r'total wave action, $N$')	
axs[2].plot( 0.5*waveaction[:,0],waveaction[:,1]/(L**2*0.5*waveaction[:,0]),linewidth=5,label=r'$Q=N/t$')	

axs[0].set_xlim(0,0.15)	
#axs[0].set_ylim(:,:)	
axs[1].set_xlim(0,0.15)	
axs[1].set_xlim(0,0.15)	
#axs[1].set_ylim(:,)	

axs[0].set_ylabel(r'total energy, $E$')

axs[0].set_xlabel(r'time')
axs[1].set_ylabel(r'total wave action, $N$')

axs[1].set_xlabel(r'time')

axs[2].set_xlabel(r'time')
axs[2].set_ylabel(r'$Q=N/t$')

plt.tight_layout()
plt.savefig('total_quantities.pdf')

fig, axs = plt.subplots(1, 2,figsize=(16,5))


print(energy[150,4])
print(waveaction[150,1])

filename = '/data2/2d_gp_selfsimilar/ensemble_1/output/dissipation.tot';
dissipation = np.loadtxt(filename)
		
axs[0].plot( 0.5*dissipation[:,0],2.0*dissipation[:,4],linewidth=5)	
axs[1].plot( 0.5*dissipation[:,0],2.dissipation[:,2],linewidth=5)	

axs[0].set_xlim(0,0.15)	
#axs[0].set_ylim(:,:)	
axs[1].set_xlim(0,0.15)	
#axs[1].set_ylim(:,)	



axs[0].set_xlabel(r'time')
axs[0].set_ylabel(r'energy dissipation rate')

axs[1].set_xlabel(r'time')
axs[1].set_ylabel(r'wave action dissipation rate')


plt.tight_layout()
plt.savefig('dissipation.pdf')

