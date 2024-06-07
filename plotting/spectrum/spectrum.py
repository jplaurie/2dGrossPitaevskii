import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)


N=300
L=2.0*np.pi
dt=0.5*0.001


for i in range(50,N+1,50):

	filename = '/data2/2d_gp_selfsimilar/ensemble_1/output/spectrum.%.5d' % i;
	print(filename)
	spec = np.loadtxt(filename)
	l1, = plt.plot(0.5*spec[:,0]**2, spec[:,1]/spec[:,0], linewidth=3,label=r't=%.3f' % round(dt*i,4))

l8, =plt.plot(0.5*spec[2:300,0]**2, 1.e-6*spec[2:300,0]/spec[2:300,0],linewidth=3,label=r'$\propto \omega^{0}$', color='black', linestyle='-.')
l9, =plt.plot(0.5*spec[10:450,0]**2, 5.e0*spec[10:450,0]**-2,linewidth=3,label=r'$\propto \omega^{-1}$', color='black', linestyle=':')

plt.xlabel(r'$\omega$')
plt.ylabel(r'Wave Action Spectrum $n_\omega$')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1,(1024*1.5)**2)
plt.ylim(1.e-12,100)
plt.legend(fontsize=12, loc='best', ncol=2)
plt.tight_layout()
plt.savefig('spectrum.pdf')








