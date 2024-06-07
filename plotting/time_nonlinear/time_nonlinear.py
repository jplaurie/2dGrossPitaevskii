import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 1.e-1
N=250
g=2000.0
Lx = 2.0*np.pi
Ly = 2.0*np.pi
Nx = 2048
Ny = 2048

N_ensemble=3
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)


fig, axs = plt.subplots()


eps = np.zeros( np.int_(2*Nx))
eps2 = np.zeros( np.int_(2*Nx))
nn = np.zeros( np.int_(2*Nx))
nn_dot = np.zeros( np.int_(2*Nx))
normalise = np.zeros( np.int_(2*Nx))


for k in range(1,N_ensemble,1):
	print('ensemble = ' + str(k))	
	filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/data/psi.00200' % (k);
	data = np.loadtxt(filename)
		

	
	psi = np.reshape(data[:,2], (Nx,Ny))
	psi_hat = np.fft.fft2(psi)/np.double(Nx*Ny)

	cubic = psi*np.power(abs(psi),2.0)

	n = np.power(abs(psi_hat),2.0) 
	n_dot = np.double(g*np.imag(np.conjugate(psi_hat)*(np.fft.fft2(cubic) / np.double(Nx*Ny))  ))


	for j in range(0,1024):
		for i in range(0,1024):
	
			w = np.power(2.0*np.pi * np.double(i) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(j) / Ly, 2.0) 

		#tau_lin = 2.0*np.pi/w
			ib = np.int_(np.sqrt(w))
			nn[ib] += n[i,j]
			nn_dot[ib] += n_dot[i,j] 
	
			normalise[ib] += 1.0
	
	
			w = np.power(2.0*np.pi * np.double(-i-1) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(j) / Ly, 2.0)
		#tau_lin = 2.0*np.pi/w
			ib = np.int_(np.sqrt(w))
			nn[ib] += n[Nx-i-1,j]
			nn_dot[ib] += n_dot[Nx-i-1,j] 
	
		#eps[ib] += tau_non[Nx-i,j]/tau_lin
			normalise[ib] += 1.0
	
			w = np.power(2.0*np.pi * np.double(i) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(-j-1) / Ly, 2.0)
		#tau_lin = 2.0*np.pi/w
			ib = np.int_(np.sqrt(w))
			nn[ib] += n[i,Ny-j-1]
			nn_dot[ib] += n_dot[i,Ny-j-1] 
	
		#eps[ib] += tau_non[i,Ny-j]/tau_lin
			normalise[ib] += 1.0
	
			w = np.power(2.0*np.pi * np.double(-i-1) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(-j-1) / Ly, 2.0)
		#tau_lin = 2.0*np.pi/w
			ib = np.int_(np.sqrt(w))
			nn[ib] += n[Nx-i-1,Ny-j-1]
			nn_dot[ib] += n_dot[Nx-i-1,Ny-j-1] 
	
		#eps[ib] += tau_non[Nx-i,Ny-j]/tau_lin
			normalise[ib] += 1.0
	
for i in range(1,2*2048):
	t_lin = 2.0*np.pi/(i*i)
	if normalise[i] !=0:
		nn[i] /= normalise[i]
		n_dot[i]/= normalise[i]

	if nn[i] !=0:
		eps[i] = t_lin*(nn_dot[i]/nn[i])
		eps2[i] = g*nn[i]/(2.0*i*i)

axs.plot(np.abs(eps[1:1024]), label=r'$\left|\frac{\tau_l}{\tau_{nl}}\right|$')
axs.plot(eps2[1:1024],label=r'$\sim \frac{E_{nl}}{E_{l}}\sim \frac{g \hat{\psi}^2_k}{2k^2 }$')
axs.set_yscale('log')
axs.set_xscale('log')
axs.set_xlabel(r'$k$')
axs.legend(fontsize=12, loc='best', ncol=1)

'''

im00 = axs[0][0].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im00, ax=axs[0,0])


data = np.reshape(psi150[:,3], (Nx,Ny))
im01= axs[0][1].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im01, ax=axs[0,1])

data = np.reshape((psi150[:,2]**2 +psi150[:,3]**2)**0.5 , (Nx,Ny))
im02= axs[0][2].imshow(data, interpolation='bilinear',cmap='cubehelix', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.5*abs(data).max(), vmin=0 )
fig.colorbar(im02, ax=axs[0,2])

data = np.reshape(psi200[:,2], (Nx,Ny))
im10=axs[1][0].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im10, ax=axs[1,0])

data = np.reshape(psi200[:,3], (Nx,Ny))
im11=axs[1][1].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im11, ax=axs[1,1])

data = np.reshape((psi200[:,2]**2 +psi200[:,3]**2)**0.5 , (Nx,Ny))
im12=axs[1][2].imshow(data, interpolation='bilinear',cmap='cubehelix', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.5*abs(data).max(), vmin=0 )
fig.colorbar(im12, ax=axs[1,2])

data = np.reshape(psi250[:,2], (Nx,Ny))
im20=axs[2][0].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im20, ax=axs[2,0])

data = np.reshape(psi250[:,3], (Nx,Ny))
im21=axs[2][1].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im21, ax=axs[2,1])

data = np.reshape((psi250[:,2]**2 +psi250[:,3]**2)**0.5 , (Nx,Ny))
im22=axs[2][2].imshow(data, interpolation='bilinear',cmap='cubehelix', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.5*abs(data).max(), vmin=0)
fig.colorbar(im22, ax=axs[2,2])

	
axs[0][0].set_xlabel(r'$x$')
axs[0][0].set_ylabel(r'$y$')
axs[1][0].set_xlabel(r'$x$')
axs[1][0].set_ylabel(r'$y$')
axs[2][0].set_xlabel(r'$x$')
axs[2][0].set_ylabel(r'$y$')


axs[0][1].set_xlabel(r'$x$')
axs[0][1].set_ylabel(r'$y$')
axs[1][1].set_xlabel(r'$x$')
axs[1][1].set_ylabel(r'$y$')
axs[2][1].set_xlabel(r'$x$')
axs[2][1].set_ylabel(r'$y$')


axs[0][2].set_xlabel(r'$x$')
axs[0][2].set_ylabel(r'$y$')
axs[1][2].set_xlabel(r'$x$')
axs[1][2].set_ylabel(r'$y$')
axs[2][2].set_xlabel(r'$x$')
axs[2][2].set_ylabel(r'$y$')

axs[0][0].set_title(r'$\textrm{Re}(\psi),\ t=0.15$')
axs[0][1].set_title(r'$\textrm{Im}(\psi),\ t=0.15$')
axs[0][2].set_title(r'$|\psi|,\ t=0.15$')
axs[1][0].set_title(r'$\textrm{Re}(\psi),\ t=0.2$')
axs[1][1].set_title(r'$\textrm{Im}(\psi),\ t=0.2$')
axs[1][2].set_title(r'$|\psi|,\ t=0.2$')
axs[2][0].set_title(r'$\textrm{Re}(\psi),\ t=0.25$')
axs[2][1].set_title(r'$\textrm{Im}(\psi),\ t=0.25$')
axs[2][2].set_title(r'$|\psi|,\ t=.025$')

'''
#plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.tight_layout()
plt.savefig('time_nonlinear.pdf')

