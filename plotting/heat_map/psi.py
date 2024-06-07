import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


dt = 1.e-1
N=250

L=2.0*np.pi
Nx=2048
Ny=2048


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)


fig, axs = plt.subplots(3, 3,figsize=(15,15))

psi150=np.loadtxt('/data2/2d_gp_selfsimilar/ensemble_1/data/psi.00150')
psi200=np.loadtxt('/data2/2d_gp_selfsimilar/ensemble_1/data/psi.00200')
psi250=np.loadtxt('/data2/2d_gp_selfsimilar/ensemble_1/data/psi.00250')

#print(psi150[2,:])
#print(psi150.shape())

data = np.reshape(psi150[:,2], (Nx,Ny))
im00 = axs[0][0].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im00, ax=axs[0,0])


data = np.reshape(psi150[:,3], (Nx,Ny))
im01= axs[0][1].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im01, ax=axs[0,1])

data = np.reshape((psi150[:,2]**2 +psi150[:,3]**2) , (Nx,Ny))
im02= axs[0][2].imshow(data, interpolation='bilinear',cmap='cubehelix', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.5*abs(data).max(), vmin=0 )
fig.colorbar(im02, ax=axs[0,2])

data = np.reshape(psi200[:,2], (Nx,Ny))
im10=axs[1][0].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im10, ax=axs[1,0])

data = np.reshape(psi200[:,3], (Nx,Ny))
im11=axs[1][1].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im11, ax=axs[1,1])

data = np.reshape((psi200[:,2]**2 +psi200[:,3]**2), (Nx,Ny))
im12=axs[1][2].imshow(data, interpolation='bilinear',cmap='cubehelix', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.5*abs(data).max(), vmin=0 )
fig.colorbar(im12, ax=axs[1,2])

data = np.reshape(psi250[:,2], (Nx,Ny))
im20=axs[2][0].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im20, ax=axs[2,0])

data = np.reshape(psi250[:,3], (Nx,Ny))
im21=axs[2][1].imshow(data, interpolation='bilinear',cmap='jet', extent=[-L/2, L/2, -L/2, L/2], vmax = 0.25*abs(data).max(), vmin=-0.25*abs(data).max() )
fig.colorbar(im21, ax=axs[2,1])

data = np.reshape((psi250[:,2]**2 +psi250[:,3]**2) , (Nx,Ny))
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
axs[0][2].set_title(r'$|\psi|^2,\ t=0.15$')
axs[1][0].set_title(r'$\textrm{Re}(\psi),\ t=0.2$')
axs[1][1].set_title(r'$\textrm{Im}(\psi),\ t=0.2$')
axs[1][2].set_title(r'$|\psi|^2,\ t=0.2$')
axs[2][0].set_title(r'$\textrm{Re}(\psi),\ t=0.25$')
axs[2][1].set_title(r'$\textrm{Im}(\psi),\ t=0.25$')
axs[2][2].set_title(r'$|\psi|^2,\ t=.025$')


#plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.tight_layout()
plt.savefig('heat_map.pdf')

fig, axs2 = plt.subplots(3, 1,figsize=(15,8))

stride = 1024
x = 2.0*np.pi*np.arange(Nx)/Nx
data = data = np.reshape(psi150[:,2]**2 + psi150[:,3]**2, (Nx,Ny))
axs2[0].plot(x,data[stride,:],linewidth=3,color='C0')

data = np.reshape(psi200[:,2]**2 + psi200[:,3]**2, (Nx,Ny))
axs2[1].plot(x,data[stride,:],linewidth=3,color='C1')

data = np.reshape(psi250[:,2]**2 + psi250[:,3]**2, (Nx,Ny))
axs2[2].plot(x,data[stride,:],linewidth=3,color='C2')



axs2[0].set_xlim(0,2.0*np.pi)
axs2[1].set_xlim(0,2.0*np.pi)
axs2[2].set_xlim(0,2.0*np.pi)
axs2[0].set_ylim(0,25)
axs2[1].set_ylim(0,25)
axs2[2].set_ylim(0,25)

axs2[0].set_xlabel(r'$x$')
axs2[1].set_xlabel(r'$x$')
axs2[2].set_xlabel(r'$x$')

axs2[0].set_ylabel(r'$|\psi|^2$')
axs2[1].set_ylabel(r'$|\psi|^2$')
axs2[2].set_ylabel(r'$|\psi|^2$')


plt.tight_layout()
plt.savefig('slice.pdf')