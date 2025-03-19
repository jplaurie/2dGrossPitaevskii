import numpy as np
import matplotlib.pyplot as plt


dt = 1.e-1

Lx = 2.0*np.pi
Ly = 2.0*np.pi
Nx = 2048
Ny = 2048


kmax = Nx
dk = 1.0

N_time = 100
start_file=100


Tmax = 1.0

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)





test = np.zeros((kmax, N_time))
omega_k = np.zeros_like(test,dtype=complex)
normalise = np.zeros(kmax)

fig, axs = plt.subplots()




for time_slice in range(0,N_time,1):
	print('loading file = psi.' + str(start_file+time_slice))	
	filename = '/data2/2d_gp_selfsimilar/ensemble_9/data/psi.%.5d' % (start_file+time_slice);
	data = np.loadtxt(filename)
	psi = data[:,2]+1j*data[:,3]
	psi = np.reshape(psi, (Nx,Ny))
	omega_k[:,time_slice] =  psi[:kmax,0]



omega_k_hat=np.fft.fft2(omega_k)


data_output = omega_k_hat

data_output[:,:(N_time//2)] = omega_k_hat[:,(N_time//2):]
data_output[:,(N_time//2):] = omega_k_hat[:,:(N_time//2)]

data_output[:kmax//2,:] = omega_k_hat[kmax//2:,:]
data_output[kmax//2:,:] = omega_k_hat[:kmax//2,:]




im11=axs.imshow(abs(data_output), interpolation='bilinear',cmap='jet', extent=[0, kmax, 0, 2.0*np.pi/Tmax], vmax = abs(omega_k_hat).max(), vmin=0 )
fig.colorbar(im11, ax=axs)

plt.tight_layout()
plt.savefig('omega_k_plot.pdf')
np.savetxt('omega_k.dat', abs(data_output), fmt='%1.4e')