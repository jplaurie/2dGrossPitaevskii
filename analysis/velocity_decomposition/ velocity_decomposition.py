import numpy as np
import matplotlib.pyplot as plt

Nx = 2048
Ny = 2048
Lx = 2.0*np.pi
Ly = 2.0*np.pi

g=2000.0
dt = 1.e-3


file_start =1
file_end = 100

def define_k():
	kx = np.zeros((Nx,Ny),dtype=complex128)
	ky = np.zeros((Nx,Ny),dtype=complex128)

	for i in range(Nx/2):
		for j in range(Ny/2):
			kx[i,j] = 1j* 2.0 * np.pi * i / Lx
			kx[Nx-i-1,j] = 1j* 2.0 * np.pi * (Nx-i) / Lx
			kx[i,Ny-j-1] = 1j* 2.0 * np.pi * i / Lx
			kx[Nx-i-1,Ny-j-1] = 1j* 2.0 * np.pi * (Nx-i) / Lx

			ky[i,j] = 1j* 2.0 * np.pi * j / Ly
			ky[Nx-i-1,j] = 1j* 2.0 * np.pi * j / Ly
			ky[i,Ny-j-1] = 1j* 2.0 * np.pi * (Ny-j) / Ly
			ky[Nx-i-1,Ny-j-1] = 1j* 2.0 * np.pi * (Ny-j) / Ly
	return kx,ky

def compute_velocity_vorticity():

	psi_hat = np.fft.fft2(psi)
	kx, ky = define_k()

	u = np.imag(np.fft.ifft2(kx*psi_hat)*np.conj(psi))
	v = np.imag(np.fft.ifft2(ky*psi_hat)*np.conj(psi))	


	return u,v, omega




kmax = Nx//8
kmin = 1.0

N_time = 200
start_file=100


for file in range(file_start,file_end+1,1):
	print('loading file', file ))	
	filename = '/data2/2d_gp_selfsimilar/ensemble_9/data/psi.%.5d' % (file);
	data = np.loadtxt(filename)
	psi = data[:,2]+1j*data[:,3]
	psi = np.reshape(psi, (Nx,Ny))
	
	u,v = compute_velocity(psi)
	
	
	
	
	psi_hat = np.fft.fft2(psi)/np.double(Nx*Ny)



#theorkmaxetical
#wmax = 2.0*np.pi/dk
#wmin = 2.0*np.pi/kmax

wmax = 2.0*np.pi/dt
wmin = 2.0*np.pi/(dt*N_time)


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
	psi_hat = np.fft.fft2(psi)/np.double(Nx*Ny)


	for j in range(0,kmax):
		for i in range(0,kmax):
	
			kk = np.power( np.power(2.0*np.pi * np.double(i) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(j) / Ly, 2.0),0.5) 
			ib = np.int_(kk/dk)
			if(ib < kmax):
				omega_k[ib,time_slice] += psi_hat[i,j]
				normalise[ib] += 1.0
	
	
			kk = np.power(np.power(2.0*np.pi * np.double(-i-1) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(j) / Ly, 2.0),0.5)
			ib = np.int_(kk/dk)
			if(ib < kmax):
				omega_k[ib,time_slice] += psi_hat[Nx-i-1,j]		
				normalise[ib] += 1.0
	
			kk = np.power(np.power(2.0*np.pi * np.double(i) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(-j-1) / Ly, 2.0),0.5)
			ib = np.int_(kk/dk)
			if(ib < kmax):
				omega_k[ib,time_slice] += psi_hat[i,Ny-j-1]
				normalise[ib] += 1.0
	
			kk = np.power(np.power(2.0*np.pi * np.double(-i-1) / Lx, 2.0)+ np.power(2.0*np.pi * np.double(-j-1) / Ly, 2.0),0.5)
			ib = np.int_(kk/dk)
			if(ib < kmax):
				omega_k[ib,time_slice] += psi_hat[Nx-i-1,Ny-j-1]
				normalise[ib] += 1.0

	for i in range(0,kmax):
		if normalise[i] !=0:
			omega_k[i,time_slice] /= normalise[i]


omega_k_hat = np.fft.fft(omega_k,axis=1)


data_output = omega_k_hat

data_output[:,:(N_time//2)] = omega_k_hat[:,(N_time//2):]
data_output[:,(N_time//2):] = omega_k_hat[:,:(N_time//2)]


im11=axs.imshow(abs(data_output), interpolation='bilinear',cmap='jet', extent=[0, kmax, wmin, wmax], vmax = abs(omega_k_hat).max(), vmin=0 )
fig.colorbar(im11, ax=axs)

plt.tight_layout()
plt.savefig('omega_k_plot.pdf')
np.savetxt('omega_k.dat', abs(data_output), fmt='%1.4e')



