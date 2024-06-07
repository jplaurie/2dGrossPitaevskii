import numpy as np
import matplotlib.pyplot as plt



g = 2000.0

Lx = 2.0*np.pi
Ly = 2.0*np.pi
Nx = 256
Ny = 256

dk = 2.0*np.pi/Lx

dt = 1.e-3
Nt = 1024
start_file = 650

#theoretical
#wmax = 2.0*np.pi/dk
#wmin = 2.0*np.pi/kmax

wmax = np.pi/dt
dw = 2.0*np.pi/(dt*Nt)


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)



psi_kx_t = np.zeros([Nx, Nt], dtype=complex)
psi_ky_t = np.zeros([Nx, Nt], dtype=complex)
data_output = np.zeros([Nx//2+1, Nt])
abs_data_output = np.zeros([Nx//2+1, Nt//2+1])


fig, axs = plt.subplots()




for time_slice in range(0,Nt,1):
	print('loading file = psi.' + str(start_file+time_slice))	
	filename = '../../src/data/psi.%.5d' % (start_file+time_slice);
	data = np.loadtxt(filename)
	psi_data = data[:,2]+1j*data[:,3]
	psi = np.reshape(psi_data, (Nx,Ny))
	psi_hat = np.fft.fft2(psi)/np.double(Nx*Ny)

	psi_kx_t[:,time_slice] = psi_hat[:,0]
	psi_ky_t[:,time_slice] = psi_hat[0,:]


psi_kx_w = np.fft.ifft(psi_kx_t,axis=1)/np.double(Nt)
psi_ky_w = np.fft.ifft(psi_ky_t,axis=1)/np.double(Nt)



for i in range(1,(Nx//2)):
	for j in range(0,Nt//2):
		data_output[i,j] = np.power(abs(psi_kx_w[i,(Nt//2) + j]),2.0) + np.power(abs(psi_kx_w[Nx-i,(Nt//2) + j]),2.0)
		data_output[i,j] += np.power(abs(psi_ky_w[i,(Nt//2) + j]),2.0) + np.power(abs(psi_ky_w[Nx-i,(Nt//2) + j]),2.0)

		data_output[i,(Nt//2)+j] +=  np.power(abs(psi_kx_w[i, j]),2.0) + np.power(abs(psi_kx_w[Nx-i,j]),2.0)
		data_output[i,(Nt//2)+j] += np.power(abs(psi_ky_w[i, j]),2.0) + np.power(abs(psi_ky_w[Nx-i,j]),2.0)



for j in range(0,Nt//2):
	data_output[0,j] = np.power(abs(psi_kx_w[0,(Nt//2) + j]),2.0) 
	data_output[0,(Nt//2)+j] += np.power(abs(psi_kx_w[0, j]),2.0) 

	data_output[(Nx//2),j] = np.power(abs(psi_kx_w[(Nx//2),(Nt//2) + j]),2.0) 
	data_output[(Nx//2),(Nt//2)+j] += np.power(abs(psi_kx_w[(Nx//2), j]),2.0) 



abs_data_output[:,:] = 0.5*data_output[:,((Nt//2)-1):]
abs_data_output[:,1:((Nt//2)+1)] += 0.5*data_output[:,(Nt//2):0:-1]

print("dw = " +  str(dw))
print("wmax = " +  str(wmax))
#im11=axs.imshow(abs(data_output), interpolation='bilinear',cmap='jet', extent=[0, kmax, dw, wmax], vmax = abs(data_output).max(), vmin=0 )
#fig.colorbar(im11, ax=axs)

#plt.tight_layout()
#plt.savefig('omega_k_plot.pdf')
np.savetxt('abs_frequency_spectrum.dat', abs(abs_data_output), fmt='%1.4e')
np.savetxt('frequency_spectrum.dat', abs(data_output), fmt='%1.4e')