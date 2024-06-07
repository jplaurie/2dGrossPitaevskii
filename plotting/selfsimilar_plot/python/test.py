

#%% Load Modules
startup_file='init_python.py'
exec(open(startup_file).read())

#number of ensembles
N_ensemble=10

#define the start/end files for plotting
file_start = 150
file_end = 300

#time between files (0.5 factor due to change in equation)
dt = 0.5*1.e-3

#forcing scale
kf=512
wf=kf**2.0
L = 2.0*np.pi
N=1024

#total_wave_action = 168.9616
#total_energy = 1.50605e7

#fluxes defined from initial forcing
wave_action_flux = 2.0 * 32.12
energy_flux = 2.0 * 8.41969e6

#fit regions for T/mu
kmin_low=1
kmax_low = 6

#fit regions for T/\omega
kmin_high=200
kmax_high=400

#if true then fit spectum by max. Default is average across modes
simon_fitting = TRUE


def scaling_fit(x,a1,b1):
    return a1*(x**b1)
	
def scaling_constant_fit(x,a2):
    return a2

def Gibbs_function(x,a):
    return 1.0/( (x/a) + 1.0)


filename = '/data2/2d_gp_selfsimilar/ensemble_0/output/spectrum.00125'
spec_avg = np.loadtxt(filename)
spec_avg = 0.0*spec_avg

n_omega = np.zeros(N)
k = np.range(N)
omega = k**2
# T and mu vs times


T_mu_data=np.zeros((file_end+1,3))
index = 0
for i in range(file_start,file_end+1,1):

	time = i*dt
	spec_avg *= 0.0

	#average over ensembles
	for j in range(0,N_ensemble,1):
		filename = '/data2/2d_gp_selfsimilar/ensemble_%.1d/output/spectrum.%.5d' % (j,i);
		spectrum = np.loadtxt(filename)
		n_omega += spectrum[:,0]
        ]
	
	n_omega /= N_ensemble
	
	
	if i==file_start:
		total_wave_action = np.sum(n_omega[:,1])
		total_energy = 0.5*np.sum(n_omega[:,0]**2.0 * n_omega[:,1])

	w = spec[1:,0]**2; # w = k^2
	n_1d = spec[1:,1]/(2.0*spec[1:,0]); # n_2d/2k = 1d sepctrum in freq space

	if simon_fitting:
		maxN= np.max(n_1d[0:(kf//2)])
		maxN_index = np.argmax(n_1d[0:(kf//2)])
		n_1d /= maxN
		w_fit = w[maxN_index:kf]
		n_fit = n_1d[maxN_index:kf]
    else:




	opt, pcov = curve_fit(Gibbs_function, w_fit, n_fit)

	mu = opt
	T = opt * maxN
	print(time, mu, T)

	T_mu_data[index,0] = i*dt
	T_mu_data[index,1] = T
	T_mu_data[index,2] = mu
	
	index=index+1


fig1, axs1 = newfig(1, 2)

mu0=T_mu_data[0,2]
T0=T_mu_data[0,1]
t0=T_mu_data[0,0]


axs1[0].plot( T_mu_data[:index,0],T_mu_data[:index,1], linestyle='-',linewidth=5)#,label=r'Numerics')	
axs1[1].plot( T_mu_data[:index,0],T_mu_data[:index,2], linestyle='-',linewidth=5)#,label=r'Numerics')	
axs1[1].plot( T_mu_data[:index,0], mu0*np.exp(-wave_action_flux*(T_mu_data[:index,0]-t0)/T0), linestyle='dashed',linewidth=3, color='black')

axs1[0].set_xlim(file_start*dt,file_end*dt)	
axs1[1].set_xlim(file_start*dt,file_end*dt)		

#axs1[0].set_ylim(0,2)	

	
axs1[0].set_xlabel(r'$t$')
axs1[0].set_ylabel(r'$T(t)$')
axs1[1].set_xlabel(r'$t$')
axs1[1].set_ylabel(r'$\mu(t)$')


axs1[1].set_yscale('log')

fig1.tight_layout()

fig1.savefig('test.pdf')
