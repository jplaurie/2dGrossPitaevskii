#!/home/sthalabard/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:21:33 2022

Analysis of forced and DECAYING 2d NLS

under the Gibbsian hypothesis

Cf the paper Stairway to UV catastrophe

@author: sthalabard
"""


#%% Load Modules
startup_file='init_python.py'
exec(open(startup_file).read())



#%% LOAD DATA
nu=0
FORCING=True
if FORCING:
    suffix='forced_%0.2e' %(nu,)
else:
    suffix='free_%0.2e' %(nu,)
IO='tmp'
os.makedirs(IO,exist_ok=True)
name=os.path.join(IO,'data_%s.dill' %(suffix,))
evo=load_dill(name)

#%% Miscellaneous functions
def func(x, x0):
        return 1/(x/x0+1)

def gibbs(x,T,mu):
        return T/(x+mu)

def gibbs_red(x,x0):
    return gibbs(x,x0,x0)

def gibbs_fit(n=None,k=None,tol=1e-30):
    from scipy.optimize import curve_fit

    x,y=k.copy(),n.copy()
    ymax=y.max()
    imax=y.argmax()
    ix=np.flatnonzero(y/ymax>tol)

    res=dic2struc()
    res.wm=x[ix[0]]
    res.wp=x[ix[-1]]
    res.wmax=x[imax]

    xdata=x[ix]
    ydata=y[ix]
    popt, pcov = curve_fit(gibbs_red, xdata, ydata/ydata.max(),bounds=(res.wmax, res.wp))

    res.x0=popt[0]
    res.mu=res.x0
    res.T=ymax*res.x0
    return res

def compute_logsum(field=None,k=None):
     deta=np.log(k[1]/k[0])
     return (field*k*deta).sum()


#%% Fig 1
fig,ax=newfig(1,1,num='Fig 1 :' + suffix )

ax.grid()
ax.set_xlabel('$\\omega$')
ax.set_ylabel('$N(\\omega,t)$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks(10**np.arange(0,16,2))
ax.tick_params(which='both', width=1,direction='in')
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4, color='k')
ax.set_ylim(1e-10,1e0)
ax.set_xlim(1e0,1e16)

if FORCING:
    ax.annotate(text='$\leftarrow f_0$',xy=(1e13,1e-4),c='firebrick')
    ax.axvspan(0.9e13,1.1e13,color='firebrick',alpha=0.2)
    ax.annotate(text='Gibbs',xy=(1,5e-2),c='k')

if not FORCING:
    ax.annotate(text='$t= 0$',xy=(1e13,4e-5),c='firebrick')
    ax.annotate(text='Gibbs',xy=(1,5e-2),c='k')

f=fig.canvas
data,=ax.plot([],[],c='royalblue',lw=2)

for i in range(evo.M):
    k,mu,T,n=evo.k,evo.mu[i],evo.T[i],evo.n[i,:]
    data.set_data(k,n)
    ax.plot(k,n,lw=0.5,c='royalblue')

ax.plot(evo.k,T/(evo.k+mu),lw=2,c='k',ls='--')

name=os.path.join(IO,'1_%s.png' %(suffix,))
fig.savefig(name)


#%% Fig 2 Rescaled spectrum Middle with Energy Inset
fig,ax=newfig(1,1,num='Fig 2 :' + suffix )
axz = inset_axes(ax, width="50%", height="45%", loc='lower left',borderpad=3)
axs=[ax,axz]
for ax in axs:
    ax.grid()

ax=axs[0]
ax.set_ylabel('$\\phi:=(\\mu/T) N(\\omega,t)$')
ax.set_xlabel('$\\eta:=\\omega/\\mu(t)$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-8,1e1)
ax.set_xticks(10**(1.*np.arange(-8,8,2)))
ax.set_xlim(1e-7,1e7)
ax.annotate(text='Gibbs: $\\frac{1}{\\eta+1}$',xy=(2e-7,2))

ax=axs[1]
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylim(0,1.1)
ax.set_xticks(10**(1.*np.arange(-8,8,2)))
ax.set_xlim(1e-7,1e7)
ax.tick_params(axis='both', which='major', labelsize=25)

f=fig.canvas

data0,=axs[0].plot([],[],c='royalblue',lw=2)
data1,=axs[1].plot([],[],c='royalblue',lw=2)
data2,=axs[1].plot([],[],c='royalblue',lw=2)

ax=axs[0]
ax.plot(evo.k,1/(evo.k+1),lw=3,c='k',ls='--',label='$\\frac{1}{\\eta+1}$')

ax=axs[1]
ax.plot(evo.k,1/(evo.k+1),lw=2,c='k',ls='--',label='$\\frac{1}{\\eta+1}$')
ax.plot(evo.k,evo.k/(evo.k+1),lw=2,c='k',ls='--',label='$\\frac{1}{\\eta+1}$')

for i in range(evo.M):
    k,mu,T,n=evo.k,evo.mu[i],evo.T[i],evo.n[i,:]
    etam=evo.target[i]/mu #left front
    data0.set_data(k/mu,(mu/T)*n)
    data1.set_data(k/mu,(mu/T)*n)
    data2.set_data(k/mu,(mu/T)*n*(k/mu))
    f.draw()    
    f.flush_events()
    if etam>1e-3:continue #only show profiles if eta_- <1e-3
    ax=axs[0]
    ax.plot(k/mu,(mu/T)*n,lw=0.4,c='royalblue')
    ax=axs[1]
    ax.plot(k/mu,(mu/T)*n,lw=0.4,c='royalblue')
    ax.plot(k/mu,(mu/T)*n*(k/mu),lw=0.4,c='royalblue')

ax=axs[1]
ax.annotate(text='$\\phi$',xy=(2e-5,0.5))
ax.annotate(text='$\\eta \\phi$',xy=(1e3,0.5))

name=os.path.join(IO,'2_%s.png' %(suffix,))
fig.savefig(name)

#%% Fig 3a: EXPONENTIAL Ansatz 
fig,axs=newfig(1,2,num='3a')
axz = inset_axes(axs[0], width="35%", height="35%", loc='center right',borderpad=0.6)
axs=[axs[0],axs[1],axz]

ax=axs[2]
ax.grid()
ax.set_xlabel('$t$',labelpad=-26)
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_xlim(0,evo.t[-1])
ax.set_ylim(0,evo.t[-1])

ax.set_xscale('linear')
ax.set_yscale('linear')

wf=1e13
if FORCING:    
    tmp=evo.param_forcing['ampli']*np.exp(-0.5*(np.log(evo.k)-np.log(evo.param_forcing['k0']))**2/evo.param_forcing['sig']**2)
#    tmp=evo.param_forcing['ampli']
    QN=compute_logsum(field=tmp,k=evo.k)
    QE=compute_logsum(field=evo.k*tmp,k=evo.k)
    ax.plot(evo.t,evo.t,'--',c='k',lw=2)
    ax.plot(evo.t[::2],evo.N[::2]/QN,'x',c='royalblue',ms=10,label='$N/I$')
    ax.plot(evo.t[1::2],evo.E[1::2]/QE,'o',c='royalblue',ms=10,label='$E/(\\omega_0 I )$')
    coeffT=QN
    labelT='$T/I$'
ax.legend(labelspacing=0.5,fontsize=20)

ax=axs[0]
ax.grid()
ax.set_ylabel('$\\eta_\\pm$')
ax.set_xlabel('$t$')
ax.set_yscale('log')
ax.set_xlim(0,1.1*evo.t[-1])
ax.set_ylim(1e-7,1e7)
ax.set_yticks(10**(1.*np.arange(-6,8,2)))

ax.set_yscale('log')
tfit=16
ix=np.nonzero(evo.t>tfit)

p=np.polyfit(evo.t[ix],np.log(evo.wp[ix]/evo.mu[ix]),1)
c=p[0]
Ap=np.exp((np.log(evo.wp[ix]/evo.mu[ix])-c*evo.t[ix]).mean())
Am=np.exp((np.log(evo.wm[ix]/evo.mu[ix])+c*evo.t[ix]).mean())

t0=np.linspace(0,30,1000)
fitp=lambda t: Ap*np.exp(c*t)
fitm=lambda t: Am*np.exp(-c*t)
ax.plot(t0,fitp(t0),ls='--',lw=2,c='k')
ax.plot(evo.t,evo.wp/evo.mu,'o',c='darkviolet',ms=16,label='$\\eta_+$',markeredgecolor='k')

ax.plot(t0,fitm(t0),ls='--',lw=2,c='k')
ax.plot(evo.t,evo.wm/evo.mu,'o',c='firebrick',ms=16,label='$\\eta_-$',markeredgecolor='w')
ax.legend(ncol=1,loc='center left',labelspacing=2,bbox_to_anchor=(0.2,0.5))
ax.annotate('$A_+e^{ct}$',xy=(3,1000))
ax.annotate('$A_-e^{-ct}$',xy=(3,1e-4))

ax.set_yticks(10**(1.*np.arange(-6,8,2)))

#%
ax=axs[1]
ax.grid()
ax.set_ylabel('$\\mu/\\omega_0\;\;\;T/T_\\infty$')
ax.set_xlabel('$t$')
ax.set_yscale('log')
ax.set_xscale('linear')

ax.set_xlim(0,1.1*evo.t[-1])
ax.set_ylim(1e-5,5e0)
ax.set_yticks(10**(1.*np.arange(-5,1)))

if FORCING:
    Too=QN/c
else:
    Too=1
ax.plot(evo.t,evo.T/Too,'o',c='royalblue',ms=16,markeredgecolor='k',label='$T/T_\\infty$')
ax.axhline(1,lw=2,c='k',ls='--')
ax.annotate('$T_\\infty=I/c $', xy=(5,1.3))
ax.plot(evo.t,evo.mu/wf,'P',c='royalblue',ms=16,label='$\\mu/\\omega_0$',markeredgecolor='k')
gam=c/Ap
gamf=np.exp((np.log(evo.mu[ix]/wf)+c*evo.t[ix]-np.log(evo.t[ix])).mean())
fitmu=lambda t: gamf*t*np.exp(-c*t)
ax.plot(t0,fitmu(t0),ls='--',lw=2,c='k')
ax.annotate('$\\gamma t e^{-c t}$ \n\n $\\gamma= %0.1f \, c/A_+$'%(gamf/gam,),xy=(5,2e-4))
ax.legend(loc='center right',labelspacing=3)

if FORCING:
    name=os.path.join(IO,'3a_%s.png' %(suffix,))
    fig.savefig(name)


#%% Fig 3b: Stretched Exponential

fig,axs=newfig(1,2,num='3b')
axz = inset_axes(axs[0], width="35%", height="35%", loc='center right',borderpad=0.6)
axs=[axs[0],axs[1],axz]

ax=axs[2]
ax.grid()
ax.set_xlabel('$t$',labelpad=-5)
#ax.set_ylabel('$N,E$',labelpad=-26)
ax.tick_params(axis='both', which='major', labelsize=25)
ax.set_xlim(0,evo.t[-1])
ax.set_ylim(0,1.2)
if not FORCING:    
    wf=1e13#QE/QN
    ax.plot(evo.t[::2],evo.N[::2]/evo.N[0],'x',c='royalblue',ms=10,label='$N/N_0$')
    ax.plot(evo.t[1::2],evo.E[1::2]/evo.E[0],'o',c='royalblue',ms=10,label='$E/(\\omega_0 N_0 )$')

ax.set_xscale('linear')
ax.set_yscale('linear')
ax.legend(labelspacing=0.5,fontsize=20)
ax=axs[0]
ax.grid()
ax.set_ylabel('$\\eta_\\pm$')
ax.set_xlabel('$t$')
ax.set_yscale('log')
ax.set_xscale('linear')
ax.set_ylim(1e-7,1e7)
ax.set_xlim(0,evo.t[-1])
ax.set_yscale('log')
tfit=12
ix=np.nonzero(evo.t>tfit)
alpha=1/3

ydata=np.log(evo.wp[ix]/evo.mu[ix])
xdata=evo.t[ix]**alpha
p=np.polyfit(xdata,ydata,1)
Bp=np.exp(p[1])
b=p[0]
fitp=lambda t: Bp*np.exp(b*t**alpha)

ydata=np.log(evo.wm[ix]/evo.mu[ix])
xdata=evo.t[ix]**alpha

Bm=np.exp((ydata+b*xdata).mean())
fitm=lambda t: Bm*np.exp(-b*t**alpha)

t0=np.linspace(1e-10,30,10000)

ax.plot(t0,fitp(t0),ls='--',lw=2,c='k')
ax.plot(evo.t,evo.wp/evo.mu,'o',c='darkviolet',ms=16,label='$\\eta_+$',markeredgecolor='k')
ax.plot(t0,fitp(t0),ls='--',lw=2,c='k')

ax.plot(t0,fitm(t0),ls='--',lw=2,c='k')
ax.plot(evo.t,evo.wm/evo.mu,'o',c='firebrick',ms=16,label='$\\eta_-$',markeredgecolor='w')
ax.legend(ncol=1,loc='center left',labelspacing=2,bbox_to_anchor=(0.2,0.5))
ax.annotate('$B_+e^{bt^{\\alpha}}$',xy=(3,1e6))
ax.annotate('$B_-e^{-bt^{\\alpha}}$',xy=(3,1e-6))

ax.set_yticks(10**(1.*np.arange(-6,7,2)))

#%
ax=axs[1]
ax.grid()
ax.set_ylabel('$\\mu/\\omega_0\;\;\;T/T_\\infty^\\prime$')
ax.set_xlabel('$t$')
ax.set_yscale('log')
ax.set_xscale('linear')

ax.set_xlim(0,1.1*evo.t[-1])
ax.set_ylim(1e-5,5e0)
ax.set_yticks(10**(1.*np.arange(-5,1)))
Too=evo.N[0]/b
ax.plot(evo.t,evo.T/Too,'o',c='royalblue',ms=16,markeredgecolor='k',label='$T/T_\\infty^\\prime$')
ax.plot(t0,t0**(-alpha),lw=2,c='k',ls='--')
ax.annotate('$T_\\infty^\\prime=N/b $', xy=(15,1.3))
ax.annotate('$t^{-\\alpha} $', xy=(5,1.3))

ax.plot(evo.t,evo.mu/wf,'P',c='royalblue',ms=16,label='$\\mu/\\omega_0$',markeredgecolor='k')
fitmu=lambda t,g: g*(b/Bp)*t**alpha*np.exp(-b*t**alpha)
gamf=np.exp((np.log(evo.mu[ix]/wf)+b*evo.t[ix]**alpha-alpha*np.log(evo.t[ix])).mean())

ax.plot(t0,fitmu(t0,5.5),ls='--',lw=2,c='k')
ax.annotate('$\\gamma^\\prime t^{\\alpha} e^{-b t^\\alpha}$ \n\n $\\gamma^\\prime= %0.1f \, b/B_+$' %(5.5,),xy=(3,2e-5))

ax.legend(loc='center right',labelspacing=3)

if not FORCING:
    name=os.path.join(IO,'3b_%s.png' %(suffix,))
    fig.savefig(name)
