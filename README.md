
<a name="top"></a>

![OS](https://img.shields.io/badge/OS-linux%2C%20macOS-0078D4) 
![language](https://img.shields.io/badge/Language-C%2B%2B-yellow?)

## Table of Contents
- [About](#about)
- [How to Run the Code](#how-to-run-the-code)

## About

This repository contains the numerical code for solving the two-dimensional Gross-Pitaevskii equation for the complex wave function $\psi(\mathbf{x},t)$ in periodic boundary conditions of length $[0,L_x)\times[0,L_y)$. The code solves the equation: 

```math
i\frac{\partial \psi}{\partial t}=-c\nabla^2\psi - \frac{\mu}{2}\psi +\frac{g}{2}|\psi|^2\psi - i\nu(-\nabla^2)^p \psi - i\alpha(-\nabla^2)^{-q} \psi+F,
```

using a pseudo-spectral discretisation method, fully dealiased using the 3/2-rule. Time stepping is performed using the fourth-order exponential time differencing Runge-Kutta method.

### Parameters 
$c$ : speed of sound  
$\mu$ : chemical potential  
$g$ : nonlinearity parameter  
$\nu$ : hyper-viscosity coefficient  
$p$ : power of hyper-viscosity ($p>0$)  
$\alpha$ : hypo-viscosity coefficient  
$q$ : power of hypo-viscosity ($q\geq0$)  
$F$ : additive forcing

### Libraries

The code uses the following libraries that will need to be installed by the user

• gcc (compiler)  
• ﬀtw-3 (fast Fourier transform library)  
• armadillo (C++ linear algebra library)  
• openBLAS (linear algebra library that is used by armadillo)  

All of the above can be installed via the local repositories on Linux, or via macports (and possibly Homebrew) on
macOS.

## How to Run the Code

The numerical code contains a */src/Makefile* that will compile the code using the gcc compiler and create an executable /src/gp2d. At the very top are paths to the default locations of the linked libraries for both Linux and MacOS. Simply comment out the one not in use. Bear in mind that the location of the libraries may diﬀer on your local machine.
All global parameters are stored in the header file /src/const.h. Any changes in /src/const.h will require a recompilation of the numerical code before running. Once /src/gp2d is created, it can be moved to any directory one wishes, but the /data/ folder must be present
in the running directory.


### The /data/ Folder
To run the code, the user must create a /src/data/ directory in the same folder as the created executable. In
the /data/ folder the user must create a file called /data/curframe.dat that will contain two numbers separated by
a tab. The numbers represent the current time and the label of the last generated file number. Each time the code
outputs a new file for the data, /data/curframe.dat will be automatically updated with the new time and label of
the last data file. The code uses this information to when restarting from the last generated file.

| current time | file number |
|:------------:|:-----------:|
|       z      |    XXXXXX   |

The folder ./data/ will also contain all the data for the wave function $\psi$. They are stored in files called
/data/psi.XXXXXX where XXXXXX is the file number padded by zeros. Files /data/psi.XXXXXX are core data for
the numerical code:

| $z$ | $Re(\psi)$ | $Im(\psi)$ |
|:---:|:-------:|:-------:|
| .|. |. |
| .|. |. |
| .|. |. |

### The ./output/ Folder
The code will automatically generate (if not already present) another folder called ./output/ that will contain all
necessary in-code post-processing of the data. Usually, this includes the data used for debugging and verification.
The code will output files at the same intervals as the main data /data/psi.XXXXXX. /output/energy.XXXXXX records the Hamiltonian or energy of the equation to ensure that it is properly conserved by the dynamics.

```math
H = H_{lin} + H_{pot} + H_{non}= \int c\left| \nabla \psi \right|^2 - \frac{\mu}{2} |\psi|^2 + \frac{g}{4}|\psi|^4 d{\bf x}.
```

 It is recorded in the form of a row of four numbers:

| current time | linear energy | potential energy | nonlinear energy |
|:---:|:-------:|:-------:|:---:|
| $z$ |  $H_{lin}$ |  $H_{pot}$ | $H_{non}$ |


/output/wave.XXXXXX records the wave action of the equation to ensure that it is properly conserved by the dynamics:

```math
N = \int  |\psi|^2 d{\bf x}.
```

It is recorded in the form:

| current time |  total wave action | zeroth mode of wave action | 
|:-:|:-:|:-:|
| $z$ |  $N$  | ${\rm abs}\left(\hat{\psi}_{k=0}\right)^2$ | 


/output/spec.XXXXXX records the Fourier wave action spectrum of $|\hat{\psi}_{k}|^2$ and wave energy spectrum $k^2\left|\hat{\psi}_{k}\right|^2$ verses $k$. It is recorded in the form of five columns of numbers:

| $k$ | $abs\left(\hat{\psi}_k\right)^2$ | $c k^2 abs\left(\hat{\psi}_k\right)^2$ | $\langle abs\left(\hat{\psi}_k\right)^2\rangle$ | $\langle ck^2 abs\left(\hat{\psi}_k\right)^2\rangle$ |
|:-:|:-:|:-:|:-:|:-:|
| . | . | . | . | . |

/output/flux.XXXXXX records spectral Fourier flux vs k for both the wave action and energy (instantaneous as
well as time averaged over data since start of code). It is recorded in the form of five columns:


| $k$ | wave action flux $\eta_k$ |  energy flux $\epsilon_k$ | $\langle \eta_k\rangle$ | $\langle \epsilon_k\rangle$ | 
|:-:|:-:|:-:|:-:|:-:|
| . | . | . | . | . | 

/output/dis.XXXXXX records the total wave action and energy dissipated via the terms that constitute D at the time of output. It is recorded in the form of five columns:

| current time | wave action dissipated by hypoviscosity |  wave action dissipated by hyperviscosity |  linear energy dissipated by hypoviscosity |  energy dissipated by hyperviscosity |
|:-:|:-:|:-:|:-:|:-:|
| . | . | . | . | . | 

Using terminal command cat energy.* > energy.tot will append all individual output files into one single file for plotting.


### The ./initial/ Folder
As probably surmised, the code initializes from the data file /data/psi.XXXXXX stored in /data/ to which the file number in /data/curframe.dat points. Therefore, to initialise the code when not restarting one needs to include an initial data file, normally /src/data/psi.000000 (it requires for /data/curframe.dat to be set to 0 0). If you require a zero-state initial condition because you are running a forced/dissipated simulation, you can set XXXXXX in /data/curframe.dat to be any negative integer and the code will interpret this as you wanting to create a zero-state initial condition and will generate this automatically. If you want to generate a custom (bespoke) initial condition, then next to the /src/ folder is the /initial/ folder that contains a simple .cpp file for generating a custom initial state that can be copied into /data/.