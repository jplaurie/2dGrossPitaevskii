# 2dGrossPitaevskii

Numerical code for  2d Gross Pitaevskii equation

```math
i\frac{\partial \psi}{\partial t}=-\frac{1}{2}\Delta\psi - \frac{\mu}{2}\psi +\frac{g}{2}|\psi|^2\psi - i\nu(-\Delta)^p \psi - i\alpha(-\Delta)^{-q} \psi+F
```

**Paramters:** Chemical potential $`\mu`$, nonlinearity $`g`$.

**Boundary Conditions:** Periodic boundary conditions.

**Spatial Discretization:** Pseudo-spectral method with full dealiasing with the 3/2-rule.

**Temporal Discretization:** 4th-order exponential time differencing Runge-Kutta method.

**Data Output:** wavefunction, energy, wave intensity, wave action spectra.

**C++ Libraries** FFTW, armadillo
