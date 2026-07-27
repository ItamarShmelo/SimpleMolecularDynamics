# Project A: Calculation of Material Properties Using Molecular Dynamics (MD)

A material is described by a collection of particles, either molecules or atoms. Assume that the interaction between the particles is the sum of the separate interactions of each particle with every other particle. In other words, the interactions are pairwise.

We will use the Lennard-Jones potential:

$$
\phi(r)=4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right],
$$

where $r$ is the distance between the particles, $\sigma$ is the distance at which the potential is zero, and $\varepsilon$ is the depth of the potential minimum.

The dynamics of the particles are described by solving Newton's equation. Suppose that particle $i$ is located at coordinate $\bar r_i$, and that there are $N$ particles, all with the same mass $m$. The equation of motion is

$$
m\ddot{\bar r}_i = \bar F_i = \sum_{j\ne i}\bar F_{i,j},
$$

where $\bar F_{i,j}$ is the force acting on particle $i$ due to its interaction with particle $j$. This force is calculated from the gradient of the potential:

$$
\bar F_{i,j} = -\bar\nabla_i \phi\left(\left|\bar r_j-\bar r_i\right|\right).
$$

The equation of motion can be solved using finite differences with the following formula, which is simple, stable, and centered in time:

$$
\forall i,\qquad m\frac{\bar r_i^{n+1} -2\bar r_i^n +\bar r_i^{n-1}}{dt^2} = \bar F_i^n.
$$

The superscript denotes the time-step number, $dt$ is the time-step size, and the force on the right-hand side is calculated at step $n$. The unknown in this equation is the position at the next step, $\bar r_i^{n+1}$, which is calculated from the positions at the current step $n$ and the previous step $n-1$.

Let the energy of particle $i$ at step $n$ be denoted by $E_i^n$. The energy is the sum of the potential energy $\phi_i^n$ and kinetic energy $K_i^n$, where

$$
\phi_i^n = \frac{1}{2}\sum_{j\ne i}\phi_{i,j}^n = \frac{1}{2}\sum_{j\ne i}\phi\left(\left|\bar r_j^n-\bar r_i^n\right|\right),
$$

$$
K_i^n = \frac{m}{2}\left|\bar v_i^n\right|^2, \qquad \bar v_i^n = \frac{\bar r_i^{n+1} - \bar r_i^{n-1}}{2dt}.
$$

Suppose that $N$ particles occupy a volume $V$, so that the mass density

$$
\rho=\frac{mN}{V}
$$

is known, and assume periodic boundary conditions. Given an initial state characterized by the position and velocity of every particle, the total energy of the system,

$$
E=\sum_i E_i,
$$

is constant in time.

The fundamental principle of molecular-dynamics calculations is that, if we follow the evolution of a system over a sufficiently long period, time averages of physical quantities describe the ensemble average of many equivalent systems. This is the ergodic hypothesis.

Assuming isotropy, define the instantaneous pressure as

$$
p^n = \frac{1}{D}\frac{1}{V}\left(\sum_{j>i}\bar r_{i,j}^{,n}\cdot\bar F_{i,j}^{,n} + m\sum_i\left(v_i^n\right)^2\right),
$$

where $D$ is the dimensionality of the system, meaning two in a two-dimensional system and three in a three-dimensional system, and

$$
\bar r_{i,j}=\bar r_j-\bar r_i
$$

is the displacement between the pair of particles.

The first term contains a sum over particle pairs and represents the contribution of the interparticle forces to the pressure. The second term represents the contribution of the particle velocities.

By averaging the instantaneous pressure values, we obtain the mean pressure

$$
p=\left\langle p^n\right\rangle_n.
$$

This pressure can be associated with a state having mass density $\rho$, as defined above, and specific total energy, meaning energy per unit mass,

$$
\epsilon=\frac{E}{mN}.
$$

We have therefore obtained a relation between the thermodynamic quantities, namely the equation of state

$$
p=p(\rho,\epsilon).
$$

## Temperature

The temperature can be defined using the mean kinetic energy:

$$
\frac{D}{2}Nk_B T = \left\langle K^n\right\rangle_n = \left\langle\sum_i K_i^n\right\rangle_n,
$$

where $k_B$ is Boltzmann's constant. This expression is valid when the mean velocity in every direction is zero.

## Two-dimensional periodic domain

We will consider a calculation in which the particles are located inside a square in a two-dimensional system. Assume that the $x$ and $y$ coordinates lie in the interval

$$
[0,L].
$$

Periodic boundary conditions are implemented by examining the coordinates of every particle after calculating its new position. For example, for the $x$-coordinate of particle $i$,

$$
x_i>L \quad\Longrightarrow\quad x_i\rightarrow x_i-L,
$$

and

$$
x_i<0 \quad\Longrightarrow\quad x_i\rightarrow x_i+L.
$$

The boundary conditions must also be taken into account when calculating the shortest distance between particles $i$ and $j$. For each component of the vector $\bar r_{i,j}$, for example $x_{i,j}$,

$$
x_{i,j}>\frac{L}{2} \quad\Longrightarrow\quad x_{i,j}\rightarrow x_{i,j}-L,
$$

and

$$
x_{i,j}<-\frac{L}{2} \quad\Longrightarrow\quad x_{i,j}\rightarrow x_{i,j}+L.
$$

When calculating the velocity from the difference between $\bar r_i^{n+1}$ and $\bar r_i^{n-1}$, make sure that the position difference does not contain an artificial displacement of one coordinate by $L$.

For example, if a particle passes from one side of the square to the other, its velocity must be calculated as though it remained beyond the system boundary on its original side. You must also make sure that the equation of motion does not contain artificial position differences caused by the implementation of the periodic boundary conditions.

## Dimensionless quantities

The calculations may be performed using dimensionless quantities. The dimensionless density, corresponding to the number of particles in a square whose side length is $\sigma$, is

$$
\rho^*=\frac{\rho\sigma^2}{m}.
$$

The average energy per particle, in units of $\varepsilon$, is

$$
\epsilon^*=\frac{\epsilon m}{\varepsilon}.
$$

The density, total number of particles, and system dimensions are related by

$$
L = \sqrt{\frac{Nm}{\rho}} = \sigma\sqrt{\frac{N}{\rho^*}}.
$$

The dimensionless time step is

$$
dt^* = \frac{dt}{\sigma}\sqrt{\frac{\varepsilon}{m}}.
$$

In the example considered here, the particles have mass

$$
m=8.3\times10^{-23}\ \mathrm{g},
$$

corresponding to an atomic weight of approximately 50.

The distance at which the potential vanishes is two angstroms,

$$
\sigma=2\times10^{-8}\ \mathrm{cm},
$$

and the depth of the potential well is

$$
\varepsilon=6.9\times10^{-13}\ \mathrm{erg},
$$

corresponding to a temperature of approximately $5000\ \mathrm{K}$.

Because we calculate the interaction between every pair of particles, rather than using neighbor tracking as is done in professional MD simulations, assume a relatively small number of particles:

$$
N=30\ldots40.
$$

A typical dimensionless density is

$$
\rho^*=1\ldots1.5,
$$

which gives a square of approximate side length

$$
L\approx5\sigma.
$$

A typical dimensionless time step for which the changes remain small is

$$
dt^*=10^{-3}.
$$

## Modified Lennard-Jones potential

For computational convenience, we will modify the Lennard-Jones potential so that it becomes exactly zero at a sufficiently large distance and so that the repulsion at short distances is moderated.

Let

$$
x=\frac{r}{\sigma}.
$$

For

$$
x\leq0.85,
$$

use

$$
\phi = \varepsilon\left(-39170.516264x^3 +102795.709007x^2 -90115.812549x +26401.649473\right).
$$

For

$$
0.85<x\leq1.24445506,
$$

use

$$
\phi = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right].
$$

For

$$
1.24445506<x\leq1.737051787=r_m,
$$

use

$$
\phi = \varepsilon\left[3.2920028(r_m-x)^3 - 4.86489008(r_m-x)^2\right].
$$

For

$$
x>r_m,
$$

use

$$
\phi=0.
$$

# Assignment A

Determine the zero-temperature energy and zero-temperature pressure as functions of density. In other words, calculate the energy and pressure at a specified density and at a temperature of zero kelvin, for

$$
1.0\leq\rho^*\leq1.5
$$

in increments of $0.05$.

## Method

After choosing $\rho$, $L$, and $N$, find the zero-temperature configuration, meaning the particle positions at zero temperature.

This can be done by minimizing the energy. The simplest method is a relaxation method. Generate random initial particle positions and perform time steps according to Newton's formula, but set the velocities to zero after every step.

Because the velocity does not explicitly appear in the equation of motion, this means that, when calculating the position at step $n+1$, assume that the position at step $n-1$ is equal to the position at step $n$.

Stop the process when the particle displacements become negligible. To calculate the mean pressure, perform additional time steps in the final state until the pressure stabilizes.

Present the following:

* A table containing the energy and pressure obtained in the different calculations.
* A graph of $p(\rho,T=0)$ in c.g.s. units.
* The particle positions obtained for several of the densities.
* An explanation of the results.

# Assignment B

Determine the dependence of pressure and energy on temperature at the fixed density

$$
\rho^*=1.1.
$$

Perform calculations for ten energy values. Choose five values between the minimum energy found in Assignment A and

$$
\epsilon^*=10,
$$

and another five values in the range

$$
10<\epsilon^*\leq100.
$$

## Method

After choosing $\rho$, $L$, and $N$, find the zero-temperature configuration, or use the particle positions calculated in Assignment A.

Next, generate random particle velocities such that the system has the specified total energy.

Perform time steps to allow the velocities to relax and to obtain a velocity field distributed according to the Maxwell-Boltzmann distribution. During this stage, evolve the system using Newton's equation without resetting the velocities.

After the velocity-relaxation stage, continue performing time steps to calculate the mean pressure and temperature, until these quantities stabilize.

Present the following:

* A table containing the energy, pressure, and temperature obtained in the different calculations.
* A log-log graph of $p(\rho^*=1.1,T)$ in c.g.s. units.
* The particle positions at several different times during the averaging stage.
* An explanation of the results.

# General Notes

* It is recommended to test the sensitivity of the results by carrying out several calculations for a given $\rho^*$, while changing $N$ and $L$ appropriately.
* A reasonable number of time steps for the different calculations is between $20{,}000$ and $200{,}000$.
* When calculating the kinetic energy for the purpose of determining temperature or pressure, subtract the mean velocity. Alternatively, subtract the mean velocity from all particle velocities after generating the initial velocities. Newton's equation conserves momentum. When necessary, the mean velocity may also be removed during the calculation if numerical-error accumulation produces a drift in the mean velocity.
