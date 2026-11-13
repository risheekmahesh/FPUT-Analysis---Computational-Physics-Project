 Numerical Investigation of Energy Transfer in a Nonlinear Oscillator Chain  
Group Members 
Adhvaith Sandeep    -     EP25BTECH11002
Risheek Mahesh       -     EP25BTECH11024
Akshat Verma           -     EP25BTECH11004

---
# Abstract
The Fermi-Pasta-Ulam-Tsingou (FPUT) problem is a classical problem in nonlinear physics that investigates the motion of a one-dimensional chain of coupled oscillators. It was originally studied to understand whether a nonlinear system of interacting particles would eventually distribute its energy equally among all its normal modes.

In this project, we numerically investigate the dynamics of a fixed-end oscillator chain with nonlinear spring interactions. The system is modelled using the FPUT potential, which includes harmonic, cubic, and quartic contributions. Since the resulting equations of motion are coupled and nonlinear, obtaining their complete time-dependent solutions analytically is difficult. We therefore implement a numerical simulation using the Velocity Verlet integration method.

The initial condition excites the first normal mode of the chain. We investigate how the energy initially concentrated in this mode is distributed among other modes as the system evolves. We also examine the conservation of total energy, calculate normal-mode energy fractions, and use spectral entropy to quantify the distribution of energy among the modes.

Finally, we compare the harmonic, FPUT-alpha, FPUT-beta, and combined nonlinear models to investigate how different nonlinear interactions influence energy transfer in the system.

---

# 1. Introduction and Description of the Physics Problem

## Background

The Fermi-Pasta-Ulam-Tsingou (FPUT) problem is a problem in computational physics and nonlinear dynamics. It started from an experiment done by Enrico Fermi, John Pasta, Stanislaw Ulam and Mary Tsingou in the 1950s.

The original goal was to see how a system of many interacting particles moves toward thermal equilibrium. According to mechanics, a system that interacts enough should eventually spread its energy among all its available degrees of freedom.

To study this, the researchers used a chain of particles linked by nonlinear springs. They began by putting energy into one mode and then used numbers to follow how the energy changed over time.

An unexpected result was observed. Instead of immediately spreading across all the normal modes, the energy showed a tendency to return to the initially excited mode. This behaviour became known as the FPUT recurrence phenomenon.

The FPUT problem shows that even a simple set of coupled oscillators can have dynamic behaviour when nonlinear interactions are added.

## Physical Model

We consider a one-dimensional chain consisting of (N) identical moving particles connected by (N+1) springs.
Each particle has mass (m), and the springs have a linear spring constant (K). The two endpoints of the chain are fixed, while the particles between them are allowed to move.
Let $(x_i)$ denote the displacement of the (i)-th particle from its equilibrium position.
The fixed boundary conditions are
$$  x_0=0,\qquad x_{N+1}=0$$
The extension of the spring between particles (i) and (i+1) is
$$  
r_i=x_{i+1}-x_i  
$$
Therefore, the system contains (N+1) spring extensions:
$$r_0,r_1,\ldots,r_N$$
The potential energy of an individual spring is described by
$$V(r_i)= 
\frac{K}{2}r_i^2+  
\frac{\alpha}{3}r_i^3+
\frac{\beta}{4}r_i^4$$
The three terms represent:
- Harmonic contribution: $(\frac{K}{2}r_i^2)$
- Cubic nonlinear contribution: $(\frac{\alpha}{3}r_i^3)$
- Quartic nonlinear contribution: $(\frac{\beta}{4}r_i^4)$

Here, $(\alpha)$ and $(\beta)$ determine the strengths of the nonlinear interactions.
We can obtain different models by choosing these coefficients appropriately:

| Model          | Cubic coefficient | Quartic coefficient |
| -------------- | ----------------- | ------------------- |
| Harmonic chain | $(\alpha=0)$      | $(\beta=0)$         |
| FPUT-alpha     | $(\alpha\ne 0)$   | $(\beta=0)$         |
| FPUT-beta      | $(\alpha=0)$      | $(\beta\ne 0)$      |
| Combined model | $(\alpha\ne 0)$   | $(\beta\ne 0)$      |

## Equations of Motion
The total potential energy of the chain is the sum of the potential energies of all the springs:
$U=\sum_{i=0}^{N}V(r_i).$
Substituting the spring potential gives
$$U=  
\sum_{i=0}^{N} 
[ 
\frac{K}{2}r_i^2+  
\frac{\alpha}{3}r_i^3+ 
\frac{\beta}{4}r_i^4] 
$$
The kinetic energy of the moving particles is
$$T=\frac{1}{2}m\sum_{i=1}^{N}v_i^2 $$
where $(v_i=\dot{x}_i)$ is the velocity of the $(i)$-th particle.
The total Hamiltonian, or total mechanical energy, is $H=T+U$ 
Using Newton's second law $m\ddot{x}_i=F_i$  
where the force on the $(i)$-th particle is obtained from the potential energy: $F_i=-\frac{\partial U}{\partial x_i}$  
Since each particle is connected to a spring on its left and a spring on its right, the resulting force is
$$F_i=  
K(r_i-r_{i-1})  
+\alpha(r_i^2-r_{i-1}^2)  
+\beta(r_i^3-r_{i-1}^3) $$
Therefore, the equation of motion becomes
$$  
\boxed{ 
m\ddot{x}_i= 
K(r_i-r_{i-1})  
+\alpha(r_i^2-r_{i-1}^2)  
+\beta(r_i^3-r_{i-1}^3)  
}
$$
for $(i=1,2,\ldots,N)$
This equation describes the coupled motion of all the particles.
## Initial Conditions

To investigate energy transfer between normal modes, we initially excite only the first normal mode of the chain.
The initial displacement of the (i)-th particle is
$$  
\boxed{  
x_i(0)=A\sin\left(\frac{i\pi}{N+1}\right)  
}  
$$
where (A) is the initial displacement amplitude.
All particles initially have zero velocity:
$$  
\boxed{v_i(0)=0}  
$$
This choice ensures that the initial energy is concentrated in the first normal mode of the corresponding harmonic chain.
The main simulation uses the following initial parameters:

| Parameter                       | Value |
| ------------------------------- | ----- |
| Number of moving particles, (N) | 32    |
| Particle mass, (m)              | 1.0   |
| Linear spring constant, (K)     | 1.0   |
| Initial amplitude, (A)          | 1.0   |
| Cubic coefficient, $(\alpha)$   | 0.25  |
| Quartic coefficient, $(\beta)$  | 0.0   |
| Initial velocity                | 0     |

The endpoints remain fixed throughout the simulation.

---
# 2. Motivation for Computational Methods
## Why Is an Analytical Solution Difficult?
To understand why computational methods are required for the FPUT problem, it is useful to first consider the simpler case of a purely harmonic chain.
### Harmonic Chain: The Linear Case
When both nonlinear coefficients are zero, i.e. ($\alpha$=$\beta$=0), the potential energy of the system is
$$  
U=\frac{K}{2}\sum_{i=0}^{N}(x_{i+1}-x_i)^2  
$$
Using Newton's second law, the equation of motion for the $(i)$-th particle becomes
$$
m\ddot{x}_i=  
K(r_i-r_{i-1}) $$
$$  
r_i=x_{i+1}-x_i  
$$
$$  
m\ddot{x}_i=K(x_{i+1}-2x_i+x_{i-1})  
$$
This is a system of linear, coupled differential equations. In matrix form, it can be written as
$$  
m\ddot{\mathbf{x}}=-K D\mathbf{x},  
$$
where (D) is the tridiagonal matrix

$$
D =
\begin{pmatrix}
2 & -1 & 0 & \cdots & 0 \\
-1 & 2 & -1 & \cdots & 0 \\
0 & -1 & 2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & -1 \\
0 & 0 & 0 & -1 & 2
\end{pmatrix}
$$

Since this matrix is symmetric, it can be diagonalized using its eigenvectors. These eigenvectors represent the normal modes of the chain.

Introducing the normal-mode coordinates ($Q_k$), the equations become
$$  
\ddot{Q}_k+\omega_k^2Q_k=0  
$$
where the normal-mode frequencies are

$$ 
\omega_k=2\sqrt{\frac{K}{m}}  
\sin\left(\frac{k\pi}{2(N+1)}\right),  
\qquad k=1,2,\ldots,N.  
$$

Each equation now describes an independent simple harmonic oscillator. Its general solution is
$$
Q_k(t)=A_k\cos(\omega_k t)+B_k\sin(\omega_k t),  
$$
where $(A_k)$ and $(B_k)$ are determined by the initial conditions.

Thus, the harmonic chain can be solved analytically by transforming the particle coordinates into normal-mode coordinates. The motion of every particle can then be reconstructed from these solutions.
### Nonlinear Chain

The situation changes when nonlinear terms are introduced into the spring potential:
$$
V(r_i)=\frac{K}{2}r_i^2+  
\frac{\alpha}{3}r_i^3+  
\frac{\beta}{4}r_i^4.  
$$
The resulting equation of motion is
$$
m\ddot{x}_i=  
K(r_i-r_{i-1})  
+\alpha(r_i^2-r_{i-1}^2)  
+\beta(r_i^3-r_{i-1}^3).  
$$

Unlike the harmonic case, this equation contains quadratic and cubic terms in the spring extensions.

For the FPUT-alpha model, the quadratic terms are present, while for the FPUT-beta model, the cubic terms are present. When both coefficients are nonzero, both kinds of interactions affect the motion.

The main challenge becomes clear when we transform the equations into normal-mode coordinates. In this case, the transformation leads to independent equations for each mode. In the nonlinear case, products of normal-mode coordinates show up in the equations.

As a result the equation for one mode depends on the amplitudes of modes. The modes are no longer independent and energy can move from one mode to another.
So even though the nonlinear equations can be written down mathematically, finding an analytical expression, for how all particles and modes evolve over time is much more difficult. Analytical approximations might work for nonlinearities or certain initial conditions but they usually do not capture the full long-term behaviour of the system.
##  Why Is Numerical Simulation Necessary?
Numerical simulation allows us to calculate the positions and velocities of all particles at successive time steps. Starting from the initial conditions, we calculate the forces acting on each particle, determine its acceleration, and update its position and velocity.
By repeating this process, we can follow the evolution of the entire chain over a long period.
Computational methods allow us to investigate the time evolution of particle displacements, conservation of total mechanical energy, transfer of energy between normal modes, effect of nonlinear interactions on the dynamics, possibility of energy recurrence, influence of changing the nonlinear coefficients.
These quantities are difficult to investigate using analytical methods alone, particularly over long time intervals.

---
# 3. Computational Method and Implementation

## a. Overview of the Computational Approach
We have made this project using Python, primarily with the NumPy and Matplotlib libraries.
NumPy is used for numerical calculations involving particle positions, velocities, forces, and normal-mode transformations. Matplotlib is used to visualize the simulation results.
The computational procedure consists of the following stages:
1. Define the physical and numerical parameters.
2. Generate the initial positions and velocities.
3. Calculate the spring extensions.
4. Calculate the potential energy, forces, and accelerations.
5. Numerically integrate the equations of motion.
6. Calculate the total energy at selected time intervals.
7. Transform the particle coordinates into normal-mode coordinates.
8. Calculate modal energies and spectral entropy.
9. Compare the results for different nonlinear models.
## b. Calculation of Spring Extensions
The first step is to calculate the extension of every spring.
For the fixed-end chain, we introduce an array containing the two fixed endpoints along with the moving particles:
$$  
\mathbf{x}_{\text{full}}=  
[0,x_1,x_2,\ldots,x_N,0]  
$$
The spring extensions are calculated using
$$  
r_i=x_{i+1}-x_i.  
$$
In Python, the implementation is:
```
def calculate_spring_extensions(x):
    x_full = np.zeros(N + 2)
    x_full[1:N+1] = x
    extensions = np.diff(x_full)
    return extensions
```
The function returns an array containing the extensions of all (N+1) springs.
Including the fixed endpoints ensures that the boundary conditions are maintained correctly.
## c. Calculation of Potential Energy
The total potential energy is calculated using
$$ 
U=  
\sum_{i=0}^{N}  
\left[  
\frac{K}{2}r_i^2+  
\frac{\alpha}{3}r_i^3+  
\frac{\beta}{4}r_i^4  
\right].  
$$
The three contributions are calculated separately:
$$ 
U_{\text{harmonic}}=  
\frac{K}{2}\sum_{i=0}^{N}r_i^2  
$$
$$ 
U_{\alpha}=  
\frac{\alpha}{3}\sum_{i=0}^{N}r_i^3  
$$
$$  
U_{\beta}=  
\frac{\beta}{4}\sum_{i=0}^{N}r_i^4  
$$
Therefore,
$$
\boxed{U=U_{\text{harmonic}}+U_{\alpha}+U_{\beta}}  
$$

The implementation is:

```
def calculate_potential_energy(x, alpha, beta):
    extensions = calculate_spring_extensions(x)
    harmonic_energy = 0.5 * K * np.sum(extensions**2)
    alpha_energy = (alpha / 3.0) * np.sum(extensions**3)
    beta_energy = (beta / 4.0) * np.sum(extensions**4)
    total_potential_energy = (
        harmonic_energy
        + alpha_energy
        + beta_energy
    )
    return total_potential_energy
```
## d. Calculation of Forces and Accelerations
The force acting on each particle is obtained by differentiating the potential energy with respect to its displacement:
$$
F_i=-\frac{\partial U}{\partial x_i}$$
The resulting expression is
$$
F_i=  
K(r_i-r_{i-1})  
+\alpha(r_i^2-r_{i-1}^2)  
+\beta(r_i^3-r_{i-1}^3).  
$$
The force is calculated by considering the contributions from the two neighbouring springs.
The acceleration is then obtained using Newton's second law:
$$
\boxed{a_i=\frac{F_i}{m}}.  
$$
The corresponding Python functions are:
```
def calculate_force(x, alpha, beta):
    extensions = calculate_spring_extensions(x)
    left_extensions = extensions[:-1]
    right_extensions = extensions[1:]
    harmonic_force = K * (right_extensions - left_extensions)
    alpha_force = alpha * (right_extensions**2 - left_extensions**2)
    beta_force = beta * (right_extensions**3 - left_extensions**3)
    force = harmonic_force + alpha_force + beta_force
    return force
```

```
def calculate_acceleration(x, alpha, beta):
    force = calculate_force(x, alpha, beta)
    acceleration = force / m
    return acceleration
```

The force calculation is also checked numerically by comparing the analytical force with the negative central finite difference of the potential energy:
$$  
F_i\approx  
-\frac{U(x_i+\epsilon)-U(x_i-\epsilon)}  
{2\epsilon}.  
$$
## e. Numerical Integration: Velocity Verlet Method
The equations of motion are second-order differential equations:
$$ 
\ddot{x}_i=a_i(\mathbf{x})  
$$
To solve them numerically, we use the **Velocity Verlet integration method**.
	The method is particularly suitable for mechanical systems because it is a second-order, time-reversible, symplectic integration method. A symplectic method is  useful for Hamiltonian systems because it generally provides good long-term energy behaviour.
Let the time step be $(\Delta t)$. At time step (n), the position, velocity, and acceleration are represented by $(x_n)$, $(v_n)$, and $(a_n)$.
The algorithm consists of three main operations.
**Step 1: First half-step velocity update**
$$
v_{n+\frac12} = v_n+\frac{\Delta t}{2}a_n  
$$
The velocity is updated by half a time step using the current acceleration.
**Step 2: Position update**
$$
\boxed{  
x_{n+1}=x_n+\Delta t*v_{n+\frac12}  
}  
$$
The new position is calculated using the intermediate velocity.
**Step 3: Calculate new acceleration and complete velocity update**
First, the acceleration at the new position is calculated:
$$ 
a_{n+1}=a(x_{n+1}).  
$$
Then the velocity is updated:
$$ v_{n+1}=v_{n+\frac12} +\frac{\Delta t}{2}a_{n+1}
$$
$$
\boxed{v_{n+1}=v_n+ \frac{\Delta t}{2} \left[a(x_n)+a(x_{n+1})\right]}
$$
These operations are repeated until the desired simulation time is reached.

The corresponding Python implementation is:
```
for step in range(1, num_steps + 1):
    # First half-step velocity update
    v_half = v_current + 0.5 * dt * a_current
    
    # Position update
    x_new = x_current + dt * v_half
    
    # Acceleration at the new position
    a_new = calculate_acceleration(x_new, alpha, beta)
    
    # Complete velocity update
    v_new = v_half + 0.5 * dt * a_new
    
    # Update current values
    x_current = x_new
    v_current = v_new
    a_current = a_new
```

The simulation saves the particle positions, velocities, and total energy after every 400 integration steps, as well as at the final time.
The main numerical parameters are:

| Parameter                 | Value     |
| ------------------------- | --------- |
| Time step, $(\Delta t)$   | 0.05      |
| Total simulation time     | 40000     |
| Total integration steps   | 800000    |
| Saving interval           | 400 steps |
| Time between saved states | 20        |

The small time step allows the continuous equations of motion to be approximated by successive numerical updates.

## f. Calculation of Total Energy

The kinetic energy is
$$  
T=\frac12m\sum_{i=1}^{N}v_i^2  
$$
The total potential energy is calculated from the spring extensions.
Therefore, the total Hamiltonian is
$$
\boxed{H=T+U}
$$
Since the system is isolated and has no damping or external driving force, its total mechanical energy should remain constant.
To examine the numerical accuracy, we calculate the relative energy error:
$$ 
\boxed{  
\epsilon_H(t)=  
\frac{H(t)-H(0)}{H(0)}  
}  
$$
A small relative energy error indicates that the numerical integration is maintaining energy conservation reasonably well.
The energy error is also useful for identifying whether the chosen time step is sufficiently small.

## g. Normal-Mode Analysis

Normal modes are collective patterns of motion in which all particles oscillate with the same frequency, with fixed relative amplitudes.
For the nonlinear chain, these are harmonic-basis modal energies. They are useful for diagnosing energy transfer, but individual modal energies are not conserved, and their sum is not necessarily equal to the full nonlinear Hamiltonian.

For a linear chain of $N$ particles with fixed ends, the equation of motion is $$ m\ddot{x}_i=K(x_{i+1}-2x_i+x_{i-1}), \qquad x_0=x_{N+1}=0.$$

Assume harmonic motion of the form $\qquad\qquad x_i(t)=\phi_i\cos(\omega t).$

Substituting into the equation of motion gives $\qquad\qquad 2\phi_i-\phi_{i+1}-\phi_{i-1}=\lambda\phi_i, \qquad \lambda=\frac{m\omega^2}{K}.$

Assume a sine-shaped mode: $\qquad\qquad \phi_i=\sin(i\theta).$

Using the trigonometric identity $$\qquad\qquad \sin((i+1)\theta)+\sin((i-1)\theta)=2\sin(i\theta)\cos\theta,$ we obtain $\qquad\qquad \lambda=2(1-\cos\theta)=4\sin^2\left(\frac{\theta}{2}\right).$$

The fixed-end boundary condition $(\phi_{N+1}=0)$ requires $\qquad\qquad \sin((N+1)\theta)=0.$

Therefore, $\qquad\qquad \theta=\frac{k\pi}{N+1}, \qquad k=1,2,\ldots,N.$

Hence, the normal-mode shape is $\qquad\qquad \boxed{\phi_k(i)=\sin\left(\frac{ik\pi}{N+1}\right)}.$

To normalize the mode shapes, we use $\qquad\qquad \sum_{i=1}^{N}\sin^2\left(\frac{ik\pi}{N+1}\right)=\frac{N+1}{2}.$

This gives the normalized mode matrix: $$\boxed{\Phi_{ik}=\sqrt{\frac{2}{N+1}}\sin\left(\frac{ik\pi}{N+1}\right)}.$$

The columns of this matrix are orthonormal, so $\qquad\qquad \Phi^T\Phi=I.$

We can transform the particle displacement vector ($x$) into normal coordinates ($Q$) using $\qquad\qquad \boxed{Q=\Phi^Tx}.$

The inverse transformation is $\qquad\qquad \boxed{x=\Phi Q}.$

Similarly, transforming the particle velocity vector gives $\qquad\qquad \boxed{\dot Q=\Phi^T\dot x}.$

Here, $(Q)_k$ and $(\dot Q_k)$ represent the displacement amplitude and velocity of the $(k)$-th normal mode.
The corresponding normal-mode frequencies are
$$
\boxed{
\omega_k=
2\sqrt{\frac Km}
\sin\left(\frac{k\pi}{2(N+1)}\right)
}.
$$
For the simulation, where (m=K=1),
$$
\boxed{
\omega_k=
2\sin\left(\frac{k\pi}{2(N+1)}\right)
}.
$$


In our implementation, the position and velocity arrays are stored with time as the first dimension. Therefore, the transformation is performed using matrix multiplication:

```
Q = position_values @ Phi
Q_dot = velocity_values @ Phi
```

The harmonic energy associated with mode (k) is
$$ 
\boxed{  
E_k=  
\frac12m\dot Q_k^2+  
\frac12m\omega_k^2Q_k^2  
}  
$$
The fraction of the total modal energy contained in mode (k) is
$$
\boxed{  
f_k=\frac{E_k}{\sum_{j=1}^{N}E_j}  
}  
$$
These fractions allow us to study how energy initially concentrated in the first mode is distributed among the other modes.
The normal coordinates allow us to study how the system's energy is distributed among the different modes. In the nonlinear FPUT chain, these linear modes are used as a reference basis to observe energy transfer between modes.

---
# 4. Results and Discussion
## Initial Configuration and Normal Modes

The initial displacement is chosen according to the first normal-mode shape:
$$
x_i(0)=A\sin\left(\frac{i\pi}{N+1}\right).  
$$
The initial configuration therefore has a sinusoidal shape, with the displacement becoming zero at both fixed endpoints.
The first six normal-mode shapes are also plotted to illustrate the possible patterns of oscillation of the chain.
 ##### **Figure 1: Initial configuration of the FPUT chain**
![[Pasted image 20261101123833.png|418]]

**Figure 2: First six normal modes of the fixed-end chain**
![[Pasted image 20261101124037.png|419]]
## Energy Conservation

The total energy of the system is calculated throughout the simulation.
The relative energy error is plotted as
$$
\epsilon_H(t)=\frac{H(t)-H(0)}{H(0)}.  
$$
 **Figure 3: Relative energy error as a function of time**
![[Pasted image 20261101124357.png|421]]

The energy error should remain small compared with the total initial energy if the numerical integration is sufficiently accurate.
The plot can be used to examine whether the Velocity Verlet method maintains the energy of the system over the long simulation period.

The actual maximum absolute relative energy error obtained from the simulation is:
$$
\boxed{\epsilon_{H,\max}={[3.5138596188914705e-05]}}  
$$
This provides a quantitative measure of the numerical accuracy of the simulation.
## ==🔵Energy Transfer Among Normal Modes==

The normal-mode energy fractions are calculated using the harmonic normal-mode basis.
Initially, the first mode contains nearly all the modal energy because the initial displacement is chosen to match its mode shape and all initial velocities are zero.
As the nonlinear system evolves, interactions between the particles can transfer energy between different normal modes.

**Figure 4: Energy fractions of the first five normal modes as a function of time**
![[Pasted image 20261101124820.png|572]]
The evolution of the modal energy fractions provides information about the extent of energy transfer.
In the harmonic chain, the normal modes are independent oscillators. Therefore, energy initially placed in one mode remains in that mode.

**In contrast, nonlinear interactions introduce coupling between the modes. Energy can consequently move from the initially excited mode to other modes.**

The plotted modal fractions allow this difference to be examined directly.
## Comparison Between Harmonic and FPUT-alpha Chains

To investigate the effect of nonlinear interactions, we compare two simulations:
1. Harmonic chain: $(\alpha=0,\ \beta=0)$
2. FPUT-alpha chain: $(\alpha=0.25,\ \beta=0)$
Both simulations use the same initial conditions and numerical parameters.

**Figure 5: Energy fraction in the first mode for the harmonic and FPUT-alpha chains.**
![[Pasted image 20261101124952.png|568]]
For the harmonic chain, the absence of nonlinear interactions means that the normal modes evolve independently. Therefore, the energy fraction in the first mode should remain close to its initial value.
For the FPUT-alpha chain, the cubic nonlinear interaction couples the modes. The first-mode energy fraction can therefore change over time as energy is transferred to other modes.
This comparison demonstrates the role of nonlinear interactions in the dynamics of the oscillator chain.

---
# 5. Extension and Parameter Study

## Comparison of FPUT-alpha, FPUT-beta, and Combined Models

The main extension in this project is the comparison of the FPUT-alpha, FPUT-beta, and combined nonlinear models.
The purpose is to investigate how the form of the nonlinear spring potential affects the energy transfer between normal modes.
The FPUT-alpha model contains a cubic contribution:
$$
V_\alpha(r)=\frac K2r^2+\frac{\alpha}{3}r^3.  
$$
The FPUT-beta model contains a quartic contribution:
$$
V_\beta(r)=\frac K2r^2+\frac{\beta}{4}r^4.  
$$
The combined model contains both:
$$
V_{\alpha\beta}(r)=  
\frac K2r^2+  
\frac{\alpha}{3}r^3+  
\frac{\beta}{4}r^4.  
$$
The parameters used for the comparison are
$$
\alpha=0.25,\qquad \beta=5.  
$$
The comparison is performed by examining the first-mode energy fraction and spectral entropy for each model. This extension demonstrates that the nonlinear dynamics depend not only on the presence of nonlinear interactions but also on the mathematical form and strength of those interactions.

| Model          | $(\alpha)$ | $(\beta)$ |
| -------------- | ---------- | --------- |
| FPUT-alpha     | 0.25       | 0         |
| FPUT-beta      | 0          | 5         |
| Combined model | 0.25       | 5         |

All three simulations use the same initial conditions, time step, and total simulation time.

**Figure 6: Comparison of the first-mode energy fraction for the three nonlinear models.**
![[Pasted image 20261101125126.png|477]]
The comparison allows us to examine how the cubic and quartic nonlinear terms influence the transfer of energy away from the first mode.
The cubic interaction and quartic interaction have different mathematical forms, so they can produce different dynamical behaviour.
The combined model includes both contributions and allows their effects to be investigated together.
The exact differences between the three curves should be discussed using the numerical results obtained from the simulation.
### Limitations of the Parameter Study
The comparison uses one selected value of $(\alpha)$ and one selected value of $(\beta)$. Therefore, it does not establish how the dynamics change over the entire range of possible nonlinear coefficients.
A more extensive parameter study could investigate several values of $(\alpha)$, $(\beta)$, and the initial amplitude (A). The present comparison is intended to demonstrate the influence of the different nonlinear terms rather than provide a complete parameter-space analysis.

## Spectral Entropy
To quantify how widely the energy is distributed among the normal modes, we calculate the spectral entropy:
$$
\boxed{  
S=-\sum_{k=1}^{N}f_k\ln f_k  
}  
$$
where $(f_k)$ is the modal energy fraction.
The entropy provides a measure of the spread of energy across the modes.
If all the energy is concentrated in a single mode, then
$$
f_1=1,\qquad f_{k\ne1}=0  , \qquad S=0
$$

If the energy is equally distributed among all (N) modes,
$$ 
f_k=\frac1N,  
$$
then
$$  
S=\ln N.  
$$
Therefore,
$$  
\boxed{0\leq S\leq\ln N}.  
$$
A low spectral entropy indicates that energy is concentrated in relatively few modes, whereas a higher entropy indicates that energy is spread across more modes.
Spectral entropy is used here as a diagnostic of modal energy distribution. A high value alone does not prove that the system is chaotic or fully thermalized.
### Spectral Entropy of the FPUT-alpha Chain
The spectral entropy is calculated using the modal energy fractions.

**Figure 7: Spectral entropy of the FPUT-alpha chain as a function of time.**
![[Pasted image 20261101125039.png|502]]
The entropy provides an overall measure of how broadly the energy is distributed among the normal modes.
When most of the energy remains in the first mode, the entropy is relatively low.
When energy spreads across a larger number of modes, the entropy increases.
The maximum possible entropy for the system is
$$
S_{\max}=\ln N=\ln 32.  
$$
This value corresponds to equal modal energy fractions.
The evolution of spectral entropy therefore provides a compact way of studying the changing distribution of energy in the nonlinear chain.
### Comparison of Spectral Entropy for Different Nonlinear Models
We also compare the spectral entropy of the FPUT-alpha, FPUT-beta, and combined models.

**Figure 8: Spectral entropy comparison of the three nonlinear models.**
![[Pasted image 20261101125236.png|510]]
This plot provides an overall comparison of the modal energy distribution. A model with a higher spectral entropy at a particular time has its energy distributed across more modes at that time.
Comparing the entropy curves helps identify differences in how the three nonlinear models redistribute energy. The entropy should be interpreted together with the modal energy fraction plots because it does not identify which particular modes contain the energy.

---
# Summary of Results
The simulations allow us to investigate the following physical features:
- The initial energy is concentrated in the first normal mode.
- Nonlinear interactions introduce coupling between the normal modes.
- The modal energy fractions describe how energy is redistributed over time.
- Spectral entropy provides a measure of the spread of energy among the modes.
- Comparing the alpha, beta, and combined models demonstrates the effect of different nonlinear interactions.
---
# 6. Strengths and Limitations of the Computational Method
## Strengths
The Velocity Verlet method has several advantages for this problem.
1. **Second-order accuracy:** The method has a global error of order $(O(\Delta t^2))$, making it more accurate than the basic Euler method for a comparable time step.
2. **Symplectic nature:** It is a symplectic integrator, which is particularly useful for Hamiltonian systems because it generally provides good long-term energy behaviour.
3. **Time reversibility:** The method has a time-reversible update structure, which is desirable for conservative mechanical systems.
4. **Computational efficiency:** The method requires force and acceleration calculations at each time step and is straightforward to implement using NumPy.
5. **Physical interpretation:** The numerical results can be transformed into normal-mode coordinates, allowing the transfer of energy between different modes to be studied.
## Limitations
Despite its advantages, the method has some limitations.
1. **Finite time-step error:** The numerical solution is an approximation to the exact continuous dynamics. A time step that is too large can produce inaccurate results.
2. **Long simulation time:** The investigation of long-term energy transfer requires a large number of integration steps, increasing computational cost.
3. **Numerical energy error:** Velocity Verlet does not conserve the exact Hamiltonian at every step. The energy error must therefore be monitored.
4. **Limited parameter range:** Only selected values of the nonlinear coefficients are investigated in the present comparison.
5. **Harmonic modal basis:** The modal energies are calculated using the normal modes of the harmonic chain. In the nonlinear system, these modal energies are diagnostic quantities and are not individually conserved.
6. **Interpretation of entropy:** Spectral entropy measures the spread of energy among modes, but it does not by itself establish chaos, or thermal equilibrium.
---
# 7. Conclusion
In this project, we numerically investigated the Fermi-Pasta-Ulam-Tsingou problem using a one-dimensional chain of coupled oscillators with fixed endpoints. The system was modelled using a nonlinear spring potential containing harmonic, cubic, and quartic contributions. Since the resulting equations of motion are coupled and nonlinear, numerical integration was used to study their time evolution. The Velocity Verlet method was implemented to calculate the particle positions and velocities at successive time steps. The total Hamiltonian was monitored to examine the numerical conservation of energy. Normal-mode analysis was then used to investigate the distribution of energy among the modes. By initially exciting only the first normal mode, we could study how nonlinear interactions influence the transfer of energy to other modes. The spectral entropy provided an additional measure of how widely the energy was distributed across the normal modes. Finally, the harmonic, FPUT-alpha, FPUT-beta, and combined nonlinear models were compared to investigate the influence of different nonlinear interactions. Overall, the project demonstrates how computational techniques can be used to investigate the dynamics of a nonlinear physical system that is difficult to solve analytically. It also illustrates the importance of numerical accuracy, energy conservation, and suitable physical diagnostics when studying long-time behaviour in Hamiltonian systems.

---
# 8. References
1. E. Fermi, J. Pasta, S. Ulam, and M. Tsingou, _Studies of Nonlinear Problems_, Los Alamos Scientific Laboratory Report LA-1940 (1955).
2. S. Paleari and T. Penati, _Fermi-Pasta-Ulam Problem_, 2008.
3. Computational Physics project implementation using Python, NumPy, and Matplotlib.



