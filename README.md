# Finite Difference Projects & CPSO Optimization in MATLAB

## Overview

This repository is dedicated to exploring numerical methods for partial differential equations (PDEs) using **finite difference schemes**, combined with **optimization techniques** such as **Continuous Particle Swarm Optimization (CPSO)**. The project is both a study of classical numerical analysis methods and a demonstration of advanced control optimization techniques in MATLAB.  

The main goals of the repository are:

1. **Develop robust finite difference solvers** for canonical PDEs (heat, Poisson, wave equations).
2. **Analyze convergence, stability, and accuracy** of these solvers.
3. **Integrate control and optimization strategies** to match a desired system trajectory.
4. Provide **clean, reusable MATLAB code** that demonstrates good practices, clarity, and efficiency.

---

## Repository Structure

```bash
├── 📁 Example
│   ├── 📄 example_1.m
│   ├── 📄 example_2.m
│   ├── 📄 example_3.m
│   └── 📄 example_4.m
│
├── 📁 project_oscillator
│   ├── 📄 build_control.m
│   ├── 📄 forward_solver.m
│   ├── 📄 main.m
│   ├── 📄 objective_cpso.m
│   ├── 📄 params.m
│   └── 📄 postprocess.m
│
├── ⚙️ .gitignore
│
└── 📝 README.md
```

---


## 1. Finite Difference Examples

### 1.1 2D Heat Equation

We solve the **2D heat equation**:

$$
u_t = \alpha (u_{xx} + u_{yy})
$$

with **Dirichlet boundary conditions** on a unit square \([0,1]^2\) and an initial Gaussian temperature distribution.  

**Highlights:**
- Explicit finite difference scheme.
- CFL-based timestep stability check.
- Real-time 2D visualization of the temperature field.
- Energy monitoring to check numerical dissipation.

```bash
# Explicit update step
u = u + alpha * dt * del2(u, dx, dy);
# Apply Dirichlet boundary conditions
u(1,:) = 0; u(end,:) = 0;
u(:,1) = 0; u(:,end) = 0;
```

--- 

### 1.2 Poisson Equation with Convergence Study

We solve the **2D Poisson equation** on a square domain \([0,1]^2\) with homogeneous Dirichlet boundary conditions:

$$
(u_{xx} + u_{yy}) = f(x,y), \quad u|_{\partial \Omega} = 0
$$

where we choose a **manufactured solution**:

$$
u_{\text{exact}}(x,y) = \sin(\pi x)\sin(\pi y)
$$

which leads to the corresponding source term:

$$
f(x,y) = 2 \pi^2 \sin(\pi x) \sin(\pi y)
$$

#### Numerical Approach

- **Spatial Discretization:** Uniform grid of \(N \times N\) internal points.
- **Finite Difference Scheme:** Standard 5-point Laplacian stencil.
- **Matrix Formulation:** Sparse matrices are used with Kronecker products to efficiently construct the 2D Laplacian:

```bash
e = ones(N,1);
L1 = spdiags([e -2*e e], -1:1, N, N) / h^2;
I = speye(N);
L2 = kron(I,L1) + kron(L1,I);
```

- Boundary Conditions: Homogeneous Dirichlet enforced by zero-padding around the internal grid.

#### Convergence Analysis

To rigorously evaluate the accuracy of the 2D Poisson solver, we perform a **grid refinement study**. This allows us to verify that the numerical solution converges to the exact solution with the expected order of accuracy.

**Procedure:**

1. Solve the Poisson equation for a sequence of grid resolutions:

$$
N = [10, 20, 40, 80]
$$

where \(N\) represents the number of internal points along each spatial direction.

2. Compute the **grid spacing**:

$$
h = \frac{1}{N+1}
$$

3. Evaluate the **discrete \(L^2\) error** between the numerical solution \(u_{\text{num}}\) and the exact solution \(u_{\text{exact}}\):

$$
\text{err}_{L^2} = \sqrt{ h^2 \sum_{i,j} \left( u_{i,j}^{\text{num}} - u_{i,j}^{\text{exact}} \right)^2 }
$$

4. Estimate the **empirical order of convergence** using successive grid refinements:

$$
\text{rate} = \log_2 \frac{\text{err}_{L^2}(N)}{\text{err}_{L^2}(2N)}
$$

--- 

### 1.3 Wave Equation 1D

The one-dimensional wave equation models the propagation of waves along a string or medium:

$$
u_{tt} = c^2 \, u_{xx}, \quad x \in [0,1], \ t>0
$$

with **Dirichlet boundary conditions**:

$$
u(0,t) = u(1,t) = 0
$$

and **initial conditions**:

$$
u(x,0) = \sin(\pi x), \quad u_t(x,0) = 0
$$

The exact solution for these initial and boundary conditions is:

$$
u(x,t) = \sin(\pi x) \cos(\pi c t)
$$

---

#### Numerical Discretization

We discretize the spatial domain using a uniform grid with \(Nx\) points:

$$
x_i = i \, \Delta x, \quad i = 0, 1, \dots, Nx-1
$$

where \(\Delta x = L_x / (Nx-1)\).  

Time is discretized with a **central difference scheme** in time and space (explicit second-order accurate):

$$
u_i^{n+1} = 2 u_i^n - u_i^{n-1} + \left( \frac{c \, \Delta t}{\Delta x} \right)^2 \left( u_{i+1}^n - 2 u_i^n + u_{i-1}^n \right)
$$

with CFL condition for stability:

$$
\frac{c \, \Delta t}{\Delta x} \le 1
$$

The first time step is initialized using the known initial displacement and velocity:

$$
u_i^1 = u_i^0 + \Delta t \, u_t(x_i,0) + \frac{1}{2} (\Delta t)^2 \, u_{tt}(x_i,0)
$$

---

#### MATLAB Implementation

```bash
Nx = 200; Lx = 1; c = 1;
dx = Lx / (Nx-1);
x = linspace(0,Lx,Nx);
dt = 0.9 * dx / c;   % CFL < 1
Tmax = 2;
Nt = ceil(Tmax/dt);

# Initial conditions
u0 = sin(pi*x);   % u(x,0)
u_prev = u0;      % u^{n-1}
u = u0;           % u^n

# Time-stepping loop
for n = 1:Nt
    t = n*dt;
    u_next = 2*u - u_prev + (c*dt/dx)^2 * ([u(2:end) 0] - 2*u + [0 u(1:end-1)]);
    u_next(1) = 0; u_next(end) = 0; % Boundary conditions
    u_prev = u;
    u = u_next;
end

# Compute final L2 error
u_exact = sin(pi*x) * cos(pi*c*Tmax);
err_L2 = sqrt(dx * sum((u - u_exact).^2));
```

---

---

### 2 CPSO-Controlled Oscillator

In this part of the project, the focus is on controlling a **nonlinear oscillator** using **Continuous Particle Swarm Optimization (CPSO)**. The goal is to reconstruct a **control input** that drives the system toward a predefined **target trajectory**.

---

#### 2.1 Dynamical System

We consider a **first-order dynamical system** of the form:

$$
\dot{u}(t) = \alpha u(t) + \beta g(t), \quad u(0) = u_0
$$

where:

- \(u(t)\) is the state of the system (oscillator amplitude),
- \(g(t)\) is the control input to be optimized,
- \(\alpha, \beta\) are system parameters defining the oscillator dynamics.

The **target trajectory** \(u_{\text{target}}(t)\) is synthetically generated using a known control \(p_{\text{true}}\), which allows evaluating the quality of the optimized control.

---

#### 2.2 Control Parametrization

The control signal \(g(t)\) is **parametrized as a vector** \(p \in \mathbb{R}^{M}\), where \(M\) is the number of control parameters:

- **Piecewise constant**: \(g(t)\) takes constant values over \(M\) intervals.
- **Spline interpolation**: \(g(t)\) is reconstructed as a smooth function passing through the \(M\) control points.

The choice of parametrization allows **reducing the dimensionality** of the optimization problem while maintaining flexibility in the control signal.

---

#### 2.3 Objective Function

The CPSO aims to find \(p\) that **minimizes the difference** between the simulated trajectory \(u(t;p)\) and the target trajectory \(u_{\text{target}}(t)\). The cost function is defined as:

$$
J(p) = \frac{1}{N_t} \sum_{i=1}^{N_t} \big(u_i(p) - u_{\text{target},i}\big)^2
$$

where \(N_t\) is the number of discrete time points.

---

#### 2.4 Forward Solver

Given a candidate control vector \(p\), the **forward solver** integrates the oscillator ODE:

```bash
[t, u] = forward_solver(p, params);
```

- Inputs: the control vector `p` and system parameters `params` (including time span, number of time steps, and system constants)

- Outputs:
    - `t`: time vector
    - `u`: state trajectory of the oscillator

The solver reconstructs the continuous control signal `g(t)` from the discrete vector `p` using the chosen parametrization (piecewise or spline). The oscillator dynamics are then integrated over the full time span using standard numerical integration methods.

---

> [!CAUTION]
> Note:
Il progetto è ancora in fase di sviluppo.

--- 

<!--───────────────────────────────────────────────-->
<!--                   AUTORE                     -->
<!--───────────────────────────────────────────────-->

<h2 align="center">✨ Autore</h2>

<p align="center">
  <strong>Giovanni Giuseppe Iacuzzo</strong><br>
  <em>Studente di Ingegneria Dell'IA e della CyberSecurity · Università degli Studi Kore di Enna</em>
</p>

<p align="center">
  <a href="https://github.com/giovanniIacuzzo" target="_blank">
    <img src="https://img.shields.io/badge/GitHub-%40giovanniIacuzzo-181717?style=for-the-badge&logo=github" alt="GitHub"/>
  </a>
  <a href="mailto:giovanni.iacuzzo@unikorestudent.com">
    <img src="https://img.shields.io/badge/Email-Contattami-blue?style=for-the-badge&logo=gmail" alt="Email"/>
  </a>
</p>