
# Chebyshev Spectral–Galerkin Method for the Helmholtz Equation

This project presents the theoretical derivation and numerical implementation of a Chebyshev spectral–Galerkin method for solving the two-dimensional Helmholtz equation with homogeneous Dirichlet boundary conditions.

The implementation is written in Python and includes matrix assembly, Kronecker-product structure exploitation, and spectral convergence analysis.

---

## Problem Statement

We consider the Helmholtz equation in a bounded domain Ω ⊂ ℝ²:

Δu + k²u = f  in Ω  
u = 0         on ∂Ω  

where:
- Δ is the Laplacian operator
- k is the wave number
- f is a forcing term

---

## Methodology

### 1. Weak Formulation

The problem is formulated in the Sobolev space:

H₀¹(Ω)

The variational formulation reads:

Find u ∈ H₀¹(Ω) such that

∫Ω ∇u · ∇v − k² u v dx = ∫Ω f v dx  
for all test functions v ∈ H₀¹(Ω).

---

### 2. Spectral–Galerkin Discretization

- Basis: Chebyshev polynomials
- Boundary conditions enforced via modified basis functions
- Tensor-product structure for 2D construction
- Global spectral approximation

The discrete system is assembled using Kronecker products, yielding a structured linear system of the form:

A u = b

---

### 3. Numerical Implementation

The implementation includes:

- Construction of Chebyshev differentiation matrices
- Assembly of stiffness and mass matrices
- Kronecker-product formulation for 2D operators
- Direct or iterative linear solver
- Error computation against analytical solution
- Spectral convergence study

---

---

## Results

The method exhibits spectral (exponential) convergence for smooth solutions, as expected from Chebyshev polynomial approximations.

Error decay is analyzed in terms of:
- L² norm
- Maximum norm

Numerical experiments confirm theoretical convergence rates.

---

## Technologies Used

- Python
- NumPy
- SciPy
- Matplotlib

---

## Key Numerical Concepts Demonstrated

- Spectral methods
- Galerkin formulation
- Functional analysis (Sobolev spaces)
- Kronecker-product matrix assembly
- Structured linear systems
- Convergence analysis
- Numerical stability considerations

---


---

## Author

A. S. Amari.

Project developed as part of the evaluation of the subject Numerical Analysis of Pdes and Approximation in Master's Degree in Physics and Mathematics - Fisymat, University of Granada.



