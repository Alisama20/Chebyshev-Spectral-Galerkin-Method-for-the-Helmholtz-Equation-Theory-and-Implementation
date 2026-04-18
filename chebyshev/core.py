"""Chebyshev spectral-Galerkin solver for the 2D Helmholtz equation.

Solves
    -Delta u + k^2 u = f          in  Omega = (-1,1)^2,
                  u = 0           on  dOmega,

with the separable boundary-adapted Chebyshev basis

    phi_n(x) = T_{n+2}(x) - T_n(x),     n = 0, 1, ..., N-1,

which vanishes at x = +/- 1.  Integrals are evaluated with Gauss-Legendre
quadrature so that the unweighted L^2 inner products of the classical
weak formulation are represented exactly.

Matrix assembly uses the tensor structure of the basis: the 2D stiffness
and mass matrices are Kronecker products of the 1D ones,

    A = K x M + M x K + k^2 (M x M).
"""

import numpy as np
import scipy.linalg as la


# ----------------------------------------------------------------------
# 1D building blocks
# ----------------------------------------------------------------------

def chebyshev_T(n, x):
    """Chebyshev polynomial of the first kind, T_n(x), for x in [-1, 1]."""
    return np.cos(n * np.arccos(np.clip(x, -1.0, 1.0)))


def chebyshev_T_derivative(n, x):
    """Derivative T_n'(x) via the identity
         (x^2 - 1) T_n'(x) = n (x T_n(x) - T_{n-1}(x)).

    The removable singularities at x = +/- 1 are patched analytically."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    if n == 0:
        return out
    # Interior values
    mask = np.abs(np.abs(x) - 1.0) > 1e-12
    xi = x[mask]
    out[mask] = n * (xi * chebyshev_T(n, xi) - chebyshev_T(n - 1, xi)) \
                / (xi ** 2 - 1.0)
    # Boundary values: T_n'(+/-1) = (+/-1)^(n-1) * n^2
    out[~mask] = (np.sign(x[~mask]) ** (n - 1)) * (n ** 2)
    return out


def phi(n, x):
    """Basis function phi_n(x) = T_{n+2}(x) - T_n(x) (vanishes at +/- 1)."""
    return chebyshev_T(n + 2, x) - chebyshev_T(n, x)


def phi_derivative(n, x):
    return chebyshev_T_derivative(n + 2, x) - chebyshev_T_derivative(n, x)


# ----------------------------------------------------------------------
# Quadrature and tabulation of the basis on the quadrature nodes
# ----------------------------------------------------------------------

def gauss_legendre_quadrature(Nq):
    """Gauss-Legendre nodes and weights on [-1, 1].

    Integrates polynomials of degree up to 2*Nq - 1 exactly against the
    unweighted Lebesgue measure.
    """
    x, w = np.polynomial.legendre.leggauss(Nq)
    return x, w


def build_phi_tables(N, x):
    """Return Phi[n, q] = phi_n(x[q]) and dPhi[n, q] = phi_n'(x[q])."""
    Phi = np.vstack([phi(n, x) for n in range(N)])
    dPhi = np.vstack([phi_derivative(n, x) for n in range(N)])
    return Phi, dPhi


# ----------------------------------------------------------------------
# 1D and 2D Galerkin matrices
# ----------------------------------------------------------------------

def build_1D_matrices(N, Nq):
    """Assemble the 1D mass and stiffness matrices in the unweighted L^2
    inner product,

        M_{ij} = int phi_i(x) phi_j(x) dx,
        K_{ij} = int phi_i'(x) phi_j'(x) dx,

    using Gauss-Legendre quadrature."""
    x, w = gauss_legendre_quadrature(Nq)
    Phi, dPhi = build_phi_tables(N, x)
    # M_ij = sum_q w_q Phi[i, q] Phi[j, q]
    M = (Phi * w) @ Phi.T
    K = (dPhi * w) @ dPhi.T
    return M, K


def build_2D_matrix(M, K, k):
    """2D Galerkin matrix A = K x M + M x K + k^2 (M x M)."""
    return np.kron(K, M) + np.kron(M, K) + (k ** 2) * np.kron(M, M)


# ----------------------------------------------------------------------
# Problem data: manufactured solution u_exact = (1-x^2)(1-y^2) exp(-x-y)
# ----------------------------------------------------------------------

def u_exact(X, Y):
    return (1 - X ** 2) * (1 - Y ** 2) * np.exp(-X - Y)


def f_source(X, Y, k):
    """Source term matching u_exact via  f = -Delta u + k^2 u.

    With  u = (1-x^2)(1-y^2) e^{-x-y},
          u_xx = (1-y^2) e^{-x-y} (-x^2 + 4x - 1),
          u_yy = (1-x^2) e^{-x-y} (-y^2 + 4y - 1),
    so -u_xx - u_yy + k^2 u = e^{-x-y} * [...].
    """
    t1 = (1 - Y ** 2) * (X ** 2 - 4 * X + 1)
    t2 = (1 - X ** 2) * (Y ** 2 - 4 * Y + 1)
    t3 = (k ** 2) * (1 - X ** 2) * (1 - Y ** 2)
    return np.exp(-X - Y) * (t1 + t2 + t3)


def build_rhs(N, Nq, k):
    """Vectorised assembly of  b_{mn} = int f(x,y) phi_m(x) phi_n(y) dx dy.

    With Gauss-Legendre quadrature,
         b_{mn} = (Phi * wx) F (Phi * wy)^T  at (m, n),
    where F_{ij} = f(x_i, x_j).  Flattened row-major to match the Kronecker
    ordering of build_2D_matrix.
    """
    x, w = gauss_legendre_quadrature(Nq)
    Phi, _ = build_phi_tables(N, x)
    X, Y = np.meshgrid(x, x, indexing="ij")
    F = f_source(X, Y, k)
    # Apply weights on both sides: sum_{i,j} w_i w_j F[i,j] Phi[m,i] Phi[n,j]
    B = (Phi * w) @ F @ (Phi * w).T
    return B.flatten()


# ----------------------------------------------------------------------
# Driver + post-processing
# ----------------------------------------------------------------------

def solve_helmholtz(N, Nq, k):
    """Solve the Helmholtz problem and return the coefficient vector u_{ij}
    (flattened row-major, length N*N)."""
    M, K = build_1D_matrices(N, Nq)
    A = build_2D_matrix(M, K, k)
    b = build_rhs(N, Nq, k)
    return la.solve(A, b)


def evaluate_solution(u_vec, N, X, Y):
    """Evaluate  u_N(x,y) = sum_{i,j} u_{ij} phi_i(x) phi_j(y)  on a grid."""
    U = np.zeros_like(X, dtype=float)
    Phi_X = np.array([phi(i, X) for i in range(N)])   # (N, nx, ny)
    Phi_Y = np.array([phi(j, Y) for j in range(N)])
    for i in range(N):
        for j in range(N):
            U += u_vec[i * N + j] * Phi_X[i] * Phi_Y[j]
    return U


def max_error(N, Nq, k, n_plot=150):
    """Solve with N basis functions and return the L-infinity error against
    the manufactured solution on an (n_plot x n_plot) uniform grid."""
    u_vec = solve_helmholtz(N, Nq, k)
    x = np.linspace(-1, 1, n_plot)
    X, Y = np.meshgrid(x, x)
    U_num = evaluate_solution(u_vec, N, X, Y)
    return float(np.max(np.abs(U_num - u_exact(X, Y))))
