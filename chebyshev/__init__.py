"""Chebyshev spectral-Galerkin package for the 2D Helmholtz equation."""

from .core import (
    chebyshev_T,
    chebyshev_T_derivative,
    phi,
    phi_derivative,
    gauss_legendre_quadrature,
    build_phi_tables,
    build_1D_matrices,
    build_2D_matrix,
    u_exact,
    f_source,
    build_rhs,
    solve_helmholtz,
    evaluate_solution,
    max_error,
)

__all__ = [
    "chebyshev_T",
    "chebyshev_T_derivative",
    "phi",
    "phi_derivative",
    "gauss_legendre_quadrature",
    "build_phi_tables",
    "build_1D_matrices",
    "build_2D_matrix",
    "u_exact",
    "f_source",
    "build_rhs",
    "solve_helmholtz",
    "evaluate_solution",
    "max_error",
]
