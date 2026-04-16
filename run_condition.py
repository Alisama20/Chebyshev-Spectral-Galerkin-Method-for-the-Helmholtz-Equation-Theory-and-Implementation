"""Condition number of the Galerkin matrix A as a function of N.

Spectral methods typically produce dense, algebraically ill-conditioned
matrices; for Chebyshev + second-order Dirichlet basis one expects a
growth roughly like N^4 for the 2D Laplacian.
"""

import numpy as np
import matplotlib.pyplot as plt
from chebyshev_core import build_1D_matrices, build_2D_matrix

# ------------------------------
# Parameters
# ------------------------------

N_values = np.arange(4, 22, 2)
Nq = 40
k = 5.0

# ------------------------------
# Sweep over N
# ------------------------------

cond_numbers = []
for N in N_values:
    M, K = build_1D_matrices(N, Nq)
    A = build_2D_matrix(M, K, k)
    cA = np.linalg.cond(A)
    cond_numbers.append(cA)
    print(f"N = {N:3d}   cond(A) = {cA:.3e}")

cond_numbers = np.array(cond_numbers)

# ------------------------------
# Plot
# ------------------------------

plt.figure(figsize=(8, 5))
plt.semilogy(N_values, cond_numbers, "o-", linewidth=2, label=r"$\kappa(A)$")
# Reference N^4 slope
ref = cond_numbers[0] * (N_values / N_values[0]) ** 4
plt.semilogy(N_values, ref, "--", alpha=0.6, label=r"$\sim N^4$")
plt.xlabel(r"$N$")
plt.ylabel("Condition number")
plt.title(rf"Conditioning of the Galerkin matrix (k = {k})")
plt.grid(True, which="both", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.savefig("figures/condition.png", dpi=150)
plt.show()
