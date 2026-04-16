"""Influence of the wavenumber k on the accuracy and on the conditioning.

The reaction term k^2 u reinforces coercivity, so for this smooth
non-oscillatory manufactured solution both the error and the condition
number decrease as k grows (the matrix becomes increasingly
mass-matrix-dominated).
"""

import numpy as np
import matplotlib.pyplot as plt
from chebyshev_core import (
    build_1D_matrices,
    build_2D_matrix,
    max_error,
)

# ------------------------------
# Parameters
# ------------------------------

k_values = np.linspace(0.5, 10.0, 20)
N = 12
Nq = 40

# ------------------------------
# Sweep
# ------------------------------

errors = []
conds = []
for k in k_values:
    errors.append(max_error(N, Nq, k))
    M, K = build_1D_matrices(N, Nq)
    A = build_2D_matrix(M, K, k)
    conds.append(np.linalg.cond(A))
    print(f"k = {k:5.2f}  error = {errors[-1]:.3e}   cond = {conds[-1]:.3e}")

errors = np.array(errors)
conds = np.array(conds)

# ------------------------------
# Plot
# ------------------------------

fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
fig.suptitle(rf"Influence of $k$  ($N = {N}$)", fontsize=14)

axs[0].semilogy(k_values, errors, "o-", linewidth=2)
axs[0].set_xlabel(r"$k$")
axs[0].set_ylabel("Max error")
axs[0].set_title("Accuracy")
axs[0].grid(True, which="both", alpha=0.4)

axs[1].semilogy(k_values, conds, "s-", linewidth=2, color="tab:orange")
axs[1].set_xlabel(r"$k$")
axs[1].set_ylabel(r"$\kappa(A)$")
axs[1].set_title("Conditioning")
axs[1].grid(True, which="both", alpha=0.4)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("figures/k_influence.png", dpi=150)
plt.show()
