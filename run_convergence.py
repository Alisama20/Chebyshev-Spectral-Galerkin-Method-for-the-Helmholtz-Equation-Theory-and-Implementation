"""Spectral convergence study.

For a smooth manufactured solution the Chebyshev spectral-Galerkin error
is expected to decay exponentially (faster than any algebraic rate) with
the number of basis functions N.
"""

import numpy as np
import matplotlib.pyplot as plt
from chebyshev_core import max_error

# ------------------------------
# Parameters
# ------------------------------

N_values = np.arange(4, 22, 2)
Nq = 40
k = 5.0

# ------------------------------
# Sweep over N
# ------------------------------

errors = []
for N in N_values:
    e = max_error(N, Nq, k)
    errors.append(e)
    print(f"N = {N:3d}  max-error = {e:.3e}")

errors = np.array(errors)

# ------------------------------
# Plot
# ------------------------------

plt.figure(figsize=(8, 5))
plt.semilogy(N_values, errors, "o-", linewidth=2, label=r"$\|u_N-u\|_{\infty}$")
plt.xlabel(r"$N$")
plt.ylabel("Max error")
plt.title(rf"Spectral convergence (k = {k})")
plt.grid(True, which="both", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.savefig("figures/convergence.png", dpi=150)
plt.show()
