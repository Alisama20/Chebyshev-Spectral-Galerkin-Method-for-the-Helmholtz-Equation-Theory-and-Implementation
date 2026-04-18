"""Solve -Delta u + k^2 u = f on (-1,1)^2 with Dirichlet BCs and plot the
numerical solution against the manufactured analytical one."""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import matplotlib.pyplot as plt
from chebyshev.core import solve_helmholtz, evaluate_solution, u_exact

FIGDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

# ------------------------------
# Parameters
# ------------------------------

N = 12           # basis size
Nq = 40          # quadrature nodes
k = 5.0

# ------------------------------
# Solve
# ------------------------------

u_vec = solve_helmholtz(N, Nq, k)

# Evaluation grid
x_plot = np.linspace(-1, 1, 150)
X, Y = np.meshgrid(x_plot, x_plot)

U_num = evaluate_solution(u_vec, N, X, Y)
U_ex = u_exact(X, Y)
err = np.abs(U_num - U_ex)

print(f"Max absolute error : {err.max():.3e}")
print(f"L2 relative error  : {np.linalg.norm(U_num - U_ex) / np.linalg.norm(U_ex):.3e}")

# ------------------------------
# Plot
# ------------------------------

fig = plt.figure(figsize=(14, 5))
fig.suptitle(
    rf"Chebyshev spectral-Galerkin, $N={N}$, $k={k}$",
    fontsize=14,
)

ax1 = fig.add_subplot(1, 3, 1, projection="3d")
ax1.plot_surface(X, Y, U_num, cmap="viridis")
ax1.set_title("Numerical solution")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

ax2 = fig.add_subplot(1, 3, 2, projection="3d")
ax2.plot_surface(X, Y, U_ex, cmap="plasma")
ax2.set_title("Exact solution")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

ax3 = fig.add_subplot(1, 3, 3)
im = ax3.pcolormesh(X, Y, err, shading="auto", cmap="magma")
ax3.set_title(rf"Absolute error (max = {err.max():.2e})")
ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_aspect("equal")
plt.colorbar(im, ax=ax3, fraction=0.045)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(os.path.join(FIGDIR, "solution.png"), dpi=150)
plt.close()
print("saved figures/solution.png")
