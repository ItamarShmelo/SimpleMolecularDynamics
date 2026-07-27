#!/usr/bin/env python3
"""
Run Assignments A and B of the MD simulation project.
Generates all required plots and saves them in figs/.
"""

import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "compiled"))
sys.path.insert(0, os.path.join(ROOT, "src", "python"))
sys.path.insert(0, ROOT)
import simple_md

import numpy as np

SEED = 42
np.random.seed(SEED)

from simulation import (
    run_zero_temperature,
    run_finite_temperature,
)
from plots.plotting import (
    plot_positions,
    plot_pressure_vs_density,
    plot_pressure_vs_temperature,
    plot_energy_vs_temperature,
)

# Physical constants (CGS)
SIGMA = 2e-8        # cm
EPSILON = 6.9e-13   # erg
MASS = 8.3e-23      # g
K_B = 1.38e-16      # erg/K

# Simulation parameters
N = 100
DT_STAR = 1e-3
ZERO_T_STEPS = 100_000
RELAX_STEPS = 50_000
AVG_STEPS = 50_000

FIGS = os.path.join(ROOT, "figs")
os.makedirs(FIGS, exist_ok=True)

# ===== Assignment A =====

rho_star_values = np.arange(1.0, 1.55, 0.05)
results_A = []

for rho_star in rho_star_values:
    rho_star = round(float(rho_star), 2)
    result = run_zero_temperature(
        rho_star, N, DT_STAR, ZERO_T_STEPS, SIGMA, EPSILON, MASS,
    )
    results_A.append(result)

print("\n\n" + "=" * 70)
print("  Assignment A - Summary: Density, Energy, Pressure")
print("=" * 70)
print(f"  {'rho [g/cm^2]':>14}  {'E/N [erg]':>14}  {'Pressure [dyn/cm]':>18}")
print("  " + "-" * 52)
for r in results_A:
    print(
        f"  {r['density_cgs']:14.6e}  {r['energy_cgs']:14.6e}  "
        f"{r['pressure_cgs']:18.6e}"
    )

# ---- plots ----
plot_pressure_vs_density(
    results_A,
    save_path=os.path.join(FIGS, "assignment_a_pressure_vs_density.png"),
)

for idx in (0, len(results_A) // 2, len(results_A) - 1):
    r = results_A[idx]
    plot_positions(
        r["positions"],
        r["L"],
        SIGMA,
        title=rf"$\rho^* = {r['rho_star']:.2f}$ (T = 0)",
        save_path=os.path.join(
            FIGS, f"assignment_a_positions_rho{r['rho_star']:.2f}.png"
        ),
    )

# ===== Assignment B =====

result_rho11 = next(r for r in results_A if r["rho_star"] == 1.1)
E_min = result_rho11["energy_star"]
eq_positions_np = result_rho11["positions"]
L_B = result_rho11["L"]

eq_positions = [
    simple_md.Vector2D(float(eq_positions_np[i, 0]), float(eq_positions_np[i, 1]))
    for i in range(N)
]

# 5 extra low-energy points targeting T ~ 100..4000 K
TEMP_UNIT = EPSILON / K_B
target_temps_K = np.array([100, 500, 1000, 2000, 4000])
eps_extra_low = E_min + target_temps_K / TEMP_UNIT

eps_low = np.linspace(E_min, 10, 7)[1:6]
eps_high = np.geomspace(10, 100, 7)[1:6]
epsilon_stars = np.concatenate([eps_extra_low, eps_low, eps_high])

results_B = []
for eps_star in epsilon_stars:
    eps_star = float(eps_star)
    result = run_finite_temperature(
        1.1, eps_star, N, DT_STAR,
        RELAX_STEPS, AVG_STEPS,
        SIGMA, EPSILON, MASS,
        equilibrium_positions=eq_positions, L=L_B,
    )
    results_B.append(result)

print("\n\nAssignment B - Summary: Energy, Pressure, Temperature")
print(f"  {'E/N [erg]':>14}  {'Pressure [dyn/cm]':>18}  {'Temperature [K]':>16}")
print("  " + "-" * 54)
for r in results_B:
    print(
        f"  {r['energy_cgs']:14.6e}  {r['pressure_cgs']:18.6e}  "
        f"{r['temperature_K']:16.1f}"
    )

# ---- plots ----
plot_pressure_vs_temperature(
    results_B,
    save_path=os.path.join(FIGS, "assignment_b_pressure_vs_temperature.png"),
)

plot_energy_vs_temperature(
    results_B,
    save_path=os.path.join(FIGS, "assignment_b_energy_vs_temperature.png"),
)

mid = len(results_B) // 2
snap_result = results_B[mid]
for i, (pos, t) in enumerate(
    zip(snap_result["snapshot_positions"], snap_result["snapshot_times"])
):
    plot_positions(
        pos,
        snap_result["L"],
        SIGMA,
        title=rf"$\varepsilon^* = {snap_result['epsilon_star']:.1f}$, step {t}",
        save_path=os.path.join(FIGS, f"assignment_b_positions_t{i}.png"),
    )

print(f"\nAll plots saved to {FIGS}/")
