"""
Plotting utilities for the MD simulation results.

All functions receive plain data (arrays, dicts) — no dependency on the
simple_md C++ module.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_positions(positions, L, sigma, title="Particle positions", save_path=None):
    positions = np.asarray(positions) * sigma
    L_cgs = L * sigma
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(positions[:, 0], positions[:, 1], s=60, edgecolors="k", linewidths=0.5)
    ax.set_xlim(0, L_cgs)
    ax.set_ylim(0, L_cgs)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ [cm]")
    ax.set_ylabel(r"$y$ [cm]")
    ax.set_title(title)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def plot_pressure_vs_density(results, save_path=None):
    rho = [r["density_cgs"] for r in results]
    p = [r["pressure_cgs"] for r in results]

    fig, ax = plt.subplots()
    ax.plot(rho, p, "o-")
    ax.set_xlabel(r"$\rho$ [g / cm$^2$]")
    ax.set_ylabel(r"$p$ [dyn / cm]")
    ax.set_title(r"$p(\rho,\; T=0)$")
    ax.grid(which="both")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def plot_pressure_vs_temperature(results, save_path=None):
    T = [r["temperature_K"] for r in results]
    p = [r["pressure_cgs"] for r in results]

    fig, ax = plt.subplots()
    ax.loglog(T, p, "o-")
    ax.set_xlabel(r"$T$ [K]")
    ax.set_ylabel(r"$p$ [dyn / cm]")
    ax.set_title(r"$p(\rho^*=1.1,\; T)$ — log-log")
    ax.grid(which="both")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def plot_energy_vs_temperature(results, save_path=None):
    T = [r["temperature_K"] for r in results]
    E = [r["energy_cgs"] for r in results]

    fig, ax = plt.subplots()
    ax.semilogx(T, E, "o-")
    ax.set_xlabel(r"$T$ [K]")
    ax.set_ylabel(r"$E/N$ [erg]")
    ax.set_title(r"$E(\rho^*=1.1,\; T)$")
    ax.grid(which="both")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)
