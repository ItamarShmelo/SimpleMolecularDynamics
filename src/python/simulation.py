import math
import numpy as np
import simple_md

K_B = 1.38e-16  # Boltzmann constant [erg / K]


def energy_to_cgs(energy_star, epsilon):
    """Dimensionless energy per particle -> erg."""
    return energy_star * epsilon


def pressure_to_cgs(pressure_star, epsilon, sigma):
    """Dimensionless 2-D pressure -> cgs (dyn / cm)."""
    return pressure_star * epsilon / (sigma ** 2)


def temperature_from_kinetic(kinetic_energy_star, N, epsilon):
    """Dimensionless total kinetic energy -> temperature [K] via equipartition (D=2)."""
    return kinetic_energy_star * epsilon / (N * K_B)


def density_to_cgs(rho_star, mass, sigma):
    """Dimensionless density rho* -> g / cm^2."""
    return rho_star * mass / (sigma ** 2)


def random_positions(N, L):
    xs = np.random.uniform(0, L, N)
    ys = np.random.uniform(0, L, N)
    return [simple_md.Vector2D(float(xs[i]), float(ys[i])) for i in range(N)]


def zero_velocities(N):
    return [simple_md.Vector2D(0.0, 0.0) for _ in range(N)]


def positions_to_array(vec2d_list):
    return np.array([[v.x, v.y] for v in vec2d_list])


def single_zero_T_attempt(N, L, dt_star, num_steps, tol, check_interval):
    """Single zero-T minimization attempt. Returns (sim, energies, converged)."""
    positions = random_positions(N, L)
    velocities = zero_velocities(N)
    sim = simple_md.MDSimulation(positions, velocities, L, dt_star)

    energies = []
    prev_pos = positions_to_array(sim.get_positions())
    converged_count = 0

    for step in range(num_steps):
        sim.step()
        sim.zero_velocities()
        if step % 100 == 0:
            energies.append(sim.get_potential_energy() / N)
        if step % check_interval == 0 and step > 0:
            curr_pos = positions_to_array(sim.get_positions())
            max_disp = np.max(np.linalg.norm(curr_pos - prev_pos, axis=1))
            prev_pos = curr_pos
            if max_disp < tol:
                converged_count += 1
                if converged_count >= 2:
                    break
            else:
                converged_count = 0

    return sim, energies


def run_zero_temperature(rho_star, N, dt_star, num_steps, sigma, epsilon, mass,
                         tol=1e-8, check_interval=1000, n_restarts=3):
    L = math.sqrt(N / rho_star)

    best_sim = None
    best_energy = math.inf
    best_energies = []

    for _ in range(n_restarts):
        sim, energies = single_zero_T_attempt(N, L, dt_star, num_steps, tol, check_interval)
        e = sim.get_potential_energy() / N
        if e < best_energy:
            best_energy = e
            best_sim = sim
            best_energies = energies

    sim = best_sim

    # Average pressure over a short window for stability
    pressures = []
    for _ in range(200):
        sim.step()
        sim.zero_velocities()
        pressures.append(sim.get_pressure())

    energy_star = sim.get_potential_energy() / N
    pressure_star = float(np.mean(pressures))

    return {
        "rho_star": rho_star,
        "N": N,
        "L": L,
        "energy_star": energy_star,
        "pressure_star": pressure_star,
        "energy_cgs": energy_to_cgs(energy_star, epsilon),
        "pressure_cgs": pressure_to_cgs(pressure_star, epsilon, sigma),
        "density_cgs": density_to_cgs(rho_star, mass, sigma),
        "positions": positions_to_array(sim.get_positions()),
        "energies": best_energies,
    }


def run_finite_temperature(rho_star, epsilon_star, N, dt_star,
                           relax_steps, avg_steps,
                           sigma, epsilon, mass,
                           equilibrium_positions=None, L=None):
    if L is None:
        L = math.sqrt(N / rho_star)

    if equilibrium_positions is None:
        pos0 = random_positions(N, L)
        vel0 = zero_velocities(N)
        sim0 = simple_md.MDSimulation(pos0, vel0, L, dt_star)
        for _ in range(50000):
            sim0.step()
            sim0.zero_velocities()
        equilibrium_positions = sim0.get_positions()

    temp_sim = simple_md.MDSimulation(
        equilibrium_positions, zero_velocities(N), L, dt_star,
    )
    potential_energy = temp_sim.get_potential_energy()
    total_energy_target = epsilon_star * N
    kinetic_energy_target = total_energy_target - potential_energy

    if kinetic_energy_target <= 0:
        raise ValueError(
            f"Cannot achieve epsilon*={epsilon_star:.4f}; minimum energy per "
            f"particle is {potential_energy / N:.4f}"
        )

    # --- generate scaled random velocities ------------------------------------
    vx = np.random.randn(N)
    vy = np.random.randn(N)
    vx -= vx.mean()
    vy -= vy.mean()

    current_ke = 0.5 * (np.sum(vx ** 2) + np.sum(vy ** 2))
    scale = math.sqrt(kinetic_energy_target / current_ke)
    vx *= scale
    vy *= scale

    vels = [simple_md.Vector2D(float(vx[i]), float(vy[i])) for i in range(N)]

    sim = simple_md.MDSimulation(equilibrium_positions, vels, L, dt_star)

    therm_window = 1000
    therm_rtol = 0.02       
    prev_mean_T = None
    stable_count = 0

    for step in range(relax_steps):
        sim.step()
        if (step + 1) % therm_window == 0:
            window_temps = []
            for _ in range(therm_window // 10):
                sim.step()
                window_temps.append(sim.get_temperature())
            window_mean = float(np.mean(window_temps))
            if prev_mean_T is not None and prev_mean_T > 0:
                rel_change = abs(window_mean - prev_mean_T) / prev_mean_T
                if rel_change < therm_rtol:
                    stable_count += 1
                    if stable_count >= 3:
                        break
                else:
                    stable_count = 0
            prev_mean_T = window_mean

    # --- averaging ------------------------------------------------------------
    pressures = []
    temperatures_K = []
    total_energies = []
    snapshot_positions = []
    snapshot_times = []
    snapshot_interval = max(1, avg_steps // 5)

    for step in range(avg_steps):
        sim.step()

        temp_star = sim.get_temperature()
        pressures.append(sim.get_pressure())
        temperatures_K.append(temp_star * epsilon / K_B)
        total_energies.append(sim.get_total_energy() / N)

        if step % snapshot_interval == 0:
            snapshot_positions.append(positions_to_array(sim.get_positions()))
            snapshot_times.append(sim.step_count)

    mean_pressure_star = float(np.mean(pressures))
    mean_temperature_K = float(np.mean(temperatures_K))
    mean_energy_star = float(np.mean(total_energies))

    return {
        "rho_star": rho_star,
        "epsilon_star": epsilon_star,
        "N": N,
        "L": L,
        "energy_star": mean_energy_star,
        "energy_cgs": energy_to_cgs(mean_energy_star, epsilon),
        "pressure_star": mean_pressure_star,
        "pressure_cgs": pressure_to_cgs(mean_pressure_star, epsilon, sigma),
        "temperature_K": mean_temperature_K,
        "snapshot_positions": snapshot_positions,
        "snapshot_times": snapshot_times,
    }
