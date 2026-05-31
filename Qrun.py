#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, Boltzmann as k_B, atomic_mass, e as e_charge

from Qsim import QuantumSimulation, analytical_energy, choose_P

# ============================================================
# Physical constants and system parameters
# ============================================================

m_Ar   = 39.948 * atomic_mass          # kg
hw_eV  = 50e-3                          # eV  (hbar * omega = 50 meV)
hw_J   = hw_eV * e_charge              # J
omega  = hw_J / hbar                   # rad/s

E0_meV = hw_eV * 1e3 / 2.0            # ground state energy = 25 meV
E0_J   = hw_J / 2.0                   # ground state energy in J

print(f"System: Ar atom, hbar*omega = {hw_eV*1e3:.1f} meV")
print(f"Ground state energy E0 = {E0_meV:.2f} meV")
print(f"Crossover temperature  = {hw_J/k_B:.1f} K\n")

# ============================================================
# Simulation settings
# ============================================================

# Dimensionless beta values 
bhw_list  = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])  # beta * hbar * omega

# Convert to physical beta (J^-1) and temperature (K)
beta_list = bhw_list / hw_J
T_list    = 1.0 / (k_B * beta_list)

Nsteps     = 50000       # TODO: change back to 50000 for production
burn_frac  = 0.30       # discard first 30 % as equilibration
printfreq  = 100        # store energy every printfreq steps
N_indep    = 5          # TODO: change back to 5 for production
P_per_unit = 20         # beads per unit of beta*hbar*omega
drmax_m    = 0.5e-10    # initial step size (m) - 0.5 Angstrom
x0_m       = 0.0        # initial bead positions

# ============================================================
# Helper: run N_indep independent simulations at one temperature
# ============================================================

def run_temperature(temp, P, Nsteps, drmax, burn_frac, N_indep, base_seed=937142):
    """
    Run N_indep independent PIMC simulations at the given temperature
    with P beads. Returns an array of mean energies (J), one per run.
    """
    mean_Es = []

    for i in range(N_indep):
        seed = base_seed + i * 1000   

        # Scale initial drmax with P (stiffer springs need smaller steps)
        drmax_scaled = drmax / np.sqrt(P)

        # Large P rings need more steps to equilibrate from x0=0
        Nsteps_scaled = max(Nsteps, int(Nsteps * P / 20))

        sim = QuantumSimulation(
            mass      = m_Ar,
            omega     = omega,
            temp      = temp,
            P         = P,
            Nsteps    = Nsteps_scaled,
            drmax     = drmax_scaled,
            seed      = seed,
            printfreq = printfreq,
            x0        = x0_m,
        )

        sim.run()

        # print acceptance ratio
        acceptance_ratio = sim.total_accept / sim.total_moves
        print(f"    run {i+1}: acceptance ratio = {acceptance_ratio:.3f}")

        # Larger P needs longer equilibration - scale burn fraction with P
        burn_frac_scaled = min(0.5, burn_frac * np.sqrt(P / 20))
        mean_E, _ = sim.mean_energy(burn_frac=burn_frac_scaled)
        mean_Es.append(mean_E)

    return np.array(mean_Es)


# ============================================================
# Part 1: Energy vs beta and temperature
# ============================================================

print("=" * 60)
print("Part 1: Mean energy vs beta*hbar*omega")
print("=" * 60)

P_list = np.array([choose_P(b, omega, P_per_unit) for b in beta_list])

print(f"\n{'beta*hbar*om':>13}  {'T (K)':>8}  {'P (beads)':>10}")
print("-" * 36)
for bhw, T, P in zip(bhw_list, T_list, P_list):
    print(f"{bhw:13.1f}  {T:8.1f}  {P:10d}")
print()

# Storage
all_mean_E = []
all_err_E  = []

for T, P, bhw in zip(T_list, P_list, bhw_list):
    print(f"beta*hbar*omega = {bhw:.1f},  T = {T:.1f} K,  P = {P} beads ...",
          flush=True)

    means = run_temperature(T, P, Nsteps, drmax_m, burn_frac, N_indep)

    all_mean_E.append(np.mean(means))
    all_err_E.append(np.std(means))

    print(f"  E = {np.mean(means)/e_charge*1e3:.3f} +/- "
          f"{np.std(means)/e_charge*1e3:.3f} meV")

all_mean_E   = np.array(all_mean_E)
all_err_E    = np.array(all_err_E)
E_exact_at_T = analytical_energy(beta_list, omega)

# Fine grids for reference curves
bhw_fine     = np.linspace(0.1, 7.0, 500)
beta_fine    = bhw_fine / hw_J
T_fine       = 1.0 / (k_B * beta_fine)
E_exact_fine = analytical_energy(beta_fine, omega)
E_class_fine = k_B * T_fine              # classical equipartition: E = k_B T

# Figure 1: Energy vs beta*hbar*omega 
fig1, ax1 = plt.subplots(figsize=(9, 6))

ax1.fill_between(bhw_list,
                 (all_mean_E - all_err_E) / e_charge * 1e3,
                 (all_mean_E + all_err_E) / e_charge * 1e3,
                 alpha=0.3, color="tab:blue",
                 label="PIMC statistical uncertainty")

ax1.plot(bhw_list, all_mean_E / e_charge * 1e3,
         "o-", color="tab:blue", label="Path Integral MC simulation", zorder=5)

ax1.plot(bhw_fine, E_exact_fine / e_charge * 1e3,
         "k-", label="Exact quantum result", linewidth=2)

ax1.plot(bhw_fine, E_class_fine / e_charge * 1e3,
         "r--", label=r"Classical limit ($k_BT$)", linewidth=1.5)

ax1.axhline(E0_meV, color="gray", linestyle=":", linewidth=1.5,
            label=f"Ground state energy $E_0$ = {E0_meV:.1f} meV")

ax1.set_xlabel(r"$\beta\hbar\omega$ (dimensionless)", fontsize=13)
ax1.set_ylabel("Mean energy (meV)", fontsize=13)
ax1.set_title("Quantum harmonic oscillator: PIMC vs exact result", fontsize=13)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("quantum_fig1_energy_vs_beta.png", dpi=200)
print("\nSaved: quantum_fig1_energy_vs_beta.png")

# Figure 2: Energy vs temperature 
fig2, ax2 = plt.subplots(figsize=(9, 6))

ax2.fill_between(T_list,
                 (all_mean_E - all_err_E) / e_charge * 1e3,
                 (all_mean_E + all_err_E) / e_charge * 1e3,
                 alpha=0.3, color="tab:blue",
                 label="PIMC statistical uncertainty")

ax2.plot(T_list, all_mean_E / e_charge * 1e3,
         "o-", color="tab:blue", label="Path Integral MC simulation", zorder=5)

ax2.plot(T_fine, E_exact_fine / e_charge * 1e3,
         "k-", label="Exact quantum result", linewidth=2)

ax2.plot(T_fine, E_class_fine / e_charge * 1e3,
         "r--", label=r"Classical limit ($k_BT$)", linewidth=1.5)

ax2.axhline(E0_meV, color="gray", linestyle=":", linewidth=1.5,
            label=f"Ground state energy $E_0$ = {E0_meV:.1f} meV")

ax2.set_xlabel("Temperature (K)", fontsize=13)
ax2.set_ylabel("Mean energy (meV)", fontsize=13)
ax2.set_title("Quantum harmonic oscillator: PIMC vs exact result", fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("quantum_fig2_energy_vs_T.png", dpi=200)
print("Saved: quantum_fig2_energy_vs_T.png")

# Zoomed temperature grid for Figure 2 

T_zoom_min = 0.0
T_zoom_max = max(T_list) * 1.05

T_zoom = np.linspace(max(1e-6, T_zoom_min + 1e-6), T_zoom_max, 500)
beta_zoom = 1.0 / (k_B * T_zoom)
E_exact_zoom = analytical_energy(beta_zoom, omega)
E_class_zoom = k_B * T_zoom

# Figure 2: Energy vs temperature 

fig2, ax2 = plt.subplots(figsize=(9, 6))
 
ax2.fill_between(T_list,
                 (all_mean_E - all_err_E) / e_charge * 1e3,
                 (all_mean_E + all_err_E) / e_charge * 1e3,
                 alpha=0.3, color="tab:blue",
                 label="PIMC statistical uncertainty")
 
ax2.plot(T_list, all_mean_E / e_charge * 1e3,
         "o-", color="tab:blue", label="Path Integral MC simulation", zorder=5)
 
ax2.plot(T_zoom, E_exact_zoom / e_charge * 1e3,
         "k-", label="Exact quantum result", linewidth=2)
 
ax2.plot(T_zoom, E_class_zoom / e_charge * 1e3,
         "r--", label=r"Classical limit ($k_BT$)", linewidth=1.5)
 
ax2.axhline(E0_meV, color="gray", linestyle=":", linewidth=1.5,
            label=f"Ground state energy $E_0$ = {E0_meV:.1f} meV")
 
ax2.set_xlim(T_zoom_min, T_zoom_max)
ax2.set_xlabel("Temperature (K)", fontsize=13)
ax2.set_ylabel("Mean energy (meV)", fontsize=13)
ax2.set_title("Quantum harmonic oscillator: PIMC vs exact result, zoomed in", fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("quantum_fig2_energy_vs_T_zoomed.png", dpi=200)
print("Saved: quantum_fig2_energy_vs_T_zoomed.png")


# ============================================================
# Part 2: Bead convergence at the highest beta (lowest temperature)
# ============================================================

print("\n" + "=" * 60)
print(f"Part 2: Bead convergence at beta*hbar*omega = {bhw_list[-1]:.1f} "
      f"(T = {T_list[-1]:.1f} K)")
print("=" * 60)

T_low       = T_list[-1]
beta_low    = beta_list[-1]
E_exact_low = analytical_energy(beta_low, omega)

P_ref  = choose_P(beta_low, omega, P_per_unit)
P_vals = sorted(set([
    4,
    8,
    12,
    max(4, P_ref // 4),
    max(4, P_ref // 2),
    P_ref,
    int(1.5 * P_ref),
    2 * P_ref,
]))

print(f"\nReference P = {P_ref}, testing P = {P_vals}\n")

conv_mean = []
conv_err  = []

for P in P_vals:
    print(f"  P = {P} ...", flush=True)
    means = run_temperature(T_low, P, Nsteps, drmax_m, burn_frac, N_indep)
    conv_mean.append(np.mean(means))
    conv_err.append(np.std(means))
    print(f"    E = {np.mean(means)/e_charge*1e3:.3f} +/- "
          f"{np.std(means)/e_charge*1e3:.3f} meV")

conv_mean = np.array(conv_mean)
conv_err  = np.array(conv_err)

# Figure 3: Bead convergence 
fig3, ax3 = plt.subplots(figsize=(9, 6))

ax3.fill_between(P_vals,
                 (conv_mean - conv_err) / e_charge * 1e3,
                 (conv_mean + conv_err) / e_charge * 1e3,
                 alpha=0.3, color="tab:green",
                 label="PIMC statistical uncertainty")

ax3.plot(P_vals, conv_mean / e_charge * 1e3,
         "s-", color="tab:green", label="Path Integral MC simulation", zorder=5)

ax3.axhline(E_exact_low / e_charge * 1e3, color="k", linestyle="-",
            linewidth=2,
            label=f"Exact quantum result: {E_exact_low/e_charge*1e3:.3f} meV")

ax3.axhline(E0_meV, color="gray", linestyle=":", linewidth=1.5,
            label=f"Ground state energy $E_0$ = {E0_meV:.1f} meV")

ax3.set_xlabel("Number of beads $P$", fontsize=13)
ax3.set_ylabel("Mean energy (meV)", fontsize=13)
ax3.set_title(
    f"Convergence with number of beads at "
    r"$\beta\hbar\omega$" + f" = {bhw_list[-1]:.1f}  (T = {T_low:.1f} K)",
    fontsize=13)
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("quantum_fig3_bead_convergence.png", dpi=200)
print("\nSaved: quantum_fig3_bead_convergence.png")


# ============================================================
# Part 3: Ring polymer snapshots at selected beta values
# ============================================================

print("\n" + "=" * 60)
print("Part 3: Ring polymer snapshots")
print("=" * 60)

snapshot_bhw = [1.0, 3.0, 6.0]
colors        = ["tab:red", "tab:orange", "tab:blue"]

fig4, axes = plt.subplots(1, 3, figsize=(14, 5))

for ax, bhw, col in zip(axes, snapshot_bhw, colors):
    beta_snap = bhw / hw_J
    T_snap    = 1.0 / (k_B * beta_snap)
    P_snap    = choose_P(beta_snap, omega, P_per_unit)

    sim_snap = QuantumSimulation(
        mass      = m_Ar,
        omega     = omega,
        temp      = T_snap,
        P         = P_snap,
        Nsteps    = Nsteps,
        drmax     = drmax_m / np.sqrt(P_snap),
        seed      = 42,
        printfreq = printfreq,
        x0        = x0_m,
    )
    sim_snap.run()

    # Plot bead positions along imaginary time axis
    x_beads  = sim_snap.x / 1e-10        # convert to Angstrom
    tau_axis = np.arange(P_snap)

    ax.plot(tau_axis, x_beads, "o-", color=col, markersize=3, linewidth=1)

    # Close the ring with a dashed line
    ax.plot([tau_axis[-1], tau_axis[0] + P_snap],
            [x_beads[-1], x_beads[0]],
            "--", color=col, linewidth=1, alpha=0.5)

    ax.set_xlabel("Bead index (imaginary time)", fontsize=11)
    ax.set_ylabel("x (Å)", fontsize=11)
    ax.set_title(
        r"$\beta\hbar\omega$" + f" = {bhw:.1f},  T = {T_snap:.0f} K,  P = {P_snap}",
        fontsize=11)
    ax.grid(True, alpha=0.3)

fig4.suptitle("Ring polymer bead positions at selected temperatures", fontsize=13)
plt.tight_layout()
plt.savefig("quantum_fig4_ring_polymer_snapshots.png", dpi=200)
print("Saved: quantum_fig4_ring_polymer_snapshots.png")


# ============================================================
# Summary table
# ============================================================

print("\n" + "=" * 70)
print(f"{'beta*hbar*om':>13}  {'T (K)':>8}  {'P':>5}  "
      f"{'E_PIMC (meV)':>14}  {'err (meV)':>10}  {'E_exact (meV)':>14}")
print("-" * 70)

for bhw, T, P, E_mc, err, E_ex in zip(bhw_list, T_list, P_list,
                                        all_mean_E, all_err_E,
                                        E_exact_at_T):
    print(f"{bhw:13.1f}  {T:8.1f}  {P:5d}  "
          f"{E_mc/e_charge*1e3:14.3f}  "
          f"{err/e_charge*1e3:10.3f}  "
          f"{E_ex/e_charge*1e3:14.3f}")

print("\nAll figures saved:")
print("  quantum_fig1_energy_vs_beta.png")
print("  quantum_fig2_energy_vs_T.png")
print("  quantum_fig2_energy_vs_T_zoom.png")
print("  quantum_fig3_bead_convergence.png")
print("  quantum_fig4_ring_polymer_snapshots.png")

plt.show()