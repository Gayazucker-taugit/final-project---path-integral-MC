#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quantum Path-Integral Monte Carlo simulation of a single particle
in a one-dimensional harmonic trap, using the ring-polymer mapping
of Feynman's path integral.

The quantum partition function of a particle in a potential V(x) at
inverse temperature beta = 1/(k_B T) is mapped onto the classical
partition function of a ring polymer of P beads:

    Z_P = (m P / (2 pi beta hbar^2))^(P/2)
          x integral dx_0...dx_{P-1}
          x exp( -beta * U_eff({x_k}) )

where the effective classical potential is:

    U_eff = sum_{k=0}^{P-1} [
              (1/2) * m * omega_P^2 * (x_{k+1} - x_k)^2   
            + (1/P) * V(x_k)                                
            ]

with omega_P = sqrt(P) / (beta * hbar)  (the inter-bead spring frequency)
and periodic boundary conditions on the ring: x_P = x_0.

The total energy is estimated using the thermodynamic (primitive) estimator:

    <E> = P/(2*beta) - (1/2)*m*omega_P^2 * sum_k (x_{k+1}-x_k)^2
          + (1/P) * sum_k V(x_k)

For the harmonic trap V(x) = (1/2)*m*omega^2*x^2 the analytical result is:

    E_exact(beta) = (hbar*omega/2) * coth(beta*hbar*omega/2)

"""

import numpy as np
from scipy.constants import hbar, Boltzmann as k_B


class QuantumSimulation:

    def __init__(self, mass, omega, temp, P, Nsteps,
                 drmax, seed=937142, printfreq=100, x0=0.0):

        # Physical parameters
        self.mass  = mass
        self.omega = omega
        self.temp  = temp
        self.beta  = 1.0 / (k_B * temp)
        self.P     = P

        # omega_P = sqrt(P) / (beta * hbar)
        self.omega_P = np.sqrt(P) / (self.beta * hbar)

        # Spring constant between beads: k_spring = m * omega_P^2
        self.k_spring = mass * self.omega_P ** 2

        # MC parameters
        self.Nsteps    = Nsteps
        self.drmax     = drmax
        self.printfreq = printfreq

        # RNG
        np.random.seed(seed)

       # Bead positions: shape (P,), 1D ring polymer
        # Initialize with small random displacements drawn from the
        # approximate ground-state width sigma = sqrt(hbar/(2*m*omega))
        # so beads are not stuck at zero for large P.
        sigma_init = np.sqrt(hbar / (2.0 * self.mass * self.omega))
        self.x = np.full(P, x0, dtype=float)

        # Energy accumulators
        self.U_spring   = 0.0   # spring (kinetic) contribution to U_eff
        self.U_ext      = 0.0   # external potential contribution to U_eff
        self.U_eff      = 0.0   # total effective potential
        self.E_thermo   = 0.0   # thermodynamic energy estimator

        # Compute initial energies
        self._eval_energies()

        # Step counter and storage
        self.step      = 0
        self.E_history = []   # stores (step, E_thermo) tuples
        self.accept    = 0    # accepted moves in last MCstep call

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _spring_energy(self, x):
       
        # Compute the total spring energy of the ring polymer.
        dx = np.roll(x, -1) - x          # x_{k+1} - x_k, periodic
        return 0.5 * self.k_spring * np.sum(dx ** 2)

    def _ext_energy(self, x):
      
        # Compute the external (harmonic trap) contribution to U_eff.
        return (1.0 / self.P) * 0.5 * self.mass * self.omega ** 2 * np.sum(x ** 2)

    def _eval_energies(self):
    
        # Compute U_spring, U_ext, U_eff and the thermodynamic energy
        # estimator from the current bead positions self.x.
        # Uses the thermodynamic (primitive) estimator,
        # which is exact in the limit P -> infinity.
    
        self.U_spring = self._spring_energy(self.x)
        self.U_ext    = self._ext_energy(self.x)
        self.U_eff    = self.U_spring + self.U_ext

        dx         = np.roll(self.x, -1) - self.x
        spring_sum = np.sum(dx ** 2)

        # V_avg = (1/P) * sum_k (1/2)*m*omega^2*x_k^2
        V_avg = 0.5 * self.mass * self.omega**2 * np.mean(self.x**2)

        # Primitive estimator
        self.E_thermo = (self.P / (2.0 * self.beta)
                         - 0.5 * self.k_spring * spring_sum
                         + V_avg)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def MCstep(self):
        """
        Perform one MC sweep: attempt to move each of the P beads once.

        For each bead k, a trial displacement is drawn uniformly from
        [-drmax, drmax]. The move is accepted or rejected using the
        Metropolis criterion applied to the change in U_eff.

        Only the two spring terms connecting bead k to its neighbours
        k-1 and k+1 change when bead k moves, plus the external term
        for bead k. This allows an O(1) delta-U calculation per bead.

        During the first quarter of the run, drmax is adaptively tuned
        every 100 sweeps to keep the acceptance ratio near 0.4-0.6.
        
        """
        accept = 0
        x      = self.x          # reference, not a copy
        P      = self.P
        beta   = self.beta
        ks     = self.k_spring
        m      = self.mass
        om2    = self.omega ** 2

        for k in range(P):
            xk_old = x[k]
            xk_new = xk_old + np.random.uniform(-self.drmax, self.drmax)

            # Neighbours with periodic boundary conditions
            xk_prev = x[(k - 1) % P]
            xk_next = x[(k + 1) % P]

            # Delta spring energy (only two segments change)
            dU_spring = (0.5 * ks * (
                (xk_new - xk_prev)**2 + (xk_next - xk_new)**2
              - (xk_old - xk_prev)**2 - (xk_next - xk_old)**2
            ))

            # Delta external energy
            dU_ext = (1.0 / P) * 0.5 * m * om2 * (xk_new**2 - xk_old**2)

            dU = dU_spring + dU_ext

            # Metropolis criterion
            if dU <= 0.0 or np.log(np.random.rand()) < -beta * dU:
                x[k]          = xk_new
                self.U_eff   += dU
                accept        += 1

        self.accept = accept

        # Adaptive step size: tune drmax during first quarter of run
        # to keep acceptance ratio near 0.4-0.6
        if self.step < self.Nsteps // 4 and self.step % 100 == 0 and self.step > 0:
            ratio = self.accept / self.P
            if ratio > 0.6:
                self.drmax *= 1.1
            elif ratio < 0.4:
                self.drmax *= 0.9

        # Recompute thermodynamic estimator from scratch 
        self._eval_energies()

    def run(self):
        
        # Run the full MC simulation for self.Nsteps sweeps.

        self.E_history = []

        self.total_accept = 0
        self.total_moves  = 0

        for step in range(self.Nsteps):
            self.MCstep()
            self.step += 1
            self.total_accept += self.accept   # FIX: was commented out
            self.total_moves  += self.P        # FIX: was commented out

            # FIX: storage line was outside the loop (wrong indentation)
            if self.step % self.printfreq == 0:
                self.E_history.append((self.step, self.E_thermo))

    def mean_energy(self, burn_frac=0.2):
        
        # Return the mean thermodynamic energy after discarding the first
        # burn_frac fraction of the stored history as equilibration.
      
        if len(self.E_history) == 0:
            raise RuntimeError("No energy history found. Did you call run()?")

        energies = np.array([e for _, e in self.E_history])
        n_burn   = max(1, int(burn_frac * len(energies)))
        energies = energies[n_burn:]

        return float(np.mean(energies)), float(np.std(energies))


def analytical_energy(beta, omega):

    # Exact thermal average energy of a quantum harmonic oscillator:

    x = 0.5 * beta * hbar * omega
    return 0.5 * hbar * omega / np.tanh(x)


def choose_P(beta, omega, P_per_unit=20):
    
    # Number of beads
    # P = max(4, round(P_per_unit * beta * hbar * omega))
  
    P = int(round(P_per_unit * beta * hbar * omega))
    return max(4, P)