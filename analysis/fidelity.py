# analysis/fidelity.py
"""Fidelity-Limit Theorem evaluation (spec §15, target T6)."""
import numpy as np
from dataclasses import replace
from physics.integrators import leapfrog_step
from physics.forces import newtonian_acceleration


def horizon(lam, delta_tol, delta0, C, h, p) -> float:
    """t_lim = (1/lam) * ln( delta_tol / (delta0 + C h^p / lam) )."""
    seed = delta0 + C * h**p / lam
    return float((1.0 / lam) * np.log(delta_tol / seed))


def empirical_divergence_rate(s0, h, n_steps, delta0, accel=newtonian_acceleration) -> float:
    """Fit log(separation) vs time for a lightly-perturbed trajectory; slope ~ lambda."""
    ref = s0
    pos = s0.positions.copy(); pos[0, 0] += delta0
    sh = replace(s0, positions=pos)
    times, logsep = [], []
    for step in range(n_steps):
        ref = leapfrog_step(ref, h, accel)
        sh = leapfrog_step(sh, h, accel)
        dp = (ref.positions - sh.positions).ravel()
        dv = (ref.velocities - sh.velocities).ravel()
        d = np.sqrt(np.dot(dp, dp) + np.dot(dv, dv))
        if d > 0:
            times.append((step + 1) * h)
            logsep.append(np.log(d))
    times = np.array(times); logsep = np.array(logsep)
    # Fit only the initial linear-growth window before saturation.
    k = max(10, len(times) // 3)
    slope = np.polyfit(times[:k], logsep[:k], 1)[0]
    return float(slope)
