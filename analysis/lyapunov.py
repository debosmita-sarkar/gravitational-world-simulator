# analysis/lyapunov.py
"""Largest Lyapunov exponent via Benettin shadow-trajectory renormalisation
(spec §11, target T5)."""
import numpy as np
from dataclasses import replace
from physics.integrators import leapfrog_step
from physics.forces import newtonian_acceleration


def _perturb(s, delta0):
    # Perturb the first body's x-position by delta0.
    pos = s.positions.copy()
    pos[0, 0] += delta0
    return replace(s, positions=pos)


def _separation(a, b):
    dp = (a.positions - b.positions).ravel()
    dv = (a.velocities - b.velocities).ravel()
    return np.sqrt(np.dot(dp, dp) + np.dot(dv, dv))


def benettin_lyapunov(s0, h, n_steps, renorm_every, delta0, accel=newtonian_acceleration) -> float:
    """lambda = (1/total_time) * sum log(d_i / delta0), renormalising the shadow
    trajectory back to separation delta0 every `renorm_every` steps."""
    ref = s0
    sh = _perturb(s0, delta0)
    log_sum = 0.0
    total_time = 0.0
    for step in range(n_steps):
        ref = leapfrog_step(ref, h, accel)
        sh = leapfrog_step(sh, h, accel)
        if (step + 1) % renorm_every == 0:
            d = _separation(ref, sh)
            if d > 0:
                log_sum += np.log(d / delta0)
                # rescale shadow back toward ref to keep separation ~ delta0
                scale = delta0 / d
                new_pos = ref.positions + (sh.positions - ref.positions) * scale
                new_vel = ref.velocities + (sh.velocities - ref.velocities) * scale
                sh = replace(sh, positions=new_pos, velocities=new_vel)
            total_time += renorm_every * h
    return float(log_sum / total_time) if total_time > 0 else 0.0


def lyapunov_time(*args, **kwargs) -> float:
    lam = benettin_lyapunov(*args, **kwargs)
    return float("inf") if lam <= 0 else 1.0 / lam
