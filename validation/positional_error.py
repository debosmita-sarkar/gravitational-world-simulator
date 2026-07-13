# validation/positional_error.py
"""Integrate from a real DE440 start state and measure positional error vs
DE440 at the end epoch (spec §14, target T3)."""
import numpy as np
from physics.integrators import leapfrog_step, integrate
from physics.forces import newtonian_acceleration, acceleration_with_gr
from .horizons import build_solar_system, fetch_state_si, HORIZONS_IDS


def validate_against_de440(bodies, jd_start, years, steps_per_year, use_gr=True):
    s0 = build_solar_system(jd_start, bodies)
    n_steps = int(round(years * steps_per_year))
    h = (years * 365.25 * 86400.0) / n_steps
    accel = acceleration_with_gr if use_gr else newtonian_acceleration
    traj = integrate(s0, h=h, n_steps=n_steps, stepper=leapfrog_step, accel=accel)
    final = traj[-1]
    jd_end = jd_start + years * 365.25
    errors = {}
    for i, name in enumerate(bodies):
        pos_true, _ = fetch_state_si(HORIZONS_IDS[name], jd_end)
        errors[name] = float(np.linalg.norm(final.positions[i] - pos_true))
    return {
        "errors_m": errors,
        "max_error_m": max(errors.values()),
        "years": years, "steps": n_steps, "h_s": h, "use_gr": use_gr,
    }
