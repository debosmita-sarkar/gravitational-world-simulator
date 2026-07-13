# analysis/precession.py
"""Perihelion precession via the Laplace-Runge-Lenz vector (spec §10, target T4).

Key method: the fixed-step integrator introduces its own spurious numerical
precession. Running Newtonian and Newtonian+1PN with IDENTICAL settings and
subtracting cancels that common-mode error, isolating the pure GR advance.
"""
import numpy as np
from physics.state import State
from physics.integrators import leapfrog_step, integrate
from physics.forces import newtonian_acceleration, acceleration_with_gr

_ARCSEC = 206264.806247   # radians -> arcseconds
_DAYS_PER_CENTURY = 36525.0


def lrl_perihelion_angle(pos, vel, GM) -> float:
    """Longitude (rad) of the Laplace-Runge-Lenz (eccentricity) vector."""
    r = np.linalg.norm(pos)
    L = np.cross(pos, vel)
    e_vec = np.cross(vel, L) / GM - pos / r
    return float(np.arctan2(e_vec[1], e_vec[0]))


def _two_body_state(a, e, GM_sun, m_planet):
    # Start at perihelion on +x. Sun at origin (dominant mass); barycentric drift
    # is negligible for precession and cancels in the GR-minus-Newtonian subtraction.
    from physics.constants import G
    m_sun = GM_sun / G
    r_peri = a * (1 - e)
    v_peri = np.sqrt(GM_sun * (1 + e) / (a * (1 - e)))
    masses = np.array([m_sun, m_planet])
    positions = np.array([[0.0, 0.0, 0.0], [r_peri, 0.0, 0.0]])
    velocities = np.array([[0.0, 0.0, 0.0], [0.0, v_peri, 0.0]])
    return State(masses, positions, velocities)


def _unwrap_angle_series(angles):
    return np.unwrap(np.array(angles))


def measure_precession(a, e, GM_sun, m_planet, orbits, steps_per_orbit, use_gr) -> float:
    """Precession rate in rad/orbit from a linear fit of the LRL angle vs orbit."""
    s = _two_body_state(a, e, GM_sun, m_planet)
    T = 2 * np.pi * np.sqrt(a**3 / GM_sun)   # ~Keplerian period
    h = T / steps_per_orbit
    accel = acceleration_with_gr if use_gr else newtonian_acceleration
    n_steps = orbits * steps_per_orbit
    traj = integrate(s, h=h, n_steps=n_steps, stepper=leapfrog_step, accel=accel)
    # Sample once per orbit.
    angles, orbit_index = [], []
    for k in range(0, len(traj), steps_per_orbit):
        st = traj[k]
        rel_pos = st.positions[1] - st.positions[0]
        rel_vel = st.velocities[1] - st.velocities[0]
        angles.append(lrl_perihelion_angle(rel_pos, rel_vel, GM_sun))
        orbit_index.append(k / steps_per_orbit)
    ang = _unwrap_angle_series(angles)
    slope = np.polyfit(np.array(orbit_index), ang, 1)[0]  # rad per orbit
    return float(slope)


def gr_precession_arcsec_per_century(a, e, GM_sun, m_planet, period_days,
                                     orbits, steps_per_orbit) -> dict:
    rate_gr = measure_precession(a, e, GM_sun, m_planet, orbits, steps_per_orbit, True)
    rate_newt = measure_precession(a, e, GM_sun, m_planet, orbits, steps_per_orbit, False)
    rate_diff = rate_gr - rate_newt            # rad/orbit, GR-only (common-mode removed)
    orbits_per_century = _DAYS_PER_CENTURY / period_days
    arcsec_per_century = rate_diff * orbits_per_century * _ARCSEC
    return {
        "rad_per_orbit_gr": rate_gr,
        "rad_per_orbit_newt": rate_newt,
        "rad_per_orbit_diff": rate_diff,
        "arcsec_per_century": arcsec_per_century,
        "orbits_per_century": orbits_per_century,
    }
