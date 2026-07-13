# tests/test_kepler_oracle.py
import numpy as np
from physics.kepler import circular_two_body
from physics.integrators import leapfrog_step, integrate
from physics.state import total_energy


def test_period_formula_matches_keplers_third_law():
    from physics.constants import G
    s, T = circular_two_body(m1=1.0e24, m2=1.0e22, separation=1.0e8)
    d = 1.0e8
    T_expected = 2 * np.pi * np.sqrt(d**3 / (G * (1.0e24 + 1.0e22)))
    assert np.isclose(T, T_expected, rtol=1e-12)


def test_body_returns_near_start_after_one_period():
    s, T = circular_two_body(m1=1.0e24, m2=1.0e22, separation=1.0e8)
    h = T / 20000.0
    traj = integrate(s, h=h, n_steps=20000, stepper=leapfrog_step)
    start, end = traj[0], traj[-1]
    sep0 = np.linalg.norm(start.positions[1] - start.positions[0])
    drift = np.linalg.norm(end.positions[0] - start.positions[0])
    # After one period the leapfrog orbit closes to well under 0.1% of the orbit size.
    assert drift / sep0 < 1e-3


def test_energy_conserved_over_one_period():
    s, T = circular_two_body(m1=1.0e24, m2=1.0e22, separation=1.0e8)
    h = T / 20000.0
    traj = integrate(s, h=h, n_steps=20000, stepper=leapfrog_step)
    e0 = total_energy(traj[0])
    e1 = total_energy(traj[-1])
    assert abs((e1 - e0) / e0) < 1e-6
