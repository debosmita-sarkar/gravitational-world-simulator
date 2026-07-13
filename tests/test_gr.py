# tests/test_gr.py
import numpy as np
from physics.state import State
from physics.forces import newtonian_acceleration, pn1_acceleration, acceleration_with_gr
from physics.constants import c


def _sun_planet():
    m_sun = 1.98892e30
    m_planet = 3.301e23   # Mercury-ish
    a = 5.79e10
    v = np.sqrt(6.674e-11 * m_sun / a)  # near-circular speed
    masses = np.array([m_sun, m_planet])
    positions = np.array([[0.0, 0.0, 0.0], [a, 0.0, 0.0]])
    velocities = np.array([[0.0, 0.0, 0.0], [0.0, v, 0.0]])
    return State(masses, positions, velocities)


def test_gr_term_is_small_correction_to_newtonian():
    s = _sun_planet()
    a_newt = np.linalg.norm(newtonian_acceleration(s)[1])
    a_gr = np.linalg.norm(pn1_acceleration(s)[1])
    ratio = a_gr / a_newt
    # 1PN correction is ~ (v/c)^2 ~ 1e-8 for Mercury: small but nonzero.
    assert 1e-9 < ratio < 1e-6


def test_central_body_gets_no_pn_term():
    s = _sun_planet()
    assert np.allclose(pn1_acceleration(s)[0], 0.0)


def test_acceleration_with_gr_equals_newtonian_plus_pn():
    s = _sun_planet()
    combined = acceleration_with_gr(s)
    expected = newtonian_acceleration(s) + pn1_acceleration(s)
    assert np.allclose(combined, expected)
