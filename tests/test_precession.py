# tests/test_precession.py
import numpy as np
from analysis.precession import (
    lrl_perihelion_angle, measure_precession, gr_precession_arcsec_per_century,
)


def test_lrl_angle_points_to_perihelion_on_x_axis():
    # At perihelion on +x with velocity +y (bound orbit), LRL points +x -> angle 0.
    GM = 1.32712440018e20
    a, e = 5.79e10, 0.2056
    r_peri = a * (1 - e)
    v_peri = np.sqrt(GM * (1 + e) / (a * (1 - e)))
    pos = np.array([r_peri, 0.0, 0.0])
    vel = np.array([0.0, v_peri, 0.0])
    ang = lrl_perihelion_angle(pos, vel, GM)
    assert abs(ang) < 1e-6


def test_newtonian_two_body_has_negligible_physical_precession():
    # Pure Newtonian ellipse does not precess (only tiny numerical drift).
    # Leapfrog's spurious precession is O(h^2); at steps_per_orbit=4000 it is
    # ~4.5e-6 rad/orbit, above the 1e-6 bar. Resolution bump (Task 4 sanctioned)
    # to 12000 steps/orbit reduces the numerical artifact to ~5e-7 rad/orbit.
    # The tolerance is NOT widened; only the integrator resolution is increased.
    rate = measure_precession(a=5.79e10, e=0.2056, GM_sun=1.32712440018e20,
                              m_planet=3.301e23, orbits=50,
                              steps_per_orbit=12000, use_gr=False)
    # rad/orbit numerical precession should be very small.
    assert abs(rate) < 1e-6


def test_mercury_gr_precession_matches_43_arcsec():
    out = gr_precession_arcsec_per_century(
        a=5.79e10, e=0.2056, GM_sun=1.32712440018e20, m_planet=3.301e23,
        period_days=87.969, orbits=60, steps_per_orbit=4000)
    # The GR-minus-Newtonian advance should land near 43"/century.
    assert 38.0 < out["arcsec_per_century"] < 48.0
