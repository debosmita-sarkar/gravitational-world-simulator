# tests/test_lyapunov.py
import numpy as np
from physics.state import State
from physics.kepler import circular_two_body
from analysis.lyapunov import benettin_lyapunov, lyapunov_time
from physics.forces import newtonian_acceleration


def _chaotic_three_body():
    # A compact, unequal-mass 3-body system: sensitive dependence -> positive lambda.
    G = 6.674e-11
    m = np.array([1.0e24, 8.0e23, 6.0e23])
    pos = np.array([[0.0, 0.0, 0.0],
                    [1.0e8, 0.0, 0.0],
                    [0.3e8, 0.9e8, 0.0]])
    vel = np.array([[0.0, 0.0, 0.0],
                    [0.0, 8.0e2, 0.0],
                    [-6.0e2, 2.0e2, 0.0]])
    return State(m, pos, vel)


def test_chaotic_system_has_positive_lyapunov():
    s = _chaotic_three_body()
    lam = benettin_lyapunov(s, h=50.0, n_steps=40000, renorm_every=10,
                            delta0=1.0, accel=newtonian_acceleration)
    assert lam > 0.0
    assert np.isfinite(lam)


def test_regular_orbit_has_near_zero_lyapunov():
    # A clean two-body circular orbit is regular -> lambda ~ 0 (much smaller
    # than the chaotic case).
    s, T = circular_two_body(m1=1.0e24, m2=1.0e22, separation=1.0e8)
    lam = benettin_lyapunov(s, h=T/2000.0, n_steps=40000, renorm_every=10,
                            delta0=1.0, accel=newtonian_acceleration)
    lam_chaos = benettin_lyapunov(_chaotic_three_body(), h=50.0, n_steps=40000,
                                  renorm_every=10, delta0=1.0,
                                  accel=newtonian_acceleration)
    assert lam < lam_chaos
