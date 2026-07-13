# tests/test_forces.py
import numpy as np
from physics.state import State
from physics.forces import newtonian_acceleration
from physics.constants import G


def test_two_body_acceleration_points_inward_and_matches_formula():
    masses = np.array([2.0, 3.0])
    positions = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    velocities = np.zeros((2, 3))
    s = State(masses, positions, velocities)
    a = newtonian_acceleration(s)
    # body 0 pulled toward +x by mass 3 at distance 10: a0 = G*3/100 in +x
    assert np.isclose(a[0, 0], G * 3.0 / 100.0)
    assert np.allclose(a[0, 1:], 0.0)
    # body 1 pulled toward -x by mass 2 at distance 10: a1 = -G*2/100 in x
    assert np.isclose(a[1, 0], -G * 2.0 / 100.0)


def test_isolated_body_has_zero_acceleration():
    s = State(np.array([1.0]), np.zeros((1, 3)), np.zeros((1, 3)))
    assert np.allclose(newtonian_acceleration(s), 0.0)
