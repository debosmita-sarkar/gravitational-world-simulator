# tests/test_leapfrog.py
import numpy as np
from physics.state import State
from physics.integrators import leapfrog_step, integrate


def test_drift_only_when_no_force():
    # Single body, no force -> straight-line motion.
    s = State(np.array([1.0]), np.zeros((1, 3)), np.array([[1.0, 0.0, 0.0]]))
    s2 = leapfrog_step(s, h=2.0)
    assert np.allclose(s2.positions[0], [2.0, 0.0, 0.0])
    assert np.allclose(s2.velocities[0], [1.0, 0.0, 0.0])


def test_integrate_returns_trajectory_of_expected_length():
    s = State(np.array([1.0]), np.zeros((1, 3)), np.array([[1.0, 0.0, 0.0]]))
    traj = integrate(s, h=1.0, n_steps=5, stepper=leapfrog_step)
    assert len(traj) == 6
    assert np.allclose(traj[-1].positions[0], [5.0, 0.0, 0.0])
