# tests/test_rk4.py
import numpy as np
from physics.state import State
from physics.integrators import rk4_step, integrate


def test_rk4_drift_only_when_no_force():
    s = State(np.array([1.0]), np.zeros((1, 3)), np.array([[1.0, 0.0, 0.0]]))
    s2 = rk4_step(s, h=2.0)
    assert np.allclose(s2.positions[0], [2.0, 0.0, 0.0])
    assert np.allclose(s2.velocities[0], [1.0, 0.0, 0.0])


def test_rk4_interchangeable_in_integrate():
    s = State(np.array([1.0]), np.zeros((1, 3)), np.array([[1.0, 0.0, 0.0]]))
    traj = integrate(s, h=1.0, n_steps=3, stepper=rk4_step)
    assert len(traj) == 4
    assert np.allclose(traj[-1].positions[0], [3.0, 0.0, 0.0])
