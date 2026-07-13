# tests/test_state.py
import numpy as np
from physics.state import (
    State, kinetic_energy, potential_energy, total_energy,
    linear_momentum, angular_momentum,
)
from physics.constants import G


def _two_body():
    # Unit-ish system: two equal masses, symmetric about origin.
    masses = np.array([1.0, 1.0])
    positions = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    velocities = np.array([[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]])
    return State(masses, positions, velocities)


def test_kinetic_energy():
    s = _two_body()
    # KE = 0.5*(1*0.25) + 0.5*(1*0.25) = 0.25
    assert kinetic_energy(s) == 0.25


def test_potential_energy():
    s = _two_body()
    # separation = 2.0; PE = -G*1*1/2
    assert np.isclose(potential_energy(s), -G * 1.0 * 1.0 / 2.0)


def test_total_energy_is_sum():
    s = _two_body()
    assert np.isclose(total_energy(s), kinetic_energy(s) + potential_energy(s))


def test_linear_momentum_zero_for_symmetric_system():
    s = _two_body()
    assert np.allclose(linear_momentum(s), np.zeros(3))


def test_angular_momentum_along_z():
    s = _two_body()
    L = angular_momentum(s)
    # r x (m v): each body contributes +z; total Lz = 1*(-1*-0.5 - 0) + 1*(1*0.5 - 0) = 1.0
    assert np.allclose(L[:2], 0.0)
    assert np.isclose(L[2], 1.0)
