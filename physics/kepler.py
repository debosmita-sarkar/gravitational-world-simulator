# physics/kepler.py
"""Analytic two-body oracle (spec §7, Appendix B.1)."""
import numpy as np
from .state import State
from .constants import G


def circular_two_body(m1: float, m2: float, separation: float):
    """Two masses on a circular orbit about their barycentre.

    Returns (State, period). Barycentric: total momentum is zero.
    """
    M = m1 + m2
    d = separation
    # Distances from barycentre.
    r1 = d * m2 / M
    r2 = d * m1 / M
    # Relative circular speed, then split by mass ratio.
    v_rel = np.sqrt(G * M / d)
    v1 = v_rel * m2 / M
    v2 = v_rel * m1 / M
    positions = np.array([[-r1, 0.0, 0.0], [r2, 0.0, 0.0]])
    velocities = np.array([[0.0, -v1, 0.0], [0.0, v2, 0.0]])
    masses = np.array([m1, m2])
    period = 2 * np.pi * np.sqrt(d**3 / (G * M))
    return State(masses, positions, velocities), period
