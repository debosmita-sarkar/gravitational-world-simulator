# physics/forces.py
"""Gravitational accelerations (spec §5, §10)."""
import numpy as np
from .state import State
from .constants import G, c


def newtonian_acceleration(s: State) -> np.ndarray:
    """a_i = G * sum_{j!=i} m_j (r_j - r_i) / |r_j - r_i|^3."""
    pos = s.positions
    n = len(s.masses)
    acc = np.zeros_like(pos)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = pos[j] - pos[i]
            r = np.linalg.norm(d)
            acc[i] += G * s.masses[j] * d / r**3
    return acc


def pn1_acceleration(s: State, central_index: int = 0) -> np.ndarray:
    """Leading Schwarzschild 1PN correction about a dominant central mass (spec §10).

    a_GR = (GM / (c^2 r^2)) * [ (4GM/r - v^2) r_hat + 4 (v . r_hat) v ]
    Applied to each non-central body relative to the central body.
    """
    pos, vel = s.positions, s.velocities
    n = len(s.masses)
    acc = np.zeros_like(pos)
    GM = G * s.masses[central_index]
    rc = pos[central_index]
    vc = vel[central_index]
    for i in range(n):
        if i == central_index:
            continue
        d = pos[i] - rc
        r = np.linalg.norm(d)
        rhat = d / r
        v = vel[i] - vc
        v2 = float(np.dot(v, v))
        vr = float(np.dot(v, rhat))
        acc[i] = (GM / (c**2 * r**2)) * ((4.0 * GM / r - v2) * rhat + 4.0 * vr * v)
    return acc


def acceleration_with_gr(s: State, central_index: int = 0) -> np.ndarray:
    return newtonian_acceleration(s) + pn1_acceleration(s, central_index)
