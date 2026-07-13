# physics/state.py
"""Barycentric N-body state and conserved-quantity evaluators (spec §6)."""
from dataclasses import dataclass
import numpy as np
from .constants import G


@dataclass(frozen=True)
class State:
    masses: np.ndarray      # (N,)
    positions: np.ndarray   # (N, 3)
    velocities: np.ndarray  # (N, 3)


def kinetic_energy(s: State) -> float:
    v2 = np.einsum("ij,ij->i", s.velocities, s.velocities)
    return float(0.5 * np.sum(s.masses * v2))


def potential_energy(s: State) -> float:
    total = 0.0
    n = len(s.masses)
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(s.positions[j] - s.positions[i])
            total -= G * s.masses[i] * s.masses[j] / d
    return float(total)


def total_energy(s: State) -> float:
    return kinetic_energy(s) + potential_energy(s)


def linear_momentum(s: State) -> np.ndarray:
    return np.sum(s.masses[:, None] * s.velocities, axis=0)


def angular_momentum(s: State) -> np.ndarray:
    return np.sum(s.masses[:, None] * np.cross(s.positions, s.velocities), axis=0)
