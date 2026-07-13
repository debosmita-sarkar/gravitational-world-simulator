# physics/integrators.py
"""Time integrators (spec §9). Leapfrog is symplectic; RK4 (Task 5) is not."""
from dataclasses import replace
import numpy as np
from .state import State
from .forces import newtonian_acceleration


def leapfrog_step(s: State, h: float, accel=newtonian_acceleration) -> State:
    """Kick-drift-kick velocity Verlet — symplectic, 2nd order."""
    a = accel(s)
    v_half = s.velocities + 0.5 * h * a
    new_pos = s.positions + h * v_half
    s_mid = replace(s, positions=new_pos)
    a_new = accel(s_mid)
    new_vel = v_half + 0.5 * h * a_new
    return State(s.masses, new_pos, new_vel)


def integrate(s: State, h: float, n_steps: int, stepper,
              accel=newtonian_acceleration) -> list:
    traj = [s]
    cur = s
    for _ in range(n_steps):
        cur = stepper(cur, h, accel)
        traj.append(cur)
    return traj


def rk4_step(s: State, h: float, accel=newtonian_acceleration) -> State:
    """Classical RK4 on (position, velocity). NOT symplectic — energy drifts."""
    def deriv(pos, vel):
        a = accel(State(s.masses, pos, vel))
        return vel, a

    p0, v0 = s.positions, s.velocities
    k1p, k1v = deriv(p0, v0)
    k2p, k2v = deriv(p0 + 0.5 * h * k1p, v0 + 0.5 * h * k1v)
    k3p, k3v = deriv(p0 + 0.5 * h * k2p, v0 + 0.5 * h * k2v)
    k4p, k4v = deriv(p0 + h * k3p, v0 + h * k3v)
    new_pos = p0 + (h / 6.0) * (k1p + 2 * k2p + 2 * k3p + k4p)
    new_vel = v0 + (h / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
    return State(s.masses, new_pos, new_vel)
