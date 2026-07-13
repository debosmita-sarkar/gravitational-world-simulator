# tests/test_energy_behaviour.py
#
# Target T2 (spec §9): a symplectic integrator's energy error is BOUNDED and
# OSCILLATORY (no secular growth), while a non-symplectic RK scheme drifts
# SECULARLY (monotone growth).
#
# Design note: we deliberately do NOT test a `late/early` magnitude ratio.
# For the symplectic integrator the energy error sits at the floating-point
# roundoff floor (~1e-13 here) and oscillates; two samples of that floor give a
# meaningless ratio (it conflates the harmless Brouwer sqrt(t) roundoff
# random-walk with genuine truncation drift). The physically correct
# discriminator is the SHAPE of the error series: symplectic = oscillatory
# (many sign changes in its increments), non-symplectic = monotone (≈0).
import numpy as np
from physics.kepler import circular_two_body
from physics.integrators import leapfrog_step, rk4_step, integrate
from physics.state import total_energy


def _energy_error_series(stepper, n_orbits=200, steps_per_orbit=500):
    s, T = circular_two_body(m1=1.0e24, m2=1.0e22, separation=1.0e8)
    h = T / steps_per_orbit
    traj = integrate(s, h=h, n_steps=n_orbits * steps_per_orbit, stepper=stepper)
    e0 = total_energy(traj[0])
    errs = np.array([abs((total_energy(st) - e0) / e0) for st in traj[::steps_per_orbit]])
    return errs


def _sign_changes(errs):
    """Number of sign changes in the increments of the error series.

    High  -> error oscillates (bounded, symplectic).
    ~zero -> error grows monotonically (secular, non-symplectic).
    """
    d = np.diff(errs[1:])           # drop t=0 (error is exactly 0 there)
    return int(np.sum(np.diff(np.sign(d)) != 0))


def test_leapfrog_energy_is_bounded_and_oscillatory():
    errs = _energy_error_series(leapfrog_step)
    # Bounded: never approaches any physically meaningful magnitude.
    assert errs.max() < 1e-4
    # Non-secular: the error oscillates (returns toward zero) rather than
    # growing monotonically. Many sign changes prove a bounded envelope.
    n = len(errs) - 2               # number of increment-pairs available
    assert _sign_changes(errs) > 0.2 * n


def test_rk4_energy_drifts_secularly():
    errs = _energy_error_series(rk4_step)
    # Monotone secular drift: almost no sign changes, and the late error
    # greatly exceeds the early error.
    n = len(errs) - 2
    assert _sign_changes(errs) < 0.05 * n
    assert errs[-10:].mean() > 5 * (errs[1:10].mean() + 1e-15)


def test_symplectic_beats_rk4_on_boundedness():
    """Direct contrast: leapfrog stays orders of magnitude tighter than RK4."""
    leap = _energy_error_series(leapfrog_step)
    rk = _energy_error_series(rk4_step)
    assert leap.max() < rk.max()
