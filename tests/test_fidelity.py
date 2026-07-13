# tests/test_fidelity.py
import numpy as np
from analysis.fidelity import horizon, empirical_divergence_rate
from analysis.lyapunov import benettin_lyapunov
from tests.test_lyapunov import _chaotic_three_body
from physics.forces import newtonian_acceleration


def test_horizon_formula_matches_hand_calculation():
    # Measurement-dominated: t_lim ~ (1/lam) ln(tol/delta0).
    lam = 1.0 / 5.0e6            # 1/yr-ish scale placeholder
    t = horizon(lam=lam, delta_tol=1.496e11, delta0=15.0, C=0.0, h=1.0, p=2)
    expected = (1.0 / lam) * np.log(1.496e11 / 15.0)
    assert np.isclose(t, expected, rtol=1e-9)


def test_smaller_step_extends_horizon_only_logarithmically():
    lam = 1e-3
    t_coarse = horizon(lam=lam, delta_tol=1e6, delta0=1e-9, C=1.0, h=1.0, p=2)
    t_fine = horizon(lam=lam, delta_tol=1e6, delta0=1e-9, C=1.0, h=0.5, p=2)
    gain = t_fine - t_coarse
    # halving h adds ~ p*ln2/lam (Corollary 1), NOT a large multiplicative gain
    assert np.isclose(gain, 2 * np.log(2) / lam, rtol=0.2)


def test_empirical_divergence_matches_lyapunov():
    s = _chaotic_three_body()
    lam = benettin_lyapunov(s, h=50.0, n_steps=40000, renorm_every=10,
                            delta0=1.0, accel=newtonian_acceleration)
    rate = empirical_divergence_rate(s, h=50.0, n_steps=3000, delta0=1.0,
                                     accel=newtonian_acceleration)
    # Same order of magnitude as the Benettin estimate.
    assert rate > 0
    assert 0.1 < rate / lam < 10.0
