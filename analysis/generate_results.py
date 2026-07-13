# analysis/generate_results.py
"""Reproducible generation of all Plan 2 Results numbers. Writes results.json.

Run: python -m analysis.generate_results
Requires network on first run (Horizons); cached thereafter.
"""
import json
import numpy as np
from physics.constants import GM_BODIES
from analysis.precession import gr_precession_arcsec_per_century
from analysis.lyapunov import benettin_lyapunov, lyapunov_time
from analysis.fidelity import horizon
from physics.forces import newtonian_acceleration
from tests.test_lyapunov import _chaotic_three_body

out = {}

# T4: Mercury GR precession (no network needed).
out["T4_precession"] = gr_precession_arcsec_per_century(
    a=5.79e10, e=0.2056, GM_sun=GM_BODIES["sun"], m_planet=3.301e23,
    period_days=87.969, orbits=80, steps_per_orbit=6000)

# T5: Lyapunov estimator validated on controlled chaos.
s = _chaotic_three_body()
lam = benettin_lyapunov(s, h=50.0, n_steps=60000, renorm_every=10, delta0=1.0,
                        accel=newtonian_acceleration)
out["T5_lyapunov"] = {"lambda_per_s": lam, "lyapunov_time_s": (1.0/lam if lam>0 else None)}

# T6: horizon formula + regime crossover, using the measured lambda above.
sweep = []
for h in [1.0, 0.5, 0.25, 0.125]:
    t = horizon(lam=lam, delta_tol=1e6, delta0=1e-6, C=1.0, h=h, p=2)
    sweep.append({"h": h, "t_lim": t})
out["T6_horizon_sweep"] = sweep

# T3: positional error vs DE440 (network). Wrapped so offline still writes T4-6.
try:
    from validation.positional_error import validate_against_de440
    out["T3_positional_error"] = validate_against_de440(
        bodies=["sun", "mercury", "venus", "emb", "mars", "jupiter"],
        jd_start=2451545.0, years=10.0, steps_per_year=8760, use_gr=True)
except Exception as e:
    out["T3_positional_error"] = {"error": f"{type(e).__name__}: {e}"}

with open("results.json", "w") as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
