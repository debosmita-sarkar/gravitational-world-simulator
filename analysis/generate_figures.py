"""Generate all manuscript figures from real simulator output.

Run: python -m analysis.generate_figures
Writes PNGs into paper/figures/. Uses the non-interactive Agg backend.
Every figure is produced from a real integration or from results.json — no
synthetic or illustrative data.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from physics.kepler import circular_two_body
from physics.integrators import leapfrog_step, rk4_step, integrate
from physics.state import total_energy
from physics.constants import GM_BODIES
from analysis.precession import lrl_perihelion_angle, _two_body_state
from physics.integrators import integrate as _integrate
from physics.forces import newtonian_acceleration, acceleration_with_gr

FIGDIR = os.path.join(os.path.dirname(__file__), os.pardir, "paper", "figures")
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150, "font.size": 10, "axes.grid": True,
    "grid.alpha": 0.3, "savefig.bbox": "tight",
})


def fig_energy_behaviour():
    """Fig 1 (T2): leapfrog bounded/oscillatory vs RK4 secular energy drift."""
    s, T = circular_two_body(m1=1.0e24, m2=1.0e22, separation=1.0e8)
    n_orbits, spo = 200, 500
    h = T / spo
    out = {}
    for name, step in [("Leapfrog (symplectic)", leapfrog_step),
                       ("RK4 (non-symplectic)", rk4_step)]:
        traj = integrate(s, h=h, n_steps=n_orbits * spo, stepper=step)
        e0 = total_energy(traj[0])
        errs = np.array([abs((total_energy(st) - e0) / e0) for st in traj[::spo]])
        out[name] = errs
    orbits = np.arange(len(next(iter(out.values()))))
    plt.figure(figsize=(6, 4))
    for name, errs in out.items():
        plt.semilogy(orbits, np.maximum(errs, 1e-18), label=name, lw=1.2)
    plt.xlabel("Orbit number")
    plt.ylabel(r"Relative energy error $|\Delta E / E_0|$")
    plt.title("Energy conservation: symplectic vs non-symplectic")
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(FIGDIR, "fig_energy_behaviour.png"))
    plt.close()
    print("wrote fig_energy_behaviour.png")


def fig_precession():
    """Fig 2 (T4): perihelion longitude accumulation, GR vs Newtonian."""
    a, e = 5.79e10, 0.2056
    GM = GM_BODIES["sun"]
    m_planet = 3.301e23
    orbits, spo = 80, 4000
    T = 2 * np.pi * np.sqrt(a**3 / GM)
    h = T / spo
    plt.figure(figsize=(6, 4))
    for label, accel, color in [("Newtonian + 1PN (GR)", acceleration_with_gr, "C0"),
                                ("Newtonian only", newtonian_acceleration, "C1")]:
        s = _two_body_state(a, e, GM, m_planet)
        traj = _integrate(s, h=h, n_steps=orbits * spo, stepper=leapfrog_step, accel=accel)
        angs, idx = [], []
        for k in range(0, len(traj), spo):
            st = traj[k]
            rp = st.positions[1] - st.positions[0]
            rv = st.velocities[1] - st.velocities[0]
            angs.append(lrl_perihelion_angle(rp, rv, GM))
            idx.append(k / spo)
        ang = np.unwrap(np.array(angs))
        # convert to arcsec, relative to start
        arcsec = (ang - ang[0]) * 206264.806247
        plt.plot(idx, arcsec, label=label, color=color, lw=1.3)
    plt.xlabel("Orbit number")
    plt.ylabel(r"Perihelion longitude shift (arcsec)")
    plt.title("Perihelion precession: GR isolates the secular advance")
    plt.legend(loc="best")
    plt.savefig(os.path.join(FIGDIR, "fig_precession.png"))
    plt.close()
    print("wrote fig_precession.png")


def fig_horizon_sweep(results):
    """Fig 3 (T6): predictability horizon vs step size (logarithmic return)."""
    sweep = results["T6_horizon_sweep"]
    hs = np.array([d["h"] for d in sweep])
    tlim = np.array([d["t_lim"] for d in sweep])
    halvings = -np.log2(hs / hs[0])  # 0,1,2,3
    plt.figure(figsize=(6, 4))
    plt.plot(halvings, tlim, "o-", lw=1.3)
    # annotate the constant per-halving gain
    gain = np.diff(tlim).mean()
    plt.xlabel(r"Number of step-size halvings ($h \to h/2$)")
    plt.ylabel(r"Predictability horizon $t_{\mathrm{lim}}$ (s)")
    plt.title(f"Logarithmic return: each halving adds $\\approx${gain:.3g}s")
    plt.savefig(os.path.join(FIGDIR, "fig_horizon_sweep.png"))
    plt.close()
    print("wrote fig_horizon_sweep.png")


def fig_de440_error(results):
    """Fig 4 (T3): positional error vs DE440 per body (common-mode clustering)."""
    errs = results["T3_positional_error"]["errors_m"]
    names = list(errs.keys())
    vals_km = [errs[n] / 1e3 for n in names]
    plt.figure(figsize=(6, 4))
    bars = plt.bar(names, vals_km, color="C0")
    bars[int(np.argmax(vals_km))].set_color("C3")
    plt.ylabel("Positional error vs DE440 (km)")
    plt.title(f"Error after {results['T3_positional_error']['years']:.0f} yr "
              f"(inner bodies cluster $\\Rightarrow$ common-mode drift)")
    plt.xticks(rotation=30, ha="right")
    plt.savefig(os.path.join(FIGDIR, "fig_de440_error.png"))
    plt.close()
    print("wrote fig_de440_error.png")


if __name__ == "__main__":
    with open(os.path.join(os.path.dirname(__file__), os.pardir, "results.json")) as f:
        results = json.load(f)
    fig_energy_behaviour()
    fig_precession()
    fig_horizon_sweep(results)
    fig_de440_error(results)
    print("all figures written to", os.path.abspath(FIGDIR))
