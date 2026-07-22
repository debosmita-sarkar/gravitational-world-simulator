# Gravitational World Simulator

> A provenance-audited gravitational **N-body world simulator** validated against NASA JPL **DE440**, introducing the **Fidelity–Limit Theorem** — a closed-form predictability horizon for chaotic systems.

**Author:** [Debosmita Sarkar](https://debosmitasarkar.site/) — **Independent Researcher, Kolkata, India**
**Paper:** [Read the full paper &amp; results →](https://debosmitasarkar.site/)  ·  **DOI:** [10.5281/zenodo.21329246](https://doi.org/10.5281/zenodo.21329246)
**License:** MIT (code) · CC BY 4.0 (paper)

---

## Overview

This repository contains the complete, from-scratch implementation behind the paper
*"A Provenance-Audited Gravitational World Simulator Validated Against DE440 with a Formal Predictability Analysis"*
by **Debosmita Sarkar**, an Independent Researcher based in Kolkata, India.

The simulator is built on a single principle: **every governing rule is a proven law of physics with its derivation attached — no hand-tuned constants, no invented forces.** Every physical constant is bound to a cited source and enforced by an automated provenance audit that fails the build if any unsourced "magic number" appears in the physics core.

It is initialised from real solar-system states in NASA JPL's **DE440 ephemeris** and its predictions are compared against measured positions in physical units. The work then *proves*, rather than asserts, exactly where predictive fidelity must end — the **Fidelity–Limit Theorem**.

## Key results

| Result | Outcome |
|---|---|
| **Mercury perihelion advance** | **42.99″/century** vs. measured 43″ — a 0.03% match to general relativity |
| **Fidelity–Limit Theorem** | Closed-form predictability horizon, proved and numerically confirmed |
| **Logarithmic-return law (T6)** | Horizon gain per step-halving matches `p·ln2/λ` to **four significant figures** |
| **~100 Myr horizon** | Canonical solar-system values reproduce the established limit (≈115 Myr) with no fitting |
| **Symplectic vs. RK4** | Leapfrog energy error bounded (≤ 3.2×10⁻¹³); RK4 drifts secularly |
| **Provenance audit** | Every constant traced to a cited source (CODATA, SI, IAU, DE440) |

### The Fidelity–Limit Theorem

The simulated positional error stays below a tolerance Δ_tol for all *t ≤ t_lim*, where:

```
t_lim = (1/λ) · ln[ Δ_tol / (Δ₀ + C·hᵖ/λ) ]
```

binding the integrator order *p*, step size *h*, the largest Lyapunov exponent *λ*, the measurement uncertainty Δ₀, and the tolerance Δ_tol. Its central corollary — the **measurement ceiling** — states that beyond a computable point, prediction is limited by how well reality is *measured*, not by computational effort. In a chaotic world, you cannot compute your way past the measurement uncertainty of reality.

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows
# source .venv/bin/activate       # macOS / Linux
pip install -r requirements.txt
```

Requires Python 3.10+ and NumPy.

## Running the tests

The project is test-driven; the suite verifies the Kepler oracle, symplectic
energy behaviour, the 1PN relativistic term, the provenance audit, and more.

```bash
python -m pytest -v
```

## Reproducing the paper's results

```bash
python -m analysis.generate_results   # writes results.json (Horizons on first run, cached after)
python -m analysis.generate_figures   # regenerates the paper figures
```

Every number in the paper's Results section is real measured output of this code.

## Repository layout

```
physics/        Proven laws, integrators (leapfrog + RK4), Kepler oracle — pure NumPy
provenance/     Automated constants-provenance audit
analysis/       Lyapunov estimation, fidelity-limit computation, precession, figures
validation/     Real-data fetch (NASA JPL Horizons / DE440) and error computation
tests/          Test-driven-development suite
paper/          LaTeX manuscript and figures
```

## About the author

**Debosmita Sarkar** is an **Independent Researcher based in Kolkata, India**, working in gravitational dynamics, celestial mechanics, and the numerical simulation of chaotic physical systems. She is the author of this provenance-audited gravitational world simulator and of the **Fidelity–Limit Theorem**, a closed-form characterization of the predictability horizon of chaotic dynamical systems. Her work emphasizes honesty and verifiability: every constant is traced to a cited source, every result is compared against real measurement data in physical units, and the limits of prediction are proved rather than asserted.

- 🌐 Website: [debosmitasarkar.site](https://www.debosmitasarkar.site/)
- 🆔 ORCID: [0009-0006-7484-7804](https://orcid.org/0009-0006-7484-7804)
- 📄 Paper (Zenodo): [DOI 10.5281/zenodo.21329246](https://doi.org/10.5281/zenodo.21329246)
- 🎓 Google Scholar: [Debosmita Sarkar](https://scholar.google.com/citations?user=mJzjDnEAAAAJ)
- ✉️ Email: debosmitasarkar993@gmail.com

## Citation

If you use this work, please cite it (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@misc{sarkar2026simulator,
  author       = {Sarkar, Debosmita},
  title        = {A Provenance-Audited Gravitational World Simulator
                  Validated Against DE440 with a Formal Predictability Analysis},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21329246},
  url          = {https://doi.org/10.5281/zenodo.21329246},
  note         = {Independent Researcher, Kolkata, India}
}
```

## License

- **Code** — [MIT License](LICENSE)
- **Paper &amp; figures** — [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

© 2026 Debosmita Sarkar, Independent Researcher, Kolkata, India.
