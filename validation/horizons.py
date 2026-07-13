# validation/horizons.py
"""Fetch real barycentric ICRF/TDB states from JPL Horizons, cached to disk.

Horizons vector table returns positions in AU and velocities in AU/day; we
convert to SI on ingest. location='@0' is the Solar-System barycentre.
"""
import json
import os
import numpy as np
from astroquery.jplhorizons import Horizons
from physics.constants import AU, G, GM_BODIES
from physics.state import State

_DAY = 86400.0
_CACHE = os.path.join(os.path.dirname(__file__), "cache")

HORIZONS_IDS = {
    "sun": "10", "mercury": "1", "venus": "2", "emb": "3", "mars": "4",
    "jupiter": "5", "saturn": "6", "uranus": "7", "neptune": "8",
}


def _cache_path(body_id: str, jd_tdb: float) -> str:
    os.makedirs(_CACHE, exist_ok=True)
    return os.path.join(_CACHE, f"{body_id}_{jd_tdb:.6f}.json")


def fetch_state_si(body_id: str, jd_tdb: float):
    """Return (position_m (3,), velocity_m_s (3,)) barycentric ICRF."""
    path = _cache_path(body_id, jd_tdb)
    if os.path.exists(path):
        with open(path) as f:
            d = json.load(f)
        return np.array(d["pos"]), np.array(d["vel"])
    obj = Horizons(id=body_id, location="@0", epochs=jd_tdb, id_type=None)
    v = obj.vectors(refplane="earth")
    pos = np.array([float(v["x"][0]), float(v["y"][0]), float(v["z"][0])]) * AU
    vel = np.array([float(v["vx"][0]), float(v["vy"][0]), float(v["vz"][0])]) * AU / _DAY
    with open(path, "w") as f:
        json.dump({"pos": pos.tolist(), "vel": vel.tolist()}, f)
    return pos, vel


def build_solar_system(jd_tdb: float, bodies) -> State:
    masses, positions, velocities = [], [], []
    for name in bodies:
        pos, vel = fetch_state_si(HORIZONS_IDS[name], jd_tdb)
        masses.append(GM_BODIES[name] / G)
        positions.append(pos)
        velocities.append(vel)
    return State(np.array(masses), np.array(positions), np.array(velocities))
