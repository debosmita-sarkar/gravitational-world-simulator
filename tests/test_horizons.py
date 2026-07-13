# tests/test_horizons.py
import numpy as np
import pytest
from physics.constants import AU

horizons = pytest.importorskip("validation.horizons")


HORIZONS_IDS = {  # our body name -> Horizons id (planet-system barycentres)
    "sun": "10", "mercury": "1", "venus": "2", "emb": "3", "mars": "4",
    "jupiter": "5", "saturn": "6", "uranus": "7", "neptune": "8",
}


@pytest.mark.network
def test_fetch_mercury_state_reasonable():
    # J2000.0 = JD 2451545.0 TDB. Mercury ~0.31-0.47 AU from barycentre.
    pos, vel = horizons.fetch_state_si("1", 2451545.0)
    r_au = np.linalg.norm(pos) / AU
    assert 0.2 < r_au < 0.6
    # orbital speed ~ 40-60 km/s
    assert 3e4 < np.linalg.norm(vel) < 7e4


@pytest.mark.network
def test_build_solar_system_has_expected_bodies():
    s = horizons.build_solar_system(2451545.0, ["sun", "emb", "jupiter"])
    assert s.positions.shape == (3, 3)
    assert s.masses[0] > s.masses[2]  # Sun heaviest
