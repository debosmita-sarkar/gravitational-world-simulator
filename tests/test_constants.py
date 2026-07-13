# tests/test_constants.py
from physics import constants as k


def test_constant_values_match_spec():
    assert k.G == 6.67430e-11
    assert k.c == 299792458.0
    assert k.GM_SUN == 1.32712440018e20
    assert k.AU == 1.495978707e11


def test_every_constant_has_a_provenance_row():
    for symbol in ("G", "c", "GM_SUN", "AU"):
        row = k.PROVENANCE[symbol]
        assert row["value"] == getattr(k, symbol)
        assert row["unit"]
        assert row["source"]
        assert row["source_type"] in {"CODATA", "SI", "IAU", "DE440", "derived"}
