# tests/test_planet_constants.py
from physics import constants as k
from provenance.audit import audit_constants


def test_gm_bodies_present_and_sourced():
    expected = {"sun", "mercury", "venus", "emb", "mars",
                "jupiter", "saturn", "uranus", "neptune"}
    assert expected <= set(k.GM_BODIES)
    # Jupiter is the largest planet GM; sanity magnitude check.
    assert 1e17 < k.GM_BODIES["jupiter"] < 1.5e17
    # Sun dominates.
    assert k.GM_BODIES["sun"] > 1e20


def test_audit_still_passes():
    assert audit_constants(k) == []
