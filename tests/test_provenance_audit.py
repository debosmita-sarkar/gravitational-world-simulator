# tests/test_provenance_audit.py
from provenance.audit import audit_constants
from physics import constants


def test_all_physics_constants_are_sourced():
    missing = audit_constants(constants)
    assert missing == [], f"Unsourced constants: {missing}"


def test_audit_flags_an_unsourced_constant():
    import types
    fake = types.SimpleNamespace(MYSTERY=42.0, PROVENANCE={})
    assert audit_constants(fake) == ["MYSTERY"]
