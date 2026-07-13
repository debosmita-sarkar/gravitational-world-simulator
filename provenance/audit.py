# provenance/audit.py
"""Enforces spec §12: every physics constant must be sourced.

Walks a module's public float attributes and returns the names of any that
lack a matching row in the module's PROVENANCE dict. Empty list == clean.
"""


def audit_constants(module) -> list:
    prov = getattr(module, "PROVENANCE", {})
    missing = []
    for name in dir(module):
        if name.startswith("_"):
            continue
        val = getattr(module, name)
        if isinstance(val, float) and name not in prov:
            missing.append(name)
    return sorted(missing)
