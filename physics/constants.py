# physics/constants.py
"""Physical constants. Every entry is mirrored in PROVENANCE (spec §12)."""

G = 6.67430e-11          # m^3 kg^-1 s^-2
c = 299792458.0          # m s^-1
GM_SUN = 1.32712440018e20  # m^3 s^-2  (carried as a product; see spec §12)
AU = 1.495978707e11      # m

PROVENANCE = {
    "G": {"value": G, "unit": "m^3 kg^-1 s^-2", "source": "CODATA 2018", "source_type": "CODATA"},
    "c": {"value": c, "unit": "m s^-1", "source": "SI base definition (2019)", "source_type": "SI"},
    "GM_SUN": {"value": GM_SUN, "unit": "m^3 s^-2", "source": "IAU 2015 / DE440", "source_type": "IAU"},
    "AU": {"value": AU, "unit": "m", "source": "IAU 2012 Resolution B2", "source_type": "IAU"},
}

# Planet-system barycentre gravitational parameters GM (SI, m^3 s^-2).
# Values from JPL DE440 / IAU 2015; carried as GM products (measured to higher
# precision than G*M separately). "emb" = Earth-Moon barycentre.
GM_BODIES = {
    "sun":     1.32712440018e20,
    "mercury": 2.2031868551e13,
    "venus":   3.2485859200e14,
    "emb":     4.0350323562e14,
    "mars":    4.2828375816e13,
    "jupiter": 1.26712764100e17,
    "saturn":  3.7940584841e16,
    "uranus":  5.7945564000e15,
    "neptune": 6.8365271000e15,
}

for _name, _gm in GM_BODIES.items():
    PROVENANCE[f"GM_{_name}"] = {
        "value": _gm, "unit": "m^3 s^-2",
        "source": "JPL DE440 / IAU 2015", "source_type": "DE440",
    }
