# tests/test_positional_error.py
import pytest
horizons = pytest.importorskip("validation.horizons")
pe = pytest.importorskip("validation.positional_error")


@pytest.mark.network
def test_short_interval_error_is_bounded():
    # Over 1 year, Sun+inner planets, error should be far below the orbit scale.
    out = pe.validate_against_de440(
        bodies=["sun", "mercury", "venus", "emb", "mars"],
        jd_start=2451545.0, years=1.0, steps_per_year=8760, use_gr=True)
    # Sanity: max positional error is well under 1 AU (1.496e11 m).
    assert out["max_error_m"] < 1e11
    assert out["max_error_m"] > 0.0  # it's a real measurement, not identically zero
