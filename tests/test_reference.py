"""
Checks of gravpy's detector models against independent references.

* Noise curves are compared with LALSimulation (``pip install gravpy[reference]``).
* LISA is re-implemented from the equations in arXiv:1803.01944.
* The antenna pattern is compared with the textbook expressions
  (e.g. Sathyaprakash & Schutz 2009, eq. 100; Jaranowski, Krolak & Schutz 1998).
* Finally, every curve is pinned at a few frequencies so that unintended changes to
  the numbers are caught while the code is being refactored.
"""
import astropy.units as u
import numpy as np
import pytest

import gravpy.interferometers as ifo

FREQS = np.linspace(10, 3000, 3000)  # Hz, 1 Hz resolution


def lal_psd(name, flow):
    """Evaluate a LALSimulation noise PSD of the 'REAL8FrequencySeries, flow' kind."""
    lal = pytest.importorskip("lal")
    lalsim = pytest.importorskip("lalsimulation")
    series = lal.CreateREAL8FrequencySeries(
        "psd", lal.LIGOTimeGPS(0), FREQS[0], FREQS[1] - FREQS[0], lal.HertzUnit, len(FREQS)
    )
    getattr(lalsim, name)(series, flow)
    return series.data.data


def lal_psd_scalar(name):
    lalsim = pytest.importorskip("lalsimulation")
    return np.array([getattr(lalsim, name)(f) for f in FREQS])


def gravpy_psd(det):
    det.frequencies = FREQS * u.hertz
    return det.psd(FREQS * u.hertz).to_value(1 / u.hertz)


def ratio_in_band(ours, theirs, lo, hi):
    band = (FREQS >= lo) & (FREQS <= hi)
    return ours[band] / theirs[band]


# ---------------------------------------------------------------------------
# Noise curves vs LAL
# ---------------------------------------------------------------------------

@pytest.mark.reference
def test_einstein_telescope_d_matches_lal():
    """ET-D (sum of the three interferometers) is tabulated data: should agree closely."""
    ours = gravpy_psd(ifo.EinsteinTelescope())
    theirs = lal_psd("SimNoisePSDEinsteinTelescopeP1600143", 1.0)
    r = ratio_in_band(ours, theirs, 45, 2500)
    assert np.median(np.abs(r - 1)) < 1e-4
    assert np.mean(np.abs(r - 1) < 1e-2) > 0.97


@pytest.mark.reference
def test_einstein_telescope_d_does_not_ring_near_lines():
    """A cubic spline through the tabulated ASD overshot by up to ~6x near narrow lines (~420 Hz)."""
    ours = gravpy_psd(ifo.EinsteinTelescope())
    theirs = lal_psd("SimNoisePSDEinsteinTelescopeP1600143", 1.0)
    r = ratio_in_band(ours, theirs, 45, 2500)
    # LAL interpolates differently on the line itself, so allow a factor of 2 there
    assert 0.5 < r.min() and r.max() < 2.0, (r.min(), r.max())


@pytest.mark.reference
def test_aplus_is_close_to_lal_t1800042():
    """A+ is tabulated, but from a different release of the design curve: ~25% in PSD."""
    ours = gravpy_psd(ifo.AdvancedLIGO(configuration="A+"))
    theirs = lal_psd("SimNoisePSDaLIGOAPlusDesignSensitivityT1800042", 10.0)
    r = ratio_in_band(ours, theirs, 45, 2500)
    assert 0.7 < r.min() and r.max() < 1.3, (r.min(), r.max())


@pytest.mark.reference
def test_initial_ligo_is_close_to_lal_srd():
    """Analytic fit vs the iLIGO SRD curve: agree to ~20% from 50 Hz."""
    r = ratio_in_band(gravpy_psd(ifo.InitialLIGO()), lal_psd_scalar("SimNoisePSDiLIGOSRD"), 50, 3000)
    assert 0.75 < r.min() and r.max() < 1.3, (r.min(), r.max())


@pytest.mark.reference
def test_geo_matches_lal():
    r = ratio_in_band(gravpy_psd(ifo.GEO()), lal_psd_scalar("SimNoisePSDGEO"), 50, 3000)
    np.testing.assert_allclose(r, 1.0, rtol=0.05)


# TAMA, Virgo and the analytic aLIGO fit are *not* compared with LAL: gravpy implements the
# Sathyaprakash & Schutz (2009) Table 1 fits (checked below), while LAL's TAMA is 10x higher,
# its Virgo is Advanced Virgo and its aLIGO is the zero-detuning/high-power design curve.


def test_einstein_telescope_is_finite_at_10hz():
    f = np.array([10.0]) * u.hertz
    det = ifo.EinsteinTelescope()
    det.frequencies = f
    assert np.isfinite(det.psd(f)[0])


# ---------------------------------------------------------------------------
# Analytic fits vs Table 1 of Sathyaprakash & Schutz, arXiv:0903.0338 (typed in independently)
# ---------------------------------------------------------------------------

TABLE1 = {
    # class: (fs/Hz, f0/Hz, S0/Hz^-1, S(x)/S0)
    "GEO": (40, 150, 1.0e-46, lambda x: (3.4*x)**-30 + 34/x + 20*(1 - x**2 + 0.5*x**4)/(1 + 0.5*x**2)),
    "InitialLIGO": (40, 150, 9.0e-46, lambda x: (4.49*x)**-56 + 0.16*x**-4.52 + 0.52 + 0.32*x**2),
    "TAMA": (75, 400, 7.5e-46, lambda x: x**-5 + 13/x + 9*(1 + x**2)),
    "Virgo": (20, 500, 3.2e-46, lambda x: (7.8*x)**-5 + 2/x + 0.63 + x**2),
    "AdvancedLIGO": (20, 215, 1.0e-49, lambda x: x**-4.14 - 5*x**-2 + 111*(1 - x**2 + 0.5*x**4)/(1 + 0.5*x**2)),
}


@pytest.mark.parametrize("name", list(TABLE1))
def test_analytic_fit_matches_table1(name):
    fs, f0, s0, shape = TABLE1[name]
    f = np.logspace(np.log10(fs), 3.5, 100) * 1.0000001
    det = getattr(ifo, name)()
    det.frequencies = f * u.hertz
    np.testing.assert_allclose(det.psd(f * u.hertz).to_value(1 / u.hertz), s0 * shape(f / f0), rtol=1e-10)
    assert det.fs.to_value(u.hertz) == fs and det.f0.to_value(u.hertz) == f0


# ---------------------------------------------------------------------------
# LISA vs arXiv:1803.01944 (plain numpy, no units)
# ---------------------------------------------------------------------------

def lisa_robson_2019(f, L=2.5e9, fstar=19.09e-3):
    p_oms = (1.5e-11) ** 2 * (1 + (2e-3 / f) ** 4)
    p_acc = (3e-15) ** 2 * (1 + (0.4e-3 / f) ** 2) * (1 + (f / 8e-3) ** 4)
    return 10 / (3 * L**2) * (p_oms + 4 * p_acc / (2 * np.pi * f) ** 4) * (1 + 0.6 * (f / fstar) ** 2)


def test_lisa_instrument_noise_matches_paper():
    """Subtract the galactic confusion noise, which is not part of the instrument curve."""
    f = np.logspace(-4, -1, 200) * u.hertz
    lisa = ifo.LISA()
    total = lisa.psd(f).to_value(1 / u.hertz)
    confusion = lisa.confusion_noise(f).to_value(1 / u.hertz)
    # gravpy uses f* = 19.08 mHz (paper: 19.09 mHz): a 0.05% difference in the roll-up term
    np.testing.assert_allclose(total - confusion, lisa_robson_2019(f.value, fstar=19.08e-3), rtol=1e-6)
    np.testing.assert_allclose(total - confusion, lisa_robson_2019(f.value), rtol=5e-3)


def test_lisa_best_sensitivity_is_in_the_millihertz_band():
    f = np.logspace(-4, -1, 2000) * u.hertz
    lisa = ifo.LISA()
    asd = np.sqrt((lisa.psd(f) - lisa.confusion_noise(f)).to_value(1 / u.hertz))
    assert 3e-3 < f[np.argmin(asd)].value < 3e-2
    assert 1e-21 < asd.min() < 1e-19


# ---------------------------------------------------------------------------
# Antenna pattern vs textbook
# ---------------------------------------------------------------------------

def textbook_response(theta, phi, psi):
    """L-shaped detector, arms along x and y, source at polar angle theta, azimuth phi."""
    c = np.cos(theta)
    fp = 0.5 * (1 + c**2) * np.cos(2 * phi) * np.cos(2 * psi) - c * np.sin(2 * phi) * np.sin(2 * psi)
    fx = 0.5 * (1 + c**2) * np.cos(2 * phi) * np.sin(2 * psi) + c * np.sin(2 * phi) * np.cos(2 * psi)
    return fp, fx


ANGLES = [
    (0.3, 0.5, 0.1),
    (1.0, 2.0, 0.7),
    (2.2, 4.0, 1.3),
    (np.pi / 2, 0.0, 0.0),
    (0.0, 0.0, np.pi / 4),
]


@pytest.mark.parametrize("theta,phi,psi", ANGLES)
def test_antenna_pattern_total_response_is_psi_independent(theta, phi, psi):
    """True of any correct pattern; also true of today's (which ignores psi altogether)."""
    det = ifo.AdvancedLIGO()
    ref = float(det.antenna_pattern(theta, phi, 0.0)[2])
    assert float(det.antenna_pattern(theta, phi, psi)[2]) == pytest.approx(ref)


def test_antenna_pattern_is_normalised():
    det = ifo.AdvancedLIGO()
    assert all(float(det.antenna_pattern(*a)[2]) <= 1.0 + 1e-12 for a in ANGLES)


def test_antenna_pattern_depends_on_psi():
    det = ifo.AdvancedLIGO()
    fp0, fx0, _ = det.antenna_pattern(0.0, 0.0, 0.0)
    fp1, fx1, _ = det.antenna_pattern(0.0, 0.0, np.pi / 4)
    # Overhead, the +/x responses swap under psi -> psi + pi/4
    assert float(fp0) == pytest.approx(1.0) and float(fx1) == pytest.approx(1.0)
    assert float(fx0) == pytest.approx(0.0, abs=1e-12) and float(fp1) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("theta,phi,psi", ANGLES)
def test_antenna_pattern_matches_textbook(theta, phi, psi):
    fp, fx, _ = ifo.AdvancedLIGO().antenna_pattern(theta, phi, psi)
    tp, tx = textbook_response(theta, phi, psi)
    assert (float(fp), float(fx)) == pytest.approx((abs(tp), abs(tx)), abs=1e-9)


# ---------------------------------------------------------------------------
# Regression pins (PSD in 1/Hz at 30, 100, 300, 1000 Hz).  NOT evidence of correctness:
# they only record what the code returned when the test-suite was written.
# ---------------------------------------------------------------------------

PIN_F = np.array([30, 100, 300, 1000.]) * u.hertz
PINS = {
    ("AdvancedLIGO", None): [3.326472e-46, 8.151228e-48, 5.10269e-48, 2.004036e-46],
    ("AdvancedLIGO", "O1"): [9.118245e-45, 1.2492e-46, 7.073643e-47, 3.251608e-46],
    ("AdvancedLIGO", "A+"): [3.233222e-47, 3.483091e-48, 2.385185e-48, 7.386551e-48],
    ("InitialLIGO", None): [np.nan, 1.496109e-45, 1.626276e-45, 1.326803e-44],
    ("Virgo", None): [2.512289e-44, 3.449036e-45, 1.383609e-45, 1.8016e-45],
    ("GEO", None): [np.nan, 6.170707e-45, 5.033333e-45, 8.182951e-44],
    ("TAMA", None): [np.nan, 8.141719e-43, 2.670737e-44, 5.284518e-44],
    ("EinsteinTelescope", None): [8.205109e-49, 1.546914e-49, 1.024328e-49, 3.3303e-49],
}


@pytest.mark.parametrize("name,config", list(PINS))
def test_ground_psd_regression(name, config):
    det = getattr(ifo, name)(configuration=config) if config else getattr(ifo, name)()
    det.frequencies = PIN_F
    np.testing.assert_allclose(det.psd(PIN_F).to_value(1 / u.hertz), PINS[(name, config)], rtol=1e-5)


SPACE_F = np.array([1e-3, 1e-2, 1e-1]) * u.hertz
SPACE_PINS = {
    "LISA": [1.444712e-37, 1.449476e-40, 2.150351e-39],
    "EvolvedLISA": [1.694643e-35, 4.752545e-39, 1.121569e-38],
    "BDecigo": [6.39936e-36, 6.399364e-40, 6.439761e-44],
    "Decigo": [5.333e-39, 5.333062e-43, 6.037244e-47],
    "BigBangObservatory": [1.26e-39, 1.260005e-43, 1.306e-47],
}


@pytest.mark.parametrize("name", list(SPACE_PINS))
def test_space_psd_regression(name):
    psd = getattr(ifo, name)().psd(SPACE_F).to_value(1 / u.hertz)
    np.testing.assert_allclose(psd, SPACE_PINS[name], rtol=1e-5)
