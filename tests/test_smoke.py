"""
Characterisation tests: what works *today*.

Anything that does not work is recorded in test_known_bugs.py as a strict xfail,
so that fixing it forces the marker to be removed.
"""
import importlib

import astropy.units as u
import numpy as np
import pytest

import gravpy.interferometers as ifo

FREQS = np.logspace(1, 3.5, 50) * u.hertz

# (class name, constructor kwargs)
DETECTORS = [
    ("AdvancedLIGO", {}),
    ("AdvancedLIGO", {"configuration": "O1"}),
    ("AdvancedLIGO", {"configuration": "A+"}),
    ("InitialLIGO", {}),
    ("Virgo", {}),
    ("GEO", {}),
    ("TAMA", {}),
    ("EinsteinTelescope", {}),
]
SPACE_DETECTORS = ["LISA", "EvolvedLISA", "BDecigo", "Decigo", "BigBangObservatory"]


@pytest.mark.parametrize(
    "module",
    ["gravpy", "gravpy.gravpy", "gravpy.general", "gravpy.noise", "gravpy.plotting",
     "gravpy.interferometers", "gravpy.sources"],
)
def test_module_imports(module):
    importlib.import_module(module)


@pytest.mark.parametrize("name,kwargs", DETECTORS)
def test_ground_detector_psd(name, kwargs):
    det = getattr(ifo, name)(frequencies=FREQS, **kwargs)
    psd = det.psd(FREQS)
    assert psd.unit.is_equivalent(1 / u.hertz)
    good = psd[np.isfinite(psd)]
    assert len(good) > 0
    assert np.all(good.value > 0)
    # Every ground-based curve lies between 1e-52 and 1e-38 /Hz over 10 Hz - 3 kHz
    assert np.all((good.value > 1e-52) & (good.value < 1e-38))


@pytest.mark.parametrize("name", SPACE_DETECTORS)
def test_space_detector_psd(name):
    cls = getattr(ifo, name)
    psd = cls().psd(cls.frequencies)
    assert psd.shape == cls.frequencies.shape
    assert np.all(psd.value[np.isfinite(psd.value)] > 0)


def test_psd_is_nan_below_low_frequency_cutoff():
    det = ifo.AdvancedLIGO(frequencies=FREQS)
    psd = det.psd(FREQS)
    assert np.all(np.isnan(psd[FREQS < det.fs]))
    assert np.all(np.isfinite(psd[FREQS >= det.fs]))


def test_aligo_has_minimum_near_f0():
    det = ifo.AdvancedLIGO()
    f = np.linspace(30, 4000, 4000) * u.hertz
    fmin = f[np.argmin(det.psd(f))]
    assert 150 * u.hertz < fmin < 400 * u.hertz


def test_noise_amplitude_and_srpsd_consistent():
    det = ifo.AdvancedLIGO()
    # Default (None) frequencies only: passing an array is broken, see test_known_bugs.py
    ok = np.isfinite(det.psd(det.frequencies))
    np.testing.assert_allclose(det.srpsd()[ok].value ** 2, det.psd(det.frequencies)[ok].value)
    np.testing.assert_allclose(
        det.noise_amplitude()[ok].value ** 2,
        (det.frequencies * det.psd(det.frequencies))[ok].value,
    )


def test_plot_returns_line():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    lines = ifo.AdvancedLIGO().plot(ax)
    assert len(lines) == 1
    plt.close(fig)


def test_configuration_changes_name():
    assert "A+" in ifo.AdvancedLIGO(configuration="A+").name


def test_antenna_pattern_runs_and_is_bounded():
    det = ifo.AdvancedLIGO()
    fp, fx, f = det.antenna_pattern(0.3, 0.5, 0.1)
    # Characterisation only: correctness is checked in test_reference.py
    assert np.isfinite(float(fp)) and np.isfinite(float(fx))
    assert float(f) == pytest.approx(np.hypot(float(fp), float(fx)))


def test_rot_x_rot_z_are_rotations():
    for rot in (ifo.rot_x, ifo.rot_z):
        r = rot(0.7)
        np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(r) == pytest.approx(1.0)
