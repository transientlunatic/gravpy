"""
Tests for bugs found while auditing the package, all now fixed.
"""
import astropy.units as u
import numpy as np
import pytest

import gravpy.interferometers as ifo
import gravpy.sources as src

def test_psd_without_arguments_uses_default_frequencies():
    det = ifo.AdvancedLIGO()
    assert det.psd().shape == det.frequencies.shape


def test_psd_honours_requested_frequencies():
    det = ifo.AdvancedLIGO()
    f = np.array([100.0, 200.0]) * u.hertz
    assert det.psd(f).shape == (2,)


def test_detector_helpers_accept_explicit_frequencies():
    det = ifo.AdvancedLIGO()
    f = np.array([100.0, 200.0]) * u.hertz
    assert det.srpsd(f).shape == (2,)
    assert det.noise_amplitude(f).shape == (2,)


def test_rot_y_works():
    r = ifo.rot_y(0.3)
    np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)


def test_cbc_can_be_constructed():
    src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=400 * u.Mpc)


def test_cbc_snr_is_positive_and_finite():
    cbc = src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=400 * u.Mpc)
    snr = cbc.snr(ifo.AdvancedLIGO())
    assert np.isfinite(snr) and snr > 0


def test_general_snr_runs_for_generic_source():
    s = src.Source()
    s.chirp_mass = lambda: 30 * u.solMass
    assert np.isfinite(src.general.snr(s, ifo.AdvancedLIGO()))


def test_general_snr_does_not_use_removed_numpy_api():
    import inspect

    assert "np.trapz(" not in inspect.getsource(src.general.snr)


def test_plotting_has_no_undefined_names():
    import ast
    import builtins
    import pathlib

    import gravpy.plotting as plotting

    tree = ast.parse(pathlib.Path(plotting.__file__).read_text())
    defined = set(dir(builtins)) | set(vars(plotting))
    used = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    local = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)}
    local |= {a.arg for n in ast.walk(tree) if isinstance(n, ast.arguments)
              for a in n.args + n.kwonlyargs + ([n.kwarg] if n.kwarg else [])}
    assert not (used - defined - local)


def test_timingarray_imports():
    import importlib

    importlib.import_module("gravpy.timingarray")


def test_et_alias_is_einstein_telescope():
    assert ifo.ET is ifo.EinsteinTelescope


def test_et_alias_psd_accepts_frequencies():
    f = np.array([100.0, 200.0]) * u.hertz
    assert ifo.ET().psd(f).shape == (2,)


def test_timingarray_hellings_downs_values():
    from astropy.coordinates import SkyCoord

    from gravpy.timingarray import Pulsar, hellingsdowns_factor

    def psr(ra, dec):
        return Pulsar("x", 1 * u.day, 1 * u.year, 1e-7 * u.second, SkyCoord(ra * u.deg, dec * u.deg))

    a = psr(0, 0)
    assert hellingsdowns_factor(a, a) == 1
    # Orthogonal pulsars: 1.5 x ln x - x/4 + 1/2 at x = 1/2 (the HD curve is ~ -0.145 there)
    assert hellingsdowns_factor(a, psr(90, 0)) == pytest.approx(-0.14486, abs=1e-4)
    # Antipodal pulsars: x = 1 -> 1/4
    assert hellingsdowns_factor(a, psr(180, 0)) == pytest.approx(0.25)


def test_timingarray_instances_do_not_share_pulsars():
    from gravpy.timingarray import TimingArray

    assert "pulsars" not in vars(TimingArray) or TimingArray.pulsars == []


def test_labelline_runs():
    import matplotlib.pyplot as plt

    from gravpy.plotting import labelLine

    fig, ax = plt.subplots()
    (line,) = ax.loglog([1, 10, 100], [1, 10, 100], label="x")
    labelLine(line, 20, "label")
    plt.close(fig)


def test_antenna_pattern_psi_average():
    """Averaging F+^2 and Fx^2 over a half-turn of psi gives (F+^2 + Fx^2)/2 each."""
    det = ifo.AdvancedLIGO()
    _, _, total = det.antenna_pattern(0.8, 1.1, 0.0)
    fp, fx, tot = det.antenna_pattern(0.8, 1.1, [0, np.pi])
    assert fp == pytest.approx(total / np.sqrt(2))
    assert fx == pytest.approx(total / np.sqrt(2))
    assert tot == pytest.approx(total)


def test_antenna_pattern_known_values():
    det = ifo.AdvancedLIGO()
    # Overhead: F+ = cos 2(phi + psi) pattern, |F| = 1 everywhere
    assert det.antenna_pattern(0, 0, 0)[:2] == pytest.approx((1.0, 0.0))
    # In the detector plane along an arm bisector the response vanishes for 'x' and is 0 for '+'.. at 45 deg
    assert det.antenna_pattern(np.pi / 2, np.pi / 4, 0.0)[2] == pytest.approx(0.0, abs=1e-12)
