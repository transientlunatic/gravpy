"""
Bugs found while auditing the package.

Each is a strict xfail: once the underlying bug is fixed the test starts passing,
pytest reports XPASS(strict) as a failure, and the marker must be removed.
"""
import astropy.units as u
import numpy as np
import pytest

import gravpy.interferometers as ifo
import gravpy.sources as src

bug = lambda reason: pytest.mark.xfail(reason=reason, strict=True)  # noqa: E731


@bug("psd() has an inverted `isinstance(frequencies, NoneType)` check")
def test_psd_without_arguments_uses_default_frequencies():
    det = ifo.AdvancedLIGO()
    assert det.psd().shape == det.frequencies.shape


@bug("psd() overwrites the requested frequencies with self.frequencies")
def test_psd_honours_requested_frequencies():
    det = ifo.AdvancedLIGO()
    f = np.array([100.0, 200.0]) * u.hertz
    assert det.psd(f).shape == (2,)


@bug("`if not frequencies:` on an array/Quantity raises (also in srpsd, energy_density)")
def test_detector_helpers_accept_explicit_frequencies():
    det = ifo.AdvancedLIGO()
    f = np.array([100.0, 200.0]) * u.hertz
    assert det.srpsd(f).shape == (2,)
    assert det.noise_amplitude(f).shape == (2,)


@bug("rot_y uses undefined name `phi` instead of `psi`")
def test_rot_y_works():
    r = ifo.rot_y(0.3)
    np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)


@bug("`if r:` on an astropy Quantity raises ValueError on modern astropy")
def test_cbc_can_be_constructed():
    src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=400 * u.Mpc)


@bug("CBC cannot be constructed (see above), so snr is unreachable")
def test_cbc_snr_is_positive_and_finite():
    cbc = src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=400 * u.Mpc)
    snr = cbc.snr(ifo.AdvancedLIGO())
    assert np.isfinite(snr) and snr > 0


@bug("Source.raw_strain does `if not frequencies:` on an array, so general.snr cannot run")
def test_general_snr_runs_for_generic_source():
    s = src.Source()
    s.chirp_mass = lambda: 30 * u.solMass
    assert np.isfinite(src.general.snr(s, ifo.AdvancedLIGO()))


@bug("general.snr calls np.trapz, which NumPy >= 2.4 no longer provides (use np.trapezoid)")
def test_general_snr_does_not_use_removed_numpy_api():
    import inspect

    assert "np.trapz(" not in inspect.getsource(src.general.snr)


@bug("plotting.labelLine uses undefined names (degrees, atan2, path_effects, lato_small)")
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


@bug("timingarray.py is Python-2 only: `import data.atnf`, functools32, implicit relative import")
def test_timingarray_imports():
    import importlib

    importlib.import_module("gravpy.timingarray")


@bug("EinsteinTelescope is defined twice; the `ET` alias (bound between them) points at the first, "
     "shadowed class, which has the 'ET-D-Sum' configuration, so ET and EinsteinTelescope differ")
def test_et_alias_is_einstein_telescope():
    assert ifo.ET is ifo.EinsteinTelescope


@bug("the shadowed first ET class: `if frequencies:` / `if not frequencies:` fail on Quantity arrays")
def test_et_alias_psd_accepts_frequencies():
    f = np.array([100.0, 200.0]) * u.hertz
    assert ifo.ET().psd(f).shape == (2,)
