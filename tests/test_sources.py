"""
Physics checks of the compact-binary source model and the optimal-SNR calculation.
"""
import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest

import gravpy.general as general
import gravpy.interferometers as ifo
import gravpy.sources as src


@pytest.fixture
def cbc():
    return src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=400 * u.Mpc)


def test_chirp_mass(cbc):
    # Mc = (m1 m2)^(3/5) / (m1 + m2)^(1/5) = 30 * 2^(-1/5) Msun for equal masses
    assert cbc.chirp_mass().to_value(u.solMass) == pytest.approx(30 * 2 ** -0.2, rel=1e-9)


def test_fisco_is_gravitational_wave_frequency(cbc):
    # f_gw,ISCO = c^3 / (6^(3/2) pi G M) = 4397 Hz (Msun / M)
    assert cbc.fisco().to_value(u.hertz) == pytest.approx(4396.8 / 60, rel=1e-3)


def test_strain_stops_at_isco(cbc):
    f = np.array([0.5, 0.99, 1.01, 2.0]) * cbc.fisco()
    h = cbc.raw_strain(f)
    assert np.all(np.isfinite(h[:2])) and np.all(np.isnan(h[2:]))


def test_strain_scales_as_f_to_minus_7_6_and_inverse_distance(cbc):
    f = np.array([20.0, 40.0]) * u.hertz
    h = cbc.raw_strain(f)
    assert (h[1] / h[0]).to_value(u.dimensionless_unscaled) == pytest.approx(2 ** (-7 / 6))
    far = src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=800 * u.Mpc)
    assert (cbc.raw_strain(f) / far.raw_strain(f)).to_value(u.dimensionless_unscaled) == pytest.approx(2.0)


@pytest.mark.reference
def test_strain_matches_lal_taylorf2(cbc):
    """|h(f)| of an optimally oriented (face-on, overhead) inspiral."""
    lal = pytest.importorskip("lal")
    lalsim = pytest.importorskip("lalsimulation")
    df = 0.25
    hp, _ = lalsim.SimInspiralChooseFDWaveform(
        30 * lal.MSUN_SI, 30 * lal.MSUN_SI, 0, 0, 0, 0, 0, 0, 400e6 * lal.PC_SI,
        0.0, 0.0, 0.0, 0.0, 0.0, df, 10.0, 1024.0, 10.0, None, lalsim.TaylorF2,
    )
    f = np.arange(len(hp.data.data)) * df
    band = (f >= 20) & (f <= 70)
    ours = cbc.raw_strain(f[band] * u.hertz).to_value(1 / u.hertz)
    np.testing.assert_allclose(ours, np.abs(hp.data.data[band]), rtol=1e-6)


@pytest.mark.reference
def test_snr_matches_independent_integral_with_lal_waveform(cbc):
    """rho^2 = 4 int |h|^2 / S df, with the LAL amplitude and gravpy's PSD, up to the ISCO."""
    lal = pytest.importorskip("lal")
    lalsim = pytest.importorskip("lalsimulation")
    det = ifo.AdvancedLIGO()
    f = det.frequencies
    df = 0.25
    hp, _ = lalsim.SimInspiralChooseFDWaveform(
        30 * lal.MSUN_SI, 30 * lal.MSUN_SI, 0, 0, 0, 0, 0, 0, 400e6 * lal.PC_SI,
        0.0, 0.0, 0.0, 0.0, 0.0, df, 10.0, 1024.0, 10.0, None, lalsim.TaylorF2,
    )
    grid = np.arange(len(hp.data.data)) * df
    h = np.interp(f.value, grid, np.abs(hp.data.data))
    band = (f <= cbc.fisco()) & np.isfinite(det.psd(f))
    rho2 = 4 * np.trapezoid(h[band] ** 2 / det.psd(f)[band].to_value(1 / u.hertz), f[band].value)
    assert general.snr(cbc, det) == pytest.approx(np.sqrt(rho2), rel=0.02)


def test_snr_scales_inversely_with_distance(cbc):
    det = ifo.AdvancedLIGO()
    far = src.CBC(m1=30 * u.solMass, m2=30 * u.solMass, r=800 * u.Mpc)
    assert general.snr(cbc, det) == pytest.approx(2 * general.snr(far, det), rel=1e-9)


def test_snr_is_a_float_in_a_physically_sensible_range(cbc):
    snr = general.snr(cbc, ifo.AdvancedLIGO())
    assert isinstance(snr, float)
    # Optimal (face-on, overhead) inspiral-only SNR of 30+30 Msun at 400 Mpc in aLIGO design:
    # ~85 using LAL's ZDHP curve from 10 Hz; gravpy's fit only starts at 30 Hz so is lower.
    assert 40 < snr < 120


def test_snr_equals_integral_of_characteristic_strain_over_ln_f(cbc):
    """rho^2 = int (h_c / h_n)^2 dln f with h_c = 2 f |h| and h_n^2 = f S_n."""
    det = ifo.AdvancedLIGO()
    f = det.frequencies
    hc = cbc.characteristic_strain(f).to_value(u.dimensionless_unscaled)
    hn2 = (f * det.psd(f)).to_value(u.dimensionless_unscaled)
    integrand = np.nan_to_num(hc ** 2 / hn2, nan=0.0)
    rho2 = np.trapezoid(integrand, np.log(f.value))
    assert general.snr(cbc, det) == pytest.approx(np.sqrt(rho2), rel=2e-3)


def test_ncycles_integrates_to_the_textbook_cycle_count(cbc):
    f1, f2 = 20 * u.hertz, 60 * u.hertz
    f = np.geomspace(f1.value, f2.value, 4000) * u.hertz
    n = cbc.ncycles(f).to_value(u.dimensionless_unscaled)
    numeric = np.trapezoid(n, np.log(f.value))
    tc = (const.G * cbc.chirp_mass() / const.c ** 3).to(u.second)
    exact = (1 / (32 * np.pi ** (8 / 3)) * tc ** (-5 / 3) * (f1 ** (-5 / 3) - f2 ** (-5 / 3))).to_value(
        u.dimensionless_unscaled
    )
    assert numeric == pytest.approx(exact, rel=1e-4)
