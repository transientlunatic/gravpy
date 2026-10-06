import numpy as np
import astropy.units as u
def snr(signal, detector):
    """
    Calculate the optimal (matched-filter) SNR of a signal in a given detector,

    .. math::
       \\rho^2 = 4 \\int \\frac{|\\tilde{h}(f)|^2}{S_n(f)} \\mathrm{d}f ,

    where :math:`\\tilde{h}` is the Fourier transform of the signal as returned by
    ``signal.raw_strain`` and :math:`S_n` is the one-sided PSD of the detector.  The
    time the signal spends at each frequency is already contained in
    :math:`|\\tilde{h}(f)|`, so no additional "number of cycles" factor is needed.
    Frequencies where either the signal or the PSD is undefined (NaN) do not contribute.

    The signal is assumed to be optimally oriented and located (and so is the
    antenna response); see arxiv.org/abs/1408.0740.

    Parameters
    ----------
    signal : Source
        A Source object which describes the source producing the
        signal, e.g. a CBC.

    detector : Detector
        A Detector object describing the instrument making the observation
        e.g. aLIGO.

    Returns
    -------
    SNR : float
        The signal-to-noise ratio of the signal in the detector.
    """
    frequencies = detector.frequencies
    noise = detector.psd(frequencies)
    ampli = signal.raw_strain(frequencies)
    # |h|^2 has units of 1/Hz^2 and the PSD 1/Hz, so the integrand is in 1/Hz
    integrand = (4 * np.abs(ampli)**2 / noise).to_value(1 / u.hertz)
    integrand = np.nan_to_num(integrand, nan=0.0)
    return float(np.sqrt(np.trapezoid(integrand, x=frequencies.to_value(u.hertz))))
