=======
History
=======


Unreleased
----------

* Requires Python >= 3.10 and NumPy >= 2.0; packaging moved to ``pyproject.toml``.
* Fixed passing explicit ``frequencies`` to ``psd``, ``srpsd``, ``noise_amplitude`` and ``energy_density``.
* Fixed ``antenna_pattern``: responses were 2x too large and ignored the polarisation angle;
  ``psi=[lo, hi]`` now gives the RMS over that range.
* Fixed the GEO600 noise curve (``0.4 x^4`` should be ``0.5 x^4``, per Sathyaprakash & Schutz 2009).
* ``EinsteinTelescope`` is now defined once (``ET`` is an alias), is valid down to 1 Hz,
  and the tabulated curves are interpolated in log-log space (PCHIP) instead of with a spline that rang near lines.
* Fixed ``general.snr`` and ``CBC``: the optimal SNR is now ``sqrt(4 int |h|^2/S df)`` (an extra ``sqrt(2 n_cycles)``
  factor double-counted the time spent at each frequency, and the waveform was cut off at twice the ISCO frequency);
  ``CBC.characteristic_strain`` is now ``2 f |h|`` like every other source; ``CBC.ncycles`` is ``f^2/fdot`` and dimensionless.
* Fixed constructors of ``Source`` and its subclasses with astropy Quantities, ``rot_y``,
  ``gravpy.plotting`` and ``gravpy.timingarray`` (now Python 3).

0.2.1 (2019-07-16)
------------------

* Added the noise curve for the Einstein Telescope (configuration D only).

0.1.0 (2016-03-15)
------------------

* First release on PyPI.
