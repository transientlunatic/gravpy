import numpy as np
import astropy.units as u
def snr(signal, detector):
    """
    Calculate the SNR of a signal in a given detector,
    assuming that it has been detected with an optimal filter.
    See e.g. arxiv.org/abs/1408.0740
    
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
    if signal.ncycles(): 
        ncycles = np.sqrt(2*signal.ncycles(detector.frequencies))
    else: 
        ncycles = 1
    noise = detector.psd(detector.frequencies)
    ampli = signal.raw_strain(detector.frequencies) * ncycles
    fraction = 4*(np.abs(ampli)**2 / noise)
    fraction[np.isnan(fraction)]=0
    return np.sqrt(np.trapz(fraction, x=detector.frequencies, dx=0.01*u.hertz))

