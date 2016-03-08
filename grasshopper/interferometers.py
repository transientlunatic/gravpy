import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
import numpy as np

class Detector():
    """
    This is the base class for all types of detectors, and 
    contains the conversion methods between the various 
    different ways of expressing the noise levels (sensitivity)
    of any detector.
    """
    def noise_amplitude(self, frequencies=None):
        """
        The noise amplitude for a detector is defined as
        $h^2_n(f) = f S_n(f)$
        and is designed to incorporate the effect of integrating 
        an inspiralling signal.
        
        Parameters
        ----------
        frequencies : ndarray
            An array of frequencies, in units of Hz
            
        Returns
        -------
        noise_amplitude : ndarray
            An array of the noise amplitudes correcsponding 
            to the input frequency values
        """
        if not frequencies: frequencies = self.frequencies
        return np.sqrt(self.frequencies*self.psd(frequencies))
    
    def energy_density(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        bigH = (2*np.pi**2)/3 * frequencies**3 * self.psd(frequencies)
        littleh = bigH / ((100*u.kilometer / u.second / u.megaparsec).to(u.hertz))**2
        return littleh
    
    def srpsd(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return np.sqrt(self.psd(frequencies))
    
    def plot(self, axis=None):
        """
        Plot the noise curve for this detector.
        """
        if axis: 
            axis.loglog(self.frequencies, self.noise_amplitude(self.frequencies), label=self.name, lw=2)
            axis.set_xlabel('Frequency [Hz]')
            axis.set_ylabel('Characteristic Strain')
            axis.legend()

class Interferometer(Detector):
    """
    The base class to describe an interferometer.
    """
    name = "Generic Interferometer"
    f0 = 150 * u.hertz
    fs = 40 * u.hertz
    S0 = 1e-46 / u.hertz
    frequencies =  np.logspace(1, 5, 10000) * u.hertz
    
    def __init__(self, frequencies=None):
        if frequencies: self.frequencies = frequencies
    
    def noise_spectrum(self, x):
        return (3.4*x)**(-30) + 34*x**(-1) + (20 * (1 - x**2 + 0.4*x**4))/(1 + 0.5*x**2)
    
    def psd(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        x = frequencies / self.f0
        xs = self.fs / self.f0
        sh = self.noise_spectrum(x)
        sh[frequencies<self.fs]=np.nan
        return sh * self.S0

    
    

class TimingArray(Detector):
    """
    A class to represent a pulsar timing array.
    """
    name = "Generic PTA"
    dt = 14*u.day  # the sampling interval
    T  = 15*u.year # the observation time
    sigma = 100 * u.nanosecond # the timing uncertainty of each observation

    frequencies =  np.logspace(-10, -6, 1000) * u.hertz
    n = 20
    zeta_sum = 4.74
    
    def Pn(self, frequencies):
        dt = self.dt.to(u.second)
        sigma = self.sigma.to(u.second)
        return 2 * dt * sigma**2
        
    def Sn(self, frequencies):
        return 12 * np.pi**2 * frequencies**2 * self.Pn(frequencies)
    
    def noise_spectrum(self, frequencies):
        return self.Sn(frequencies)*self.zeta_sum**(-0.5)
    
    def psd(self, frequencies):
        # We're currently over-estimating the sensitivity, 
        # we can get around this using 
        # http://iopscience.iop.org/article/10.1088/0264-9381/30/22/224015/pdf
        lower = 1 / self.T
        upper = 1 / self.dt
        sh = self.noise_spectrum(frequencies)
        sh[frequencies<lower]=np.nan
        sh[frequencies>upper]=np.nan
        return sh 

    
            
class AdvancedLIGO(Interferometer):
    """
    The aLIGO Interferometer
    """
    name = "aLIGO"
    f0 = 215 * u.hertz
    fs = 20 * u.hertz
    S0 = 1.0e-49 / u.hertz
    
    def noise_spectrum(self, x):
        return (x)**(-4.14) -5*x**(-2) + ((111 * (1-x**2 +0.5*x**4))/(1+0.5*x**2))
    
class GEO(Interferometer):
    """
    The GEO600 Interferometer
    """
    name = "GEO600"
    f0 = 150 * u.hertz
    fs = 40 * u.hertz
    S0 = 1e-46 / u.hertz
    
    def noise_spectrum(self, x):
        return (3.4*x)**(-30) + 34*x**(-1) + (20 * (1 - x**2 + 0.4*x**4))/(1 + 0.5*x**2)
    
class InitialLIGO(Interferometer):
    """
    The iLIGO Interferometer
    """
    name = "Initial LIGO"
    f0 = 150 * u.hertz
    fs = 40 * u.hertz
    S0 = 9e-46 / u.hertz
    
    def noise_spectrum(self, x):
        return (4.49*x)**(-56) + 0.16*x**(-4.52) + 0.52 + 0.32*x**2
    
class TAMA(Interferometer):
    """
    The TAMA Interferometer
    """
    name = "TAMA"
    f0 = 400 * u.hertz
    fs = 75 * u.hertz
    S0 = 7.5e-46 / u.hertz
    
    def noise_spectrum(self, x):
        return x**(-5) + 13*x**-1 + 9*(1+x**2)
    
class VIRGO(Interferometer):
    """
    The VIRGO Interferometer
    """
    name = "VIRGO"
    f0 = 500 * u.hertz
    fs = 20 * u.hertz
    S0 = 3.2e-46 / u.hertz
    
    def noise_spectrum(self, x):
        return (7.8*x)**(-5) + 2*x**(-1) + 0.63 + x**2
    
class EvolvedLISA(Interferometer):
    """
    The eLISA Interferometer
    """
    name = "eLISA"
    frequencies =  np.logspace(-6, 0, 10000) * u.hertz
    L = 1e9*u.meter
    fs = 3e-5 * u.hertz
    def psd(self, frequencies):
        #residual acceleration noise
        sacc = 9e-28 * (1*u.hertz)**4 * (2*np.pi*frequencies)**-4 * (1+(1e-4*u.hertz)/frequencies) * u.meter**2 * u.hertz**-1 # * u.second**-4 
        # shot noise
        ssn = 5.25e-23 * u.meter**2 / u.hertz
        # other measurement noise
        son = 6.28e-23 * u.meter**2 / u.hertz
        #
        s  =(20./3) * (4*(sacc + ssn + son) / self.L**2) * ( 1+ (frequencies/(0.41 * (c.c/(2*self.L))))**2)
        s[frequencies<self.fs]=np.nan
        return s