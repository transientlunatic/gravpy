import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from . import general
from .data import atnf as atnf


class Source():
    """
    The base class for a gravitational wave source.
    """
    name = "Generic Source"
    frequencies =  np.logspace(-5, 5, 1000) * u.hertz
    M = 30 * u.solMass
    r = 300 * u.parsec
    
    def __init__(self, frequencies=None, M=None, r=None):
        if frequencies: self.frequencies = frequencies
        if r: self.r = r
        if M: self.M = M
      
    def raw_strain(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return ((1./self.r) * ((5*np.pi)/(24*c.c**3))**(0.5) * (c.G * self.chirp_mass())**(5./6) * (np.pi*frequencies)**(-7./6)).to(1/u.hertz)
    
    def psd(self, frequencies=None):
        """
        The one-sided power spectral density
        
        Parameters
        ----------
        frequencies : ndarray
            An array of frequencies where the PSD should be calculated.
            
        Returns : ndarray
            An array of the PSDs at the given frequencies for this source.
        """
        if not frequencies: frequencies = self.frequencies
        return 2 * (frequencies**0.5) * np.abs(self.raw_strain(frequencies))
    
    def srpsd(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return np.sqrt(self.psd(frequencies)) 
        
    def characteristic_strain(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return np.sqrt(4 * frequencies**2 * np.abs(self.raw_strain(frequencies))**2)
    
    def energy_density(frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return (2*pi**2)/3 * frequencies**3 * self.psd(frequencies)
    
    def plot(self, axis, label=None):
        if axis:
            if not label:
                label = self.name
            axis.loglog(self.frequencies, self.characteristic_strain(self.frequencies), label=label, lw=2)
            axis.set_xlabel('Frequency [Hz]')
            #axis.set_ylabel('Root Noise Power spectral density')
            axis.legend()
            
    def snr(self, detector):
        return general.snr(self, detector)

class Pulsar(Source):
    """
    A gravitational-wave pulsar.
    """
    name = "Pulsar"

    def __init__(self, psrj,  Izz=1e-5 * 10**38 * u.kilogram * u.meter**2):
        """
        Object representing a pulsar.
        
        Parameters
        ----------
        prsj : str
            The Julian (J) name of the pulsar.
        Izz : float
           The magnitude of the zz component of the moment of inertia tensor.

        """
        self.Izz = Izz
        catalogue = atnf.get_atnf()
        rowdata = catalogue.loc['PSRJ', psrj]
        self.data = rowdata
        self.name = psrj

    def raw_strain(self, frequencies = None):
        """Calculate the raw strain which the pulsar should produce. Note
        that unlike other sources this will be at a single frequency,
        since pulsars are not broadband emitters.

        Parameters
        ----------
        
        """
        if not frequencies: frequencies = self.frequencies
        response = np.ones(len(frequencies)) * np.nan
        def find_nearest(array,value):
            idx = (np.abs(array-value)).argmin()
            return idx
        response[find_nearest(frequencies, 2*self.data['F0']*u.hertz)] = 1
        distance = self.data['DIST'] * 1000 * u.parsec
        f = 2*self.data['F0'] * u.hertz
        fdot = self.data['F1']
        fratio = fdot / f
        GoC = c.G / c.c**3
        rational = - (5.0/4.0) * GoC * self.Izz * fratio
        return response * (1/distance) * np.sqrt(rational)
    
    def plot(self, axis):
        if axis: 
            axis.loglog(self.frequencies, self.characteristic_strain(self.frequencies), 'o', label=self.name,)
            axis.set_xlabel('Frequency [Hz]')
            #axis.set_ylabel('Root Noise Power spectral density')
            axis.legend()
    
class Type1ASupernova(Source):
    """
    A Type-1A supernova source. Based on https://arxiv.org/abs/1511.02542.
    """
    name = "Type Ia SN"
    r = 10 * 1000 * u.parsec
    
    def __init__(self, frequencies = None, r = None):
        if frequencies: self.frequencies = frequencies
        if r: self.r = r

    def raw_strain(self, frequencies = None):
        if not frequencies: frequencies = self.frequencies
        response = np.ones(len(frequencies))
        response *= ((9e-21) *   (10 * 1000 * u.parsec) / self.r)
        response[frequencies < 0.25 * u.hertz ] = np.nan
        response[frequencies > 1.5 * u.hertz ] = np.nan
        
        return response

class CoreCollapseSupernova(Source):
    """
    A core-collapse supernova source. Based on Dimmelmeier.
    """
    name = "CCSN"
    r = 10 * 1000 * u.parsec
    
    def __init__(self, frequencies = None, r = None):
        if frequencies: self.frequencies = frequencies
        if r: self.r = r

    def raw_strain(self, frequencies = None):
        if not frequencies: frequencies = self.frequencies
        return np.ones(len(frequencies)) * ((8.9e-21) * (10 * 1000 * u.parsec) / self.r)
    
    
class CBC(Source):
    """
    A compact binary coallescence source
    """
    name = "CBC"
    M = 30 * u.solMass
    r = 300 * u.parsec
    
    def __init__(self, frequencies=None, m1=None, m2=None, r=None):
        if frequencies: self.frequencies = frequencies
        if r: self.r = r
        if m1: self.m1 = m1
        if m2: self.m2 = m2
        self.M = self.chirp_mass()
        
    def fdot(self, frequencies=None, M=None):
        """
        Calculate the first time derivative of the CBC's frequency.
        
        Parameters
        ---------
        frequencies : ndarray
            The frequencies at which the number of cycles need to be found.
            
        M : float
            The chirp mass of the CBC.
            
        Returns
        -------
        fdot : ndarray
            The df/dt of each frequency.
        """
        if not frequencies: frequencies = 0.5*self.frequencies
        if not M: M = self.chirp_mass()
        return (((96*np.pi**(8./3)) / (5 * c.c**5)) * (c.G*M)**(5./3) * frequencies**(11./3))#.to(u.hertz**2)

    def ncycles(self, frequencies=None, M=None):
        """
        Calculate the number of cycles that the CBC spends in each frequency bin.
        
        Parameters
        ---------
        frequencies : ndarray
            The frequencies at which the number of cycles need to be found.
            
        M : float
            The chirp mass of the CBC.
            
        Returns
        -------
        ncycles : ndarray
            The number of cycles in each frequency bin.
        """
        if not frequencies: frequencies = 0.5*self.frequencies
        if not M: M = self.chirp_mass()
        return np.sqrt(frequencies**2/ self.fdot(frequencies, M))#.to(1)
    
    def characteristic_strain(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return np.sqrt(2*self.ncycles())*np.sqrt(4 * frequencies**2 * np.abs(self.raw_strain())**2)
    
    def chirp_mass(self):
        return ((self.m1*self.m2)**(3./5) / (self.m1 + self.m2)**(1./5)).to(u.kilogram)
    
    def fisco(self):
        return ((c.c**3) / (np.pi*c.G*(self.m1+self.m2)*6*6**0.5 )).to(u.hertz)
    
    def raw_strain(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        h = ((1./self.r) * ((5*np.pi)/(24*c.c**3))**(0.5) * (c.G * self.M)**(5./6) * (np.pi*frequencies)**(-7./6)).to(1/u.hertz)
        h[frequencies>2*self.fisco()] = np.nan
        return h



class IMR(Source):
    """
    An inspiral, merger, ringdown frequency spectrum.

    Modelled on IMRPhenomA, and does not include contributions from spin.
    """

    def __init__(self, frequencies=None, m1=None, m2=None, r=None):
        if frequencies: self.frequencies = frequencies
        self.distance = r.to(u.meter)
        self.mass1 = m1.to(u.kilogram)
        self.mass2 = m2.to(u.kilogram)
        
    @property
    def eta(self):
        """
        The symmetric mass ratio of the CBC system.
        """
        eta = (self.mass1 * self.mass2) / (self.mass1 + self.mass2)**2
        return eta
    
    def fk(self, k):

        # The various transition frequencies.
        # Broadly
        # 0 is the merger,
        # 1 is the ringdown
        # 2 decay width
        # 3 cut-off frequency
        a = [2.9740e-1, 5.9411e-1, 5.0801e-1, 8.4845e-1]
        b = [4.4810e-2, 8.9794e-2, 7.7515e-2, 1.2848e-1]
        d = [9.5560e-2, 1.9111e-1, 2.2369e-2, 2.7299e-1]
        
        top = a[k] * self.eta**2 + b[k] * self.eta + d[k]
        bot = np.pi * (c.G*(self.mass1+self.mass2) / c.c**3)
        return top / bot
    
    @property
    def chirp_mass(self):
        return ((self.mass1*self.mass2)**(3./5) / (self.mass1 + self.mass2)**(1./5)).to(u.kilogram)
    
    
    @property
    def w(self):
        first = (np.pi * self.fk(2)/2)
        second = (self.fk(0) / self.fk(1))**(2./3)

        return first * second

    def L(self, f):
        first = (1/(2*np.pi))
        second = (self.fk(2)/((f - self.fk(1))**2 + self.fk(2)**2/4.))

        return first * second

    def amplitude(self, f):
        first = np.sqrt(5./24)
        second = (c.G * self.chirp_mass / c.c**3)**(5./6) * (self.fk(0))**(-7./6)
        third = (np.pi**(2/3.) * (self.distance / c.c))

        tail = np.ones(len(f))*np.nan
        tail[f<self.fk(0)] = (f[f<self.fk(0)]/self.fk(0))**(-7./6)
        tail[(self.fk(0)<f) & (f<self.fk(1))] = (f[(self.fk(0)<f) & (f<self.fk(1))] / self.fk(0))**(-2/3.)
        tail[(self.fk(1)<f) & (f<self.fk(3))] = self.w * self.L(f[(self.fk(1)<f) & (f<self.fk(3))])

        return first * (second/third) * tail

    def raw_strain(self, frequencies):
        return self.amplitude(frequencies)
