#
import numpy as np
import matplotlib.pyplot as plt
import data.atnf
import scipy.linalg as la
from astropy.io import ascii
from astropy.table import Table, join
import astropy.units as u

from astropy.coordinates import SkyCoord
import functools32
import psrqpy

from interferometers import Detector

from pkg_resources import resource_string, resource_stream, resource_filename

def hellingsdowns_factor(pulsar1, pulsar2):
    """
    The Hellings-Downs factor for two pulsars, normalised
    so that if pulsar1==pulsar2 the returned valueis 1. 
    See equation 20 in http://arxiv.org/pdf/1310.5300v2.pdf
    """
    sep = pulsar1.location.separation(pulsar2.location).to('radian')
    first = (1 - np.cos(sep))/2
    first = first * np.log(first)
    second = (1 - np.cos(sep))/2
    last = 0
    if pulsar1 == pulsar2: return 1 #last = 0.5
    return 1.5*first - 0.25 * second + 0.5 

class Pulsar(object):
    def __init__(self, psrj, cadence, obstime, rms, position):
        """
        Object representing a pulsar.
        
        Parameters
        ----------
        prsj : str
            The Julian (J) name of the pulsar.
        cadence : float
            The observation cadence of the pulsar, in seconds.
        obstime : float
            The observation time of the pulsar.
        rms : float
            The root mean square of the timing residuals for the pulsar.
        position : `astropy.SkyCoord`
            The location of the pulsar in the Sky
        """
        #catalogue = data.atnf.get_atnf()
        #rowdata = catalogue.loc['PSRJ', psrj]
        #self.location = rowdata['POS']
        self.cadence = cadence
        self.rms = rms
        self.obstime = obstime
        self.location = position
        
        
    @property
    def p_vector(self):
        return self.location.cartesian.xyz
    
    def psd(self, frequency):
        """
        Return the power spectral density of the pulsar.
        
        Parameters
        ----------
        frequency : ndarray
            The frequencies at which the PSD should be evaluated.
        
        Returns
        -------
        psd : float
            The pulsar's one-sided PSD.
        """
        #if frequency < 1 / self.obstime: return np.nan
        #if frequency > 1 / self.cadence: return np.nan
        outs = np.ones(len(frequency))
        frequency = frequency * u.hertz
        outs[frequency < 1./self.obstime] = np.nan
        outs[frequency > 1./self.cadence] = np.nan
        return (2 * 1./self.cadence *  self.rms**2)*outs

class TimingArray(Detector):
    pulsars = []
    frequencies = np.logspace(-9, -6)
    def __init__(self, pulsars):
        pulsar_list =  ascii.read(pulsars, delimiter=" ", guess=False)
        pulsar_list.add_index('Name')

        params = ["RAJ", "DECJ"]
        
        print(pulsar_list['Name'])

        query = psrqpy.QueryATNF(params, psrs=pulsar_list['Name'])
        pq = query.get_pulsars()
        
        for pulsar in pulsar_list:
            cadence = pulsar['Cadence']*u.day
            obstime = pulsar['Timespan']*u.year
            pos = SkyCoord(getattr(pq[pulsar], 'RAJ'), getattr(pq[pulsar],'DECJ'), unit=(u.deg, u.deg))
            self.pulsars.append(Pulsar(pulsar['Name'], cadence.to(u.second), obstime.to(u.second), pulsar['RMSRes'], position=pos ))
            
    def plot_array(self):
        """
        Plot the skymap of the array.
        
        Returns
        -------
        matplotlib figure
        """
        locations = np.array([getattr(pulsar, 'location') for pulsar in self.pulsars])
        fig = plt.figure()
        ax = plt.subplot(111, projection="hammer")
        for location in locations:
            plt.plot(location.ra, location.dec, '.', color='b')
        return fig
            
    @functools32.lru_cache(maxsize=1)
    def hdmatrix(self):
        """
        Compute the Hellings-Down matrix for the entire array.
        
        Returns
        -------
        hdmatrix : ndarray
            The matrix of the entire HD matrix
            
        References
        ----------
        [1] 10.1103/PhysRevD.88.124032
        """
        hdmat = np.zeros((len(self.pulsars), len(self.pulsars)))

        for i,pulsar1 in enumerate(self.pulsars):
            for j,pulsar2 in enumerate(self.pulsars):
                hdmat[i,j] = hellingsdowns_factor(pulsar1, pulsar2)
                hdmat[j,i] = hdmat[i,j]
                self.hdm = hdmat
        return hdmat

        
    def effective_pairs(self):
        """
        Compute the "effective number of pairs", that is, the sum over all of the
        Hellings-Downs factors.
        
        Returns
        -------
        effective pairs : float
            The effective number of pairs in the array.
            
        References
        ----------
        [1] 10.1103/PhysRevD.88.124032
        """
        out = 0
        hdmat = self.hdmatrix()
        for i in xrange(len(hdmat[0])):
            for j in xrange(i+1, len(hdmat[0])): 
                out += hdmat[i,j]**2
        return out
    
    def psd(self, frequency=None):
        """
        Compute the effective power spectral density of the array.
        
        Parameters
        ----------
        frequency : float
            The frequency at which to evaluate the PSD. Defaults to the
            frequencies set for the object.
        
        Returns
        -------
        effective pairs : float
            The effective psd of the array
            
        References
        ----------
        [1] 10.1103/PhysRevD.88.124032
        
        """
        if frequency == None: frequency = self.frequencies
        out = 0
        hdmat = self.hdmatrix()
        for i in xrange(len(hdmat[0])):
            for j in xrange(i+1, len(hdmat[0])): 
                out += hdmat[i,j]**2 / (self.pulsars[i].psd(frequency) * self.pulsars[j].psd(frequency))
        eff_psd = (12*np.pi*frequency**2) *out**(-0.5)
        return eff_psd

##

class IPTA(TimingArray):
    """
    The international pulsar timing array.
    """
    name = "IPTA"
    def __init__(self):
        data_file = resource_filename(__name__, "data/IPTA-pulsars.dat")
        pulsar_list =  ascii.read(data_file, delimiter=" ", guess=False)
        pulsar_list.add_index('Name')

        params = ["RAJ", "DECJ"]

        query = psrqpy.QueryATNF(params, psrs=list(pulsar_list['Name']))
        pq = query.get_pulsars()
        
        for pulsar in pulsar_list:
            cadence = pulsar['Cadence']*u.day
            obstime = pulsar['Timespan']*u.year
            pos = SkyCoord(getattr(pq[pulsar['Name']], 'RAJ'), getattr(pq[pulsar['Name']],'DECJ'), unit=(u.deg, u.deg))
            self.pulsars.append(Pulsar(pulsar['Name'], cadence.to(u.second), obstime.to(u.second), pulsar['RMSRes'], position=pos ))
            
        
        # for pulsar in pulsar_list:
        #     cadence = pulsar['Cadence']*u.day
        #     obstime = pulsar['Timespan']*u.year
        #     self.pulsars.append(Pulsar(pulsar['Name'], cadence.to(u.second), obstime.to(u.second), pulsar['RMSRes'] ))
    
