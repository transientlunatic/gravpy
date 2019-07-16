import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import numpy.linalg as la
import os

def rot_z(phi):
    return np.array([
            [np.cos(phi), -np.sin(phi), 0],
            [np.sin(phi), np.cos(phi), 0],
            [0,0,1]
        ])

def rot_y(psi):
    return np.array([
            [np.cos(phi), 0,  -np.sin(phi)],
            [0, 1,   0],
            [np.sin(phi),           0,             np.cos(phi)]
        ])

def rot_x(theta):
    return np.array([
            [1,             0,              0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta),  np.cos(theta)]
        ])

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
        """
        Produce the sensitivity curve of the detector in terms of the 
        energy density.
        
        Parameters
        ----------
        frequencies : ndarray
            An array of frequencies, in units of Hz
            
        Returns
        -------
        energy_density : ndarray
            An array of the dimensionless energy density of the sensitivity of
            the detector.
        """
        if not frequencies: frequencies = self.frequencies
        bigH = (2*np.pi**2)/3 * frequencies**3 * self.psd(frequencies)
        littleh = bigH / ((100*u.kilometer / u.second / u.megaparsec).to(u.hertz))**2
        return littleh
    
    def srpsd(self, frequencies=None):
        if not frequencies: frequencies = self.frequencies
        return np.sqrt(self.psd(frequencies))
    
    def plot(self, axis=None, **kwargs):
        """
        Plot the noise curve for this detector.
        """
        if "lw" not in kwargs.keys():
            lw=2
        else:
            
            lw=kwargs.pop('lw')
        
        if axis: 
            axis.loglog(self.frequencies, self.noise_amplitude(), label=self.name, lw=lw, **kwargs)
            axis.set_xlabel('Frequency [Hz]')
            axis.set_ylabel('Characteristic Strain')
        

class Interferometer(Detector):
    """
    The base class to describe an interferometer.
    """
    name = "Generic Interferometer"
    f0 = 150 * u.hertz
    fs = 40 * u.hertz
    S0 = 1e-46 / u.hertz
    frequencies =  np.logspace(1, 5, 4000) * u.hertz
    
    xhat = np.array([1,0,0])
    yhat = np.array([0,1,0])
    zhat = np.array([0,0,1])
    length = 4 * u.kilometer
    
    detector_tensor = length * (np.outer(xhat, xhat) - np.outer(yhat, yhat))
    
    def __init__(self, frequencies=None, configuration=None, obs_time=None):
        if frequencies: self.frequencies = frequencies
        self.configuration = configuration
        self.obs_time = obs_time
        
        if configuration: 
            self.name = "{} [{}]".format(self.name, configuration)
    
    def noise_spectrum(self, x):
        """
        Return the default noise spectrum for the interferometer.
        
        Parameters
        ----------
        x : float
            The rescaled frequency at which the fit should be evaluated.
            
        Returns
        -------
        float
            The sensitivity of the interferometer at the frequency provided.
        """
        return (3.4*x)**(-30) + 34*x**(-1) + (20 * (1 - x**2 + 0.4*x**4))/(1 + 0.5*x**2)
    
    def psd(self, frequencies=None):
        """
        Calculate the one-sided power spectral desnity for a detector. 
        If a particular configuration is specified then the results will be
        returned for a spline fit to that configuration's curve, if available.
        
        Parameters
        ----------
        frequencies : ndarray
            An array of frequencies where the PSD should be evaluated.
            
        configuration : str
            The configuration of the detector for which the curve should be returned.
        """
        if not frequencies: frequencies = self.frequencies
            
        if self.configuration:
            configuration = self.configuration
            d_frequencies, d_sensitivity = self.configurations[configuration]
            d_frequencies, d_sensitivity = np.genfromtxt(os.path.join(os.path.dirname(__file__),d_frequencies)), np.genfromtxt(os.path.join(os.path.dirname(__file__),d_sensitivity))
            tck = interpolate.splrep(d_frequencies, d_sensitivity, s=0)
            interp_sensitivity = interpolate.splev(frequencies, tck, der=0)
            interp_sensitivity[frequencies<self.fs]=np.nan
            return (interp_sensitivity)**2 * u.hertz**-1
            
        
        x = frequencies / self.f0
        xs = self.fs / self.f0
        sh = self.noise_spectrum(x)

        if self.obs_time:
            sh /= (self.obs_time.to(u.second))
        
        sh[frequencies<self.fs]=np.nan
        return sh * self.S0

    def antenna_pattern(self, theta, phi, psi):
        """
        Produce the antenna pattern for a detector, given its detector tensor, 
        and a set of angles.
        
        Parameters
        ----------
        theta : float
            The altitude angle.
        phi : float
            The azimuthal angle.
        psi : float or list
            The polarisation angle. If psi is a list of two angles the returned 
            antenna patterns will be the integrated response between those two 
            polsarisation angles.
            
        Returns
        -------
        F+ : float
            The antenna response to the '+' polarisation state.
        Fx : float
            The antenna response to the 'x' polsarisation state.
        |F| : float
            The combined antenna response (sqrt(F+^2 + Fx^2)).
        """
        detector = self.detector_tensor / self.length
        # The unrotated basis of the gravitational wave
        e = np.array([
            [1,0,0],
            [0,1,0],
            [0,0,1]
        ])
        # Calculate the rotated basis
        # Rotate phi about z
        # Rotate theta about x
        # Rotate psi about z
        #rot_basis = np.dot(np.dot(np.dot(np.dot(dhat,rot_x(theta)), rot_z(phi)), rot_z(psi)), e)
        rot_basis = np.dot( np.dot( rot_x(theta), rot_z(phi)), e)


        def plus_polarisation(psi, rot_basis):
            alpha, beta, _ = rot_basis
            rot_basis = np.dot(rot_basis, rot_z(psi))
            return np.outer(alpha, alpha) - np.outer(beta, beta)
        def cross_polarisation(psi, rot_basis):
            alpha, beta, _ = rot_basis
            ot_basis = np.dot(rot_basis, rot_z(psi))
            return np.outer(alpha, beta) + np.outer(beta, alpha)

        # Now the antenna pattern
        if isinstance(psi, list):
            fplus  = integrate.quad(lambda psi: (detector*plus_polarisation(psi, rot_basis)).sum(),  psi[0], psi[1])[0]
            fcross = integrate.quad(lambda psi:((detector* cross_polarisation(psi, rot_basis)).sum()),  psi[0], psi[1])[0]
        else:
            fplus = (detector*plus_polarisation(psi, rot_basis)).sum()
            fcross = (detector* cross_polarisation(psi, rot_basis)).sum()

        return np.abs(fplus), np.abs(fcross), np.sqrt(fplus**2 + fcross**2)

    def skymap(self, nx=200, ny=100, psi=[0, np.pi]):
        """
        Produce a skymap of the antenna repsonse of the interferometer.
        
        Parameters
        ----------
        nx : int
            The number of locations along the horizontal axis to produce the map at
            defaults to 200.
        ny : int
            The number of locations along the vertical axis to produce the map at
            defaults to 100
        psi : float or list
            The polarisation angle to produce the map at. If a list is given then the integrated 
            response is given between those angles.
            
        Returns
        -------
        x : ndarray
            The x values for the map
        y: ndarray
            The y values for the map
        antennap : ndarray
            The values of the sensitivity in the + polarisation
        antennax : ndarray
            The values of the sensitivity in the x polarisation
        antennac : ndarray
            The values of the combined polarisation sensitivities
        """
        
        # Note these are, confusingly, the wrong way 
        # around, and I should fix them.
        x = np.linspace(0, np.pi, ny)
        y = np.linspace(0, 2*np.pi, nx)
        xv, yv = np.meshgrid(x,y)

        H = np.zeros((nx, ny))
        A = np.zeros((nx, ny))
        B = np.zeros((nx, ny))
        
        for i in range(nx):
            for j in range(ny): 
                A[i,j], B[i,j], H[i,j] = self.antenna_pattern(xv[i,j], yv[i,j],psi)
        
        return x, y, A, B, H

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


class BDecigo(Interferometer):
    # B-DECIGO noise curve in arxiv:1802.06977
    S0 = 4.040e-46 * u.hertz**-1

    frequency_range = [1e-2, 1e2] * u.hertz
    frequencies =  np.logspace(-2, 2, 10000) * u.hertz

    def psd(self, frequencies):
        return self.S0 * (1.0 + 1.584e-2 * frequencies.value**-4 + 1.584e-3 * frequencies.value**2)

class Decigo(Interferometer):
    """
    The full, original Decigo noise curve, from  arxiv:1101.3940.
    """

    fp = 7.36 * u.hertz

    frequency_range = [1e-2, 1e2] * u.hertz
    frequencies =  np.logspace(-2, 2, 10000) * u.hertz
    
    def psd(self, frequencies):
        """
        The power spectrum density of the detector, taken from equation 5 of arxiv:1101.3940.
        """
        first_c = 7.05e-48
        second_c = 4.8e-51
        third_c = 5.33e-52

        first = first_c * (1 + (frequencies / self.fp)**2)
        second =  second_c * frequencies.value**-4 / (1+(frequencies/self.fp)**2)
        third = third_c * frequencies.value**-4
        
        return  (first + second + third) * u.hertz**-1

class BigBangObservatory(Interferometer):
    """
    The Big Bang Observatory.
    """

    frequency_range = [1e-3, 1e2] * u.hertz
    frequencies =  np.logspace(-3, 2, 10000) * u.hertz

    def psd(self, frequencies):
        """
        The power spectrum density of the detector, taken from equation 6 of arxiv:1101.3940.
        """
        first = 2.00e-49 * frequencies.value**2
        second = 4.58e-49
        third = 1.26e-51*frequencies.value**-4

        return (first + second + third)*u.hertz**-1
    
class AdvancedLIGO(Interferometer):
    """
    The aLIGO Interferometer
    """
    name = "aLIGO"
    f0 = 215 * u.hertz
    fs = 20 * u.hertz
    S0 = 1.0e-49 / u.hertz
    
    frequency_range = [30, 4e3] * u.hertz
    xhat = np.array([1,0,0])
    yhat = np.array([0,1,0])
    zhat = np.array([0,0,1])
    length = 4 * u.kilometer
    
    detector_tensor = length * (np.outer(xhat, xhat) - np.outer(yhat, yhat))
    
    configurations = {
            'O1': ['data/aligo_freqVector.txt', 'data/o1_data50Mpc_step1.txt']
                      }
    
    def noise_spectrum(self, x):
        return (x)**(-4.14) -5*x**(-2) + ((111 * (1-x**2 +0.5*x**4))/(1+0.5*x**2))

class EinsteinTelescope(Interferometer):
    """
    The Einstein Telescope.
    """
    name = "Einstein Telescope"
    f0 = 1.0 * u.hertz

    frequency_range = [f0, 1e4*u.hertz]

    frequencies =  np.logspace(0, 4, 4000) * u.hertz
    
    configurations = {
        "ET-D-Sum": "data/et-d-curve.txt",
        }

    def __init__(self, frequencies=None, configuration="ET-D-Sum", obs_time=None):
        """
        Create a new Einstein Telescope object.
        By default the ET-D configuration is used, and the PSD is the sum of the two interferometers' sensitivity curves.
        """
        
        if frequencies: self.frequencies = frequencies
        self.configuration = configuration
        self.obs_time = obs_time
        
        if configuration: 
            self.name = "{} [{}]".format(self.name, configuration)


    def psd(self, frequencies=None):
        """
        Calculate the one-sided power spectral desnity for a detector. 
        If a particular configuration is specified then the results will be
        returned for a spline fit to that configuration's curve, if available.
        
        Parameters
        ----------
        frequencies : ndarray
            An array of frequencies where the PSD should be evaluated.
            
        configuration : str
            The configuration of the detector for which the curve should be returned.
        """
        if not frequencies: frequencies = self.frequencies


        # The ET curves are all given as PSDs
        if self.configuration:
            configuration = self.configuration
            datafile = self.configurations[configuration]
            data = np.genfromtxt(os.path.join(os.path.dirname(__file__), datafile))

            d_frequencies = data[:,0]

            # This would almost definitely be better handled by splitting these curves into their own files.
            if self.configuration == "ET-D-Sum":
                col = 3
            
            d_sensitivity = data[:,col]
            
            tck = interpolate.splrep(d_frequencies, d_sensitivity, s=0)
            interp_sensitivity = interpolate.splev(frequencies, tck, der=0)
            interp_sensitivity[frequencies<self.fs]=np.nan
            return (interp_sensitivity)**2 * u.hertz**-1
            
        
        x = frequencies / self.f0
        xs = self.fs / self.f0
        sh = self.noise_spectrum(x)

        if self.obs_time:
            sh /= (self.obs_time.to(u.second))
        
        sh[frequencies<self.fs]=np.nan
        return sh * self.S0
# Make a little shim so you can call EinsteinTelscope as ET
ET = EinsteinTelescope    
    
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
    
class Virgo(Interferometer):
    """
    The Virgo Interferometer
    """
    name = "Virgo"
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

class LISA(Interferometer):
    """
    The LISA Interferometer in its mission-accepted state, as of 2018
    """
    name = "eLISA"
    frequencies =  np.logspace(-5, 0, 10000) * u.hertz
    L = 2.5e9*u.meter
    fstar = 19.08*1e-3 * u.hertz
    fs = 3e-5 * u.hertz

    def metrology_noise(self, frequencies):
        """
        Calculate the noise due to the single-link optical metrology, from arxiv:1803.01944.
        """
        first =  (1.5e-11 * u.meter)**2
        second = (1 * u.dimensionless_unscaled +(2e-3*u.hertz/frequencies)**4) * u.hertz**-1

        return first*second

    def single_mass_noise(self, frequencies):
        """
        The acceleration noise for a single test mass.
        """
        first = (3e-15 * u.meter * u.second**-2)**2
        second = (1 * u.dimensionless_unscaled +(0.4e-3*u.hertz/frequencies)**2)
        third = (1 * u.dimensionless_unscaled+(frequencies/(8e-3*u.hertz))**4)*u.hertz**-1

        return first*second*third

    def confusion_noise(self, frequencies, observation_time=0.5):
        """
        The noise created by unresolvable galactic binaries at low frequencies.
        """
        amp = 9e-45

        alpha = {0.5: 0.133, 1: 0.171, 2: 0.165, 4: 0.138}
        beta  = {0.5: 243, 1: 292, 2: 299, 4: -221}
        kappa = {0.5: 482, 1: 1020, 2: 611, 4: 521}
        gamma = {0.5: 917, 1: 1680, 2: 1340, 4: 1680}
        fk = {0.5: 0.00258, 1: 0.00215, 2: 0.00173, 4: 0.00113}

        first = amp * frequencies**(-7./3.) * np.exp((- frequencies**alpha[observation_time]).value
                                                      + (beta[observation_time] * frequencies * np.sin((kappa[observation_time] * frequencies).value)).value)
        second = (1+np.tanh((gamma[observation_time] * (fk[observation_time]*u.hertz - frequencies)).value))

        return (first * second).value * (u.hertz**-1)
        
    
    def psd(self, frequencies):
        # See https://arxiv.org/pdf/1803.01944.pdf for this

        first = (10 / 3 * self.L**-2)
        second = (self.metrology_noise(frequencies) + (4*self.single_mass_noise(frequencies)/(2*np.pi*frequencies)**4))
        third = (1 * u.dimensionless_unscaled + (6./10)*(frequencies / self.fstar)**2)
        return (first*second*third) + self.confusion_noise(frequencies)
    
