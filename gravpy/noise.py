import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import numpy.linalg as la
import os


class NoiseSource(object):
    """
    This is the base class for all noise sources which we might want to define or plot.
    """

    pass


