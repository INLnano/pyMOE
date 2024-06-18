"""
aperture.py


Definition of Class Aperture


"""

import numpy as np
from pyMOE.utils import digitize_array_to_bins
from pyMOE.utils import discretize_array
from pyMOE.sag_functions import phase2height, height2phase



class Aperture:
    """
    Class Aperture:
        Creates an Aperture object that is an homogenous array of values corresponding
        to the transfer function matrix across the aperture
    
    Args:
        :x:         Vector for the x axis
        :y:         Vector for the y axis
    
    Methods:
        :aperture:  returns the aperture
        :shape:     returns the shape of the aperture

    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.XX, self.YY = np.meshgrid(x, y)
        self.pixel_x = self.x[1]-self.x[0]
        self.pixel_y = self.y[1]-self.y[0]
        
        # XX.shape has the same shape as YY.shape
        self.aperture = np.zeros(self.XX.shape)
        self.aperture_original = None
        self.levels = None
        self.aperture_discretized = None
        self.discretized_flag = False

        self.is_height = False

    @property
    def shape(self):
        return self.aperture.shape

    def discretize(self, levels):
        """Discretizes the aperture to the number of levels"""
        if self.aperture_original is None:
            self.aperture_original = np.copy(self.aperture)
        levels, digitized = digitize_array_to_bins(self.aperture, levels)
        
        self.levels = levels
        self.aperture_discretized = digitized
        self.aperture = levels[digitized]
        self.discretized_flag=True

    def modulos(self, mod, normalize_to_max=True):
        """Discretizes the aperture to the number of levels"""
        if self.aperture_original is None:
            self.aperture_original = np.copy(self.aperture)

        aux = self.aperture        
        self.aperture = (aux-np.max(aux)) % (mod)

    

    def phase_unwrap(self):
        """Unwraps the phase of the aperture"""
        assert self.is_height is False, "Cannot unwrap height"
            
        self.aperture = np.unwrap(np.unwrap(self.aperture, axis=0), axis=1)
        # self.aperture = np.apply_over_axes(np.unwrap, self.aperture, np.arange(len(self.aperture.shape)))
        # self.aperture = np.apply_over_axes(np.unwrap, self.aperture, np.arange(len(self.aperture.shape)))


    def phase2height(self, wavelength, n1, n0=1):
        """Converts the phase to height
        Args:
            :wavelength:    Wavelength of the light
            :n1:            Refractive index of the medium where the light is propagating
            :n0:            Refractive index of the medium background"""
        if self.is_height:
            return
        self.aperture = phase2height(self.aperture, wavelength, n1)
        self.is_height = True

    def height2phase(self, wavelength, n1, n0=1):
        """Converts the height to phase
        Args:
            :wavelength:    Wavelength of the light
            :n1:            Refractive index of the medium where the light is propagating
            :n0:            Refractive index of the medium background"""
        if not self.is_height:
            return
        self.aperture = height2phase(self.aperture)
        self.is_height = False



class ApertureField:
    """
    Class Aperture:
        Creates an Aperture object that is an homogenous array of complex values corresponding
        to the transfer function matrix across the aperture
    
    Args:
        :x:         Vector for the x axis
        :y:         Vector for the y axis
    
    Methods:
        :aperture:  returns the complex aperture array
        :amplitude: sets or returns the amplitude of the aperture array
        :phase:     sets or returns the phase of the aperture array
        :unwrap:    retursn the unwrapped phase
        :shape:     returns the shape of the aperture

    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.XX, self.YY = np.meshgrid(x, y)
        self.pixel_x = self.x[1]-self.x[0]
        self.pixel_y = self.y[1]-self.y[0]
        

        self.aperture = np.ones(self.XX.shape)*np.exp(1j*np.zeros(self.XX.shape))

    @property
    def shape(self):
        return self.aperture.shape
    @property
    def amplitude(self):
        return np.abs(self.aperture)

    @amplitude.setter
    def amplitude(self, amplitude):
        assert amplitude.shape == self.shape, "Provided array shape does not match Aperture shape"
        self.aperture = amplitude*np.exp(1j*self.phase)
    
    @property
    def phase(self):
        return np.angle(self.aperture)

    @phase.setter
    def phase(self, phase):
        assert phase.shape == self.shape, "Provided array shape does not match Aperture shape"
        self.aperture = self.amplitude*np.exp(1j*phase)

    @property
    def unwrap(self):
        return np.unwrap(self.phase)


