"""
mask.py


Definition of Class Aperture


"""

import numpy as np
from pyMOE.utils import digitize_array_to_bins
from pyMOE.utils import discretize_array


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

    @property
    def shape(self):
        return self.aperture.shape

    def discretize(self, levels):
        """Discretizes the aperture to the number of levels"""
        if self.aperture_original is None:
            self.aperture_original = np.copy(self.aperture)
        levels, digitized = digitize_array_to_bins(self.aperture_original, levels)
        
        self.levels = levels
        self.aperture_discretized = digitized
        self.aperture = levels[digitized]
        self.discretized_flag=True



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


