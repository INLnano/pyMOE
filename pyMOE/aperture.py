"""
mask.py


Definition of Class Apertures


"""
import numpy as np

class Aperture:
    """
    Class Aperture:
        Creates an Aperture object that is an homogenous array of complex values corresponding
        to the transfer function matrix across the aperture
    
    Args:
        x = Vector for the x axis
        y = Vector for the y axis
    
    Methods:
        aperture: returns the complex aperture array
        amplitude: sets or returns the amplitude of the aperture array
        phase: sets or returns the phase of the aperture array
        unwrap: retursn the unwrapped phase
        shape: returns the shape of the aperture

    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.XX, self.YY = np.meshgrid(x, y)

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


