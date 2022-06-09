"""
Library module for hologram generation and visualiation
"""


from PIL import Image
import numpy as np

from pyMOE.utils import progress_bar, Timer, mean_squared_error, discretize_array

class Field():
    """ 
    Class to define complex fields with methods for amplitude and phase"""
    def __init__(self, a,b):
        self.signal = a+1j*b
    
    @property
    def phase(self):
        return np.angle(self.signal)
    
    @phase.setter
    def phase(self, phase):
        assert phase.shape == self.signal.shape, "Provided array shape must match existing signal"
        self.signal = self.amplitude*np.exp(1j*phase)

    @property
    def amplitude(self):
        return np.abs(self.signal)

    @amplitude.setter
    def amplitude(self, amplitude):
        assert amplitude.shape == self.signal.shape, "Provided array shape must match existing signal"
        self.signal = amplitude*np.exp(1j*self.phase)
        
    @property
    def intensity(self):
        return np.power(self.amplitude,2)


def algorithm_Gerchberg_Saxton(target_intensity, iterations=3, levels=None, input_phase=None, source_beam=None, verbose=True):
    """
    Creates hologram phase mask to generate target intensity using Gerchberg Saxton Algorithm.
    
    U0 - input field
    U1 - calculated far-field
    
    Args:
        target_intensity: 2D array of intensity values 
        levels: Scalar or array: levels to consider in the phase mask as physical constraint. If None, it does not discretize.
        input_phase: If given, uses this input_phase as starting phase instead of random.
    
    Returns:
        phase_mask: 2D phase mask
        error_list: list of errors measured in each iteration
    """
    shape = target_intensity.shape
        
    if levels is not None:
        if isinstance(levels, int):
            levels = np.linspace(-np.pi, np.pi, levels, endpoint=False)
            
    # Due to the way numpy fft works, we must first fftshift all fields
    target_intensity = np.fft.fftshift(target_intensity)
    
    if input_phase is not None:
        input_phase = np.fft.fftshift(input_phase)
        field_0.phase = input_phase
    else: 
        input_phase = np.random.random(shape)*2*np.pi
    if source_beam is not None:
        source_beam = np.fft.fftshift(source_beam)
    else:
        source_beam = np.ones(shape)
    
    field_0 = Field(source_beam, input_phase)
    
    field_1 = Field(np.ones(shape), input_phase)

    list_iteration_errors = []
    
    with Timer("Gerchberg Saxton Algorithm"):
        for i in range(iterations):

            # Add input beam amplitude
            field_0.amplitude = source_beam

            # Apply physical constraints

            # Discretize phase array to the levels
            if levels is not None:
                physical_phase = discretize_array(field_0.phase, levels) 
                field_0.phase = physical_phase

            # Output_phase is here
            phase_mask = field_0.phase

            # Calculate forward Fourier Transform
            field_1.signal = np.fft.fft2(field_0.signal)

            # Calculate error metrics
            # TODO
            norm_amplitude = field_1.amplitude/field_1.amplitude.max()
            error = mean_squared_error(norm_amplitude, target_intensity)

            list_iteration_errors.append(error)
            field_1.amplitude = target_intensity

            # Calculate inverse Fourier Transform
            field_0.signal = np.fft.ifft2(field_1.signal)

            if verbose:
                progress_bar(i/iterations)

        if verbose:
            progress_bar(1)
    phase_mask = np.fft.fftshift(phase_mask)
    return phase_mask, list_iteration_errors

def calculate_phase_farfield(phase, source_beam=None):
    """
    Calculates a proportional far field of a given phase aperture and source_beam
    
    TODO Fraunhofer propagation
    
    Args
        phase: phase mask
        source_beam: if not given, assumes constant intensity=1
        
    Returns
        far_field: far field amplitude
    """
    
    shape = phase.shape
    if source_beam is not None:
        source_beam = np.fft.fftshift(source_beam)
    else:
        source_beam = np.ones(shape)
    phase_mask = np.fft.fftshift(phase)
    field_0 = Field(np.ones(shape), np.zeros(shape))
    field_0.amplitude = source_beam
    field_0.phase = phase_mask
    field_1 = Field(np.ones(shape), np.zeros(shape))
    
    field_1.signal = np.fft.fft2(field_0.signal)
    field_1.signal = np.fft.fftshift(field_1.signal)
    
    return field_1.amplitude
    