"""
ensemble.py


Definition of Class Ensemble


"""

import numpy as np

import pyMOE.aperture as aperture 
import pyMOE.field as field 

from pyMOE.aperture import Aperture
from pyMOE.field import Field
from pyMOE.field import Screen

import scipy.fftpack as sfft 

class Ensemble:
    """
    Class Ensemble:
        Creates a object with an optical ensemble of apertures, fields and screens 
        
    Args:
        :aperture_array:               array of apertures (Aperture objects)
        :screen_array:                 array of screens (Screen objects)
        :propagation_methods_array:    propagators from propagate (optimize) module to use (function names)
        :wavelength_array:             array of wavelengths (scalar)
        :input_light_field:            input light (1 Field object)
        :topografies_array:            (optional) if given consider it as the topografies from where the phases can be obtained 
    
    Methods:
        :fields:         returns the fields at the apertures 
        :screens:        returns the fields stored 
        :overall_screen: returns the field stored in an overall screen 
    """
    def __init__(self, aperture_array_amp, aperture_array_phase, screen_array,  wavelength_array, input_light_field, topographies_array=None, propagation_methods=None):
        #if topographies_array==None:
            #print("Assuming the Apertures to contain already the actual phases of the objects.")
            self.aperture_array_amp = aperture_array_amp
            self.aperture_array_phase = aperture_array_phase
            self.screen_array = screen_array 
            self.propagation_methods = propagation_methods 
            self.wavelength_array = wavelength_array
            self.input_light_field = input_light_field
            
            
            self.field_array = [field.create_empty_field_from_aperture(ap) for ap in self.aperture_array_amp ] 

            self.overall_screen = np.dstack(np.array(screen_array))

        ##TODO: CHECK ALL arrays have same dimension (apertures, screens, propagation-METHODS)
        #TODO     
        #else: 
        #    #calculates the phase from index of refraction and assigns it to the Aperture objects 
        

    @property
    def fields(self):
        return self.field_array
    @property
    def screens(self):
        return self.screen_array
    @property
    def apertures_amplitudes(self):
        return self.aperture_array_amp
    @property
    def apertures_phases(self):
        return self.aperture_array_phase
    @property
    def apertures_xar(self):
        flat_apertures = np.array([phase.aperture.flatten() for phase in self.aperture_array_phase])
        return flat_apertures           
