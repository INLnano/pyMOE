"""
field.py


Definition of Class Field and related functions


"""

import numpy as np


from pyMOE.aperture import Aperture



class Field:
    """
    Class Field:
        Creates an E Field object is an of values corresponding
        to an input or modulated electric field 
        
        
    Args:
        :x:         Vector for the x axis
        :y:         Vector for the y axis
    
    Methods:
        :field:     returns the field
        :shape:     returns the shape of the field

    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.XX, self.YY = np.meshgrid(x, y) # indexing='ij')
        self.pixel_x = self.x[1]-self.x[0]
        self.pixel_y = self.y[1]-self.y[0]
    
        self.field = np.zeros(self.XX.shape, dtype=complex)
        
    @property
    def shape(self):
        return self.field.shape
    @property
    def amplitude(self):
        return np.abs(self.field)
    @property
    def phase(self):
        return np.angle(self.field)
    @property 
    def intensity(self):
        return np.abs(self.field)**2
    

    
def create_empty_field(xmin, xmax, N_x, ymin, ymax, N_y):
    """
    Creates an empty field max of the mesh dimensions provided
    
    Args: 
        :xmin, xmax:    range for x 
        :N_x:           number of x points
        :ymin, ymax:    range for y 
        :N_y:           number of y points
    
    Returns:
        :field: empty Field
    """
    x = np.linspace(xmin, xmax, N_x)
    y = np.linspace(ymin, ymax, N_y)
    
    return Field(x,y)

def create_empty_field_from_field(field):
    """
    Creates an empty field with the same spatial dimensions of the given field
    
    Args:
        :field: field
    Returns:
        :field: empty Field of same spatial dimensions
    """
    assert type(field) is Field, "aperture must be of type Field"


    return Field(field.x, field.y)

def create_empty_field_from_aperture(aperture):
    """
    Creates an empty field with the same spatial dimensions of the given aperture
    but does not modulate the field.
    
    Args:
        :aperture: aperture
    Returns:
        :field: empty Field of same spatial dimensions
    """
    assert type(aperture) is Aperture, "aperture must be of type Aperture"


    return Field(aperture.x, aperture.y)


def modulate_field(field, amplitude_mask=None, phase_mask=None):
    """
    Modulates the input field with the given amplitude or phase mask.
    If amplitude_mask is not given, it assumes an amplitude of 1 for the amplitude modulatiom.
    If phase is not given, it assumes a phase of 0 for the phase modulation
    The field dimensions and the provided masks must be the same otherwise it raises an error.
    
    Args:
        :field: Input field to be modulated
        :amplitude_mask: Mask of amplitude aperture
        :phase_mask: Mask of phase aperture
    Returns:
        :field: returns the modulated field.
    """
    assert type(field) is Field, "field must be of type Field"


    modulation_amplitude = np.ones(field.XX.shape)
    modulation_phase = np.zeros(field.XX.shape)

    if amplitude_mask is not None:
        assert type(amplitude_mask) is Aperture, "amplitude_mask must be of type Aperture"
        assert np.all(amplitude_mask.XX == field.XX) and np.all(amplitude_mask.YY == field.YY), "Spatial dimensions of field and amplitude_mask must be the same"
        modulation_amplitude = amplitude_mask.aperture
    if phase_mask is not None:
        assert type(phase_mask) is Aperture, "phase_mask must be of type Aperture"
        assert np.all(phase_mask.XX == field.XX) and np.all(phase_mask.YY == field.YY), "Spatial dimensions of field and phase_mask must be the same"
        modulation_phase = phase_mask.aperture

    # Creates a new empty field to store the modulated field
    modulated_field = create_empty_field_from_field(field)

    # Calculates the modulation field from the provided masks
    modulation = modulation_amplitude *np.exp(1.0j*modulation_phase)

    # Modulates the input field
    modulated_field.field = field.field*modulation
    

    return modulated_field


def generate_uniform_field(field, E0=1 ):
    """
    Generates a unifor, wavefront field with amplitude E
    
    Args:
        :field: Input field object to add the flat field
        :E0: Amplitude of the electric field
    Returns:
        :field: returns the field.
    """
    assert type(field) is Field, "field must be of type Field"

    field.field = np.ones(field.XX.shape)*E0
    

    return field

def generate_gaussian_field(field, E0, w0, center=(0,0) ):
    """
    Generates a Gaussian beam amplitude (no imaginary part) 
    **to be deprecated in favor of generate_gaussian_beam

    E = E0*exp(-r^2/w0^2)
    E0 - amplitude on axis
    r - radius
    w0 - radius where amplitude is 1/e of it's on axis value
    
    Args:
        :field: Input field object to add the flat field
        :E: Amplitude of the electric field
    Returns:
        :field: returns the field.
    """
    assert type(field) is Field, "field must be of type Field. Field is type %s"%(type(field))


    x0,y0 = center

    field.field = E0*np.exp(-((field.XX-x0)**2+(field.YY-y0)**2)/(w0**2))


    return field



def generate_gaussian_beam(field, w0, z, wavelength, center=(0,0), E0=1 ):
    """
    Generates a Gaussian beam (with wimaginary part) https://spie.org/publications/spie-publication-resources/optipedia-free-optics-information/fg12_p18-19_gaussian_beams#_=_
    
    
    Args:
        :field: Input field object to add the flat field
        :w0: Beam waist at z=0 
        :z: Z-coordinate, with origin in minimum waist
        :wavelength: Wavelength in meters 
        :center: center in (x,y) plane, default (0,0) 
        :E0: Amplitude of the electric field
    Returns:
        :field: returns the field.
    """
    assert type(field) is Field, "field must be of type Field. Field is type %s"%(type(field))

    x0,y0 = center

    k = 2*np.pi/wavelength
    r2 = ((field.XX-x0)**2+(field.YY-y0)**2)
    
    zR = np.pi*w0**2/wavelength 
    wz = w0*np.sqrt(1+(z/zR)**2)
    Rz = z*(1+(zR/z)**2)
    phi = k*z - np.arctan(z/zR)+ k* r2/(2*Rz)
    
    Exyz = E0 * np.exp(-r2/wz**2) * np.exp(1.0j * phi)
    
    field.field = Exyz

    return field





class Screen:
    """
    Class Screen:
        Creates a Screen object that is 1D or 2D and has internal xyz coordinates
         to facilitate the propagation of a field onto a target screen
          
        
        
    Args:
        :x:         Vector for the x axis
        :y:         Vector for the y axis
        :z:         Vector for the z axis
    
    Methods:
        :screen:  returns the field on the screen
        :shape:     returns the shape of the field

    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.XX, self.YY, self.ZZ = np.meshgrid(x, y, z)#, indexing='ij')
        # self.pixel_x = self.x[1]-self.x[0]
        # self.pixel_y = self.y[1]-self.y[0]
        # self.pixel_z = self.z[1]-self.z[0]
        
        self.n = np.ones(self.XX.shape)
    
        self.screen = np.zeros(self.XX.shape, dtype=complex)
    @property
    def shape(self):
        return self.screen.shape
    @property
    def amplitude(self):
        return np.abs(self.screen)
    @property
    def phase(self):
        return np.angle(self.screen)
    @property 
    def intensity(self):
        return np.abs(self.screen)**2
    @property
    def nindex(self):
        return self.n
    

    
    
def create_screen_XY(xmin, xmax, N_x, ymin, ymax, N_y, z):
    """
    Creates an empty screen of the mesh dimensions provided
    
    Args: 
        :xmin, xmax:    range for x 
        :N_x:           number of x points
        :ymin, ymax:    range for y 
        :N_y:           number of y points
        :z:             z position of the screen plane
    
    Returns:
        :screen: empty Screen
    """
    x = np.linspace(xmin, xmax, N_x)
    y = np.linspace(ymin, ymax, N_y)
    z=z
    
    return Screen(x,y,z)



    
def create_screen_YZ(ymin, ymax, N_y, zmin, zmax, N_z, x=0):
    """
    Creates an empty screen of the mesh dimensions provided
    
    Args: 
        :ymin, ymax:    range for y 
        :N_y:           number of y points
        :zmin, zmax:    range for z
        :N_z:           number of z points
        :x:             x position of the screen plane
    
    Returns:
        :screen: empty Screen
    """
    x=x
    y = np.linspace(ymin, ymax, N_y)
    z = np.linspace(zmin, zmax, N_z)

    return Screen(x,y,z)


    
def create_screen_ZZ(zmin, zmax, N_z, x=0, y=0, log=False):
    """
    Creates an empty screen of the mesh dimensions provided
    
    Args: 
        :zmin, zmax:    range for z
        :N_z:           number of z points
        :x:             x position of the screen line
        :y:             y position of the screen line
    
    Returns:
        :screen: empty Screen
    """
    z = np.linspace(zmin, zmax, N_z)
    
    if log==True: 
        z = np.logspace(zmin, zmax, N_z)
    y=y
    x=x
    
    return Screen(x,y,z)