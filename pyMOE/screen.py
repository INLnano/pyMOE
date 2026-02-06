"""
screen.py


Definition of Class Screen and related functions


"""

import numpy as np

from pyMOE.aperture import Aperture

import scipy.fftpack as sfft 


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
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        
        self.Nx, self.Ny, self.Nz = self.x.size, self.y.size, self.z.size
        
        self.XX, self.YY, self.ZZ = np.meshgrid(x, y, z, indexing='ij')
    
        if self.Nx!=1: 
            self.pixel_x = np.diff(self.x)[0]
            self.fx = sfft.fftshift(sfft.fftfreq(self.Nx, d=self.pixel_x) )
        if self.Ny!=1:        
            self.pixel_y = np.diff(self.y)[0]
            self.fy = sfft.fftshift(sfft.fftfreq(self.Ny, d=self.pixel_y) )
        if self.Nz!=1:        
            self.pixel_z = np.diff(self.z)[0]
        
        self.n = np.ones(self.XX.shape)
        
        if (self.Nx !=1) and (self.Ny!=1):
            self.FX, self.FY = np.meshgrid(self.fx, self.fy, indexing='ij')
        
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
    
    def slice_XY(self, z):
        z_index = np.argmin(np.abs(self.z - z))
        return self.screen[:,:,z_index]
    def slice_XZ(self, y):
        y_index = np.argmin(np.abs(self.y - y))
        return self.screen[:,y_index,:]
    def slice_YZ(self, x):
        x_index = np.argmin(np.abs(self.x - x))
        return self.screen[x_index,:,:]
    
    def slice_X(self, y, z):
        y_index = np.argmin(np.abs(self.y - y))
        z_index = np.argmin(np.abs(self.z - z))
        return self.screen[:,y_index,z_index]
    def slice_Y(self, x, z):
        x_index = np.argmin(np.abs(self.x - x))
        z_index = np.argmin(np.abs(self.z - z))
        return self.screen[x_index,:,z_index]
    def slice_Z(self, x, y):
        x_index = np.argmin(np.abs(self.x - x))
        y_index = np.argmin(np.abs(self.y - y))
        return self.screen[x_index,y_index,:]

    
    
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
    z=  np.array([z])
    
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
    x = np.array([x])
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
    y = np.array([y])
    x = np.array([x])
    
    return Screen(x,y,z)