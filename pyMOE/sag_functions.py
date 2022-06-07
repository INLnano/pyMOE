"""
Module containing several sag or phase functions for various optical elements


"""

import numpy as np

def fresnel_lens_phase(XX,YY,focal_length,wavelength): 
    """
    returns the COMPLEX PHASE of a fresnel lens with input meshgrid (x,y) with center at (x0,y0)
    Args:
        XX:  x array from meshgrid 
        YY: y array from meshgrid 
        focal_length: focal distance 
        wavelength: wavelength of design
    
    Note: for angle (in rad), call numpy.angle(...)
    """

    rc = np.sqrt((XX)**2 + (YY)**2)
    fresn = np.exp(1.0j*(focal_length-np.sqrt(focal_length**2 + rc**2))*(2*np.pi)/(wavelength))
    return fresn     




def spiral(x,y,L):
    """
    returns a spiral COMPLEX PHASE with input meshgrid (x,y) with center at (x0,y0)
    Args:
        x = x array from meshgrid 
        y = y array from meshgrid 
        L = topological charge 
        
    Returns
        spiral phase
    """

    theta = np.arctan2(y, x)
    sp = np.exp(1.0j*L*theta)
    return sp


def saddle(x,y,a,b):
    """
    returns a COMPLEX function (the function to implement NEEDS to be complex to extract the phase!)
    Args:
        x = x array from meshgrid 
        y = y array from meshgrid 
        a = arbitrary parameter
        b  = arbitrary parameter 
    """

    sfunc =  (a * ((x*x - y*y)) -b) 
    func = np.exp(1.0j*sfunc)

    return func




