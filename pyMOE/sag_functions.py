"""
sag_functions.py 
Module containing several sag or phase functions for various optical elements


"""

import numpy as np
# from zernike import RZern


def fresnel_lens_phase(XX,YY,focal_length,wavelength): 
    """
    returns the COMPLEX PHASE of a fresnel lens with input meshgrid (x,y) with center at (x0,y0)
    
    Args:
        :XX:            x array from meshgrid 
        :YY:            y array from meshgrid 
        :focal_length:  focal distance 
        :wavelength:    wavelength of design
    
    Note: for angle (in rad), call numpy.angle(...)
    """

    rc = np.sqrt((XX)**2 + (YY)**2)
    fresn = np.exp(1.0j*(focal_length-np.sqrt(focal_length**2 + rc**2))*(2*np.pi)/(wavelength))
    fresn = np.angle(fresn)
    fresn = fresn-np.min(fresn)
    
    return fresn     




def spiral(x,y,L):
    """
    returns a spiral COMPLEX PHASE with input meshgrid (x,y) with center at (x0,y0)
    
    Args:
        :x:  x array from meshgrid 
        :y:  y array from meshgrid 
        :L:  topological charge 
        
    Returns
        spiral phase
    """

    theta = np.arctan2(y, x)
    sp = np.exp(1.0j*L*theta)
    sp = np.angle(sp)
    return sp


def saddle(x,y,a,b):
    """
    returns a COMPLEX PHASE saddle function 
    
    Args:
        :x:   x array from meshgrid 
        :y:   y array from meshgrid 
        :a:   arbitrary parameter
        :b:   arbitrary parameter 
    """

    sfunc =  (a * ((x*x - y*y)) -b) 
    func = np.exp(1.0j*sfunc)
    func = np.angle(func)

    return func


def monkey_saddle(x,y,a,b):
    """
    returns a COMPLEX PHASE monkey saddle function 
    
    Args:
        :x:   x array from meshgrid 
        :y:   y array from meshgrid 
        :a:   arbitrary parameter
        :b:   arbitrary parameter 
    """
    
    sfunc =  (a * ((x*x*x- 3*x*y*y)*1e6) -b) 
    func = np.exp(1.0j*sfunc)
    func = np.angle(func)

    return func


def Alvarez_phase(XX,YY, f1, f2, tuning_distance, wavelength): 
    """
    returns the phase of an Alvarez lens profile
    
    Args:
        :XX:            x array from meshgrid 
        :YY:            y array from meshgrid 
        :f1:  shorter focal distance f1<f2
        :f2:  longer focal distance f1<f2
        :tuning_distance: tuning displacement range for the Alvarez lens
        :wavelength:    wavelength of design
    
    Note: 
    """

    A = np.pi/(2*tuning_distance*wavelength)*(1/f1-1/f2)

    phase = A*(XX*np.power(YY,2)+np.power(XX,3)/3)

    return phase







