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





### Dammann gratings

def alternate_transitions(x, transitions, start_value=1):
    y = np.full_like(x, start_value)
    for i in range(len(transitions)-1):
        mask = np.logical_and(x >= transitions[i], x < transitions[i+1])
        y[mask] = start_value * (-1) ** (i+1)
    mask = x >= transitions[-1]  # Add this line to handle the last transition
    y[mask] = start_value * (-1) ** len(transitions)
    return y



def alternate_transitions_symmetric(x, transitions, start_value=1):
    assert np.all(np.diff(transitions) > 0), "transitions should be monotonic increasing"
    y = np.full_like(x, start_value)
    for i in range(len(transitions)-1):
        mask = np.logical_and(abs(x) >= transitions[i], abs(x) < transitions[i+1])
        y[mask] = start_value * (-1) ** (i+1)
    mask = abs(x) >= transitions[-1]  # Add this line to handle the last transition
    y[mask] = start_value * (-1) ** len(transitions)
    return y



    

def dammann_grating_element(x,transitions, start_value=1):
    assert np.all(np.array(transitions)<=0.5), "transitions should be <= 0.5"

    y = alternate_transitions_symmetric(x,transitions=transitions, start_value=start_value)
    return y



def dammann_grating_periodic(x,transitions, period=1, start_value=1):
    # assert np.all(np.array(transitions)<=0.5), "transitions should be <= 0.5"

    transitions = np.array(transitions)*period
    # x = clip_remainder(x, period)
    half_period = period/2
    x = np.mod(np.abs(x+half_period), period)-half_period

    y = alternate_transitions_symmetric(x,transitions=transitions, start_value=start_value)
    return y


    
def dammann_2d(XX,YY,transitions_x=None, period_x=1, transitions_y=None, period_y=1, start_value_x=1, start_value_y=1):
    x = XX[0,:]
    y = YY[:,0]

    if transitions_x is not None:
        grating_x = dammann_grating_periodic(x, transitions_x, period=period_x, start_value=start_value_x)
    else:
        grating_x = np.ones_like(x)
    
    if transitions_y is not None:
        grating_y = dammann_grating_periodic(y, transitions_y, period=period_y, start_value=start_value_y)
    else:  
        grating_y = np.ones_like(y)

    grating = np.outer(grating_y, grating_x)
    return grating



