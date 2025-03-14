
"""
utils.py 
Module containing utility functions 

""" 


import time
from datetime import timedelta
import numpy as np

def progress_bar(progress, bar_length=20, bar_character='#'):
    """
    Progress bar.
    Writes a progress bar in place in the output
    
    Args:
        :progress: value between 0 and 1
        :bar_length: number of characters to consider in the bar
    """
    
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
    block = int(round(bar_length * progress))
#     clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( bar_character * block + "-" * (bar_length - block), progress * 100)
    if progress <1:
        end = '\r'
    else:
        end = "\n"
    print(text, end=end) 



class Timer(object):
    """
    Timer helper class to calculated elapsed time of chunk of code, from https://stackoverflow.com/a/5849861/7996766
    """
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name,)
        print('Elapsed: %s' % str(timedelta(seconds=(time.time() - self.tstart))))



def mean_squared_error(new, original):
    """ Calculates mean squared error between the new array and the original"""
    return (np.square(new-original)).mean()


def create_levels(start, end, levels):
    """ Creates linear levels list with N"""
    return np.linspace(start, end, levels, endpoint=False)

def digitize_array_to_bins(array, levels):
    """Digitizes the given array to within the number of levels provided 
    
    Args:
        :array:  input array of values
        :levels: integer number of levels to consider or array of levels
        
    Returns:
        :bins:      bins corresponding to the levels
        :digitized: digitized array
        
    To do:
        Consider the midpoint selection in the future
    """    
    assert isinstance(levels, (np.ndarray, int)), "levels must be a scalar or numpy array"
    if isinstance(levels, int):
        bins = np.linspace(np.nanmin(array), np.nanmax(array) , levels, endpoint=False) 
    else:
        bins = levels
    
    dig = np.digitize(array, bins, )
    
    # Everything below the minimum bin level is changed to the minimum level
    dig[dig==0] = 1
    dig = dig-1
    # dig[np.isnan(array)] = np.nan
    return bins, dig


def discretize_array(array, levels):
    bins, dig = digitize_array_to_bins(array, levels)
    
    return bins[dig]
    
def simpson2d(f,ax,bx,ay,by):
    """
    Implements Simpson method for calculating a double integral in 2D array f
    
    Arguments: 
        :f:         2D array to calculate integral 
        :[ax, bx]:  limits of integration in x, [lower, upper]
        :[ay, by]:  limits of integration in y, [lower, upper]
    """
    
    num = len(f)
    hx = (bx-ax)/(num-1)
    hy = (by-ay)/(num-1)
    h = hx * hy / 9

    # Simpson coefficients 
    #1 4 2 4 ...2 4 1
    sc = 2*np.ones(num)
    sc[np.arange(1,num-1,2)] = 4
    sc[0] = 1
    sc[num-1] = 1
    #print(sc)

    scx = np.meshgrid(sc,sc)[0]
    scxy = np.ones((num,num))
    scxy[np.arange(1,num-1,2,dtype=int),:] = scx[np.arange(1,num-1,2,dtype=int),:]*sc[1]
    scxy[np.arange(2,num-2,2, dtype=int),:] = scx[np.arange(2,num-2,2,dtype=int),:]*sc[2]
    scxy[0,:] = sc
    scxy[num-1,:] = sc
    
    #print(scxy)

    # integral    
    tint = h * np.sum(np.sum(scxy * f))
    
    return tint



def find_closest_indices(x, y):
    from scipy.spatial.distance import cdist

    # Reshape y to a column vector
    y = y.reshape(-1, 1)

    # Calculate the distances between each value of y and x
    distances = cdist(y, x.reshape(-1, 1))

    # Find the index of the closest value in x for each value of y
    closest_indices = np.argmin(distances, axis=1)

    return closest_indices
