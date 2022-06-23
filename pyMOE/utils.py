
import time
from datetime import timedelta
import numpy as np

def progress_bar(progress, bar_length=20, bar_character='#'):
    """
    Progress bar.
    Writes a progress bar in place in the output
    
    Args:
        progress: value between 0 and 1
        bar_length: number of characters to consider in the bar
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
        array : input array of values
        levels : integer number of levels to consider or array of levels
        
    Returns:
        bins: bins corresponding to the levels
        digitized: digitized array
        
    To do:
        Consider the midpoint selection in the future
    """    
    assert isinstance(levels, (np.ndarray, int)), "levels must be a scalar or numpy array"
    if isinstance(levels, int):
        bins = np.linspace(array.min(), array.max() , levels, endpoint=False)
    else:
        bins = levels
    
    dig = np.digitize(array, bins, )
    
    # Everything below the minimum bin level is changed to the minimum level
    dig[dig==0] = 1
    dig = dig-1
    return bins, dig


def discretize_array(array, levels):
    bins, dig = digitize_array_to_bins(array, levels)
    
    return bins[dig]