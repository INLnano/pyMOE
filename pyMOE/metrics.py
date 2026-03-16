"""
metrics.py 
Module containing some metric functions 

""" 

import numpy as np
import scipy.fftpack as sfft

import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D

from scipy.signal import find_peaks, peak_widths 
from scipy.special import j0, j1

from pyMOE.utils import simpson2d, qmc_2d


def calculate_OTF(intensity, plot_OTF=False, optimize=False):    
    """
    Calculates the OTF at an observation screen 
    
    Args: 
        :intensity: 2D intensity of the field at the observation screen 
        :plot_OTF:  plots the abs of the OTF with imshow (no freq information) for quick inspection only
        :optimize:  If True will import jax for the gradient descent based optimization routine, defaults False 
    Returns: 
        :otf:    OTF (2D complex array) of the field at the screen 
    
    """ 
    if optimize==True: 
        #print("OPTIMIZE")
        import jax.numpy.fft as sfft
        import jax.numpy as np
    else: 
        import numpy as np
        import scipy.fftpack as sfft
    
    #intensity = np.abs(Escreen)**2
    norm_intensity = intensity/np.sum(intensity)  # normalize intensity, not field

    otf = sfft.fftshift(sfft.fft2(sfft.ifftshift(norm_intensity)))

    center = (otf.shape[0] // 2, otf.shape[1] // 2)
    otf_norm =  otf/np.abs(otf[center])
    
    if plot_OTF==True: 
        fig = plt.figure() 

        plt.imshow( np.abs(otf_norm) )
        plt.colorbar()
        
        plt.show()
        
        fig.savefig("OTF.png")
    
    return otf_norm

    
def MTF_from_OTF(OTF, optimize=False):
    """
    Create the MTF from the OTF (MTF adim.)
    
    Args: 
        :OTF:  2D arrays with the (complex) OTF 
        
    Returns: 
        :MTF, MTF_x, MTF_y: 2D array with MTF, MTF slice in (positive) x axis, MTF slice in (positive) y axis
    """
    if optimize==True: 
        import jax.numpy as np
    else: 
        import numpy as np 
        
    MTF = np.abs(OTF)
    MTF = MTF/np.max(MTF)
    
    m,n = MTF.shape
    
    MTF_x = MTF[int(m/2):, int(n/2)]
    MTF_y = MTF[int(m/2), int(n/2):]
    
    return MTF, MTF_x, MTF_y



def PTF_from_OTF(OTF, optimize=False):
    """
    Create the PTF from the OTF (PTF in radians)
    
    Args: 
        :OTF:  2D arrays with the (complex) OTF 
        
    Returns: 
        :PTF, PTF_x, PTF_y: 2D array with PTF, PTF slice in (positive) x axis, PTF slice in (positive) y axis
  
    """
    if optimize==True: 
        import jax.numpy as np
    else: 
        import numpy as np 
        
    PTF = np.angle(OTF)
    
    m,n = PTF.shape
    
    PTF_x = PTF[int(m/2):, int(n/2)]
    PTF_y = PTF[int(m/2), int(n/2):]
    
    return PTF, PTF_x, PTF_y
    
    
def spatial_freq_space(obj): 
    """
    Creates a real 2D space into a 2D reciprocal space (e.g. coordinate -> spatial frequency space)
    Can be used for Screen, Field, and Aperture objects 
    
    Args: 
        :obj:   Object of type Screen, Field or Aperture 
        
    Returns: 
        :fx,fy,FX,FY: (fx,fy) arrays with freqs in (x,y) and meshgrids (FX,FY)
    """
    Nx, Ny = obj.x.size, obj.y.size
    dx, dy = obj.pixel_x, obj.pixel_y
    
    fx = sfft.fftshift(sfft.fftfreq(Nx, d=dx) )
    fy = sfft.fftshift(sfft.fftfreq(Ny, d=dy) )
    
    FX, FY = np.meshgrid(fx, fy)
    
    return fx,fy,FX,FY
    
   
def theoretical_ideal_OTF(screen, D, wavelength, z, optimize=False):
    """
    Theoretical OTF for circular aperture of diameter D 
    From Goodman (6-32) 
    
    Args: 
        :screen:     Screen object where to compute the (ideal) OTF
        :D:          Aperture diameter (larger one)
        :wavelength: Wavelength in m 
        :z:          Distance 
    
    Returns: 
        :OTF:  2D array with the OTF 
    """
    if optimize==True: 
        import jax.numpy as np
    else: 
        import numpy as np 
    
    fx, fy, FX, FY = screen.fx, screen.fy, screen.FX, screen.FY
    
    rho = np.sqrt(FX**2 + FY**2)

    #NA = np.sin(np.arctan(D/2/z)) 
    #fc = 2 * NA/wavelength  
    fc = D/(wavelength * z) 
    
    OTF = (2 / np.pi) * (np.arccos(rho / fc) - (rho / fc) * np.sqrt(1 - (rho / fc)**2))
    
    #OTF[rho > fc] = 0.0
    
    OTF = np.nan_to_num(OTF, nan=0)
    
    return OTF
 
def strehl_ratio(intensity, screen, D, wavelength, z , strehl_from_mtf=False, plot_MTF=False, plot_OTF=False, optimize=False): 
    """ 
    Calculates the Strehl ratio, based on the OTF volumes (Goodman prob. 6-9)
    Assumes Circular Aperture as theoretical. 
    
    E.g. for a Screen, strehl_ratio(np.abs(screen.screen[:,:,-1])**2, screen, max([max(screen.x), max(screen.y)])*2, wavelength, z = screen.z[-1]) 
   
    Args: 
        :intensity:       Intensity of the 2D array at the XY screen 
        :screen:          Screen object used for definition of equivalent ideal (theoretical) circular aperture
        :D:               Aperture diameter (rule of thumb, choose the larger one between x and y )
        :wavelength:      Wavelength in m 
        :z:               Distance to xy screen in m 
        :strehl_from_mtf: (default False) Calculate the Strehl from the MTF (abs(OTF)), although actual definition is from the OTF, not MTF 
        :plot_MTF:        (default False) Plot the MTF, for inspection
        :plot_OTF:        (defaul False) Plot the OTF projecction on positive fx axis   
    Returns: 
        :strehl:          Strehl ratio (value)
    """

    if optimize==True: 
        import jax.numpy as np
        from jax.scipy.integrate import trapezoid 
    else: 
        import numpy as np 
        from scipy.integrate import trapezoid 
        
    OTF_exp = calculate_OTF(intensity, optimize=optimize)
    MTFexp, MTFexpx, MTFexpy = MTF_from_OTF(OTF_exp, optimize=optimize)

    fx, fy = screen.fx, screen.fy
    
    if strehl_from_mtf==True: 
        volume_exp = trapezoid(trapezoid(np.abs(OTF_exp), fy, axis=1), fx)
    else: 
        volume_exp = trapezoid(trapezoid(OTF_exp, fy, axis=1), fx)
    
    OTF_ideal = theoretical_ideal_OTF(screen, D, wavelength, z , optimize=optimize)
    MTFid, MTFidx, MTFidy = MTF_from_OTF(OTF_ideal, optimize=optimize)
    
    if strehl_from_mtf==True: 
        volume_ideal = trapezoid(trapezoid(np.abs(OTF_ideal), fy, axis=1), fx)
    else: 
        volume_ideal = trapezoid(trapezoid((OTF_ideal), fy, axis=1), fx)

    strehl = np.abs(volume_exp)/np.abs(volume_ideal)
        
    if plot_MTF==True:     
        fig = plt.figure() 
        
        plt.plot(fx[fx>=0]*1e-3, MTFidx, '-', label="MTFx Ideal")
        plt.plot(fx[fx>=0]*1e-3, MTFexpx, label="MTFx Experimental")
        plt.plot(fy[fy>=0]*1e-3, MTFidy, ':', label="MTFy Ideal")
        plt.plot(fy[fy>=0]*1e-3, MTFexpy, ':', label="MTFy Experimental")
        plt.xlabel("Freq. (cycles/mm)")
        plt.ylabel("MTF")
        plt.title("MTF plot (Strehl Ratio ="+str(strehl)+")")

        plt.legend()
        
        fig.savefig("MTFs.png") 
    
    if plot_OTF==True:      
        fig, ax1 = plt.subplots()

        # Plot PSF
        color1 = 'C0'
        ax1.set_xlabel('Freq. (cycles/mm)')
        ax1.set_ylabel('Real(OTF)', color=color1)
        ax1.plot(fx[fx>=0]*1e-3, np.real(OTF_ideal)[fx>=0][:,int(len(fx)/2)] , color=color1, label='Ideal')
        ax1.plot(fx[fx>=0]*1e-3, np.real(OTF_exp)[fx>=0][:,int(len(fx)/2)] , '--', color=color1, label='Experimental')
        ax1.tick_params(axis='y', labelcolor=color1)

        # Create a second y-axis for the OTF
        ax2 = ax1.twinx()
        color2 = 'C1'
        ax2.set_ylabel('Imag(OTF)', color=color2)
        ax2.plot(fx[fx>=0]*1e-3, np.imag(OTF_ideal)[fx>=0][:,int(len(fx)/2)] , color=color2, label='Ideal')
        ax2.plot(fx[fx>=0]*1e-3, np.imag(OTF_exp)[fx>=0][:,int(len(fx)/2)] , '--', color=color2, label='Experimental')
        ax2.tick_params(axis='y', labelcolor=color2)
        
        colors = ['black', 'black']
        lines = [Line2D([0], [0], color='black', linestyle='-'), Line2D([0], [0], color='black', linestyle='--') ]
        labels = ['Ideal', 'Experimental']
        plt.legend(lines, labels)
        #plt.legend()

        fig.suptitle('OTF')
        fig.tight_layout()
        fig.savefig("OTFs.png") 

    return strehl


def intensity_theo_Airy(screen, D, wavelength, z = 1):
    """
    Airy pattern (theoretical solution to circular aperture) on a (x,y) screen grid 
    
    Args: 
        :screen:     Screen object for definition of equivalent ideal circular aperture
        :D:          Aperture diameter (larger one)
        :wavelength: Wavelength in m 
        :z:          Distance 
        
    Returns: 
        :I:          Intensity (2D array)
    """
    #if len(screen.z) !=1: 
    XX, YY = screen.XX[:,:,-1], screen.YY[:,:,-1]
    r = np.sqrt(XX**2 + YY**2)
    z = screen.z[-1]
    #print("Picking up the last z propagated assumed to be the screen position")
        
    
    u = np.pi*D*r/(wavelength*z)
    
    I = (2*j1(u)/u)**2
    
    return I
    


def encircled_energy(screen, p, center=(0,0), optimize=False):
    """
    Calculate the encircled energy 
    
    Args: 
        :screen:   Screen in XY, where to calculate the encircled energy 
        :p:        Radius to encircle the energy
        :center:   Center of the circle, default = (0,0)
        
    Returns: 
        :energy:   Encircled energy 
    """
    if optimize==True: 
        import jax.numpy as np
    else: 
        import numpy as np
        
    #Intensity (ideally in W/m2) 
    intensity = np.abs(screen.screen)**2 
    
    #Circular mask 
    mask= np.where(((screen.XX-center[0])**2 + (screen.YY-center[1])**2 ) <= p**2, 1, 0)
    
    #Masked intensity 
    masked = intensity*mask
    
    areapx = screen.pixel_x * screen.pixel_y 
    
    #Calculate the integral   
    energy = masked.sum()*areapx
    
    return energy
    
    
def theoretical_encircled_energy(screen, p):
    """
    Theoretical encircled energy for a circular aperture, see e.g. 2015 Andersen Eq (6)
    
    Args: 
        :screen:   Screen in XY, where to calculate the encircled energy 
        :p:        Radius from the center 
        
    Returns: 
        :energy:   Encircled energy 
    """
    
    rmax = np.max(screen.x)
    
    res = 1 - j0(p/rmax)**2 - j1(p/rmax)**2
    
    return res
    

def ensquared_energy(screen, px, py, center=(0,0), optimize=False): 
    """
    Calculate the ensquared energy 
    
    Args: 
        :screen:   Screen in XY, where to calculate the encircled energy 
        :p:        Radius to encircle the energy
        :center:   Center of the circle, default = (0,0)
        
    Returns: 
        :energy:   Encircled energy 
    """
    if optimize==True: 
        import jax.numpy as np
    else: 
        import numpy as np
        
    #Intensity (~W/m2) 
    intensity = np.abs(screen.screen)**2 
    
    #Circular mask 
    mask1= np.where((np.abs(screen.XX-center[0]) <= px ), 1, 0)
    mask2= np.where((np.abs(screen.YY-center[0]) <= py ), 1, 0)
    mask = mask1 * mask2
    
    #Masked intensity 
    masked = intensity*mask
    
    areapx = screen.pixel_x * screen.pixel_y 
    
    #Calculate the integral   
    energy = masked.sum()*areapx
    
    return energy
    
    
def focus_efficiency(screen, p, input_energy, center=(0,0), optimize=False): 
    """
    Focus efficiency 
    
    Args: 
        :screen:       Screen in XY, where to calculate the encircled energy 
        :p:            Radius to encircle the energy
        :center:       Center of the circle, default = (0,0)
        :input_energy: Energy of field previous to aperture 
        
    Returns: 
        :energy:   Encircled energy 
        
    
    Notes: For a rule of thumb, some works use 1*FWHM and 3*FWHM of the central diffraction peak radius,
    but following Engelberg et al. https://doi.org/10.1038/s41566-022-00963-7 and https://doi.org/10.1515/nanoph-2022-0196 
    the actual focus efficiency should be calculated for a radius where the PSF flattens out. 
    For Meem et al. (https://doi.org/10.1364/OPTICA.388697) this happens at 10*FWHM. 
    """
    
    encircled_energy_screen = encircled_energy(screen, p, center=(0,0), optimize=optimize)
    result = encircled_energy_screen/input_energy
    
    return result 



def find_FWHM_2d(XY, prominence=0.5): 
    """
    Find FWHM of largest (identified) peak in XY data 
    
    Args: 
        :XY:         2D mesh data e.g. XY = screen.screen[:,:,-1]
        :prominence: Prominence of the peak for find_peaks function 

    Return: 
        :largestpeak_height, largestpeak_width:  Height of the largest peak, FWHM of that peak 
    
    """ 
    all_maxpeakwidths = []
    all_maxpeakheights = []
    all_peak_positions = []

    rows, cols = XY.shape  

    for y in range(rows):
        profile = XY[y, :]  # horizontal line (x-direction at y)
        peak_ndxs, _ = find_peaks(profile, prominence=prominence)

        fwhm, heights, left_ips, right_ips = peak_widths(profile, peak_ndxs, rel_height=0.9)

        max_idx = np.argmax(heights)
        x = peak_ndxs[max_idx]

        all_maxpeakheights.append(heights[max_idx])
        all_maxpeakwidths.append(fwhm[max_idx])
        all_peak_positions.append((x, y))  

    all_maxpeakheights = np.array(all_maxpeakheights)
    all_maxpeakwidths = np.array(all_maxpeakwidths)
    all_peak_positions = np.array(all_peak_positions)
    
    max_idx = np.argmax(all_maxpeakheights)
    largestpeak_height = all_maxpeakheights[max_idx]
    largestpeak_width = all_maxpeakwidths[max_idx]
    largestpeak_position = all_peak_positions[max_idx]  # (x, y) in array indices

    return largestpeak_height, largestpeak_width, largestpeak_position



def Engelberg_OPM(screen, p, input_energy, transmission , wavelength, center=(0,0), intensity = None, z = None, D = None, optimize=False ):
    """
    Overall performance metric (OPM) from eq (2) defined at Engelberg et al. https://doi.org/10.1515/nanoph-2022-0196 
    
    Args: 
        :screen:       Screen where to calculate the encircled energy (last XY plane at z_focal=screen.z[-1])
        :p:            Radius to encircle the energy
        :input_energy: Energy of field previous to aperture 
        :transmission: Overall transmission [might be = eta, = 1]
        :wavelength:   Wavelength in m
        :center:       Center of the circle, default = (0,0)
        :intensity:    2D array with the intensity at the screen (observation) plane, by default the last XY plane screen.screen[:,:,-1]
        :z:            Distance to the observation plane, by default the last z value
        :D:            Size of the aperture, by default max([max(screen.x), max(screen.y)])
        
    Returns: 
        :OPM:          Calculated OPM
    
    """
    if optimize==True: 
        import jax.numpy as np
    else: 
        import numpy as np
        
    if intensity is None: 
        intensity = np.abs(screen.screen[:,:,-1])**2
    if z is None: 
        z = screen.z[-1]
    if D is None: 
        D = max([max(screen.x), max(screen.y)])*2
    
    #print(D,z, transmission, wavelength)
    
    eta = focus_efficiency(screen, p, input_energy, optimize=optimize)
    #print(eta)
    T = transmission
    strehl = strehl_ratio(intensity, screen, D, wavelength, z, optimize=optimize) 
    #print(strehl)
    
    OPM = eta*strehl/np.sqrt(transmission)
   
    return OPM



def Engelberg_EOPM(screen, p, input_energy, transmission , wavelength, center=(0,0), intensity = None, z = None, D = None ):
    """
    Overall performance metric (OPM) from eq (2) defined at Engelberg et al. https://doi.org/10.1515/nanoph-2022-0196 
    
    Args: 
        :screen:       Screen where to calculate the encircled energy (last XY plane at z_focal=screen.z[-1])
        :p:            Radius to encircle the energy
        :input_energy: Energy of field previous to aperture 
        :transmission: Overall transmission [might be = eta, = 1]
        :wavelength:   Wavelength in m
        :center:       Center of the circle, default = (0,0)
        :intensity:    2D array with the intensity at the screen (observation) plane, by default the last XY plane screen.screen[:,:,-1]
        :z:            Distance to the observation plane, by default the last z value
        :D:            Size of the aperture, by default max([max(screen.x), max(screen.y)])
        
    Returns: 
        :EOPM:          Calculated EOPM
    
    """

   
    return 