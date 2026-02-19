"""
propagate.py   

This module had origin in the "propopt" repository module with same name 
-> For more information, examples of use and validation tests, 
   please see link for repository : https://github.com/cunhaJ/propopt 

Here, "propopt" code has been modified and extended for use with apertures and masks of "pyMOE"
"""

import numpy as np 
import scipy 
import scipy.fft as sfft 
from scipy import integrate

import decimal

from pyMOE.utils import simpson2d, qmc_2d


import dask
import dask.array as da
from dask import delayed
import dask.bag as db
from dask.diagnostics import ProgressBar

from pyMOE.utils import progress_bar, Timer
from pyMOE.field import Field, Screen, create_screen_XY, create_screen_YZ, modulate_field

from pyMOE.plotting import plot_field 

import matplotlib.pyplot as plt 



def circ_zz24(aperture_rad, zdist, wavelength):
    """
    This function is the axial intensity for a circular aperture following 1992 JOSA 9(2) paper "Diffraction by a circular aperture: a generalization of Fresnel diffraction theory" , exp (24) 
    
    Args: 
        :aperture_rad:  Radius of the circular aperture in m 
        :zdist:         Distance of propagation
        :wavelength:    Wavelength of incoming illumination in m

    Returns: 
        Propagated intensity to zdist 

    """
    import numpy as np 

    k = 2*np.pi /wavelength 

    izz1 = 1/(1+aperture_rad**2/zdist**2)
    izz2 = 2/np.sqrt(1+aperture_rad**2/zdist**2)
    izz3 = (k * aperture_rad**2/(zdist)) /(np.sqrt(1+aperture_rad**2/zdist**2)+1)
    itot = 0.25*(1+izz1-izz2* np.cos(izz3))

    return itot
    
    

def fresnel(z, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength):
    """
    Calculate Fresnel approximation, following Goodman exp 4-17 
    
    Args: 
        :z:             distance to the observation plane in m
        :mask:          2D map (x-y plane) npixel vs npixel mask 
        :npixmask:      number of pixels of the mask 
        :pixsizemask:   size of the pixel at the mask in m
        :npixscreen:    size of the pixel at the screen 
        :dxscreen:      x-size of the screen  in m 
        :dyscreen:      y-size of the screen   in m
        :wavelength:    wavelength in m 
    
    Returns: 
        2D map (x-y plane) with abs of Electric field at distance z 
    
    """ 
    k = 2* np.pi/wavelength
    
    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    ##calculate the resolution  -
    res = wavelength *z/ (pixsizemask *npm)
    
    dmask =   res * npm
    
    if z < Fresnel_criterion(wavelength, dmask/2):
        print("The propagation distance is too short for Fresnel propagation! Propagation results might be incorrect.")
    
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1, ym1)
    
    xs1 = np.linspace(-dxscreen, dxscreen, nps)
    ys1 = np.linspace(-dyscreen, dyscreen, nps)
    (xs, ys) = np.meshgrid(xs1, ys1)
    
    Ef = fresnel_kernel(k, xm, ym, z, mask)
    
    return Ef

def fresnel_kernel(k, xm, ym, z, mask):
    #Goodman, exp 4-17
    v1  = np.exp(1.0j*k* (xm*xm + ym*ym)/ (2*z))
    intarg = v1 * mask
    Ef = sfft.fftshift(sfft.fft2(sfft.ifftshift(intarg)))

    return Ef 

def fraunhofer(z, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength):
    """
    Calculate Fraunhofer approximation, following Goodman exp 4-25
    
    Args: 
        :z:             distance to the observation plane in m
        :mask:          2D map (x-y plane) npixel vs npixel mask 
        :npixmask:      number of pixels of the mask 
        :pixsizemask:   size of the pixel at the mask in m
        :npixscreen:    size of the pixel at the screen 
        :dxscreen:      x-size of the screen  in m 
        :dyscreen:      y-size of the screen   in m
        :wavelength:    wavelength in m 

    Returns: 
        2D map (x-y plane) with abs of Electric field
    
    """
    
    k = 2* np.pi/wavelength
    #number of pixels 
    nps = npixscreen
    npm = npixmask 
    
    dmask = npixmask * npm  
    
    if z < Fraunhofer_criterion(wavelength, dmask/2):
        print("The propagation distance is too short for Fraunhofer propagation! Propagation results might be incorrect.")
    
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1, ym1)
    
    xs1 = np.linspace(-dxscreen, dxscreen, nps)
    ys1 = np.linspace(-dyscreen, dyscreen, nps)
    (xs, ys) = np.meshgrid(xs1, ys1)
    
    delta1=pixsizemask

    Ef = fraunhofer_kernel(k, xm, ym, mask, z, wavelength)

    return Ef


def fraunhofer_kernel(k, xm, ym, mask, z, wavelength): 
    #Goodman, exp 4-25
    v2  = np.exp(1.0j*k* (xs*xs + ys*ys)/ (2*z)) 
    v3  = np.exp(1.0j*k*z)/ (1.0j*wavelength*z)
    Ef = v2 * v3* sfft.fftshift(sfft.fft2(resized))
    return Ef 
    



def Fresnel_num(width, wavelength, zdist):
    """
    Calculation of Fresnel number, Goodman pag 85
    
    Args:
        :width:         size of the aperture in m 
        :wavelength:    wavelength of the aperture in m 
        :zdist:         distance to the screen in m 
    
    Returns: 
        Fresnel number 
    """
    NF = width**2 / (wavelength * zdist)
    return NF 
    

def Fraunhofer_criterion(wavelength, radius):
    """
    Calculation of "Fraunhofer distance" , Goodman 4-27  
    
    Args:
        :wavelength: wavelength of the illumination m 
        :radius: (max) radius of the aperture in m 
    
    Returns: 
        Fraunhofer distance 
    """
    z = (2*np.pi/wavelength)/2*radius**2
    return z

def Fresnel_criterion(wavelength, radius):
    """
    Calculation of "Fresnel distance" , Goodman 
    
    Args:
        :wavelength: wavelength of the illumination m 
        :radius: (max) radius of the aperture in m
    
    Returns: 
        Fresnel distance 
    """
    z = ((np.pi/(4*wavelength))*radius**4)**(1/3)
    return z




@dask.delayed
def kernel_RS(field, k, x,y,z, simp2d=False, method=False, sampler=None):
    """
    Calculates the RS kernel integral from a field input aperture, assumed to be at z=0
    and returns the calculated E field
    
    Implements the Kernel in Mahajan 2011 part II eq 1-20 

    Args:
        :field:     input field
        :k:         Calculated wavenumber k=2pi/(wl*n)
        :x,y,z:     x, y, z coordinates of the screen point being evaluated
        :simp2d:    Defaults False, if True uses the simpson2d function
        :method     Choose between "simp2d", "qmc" or "trap", or defults to scipy simpson 
        :sampler:   Integer sampler for the QMC, e.g. sampler.integers(l_bounds=0, u_bounds=n, n=n) where n is side lenth 
    Returns:
        :E:         Calculated field
    """

    z_field = 0 # the field source is assumed at z=0
    r = np.sqrt( (field.XX-x)**2 + (field.YY-y)**2 + (z_field-z)**2)
    r2 = r*r

    prop1 = np.exp(r*1.0j*k)/r2
    prop2 = z * k/(2*np.pi) *( 1/(r*k) - 1.0j)
    propE = field.field * prop1 * prop2

    # integrate over the input field and return field
    if simp2d==True or method=="simp2d": 
        Exyz = simpson2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1]) /(2*np.pi)
        #print("Simpson 2D method")
    elif method=="qmc":
        if sampler is None: 
            print("Please add a sample selection of points, e.g. Sobol.")
        Exyz = qmc_2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1], sampler) /(2*np.pi)
        #print("QMC 2D method")
    elif method=="trap":
        Exyz = np.trapz(np.trapz(propE, x=field.x, axis=0), x=field.y, axis=0)
    else: 
        Exyz = integrate.simpson(integrate.simpson(propE, field.y),field.x)/(2*np.pi) 

    return Exyz
    

def RS_integral(field, screen, wavelength, n=None, parallel_computing=True, simp2d=False, sampler=None, method=None):
    """
    Calculates the Raleyigh Sommerfeld integral in the  of the first kind (Mahajan 2011 part II eq 1-20), receiving an input field and an observation screen plane on which to 
    calculate the integral.
    
    Args: 
        
        :field:     input Field
        :screen:    Observation Screen
        :wavelength:    wavelength to consider
        :n:         refractive index of the propagation medium (default=1 for vacuum/air)
        :parallel_computing: Flag to trigger the concurrent computation of the kernels using Python Dask library
        :simp2d:    Defaults False, if True uses the simpson2d function
    Returns:
        :screen:    Returns the screen populated with the result
    """

    if (field.pixel_x > wavelength/2) or (field.pixel_y > wavelength/2):
        print("Warning: Sampling field pixel is larger than wavelength/2!")
    k = 2* np.pi/(wavelength)

    xlen,ylen,zlen = screen.XX.shape

    if parallel_computing:
        delayed_tasks = []
        # For each cell on the screen, the RS integral will be calculated based on the input field
        # this loop sets up the delayed tasks to be executed
        for x_i in range(xlen):
            for y_i in range(ylen):
                for z_i in range(zlen):

                    x = screen.XX[x_i, y_i, z_i]
                    y = screen.YY[x_i, y_i, z_i]
                    z = screen.ZZ[x_i, y_i, z_i]
                    
                    if n is not None: 
                        k = 2* np.pi*n[x_i,y_i,z_i]/(wavelength)
                    # the kernel is configured as a dask delayed task
                    result = kernel_RS(field, k ,x,y,z, simp2d, method, sampler)

                    delayed_tasks.append(result)
                    # screen.screen[x_i, y_i, z_i] = a
        
        # the dask.compute triggers the computation of the delayed tasks and stores the result
        # into the results list
        # print(delayed_tasks)
        with ProgressBar():
            results = list(dask.compute(*delayed_tasks))
        # print(results)
        # again we go through the for loop to pop the results and insert it into the screen position
        for x_i in range(xlen):
            for y_i in range(ylen):
                for z_i in range(zlen):
                    screen.screen[x_i, y_i, z_i] = results.pop(0)
    else:
        with Timer():
            for x_i in range(xlen):
                for y_i in range(ylen):
                    for z_i in range(zlen):

                        x = screen.XX[x_i, y_i, z_i]
                        y = screen.YY[x_i, y_i, z_i]
                        z = screen.ZZ[x_i, y_i, z_i]
                        
                        if n is not None: 
                            k = 2* np.pi*n[x_i,y_i,z_i]/(wavelength)
                        
                        result = kernel_RS(field, k ,x,y,z, simp2d, method, sampler).compute()

                        screen.screen[x_i, y_i, z_i] = result
                        progress_bar((x_i*zlen*ylen+y_i*zlen+z_i)/(xlen*ylen*zlen))
            progress_bar(1)

    return screen

    
    
def bluestein_czt(x, f1, f2, fs, mout):
    """
    Bluestein from Hu et al. 2020 
    x =  field.field  * F  
    
    TO COMPLETE
    """
    import numpy as np 
    
    m, n = x.shape #m = x, n = y 
    
    f11 = f1 + (mout * fs + f2 - f1) / (2 * mout)
    f22 = f2 + (mout * fs + f2 - f1) / (2 * mout)
    
    #Eq S15, complex starting point 
    a = np.exp(1j * 2 * np.pi * f11 / fs)
    #Eq S16, complex point spacing
    w = np.exp(-1j * 2 * np.pi / ((mout * fs) /((f22 - f11) )))
    
  
    #calculation of the size of the tile N
    mp = m + mout - 1
    N = 2 ** int(np.ceil(np.log2(mp)))

    #defining the number of vals where to evaluate the czt, from -m +1 to m -1  (2m-1 size)
    h = np.arange(-m + 1, max(mout, m) )

    
    # W^(m^2 /2) 
    h_vals = w ** ((h ** 2) / 2)


    # W^(-m^2 /2) #all m, but cropped to N after 
    hin = 1 / h_vals[0:mp+1]
    
    ### fft ( W^(-m^2 /2) )
    ft = sfft.fft(hin , N, norm="ortho") #if N>len(hin) zero pads the rest, increases the pad with mout 

    #b is A^(-m) * W^(m^2 /2) , positive m 
    b = a ** (-(np.arange(0,m))) * h_vals[m -1 : 2 * m-1]
    #apply the b in each tile * field 
    x_t = x * np.tile(b[:, np.newaxis], (1, n))
    ###fft ( field * A^(-m) * W^(m^2 /2)  )
    b_fft = sfft.fft(np.array(x_t), N, axis=0, norm="ortho")
    ###convolution: ifft( fft[ field * A^(-m) * W^(m^2 /2)] * fft (W^(-m^2/2)) ) * W^(m^2/2)
    b_ifft = sfft.ifft(b_fft * np.tile(ft[:, np.newaxis], (1, n)), axis=0, norm="ortho")
    bconv = b_ifft[m - 1:mp, :].T * np.tile(h_vals[m - 1:mp], (n, 1)) #positive m until mp, not 2m-1 
    
    #Eq S18, but in center of tiles 
    Mshift = (-m / 2) +  0.5 

    #Eq S19 , phase shift
    Pshift = np.exp(-1j * 2 * np.pi * (np.linspace(0, mout - 1, mout) * (f22 - f11) + mout *f11 ) * (Mshift) / (mout*fs))
    
    #tile over 
    Mshift_til = np.tile(Pshift, (n, 1))
    bout = bconv * Mshift_til 

    
    return bout
    
    
    
def scalar_bluestein(field, screen, d, wavelength):
    """
    Scalar diffraction computation using Bluestein's method.

    Args: 
            :field:         input Field
            :screen:        Observation Screen
            :d:             Observation distance z of the the Screen position 
            :wavelength:    Wavelength
            
    Returns:
            :screen:    Returns the screen populated with the result
    """
    k = 2 * np.pi / wavelength  

    xstart, xend = np.min(screen.x), np.max(screen.x)
    ystart, yend = np.min(screen.y), np.max(screen.y)
    
    x0, y0 = field.XX, field.YY
    x1, y1 = screen.XX[:,:,-1], screen.YY[:,:,-1]
    
    #Eqs 4 &5 
    F0 = np.exp(1j * k * d) / (1j * wavelength * d) * np.exp(1j * k / (2 * d) * (x1**2 + y1**2))
    F = np.exp(1j * k / (2 * d) * (x0**2 + y0**2))

    #arg of Eq 6 
    gout = field.field  * F

   
    #args for bluestein 
    fsx = wavelength * d / field.pixel_x
    fsy = wavelength * d / field.pixel_y
    
    mxout, myout = len(screen.x), len(screen.y)

    # X-direction 
    fx1 = xstart + fsx / 2
    fx2 = xend + fsx / 2
    gout = bluestein_czt(gout, fx1, fx2, fsx, mxout)
      

    #apply calc to obtain the CZT 
    # Y-direction 
    fy1 = ystart + fsy/ 2
    fy2 = yend + fsy / 2
    gout = bluestein_czt(gout, fy1, fy2, fsy, myout)

    #Eq 6, now gout is the result of the CZT
    gout = F0 * gout

    return gout
    
    

def Bluestein(field, screen, wavelength, n=None):
    """
    Implements the Bluestein propagation 
    
    Args: 
        
        :field:     input Field
        :screen:    Observation Screen
        :wavelength:    wavelength to consider
        :n:         refractive index of the propagation medium (default=1 for vacuum/air)
    Returns:
        :screen:    Returns the screen populated with the result
    """
    assert screen.x.shape != () and screen.y.shape != (), "Screen x and y must not be empty"
    xlen,ylen,zlen = screen.XX.shape
    with Timer():
        for z_i in range(zlen):
            z = screen.ZZ[:, :, z_i][-1][0]

            g1  = scalar_bluestein(field, screen, z, wavelength)
            
            screen.screen[:, :, z_i] = g1
            
            progress_bar((z_i)/(zlen))
        progress_bar(1)

    return screen
    
    
def zero_pad(arr, pad_factor=2):
    H, W = arr.shape
    H_new, W_new = pad_factor * H, pad_factor * W
    out = np.zeros((H_new, W_new), dtype=arr.dtype)
    start_H = (H_new - H) // 2
    start_W = (W_new - W) // 2
    out[start_H:start_H + H, start_W:start_W + W] = arr
    return out

def crop_to_physical_size(arr, dx_out, desired_size_meters):
    """
    Center crop an array based on physical size and pixel spacing.
    """
    N_total = arr.shape[0]
    crop_size = int(desired_size_meters / dx_out)
    if crop_size % 2 == 0:
        crop_size -= 1  # force odd for center alignment
        
    start = (N_total - crop_size) // 2
    return arr[start:start+crop_size, start:start+crop_size]

def resize_field_to_shape(field, output_shape):
    """
    Resize a 2D complex array to the desired output shape by
    resizing real and imaginary parts separately.
    """
    from skimage.transform import resize 
    
    amplitude = np.abs(field)
    phase = np.angle(field)

    amp_resized = resize(amplitude, output_shape, mode='constant', anti_aliasing=True)
    phase_resized = resize(phase, output_shape, mode='constant', anti_aliasing=True)

    return amp_resized * np.exp(1j * phase_resized)
    
    
def scalable_angular_spectrum_method(field, screen, z, wavelength, pad_factor=2, skip_final_phase=True, crop=False):
    """
    kernel based on Heintzmann et al. 2023  
    """
    
    N = field.field.shape[0]
    L = np.max(field.x) - np.min(field.x) #assuming square.. 
    
    L_new = pad_factor * L
    N_new = pad_factor * N
    psi_p = zero_pad(field.field, pad_factor=pad_factor)
    
    
    z_limit = (- 4 * L * np.sqrt(8*L**2 / N**2 + wavelength**2) * np.sqrt(L**2 * 1 / (8 * L**2 + N**2 * wavelength**2))\
               / (wavelength * (-1+2 * np.sqrt(2) * np.sqrt(L**2 * 1 / (8 * L**2 + N**2 * wavelength**2)))))
    
    z_limit = (-4 * L_new * np.sqrt(8 * L_new**2 / N_new**2 + wavelength**2)
           * np.sqrt(L_new**2 / (8 * L_new**2 + N_new**2 * wavelength**2)) /
           (wavelength * (-1 + 2 * np.sqrt(2) * np.sqrt(L_new**2 / (8 * L_new**2 + N_new**2 * wavelength**2)))))
    
    if z > z_limit:
        print("Zlimit is " +str(z_limit)+" but z is "+str(z) )
        

    k = 2 * np.pi / wavelength
    df = 1 / L_new
    Lf = N_new * df

    fx = np.fft.fftfreq(N_new, d=1/Lf)
    fy = np.fft.fftfreq(N_new, d=1/Lf)
    FX, FY = np.meshgrid(fx, fy, indexing='ij')

    x = np.fft.ifftshift(np.linspace(-L_new/2, L_new/2, N_new, endpoint=False))
    y = x.copy()
    X, Y = np.meshgrid(x, y, indexing='ij')

    cx = wavelength * FX
    cy = wavelength * FY
    tx = L_new / (2 * z) + np.abs(cx)
    ty = L_new / (2 * z) + np.abs(cy)

    W = ((cx**2 * (1 + tx**2) / tx**2 + cy**2 <= 1) *
         (cy**2 * (1 + ty**2) / ty**2 + cx**2 <= 1))

    H_AS = np.sqrt(np.clip(1 - (cx)**2 - (cy)**2, 0, None))
    H_Fr = 1 - (cx)**2 / 2 - (cy)**2 / 2
    delta_H = W * np.exp(1j * k * z * (H_AS - H_Fr))

    psi_fft = np.fft.fft2(np.fft.ifftshift(psi_p))
    psi_precomp = np.fft.ifft2(psi_fft * delta_H)

    dq = wavelength * z / L_new
    Q = dq * N * pad_factor

    q = np.fft.ifftshift(np.linspace(-Q/2, Q/2, N_new, endpoint=False))
    QX, QY = np.meshgrid(q, q, indexing='ij')

    H_1 = np.exp(1j * k / (2 * z) * (X**2 + Y**2))

    if skip_final_phase:
        psi_p_final = np.fft.fftshift(np.fft.fft2(H_1 * psi_precomp))
    else:
        H_2 = np.exp(1j * k * z) * np.exp(1j * k / (2 * z) * (QX**2 + QY**2))
        psi_p_final = np.fft.fftshift(H_2 * np.fft.fft2(H_1 * psi_precomp))

    #psi_final = zero_unpad(psi_p_final, field.field.shape, pad_factor=pad_factor)
    #psi_final = psi_p_final
    # After propagation:
    dx_out = wavelength * z / L_new  # new pixel size
    desired_output_size = np.max(screen.x) - np.min(screen.x) # e.g. 5 mm screen

    if crop ==True: 
        psi_final = crop_to_physical_size(psi_p_final, dx_out, desired_output_size)
        #print(psi_final.shape)
        psi_final = resize_field_to_shape(psi_final, (N,N) )
    else:
        psi_final = psi_p_final

    return psi_final




def SASM(field, screen, wavelength, pad_factor = 2, crop = False):
    """
    Implements the Scalable Angular Spectrum Method propagation (Heintzmann et al. )
    
    Args: 
        :field:         input Field
        :screen:        Observation Screen
        :wavelength:    wavelength to consider
        :pad_factor:    zero padding factor all around (integrer multiple of the side pixel size) 
        :crop:          boolean, if True, crops the size to the non zero pixels 
    Returns:
        :screen:    Returns the screen populated with the result
    """
    
    xlen,ylen,zlen = screen.XX.shape
    
    g0  = scalable_angular_spectrum_method(field, screen, screen.z[0], wavelength, pad_factor= pad_factor, crop=crop)
    
    #print(np.shape(g0))
    
    x1 = np.linspace(np.min(screen.x), np.max(screen.x), np.shape(g0)[0])
    y1 = np.linspace(np.min(screen.y), np.max(screen.y), np.shape(g0)[1])
    z1 = screen.z
    
    #print(z1)

    newscreen = Screen(x1,y1,z1)

    with Timer():
        for z_i in range(zlen):
            z = screen.ZZ[:, :, z_i][-1][0]

            g1  = scalable_angular_spectrum_method(field, screen, z, wavelength, pad_factor= pad_factor, crop = crop)
            
            #print(np.shape(g1))
            
            newscreen.screen[:, :, z_i] = g1
            
            progress_bar((z_i)/(zlen))
        progress_bar(1)

    return newscreen
    

def ASM_kernel(field, z, wavelength, input_extent, input_df, n = 1.0,  bandlimit = False, shift_yx = (0.0, 0.0),  kykx = (0.0, 0.0)):
    """
    Compute the ASM propagation kernel H(kx, ky)

    Args: 
        :field:        Input Field
        :z:            Distance to Screen plane 
        :wavelength:   Wavelength
        :input_extent: Spatial extent of the input field (might be already padded) 
        :n:            Refractive index of the propagation medium (default=1 for vacuum/air)
        :bandlimit:    Boolean, default False, if True enforces band limit filters akin Matushima & Shimobaba
        :shift_yx:     tuple (shift_y, shift_x) with shift at the output screen plane, default (0.,0.)  
        :kykx:         tuple (ky,kx) angular input direction, default = (0.,0.) 
    Returns:
        :H field kernel:        Returns the kernel 
        
    """
    axes = (-2, -1) 
    
    Ny, Nx = field.Nx, field.Ny
    dy, dx = field.pixel_y, field.pixel_x
    fx, fy = field.fx, field.fy 
    
    fx_grid, fy_grid = field.FX, field.FY
    f_grid = np.stack((fx_grid, fy_grid), axis=-1)
    
    kykx_arr = np.asarray(kykx) / (2 * np.pi) 
    f_shifted = f_grid - kykx_arr
    f2_shifted = np.sum(f_shifted**2, axis=-1)
    
    phase_delay = np.sqrt(np.complex128(1 - (wavelength / n)**2 * f2_shifted))
        
    # shift in output plane
    shift_yx = np.asarray(shift_yx)
    out_shift = 2 * np.pi * np.sum(f_grid * shift_yx, axis=-1)
    
    phase = (2 * np.pi * (n / wavelength) * np.abs(z) * phase_delay) + out_shift
    
    propagator_field = np.where(z >= 0, np.exp(1j * phase), np.conj(np.exp(1j * phase)))
    
    if bandlimit:
        # Bandlimit (Matsushima & Shimobaba)
        z_arr = np.array(z)
        shift_yx_grid = shift_yx[np.newaxis, np.newaxis, :]
        input_extent_grid = input_extent[np.newaxis, np.newaxis, :]
        input_df_grid = input_df[np.newaxis, np.newaxis, :]

        k_limit_p = (
            ((shift_yx_grid + 1 / (2 * input_df_grid) ) ** -2 * z_arr**2 + 1) ** (-1 / 2)
        ) / wavelength * n
        
        k_limit_n = (
            ((shift_yx_grid - 1 / (2 * input_df_grid) ) ** -2 * z_arr**2 + 1) ** (-1 / 2)
        ) / wavelength * n

        # k0: Center of the bandlimit filter 
        k0 = (1 / 2) * (
            np.sign(shift_yx_grid + input_extent_grid) * k_limit_p
            + np.sign(shift_yx_grid - input_extent_grid) * k_limit_n
        )
        
        # k_width: Width of the bandlimit filter 
        k_width = (
            np.sign(shift_yx_grid + input_extent_grid) * k_limit_p
            - np.sign(shift_yx_grid - input_extent_grid) * k_limit_n
        )
        
        k_max = k_width / 2 # Half the width
        
        # H band limit filter 
        H_filter_yx = np.abs(f_grid - k0) <= k_max
        H_filter = H_filter_yx[..., 0] * H_filter_yx[..., 1]
        
        propagator_field = propagator_field * H_filter
        
    return np.fft.ifftshift(propagator_field, axes=axes)




def ASM_propagate(field, screen, z, wavelength, pad_width, n = 1.0, mode = None, bl = True, shift = None, kykx = (0.0, 0.0)):
    """
    Angular Spectrum Method (ASM) propagation computation.

    Args: 
        :field:       Input Field
        :screen:      Observation Screen
        :z:           Distance to Screen plane 
        :wavelength:  Wavelength
        :pad_width:   Padding around area of interest (assumed constant =0 all around), in nr. of pixels  
        :n:           refractive index of the propagation medium (default=1 for vacuum/air)
        :mode:        ASM mode, options: None (default) = convetional ASM; "czt" = Chirp Z-Transform; "BLAS" = Band-Limited ASM
        :bl:          Boolean, default False, if True enforces band limit filters akin Matushima & Shimobaba
        :shift:       tuple (shift_y, shift_x) with shift at the output screen plane, default None = calculates the shift from screen limits. If mode="czt" shift is not used. 
        :kykx:        tuple (ky,kx) angular input direction, default = (0.,0.), required for "off-axis" for non-czt 
    Returns:
        :field:       Returns the calculated field
    """
    # zero pad the field  --- asssumes symmetric padding!
    padding = ((pad_width, pad_width), (pad_width, pad_width))
    padded_vals = np.pad(field.field, padding, mode='constant', constant_values=0)

    #print(padded_vals.shape)

    paddedx = np.linspace(field.x.min() - pad_width*field.pixel_x, field.x.max() + pad_width*field.pixel_x, len(field.x) + 2*pad_width )
    paddedy = np.linspace(field.y.min() - pad_width*field.pixel_y, field.y.max() + pad_width*field.pixel_y, len(field.y) + 2*pad_width )

    #print(field.x, paddedx, field.x.shape, paddedx.shape)
    #print(field.y, paddedy, field.y.shape, paddedy.shape)
    
    padded_field = Field(paddedx, paddedy)
    
    padded_field.field = padded_vals
    
    axes = (-2, -1)  #define axes as last two dims for generalization, although with 2D does not make a difference
    spatial_shape = np.array(padded_field.shape)
    
    input_dx = np.array([field.pixel_x, field.pixel_y])
    mask_field_extent = input_dx * spatial_shape 
    input_df = 1.0 / mask_field_extent 
    

    # Handle shifts from screen coordinates or given by user 
    # NOTE, czt mode does not take shifts
    if (shift is None) & (mode != "czt"):
        ymin, ymax = np.min(screen.y), np.max(screen.y)
        xmin, xmax = np.min(screen.x), np.max(screen.x)

        shift_y = (ymax + ymin)/2
        shift_x = (xmax + xmin)/2
        
        shift_yx_for_kernel = (shift_x, shift_y) 
    
    elif (shift is None) & (mode == "czt"):
        shift_yx_for_kernel = (0.,0.) 
    elif (shift is not None):
        shift_yx_for_kernel = shift
    else: 
        shift_yx_for_kernel = (0.,0.)

    # Compute kernel H(kx, ky)
    kernel_H = ASM_kernel(padded_field, z, wavelength, input_extent=mask_field_extent, \
                          input_df=input_df, shift_yx = shift_yx_for_kernel, bandlimit = bl, n = n, kykx=kykx )
    
    # APPLY KERNEL FORWARD TRANSFORM 
    field_transform = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(padded_field.field, axes=axes), axes=axes) * kernel_H, axes=axes)

    # INVERSE TRANSFORM
    if mode=="czt":
        output_shape = np.array(screen.shape)
        ymin, ymax = screen.y.min(), screen.y.max()
        xmin, xmax = screen.x.min(), screen.x.max()
        
        output_dx = (ymax - ymin) / (output_shape[0] - 1), (xmax - xmin) / (output_shape[1] - 1)
        
        # Scaling factor: alpha = output_dx / input_df
        alpha = output_dx / input_df 

        limits_min = [xmin, ymin]
        limits_max = [xmax, ymax]
        
        
        for d, axis in enumerate(axes):
            m = output_shape[d]
            
            # czt parameters, 
            # a: starting point on the circle (related to the min limit) w: angular step/ratio (related to the span/range)
            a_czt = np.exp(-1j * 2 * np.pi / mask_field_extent[d] * limits_min[d])
            w_czt = np.exp(1j * (2 * np.pi / mask_field_extent[d]) * (limits_max[d] - limits_min[d]) / (m - 1))

            # apply the czt
            field_czt = scipy.signal.czt(field_transform, m, w_czt, a_czt, axis=axis)
            
            # Ensure complex type just to be sure
            field_transform = field_czt.astype(field_transform.dtype)
            
            center = (m - 1) // 2 
            
            # phase compensation/modulation factor after the czt
            compensation = w_czt ** (-center * np.arange(m)) * (a_czt**center)
            
            field_transform = np.moveaxis(field_transform, axis, -1)
            field_transform = field_transform * compensation 
            field_transform = np.moveaxis(field_transform, -1, axis)
        
        # Apply the final scaling factor
        final_scaling = np.prod(1.0 / alpha)
        field_transform = field_transform * final_scaling
        
    elif (mode == "BLAS") & (bl==True):
        #print("BLAS")
        fx_grid, fy_grid = padded_field.FX, padded_field.FY
        f_grid = np.stack((fx_grid, fy_grid), axis=-1)
            
        output_shape = np.array(screen.shape)
        ymin, ymax = screen.y.min(), screen.y.max()
        xmin, xmax = screen.x.min(), screen.x.max()
        
        output_dx_y = (ymax - ymin) / (output_shape[0] - 1)
        output_dx_x = (xmax - xmin) / (output_shape[1] - 1)

        output_dx = np.array([output_dx_x, output_dx_y])
        
        # Scaling factor: alpha = output_dx / input_df (Eq 7)
        alpha = output_dx / input_df
        
        # Eq 9 of "Band-limited angular spectrum numerical propagation method
        # with selective scaling of observation window size and sample number"
        # (2012)
        wn = alpha * f_grid 
        
        # f = kernel for convolution, Eq 9 first term
        f = np.prod(np.exp(-1j * np.pi / alpha * wn**2), axis=-1)
        # B = modulated k-space field, Eq 9 second term
        B = field_transform * np.prod( (1 / alpha) * np.exp(1j * np.pi / alpha * wn**2), axis=-1)

        prefactor = np.prod( output_dx * np.exp(1j * np.pi / alpha * f_grid**2),axis=-1,)
        
        field_transform = prefactor * scipy.signal.fftconvolve(B, f, mode="same", axes=axes)
    
        # crop to unpadded 
        y_slice = slice(int(padded_field.Ny/2) - int(screen.Ny/2),  int(padded_field.Ny/2) - int(screen.Ny/2) + screen.Ny)
        x_slice = slice(int(padded_field.Nx/2) - int(screen.Nx/2),  int(padded_field.Nx/2) - int(screen.Nx/2) + screen.Nx)
        field_transform = field_transform[x_slice, y_slice]

    else:   
        #print("Conventional ")
        # just IFFT 
        propagated_field = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(field_transform, axes=axes), axes=axes), axes=axes)

        #fig = plt.figure() 
        #plt.imshow(np.abs(propagated_field))

        scaling_x = screen.pixel_x / padded_field.pixel_x
        scaling_y = screen.pixel_y / padded_field.pixel_y

        y_slice = slice(int(padded_field.Ny/2) - int(screen.Ny/2 *scaling_y),  \
                        int(padded_field.Ny/2) - int(screen.Ny/2 *scaling_y) + int(screen.Ny*scaling_y) )
        x_slice = slice(int(padded_field.Nx/2) - int(screen.Ny/2 *scaling_y), \
                        int(padded_field.Nx/2) - int(screen.Ny/2 *scaling_y) + int(screen.Nx*scaling_x) )
        field_transform = propagated_field[x_slice, y_slice]

        #fig = plt.figure() 
        #plt.imshow(np.abs(psi_final))
        
        if (screen.Nx != field.Nx) or (screen.Ny != field.Ny):
            field_transform = resize_field_to_shape(field_transform, (screen.Nx,screen.Ny) )
        else: 
            field_transform = resize_field_to_shape(field_transform, (field.Nx,field.Ny) )
        
    #print(field_transform.shape)
        
    return field_transform.astype(complex)


def ASM(field, screen, wavelength, pad, n = 1.0, mode = None, bl = False, shift = None, kykx = (0.,0.) ):
    """
    Implements the Angular Spectrum Method (ASM), with czt and bandlimit options.

    Args: 
        :field:       Input Field
        :screen:      Observation Screen
        :wavelength:  Wavelength
        :pad:         Padding around area of interest (assumed constant =0 all around), in nr. of pixels  
        :n:           refractive index of the propagation medium (default=1 for vacuum/air)
        :mode:        ASM mode, options: None (default) = convetional ASM; "czt" = Chirp Z-Transform; "BLAS" = Band-Limited ASM
        :bl:          Boolean, default False, if True enforces band limit filters akin Matushima & Shimobaba
        :shift:       tuple (shift_y, shift_x) with shift at the output screen plane, default None = calculates the shift from screen limits. If mode="czt" shift is not used. 
        :kykx:        tuple (ky,kx) angular input direction, default = (0.,0.) 
    Returns:
        :screen:      Returns the screen populated with the result   
    """
    xlen, ylen, zlen = screen.Nx, screen.Ny, screen.Nz

    info = "with band limit." if bl else "without band limit."
    
    if mode == "czt": 
        print("ASM-CZT mode ON, "+info)
    elif mode == "BLAS": 
        print("BLAS mode ON.")  
    else:
        print("Conventional ASM, "+info)

    
    with Timer():
        for z_i, z in enumerate(screen.z):
    
            g1 = ASM_propagate(field, screen, z, wavelength, pad, n=n, mode=mode, bl=bl, shift=shift, kykx=kykx)
            
            if screen.Ny == 1:
                screen.screen[:, :, z_i] = np.reshape(g1, (screen.Nx, screen.Ny))
            else:
                screen.screen[:, :, z_i] = g1
            
            progress_bar((z_i) / (zlen))
        progress_bar(1)

    return screen

    
def propagate_through_ensemble(ensemble,  wavelength , xar_plus_z=None, propagation_methods_array=None): 
    """
    Propagates though Ensemble object (various MOE surfaces)  
    
    Args:
        :ensemble:                   Ensemble object 
        :wavelength:                 Propagation wavelength
        :xar_plus_z:                 Array with optimization variables (and optionally z position as variable) 
        :propagation_methods_array:  Propagation method 
    Returns:
        :overall_field_at_screen:    Returns a screen concatenating both 
    """ 
    
    #initial aperture (the rest will be inside a loop)
    
    fieldi = ensemble.input_light_field 
    ensemble.field_array[0] = fieldi
    
    field0 = modulate_field(fieldi, amplitude_mask = ensemble.aperture_array_amp[0],\
                                      phase_mask=ensemble.aperture_array_phase[0])

    print("Surface #1")

    plot_field(field0)
    
    screen0 = ensemble.screen_array[0]
    #print(np.max(screen0.z))
    
    if propagation_methods_array==None:
        propagation_methods_array= ["bluestein" for i in ensemble.screen_array]
    
    
    if propagation_methods_array[0]=="bluestein":
        EXYZ0 = Bluestein(field0, screen0, wavelength)
        ensemble.field_array[1].field = EXYZ0.screen[:,:, -1]/np.max(EXYZ0.screen[:,:, 0])
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)

    elif propagation_methods_array[0]=="ASM": 
        EXYZ0 = ASM(field0, screen0, wavelength, pad=int(len(field0.x)/2 ))
        ensemble.field_array[1].field = EXYZ0.screen[:,:, -1]/np.max(EXYZ0.screen[:,:, 0])
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)
        
        
    elif propagation_methods_array[0]=="ASM-CZT": 
        EXYZ0 = ASM(field0, screen0, wavelength, pad=int(len(field0.x)/2 ), mode="czt")
        ensemble.field_array[1].field = EXYZ0.screen[:,:, -1]/np.max(EXYZ0.screen[:,:, 0])
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)

    elif propagation_methods_array[0]=="SASM": 
        EXYZ0 = SASM(field0, screen0, wavelength, pad_factor=4, crop=True)
        ensemble.field_array[1].field = EXYZ0.screen[:,:, -1]/np.max(EXYZ0.screen[:,:, 0])
        
        EXYZ = EXYZ0
        
        overall_field_at_screen = screen0
        overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)
        
    elif propagation_methods_array[0]=="rayleigh-sommerfeld":
        screen_XY = create_screen_XY(np.min(screen0.x), np.max(screen0.x), len(screen0.x), 
                                        np.min(screen0.y), np.max(screen0.y), len(screen0.y), 
                                        z=np.max(screen0.z) )

        EXYZ0 = RS_integral(field0, screen_XY, wavelength, simp2d=True)
        ensemble.field_array[1].field = EXYZ0.screen[:,:,0]/np.max(EXYZ0.screen[:,:, 0])
    
        overall_field_at_screen = screen0
        overall_field_at_screen.screen[:,:,-1] = EXYZ0.screen[:,:,0]/np.max(EXYZ0.screen[:,:, 0])
        
        print("RS only calculated correctly the XY at the end Z, the other values of the screen object (e.g. YZ plane) have been repeated")
    
    EXYZ = EXYZ0
    
    #to the rest for all 
    for i in np.arange(1,len(ensemble.aperture_array_phase)): 
        #print(i)
        field = ensemble.field_array[i]
        
        aperture_amp = ensemble.aperture_array_amp[i]
        aperture_phase = ensemble.aperture_array_phase[i]
    
        field = modulate_field(field, amplitude_mask = aperture_amp, phase_mask=aperture_phase)

        print("Surface #"+str(i+1))
        plot_field(field)

    
        ensemble.field_array[i].field = field
        
        screen = ensemble.screen_array[i]
        
        if propagation_methods_array[i]=="bluestein":
            EXYZ1 = Bluestein(field, screen, wavelength)
            EXYZ = EXYZ1
            
            if i<len(ensemble.aperture_array_phase)-1 : 
                ensemble.field_array[i+1].field = EXYZ.screen[:,:, -1]/np.max(EXYZ.screen[:,:, -1])

            overall_field_at_screen = screen
            overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)
            
        if propagation_methods_array[i]=="ASM":
            EXYZ1 = ASM(field, screen, wavelength, pad=int(len(field0.x)/2) )
                        
            EXYZ = EXYZ1
            
            if i<len(ensemble.aperture_array_phase)-1 : 
                ensemble.field_array[i+1].field = EXYZ.screen[:,:, -1]/np.max(EXYZ.screen[:,:, -1])


            overall_field_at_screen = screen
            overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)


        if propagation_methods_array[i]=="SASM":
            EXYZ1 = SASM(field, screen, wavelength, pad_factor=4, crop=True)
            
            EXYZ = EXYZ1
            
            if i<len(ensemble.aperture_array_phase)-1 : 
                ensemble.field_array[i+1].field = EXYZ.screen[:,:, -1]/np.max(EXYZ.screen[:,:, -1])
            
            overall_field_at_screen = screen
            overall_field_at_screen.screen = EXYZ.screen/np.max(EXYZ.screen)
        
        
        # NOT TESTED - TO REVIEW RS 
        elif propagation_methods_array[i]=="rayleigh-sommerfeld":
            if i==len(ensemble.aperture_array_phase)-1:
                print("YZ")
                screen_YZ = create_screen_YZ(np.min(screen.y), np.max(screen.y), len(screen.y), 
                                            np.min(screen.z), np.max(screen.z), len(screen.z), 
                                            #z=np.max(screen0.z) )
                                            x=0 )

                plot_field(field)
        
                EXYZ0 = RS_integral(field, screen_YZ, wavelength, simp2d=True)
                
                if i<len(ensemble.aperture_array_phase)-1 :
                    ensemble.field_array[i+1].field = EXYZ0.screen[:,:,-1]/np.max(EXYZ0.screen[:,:, -1])

                overall_field_at_screen = screen
                overall_field_at_screen.screen[int(len(screen0.x)/2),:,:] = EXYZ0.screen[:,0,:]
 
                EXYZ = screen

                
            else:  
                print("XY")
                screen_XY = create_screen_XY(np.min(screen.x), np.max(screen0.x), len(screen0.x), 
                                                np.min(screen0.y), np.max(screen0.y), len(screen0.y), 
                                                z=np.max(screen0.z) )
                                                
                EXYZ0 = RS_integral(field, screen_XY, wavelength, simp2d=True)
                
                print(np.shape(EXYZ0.screen[:,:,0]), np.shape(EXYZ0.screen))
                
                if i<len(ensemble.aperture_array_phase)-1 :
                    ensemble.field_array[i+1].field = EXYZ0.screen[:,:,0]/np.max(EXYZ0.screen[:,:, -1])

                overall_field_at_screen = screen
                overall_field_at_screen.screen[:,:,-1] = EXYZ0.screen[:,:,0]

                #attributing this as screen to avoid dimension errors, however ONLY the XY is meaningful!
                #might be fixed 
                EXYZ = screen

                print("RS only calculated correctly the XY at the end Z, the other values of the screen object (e.g. YZ plane) have been repeated")
        
        del field 
        
        #print(EXYZ.screen.shape)        
        #overall_field_at_screen = screen0
        #overall_field_at_screen.screen = EXYZ.screen
        
        if i == 1:
            overall_field_at_screen2 = np.concatenate( (EXYZ0.screen, EXYZ.screen), axis=2)
            #print(overall_field_at_screen2.shape)
        else:
            #print(EXYZ.screen.shape, overall_field_at_screen.shape)
            overall_field_at_screen2 = np.concatenate((overall_field_at_screen2, EXYZ.screen), axis=2)
            #print(overall_field_at_screen2.shape)
            
        
    return overall_field_at_screen2 
    
    
##################################################
##RS integral functions introduced in v1.3
    
    
def RS_intXY(zs, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength, verbose =False ): 
    """
    Calculates the RS int of the first kind
    returns Escreen (complex electric field at obs screen), Iscreen (intensity at obs screen), iplot (the actual intensity) 
    
    Args: 
        :zs:            distance to screen [m]
        :mask:          Electric field at the mask, complex valued 2D function   
        :npixmask:      number of pixels on the side of the mask 
        :pixsizemask:   size of pixel of the mask [m]
        :npixscreen:    number of pixels on the side of the screen (observation) 
        :dxscreen:      max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
        :dyscreen:      max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
        :wavelength:    wavelength of the light [m]
        :verbose:       defaults to False, if True prints 
    """
    # set the precision to double that of float64
    decimal.setcontext(decimal.Context(prec=34))

    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    #if nps <= 2*npm: #Increase sampling values 
    #    nps = npm*4 +1
    
    #size of mask 
    #dmask = pixsizemask * npm
    dmask = 2*dxscreen
    
    #prop const
    k = 2* np.pi/wavelength
    
    #define the zpos of the mask at 0 
    zm =0 
    
    #definition of the mask
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1,ym1)
    
    xs1 = np.linspace(-dxscreen, dxscreen, nps)
    ys1 = np.linspace(-dyscreen, dyscreen, nps)
    (xs, ys) = np.meshgrid(xs1,ys1)
    
    #definition of the electric field at the mask 
    E0m = mask

    ###### calculate the Rayleigh Sommerfeld integral 
    Escreen = RS_intXY_kernel(dmask, nps, npm, xs, ys,zs,xm,ym,zm,k, E0m,verbose)

    Iscreen = np.abs(Escreen)**2
    iplot = 10*Iscreen**0.2
    
    return Escreen, Iscreen, iplot 


def RS_intXY_kernel(dmask, nps, npm, xs, ys,zs,xm,ym,zm,k, E0m,verbose):
    ## definitions 
    unit = np.ones((npm,npm), dtype=complex)
    r = np.zeros((npm,npm)) 
    prop1 = np.zeros((npm,npm), dtype=complex)
    prop2 = np.zeros((npm,npm), dtype=complex)
    propE = np.zeros((npm,npm), dtype=complex)
    
    #electric field real and imaginary and total 
    rEs =np.zeros ((nps,nps))
    iEs =np.zeros ((nps,nps))
    Escreen =np.zeros ((nps,nps), dtype =complex)

    with Timer():
        #e.g. Mahajan 2011 part II eq 1-20
        for isc in np.arange(0,nps-1):
            if verbose == True: 
                progress_bar(isc/nps)
                
            for jsc in np.arange(0,nps-1): 
                r = np.sqrt((xs[isc,jsc]-xm)**2 + (ys[isc,jsc]-ym)**2 + (zs-zm)**2)
                r2 = r*r
                prop1= np.exp(r*1.0j*k)/r2
                prop2 = zs * k/(2*np.pi) *(unit / (r*k)  - 1.0j * unit)
                propE = E0m * prop1 * prop2

                #rEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.real(propE))/(2*np.pi)
                #iEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.imag(propE))/(2*np.pi)
                
                Escreen[isc,jsc] = simpson2d(propE,-dmask/2, dmask/2, -dmask/2, dmask/2) /(2*np.pi)
        progress_bar(1)
        
    return Escreen 
    
    
def RS_intZZ(zmin, zmax, nzs, xfixed, yfixed, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength, nind, verbose =False ): 
    """
    Calculates the RS_int in the  of the first kind, taking information about mask, distance to screen, and screen information
    returns Escreen (complex electric field at obs screen), Iscreen (intensity at obs screen), iplot (the actual intensity) 
        
    Args: 
        :[zmin, zmax]:      distance range limits in z [m]
        :nzs:               number of points along the optical axis 
        :xfixed, yfixed:    x and y fixed coordinates
        :mask:              Electric field at the mask, complex valued 2D function   
        :npixmask:          number of pixels on the side of the mask 
        :pixsizemask:       size of pixel of the mask [m]
        :npixscreen:        number of pixels on the side of the screen (observation) 
        :dxscreen:          max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
        :dyscreen:          max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
        :nind:              scaling refractive index 
        :wavelength:        wavelength of the light [m]
        :verbose:           defaults to False, if True prints 
    """
    # set the precision to double that of float64
    decimal.setcontext(decimal.Context(prec=34))

    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    #size of mask 
    #dmask = pixsizemask * npm
    dmask = 2*dxscreen
    
    
    #prop const
    k = 2* np.pi*nind/(wavelength)

    #define the zpos of the mask at 0 
    zm =0 

    #definition of the mask
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1,ym1)

    #definition of the electric field at the mask 
    E0m = mask

    ##zdists 
    zarray = np.linspace(zmin,zmax,nzs)
    
    ####calculation of the integral 
    Escreen = RS_intZZ_kernel(dmask, npm, nzs, xfixed, yfixed,zarray,xm,ym,zm,k, E0m, verbose)

    Iscreen = np.abs(Escreen)**2
    iplot = 10*Iscreen**0.2

    return Escreen, Iscreen, iplot 

 
def RS_intZZ_kernel(dmask, npm, nzs, xs, ys,zs,xm,ym,zm,k, E0m, verbose): 
    ## definitions 
    unit = np.ones((npm,npm), dtype=complex)
    r = np.zeros((npm,npm)) 
    prop1 = np.zeros((npm,npm), dtype=complex)
    prop2 = np.zeros((npm,npm), dtype=complex)
    propE = np.zeros((npm,npm), dtype=complex)

    #electric field real and imaginary and total 
    rEs =np.zeros(nzs)
    iEs =np.zeros(nzs)
    Escreen =np.zeros(nzs, dtype=complex) 

    
    ###### calculate the Rayleigh Sommerfeld integral 
    #e.g. Mahajan 2011 part II eq 1-20           
    with Timer():
        for jsc in np.arange(0,nzs-1): 
            if verbose == True: 
                progress_bar(jsc/nzs)
                
            r = np.sqrt((xs-xm)**2 + (ys-ym)**2 + (zs[jsc]-zm)**2)
            r2 = r*r
            prop1= np.exp(r*1.0j*k)/r2
            prop2 = zs[jsc] * k/(2*np.pi) *(unit / (r*k)  - 1.0j * unit)
            propE = E0m * prop1 * prop2

            #rEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.real(propE))/(2*np.pi)
            #iEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.imag(propE))/(2*np.pi)

            Escreen[jsc] = simpson2d(propE,-dmask/2, dmask/2, -dmask/2, dmask/2) /(2*np.pi)
        progress_bar(1)   
          
    return Escreen
 
 
def RS_intYZ(zmin, zmax, nzs, yfixed, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength, nind, verbose =False ): 
    """
    Calculates the RS_int in the  of the first kind, taking information about mask, distance to screen, and screen information
    returns Escreen (complex electric field at obs screen), Iscreen (intensity at obs screen), iplot (the actual intensity) 
    
    Args: 
        
        :[zmin, zmax]:      distance range limits in z [m]
        :nzs:               number of points along the optical axis 
        :yfixed:            y fixed coordinates
        :mask:              Electric field at the mask, complex valued 2D function   
        :npixmask:          number of pixels on the side of the mask 
        :pixsizemask:       size of pixel of the mask [m]
        :npixscreen:        number of pixels on the side of the screen (observation) 
        :dxscreen:          max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
        :dyscreen:          max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
        :nind:              scaling refractive index 
        :wavelength:        wavelength of the light [m]
        :verbose:           defaults to False, if True prints 

    """
    # set the precision to double that of float64
    decimal.setcontext(decimal.Context(prec=34))

    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    #size of mask 
    #dmask = pixsizemask * npm
    dmask = 2*dxscreen
    
    #prop const
    k = 2* np.pi*nind/(wavelength)

    #define the zpos of the mask at 0 
    zm =0 

    #definition of the mask
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1,ym1)
    
    #definition of the electric field at the mask 
    E0m = mask 

    xfixed = 0
    ys = np.linspace(-dmask/2, dmask/2, npixscreen)
    
    ##zdists 
    zarray = np.linspace(zmin,zmax,nzs)

    ####calculation of the integral 
    Escreen = RS_intYZ_kernel(dmask, npm, nzs, npixscreen, xfixed, ys,zarray,xm,ym,zm,k, E0m, verbose)

    #Escreen = rEs + 1.0j*iEs 
    Iscreen = np.abs(Escreen)**2
    iplot = 10*Iscreen**0.2
    
    return Escreen, Iscreen, iplot 


def RS_intYZ_kernel(dmask, npm, nzs, npixscreen, xfixed, ys,zarray,xm,ym,zm,k, E0m, verbose): 
    ## definitions 
    unit = np.ones((npm,npm), dtype=complex)
    r = np.zeros((npm,npm)) 
    prop1 = np.zeros((npm,npm), dtype=complex)
    prop2 = np.zeros((npm,npm), dtype=complex)
    propE = np.zeros((npm,npm), dtype=complex)

    #electric field real and imaginary and total 
    rEs =np.zeros ((npixscreen,nzs))
    iEs =np.zeros ((npixscreen,nzs))
    Escreen =np.zeros((npixscreen,nzs), dtype=complex) 

    ###### calculate the Rayleigh Sommerfeld integral 
    #e.g. Mahajan 2011 part II eq 1-20           
    with Timer():
        for isc in np.arange(0,npixscreen):
            if verbose == True: 
                progress_bar(isc/npixscreen)
                
            for jsc in np.arange(0,nzs-1): 
                r = np.sqrt((xfixed-xm)**2 + (ys[isc]-ym)**2 + (zarray[jsc]-zm)**2)
                r2 = r*r
                prop1= np.exp(r*1.0j*k)/r2
                prop2 = zarray[jsc] * k/(2*np.pi) *(unit / (r*k)  - 1.0j * unit)
                propE = E0m * prop1 * prop2
                
                
                #rEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.real(propE))/(2*np.pi)
                #iEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.imag(propE))/(2*np.pi)
                
                Escreen[isc,jsc] = simpson2d(propE,-dmask/2, dmask/2, -dmask/2, dmask/2) /(2*np.pi) 
        progress_bar(1)
        
    return Escreen