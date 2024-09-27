"""
propagate.py   

This module had origin in the "propopt" repository module with same name 
-> For more information, examples of use and validation tests, 
   please see link for repository : https://github.com/cunhaJ/propopt 

Here, "propopt" code has been modified and extended for use with apertures and masks of "pyMOE"
"""

import numpy as np
import scipy.fftpack as sfft 

import decimal

from pyMOE.utils import simpson2d 

from scipy import integrate

import dask
import dask.array as da
from dask import delayed
import dask.bag as db
from dask.diagnostics import ProgressBar

from pyMOE.utils import progress_bar, Timer



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
def kernel_RS(field, k, x,y,z, simp2d=False):
    """
    Calculates the RS kernel integral from a field input aperture, assumed to be at z=0
    and returns the calculated E field
    
    Implements the Kernel in Mahajan 2011 part II eq 1-20 

    Args:
        :field:     input field
        :k:         Calculated wavenumber k=2pi/(wl*n)
        :x,y,z:     x, y, z coordinates of the screen point being evaluated
        :simp2d:    Defaults False, if True uses the simpson2d function
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
    if simp2d==True: 
        Exyz = simpson2d(propE,field.x[0], field.x[-1], field.y[0], field.y[-1]) /(2*np.pi)
    else: 
        Exyz = integrate.simpson(integrate.simpson(propE, field.x),field.y)/(2*np.pi) 

    return Exyz
    

def RS_integral(field, screen, wavelength, n=None, parallel_computing=True, simp2d=False):
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
                    result = kernel_RS(field, k ,x,y,z, simp2d)

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
                        
                        result = kernel_RS(field, k ,x,y,z, simp2d).compute()

                        screen.screen[x_i, y_i, z_i] = result
                        progress_bar((x_i*zlen*ylen+y_i*zlen+z_i)/(xlen*ylen*zlen))
            progress_bar(1)

    return screen


    
    
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