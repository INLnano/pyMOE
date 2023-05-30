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

##We will use the simpson method for cnumerical calc of integral
from pyMOE.utils import simpson2d 
##Another integration methods are also possible 
#from pyMOE.utils import double_Integral
 
from pyMOE.utils import progress_bar, Timer

def fresnel(z, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength):
    """
    Calculate Fresnel approximation, following Goodman exp 4-17 
    
    inputs: 
        z = distance to the observation plane in m
        mask = 2D map (x-y plane) npixel vs npixel mask 
        npixmask = number of pixels of the mask 
        pixsizemask = size of the pixel at the mask in m
        npixscreen = size of the pixel at the screen 
        dxscreen = x-size of the screen  in m 
        dyscreen = y-size of the screen   in m
        wavelength = wavelength in m 
    
    returns 2D map (x-y plane) with abs of Electric field at distance z
    
    """ 
    k = 2* np.pi/wavelength
    
    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    ##calculate the resolution  -
    res = wavelength *z/ (pixsizemask *npm)
    
    dmask =   res * npm
    
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
    #v2  = np.exp(1.0j*k* (xs*xs + ys*ys)/ (2*z)) 
    #v3  = np.exp(1.0j*k*z)/ (1.0j*wavelength*z)
    intarg = v1 * mask
    Ef = sfft.fftshift(sfft.fft2(sfft.ifftshift(intarg)))

    return Ef 

def fraunhofer(z, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength):
    """
    Calculate Fraunhofer approximation, following Goodman exp 4-25
    
    inputs: 
        z = distance to the observation plane in m
        mask = 2D map (x-y plane) npixel vs npixel mask 
        npixmask = number of pixels of the mask 
        pixsizemask = size of the pixel at the mask in m
        npixscreen = size of the pixel at the screen 
        dxscreen = x-size of the screen  in m 
        dyscreen = y-size of the screen   in m
        wavelength = wavelength in m 

    returns 2D map (x-y plane) with abs of Electric field at distance z 
    
    """
    
    k = 2* np.pi/wavelength
    #number of pixels 
    nps = npixscreen
    npm = npixmask 
    
    dmask = npixmask * npm  
    
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
    

def RS_intXY(zs, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength, verbose =False ): 
    """
    Calculates the RS int of the first kind
    returns Escreen (complex electric field at obs screen), Iscreen (intensity at obs screen), iplot (the actual intensity) 
    
    Inputs: 
        zs = distance to screen [m]
        mask = Electric field at the mask, complex valued 2D function   
        npixmask = number of pixels on the side of the mask 
        pixsizemask = size of pixel of the mask [m]
        npixscreen = number of pixels on the side of the screen (observation) 
        dxscreen = max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
        dyscreen = max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
        wavelength = wavelength of the light [m]

    ------- 
    optional 
    verbose, defaults to False, if True prints 
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
                prop2 = zs * (1.0j * k  - unit/r)
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
        
    Inputs: 
        [zmin, zmax] = distance range limits in z [m]
        nzs = number of points along the optical axis 
        xfixed, yfixed = x and y fixed coordinates
        mask = Electric field at the mask, complex valued 2D function   
        npixmask = number of pixels on the side of the mask 
        pixsizemask = size of pixel of the mask [m]
        npixscreen = number of pixels on the side of the screen (observation) 
        dxscreen = max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
        dyscreen = max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
        nind = scaling refractive index 
        wavelength = wavelength of the light [m]

    ------- 
    optional 
    verbose, defaults to False, if True prints 
    """
    # set the precision to double that of float64
    decimal.setcontext(decimal.Context(prec=34))

    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    #size of mask 
    dmask = pixsizemask * npm
    
    #prop const
    k = 2* np.pi/(wavelength*nind)

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
            prop2 = zs[jsc] * (1.0j * k  - unit/r)
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
    
    Inputs: 
        [zmin, zmax] = distance range limits in z [m]
        nzs = number of points along the optical axis 
        yfixed = y fixed coordinates
        mask = Electric field at the mask, complex valued 2D function   
        npixmask = number of pixels on the side of the mask 
        pixsizemask = size of pixel of the mask [m]
        npixscreen = number of pixels on the side of the screen (observation) 
        dxscreen = max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
        dyscreen = max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
        nind = scaling refractive index 
        wavelength = wavelength of the light [m]
    
    ------- 
    optional 
    verbose, defaults to False, if True prints 
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
    k = 2* np.pi/(wavelength*nind)

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
                prop2 = zarray[jsc] * (1.0j * k - unit/r)
                propE = E0m * prop1 * prop2
                
                #rEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.real(propE))/(2*np.pi)
                #iEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*100,npm*100,np.imag(propE))/(2*np.pi)
                
                Escreen[isc,jsc] = simpson2d(propE,-dmask/2, dmask/2, -dmask/2, dmask/2) /(2*np.pi) 
        progress_bar(1)
        
    return Escreen


def Fresnel_num(width, wavelength, zdist):
    """
    Calculation of Fresnel number, Goodman pag 85
    inputs:
    width = size of the aperture in m 
    wavelength = wavelength of the aperture in m 
    zdist = distance to the screen in m 
    
    returns Fresnel number 
    """
    NF = width**2 / (wavelength * zdist)
    return NF 
    
def Fraunhofer_criterion(aperturesiz, wavelength): 
    """
    Calculation of "Fraunhofer distance" , Goodman exp 4-27  
    inputs:
    aperturesiz = size of the aperture in m 
    wavelength = wavelength of the aperture in m 
    
    returns Fraunhofer distance 
    """
    zfraun = 2 * aperturesiz**2 / wavelength 
    
    return zfraun 