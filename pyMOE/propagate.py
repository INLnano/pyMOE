####propagate.py   

####This file is "duplicate" from the propopt repository on January 2022
####For more information, examples of use and validation tests, please see
####link for repository : https://github.com/cunhaJ/propopt 

import numpy as np
import scipy.fftpack as sfft 

import decimal

import numpy as np  

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
    
    #Goodman, exp 4-17
    v1  = np.exp(1.0j*k* (xm*xm + ym*ym)/ (2*z))
    v2  = np.exp(1.0j*k* (xs*xs + ys*ys)/ (2*z)) 
    v3  = np.exp(1.0j*k*z)/ (1.0j*wavelength*z)
    intarg = v1 * mask
    Ef = sfft.fftshift(sfft.fft2(sfft.ifftshift(intarg)))

    #consider changing this into complex field 
    return abs(Ef)
    

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

    #delta2 = z*wavelength /(npm*delta1)
    #num = np.round(dxscreen/delta2)
    #print(num)
    
    #import cv2 
    #resized = cv2.resize(mask,(int(4*num),int(4*num)), interpolation = cv2.INTER_AREA)
    resized = mask
    
    #Goodman, exp 4-25
    v2  = np.exp(1.0j*k* (xs*xs + ys*ys)/ (2*z)) 
    v3  = np.exp(1.0j*k*z)/ (1.0j*wavelength*z)
    Ef = v2 * v3* sfft.fftshift(sfft.fft2(resized))
    

    #print(resized.shape)
    
    #Ef2 = cv2.resize(Ef,(int(num),int(num)), interpolation = cv2.INTER_AREA)
    
    return abs(Ef)
    

def RS_int(zs, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength, I0, verbose =False ): 
    """
    Calculates the RS_int in the  of the first kind, taking information about mask, distance to screen, and screen information
    returns Escreen (complex electric field at obs screen), Iscreen (intensity at obs screen), iplot (the actual intensity) 
    inputs: 
    zs = distance to screen [m]
    mask = image of the mask (normalized to [0,1]) - in principle could be grayscale between 0 and 1 
    npixmask = number of pixels on the side of the mask 
    pixsizemask = size of pixel of the mask [m]
    npixscreen = number of pixels on the side of the screen 
    dxscreen = max_x of the screen [m], the screen range is [-dxscreen, dxscreen]
    dyscreen = max_y of the screen [m], the screen range is [-dyscreen, dyscreen]
    wavelength = wavelength of the light [m]
    I0 = intensity of the light at the mask plane [W/m2]
    
    ------- 
    optional 
    verbose, defaults to False, if True prints 
    """
    import decimal
    import numpy as np 
    import matplotlib.pyplot as plt 
    
    # set the precision to double that of float64.. or whatever you want.
    decimal.setcontext(decimal.Context(prec=34))

    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    if nps <= 2*npm: 
        print("The number of screen is not large enough, I will resize them for you.")
        nps = npm*4 
        print("Rescaled the screen pixels to "+str(nps)+" . The computation will now proceed")
    
    #size of mask 
    dmask = pixsizemask * npm
    
    #physical constants 
    c_const = 3e8 #m/s
    eps0 = 8.85e-12 #F/m
    n = 1 #refractive index of medium  

    k = 2* np.pi/wavelength

    ## definitions 
    unit = np.ones((npm,npm), dtype=complex)
    r = np.zeros((npm,npm)) 
    prop1 = np.zeros((npm,npm))
    prop2 = np.zeros((npm,npm))
    propE = np.zeros((npm,npm))
    
    #electric field real and imaginary and total 
    rEs =np.zeros ((nps,nps))
    iEs =np.zeros ((nps,nps))
    Escreen =np.zeros ((nps,nps), dtype =complex)

    #define the zpos of the mask at 0 
    zm =0 

    #definition of the mask
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1,ym1)
    
    ##Electric field calc from intensity
    E0 = np.sqrt(2*I0/(c_const*n*eps0))
    
    #definition of the electric field at the mask 
    E0m = E0 * mask

    fig=plt.figure()
    plt.imshow(E0m, vmin=0, vmax=1, cmap=plt.get_cmap("jet"))
    plt.title("Efield at mask")
    
    xs1 = np.linspace(-dxscreen, dxscreen, nps)
    ys1 = np.linspace(-dyscreen, dyscreen, nps)
    (xs, ys) = np.meshgrid(xs1,ys1)

    ###### calculate the Rayleigh Sommerfeld integral 
    #From Oshea formulation, eq 2.8
    for isc in np.arange(0,nps-1):
        if verbose == True: 
            print(isc/nps)
            
        for jsc in np.arange(0,nps-1): 
            r = np.sqrt((xs[isc,jsc]-xm)**2 + (ys[isc,jsc]-ym)**2 + (zs-zm)**2)
            r2 = r*r
            prop1= np.exp(-r*1.0j*k)/r2
            prop2 = zs * (1.0j * k  + unit/r)
            propE = E0m * prop1 * prop2
            #here npm*400, is a guess for the number of points to calc the int
            rEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*400,npm*400,np.real(propE))/(2*np.pi)
            iEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*400,npm*400,np.imag(propE))/(2*np.pi)

    Escreen = rEs + 1.0j*iEs 
    Iscreen = (c_const*n*eps0/2) * np.abs(Escreen)**2
    iplot = 10*Iscreen**0.2
    iplotmax = np.max(iplot)
    
    return Escreen, Iscreen, iplot 
    

def RS_int_XZ2(zs, nzds, mask, npixmask, pixsizemask, npixscreen, dxscreen, dyscreen, wavelength, I0, verbose=False, logscale = False): 
    """
    Calculates the RS_int in the  of the first kind, taking information about mask, distance to screen, and screen information
    returns Escreen (complex electric field at obs screen), Iscreen (intensity at obs screen), iplot (the actual intensity) 
    inputs: 
    zs = distance to screen [m]
    nzds = number of observation planes along z 
    mask = image of the mask (normalized to [0,1]) - in principle could be grayscale between 0 and 1 
    npixmask = number of pixels on the side of the mask 
    pixsizemask = size of pixel of the mask [m]
    npixscreen = number of pixels on the side of the screen 
    dxscreen = x side of the screen [m]
    dyscreen = y side of the screen [m]
    wavelength = wavelength of the light [m]
    I0 = intensity of the light at the mask plane [W/m2]
    
    logscale makes the calculation of the points in a log scale 
    """

    # set the precision to double that of float64
    decimal.setcontext(decimal.Context(prec=34))

    #number of pixels 
    nps = npixscreen 
    npm = npixmask 
    
    if nps <= 2*npm: 
        print("The number of screen is not large enough, I will resize them for you.")
        nps = npm*4 
        print("Rescaled the screen pixels to "+str(nps)+" . The computation will now proceed")
    
    #size of mask 
    dmask = pixsizemask * npm
    
    #physical constants 
    c_const = 3e8 #m/s 
    eps0 = 8.85e-12 #F/m
    n = 1 #refractive index of medium  

    k = 2* np.pi/wavelength

    ## definitions 
    unit = np.ones((npm,npm), dtype=complex)
    r = np.zeros((npm,npm)) 
    prop1 = np.zeros((npm,npm))
    prop2 = np.zeros((npm,npm))
    propE = np.zeros((npm,npm))

    
    #electric field real and imaginary and total 
    rEs =np.zeros ((nps,nps))
    iEs =np.zeros ((nps,nps))
    Escreen =np.zeros ((nps,nps), dtype =complex)

    #define the zpos of the mask at 0 
    zm =0 

    #definition of the mask
    xm1 = np.linspace(-dmask/2, dmask/2, npm)
    ym1 = np.linspace(-dmask/2, dmask/2, npm)
    (xm, ym) = np.meshgrid(xm1,ym1)
    
    ##Electric field calc from intensity
    E0 = np.sqrt(2*I0/(c_const*n*eps0))
    
    #definition of the electric field at the mask
    ##Atention here is being done the amplitude only... what about phase
    E0m = E0 * mask

    fig=plt.figure()
    plt.imshow(E0m, vmin=0, vmax=1, cmap=plt.get_cmap("jet"))
    plt.title("Efield at mask")

 
    xs1 = np.linspace(-dxscreen, dxscreen, nps)
    ys1 = np.linspace(-dyscreen, dyscreen, nps)
    (xs, ys) = np.meshgrid(xs1,ys1)
    
    #print(ys[0,:])
    
    
    #array with the distance in z until the distance zs given as argument 
    zdistarray = np.linspace(0,zs,nzds)
    
    if logscale ==True: 
        zdistarray = np.logspace(np.log(1e-6),np.log(zs),nzds, base=np.e)
    
    inten = np.zeros((len(xs1), len(zdistarray)))
    
    print (inten)

    ###### calculate the Rayleigh Sommerfeld integral 
    for iz, zsd in enumerate(zdistarray): 
        ##### calculate the Rayleigh Sommerfeld integral 
        #From Oshea formulation, eq 2.8
        for isc in np.arange(0,nps-1):
            if verbose == True: 
                print(isc/nps)
                
            for jsc in np.arange(0,nps-1): 
                r = np.sqrt((xs[isc,jsc]-xm)**2 + (ys[isc,jsc]-ym)**2 + (zsd-zm)**2)
                r2 = r*r
                prop1= np.exp(-r*1.0j*k)/r2
                prop2 = zsd * (1.0j * k  + unit/r)
                propE = E0m * prop1 * prop2
                #here npm*400, is a guess for the number of points to calc the int
                rEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*400,npm*400,np.real(propE))/(2*np.pi)
                iEs[isc,jsc] = double_Integral(-dmask/2, dmask/2, -dmask/2, dmask/2, npm*400,npm*400,np.imag(propE))/(2*np.pi)

        Escreen = rEs + 1.0j*iEs 
        Iscreen = (c_const*n*eps0/2) * np.abs(Escreen)**2
        iplot = 10*Iscreen**0.2
        #print(iplot)
        midpoint = int(nps/2)-1
        inten[:,iz] = iplot[:,midpoint]#This is a line cut in y  for each iteration
        #inten[isc,iz] = iplot
        
        print(inten[:,iz])
        
        iplotmax = np.max(iplot)
    
    return Escreen, Iscreen, inten
   
   
def double_Integral(xmin, xmax, ymin, ymax, nx, ny, A):
    """
    #rudimentary 2D integral, following https://stackoverflow.com/questions/20668689/integrating-2d-samples-on-a-rectangular-grid-using-scipy 
    """

    dS = ((xmax-xmin)/(nx-1)) * ((ymax-ymin)/(ny-1))

    A_Internal = A[1:-1, 1:-1]

    # sides: up, down, left, right
    (A_u, A_d, A_l, A_r) = (A[0, 1:-1], A[-1, 1:-1], A[1:-1, 0], A[1:-1, -1])

    # corners
    (A_ul, A_ur, A_dl, A_dr) = (A[0, 0], A[0, -1], A[-1, 0], A[-1, -1])

    return dS * (np.sum(A_Internal)\
                + 0.5 * (np.sum(A_u) + np.sum(A_d) + np.sum(A_l) + np.sum(A_r))\
                + 0.25 * (A_ul + A_ur + A_dl + A_dr))

    
    

def circ_fraun(aperture_rad, rcoord, zdist, wavelength): 
    """
    Analytical Fraunhofer approximation for circular aperture
    implements circular aperture fraunhofer approximation analytic solution for x-y planes, Goodman exp 4-31 
    inputs: 
    aperture_rad = aperture radius in m
    rcoord = radius coordinate 2D map (x-y plane)
    zdist = distance to screen in m 
    wavelength = wavelength in m 
    
    returns intensity 2D array (x-y)
    """

    from scipy.special import j1
    ###use like jv(v,z) where v is the order and z the dist in z 
    
    area = np.pi * aperture_rad**2
    
    k = 2* np.pi/wavelength
    argument = k*aperture_rad*rcoord/zdist
    
    #function 4-31
    intensity = (area/(wavelength * zdist))**2 * (2*j1( argument)/argument)**2
    
    return intensity 



def circ_zz(aperture_rad, zdist, wavelength):
    """
    implements circular aperture fresnel approximation analytic solution for x-z planes, Sheppard  1992, exp28 
    ###This function is the axial intensity for a circular aperture following 1992 Sheppard paper, exp (28) 
    inputs: 
    aperture_rad = aperture radius in m
    zdist = distance to screen in m 
    wavelength = wavelength in m 
    
    returns intensity 2D array (x-y)
    
    """
    
    k = 2*np.pi /wavelength 
    
    izz1 = (1 + np.sqrt(1 + aperture_rad**2 /zdist**2))
    izz2 = 1 + aperture_rad**2/(2*zdist**2)
    izz3 = (k * aperture_rad**2/(2 *zdist))/(np.sqrt(1+aperture_rad**2/zdist**2)+1)
    itot = 0.25*(izz1/izz2) * np.sin(izz3)**2
    
    return itot


def circ_zz24(aperture_rad, zdist, wavelength):
    """
    implements circular aperture fresnel approximation analytic solution for x-z planes, Sheppard  1992 , exp24
    ###This function is the axial intensity for a circular aperture following 1992 Sheppard paper, exp (24) 
    inputs: 
    aperture_rad = aperture radius in m
    zdist = distance to screen in m 
    wavelength = wavelength in m 
    
    returns intensity 2D array (x-y)
    
    """
    
    k = 2*np.pi /wavelength 
    
    izz1 = 1/(1+aperture_rad**2/zdist**2)
    izz2 = 2/np.sqrt(1+aperture_rad**2/zdist**2)
    izz3 = (k * aperture_rad**2/(zdist)) /(np.sqrt(1+aperture_rad**2/zdist**2)+1)
    itot = 0.25*(1+izz1-izz2* np.cos(izz3))
    
    return itot
    
    
    
def rect_fraun(sizex, sizey, xcoord, ycoord, zd, wl): 
    """
    Analytical Fraunhofer approximation for rectangular aperture
    implements rectangular aperture fraunhofer approximation analytic solution for x-y planes, Goodman exp 4-28
    inputs: 
    aperture_rad = aperture radius in m
    rcoord = radius coordinate 2D map (x-y plane)
    zdist = distance to screen in m 
    wavelength = wavelength in m 
    
    returns intensity 2D array (x-y)
    """
    from scipy.special import sinc 
    #sinc uses sin(pi*x)/(pi*x) with x as the argument 
    
    area = sizex*sizey
    
    k= 2* np.pi/wl
    
    intensity = (area/(wl*zd))**2 * sinc(sizex*xcoord/(wl*zd))**2 * sinc(sizey * ycoord/(wl* zd))**2
    
    return intensity 
    


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