####generate.py 

#import the gds operations from the gds_klops file 

from gds_klops.py import * 

####Function that defines circular aperture mask 
def circ_mask(npix, pixsize, partial, filename, plotting=False ):
    """
    returns 2D circular aperture mask 
    npix = nr of pixels 
    pixsize = pixel size 
    partial = size of circ radius as fraction of npix [0,1]
    filename = string with image name 'image.png'
    if plotting=True, shows the mask 
    
    """
    
    from matplotlib import pyplot as plt
    import numpy as np 

    xsiz = npix
    ysiz = npix
    
    #by default centered  
    xcmm =  0.5* xsiz
    ycmm =  0.5* ysiz 

    a = partial*npix #radius of the circular aperture 
    maskcir = np.ones((npix,npix))
    xc1 = np.linspace(0, xsiz, npix)
    yc1 = np.linspace(0, ysiz, npix)
    (xc, yc) = np.meshgrid(xc1,yc1)
    
    #definition of the circular aperture 
    rc = np.sqrt((xc-xcmm)**2 + (yc-ycmm)**2)
    maskcir[np.where(rc>a)] = 0

    if plotting == True: 
        fig=plt.figure()
        plt.imshow(maskcir, vmin=0, vmax=1, cmap=plt.get_cmap("Greys"))
        plt.show()
        fig.savefig(filename)

    return maskcir 


####Function that defines rectangular aperture mask 
def rect_mask(npix, pixsize, partial, filename, plotting=False ):
    """
    returns 2D rectangular aperture mask 
    npix = nr of pixels 
    pixsize = pixel size 
    partial = size of circ radius as fraction of npix [0,1]
    filename = string with image name 'image.png'
    if plotting=True, shows the mask 
    
    """
    import numpy as np 
    from matplotlib import pyplot as plt 
    
    xsiz = npix
    ysiz = npix
     
    #by default centered  
    xcmm =  0.5* xsiz
    ycmm =  0.5* ysiz 
    
    maskrect = np.zeros((npix,npix))
    xc1 = np.linspace(0, xsiz, npix)
    yc1 = np.linspace(0, ysiz, npix)
    (xc, yc) = np.meshgrid(xc1,yc1)
    
    ##definition of the mask pixels
    maskrect[np.where((xc>(xcmm-xcmm*partial*2)) & (xc<(xcmm+xcmm*partial*2)) & \
                      (yc>(ycmm-ycmm*partial*2)) & (yc<(ycmm+ycmm*partial*2)))] = 1

    if plotting == True: 
        fig=plt.figure()
        plt.imshow(maskrect, vmin=0, vmax=1, cmap=plt.get_cmap("Greys"))
        plt.show()
        fig.savefig(filename)
    
    return maskrect
    
####Function that defines a Fresnel Zone Plate mask 
def fzp_mask(npix, foc, lda, xsiz, ysiz, filename, plotting=False ):
    """
    returns a fresnel zone plate (as a numpy 2D array)
    npix = nr of pixels 
    foc = focal length in um
    lda = wavelength in um 
    xsiz = size in x in um 
    ysiz = size in y in um 
    filename = string with mask image name 'image.png'
    if plotting=True, shows the mask 
    
    Example of use: 
    
    fzp_mask(npix = 50,\
         foc = 5000 ,\
         lda = 0.6328 ,\
         xsiz = 500, \
         ysiz = 500, \
         filename = 'fresnel2.png', \
         plotting=True )
         
    """
    from matplotlib import pyplot as plt
    import numpy as np 

    #by default centered 
    xcmm =  0.5* xsiz
    ycmm =  0.5* ysiz 

    a = 0.5 * np.min([xsiz,ysiz])  #radius of the circular aperture 
    maskfres = np.ones((npix,npix))
    xc1 = np.linspace(0, xsiz, npix)
    yc1 = np.linspace(0, ysiz, npix)
    (xc, yc) = np.meshgrid(xc1,yc1)

    #definition of the circular aperture 
    rc = np.sqrt((xc-xcmm)**2 + (yc-ycmm)**2)
    
    #definition of the phase profile 
    fzp = np.exp(-1.0j*(foc-np.sqrt(foc**2 + rc**2))*(2*np.pi)/(lda))

    #Define the zones 
    fzp[np.where((np.angle(fzp)>-np.pi/2 )& (np.angle(fzp)<np.pi/2) )] = 0 

    i,j = fzp.shape 

    #final plate
    fzp2 = np.zeros((i,j)) 

    for ie in np.arange(0,i):
        for je in np.arange(0,j):
            if ((np.angle(fzp[ie][je]) >= -np.pi/2) & (np.angle(fzp[ie][je]) <= np.pi/2)): 
                fzp2[ie][je] = 0
            else: 
                fzp2[ie][je] = 1

    fzp2[np.where(rc>a)] = 1

    if plotting == True: 
        fig=plt.figure()
        plt.imshow(fzp2, cmap=plt.get_cmap("Greys"))
        #plt.colorbar()
        plt.show()
        fig.savefig(filename)
        
    return fzp2    

##Code to create a gray scale with successive gray levels 
def create_scale(npixel, nsz, ngs): 
    """
    returns a 2D array with a scale of successive gray levels 
    npixel= nr of pixels 
    nsz = division in size 
    ngs = nr of gray levels 
    
    """
    import cv2 
    import numpy as np

    scale_img = np.zeros((npixel,npixel,3), np.uint8)

    width = npixel 
    height = npixel 
    
    nsz = npixel/ngs 
    
    xdims = np.arange(0,width, nsz)
    
    xdist = 255/ngs
    xlevs = np.arange(0,255, xdist)
    print(xlevs)

    print(xdims)

    gslev = ngs

    xdimsint = np.array(xdims, dtype=int)
    xdimsround = np.round(xdimsint)

    for iw, wd in enumerate(xlevs):
        if iw == (len(xdims)-1): 
            break 
            
        gss =255 - np.round(wd)
        
        colorgray = np.uint8([[gss,gss,gss]])
        
        graycolor = cv2.cvtColor(colorgray, cv2.COLOR_GRAY2RGB)
        
        scale_img[:,xdimsround[iw]:xdimsround[iw+1]] = (int(graycolor[0][0][0]),int(graycolor[0][0][1]),int(graycolor[0][0][2]))      # (B, G, R)
        
    
    fig1 = plt.figure()
    plt.imshow(scale_img, vmin=0, vmax=255, cmap=plt.get_cmap("Greys"))
    plt.title("scaled")
 
    return scale_img
    
    
def lensfres(x,y,x0,y0,fo,lda): 
    """
    returns the COMPLEX PHASE of a fresnel lens ith input meshgrid (x,y) with center at (x0,y0)
    x = x array from meshgrid 
    y = y array from meshgrid 
    x0 = coordinate of center of the lens 
    y0 = coordinate of center of the lens
    fo = focal distance 
    lda = wavelength 
    
    Note: for angle (in rad), call numpy.angle(...)
    """
    import numpy as np
    
    rc = np.sqrt((x-x0)**2 + (y-y0)**2)
    fresn = np.exp(1.0j*(fo-np.sqrt(fo**2 + rc**2))*(2*np.pi)/(lda))
    return fresn 
    
    
def spiral(x,y,x0,y0,L):
    """
    returns a spiral COMPLEX PHASE with input meshgrid (x,y) with center at (x0,y0)
    x = x array from meshgrid 
    y = y array from meshgrid 
    x0 = x-coordinate of center of the lens 
    y0 = y-coordinate of center of the lens
    L = topological charge 
    """
    import numpy as np 

    theta = np.arctan2((y-y0),(x-x0))
    sp = np.exp(1.0j*L*theta)
    return sp


def fresnel_phase_mask(npix, foc, lda, xsiz, ysiz,n, filename=None, plotting=False ,prec = 1e-6, mpoints = 1e9):
    """
    returns a Fresnel "phase mask" (2D array of the phase IN RADIANS)
    parameters: 
    npix = nr of pixels , by default the results 2D array is npix by npix 
    foc = focal length in um
    lda = wavelength in um 
    xsiz = size in x in um 
    ysiz = size in y in um
    n = number of gray levels 
    
    optional: 
    filename = string with mask output into GDS  (default None)
    plotting = True, shows the mask  (default False)
    prec = precision of the gdspy boolean operation  (default 1e-6)
    mpoints = max_points of the gdspy polygon (default 1e9)
    
    Example of use: 
    fresnel_phase_mask(npix = 5000, \
                   foc = 5000,\
                   lda = 0.6328 ,\
                   xsiz = 500,\
                   ysiz =500,\
                   n=10,\
                   filename='fresnel_phase_mask.gds',\
                   plotting=True )   #Should take around ~30 s 
         
    """  
    import numpy as np
    import matplotlib.pyplot as plt 
    
    #by default centered 
    xcmm =  0.5* xsiz
    ycmm =  0.5* ysiz 
    
    a = 0.5 * np.min([xsiz,ysiz])  #radius of the circular aperture 
    maskfres = np.ones((npix,npix))
    xc1 = np.linspace(0, xsiz, npix)
    yc1 = np.linspace(0, ysiz, npix)
    (xc, yc) = np.meshgrid(xc1,yc1)
    
    #definition of the circular aperture 
    rc = np.sqrt((xc-xcmm)**2 + (yc-ycmm)**2)

    #calculate the fresnel complex phase 
    ###TODO: Generalize the input function for any phase map function, given as argument to the function 
    fresarray = lensfres(xc,yc,xcmm,ycmm,foc,lda)
    
    fresarray[np.where(rc>a)] = np.pi
    fresarray_rad = np.angle(fresarray)
    
    #make array with the z plane intersections  (n gray levels)
    zlevs = np.linspace(np.min(fresarray_rad), np.max(fresarray_rad), n+1)
    #print(zlevs)

    if plotting == True: 
        plt.figure()
        plt.axis('equal')
        cs = plt.contourf(xc,yc,fresarray_rad, zlevs, cmap=plt.get_cmap("Greys"))
        plt.xlabel('x ($\mu$m)')
        plt.ylabel('y ($\mu$m)')
        plt.colorbar(label='Phase (rad)')
        plt.tight_layout()
    else: 
        cs = plt.contourf(xc,yc,fresarray_rad, zlevs, cmap=plt.get_cmap("Greys"))
      
    #possible improvement, pass this function as argument
    #lib1, cell1 = cell_wpol_gdspy_fast(cs, 'TOP', prec, mpoints)
    lib1, cell1 = cell_wpol_gdspy(cs, 'TOP', prec, mpoints)

    if filename is not None: 
        lib1.write_gds(filename)
        print("Saved the phase profile with " + str(n) +  " layers into the file " + filename)
        
    return fresarray_rad 

###ANY FUNCTION PHASE MASK 
def phase_mask(npix, xsiz, ysiz, n, fname,*args,filename=None, plotting=False ,prec = 1e-6, mpoints = 1e9 ,**kwargs):
    """
    returns a "phase mask" (2D array of the phase IN RADIANS) from COMPLEX PHASE function fname  given as argument
    
    parameters: 
    npix = nr of pixels (or points) , by default the results 2D array is npix by npix 
    xsiz = size in x in um 
    ysiz = size in y in um
    n = number of gray levels
    fname = function name (e.g. lensfres([x,y,x0,y0], args) , where [x,y,x0,y0] are fixed)
    *args = arguments fname, excluding the [x,y,x0,y0], e.g. with syntax fo = ..., lda ..., etc 
    
    optional: 
    filename = string with mask output into GDS  (default None)
    plotting = True, shows the mask  (default False)
    prec = precision of the gdspy boolean operation  (default 1e-6 um)
    mpoints = max_points of the gdspy polygon (default 1e9 points)
    
    Examples of use: #Should take around ~30 s for any of these 
    phase_mask(5000, 500,500, 10,\
           lensfres, fo=5000, lda=0.6328, \
           filename="fresnel_phase_plate.gds", plotting=True ,prec = 1e-6, mpoints = 1e9 )
           
    phase_mask(5000, 500,500, 60,\
           spiral, L=1, \
           filename="spiral_phase_plate.gds", plotting=True ,prec = 1e-12, mpoints = 1e9 )
         
    """  
    import numpy as np
    import matplotlib.pyplot as plt
    
    #by default centered 
    xcmm =  0.5* xsiz
    ycmm =  0.5* ysiz 
    
    a = 0.5 * np.min([xsiz,ysiz])  #radius of the circular aperture 
    maskfres = np.ones((npix,npix))
    xc1 = np.linspace(0, xsiz, npix)
    yc1 = np.linspace(0, ysiz, npix)
    (xc, yc) = np.meshgrid(xc1,yc1)
    
    #definition of the circular aperture 
    rc = np.sqrt((xc-xcmm)**2 + (yc-ycmm)**2)

    #calculate the complex phase  fname function  
    farray = fname(xc,yc,xcmm,ycmm,*args, **kwargs)
    
    #farray[np.where(rc>a)] = np.pi
    farray_rad = np.angle(farray)
    
    #make array with the z plane intersections  (n gray levels)
    zlevs = np.linspace(np.min(farray_rad), np.max(farray_rad), n+1)
    #print(zlevs)

    if plotting == True: 
        plt.figure()
        plt.axis('equal')
        cs = plt.contourf(xc,yc,farray_rad, zlevs, cmap=plt.get_cmap("Greys"))
        plt.xlabel('x ($\mu$m)')
        plt.ylabel('y ($\mu$m)')
        plt.colorbar(label='Phase (rad)')
        plt.tight_layout()
      
    #possible improvement, pass this function as argument
    #lib1, cell1 = cell_wpol_gdspy_fast(cs, 'TOP', prec, mpoints)
    lib1, cell1 = cell_wpol_gdspy(cs, 'TOP', prec, mpoints)

    if filename is not None: 
        lib1.write_gds(filename)
        print("Saved the phase profile with " + str(n) +  " layers into the file " + filename)
        
    return farray_rad 
