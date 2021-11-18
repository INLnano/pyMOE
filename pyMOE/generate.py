####generate.py 

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
    
    #by default the aperture is at the center of the mask 
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
     
    #by default the aperture is at the center of the mask 
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

def fzp_mask(npix, foc, lda, xsiz, ysiz, filename, plotting=False ):
    """
    returns a fresnel zone plate
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

    #by default the aperture is at the center of the mask 
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

##Code to create gray scale 
def create_scale(npixel, nsz, ngs): 
    """
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
