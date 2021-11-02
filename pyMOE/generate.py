####generate.py 

###Code to create gray scale 
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

