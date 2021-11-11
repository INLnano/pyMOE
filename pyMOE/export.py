####export.py 

###Function exports an image file (converted to gray) into a gds file 
def grayim2gds(infile, outfile, pixelx, pixely, cellname, level, layer=0, datatype=0):
    """
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 for reading the image
    'infile'  = input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
    'outfile' = output GDS file, e.g. "image.gds"
    'pixelx'  = pixel size in x, in um 
    'pixely'  = pixel size in y, in um 
    'cellname'= string with cellname, e.g. "TOP"
    'level'   = int from 0 to 255, pixels gray value to be passed to GDS 
    
    ---- 
    Usage example: 
    
    infilxe = "image.png"
    outfilxe = "image.gds"
    pixelx = 1 #um 
    pixely = 1 #um 
    grayim2gds(infilxe, outfilxe, pixelx, pixely,"TOP", 0)
    """
    import cv2
    import gdspy 
    import numpy as np 

    img = cv2.imread(infile, cv2.IMREAD_GRAYSCALE)
    
    h,w = img.shape 
    #print(h)
    #print(w)

    pols=[]
    for i in np.arange(0,h):
        #print(i/h)
        for j in np.arange(0, w):
            #print(j/w)
            #here we can also think of selectin pixels at a certain level only
            #and creating a GDS from a grayscale image 
            if img[i][j] == int(level):
                #rectangle takes the two  opposite corners 
                pols.append(gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype))

    if len(pols) !=0: 
        polygons=gdspy.boolean(pols[0], pols[1:], "or") #max_points could be used

        #define current as Gdslib with default properties unit=1e-06, precision=1e-09
        gdspy.current_library = gdspy.GdsLibrary() 
        cell = gdspy.Cell(cellname)
        cell.add(polygons)
        gdspy.write_gds(outfile)
        print("Exported the image file "+str(infile) + " into " + str(outfile))
        
    else: 
        print("There are no pixels in this gray level! Please try another gray level.")








###OTHERS 




