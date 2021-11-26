####export.py 

###Function exports an image file (converted to gray) into a gds file 
def grayim2gds(infile, outfile, pixelx, pixely, cellname, level, layer=0, datatype=0, verbose=False):
    """
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 for reading the image
    'infile'  = input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
    'outfile' = output GDS file, e.g. "image.gds"
    'pixelx'  = pixel size in x, in um 
    'pixely'  = pixel size in y, in um 
    'cellname'= string with cellname, e.g. "TOP"
    'level'   = int from 0 to 255, pixels gray value to be passed to GDS 
    
    optional:
    'layer' = gray level, defaults to 0 
    'datatype' = gds datatype, defaults to 0 
    'verbose' defaults to False, if True prints 
    
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
        if verbose == True: 
            print(i/h)
        for j in np.arange(0, w):
            #print(j/w)
            #here we can also think of selectin pixels at a certain level only
            #and creating a GDS from a grayscale image 
            if img[i][j] == int(level):
                #print(i)
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

    
def grayim2gds_writer_frac(infile, outfile, pixelx, pixely, cellname, level, nm=None, layer=0, datatype=0 , verbose=False):
    """
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 for reading the image
    If nm is given fractionates the image 
    by default adds the image to (layer, datatype) = (0,0)
    'infile'  = input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
    'outfile' = output GDS file, e.g. "image.gds"
    'pixelx'  = pixel size in x, in um 
    'pixely'  = pixel size in y, in um 
    'cellname'= string with cellname, e.g. "TOP"
    'level'   = int from 0 to 255 (0 = black, 255=white) , pixels gray value to be passed to GDS

    optional:
    'nm' = nr of pixels for each fractioned part  (should be multiple of pixel sizes), defaults to None
    'layer' = gray level, defaults to 0 
    'datatype' = gds datatype, defaults to 0 
    'verbose' defaults to False, if True prints 
    
    ---- 
    Usage example: 
    
    infilxe = "image.png"
    outfilxe = "image.gds"
    pixelx = 1 #um 
    pixely = 1 #um 
    cellname = "TOP" #name of the gds cell 
    graycolor = 0 #black pixels 
    frac = 250 #size of frac pixels in the image 

    grayim2gds_writer_frac(infilxe, outfilxe, pixelx, pixely, cellname, graycolor, frac, verbose=True) 
    """
    import cv2
    import gdspy 
    import numpy as np 

    img = cv2.imread(infile, cv2.IMREAD_GRAYSCALE)
    
    if img is not None: 
        print("Sucessfully imported img!")
        
    h,w = img.shape 
    print(h)
    print(w)
    
    if nm == None:
        nmx = w
        nmy = h
    else: 
        nmx = nm
        nmy = nm
    
    harray = np.arange(0,h+1,nmy)
    warray = np.arange(0,w+1,nmx)
    cn = 1
    print(harray)
    
    for hn, hi in enumerate(harray):
        if hn == (len(harray)-1):

            break
        print(hi)
        for hw, wi in enumerate(warray):
            if hw == (len(warray)-1): 
                break
                
            print(wi)
            
            pols=[]

            gdspy.current_library = gdspy.GdsLibrary() 
            lib = gdspy.GdsLibrary()
            outfilen = str(cn)+outfile
            writer = gdspy.GdsWriter(outfilen,unit=1.0e-6,precision=1.0e-9)
            cell = lib.new_cell(cellname)

            for i in np.arange(hi,hi+nmy):
                if verbose == True: 
                    print(i/h)
                for j in np.arange(wi, wi+nmx):
                    #print(j/w)
                    #here we can also think of selectin pixels at a certain level only
                    #and creating a GDS from a grayscale image 
                    if img[i][j] == int(level):
                        #rectangle takes the two  opposite corners 
                        #pols.append(gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype))
                        rect = gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype)
                        cell.add(rect)
                   
                writer.write_cell(cell)

            writer.close()
            cn = cn+1


    print(cn)
    print("Exported the image file "+str(infile) + " into " + str(outfile))

  
def grayim2gds_writer(infile, outfile, pixelx, pixely, cellname, level, layer=0, datatype=0 , verbose=False):
    """
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 for reading the image
    by default adds the image to (layer, datatype) = (0,0)
    'infile'  = input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
    'outfile' = output GDS file, e.g. "image.gds"
    'pixelx'  = pixel size in x, in um 
    'pixely'  = pixel size in y, in um 
    'cellname'= string with cellname, e.g. "TOP"
    'level'   = int from 0 to 255 (0 = black, 255=white) , pixels gray value to be passed to GDS 

    optional:
    'layer' = gray level, defaults to 0 
    'datatype' = gds datatype, defaults to 0 
    'verbose' defaults to False, if True prints 
    
    ---- 
    Usage example: 
    
    infilxe = "image.png"
    outfilxe = "image.gds"
    cellname = "TOP" #name of the gds cell 
    graycolor = 0 #black pixels 
    pixelx = 1 #um 
    pixely = 1 #um 

    grayim2gds_writer(infilxe, outfilxe, pixelx, pixely,cellname, graycolor, verbose=True)"""
    import cv2
    import gdspy 
    import numpy as np 

    img = cv2.imread(infile, cv2.IMREAD_GRAYSCALE)
    
    if img is not None: 
        print("Sucessfully imported img!")
        
    h,w = img.shape 
    print(h)
    print(w) 
    
    nmx = w
    nmy = h
    
    harray = np.arange(0,h+1,nmy)
    warray = np.arange(0,w+1,nmx)
    #print(harray)
    
    lib = gdspy.GdsLibrary()
    gdspy.current_library = gdspy.GdsLibrary() 

    outfilen = outfile
    writer = gdspy.GdsWriter(outfilen,unit=1.0e-6,precision=1.0e-9)
    cell = lib.new_cell(cellname)
    
    for hn, hi in enumerate(harray):
        if hn == (len(harray)-1):
            #writer.close()
            break
        #print(hi)
        for hw, wi in enumerate(warray):
            if hw == (len(warray)-1): 
                break   
            #print(wi)

            for i in np.arange(hi,hi+nmy):
                if verbose == True: 
                    print(i/h)
                
                for j in np.arange(wi, wi+nmx):
                    #print(j/w)
                    #here we can also think of selectin pixels at a certain level only
                    #and creating a GDS from a grayscale image 
                    if img[i][j] == int(level):
                        #rectangle takes the two  opposite corners 
                        #pols.append(gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype))
                        rect = gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype)
                        
                        cell.add(rect)
                
                writer.write_cell(cell)

        writer.close()

    print("Exported the image file "+str(infile) + " into " + str(outfile))




###OTHERS 




