"""
export.py 
Module containing several functions to export masks to gds. 

"""
import cv2
# import gdspy 
import numpy as np 
import gdspy

    
    
###Function exports an image file (converted to gray) into a gds file 
def grayim2gds(infile, outfile, pixelx, pixely, cellname, level, layer=0, datatype=0, verbose=False):
    """
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 
    
    Args: 
        :infile:    input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
        :outfile:   output GDS file, e.g. "image.gds"
        :pixelx:    pixel size in x, in um 
        :pixely:    pixel size in y, in um 
        :cellname:  string with cellname, e.g. "TOP"
        :level:     int from 0 to 255, pixels gray value to be passed to GDS 
        :layer:     (optional) gray level, defaults to 0 
        :datatype:  (optional) gds datatype, defaults to 0 
        :verbose:   (optional) defaults to False, if True prints 
    
    ---- 
    Usage example: 
    
    infilxe = "image.png"
    outfilxe = "image.gds"
    pixelx = 1 #um 
    pixely = 1 #um 
    grayim2gds(infilxe, outfilxe, pixelx, pixely,"TOP", 0)
    """
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
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 
    
    By default adds the image to (layer, datatype) = (0,0)
    
    Args:
        :infile:    input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
        :outfile:   output GDS file, e.g. "image.gds"
        :pixelx:    pixel size in x, in um 
        :pixely:    pixel size in y, in um 
        :cellname:  string with cellname, e.g. "TOP"
        :level:     int from 0 to 255 (0 = black, 255=white) , pixels gray value to be passed to GDS
        :nm:        (optional) If nm is given fractionates the image, nr of pixels for each fractioned part  (should be multiple of pixel sizes), defaults to None
        :layer:     (optional) gray level, defaults to 0 
        :datatype:  (optional) gds datatype, defaults to 0 
        :verbose:   (optional) defaults to False, if True prints 
    
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
    #print(harray)
    
    for hn, hi in enumerate(harray):
        if hn == (len(harray)-1):

            break
        #print(hi)
        for hw, wi in enumerate(warray):
            if hw == (len(warray)-1): 
                break
                
            #print(wi)
            
            pols=[]

            gdspy.current_library = gdspy.GdsLibrary() 
            lib = gdspy.GdsLibrary()
            outfilen = str(cn)+outfile
            writer = gdspy.GdsWriter(outfilen,unit=1.0e-6,precision=1.0e-9)
            cell = lib.new_cell(cellname)

            for i in np.arange(hi,hi+nmy):
                cell.remove_polygons(lambda pts, layer, datatype: layer == 0)
                if verbose == True: 
                    print(i/h)
                for j in np.arange(wi, wi+nmx):
                    #print(j/w)
                    #here we can also think of selectin pixels at a certain level only
                    #and creating a GDS from a grayscale image 
                    if img[i][j] == int(level):
                        #rectangle takes the two  opposite corners 
                        pols.append(gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype))

            cell.add(pols)
            writer.write_cell(cell)
            del cell
            writer.close()
            cn = cn+1

    #print(cn)
    print("Exported the image file "+str(infile) + " into " + str(outfile))

  
def grayim2gds_writer(infile, outfile, pixelx, pixely, cellname, level, layer=0, datatype=0 , verbose=False):
    """
    (void) Transforms one image (converted to grayscale) into a gds, using cv2 
    by default adds the image to (layer, datatype) = (0,0)
    
    Args:
        :infile:    input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
        :outfile:   output GDS file, e.g. "image.gds"
        :pixelx:    pixel size in x, in um 
        :pixely:    pixel size in y, in um 
        :cellname:  string with cellname, e.g. "TOP"
        :level:     int from 0 to 255 (0 = black, 255=white) , pixels gray value to be passed to GDS 
        :layer:     (optional) gray level, defaults to 0 
        :datatype:  (optional) gds datatype, defaults to 0 
        :verbose:   (optional) defaults to False, if True prints 
    
    ---- 
    Usage example: 
    
    infilxe = "image.png"
    outfilxe = "image.gds"
    cellname = "TOP" #name of the gds cell 
    graycolor = 0 #black pixels 
    pixelx = 1 #um 
    pixely = 1 #um 

    grayim2gds_writer(infilxe, outfilxe, pixelx, pixely,cellname, graycolor, verbose=True)"""


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

    pols = []
    
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
                cell.remove_polygons(lambda pts, layer, datatype: layer == 0)
                if verbose == True: 
                    print(i/h)
                
                for j in np.arange(wi, wi+nmx):
                    #print(j/w)
                    #here we can also think of selectin pixels at a certain level only
                    #and creating a GDS from a grayscale image 
                    if img[i][j] == int(level):
                        #rectangle takes the two  opposite corners 
                        pols.append(gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype))

    cell.add(pols)         
    writer.write_cell(cell)
    del cell 

    writer.close()

    print("Exported the image file "+str(infile) + " into " + str(outfile))

def grayim2gds_writer_klops(infile,  output_filename , pixelx, pixely, cellname, level, layer=0, datatype=0 , verbose=False):
    """
    (void)  Transforms one image (converted to grayscale) into a gds, using cv2 for reading the image
    by default adds the image to (layer, datatype) = (0,0)
    
    Args: 
        :infile:            input IMAGE file (on extension that cv2 can read ), e.g. "image.png"
        :output_filename:   intermediate GDS file
        :pixelx:            pixel size in x, in um 
        :pixely:            pixel size in y, in um 
        :cellname:          string with cellname, e.g. "TOP"
        :level:             int from 0 to 255 (0 = black, 255=white) , pixels gray value to be passed to GDS 
        :layer:             (optional) gray level, defaults to 0 
        :datatype:          (optional) gds datatype, defaults to 0 
        :verbose:           (optional) defaults to False, if True prints 

    """
    
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

    #this is a gds file with single square that will be instantiated 
    outfile="image.gds"
    writer = gdspy.GdsWriter(outfile,unit=1.0e-6,precision=1.0e-9)
    cell = lib.new_cell(cellname)
    
    i=0 
    j=0 
    
    rect = gdspy.Rectangle((pixelx*j,-pixely*i),(pixelx*(j+1), -pixely*(i+1)), layer, datatype)                   
    cell.add(rect)        
    writer.write_cell(cell)
    writer.close()
    
    print("Exported the image file "+str(infile) + " into " + str(outfile))
    
    #####--------------------------------
    print("Starting making instances")
    
    import pya
    
    layout = pya.Layout()
    
    #input_filename = "fresnel_phase_mask_newlayers.gds"
    cell_name = "TOP" #name of the cell for instance array cannot be the same as the input gds filename

    #create cell at top 
    top = layout.create_cell(cell_name)

    #gds files to read (could also be a list)
    gds_files = [outfile]

    for file in gds_files:
        layout.read(file) #read the files

        for top_cll in layout.top_cells():
            if (top_cll.name != cell_name): #Don't insert in the top_cell
                cell_index = top_cll.cell_index()

                for hn, hi in enumerate(harray):
                    if hn == (len(harray)-1):
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
                                    new_instance = pya.CellInstArray( cell_index, pya.Trans(pya.Vector(int(j*pixelx*1000),int(-i*pixely*1000))), pya.Vector(0, 0), pya.Vector(0, 0), 0, 0)
                                    top.insert( new_instance ) #insert the cell in the array
            else:
                if cellname==cell_name:
                    print("Cell names seem to be the same. Proceed with caution, the final file might not be complete.")
                
    #write to gds
    layout.write(output_filename)
    print("Done")
