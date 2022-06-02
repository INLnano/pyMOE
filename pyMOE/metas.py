##metas.py 

def metasurface_pillars(xsiz,ysiz, pixelx, pixely, p, dvarx,cellname, outfilen, tolerance=0.01, verbose=False ): 
    """
    (void) Transform a phase profile in 2D for a pillar based metasurface 
    'xsiz' = size in x in um 
    'ysiz' = size in y in um
    'pixelx'   = pixel size in x in um
    'pixely'   = pixel size in y in um
    'p'        = periodicity in um 
    'dvarx'    = 2D matrix with the pillar diameters  -> Can be approximated to required diameters
    'cellname' = string with name of cell, e.g. 'TOP'
    'oufilen'  = string filename of output gds
    
    'tolerance' defaults to 0.01, can be decreased to have better defined circles
    'verbose' if True, prints during execution 
    """    
    import gdspy 
    import numpy as np 

    harray = np.arange(0, ysiz, p)
    warray = np.arange(0, xsiz, p)

    lib = gdspy.GdsLibrary()
    gdspy.current_library = gdspy.GdsLibrary() 

    writer = gdspy.GdsWriter(outfilen,unit=1.0e-6,precision=1.0e-9)
    cell = lib.new_cell(cellname)

    for hn, hi in enumerate(harray):
        if verbose == True: 
            print(hn/len(harray))
            
        if hn == (len(harray)-1):
            break

        for hw, wi in enumerate(warray):
            if hw == (len(warray)-1): 
                break   

            diameter = dvarx[hw,hn]
            #print(diameter)
            
            #avoid pillars smaller than 50 nm diameter 
            if diameter < 0.05: 
                diameter = 0.05 

            #Attention, we are using default tolerance, the tolerance could be decreased 
            circle = gdspy.Round((pixelx*hn, -pixely*hw), diameter/2,tolerance = tolerance, max_points=1e6)
            cell.add(circle)

        writer.write_cell(cell)

    writer.close()

    print("Created the metasurface mask in the file "+str(outfilen))