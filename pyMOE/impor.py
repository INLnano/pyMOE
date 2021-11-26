####imports gds and transforms into gray image 

#plot in subplot
def makesubplot(x,y, *argv, **kwargs):
    import matplotlib.pyplot as plt 
    
    axes = kwargs.pop("axes", None)
    if not axes:
        fig, axes = plt.subplot()

    return axes.plot(x,y, *argv, **kwargs)

##### INSPECT GDS IN MATPLOTLIB 
def inspect_gds2(filename, colors, rescale=0, **kwargs): 
    """
    (void) plots the gds for inspection in python matplotlib, only with the borders of the features 
    Note: if plotting different gds file together, please make sure they are aligned (e.g. centered at origin) 
    'filename' = string gds filename (e.g. 'yolo.gds')
    'color'    = string color name for plotting  
    rescale    = int rescaling factor for all points in mask, by default is 0 
    
    """
    import gdspy 
    from gdspy import FlexPath
    import matplotlib.pyplot as plt 
    from shapely.geometry import Polygon
    import pickle 
    import numpy as np
    
    lib = gdspy.GdsLibrary(infile=filename)
    main_cell = lib.top_level()[0]
    
    axess = kwargs.pop("axes", None)

    #control vars 
    n=0
    ct=1

    #while it is empty read again 
    while main_cell.polygons ==[]:
        lib = gdspy.GdsLibrary(infile=filename)
        main_cell = lib.top_level()[0]
        n=n+1 #just a control to avoid infinite loop 
        if n==10: #if we got to 10 trials 
            print("Cannot read polygons in this GDS file after"+str(n)+" trials.")
            break 
    else: 
        if main_cell.polygons !=[]: 
            print(str(np.size(main_cell.polygons))+" polygons found...")
            pol_dict = main_cell.get_polygons(by_spec=False)
            #we first get all the polygons
            #print(np.size(pol_dict[:]))
            for po, pod in enumerate(pol_dict[:]):
                polygon1 = Polygon(pod)
                x1,y1 = polygon1.exterior.xy
                if rescale!=0:
                    x1= np.array(x1)*rescale #get coords in um 
                    y1= np.array(y1)*rescale #get coords in um 
                axess.fill(x1, y1, alpha=0.5, fc=colors, ec='none')
                makesubplot(x1,y1, color=colors, axes=axess)

        #now we get all the paths 
        main_cell = lib.top_level()[0]
        paths = main_cell.get_paths()
        
        layers = list(main_cell.get_layers())
        print(layers)
        datatps = list(main_cell.get_datatypes())
                
        if main_cell.get_paths() !=[]:
            print(str(np.size(main_cell.get_paths()))+" paths found...")
            for po, pod in enumerate(paths):
                #print(po)
                #to be sure we get all the polygons lets make the paths thicker 
                polypath = FlexPath(paths[po].points, width =0.1, offset=0, corners='natural', ends='flush', bend_radius=None, tolerance=0.01, precision=0.001, max_points=1000000, gdsii_path=False, width_transform=True, layer=2, datatype=0)
                polygon = polypath.get_polygons(by_spec=False)
                
                #print(polygon)
                if polygon!=[]:
                    if ct==1:
                        print("And paths are being converted to polygons...")
                        ct=0
                    polygon1 = Polygon(polygon[0])
                    #x1,y1 = polygon1.interior.xy
                    x1,y1 = polygon1.exterior.xy
                    if rescale!=0:
                        x1= np.array(x1)*rescale #get coords in um 
                        y1= np.array(y1)*rescale #get coords in um  
                    
                    #axess.fill(x1, y1, alpha=0.5, fc=colors, ec='none')
                    makesubplot(x1,y1, color=colors, axes=axess)
                    
    axess.set_aspect('equal')        
    print("Done.")
    

##########FUNCTION TO INSPECT ALL LAYERS OF GDS (not only shapes)
def inspect_gds2layers(filename, norm,verbose = False, **kwargs ): 
    """
    Returns cell, polygon dictionary, (xmin,xmax) and (ymin,ymax) for representation of the gds layers into a grayscale image  
    'filename' = string gds filename (e.g. 'yolo.gds')
    'norm' = maximum level of gray (e.g. 128)  
    'verbose' = show some information while doing the operations 
    for **kwargs, add 'axes = subplot', where subplot has been previously defined as subplot = fig.add_subplot(111) 
    with 'import matplotlib.pyplot as plt' and 'fig = plt.figure()'
    """
    from gdshelpers.geometry.chip import Cell
    import gdspy
    from shapely.geometry import MultiPolygon, Polygon 
    import numpy as np 
    
    #create cell to store the layers 
    cell_multipol = Cell('top')

    lib = gdspy.GdsLibrary(infile=filename)
    main_cell = lib.top_level()[0]
    
    axess = kwargs.pop("axes", None)
        
    #control var
    n=0

    #while it is empty read again 
    while main_cell.polygons ==[]:
        try:
            lib = gdspy.read_gds(infile=filename)
        except: 
            continue 
            
        main_cell = lib.top_level()[0]
        n=n+1 #just a control to avoid infinite loop 
        if n==10: #if we got to 10 trials 
            print("Cannot read polygons in this GDS file after "+str(n)+" trials.")
            break 
    else: 
        if main_cell.polygons !=[]: 
            print(str(np.size(main_cell.polygons))+" polygons found...")
            pol_dict = main_cell.get_polygons(by_spec=True)

    layers = list(main_cell.get_layers())
    datatps = list(main_cell.get_datatypes())
    
    if verbose == True: 
        print("Layers: "+ str(layers))
        print("Datatypes: "+ str(datatps))
    
    xmaxs =[]
    ymaxs =[]
    xmins =[]
    ymins =[]
 
    for i in layers: 
        for j in datatps: #here is fine because datatps is single element         
            if np.max(layers)>=norm: 
                gss_norm = i/np.max(layers)
                #print(gss_norm)
            else: #the levels are within 0-127 range, the gray level is norm to the available levels 
                gss_norm = i/(norm-1)
                #print(gss_norm)
            
            ### lower levels show more black 
            #(0,0,0) is black
            #(1,1,1) is white
            #we want lower levels to show more black
            
            colorx=(gss_norm,gss_norm,gss_norm)
            
            #print(i)
            for k, pols in enumerate(pol_dict[(i,j)][:]):
                pol_dict[(i,j)][k] = Polygon(pol_dict[(i,j)][k])
                x1,y1 = pol_dict[(i,j)][k].exterior.xy 
                #print(x1)
                axess.fill(x1, y1, fc=colorx)
                
                xmaxs = np.append(np.max(x1),xmaxs)
                ymaxs = np.append(np.max(y1),ymaxs)
                xmins = np.append(np.min(x1),xmins)
                ymins = np.append(np.min(y1),ymins)
                
                
            pol_dict[(i,j)] = MultiPolygon(pol_dict[(i,j)])
            polslayer = pol_dict[(i,j)]
            if polslayer.is_valid ==False: 
                if verbose == True: 
                    print(i/(norm-1))
                    print("There is an invalid polygon at "+str(i)+", but I will proceed, please check the result at the end")
            
            cell_multipol.add_to_layer(i, polslayer)

    xmx = np.max(xmaxs)
    ymx = np.max(ymaxs)
    xmn = np.min(xmins)
    ymn = np.min(ymins)
    
    print("xmax is "+ str(xmx))
    print("ymax is "+ str(ymx))
    print("xmin is "+ str(xmn))
    print("ymin is "+ str(ymn)) 
    
    axess.set_xlim([xmn, xmx])
    axess.set_ylim([ymn, ymx])
    #axess.margins('tight')
    
    
    return cell_multipol, pol_dict, xmn, xmx, ymn, ymx
    
    
##########FUNCTION TO INSPECT ALL LAYERS OF GDS (not only shapes)
def inspect_gds2layersplt(filename, norm, verbose = False, **kwargs ):  
    """
    Returns cell, polygon dictionary, (xmin,xmax) and (ymin,ymax) for representation of the gds layers into a grayscale image  
    'filename' = string gds filename (e.g. 'yolo.gds')
    'norm' = maximum level of gray (e.g. 128)  
    'verbose' = show some information while doing the operations 
    for **kwargs, add 'axes = plt' where we have previously 
    'import matplotlib.pyplot as plt' and 'fig = plt.figure()'
    """
    from gdshelpers.geometry.chip import Cell
    import gdspy
    from shapely.geometry import MultiPolygon, Polygon 
    import numpy as np 
    import matplotlib.pyplot as plt 
    
    #create cell to store the layers 
    cell_multipol = Cell('top')

    lib = gdspy.GdsLibrary(infile=filename)
    main_cell = lib.top_level()[0]
    
    axess = kwargs.pop("axes", None)
    #print(axess)
    #print(*kwargs)
        
    #control var
    n=0

    #while it is empty read again 
    while main_cell.polygons ==[]:
        #lib = gdspy.GdsLibrary(infile=filename)
        try:
            lib = gdspy.read_gds(infile=filename)
        except: 
            continue 
            
        main_cell = lib.top_level()[0]
        n=n+1 #just a control to avoid infinite loop 
        if n==10: #if we got to 10 trials 
            print("Cannot read polygons in this GDS file after "+str(n)+" trials.")
            break 
    else: 
        if main_cell.polygons !=[]: 
            print(str(np.size(main_cell.polygons))+" polygons found...")
            pol_dict = main_cell.get_polygons(by_spec=True)

    layers = list(main_cell.get_layers())
    datatps = list(main_cell.get_datatypes())
    
    if verbose == True: 
        print("Layers: "+ str(layers))
        print("Datatypes: "+ str(datatps))
    
    xmaxs =[]
    ymaxs =[]
    xmins =[]
    ymins =[]
 
    for i in layers: 
        for j in datatps: #here is fine because datatps is single element         
            if np.max(layers)>=norm: 
                gss_norm = i/np.max(layers)
                #print(gss_norm)
            else: #the levels are within 0-127 range, the gray level is norm to the available levels 
                gss_norm = i/(norm-1)
                #print(gss_norm)
            
            ### lower levels show more black 
            #(0,0,0) is black
            #(1,1,1) is white
            #we want lower levels to show more black
            
            colorx=(gss_norm,gss_norm,gss_norm)
            
            #print(i)
            for k, pols in enumerate(pol_dict[(i,j)][:]):
                pol_dict[(i,j)][k] = Polygon(pol_dict[(i,j)][k])
                x1,y1 = pol_dict[(i,j)][k].exterior.xy 
                #print(x1)
                fplot = axess.fill(x1, y1, fc=colorx)
                axess.axis('off')
                #fplot.axes.get_xaxis().set_visible(False)
                #fplot.axes.get_yaxis().set_visible(False)
                
                xmaxs = np.append(np.max(x1),xmaxs)
                ymaxs = np.append(np.max(y1),ymaxs)
                xmins = np.append(np.min(x1),xmins)
                ymins = np.append(np.min(y1),ymins)
                
                
            pol_dict[(i,j)] = MultiPolygon(pol_dict[(i,j)])
            polslayer = pol_dict[(i,j)]
            if polslayer.is_valid ==False: 
                if verbose == True: 
                    print(i/(norm-1))
                    print("There is an invalid polygon at "+str(i)+", but I will proceed, please check the result at the end")
            
            
            cell_multipol.add_to_layer(i, polslayer)
    
      
    xmx = np.max(xmaxs)
    ymx = np.max(ymaxs)
    xmn = np.min(xmins)
    ymn = np.min(ymins)  
    axess.xlim([xmn,xmx])
    axess.ylim([ymn,ymx])     
    axess.savefig("special.tiff", bbox_inches=0, pad_inches = 0) 
    
    print("xmax is "+ str(xmx))
    print("ymax is "+ str(ymx))
    print("xmin is "+ str(xmn))
    print("ymin is "+ str(ymn)) 
    
    #axess.set_xlim([xmn, xmx])
    #axess.set_ylim([ymn, ymx])
    #axess.margins('tight')
            
    return cell_multipol, pol_dict, xmn, xmx, ymn, ymx
    


def gds2img(infile,outfile,norm, verbose=False): 
    """
    (void) plots the gds for inspection in python matplotlib and saves the figure to image file  
    Note: if plotting different gds file together, please make sure they are aligned (e.g. centered at origin) 
    'infile' = string gds filename (e.g. 'yolo.gds')
    'outfile' = string img filename (e.g. 'yolo.tiff')
    'norm' = maximum level of gray (e.g. 128)  
    'verb' if True, shows verbose, defaults to False 
    """
    import cv2 
    from matplotlib import pyplot as plt 
    import numpy as np
    
    fig = plt.figure()
    cell_multipol, pol_dict, xmn, xmx, ymn, ymx = inspect_gds2layersplt(infile,norm,verbose = verbose, axes=plt)
    
    #This dpi was chosen as high enough... but can be changed 
    my_dpi=4000

    plt.axis('equal')
    plt.axis('off')
    
    fig.set_size_inches(((xmx-xmn)+4)/my_dpi, ((ymx-ymn)+4)/my_dpi)
    fig.set_dpi(my_dpi)
    
    plt.xlim([0,xmx-xmn])
    plt.ylim([0,ymx-ymn])
    
    plt.tight_layout(pad=0, w_pad=0, h_pad=0) 
    plt.subplots_adjust(hspace = 0, wspace=0)
    fig.tight_layout(w_pad=0, h_pad=0, pad =0)
    fig.savefig(outfile,bbox_inches=0, pad_inches = 0)
    
    #remove any white padding
    im = cv2.imread(outfile)
    gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    gray = 255*(gray < 128).astype(np.uint8) 
    coords = cv2.findNonZero(gray) 
    x, y, w, h = cv2.boundingRect(coords) 
    nopad = im[y:y+h, x:x+w] 

    cv2.imwrite(outfile, nopad)
    
    print("Imported file "+infile+" and exported into "+outfile+ " with size "+ str(int(xmx)) + " x " + str(int(ymx)) + " pixels.")
    
