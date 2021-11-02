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
    (void) plots the gds for inspection in python matplotlib 
    Note: if plotting different gds together, please make sure they are aligned (e.g. centered at origin) 
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
    
    #fig = plt.figure(figsize=(6,5))
    #subplot = fig.add_subplot(111)
    #subplot.title.set_text(filename)

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
def inspect_gds2layers(filename, **kwargs): 
    from gdshelpers.geometry.chip import Cell
    import gdspy
    from shapely.geometry import MultiPolygon, Polygon 
    from gdshelpers.geometry.chip import Cell
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
    print(layers)
    datatps = list(main_cell.get_datatypes())
    print(datatps)
    
    xmaxs =[]
    ymaxs =[]
    xmins =[]
    ymins =[]
 
    for i in layers: 
        for j in datatps: #here is fine because datatps is single element 
                        
            if np.max(layers)>=128: 
                gss_norm = i/np.max(layers)
                print(gss_norm)
            else: #the levels are within 0-127 range
                gss_norm = i/127
                print(gss_norm)
            
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
                print("There is an invalid layer!")
            
            
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
            
    return cell_multipol, pol_dict, xmn, xmx, ymn, ymx