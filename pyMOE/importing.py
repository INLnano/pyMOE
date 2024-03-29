"""
importing.py 
Module containing functions to import gds files, inspect them in matplotlib and save gds files as image files 


"""

import numpy as np 
import matplotlib.pyplot as plt 
import gdspy 
from gdspy import FlexPath
from shapely.geometry import MultiPolygon, Polygon
import pickle 
import cv2 

# from gdshelpers.geometry.chip import Cell

 
 
def makesubplot(x,y, *argv, **kwargs):
    """
    Makes a subplot object for plotting 
    """
    #plot in subplot
    axes = kwargs.pop("axes", None)
    if not axes:
        fig, axes = plt.subplot()

    return axes.plot(x,y, *argv, **kwargs)


def inspect_gds2(filename, colors, rescale=0, **kwargs): 
    """
    Plots the gds for inspection in python matplotlib, only with the borders of the features 
    Note: if plotting different gds file together, please make sure they are aligned (e.g. centered at origin) 
    
    Args:
        :filename:   string gds filename (e.g. 'yolo.gds')
        :color:      string color name for plotting  
        :rescale:    int rescaling factor for all points in mask, by default is 0 
    
    """

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
    

def inspect_gds2layers(filename, norm, rescale=0,verbose = False, **kwargs ): 
    """
    Returns cell, polygon dictionary, (xmin,xmax) and (ymin,ymax) for representation of the gds layers into a grayscale image  
    
    Args:
        :filename:  string gds filename (e.g. 'yolo.gds')
        :norm:      maximum level of gray (e.g. 128)  
        :verbose:   if True show some information while doing the operations. defaults false
        
        for **kwargs, add 'axes = subplot', where subplot has been previously defined as subplot = fig.add_subplot(111) 
        with 'import matplotlib.pyplot as plt' and 'fig = plt.figure()'
    """

    
    #create cell to store the layers 
    cell_multipol = gdspy.Cell('top')

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
                if rescale!=0:
                    x1= np.array(x1)*rescale #get coords in um 
                    y1= np.array(y1)*rescale #get coords in um 
                #print(x1)
                axess.fill(x1, y1, fc=colorx)
                
                xmaxs = np.append(np.max(x1),xmaxs)
                ymaxs = np.append(np.max(y1),ymaxs)
                xmins = np.append(np.min(x1),xmins)
                ymins = np.append(np.min(y1),ymins)
                
                
            pol_dict[(i,j)] = MultiPolygon(pol_dict[(i,j)])
            polslayer = gdspy.PolygonSet(pol_dict[(i,j)], layer=i)
            #if polslayer.is_valid ==False: 
            #    if verbose == True: 
            #        print(i/(norm-1))
            #        print("There is an invalid polygon at "+str(i)+", but I will proceed, please check the result at the end")
            
            
            cell_multipol.add(polslayer)
    
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
    

def inspect_gds2layersplt(filename, norm, rescale=0, verbose = False, **kwargs ):  
    """
    Returns cell, polygon dictionary, (xmin,xmax) and (ymin,ymax) for representation of ALL gds layers into a grayscale image  
    
    Args:
        :filename:  string gds filename (e.g. 'yolo.gds')
        :norm:      maximum level of gray (e.g. 128)  
        :verbose:   if True show some information while doing the operations. defaults false
        
        for **kwargs, add 'axes = subplot', where subplot has been previously defined as subplot = fig.add_subplot(111) 
        with 'import matplotlib.pyplot as plt' and 'fig = plt.figure()'
    """
    
    #create cell to store the layers 
    gdspy.current_library = gdspy.GdsLibrary()
    cell_multipol = gdspy.Cell('top')

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
                if rescale!=0:
                    x1= np.array(x1)*rescale #get coords in um 
                    y1= np.array(y1)*rescale #get coords in um 
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
            polslayer = gdspy.PolygonSet(pol_dict[(i,j)], layer=i)
            #if polslayer.is_valid ==False: 
            #    if verbose == True: 
            #        print(i/(norm-1))
            #        print("There is an invalid polygon at "+str(i)+", but I will proceed, please check the result at the end")
            
            cell_multipol.add(polslayer)
    
      
    xmx = np.max(xmaxs)
    ymx = np.max(ymaxs)
    xmn = np.min(xmins)
    ymn = np.min(ymins)  
    axess.xlim([xmn,xmx])
    axess.ylim([ymn,ymx])     
    #axess.savefig("special.tiff", bbox_inches=0, pad_inches = 0) 
    
    print("xmax is "+ str(xmx))
    print("ymax is "+ str(ymx))
    print("xmin is "+ str(xmn))
    print("ymin is "+ str(ymn)) 
            
    return cell_multipol, pol_dict, xmn, xmx, ymn, ymx
    


def gds2img(infile,outfile,norm, rescaled=0, verbose=False): 
    """
    (void) plots the gds for inspection in python matplotlib and saves the figure to image file  
    Note: if plotting different gds file together, please make sure they are aligned (e.g. centered at origin) 
    
    Args:
        :infile:    string gds filename (e.g. 'yolo.gds')
        :outfile:   string img filename (e.g. 'yolo.tiff')
        :norm:      maximum level of gray (e.g. 128)  
        :rescaled:  if !=0 will rescale the layout 
        :verb:      if True, shows verbose, defaults to False 
    """
    
    fig = plt.figure()
        
    cell_multipol, pol_dict, xmn, xmx, ymn, ymx = inspect_gds2layersplt(infile,norm,rescale=rescaled, verbose = verbose, axes=plt)
    
    ax = fig.add_subplot(1,1,1)
    
    fig.set_dpi(1)

    plt.axis('equal')
    plt.axis('off')

    plt.tight_layout(pad=0, w_pad=0, h_pad=0) 
    fig.tight_layout(w_pad=0, h_pad=0, pad =0)
    
    plt.gca().set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0) 
    plt.margins(0,0)
    
    fig.set_size_inches(((xmx-xmn)), ((ymx-ymn)))

    fig.savefig("temp.png", bbox_inches=0,pad_inches = 0, dpi=1)

    #remove any white padding
    im = cv2.imread("temp.png")

    cv2.imwrite(outfile, im)
    
    print("Imported file "+infile+" and exported into "+outfile+ " with size "+ str(int(xmx)) + " x " + str(int(ymx)) + " pixels.")
    
