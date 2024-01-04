"""
metas.py 
Module containing functions to create metasurfaces from phase masks  


"""

import cv2
import gdspy 
import numpy as np 
from pyMOE.utils import progress_bar, Timer
from pyMOE.gds_klops import rescale_layout, rotate_layout

import pya


def metasurface_from_phase(xsiz, ysiz, pixelx, pixely, p, aperture_vals, topcellname, outfilen, gdspyelements='pillar', \
                           verbose=False, rotation=None, scaling=None, grid='square', mindim = 0.05, smallerdim =0, \
                           largest_phase=None): 
    """
    Transform a 2D array (aperture_vals) representing the phase into a 2D metasurface and saves it to gds 
    
    Args: 
        :xsiz:             x size of aperture in x in um 
        :ysiz:             y size of aperture size in y in um
        :pixelx:           pixel size in x in um
        :pixely:           pixel size in y in um
        :p:                periodicity in um 
        :aperture_vals:    2D array with the phase  
        :topcellname:      string with name of top cell, e.g. 'TOP'
        :oufilen:          string filename of output gds
        :gdspyelements:    gdspy element to be used as individual meta-element (also accepts array of such elements for iteration, with same dimension as unique values in aperture_vals). If == 'pillar' (default) -> gdspy circle with 1 um diameter 
        :verbose:          if True, prints during execution 
        :rotation:         array with the rotation angles of unique meta-elements (1:1 correspondence with unique aperture_vals values!), Rotation angle is anti-clockwise in radians. If None (default), sets the rotation angle = 0 for all elements. 
        :scaling:          array with the scaling factor of unique meta-element (1:1 correspondence with unique aperture_vals values!), Scaling factor is with respect to the dimension of the individual meta-element. If None (default), sets scaling factor to 1.0 for all elements.  
        :grid:             Type of grid, options are 'square' or 'hex'. Default is 'square'. Please make sure the aperture_vals have been evaluated in an hexagonal grid, to make sure the values match. 
        :mindim:           clipping scaling factor (cannot scale below a certain value, to avoid very small elements)
        :smallerdim:       lowest scaling factor
        :largest_phase:    largest phase in the phase mask. If None, takes the maximum of aperture_vals  
    
    Returns:
        None
    """  
    from gdspy import Polygon, PolygonSet 
    from pyMOE.utils import Timer, progress_bar
    
    #total number of elements count
    tot_meta = 0
    
    #some global variables 
    tolerance, nr_points, mindim, smallerdim = 0.001, 15, 0.05, 0

    #Start the metasurface library  
    lib = gdspy.GdsLibrary()

    if largest_phase is None: 
        largest_phase = np.max(aperture_vals)
    
    #Extract unique values of phase from the aperture 2D array 
    phase_array = np.unique(aperture_vals[aperture_vals <=largest_phase])
    
    ##############################################################
    ###Various options for the metasurface currently given as args  
    if rotation is not None: 
        flag=1
        if isinstance(rotation, (int, float)): 
            rotation_array, flag = np.ones(len(phase_array)) * rotation, 0
        else: 
            assert len(rotation)==len(phase_array), "The length of unique phase values and rotation array is different." 
            rotation_array, flag = rotation , 0
        if flag: 
            print("Unsuported rotation argument!")
    else: 
        rotation_array = np.zeros(len(phase_array))  #default rotate by 0 degs 
    
    #################################
    if scaling is not None: 
        scaling_flag=0
        flag=1
        if isinstance(scaling, (int, float)): 
            scaling_array, flag = np.ones(len(phase_array)) * scaling, 0
        else: 
            assert len(scaling)==len(phase_array), "The length of unique phase values and rscaling array is different." 
            scaling_array, flag = scaling, 0
        if flag: 
            print("Unsuported scaling argument!")
    else: 
        scaling_array = np.ones(len(phase_array))  #default scale by 1 
        scaling_flag = 1 
    
    #################################
    if grid == 'square':
        xv, yv = np.meshgrid(np.arange(0, xsiz, pixelx, dtype=float), np.arange(0, ysiz, pixely, dtype=float))
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
        
    elif grid == 'hex': 
        x= np.arange(0, xsiz, pixelx, dtype=float) # arange is preferred over linspace because it keeps the pixel size! 
        y= np.arange(0, ysiz, pixely, dtype=float) # arange is preferred over linspace because it keeps the pixel size! 
        xv, yv = np.meshgrid(x,y)
        xv[::2, :] += pixelx/2
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
    else: 
        print("Unsuported grid argument!")
    
    #################################
    ###elements options:
    pflag = 3 
    if type(gdspyelements) is not str:
        print("Custom metasurface")
        
        #print(np.asarray(gdspyelements).size)
                
        if np.asarray(gdspyelements).size>1:
            assert len(gdspyelements)==len(phase_array), "The length of unique phase values and gdspyelements argument array is different." 
            pflag = 2
        elif np.asarray(gdspyelements).size==1: 
            print("Single gdspyelement element")
            pflag =0

    elif type(gdspyelements) is str: ###This corresponds to the default
        if gdspyelements=='pillar':  
            print("Pillar metasurface")
            diameter = 1 #standard 1 um 
            if scaling_flag: 
                print("By default all pillars have "+str(diameter)+" um diameter without scaling. Not sure this is was what is wanted.")

            gdspyelements = gdspy.Round((0, 0), diameter/2, tolerance = tolerance, number_of_points = nr_points, max_points=100)
            pflag = 0 
        
        else: 
            pflag = 1 
    else: 
        if gdspyelements is None: 
            pflag = 4
        else: 
            pflag = 1
        
    if pflag==1: 
        print("Unsuported gdspyelements argument!")  
    ######################################################################################################
    #####--------------------------------
    print("Building the metasurface...")
    print("Total of "+str(len(phase_array))+" layers.")
    
    with Timer():
        gdspy.current_library = gdspy.GdsLibrary() 
        writer = gdspy.GdsWriter(outfilen,unit=1.0e-6,precision=1.0e-9)
        cell = lib.new_cell(topcellname)
   
        for ids, phase in enumerate(phase_array):          
            angle  = rotation_array[ids]
            scaling_factor = scaling_array[ids]
            
            #tempcellname = "p"+str(np.round(phase,3))+"_s"+str(np.round(scaling_factor[ids],3))+"_r"+str(np.round(angle,3))

            print("Building meta-elements in layer "+str(ids)+":")
            
            if grid == 'square': 
                selection_ids = np.where(aperture_vals == phase)
                harray = selection_ids[1]*pixelx
                warray = selection_ids[0]*pixely
                
            if grid == 'hex':
                ###select the positions of each phase value in the aperture
                selection_ids = np.where(aperture_vals.ravel()==phase)
                harray = positions_xv[selection_ids]
                warray = positions_yv[selection_ids]

            with Timer():
                if (scaling_array[ids] >0) and (phase<=largest_phase):# & (scaling_array[ids] < p): 
                    for hn, (hi, wi) in enumerate(zip(harray,warray)):  
                        if verbose == True:                         
                            progress_bar(hn/len(harray))
                            
                        #avoid features with scaling smaller than mindim, setting them to smallerdim(=0) 
                        if scaling_array[ids] < mindim: 
                            scaling_array[ids] = smallerdim
 
                        newpolygon, newpolygon2 = [], []
                        cell.remove_polygons(lambda pts, layer, datatype: layer == 0)

                        if pflag==2:
                            newpolygon = gdspyelements[ids]
                        elif pflag==0:
                            newpolygon = gdspyelements

                        newpolygon2 = gdspy.copy(newpolygon)
                        newpolygon2 = newpolygon2.scale(scaling_factor)
                        newpolygon2 = newpolygon2.rotate(angle)
                        newpolygon2 = newpolygon2.translate(hi, wi) 
                        cell.add(newpolygon2)

                        tot_meta = tot_meta + 1
                        writer.write_cell(cell)
                    
                    progress_bar(1)
                    if verbose == True:    
                        print("So far "+str(tot_meta)+" elements and counting.")
                    
        writer.close() 
        
    print("\n Saved the metasurface mask with "+str(tot_meta)+" meta-elements in the file "+str(outfilen))
    
    

def metasurface_from_phase_instances (xsiz, ysiz, pixelx, pixely, p, aperture_vals, topcellname, outfilen, gdspyelements='pillar', \
                                      infile=None, verbose=False, rotation=None, scaling=None, grid='square',\
                                      mindim = 0.05, smallerdim =0, tempfile="temp.gds", largest_phase=None): 
    """
    Transform a 2D array (aperture_vals) representing the phase into a 2D metasurface using instances (from pya) package and saves it to gds 
    
    Args: 
        :xsiz:              x size of aperture in x in um 
        :ysiz:              y size of aperture size in y in um
        :pixelx:            pixel size in x in um
        :pixely:            pixel size in y in um
        :p:                 periodicity in um 
        :aperture_vals:     2D array with the phase  
        :topcellname:       string with name of top cell, e.g. 'TOP'
        :oufilen:           string filename of output gds
        :gdspyelements:     gdspy element to be used as individual meta-element (also accepts array of such elements for iteration, with same dimension as unique values in aperture_vals). If == 'pillar' (default) -> gdspy circle with 1 um diameter. If 'infile' is provided, ignores this design. 
        :infile:            string with filename to be used as meta-element
        :verbose:           if True, prints during execution 
        :rotation:          array with the rotation angles of unique meta-elements (1:1 correspondence with unique aperture_vals values!), Rotation angle is anti-clockwise in radians. If None (default), sets the rotation angle = 0 for all elements. 
        :scaling:           array with the scaling factor of unique meta-element (1:1 correspondence with unique aperture_vals values!), Scaling factor is with respect to the dimension of the individual meta-element. If None (default), sets scaling factor to 1.0 for all elements.  
        :grid:              Type of grid, options are 'square' or 'hex'. Default is 'square'. Please make sure the aperture_vals have been evaluated in an hexagonal grid, to make sure the values match. 
        :mindim:            clipping scaling factor (cannot scale below a certain value, to avoid very small elements)
        :smallerdim:        lowest scaling factor
        :tempfile:          string with name of a temporary file that will be used to have the individual elements and make the instances 
        :largest_phase:     largest phase in the phase mask. If None, takes the maximum of aperture_vals  
    
    Returns:
        None
    """   
    from gdspy import Polygon, PolygonSet 
    from pyMOE.utils import Timer, progress_bar
    
    #total number of elements count
    tot_meta = 0
    
    #some global variables 
    tolerance, nr_points, mindim, smallerdim = 0.001, 15, 0.05, 0

    #Start the metasurface library  
    lib = gdspy.GdsLibrary()
    
    if largest_phase is None: 
        largest_phase = np.max(aperture_vals)
        
    #Extract unique values of phase from the aperturea 2D array 
    phase_array = np.unique(aperture_vals[aperture_vals <=largest_phase])
    
    ##############################################################
    ###Various options for the metasurface currently given as args  
    if rotation is not None: 
        flag=1
        if isinstance(rotation, (int, float)): 
            rotation_array, flag = np.ones(len(phase_array)) * rotation, 0
        else: 
            assert len(rotation)==len(phase_array), "The length of unique phase values and rotation array is different." 
            rotation_array, flag = rotation , 0
        if flag: 
            print("Unsuported rotation argument!")
    else: 
        rotation_array = np.zeros(len(phase_array))  #default rotate by 0 degs 
    
    #################################
    if scaling is not None: 
        scaling_flag=0
        flag=1
        if isinstance(scaling, (int, float)): 
            scaling_array, flag = np.ones(len(phase_array)) * scaling, 0
        else: 
            assert len(scaling)==len(phase_array), "The length of unique phase values and rscaling array is different." 
            scaling_array, flag = scaling, 0
        if flag: 
            print("Unsuported scaling argument!")
    else: 
        scaling_array = np.ones(len(phase_array))    #default scale by 1 
        scaling_flag = 1 
    
    #################################
    if grid == 'square':
        xv, yv = np.meshgrid(np.arange(0, xsiz, pixelx, dtype=float), np.arange(0, ysiz, pixely, dtype=float))
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
        
    elif grid == 'hex': 
        x= np.arange(0, xsiz, pixelx, dtype=float) # arange is preferred over linspace because it keeps the pixel size! 
        y= np.arange(0, ysiz, pixely, dtype=float) # arange is preferred over linspace because it keeps the pixel size! 
        xv, yv = np.meshgrid(x,y)
        xv[::2, :] += pixelx/2
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
    else: 
        print("Unsuported grid argument!")
        
    #################################
    ###elements options:
    pflag = 3
    if infile is None: 
        if type(gdspyelements) is not str:
            print("Custom metasurface")
            
            #print(np.asarray(gdspyelements).size)
                    
            if np.asarray(gdspyelements).size>1:
                assert len(gdspyelements)==len(phase_array), "The length of unique phase values and gdspyelements argument array is different." 
                pflag = 2
            elif np.asarray(gdspyelements).size==1: 
                print("Single gdspyelement element")
                pflag =0

        elif type(gdspyelements) is str: ###This corresponds to the default
            if gdspyelements=='pillar':  
                print("Pillar metasurface")
                diameter = 1 #standard 1 um 
                if scaling_flag: 
                    print("By default all pillars have "+str(diameter)+" um diameter without scaling. Not sure this is was what is wanted.")

                gdspyelements = gdspy.Round((0, 0), diameter/2, tolerance = tolerance, number_of_points = nr_points, max_points=100)
                pflag = 0 
            
            else: 
                pflag = 1 
        else: 
            if gdspyelements is None: 
                pflag = 4
            else: 
                pflag = 1
    
    if pflag==1: 
        print("Unsuported gdspyelements argument!")              

    ########################################################################################################
    #####--------------------------------
    print("Building the metasurface...")
    print("Total of "+str(len(phase_array))+" layers.")
    
    with Timer():
        layout = pya.Layout()

        #create cell at top 
        top = layout.create_cell(topcellname)
   
        for ids, phase in enumerate(phase_array): 
            first = 1 
            cell_index2 = None 
            
            angle  = np.degrees(rotation_array[ids])
            scaling_factor = scaling_array[ids]
            fvalue = phase 
            
            tempcellname = "layer_"+str(ids)+"p"+str(np.round(fvalue,3))+"_s"+str(np.round(scaling_factor,3))+"_r"+str(np.round(angle,3))

            print("Building meta-elements in layer "+str(ids)+":")
               
            if grid == 'square': 
                selection_ids = np.where(aperture_vals == phase)
                harray = selection_ids[1]*pixelx
                warray = selection_ids[0]*pixely
                
            if grid == 'hex':
                ###select the positions of each phase value in the aperture
                selection_ids = np.where(aperture_vals.ravel()==phase)
                harray = positions_xv[selection_ids]
                warray = positions_yv[selection_ids]

            with Timer(): 
                if (scaling_array[ids] >0) and (phase<=largest_phase):# & (scaling_array[ids] < p):
                    for hn, (hi, wi) in enumerate(zip(harray,warray)):  
                        if verbose == True:                         
                            progress_bar(hn/len(harray))
                            
                        #avoid features with scaling smaller than mindim, setting them to smallerdim(=0) 
                        if scaling_array[ids] < mindim: 
                            scaling_array[ids] = smallerdim
 
                        if infile is None: 
                            lib = gdspy.GdsLibrary()
                            gdspy.current_library = gdspy.GdsLibrary()                                    

                            if first==1:
                                writer = gdspy.GdsWriter(tempfile, unit=1.0e-6, precision=1.0e-9) #the precision could be passed as argument if needed 

                                cell = lib.new_cell(tempcellname) 
                                newpolygon, newpolygon2 = [], []
                                cell.remove_polygons(lambda pts, layer, datatype: layer == 0)

                                if pflag==2:
                                    newpolygon = gdspyelements[ids]
                                else:
                                    newpolygon = gdspyelements

                                newpolygon2 = gdspy.copy(newpolygon)
                                cell.add(newpolygon2)

                                writer.write_cell(cell)
                                writer.close()
                                first = 0 

                                layout.read(tempfile)
                                cell_index1 = layout.cell(tempcellname).cell_index()

                                
                            tot_meta = tot_meta + 1

                            new_instance1 = pya.DCellInstArray(cell_index1 , pya.DCplxTrans(scaling_factor, angle, False, pya.DVector(float(hi),float(wi))))
                            top.insert( new_instance1 ) 

                        else: 
                            if first==1: 
                                layout.read(infile)

                                rotate_layout(infile, tempcellname, angle, tempfile, transx =0, transy=0)
                                first = 0 

                                layout.read(tempfile)
                                cell_index2 = layout.cell(tempcellname).cell_index()

                            if cell_index2 is not None:    
                                new_instance1 = pya.DCellInstArray(cell_index2 , pya.DCplxTrans(scaling_factor, angle, False, pya.DVector(float(hi),float(wi))))
                                top.insert( new_instance1 ) 

                                tot_meta = tot_meta + 1
                    
                    layout.write(outfilen)
            progress_bar(1)      
            print("So far "+str(tot_meta)+" elements and counting.")
            

    print("\n Saved the metasurface mask with "+str(tot_meta)+" meta-elements in the file "+str(outfilen))
