"""
metas.py 
Module containing functions to create metasurfaces from phase masks  


"""

import cv2
import numpy as np 
from pyMOE.utils import progress_bar, Timer
from pyMOE.gds_klops import rescale_layout, rotate_layout

import pya


def metasurface_from_phase(xsiz, ysiz, pixelx, pixely, p, aperture_vals, topcellname, outfilen, elements='pillar', \
                           verbose=False, rotation=None, scaling=None, grid='square', mindim = 0.05, smallerdim =0, \
                           largest_phase=None, dbu=1): 
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
        :elements:    ´    element (as klayout polygon) to be used as individual meta-element (also accepts array of such elements for iteration, with same dimension as unique values in aperture_vals). If == 'pillar' (default) -> circle with 1 um diameter 
        :verbose:          if True, prints during execution 
        :rotation:         array with the rotation angles of unique meta-elements (1:1 correspondence with unique aperture_vals values!), Rotation angle is anti-clockwise in degrees. If None (default), sets the rotation angle = 0 for all elements. 
        :scaling:          array with the scaling factor of unique meta-element (1:1 correspondence with unique aperture_vals values!), Scaling factor is with respect to the dimension of the individual meta-element. If None (default), sets scaling factor to 1.0 for all elements.  
        :grid:             Type of grid, options are 'square' or 'hex'. Default is 'square'. Please make sure the aperture_vals have been evaluated in an hexagonal grid, to make sure the values match, and that the pixely is given for an hexagonal lattice
        :mindim:           clipping scaling factor (cannot scale below a certain value, to avoid very small elements)
        :smallerdim:       lowest scaling factor
        :largest_phase:    largest phase in the phase mask. If None, takes the maximum of aperture_vals  
        :dbu:              scaling because of a particular dbu 
    
    Returns:
        None
    """   
    from pyMOE.utils import Timer, progress_bar
    
    #total number of elements count
    tot_meta = 0
    
    #some global variables 
    tolerance, nr_points, mindim, smallerdim = 0.001, 15, 0.05, 0

    
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
        xv, yv = np.meshgrid(np.arange(0, xsiz, pixelx, dtype=float), 
                             np.arange(0, ysiz, pixely, dtype=float))
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
        
    elif grid == 'hex': 
        x= np.arange(0, xsiz, pixelx, dtype=float) # arange is preferred over linspace because it keeps the pixel size! 
        y= np.arange(0, ysiz, pixely, dtype=float) # arange is preferred over linspace because it keeps the pixel size! 
        xv, yv = np.meshgrid(x,y)
        xv[::2, :] += pixelx/2
        #yv = yv*np.cos(np.radians(30))
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
    else: 
        print("Unsuported grid argument!")
    
    #################################
    ###elements options:
    pflag = 3 
    if type(elements) is not str:
        print("Custom metasurface")
        
        #print(np.asarray(elements).size)
                
        if np.asarray(elements).size>1:
            assert len(elements)==len(phase_array), "The length of unique phase values and elements argument array is different." 
            pflag = 2
        elif np.asarray(elements).size==1: 
            print("Single element provided.")
            pflag =0

    elif type(elements) is str: ###This corresponds to the default
        if elements=='pillar':  
            print("Pillar metasurface")
            diameter = 1 #standard 1 um 
            if scaling_flag: 
                print("By default all pillars have "+str(diameter)+" um diameter without scaling. Not sure this is was what is wanted.")

            elements = pya.DPolygon.ellipse(pya.DBox(-0.5, -0.5, 0.5, 0.5), nr_points)
            pflag = 0 
        
        else: 
            pflag = 1 
    else: 
        if elements is None: 
            pflag = 4
        else: 
            pflag = 1
        
    if pflag==1: 
        print("Unsuported elements argument!")  
    ######################################################################################################
    #####--------------------------------
    print("Building the metasurface...")
    print("Total of "+str(len(phase_array))+" layers.")
    
    with Timer():
        layout2 = pya.Layout()

        #create cell at top 
        cell = layout2.create_cell(topcellname)
        layer = layout2.layer(0,0)
        cell_index1 = layout2.cell(topcellname).cell_index()
        
        layout2.dbu = 0.001
   
        for ids, phase in enumerate(phase_array):
            first = 1 
            
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
                
            rect = pya.DPolygon()

            with Timer():
                if (scaling_array[ids] >0) and (phase<=largest_phase):# & (scaling_array[ids] < p): 
                    for hn, (hi, wi) in enumerate(zip(harray,warray)):  
                        #print(hi,wi)
                        if first==1: 
                            if verbose == True:                         
                                progress_bar(hn/len(harray))

                            #avoid features with scaling smaller than mindim, setting them to smallerdim(=0) 
                            if scaling_array[ids] < mindim: 
                                scaling_array[ids] = smallerdim

                            newpolygon, newpolygon2 = pya.DPolygon(), pya.DPolygon()

                            if pflag==2:
                                newpolygon = elements[ids]
                            elif pflag==0:
                                newpolygon = elements

                            t = pya.DCplxTrans(float(scaling_factor), float(angle), False, pya.DVector(float(hi/dbu),float(wi/dbu)))
                            
                            newpolygon2 = newpolygon.dup()
                            newpolygon2.transform(t)

                            cell.shapes(cell_index1).insert(newpolygon2)
                            
                        else: 
                            t = pya.DCplxTrans(float(scaling_factor), float(angle), False, pya.DVector(float(hi/dbu),float(wi/dbu)))

                            newpolygon2 = newpolygon.dup()
                            newpolygon2.transform(t)

                            cell.shapes(cell_index1).insert(newpolygon2)

                        tot_meta = tot_meta + 1
            
                        layout2.write(outfilen)
                    
                    progress_bar(1)
                    if verbose == True:    
                        print("So far "+str(tot_meta)+" elements and counting.") 
        
    print("\n Saved the metasurface mask with "+str(tot_meta)+" meta-elements in the file "+str(outfilen))
    
    
      

def metasurface_from_phase_instances (xsiz, ysiz, pixelx, pixely, p, aperture_vals, topcellname, outfilen, elements='pillar', \
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
        :elements:    ´     element (as klayout polygon) to be used as individual meta-element (also accepts array of such elements for iteration, with same dimension as unique values in aperture_vals). If == 'pillar' (default) -> circle with 1 um diameter  
        :infile:            string with filename to be used as meta-element
        :verbose:           if True, prints during execution 
        :rotation:          array with the rotation angles of unique meta-elements (1:1 correspondence with unique aperture_vals values!), Rotation angle is anti-clockwise in degrees. If None (default), sets the rotation angle = 0 for all elements. 
        :scaling:           array with the scaling factor of unique meta-element (1:1 correspondence with unique aperture_vals values!), Scaling factor is with respect to the dimension of the individual meta-element. If None (default), sets scaling factor to 1.0 for all elements.  
        :grid:             Type of grid, options are 'square' or 'hex'. Default is 'square'. Please make sure the aperture_vals have been evaluated in an hexagonal grid, to make sure the values match, and that the pixely is given for an hexagonal lattice
        :mindim:            clipping scaling factor (cannot scale below a certain value, to avoid very small elements)
        :smallerdim:        lowest scaling factor
        :tempfile:          string with name of a temporary file that will be used to have the individual elements and make the instances 
        :largest_phase:     largest phase in the phase mask. If None, takes the maximum of aperture_vals  
    
    Returns:
        None
    """   
    from pyMOE.utils import Timer, progress_bar
    
    #total number of elements count
    tot_meta = 0
    
    #some global variables 
    tolerance, nr_points, mindim, smallerdim = 0.001, 15, 0.05, 0

    
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
        #yv = yv*np.cos(np.radians(30))
        positions = np.array([xv.ravel(), yv.ravel()])
        positions_xv = positions[0]
        positions_yv = positions[1]
    else: 
        print("Unsuported grid argument!")
        
    #################################
    ###elements options:
    pflag = 3
    if infile is None: 
        if type(elements) is not str:
            print("Custom metasurface")
            
            #print(np.asarray(elements).size)
                    
            if np.asarray(elements).size>1:
                assert len(elements)==len(phase_array), "The length of unique phase values and elements argument array is different." 
                pflag = 2
            elif np.asarray(elements).size==1: 
                print("Single element provided.")
                pflag =0

        elif type(elements) is str: ###This corresponds to the default
            if elements=='pillar':  
                print("Pillar metasurface")
                diameter = 1 #standard 1 um 
                if scaling_flag: 
                    print("By default all pillars have "+str(diameter)+" um diameter without scaling. Not sure this is was what is wanted.")

                elements = pya.DPolygon.ellipse(pya.DBox(-0.5, -0.5, 0.5, 0.5), nr_points)
                
                pflag = 0 
            
            else: 
                pflag = 1 
        else: 
            if elements is None: 
                pflag = 4
            else: 
                pflag = 1
    
    if pflag==1: 
        print("Unsuported elements argument!")              

    ########################################################################################################
    #####--------------------------------
    print("Building the metasurface...")
    print("Total of "+str(len(phase_array))+" layers.")
    
    with Timer():
        layout2 = pya.Layout()

        #create cell at top 
        top = layout2.create_cell(topcellname)
   
        for ids, phase in enumerate(phase_array): 
            first = 1 
            
            angle  = rotation_array[ids]
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
                            layout3 = pya.Layout()

                            if first==1:
                                cell = layout3.create_cell(tempcellname)
                                layer = layout3.layer(0,0)
                                cell_index1 = layout3.cell(tempcellname).cell_index()

                                newpolygon, newpolygon2 = [], []

                                if pflag==2:
                                    newpolygon = elements[ids]
                                else:
                                    newpolygon = elements

                                newpolygon2 = newpolygon
                                
                                cell.shapes(cell_index1).insert(newpolygon2)

                                cell.write(tempfile)

                                first = 0 
                                
                                #define load opt for the next gds file with fst layer map (better to have...)
                                load_layout_options = pya.LoadLayoutOptions()

                                layout2.read(tempfile, load_layout_options)
                                layout2.layer(0,0)
                                cell_index1 = layout2.cell(tempcellname).cell_index()

                                
                            tot_meta = tot_meta + 1

                            new_instance1 = pya.DCellInstArray(cell_index1 , pya.DCplxTrans(float(scaling_factor),  \
                            float(angle), False, pya.DVector(float(hi),float(wi))))
                            top.insert( new_instance1 ) 

                        else: 
                            if first==1: 
                                layout2.read(infile)

                                rotate_layout(infile, tempcellname, angle, tempfile, transx =0, transy=0)
                                first = 0 

                                layout2.read(tempfile)
                                cell_index2 = layout2.cell(tempcellname).cell_index()

                            if cell_index2 is not None:    
                                new_instance1 = pya.DCellInstArray(cell_index2 , pya.DCplxTrans(float(scaling_factor), float(angle), False, pya.DVector(float(hi),float(wi))))
                                top.insert( new_instance1 ) 

                                tot_meta = tot_meta + 1
                    
                    layout2.write(outfilen)
            progress_bar(1)      
            print("So far "+str(tot_meta)+" elements and counting.")
            

    print("\n Saved the metasurface mask with "+str(tot_meta)+" meta-elements in the file "+str(outfilen))