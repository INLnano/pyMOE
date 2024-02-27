
"""
gds_klops.py 
Module containing several functions for operations with/within mask files (e.g. gds) 

These functions might require pya, gdspy, gdstk or gdshelpers libraries. 
Functions using gdspy and/or pya are exemplified within the notebooks.  

####DISCLAIMER: these functions were written with klayout-0.26.11 version
####There might be some inconsistencies for later versions of pya library
####Nevertheless, tests were made for klayout-0.28.17 and no apparent conflict was found. 

"""


import pya 


####### #MERGE FUNCTION 
def merge_layer(readfile, cellname, layer_nr, datatype_nr, outputfile): 
    """
    (void) Merges all shapes within gds (only tested gds with single layer)
    
    Args:
        :readfile:     string filename of input gds
        :cellname:     string name of cell 
        :layer_nr:     int    layer number 
        :datatype_nr:  int    datatype number
        :outputfile:   string filename of output gds
    """
    layoutor = pya.Layout()

    lmap = layoutor.read(readfile)
    cell = layoutor.cell(cellname)
    cell.flatten(1)
    layer = layoutor.layer(layer_nr,datatype_nr)

    region = pya.Region(cell.shapes(layer))
    region.merge()
    
    cell.layout().clear_layer(layer)
    cell.shapes(layer).insert(region)
    
    layoutor.write(outputfile)
    
    print("Merged layers in " + outputfile)
    
########IMPORT FUNCTION 
def import_gds(fstgds_filename, fst_cellname, fst_layer_nr, fst_datatype_nr, \
               sndgds_filename, snd_cellname, snd_layer_nr, snd_datatype_nr, \
               output_filename, clear_gds = True ):
    """
    (void) Imports shapes from fst gds file INTO snd gds file 
    
    Args:
        :Ngds_filename:   string filename of N(=fst,snd) gds
        :N_cellname:      string name of cell in N(=fst,snd) gds 
        :N_layer_nr:      int layer number in N(=fst,snd) gds  
        :N_datatype_nr:   int datatype number in N(=fst,snd) gds  
        :output_filename: string filename of output gds
        :clear_gds:       clear gds, before inserting shapes, defaults to True 
    """
    
    #ly1 with 1st gds file
    ly1 = pya.Layout()
    lmap1 = ly1.read(fstgds_filename)
    cll1 = ly1.cell(fst_cellname)
    lyr1 = ly1.layer(fst_layer_nr,fst_datatype_nr)
    region1 = pya.Region(cll1.shapes(lyr1)) #define region1 as shapes from ly1-lyr1
    if clear_gds == True: 
        cll1.layout().clear_layer(lyr1) #clear ly1-lyr1

    #define load opt for the next gds file with fst layer map (better to have...)
    load_layout_options = pya.LoadLayoutOptions()
    load_layout_options.set_layer_map(lmap1,True)

    #ly2 with 2nd gds file
    ly2 = pya.Layout()
    lmap2 = ly2.read(sndgds_filename,load_layout_options)
    cll2 = ly2.cell(snd_cellname)
    lyr2 = ly2.layer(snd_layer_nr,snd_datatype_nr)
    region2 = pya.Region(cll2.shapes(lyr2)) #define region2 as shapes from ly2-lyr2
    if clear_gds == True: 
        cll2.layout().clear_layer(lyr2) #clear ly2-lyr2

    #create new layer in ly1 to receive shapes from ly2-lyr2 
    ly1.insert_layer(pya.LayerInfo(snd_layer_nr, snd_datatype_nr)) #create with same info
    lyr12 = ly1.layer(snd_layer_nr,snd_datatype_nr) #define the lyr2 in ly1 
    #USE ONLY ONE CELL, here the existing one in ly1
    cll1.shapes(lyr12).insert(region2) #cp shapes from ly2 region2 to ly1 lyr12
    cll1.shapes(lyr1).insert(region1)

    ly1.write(output_filename)
    
    print("Imported "+sndgds_filename+" layer " + str(snd_layer_nr)+ " into "+\
    fstgds_filename+" layer "+str(fst_layer_nr)+". Output file "+output_filename+" .")

######WRITE .DXF FILES 
def gds_to_dxf(inputfilename_gds, outputfilename_dxf):
    """
    (void) writes to file (it is named to be dxf... other extensions also work)
    """
    layout=pya.Layout()
    layout.read(inputfilename_gds)
    layout.write(outputfilename_dxf) 
    
def correct_gds(inputfilename_gds, outputfilename_dxf): 
    """
    (void) corrects gds to be able to read/write gds files 
    """
    gds_to_dxf(inputfilename_gds, outputfilename_dxf)
    gds_to_dxf(outputfilename_dxf, inputfilename_gds)
    
######CREATE INSTANCE ARRAY 
def instance_array(cell_name, input_filename, transx, transy, nr_inst_X, nr_inst_Y, pitx, pity, output_filename):
    """
    (void) Makes array of existing object(s) in the gds
    
    Args:
        :cell_name:        string name of cell
        :input_filename:   string name of the input gds 
        :transx:           translation vector in x in um 
        :transy:           translation vector in y in um
        :nr_inst_X:        int number of repeated objects in x
        :nr_inst_Y:        int number of repeated objects in y 
        :pitx:             int pitch in x in um
        :pity:             int pitch in y in um 
        :output_filename:  string filename of output gds
    """

    layout = pya.Layout()

    #pitches in x and y
    pitchx = pitx*1000 # pitx in um
    pitchy = pity*1000 # pit um

    #create cell at top 
    top = layout.create_cell(cell_name)

    #gds files to read (could also be a list)
    gds_files = [input_filename]
    cnt = 0

    for file in gds_files:
        layout.read(file) #read the files

        for top_cll in layout.top_cells():
            if (top_cll.name != cell_name): #Don't insert in the top_cell
                cell_index = top_cll.cell_index()
                #print(cell_index)
                #define new origin point 
                #newox = -pitx*(nr_inst_X/2-0.5) 
                #newoy = -pity*(nr_inst_Y/2-0.5)
                #print(newox)
                #print(newoy)
                
                new_instance = pya.CellInstArray( cell_index, pya.Trans(pya.Vector(transx*1000,transy*1000)), pya.Vector(pitchx, 0), pya.Vector(0, pitchy), nr_inst_X, nr_inst_Y)
                top.insert( new_instance ) #insert the cell in the array
                cnt = cnt +1 

        if cnt==0:
            print("Instantiation was unsuccessful. The cell_name needs to be different than the top cell name in "+str(input_filename)+ ". ")

    #write to gds
    layout.write(output_filename)

########RESET DATATYPES 
def reset_datatypes(fstgds_filename, fst_cellname, fst_layer_nr, fst_datatype_nr, snd_layer_nr, snd_datatype_nr, output_filename):
    """
    (void) resets datatypes 
    
    Args:
        :Ngds_filename:    string filename of N(=fst,snd) gds
        :N_cellname:       string name of cell in N(=fst,snd) gds 
        :N_layer_nr:       int layer number in N(=fst,snd) gds  
        :N_datatype_nr:    int datatype number in N(=fst,snd) gds  
        :output_filename:  string filename of output gds
    """
    
    #ly1 with 1st gds file
    ly1 = pya.Layout()
    lmap1 = ly1.read(fstgds_filename)
    cll1 = ly1.cell(fst_cellname)
    lyr1 = ly1.layer(fst_layer_nr,fst_datatype_nr)
    region1 = pya.Region(cll1.shapes(lyr1)) #define region1 as shapes from ly1-lyr1
    
    #create new layer with correct datatype
    ly1.insert_layer(pya.LayerInfo(snd_layer_nr, snd_datatype_nr)) #create with same info
    lyr12 = ly1.layer(snd_layer_nr,snd_datatype_nr) #define the lyr2 in ly1 
    
    cll1.shapes(lyr12).insert(region1)
    cll1.layout().clear_layer(lyr1) #clear ly1-lyr1

    ly1.write(output_filename)
    
########CREATES A CELL WITH THE POLYGONS 
def cell_wpol(cs, cellname):
    """
    Transforms contourf into GDS2 using gdshelpers 
    'cs'       = contours FROM matplotlif contourf function 
    'cellname' = string cellname, e.g. 'TOP' 
    Returns: cell with polygon, multipolygon array  
    """
    
    from gdshelpers.geometry.chip import Cell

    from shapely.geometry import Polygon
    import pickle 
    import numpy as np 

    multipoly=[]

    cell = Cell(cellname)
    
    print(len(cs.collections))

    #EXTRACT THE CONTOURLEVELS AS POLYGONS, AND ADD POLYGONS TO A GDS CELL 
    for ncol,col in enumerate(cs.collections):
        print(ncol)
        # Loop through all polygons that have the same intensity level
        for cps,contour_path in enumerate(col.get_paths()): 
            #print(cps)
            # Create the polygon for this intensity level
            # The first polygon in the path is the main one, the following ones are "holes"
            for ncp,cp in enumerate(contour_path.to_polygons()):
                #print(ncp)
                x = cp[:,0]
                y = cp[:,1]
                new_shape = Polygon([(i[0], i[1]) for i in zip(x,y)])
                new_shape = new_shape.buffer(0)
                print(new_shape.is_valid)
                if ncp == 0:
                    poly = new_shape
                    poly = poly.buffer(0)
                    
                else:
                    # Remove the holes if there are any
                    poly = poly.difference(new_shape)
                
                #Add the Polygon to the layer (ncol+1)
                #NOTE!!! IT LEAVES THE DATATYPE UNDEFINED
                cell.add_to_layer(ncol, poly)
        
            #if cps == 0:
            #    poly1 = new_shape
        
            # Simplify the shapes in the previous layer 
            cell.get_reduced_layer(ncol-1)
        
            # For debug add the polygons to a multipoly array 
            multipoly = np.append(multipoly, poly)
            
    return cell, multipoly
    
##########FUNCTION TO GET LAYERS OF GDS INTO GDSTK POlygon
def inspect_gds2layers_gdstk(filename): 
    """
    Get all layer polygons from gds filename ='filename'
    
    Returns: 
        [0] gstk lib object with all information of cell, layers and shapes 
        [1] pol_dict, dictionary of gdstk Polygons {(layer_nr, datatype_nr): [Polygon1, Polygon2,...]}
        [2] nmpy array with (layer_nr, datatype_nr)
        [3] list with (layer_nr, datatype_nr)
    """
    from gdshelpers.geometry.chip import Cell
    import gdstk
    from shapely.geometry import MultiPolygon, Polygon 

    import numpy as np 
    
    #create cell to store the layers 
    cell_multipol = Cell('TOP')
    
    gdstk.current_library = gdstk.Library()

    lib = gdstk.Library()
    lib = gdstk.read_gds(infile=filename)
    main_cell = lib.top_level()[0]

    #return the layers 
    layerspol = [main_cell.polygons[i-1].layer for i in np.arange(1,np.size(main_cell.polygons)+1)]
    #print(layerspol)
    layers = np.unique(layerspol)
    #print(layers)

    #retun the index where they are unique 
    indexespol = [int(list(layerspol).index(x)) for x in set(list(layerspol))]
    #print(indexespol)
    indexespol = np.array(indexespol)
    indexespol = indexespol.astype(int)
    #print(indexespol.astype(int))

    #this returns the datatypes in the index that are unique ones 
    datatypespol = [main_cell.polygons[i-1].datatype for i in np.arange(1,np.size(main_cell.polygons)+1)]
    datatypespol = np.array(datatypespol)
    datatypes = datatypespol[indexespol]
    #print(datatypes)
    
    #control var
    n=0

    #while it is empty read again 
    while main_cell.polygons ==[]:
        lib = gdstk.read_gds(infile=filename)
        main_cell = lib.top_level()[0]
        n=n+1 #a control to avoid infinite loop 
        if n==10: #if we got to 10 trials 
            print("Cannot read polygons in this GDS file after "+str(n)+" trials.")
            break 
    else: 
        if main_cell.polygons !=[]: 
            print(str(np.size(main_cell.polygons))+" polygons found...")
            pol_dict ={}
            dkeys = [(i,j) for i,j in zip(layers, datatypes)]
            for d in dkeys: 
                pol_dict[d] = []
                pol_dict[d] = [pols for pols in main_cell.polygons if (pols.layer==d[0] and pols.datatype==d[1])]

    for i, el in enumerate(dkeys):
        print(el)
        for k, pols in enumerate(pol_dict[el]):
            print(pol_dict[el][k].layer)
            pol_dict[el][k].layer = int(el[0])
            
    ldtpsa = dkeys 
    ldtps = np.array(ldtpsa)
            
    #pol_dict[1,1]
    return lib, pol_dict, ldtpsa, ldtps
    
  
########CREATES A GDSTK CELL WITH THE POLYGONS USING GDSTK 
def cell_wpol_gdstk(cs, cellname):
    """
    Transforms contourf into GDS2 using gdstk
    
    Args:
        :cs:        contours FROM matplotlif contourf function 
        :cellname:  string cellname, e.g. 'TOP' 
    
    Returns:
        [0] gdstk library
        [1] cell with polygons  
    """
    
    import gdstk

    import pickle 
    import numpy as np 
    
    gdstk.current_library = gdstk.Library()

    lib = gdstk.Library()
    
    main_cell = lib.new_cell(cellname)
    
    
    #EXTRACT THE CONTOURLEVELS AS POLYGONS, AND ADD POLYGONS TO A GDS CELL 
    for ncol,col in enumerate(cs.collections):
        #print(ncol)
        # Loop through all polygons that have the same intensity level
        for cps,contour_path in enumerate(col.get_paths()): 
            #print(cps)
            # Create the polygon for this intensity level
            # The first polygon in the path is the main one, the following ones are "holes"
            for ncp,cp in enumerate(contour_path.to_polygons()):
                #print(ncp)
                x = cp[:,0]
                y = cp[:,1]
                new_shape = gdstk.Polygon([(i[0], i[1]) for i in zip(x,y)], layer=int(ncol+1),datatype= 0)
 
                #Add the Polygon to the layer (ncol+1)
                main_cell.add(new_shape)
        
            # Simplify the shapes in the previous layer 
            #cell.get_reduced_layer(ncol-1)
        
            # For debug add the polygons to a multipoly array 
            #multipoly = np.append(multipoly, poly)
            
    return lib, main_cell
    

##########CHANGE LAYERS ON THE GDS, FROM layerspol TO new_layers
###THIS IS USING DATATYPE AS ZERO by default, can be changed 
def change_layers(fstgds_filename, fst_cellname, layerspol,\
                  #fst_layer_nr, fst_datatype_nr, \
                  new_layers, output_filename):
    """
    (void) Transforms layers from the source layer into the destination
    #by default considers datatypes are int(0), set datatypes to 0 function can be used before
    
    Args:
        :fstgds_filename:    string filename of gds to read
        :fst_cellname:       string name of cell in the gds 
        :layerspol:          array of the layers of the gds file 
        :new_layers:         array of destination layers - MUST HAVE THE SAME CORRESPONDENCE 
        :output_filename:    string filename of output gds
    """
    
    import pya
    
    #ly1 with 1st gds file
    ly1 = pya.Layout()
    lmap1 = ly1.read(fstgds_filename)
    cll1 = ly1.cell(fst_cellname)

    #ly2 with copy of gds file - it will be cleared and written in the correct layers 
    ly2 = pya.Layout()
    lmap2 = ly2.read(fstgds_filename)
    cll2 = ly2.cell(fst_cellname)

    
    #clear the layers in the destination cell 
    for li, lyr in enumerate(layerspol):
        try: 
            cll2.layout().clear_layer(int(li))
        except: 
            print("Could not clear layers... check layer names/nrs in gds file.")
    
    #ntot = np.size(layerspol)
    
    for li, lyr in enumerate(layerspol):
        print(li)
        #select the layer lyr1
        lyr1 = ly1.layer(int(lyr),int(0))
        #print(lyr)
        #define region1 as shapes from lyr1
        region1 = pya.Region(cll1.shapes(lyr1)) 

        #change the shapes from one layer to the other 
        ly2.insert_layer(pya.LayerInfo(int(new_layers[li]), int(0)))
        #select layer in the destination gds, corresponds to layerspol one to one 
        lyr12 = ly2.layer(int(new_layers[li]), int(0))
        #insert the region1 in the selected layer in the destination gds
        cll2.shapes(lyr12).insert(region1)
        
        print("Changed the shapes in layer "+str(lyr)+" into "+str(new_layers[li]))
        
    ly2.write(output_filename)
    
    print("Changed layers - wrote result to " +str(output_filename))
    
    
    
def change_layers_gdspy(fstgds_filename, new_cellname, layerspol, new_layers, output_filename):
    """
    Transforms layers from the source layer (layerpol) into the destination (new_layers) 
    By default considers datatypes are int(0), set datatypes to 0 function can be used before
    Assumes that we have the polygons in the top level of the input gds 
    
    Args: 
        :fstgds_filename:    string filename of gds to read
        :new_cellname:       string name of cell in the gds 
        :layerspol:          array of the layers of the gds file, if it is not the same, leaves the absent layers untouched 
        :new_layers:         array of destination layers - MUST HAVE THE SAME CORRESPONDENCE 
        :output_filename:    string filename of output gds
        

    """
    import gdspy
    import numpy as np
    
    lib = gdspy.GdsLibrary()
    
    #open the inout gds file
    lib.read_gds(fstgds_filename)
    gdspy.current_library = lib

    #get the top cell of the  input gds file 
    currentcell = lib.top_level()[0]

    #get all polygons within the input gds file 
    polygons_dict= currentcell.get_polygons(by_spec=True)
    
    #for info get the layers in the current file 
    listlayers = currentcell.get_layers() 
    
    #make sure the arrays are both int 
    filelayers = np.array(list(listlayers), dtype = int)
    layerpols = np.array(layerspol, dtype=int)

    #new library with the new cell 
    lib2 = gdspy.GdsLibrary() 
    newcell = lib2.new_cell(new_cellname)

    #Check if given array corresponds to the layers within file 
    comp = np.array_equal(filelayers, layerpols)
    if comp is False:
        print("Attention: The layers given " + str(layerpols) + " are NOT the same as the layers in the file " + str(filelayers))
    
    #change the layers
    for ips, ids in zip(layerpols, new_layers): 
        newpols = gdspy.PolygonSet(polygons_dict[(ips, 0)],layer=ids, datatype=0)
        currentcell.remove_polygons(lambda pts, layer, datatype: layer == ips)
        newcell.add(newpols)
        
        print("Changed the shapes in layer "+str(ips)+" into "+str(ids)) 
    
    #layers that are not within the layers list, remain the same 
    for ips in filelayers:
        if ips not in layerpols:
            newpols = gdspy.PolygonSet(polygons_dict[(ips, 0)],layer=ips, datatype=0)
            currentcell.remove_polygons(lambda pts, layer, datatype: layer == ips)
            newcell.add(newpols)

    lib.remove(currentcell)    
    lib2.write_gds(output_filename)
    
    print("Changed layers - wrote result to " +str(output_filename))

     
    
####FUNCTION USING KLAYOUT PYTHON LIB TO RESCALE THE WHOLE LAYOUT
def rescale_layout(readfile, cellname, factor, outfile, divfactor=1, newcellname=None, verbose=True): 
    """
    (void) Rescales the layout
    
    Args:
        :readfile:     string filename of input gds
        :cellname:     string name of cell in the readfile
        :factor:       int with scaling factor multiply   
        :outfile:      string filename of output gds
        :divfactor:    int with scaling factor division, defaults to 1   
        :newcellname:  string name of the cell in the outfile, defaults to cellname 
    """
    import pya
    
    #define layout and read layout from file
    layoutor = pya.Layout()
    lmap = layoutor.read(readfile)
    cell = layoutor.cell(cellname)
    cell_index = layoutor.cell(cellname).cell_index()
    
    if newcellname is not None: 
        newcell = layoutor.rename_cell(cell_index, newcellname)
    
    cell.layout().scale_and_snap(cell, int(1), int(factor), int(divfactor))
    
    #https://www.klayout.de/doc/code/class_Layout.html#m_scale_and_snap
    #void scale_and_snap (Cell cell, int grid, int mult, int div)
    #This method is useful to scale a layout by a non-integer factor. 
    #The scale factor is given by the rational number mult / div. 
    #After scaling, the layout will be snapped to the given grid.
    
    layoutor.write(outfile)
    
    if verbose == True: 
        print("Rescaled "+str(readfile)+  "by a factor of " +str(factor/divfactor))
        print("Saved the result to "+str(outfile))
        
        

def rotate_layout(readfile, cellname, angle, outputfile, transx =0, transy=0): 
    """
    rotate layout by angle, inspired in 
    https://www.klayout.de/forum/discussion/2143/import-a-gds-then-make-into-a-cell-to-rotate 
    
    + 
    translation with vector (transx, transy), default = (0,0)
    """
    ly = pya.Layout()
    ly.read(readfile)

    org_top = ly.top_cell()
    new_top = ly.create_cell("TEMP")
    new_top.insert(pya.DCellInstArray(org_top.cell_index(), pya.DCplxTrans(1.0, angle, False, pya.DVector(transx, transy))))
    
    new_top.flatten(1)
    ly.rename_cell(new_top.cell_index(), cellname)
    
    ly.write(outputfile)
    
    print("Rotated " +readfile + " by " +str(angle)+ " degrees. Saved in " + outputfile)
        
    
####FUNCTION TO MAKE THE DIFFS BETWEEN THE LAYERS IN THE LAYERS ARRAY 
def diffs_layers_arrays(readfile, cellname, layerspol1, datatypes1, layerspol2, datatypes2, outfile): 
    """
    (void) Sequentially makes the difference from one layer to the other, from layerspol1 to layerspol2
    
    Args: 
        :readfile:     string filename of input gds
        :cellname:     string name of cell 
        :layerspol1:   numpy array with all layer that will be substracted
        :datatypes1:   numpy array with all datatypes   that will be subtracted 
        :layerspol2:   numpy array with all layer where we will subtract
        :datatypes2:   numpy array with all datatypes where we will substract
        :outfile:      string filename of output gds
    """
    import pya
    import numpy as np 

    #define layout and read layout from file
    layoutor = pya.Layout()
    lmap = layoutor.read(readfile)
    cell = layoutor.cell(cellname)

    #for all the layers in layerspol array 
    for lyrs1, lyrs2, dtps1,dtps2 in zip(layerspol1, layerspol2, datatypes1, datatypes2):
        #print(lyrs)
        #print(datatypes[lyrs])
        
        npd = int(lyrs1)
        dts = int(dtps1)
        minval = np.min(layerspol1)
        maxval = np.max(layerspol2)
        
        layer1 = layoutor.layer(npd,dts)
        region1 = pya.Region(cell.shapes(layer1))
        
        #print(npd)
        #make the collection from all areas that already have a shape
        if npd > 0: 
            #print(npd)
            nar = np.flip(np.arange(minval,npd))
            print(nar)
            dtar = np.zeros(len(nar))
             
            for nps, dtps in zip(nar, dtar):
                print(nps)
                layer1s = layoutor.layer(int(nps),int(dtps))
                #region with all the shapes within the TOP cell
                region1s = pya.Region(cell.shapes(layer1s))
                region1 = region1 + region1s 
        
        #print(npd)
        #if npd>2: 
        #    print(np.arange(0,float(npd)))
    
        np1 = int(lyrs2)
        dts1 = int(dtps2)
        layer2 = layoutor.layer(np1,dts1)
        #region with all the shapes within the TOP cell
        region2 = pya.Region(cell.shapes(layer2))
    
        #make the difference of the layer2 on layer1
        result = region2-region1 
    
        #temp save of the resull layer (to be cleared)
        resultlay = layer2
        cell.layout().clear_layer(layer2) #clear the results layer
    
        #insert shapes from boolean into the resultslayer 
        cell.shapes(resultlay).insert(result)
    
    layoutor.write(outfile)
    
    print("Substracted "+str(layerspol1)+" in " + str(layerspol2)+" of the file "+str(readfile))
    print("Saved the result to "+str(outfile))
    

def cell_wpol_gdspy(cs, cellname, prec=1e-6, mpoints=1e9):
    """
    Cell made with cut polygons from the z profile

    Args:
        :cs:        contours FROM matplotlif contourf function 
        :cellname:  string cellname, e.g. 'TOP' 
    
    Returns:
        [0] gdspy library
        [1] cell with polygons  
  
    ##By default the levels start at 0 and go to the number of levels in the contourplot 
    """
    
    from shapely import geometry
    import numpy as np
    import gdspy

    #get collections from contour 
    collecs = cs.collections 

    # lib of the gdsii file 
    lib = gdspy.GdsLibrary()
    gdspy.current_library = gdspy.GdsLibrary()

    cell = lib.new_cell(cellname)

    ncolec = len(collecs)
    print("Passing contours into GDS. ")

    #go through all the elements in collections (gray levels)
    for ncol,col in enumerate(collecs):
        print(ncol)
        # Loop through all contours that have the same gray level
        paths = col.get_paths() 
        #lenpat = len(paths)

        #go through the paths of the contours 
        for ec, contour_path in enumerate(paths): 
            #get the polygons at certain gray level 
            arr = contour_path.to_polygons()
            polset = gdspy.PolygonSet(arr, layer=int(ncol),datatype=int(0))
            #print("pols "+str(len(polset.polygons)))
            
            #layers and datatypes of polygon set 
            layes = polset.layers[0]
            dts  = polset.datatypes[0]

            #print(polset.polygons)
            if len(polset.polygons)==1:
                pols = gdspy.Polygon(polset.polygons[0], layer=int(layes), datatype=int(dts))
                #print(pols)
                cell.add(pols)

            for ncp,cp in enumerate(arr):
                #print(ncp)
                x = cp[:,0]
                y = cp[:,1]
                new_shape = gdspy.Polygon([(i[0], i[1]) for i in zip(x,y)])
                
                if ncp == 0:
                    poly = new_shape
                elif ncp==1:
                    de   = gdspy.boolean(poly,new_shape, "not", precision=prec, max_points=mpoints, layer=int(ncol), datatype=int(0)) 
                    cell.add(de)
                    
               
    return lib, cell 
    
 

def cell_wpol_gdspy_fast(cs, cellname, prec=1e-6, mpoints=1e9):
    """
    Cell made with cut polygons from the z profile 

    Args:
        :cs:        contours FROM matplotlif contourf function 
        :cellname:  string cellname, e.g. 'TOP' 
    
    Returns:
        [0] gdspy library
        [1] cell with polygons 
    ##By default the levels start at 0 and go to the number of levels in the contourplot 
    #CAREFUL, this is better for estimation ONLY as we need enough resolution in the contourplot for this approach to work properly 
    """
    
    from shapely import geometry
    import numpy as np
    import gdspy

    #get collections from contour 
    collecs = cs.collections 

    # lib of the gdsii file 
    lib = gdspy.GdsLibrary()
    gdspy.current_library = gdspy.GdsLibrary()

    cell = lib.new_cell(cellname)

    ncolec = len(collecs)
    print("Passing contours into GDS. ")

    #go through all the elements in collections (gray levels)
    for ncol,col in enumerate(collecs):
        # Loop through all contours that have the same gray level
        paths = col.get_paths() 
        lenpat = len(paths)

        #go through the paths of the contours 
        for ec, contour_path in enumerate(paths): 
            #get the polygons at certain gray level 
            arr = contour_path.to_polygons()
            polset = gdspy.PolygonSet(arr, layer=int(ncol),datatype=int(0))
            #print("pols "+str(len(polset.polygons)))
            #layers and datatypes of polygon set 
            layes = polset.layers[0]
            dts  = polset.datatypes[0]

            #print(polset.polygons)
            if len(polset.polygons)==1:
                pols = gdspy.Polygon(polset.polygons[0], layer=int(layes), datatype=int(dts))
                #print(pols)
                cell.add(pols)

            #if there is also another polygons, subtract 
            if len(polset.polygons)==2: 
                #print(len(polset.polygons))
                for ips in np.arange(0,len(polset.polygons)-1):
                    pol1 = gdspy.Polygon(polset.polygons[ips])
                    pol2 = gdspy.Polygon(polset.polygons[ips+1])
                    de   = gdspy.boolean(pol1,pol2, "not", precision=prec, max_points=mpoints, layer=int(layes), datatype=int(dts)) 
                    cell.add(de)
            
    return lib, cell