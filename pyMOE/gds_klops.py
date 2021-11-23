####IMPORTANT: these functions work on the klayout-0.26.11 version
####There might be some inconsistencies for later versions 
import pya 


####### #MERGE FUNCTION 
def merge_layer(readfile, cellname, layer_nr, datatype_nr, outputfile): 
    """
    (void) Merges all shapes within gds (only tested gds with single layer)
    'readfile'    = string filename of input gds
    'cellname'    = string name of cell 
    'layer_nr'    = int    layer number 
    'datatype_nr' = int    datatype number
    'outputfile'  = string filename of output gds
    """
    layoutor = pya.Layout()
    lmap = layoutor.read(readfile)
    cell = layoutor.cell(cellname)
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
               output_filename):
    """
    (void) Imports shapes from fst gds file INTO snd gds file 
    'Ngds_filename'   = string filename of N(=fst,snd) gds
    'N_cellname'      = string name of cell in N(=fst,snd) gds 
    'N_layer_nr'      = int layer number in N(=fst,snd) gds  
    'N_datatype_nr'   = int datatype number in N(=fst,snd) gds  
    'output_filename' = string filename of output gds
    """
    
    #ly1 with 1st gds file
    ly1 = pya.Layout()
    lmap1 = ly1.read(fstgds_filename)
    cll1 = ly1.cell(fst_cellname)
    lyr1 = ly1.layer(fst_layer_nr,fst_datatype_nr)
    region1 = pya.Region(cll1.shapes(lyr1)) #define region1 as shapes from ly1-lyr1
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
    cll2.layout().clear_layer(lyr2) #clear ly2-lyr2

    #create new layer in ly1 to receive shapes from ly2-lyr2 
    ly1.insert_layer(pya.LayerInfo(snd_layer_nr, snd_datatype_nr)) #create with same info
    lyr12 = ly1.layer(snd_layer_nr,snd_datatype_nr) #define the lyr2 in ly1 
    #USE ONLY ONE CELL, here the existing one in ly1
    cll1.shapes(lyr12).insert(region2) #cp shapes from ly2 region2 to ly1 lyr12
    cll1.shapes(lyr1).insert(region1)

    ly1.write(output_filename)
    
    print("Imported "+sndgds_filename+" into "+fstgds_filename+". Output file "+output_filename+" .")

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
    'cell_name'       = string name of cell
    'input_filename'  = string name of the input gds 
    'nr_inst_X'       = int number of repeated objects in x
    'nr_inst_Y'       = int number of repeated objects in y 
    'pitx'            = int pitch in x in um
    'pity'            = int pitch in y in um 
    'output_filename' = string filename of output gds
    """

    layout = pya.Layout()

    #pitches in x and y
    pitchx = pitx*1000 # pitx in um
    pitchy = pity*1000 # pit um

    #create cell at top 
    top = layout.create_cell(cell_name)

    #gds files to read (could also be a list)
    gds_files = [input_filename]

    for file in gds_files:
        layout.read(file) #read the files

        for top_cll in layout.top_cells():
            if (top_cll.name != cell_name): #Don't insert in the top_cell
                cell_index = top_cll.cell_index()
                #define new origin point 
                #newox = -pitx*(nr_inst_X/2-0.5) #-18498*1000
                #newoy = -pity*(nr_inst_Y/2-0.5)
                #print(newox)
                #print(newoy)
                
                new_instance = pya.CellInstArray( cell_index, pya.Trans(pya.Vector(transx*1000,transy*1000)), pya.Vector(pitchx, 0), pya.Vector(0, pitchy), nr_inst_X, nr_inst_Y)
                top.insert( new_instance ) #insert the cell in the array

    #write to gds
    layout.write(output_filename)

########RESET DATATYPES 
def reset_datatypes(fstgds_filename, fst_cellname, fst_layer_nr, fst_datatype_nr, snd_layer_nr, snd_datatype_nr, output_filename):
    """
    (void) resets datatypes 
    'Ngds_filename'   = string filename of N(=fst,snd) gds
    'N_cellname'      = string name of cell in N(=fst,snd) gds 
    'N_layer_nr'      = int layer number in N(=fst,snd) gds  
    'N_datatype_nr'   = int datatype number in N(=fst,snd) gds  
    'output_filename' = string filename of output gds
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
    Cell made with cut polygons 
    'cs' = contours FROM matplotlif 'contourf' 
    'cellname' = cellname 
    Returns: cell with polygon, multipolygon array  
    """
    
    from gdshelpers.geometry.chip import Cell

    from shapely.geometry import Polygon
    import pickle 
    import numpy as np 

    multipoly=[]

    cell = Cell(cellname)

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
                new_shape = Polygon([(i[0], i[1]) for i in zip(x,y)])
                if ncp == 0:
                    poly = new_shape
                else:
                    # Remove the holes if there are any
                    poly = poly.difference(new_shape)
                
                #Add the Polygon to the layer (ncol+1)
                #NOTE!!! IT LEAVES THE DATATYPE UNDEFINED
                cell.add_to_layer(ncol+1, poly)
        
            if cps == 0:
                poly1 = new_shape
        
            # Simplify the shapes in the previous layer 
            cell.get_reduced_layer(ncol-1)
        
            # For debug add the polygons to a multipoly array 
            multipoly = np.append(multipoly, poly)
            
    return cell, multipoly
    
##########FUNCTION TO GET LAYERS OF GDS INTO GDSTK POlygon
def inspect_gds2layers_gdstk(filename): 
    """
    Get all layer polygons from gds filename ='filename'
    returns 
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
    
  
########CREATES A GDSTK CELL WITH THE POLYGONS 
def cell_wpol_gdstk(cs, cellname):
    """
    Cell made with cut polygons 
    'cs'       = contours FROM matplotlif 'contourf' 
    'cellname' = cellname 
    Returns: cell with polygon, multipolygon array  
    """
    
    from gdshelpers.geometry.chip import Cell
    import gdstk
    from shapely.geometry import Polygon
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
                #NOTE!!! IT LEAVES THE DATATYPE UNDEFINED
                main_cell.add(new_shape)
                #cell.add_to_layer(ncol+1, poly)
        
        
            # Simplify the shapes in the previous layer 
            #cell.get_reduced_layer(ncol-1)
        
            # For debug add the polygons to a multipoly array 
            #multipoly = np.append(multipoly, poly)
            
    return lib, main_cell
    
##########CHANGE LAYERS ON THE GDS, FROM layerspol TO gvts
###THIS IS USING DATATYPE AS ZERO
def change_layers(fstgds_filename, fst_cellname, layerspol,\
                  #fst_layer_nr, fst_datatype_nr, \
                  gvts, output_filename):
    """
    (void) Transforms layers from the source layer into the destination
    'fstgds_filename'   = string filename of gds to read
    'fst_cellname'      = string name of cell in the gds 
    'layerspol'         = array of the layers of the gds file 
    'gvts'              = array of destination layers - MUST HAVE THE SAME CORRESPONDENCE 
    'output_filename'   = string filename of output gds
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
        cll2.layout().clear_layer(int(li))
    
    #ntot = np.size(layerspol)
    
    for li, lyr in enumerate(layerspol):
        #select the layer lyr1
        lyr1 = ly1.layer(int(lyr),int(0))
        print(lyr)
        #define region1 as shapes from lyr1
        region1 = pya.Region(cll1.shapes(lyr1)) 

        #change the shapes from one layer to the other 
        ly2.insert_layer(pya.LayerInfo(int(gvts[li]), int(0)))
        #select layer in the destination gds, corresponds to layerspol one to one 
        lyr12 = ly2.layer(int(gvts[li]), int(0))
        #insert the region1 in the selected layer in the destination gds
        cll2.shapes(lyr12).insert(region1)
        
        print("Changed the shapes in layer "+str(lyr)+" into "+str(gvts[li]))
        
    ly2.write(output_filename)
    
    print("Changed layerspol layer to gvts - wrote result to " +str(output_filename))
    
####FUNCTION USING KLAYOUT PYTHON LIB TO RESCALE THE WHOLE LAYOUT
def rescale_layout(readfile, cellname, factor, outfile, divfactor=1): 
    """
    (void) Rescales the layout
    'readfile'    = string filename of input gds
    'cellname'    = string name of cell 
    'factor'      = int with scaling factor multiply   
    'outfile'     = string filename of output gds
    'divfactor'   = int with scaling factor division, defaults to 1   
    """
    import pya

    #define layout and read layout from file
    layoutor = pya.Layout()
    lmap = layoutor.read(readfile)
    cell = layoutor.cell(cellname)
    
    cell.layout().scale_and_snap(cell, int(1), int(factor), int(divfactor))
    
    #https://www.klayout.de/doc/code/class_Layout.html#m_scale_and_snap
    #void scale_and_snap (Cell cell, int grid, int mult, int div)
    #This method is useful to scale a layout by a non-integer factor. 
    #The scale factor is given by the rational number mult / div. 
    #After scaling, the layout will be snapped to the given grid.
    
    layoutor.write(outfile)
    
    print("Rescaled "+str(readfile)+  "by a factor of " +str(factor/divfactor))
    print("Saved the result to "+str(outfile))
    
####FUNCTION TO MAKE THE DIFFS BETWEEN THE LAYERS IN THE LAYERS ARRAY 
def diffs_layers_arrays(readfile, cellname, layerspol1, datatypes1, layerspol2, datatypes2, outfile): 
    """
    (void) Sequentially makes the difference from one layer to the other, from layerspol[1] to layerspol[0]
    'readfile'    = string filename of input gds
    'cellname'    = string name of cell 
    'layerspol1'   = numpy array with all layer that will be substracted
    'datatypes1'   = numpy array with all datatypes   that will be subtracted 
    'layerspol2'   = numpy array with all layer where we will subtract
    'datatypes2'   = numpy array with all datatypes where we will substract
    'outfile'     = string filename of output gds
    """
    import pya

    #define layout and read layout from file
    layoutor = pya.Layout()
    lmap = layoutor.read(readfile)
    cell = layoutor.cell(cellname)

    #for all the layers in layerspol array 
    for lyrs1, lyrs2, dtps1,dtps2 in zip(layerspol1, layerspol2, datatypes1, datatypes2):
        #print(lyrs)
        #print(datatypes[lyrs])
        
        np = int(lyrs1)
        dts = int(dtps1)
        layer1 = layoutor.layer(np,dts)
    
        #region with all the shapes within the TOP cell
        region1 = pya.Region(cell.shapes(layer1))
    
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
