"""
GDS converter module

"""
import gdspy
import numpy as np

# import pyMOE as moe
from pyMOE.aperture import Aperture
from pyMOE.utils import progress_bar, Timer

import matplotlib.pyplot as plt 


def count_vertices(pols):
    """ Counts vertices of polygons
    
    Args: 
        polygon or polygon set 
    
    Returns: 
        number of vertices 
    """
    sum = 0 
    for p in pols:
        sum += p.shape[0]
    return sum

def merge_polygons(polygons, layer=0, assume_non_overlap=True, break_vertices=250, verbose=True, ):
    """
    Merge polygons function receives a list of polygons or polygon set and 
    will iteratively merge consecutive polygons, assuming to be in the same layer. If we assume that polygons are
    sequentially located in the matrix, then it is expected that consecutive polygons
    can be merged, while polygons far apart are probably not in the same cluster.
    
    The boolean operation compares the new polygon with the existing operand, and it
    becomes slower with more polyons/vertices already included in the merged set.
    To optimize this, we consider a break_vertices threshold where the merged set of 
    polygons is broken into a separate list to continue the merging. This will result 
    in a merged set that could be smaller, but is much faster to execute, scaling with 
    N instead of N^2.
    
    Args:
        polygons: must be a list of polygons, polygonset, rectangles etc that the boolean operation accepts
        layer: default 0 and ignored. The merged list of merged polygons will have the same layer as the input polygons (assumed the same for all)
        assume_non_overlap: default True, if True the function will break the merged list of polygons when reaching break_vertices
        break_vertices: approximate number of vertices to consider in the boolean operation before breaking the list.
        verbose: default True. Prints the progress bar.
        
    Returns:
        list of merged polygons or polygonsets.
    """
    
    # TODO assert that polygons are list of polygons or polygonset etc...
    
    total_polygons = len(polygons)
    merged = gdspy.PolygonSet([])
    list_polygonsets = []
    
    restart_merge = True
    for i,pol in enumerate(polygons):

        # pops one polygon to be considered in the next boolean operation
        if restart_merge:
            merged = pol
            if verbose:
                progress_bar(i/total_polygons)
            restart_merge = False
            continue
        
        layer = pol.layers[0]
        merged = gdspy.boolean(merged, pol, "or", layer=layer)
        
        if assume_non_overlap:
            # Breaks the polygons considered for boolean operation
            if (count_vertices(merged.polygons)>=break_vertices) & (i <total_polygons-2):
                list_polygonsets.append(merged)
                restart_merge = True
    list_polygonsets.append(merged)
    progress_bar(1)
    
    return list_polygonsets
    
##############################################################################################################################
    
def change_layers_gdspy(fstgds_filename, new_cellname, layerspol, new_layers, output_filename):
    """
    Transforms layers from the source layer (layerpol) into the destination (new_layers) 
    By default considers datatypes are int(0), set datatypes to 0 function can be used before
    Assumes that we have the polygons in the top level of the input gds 
    
    Args: 
        'fstgds_filename'   : string filename of gds to read
        'new_cellname'      : string name of cell in the gds 
        'layerspol'         : array of the layers of the gds file, if it is not the same, leaves the absent layers untouched 
        'new_layers'              : array of destination layers - MUST HAVE THE SAME CORRESPONDENCE 
        'output_filename'   : string filename of output gds
        

    """
    
    
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
    comp = np.array_equal(filelayers, layerspol)
    if comp is False:
        print("Attention: The layers given " + str(layerspol) + " are NOT the same as the layers in the file " + str(filelayers))
    
    #change the layers
    for ips, ids in zip(layerspol, new_layers): 
        newpols = gdspy.PolygonSet(polygons_dict[(ips, 0)],layer=ids, datatype=0)
        currentcell.remove_polygons(lambda pts, layer, datatype: layer == ips)
        newcell.add(newpols)
        
        print("Changed the shapes in layer "+str(ips)+" into "+str(ids)) 
    
    #layers that are not within the layers list, remain the same 
    for ips in filelayers:
        if ips not in layerspol:
            newpols = gdspy.PolygonSet(polygons_dict[(ips, 0)],layer=ips, datatype=0)
            currentcell.remove_polygons(lambda pts, layer, datatype: layer == ips)
            newcell.add(newpols)

    lib.remove(currentcell)    
    lib.write_gds(output_filename)
    
    print("Changed layers - wrote result to " +str(output_filename))



def cell_wpol_gdspy(cs, cellname, prec=1e-6, mpoints=1e9):
    """
    ###Initially coming from gds_klops 
    ###TODO: migrate functions from the gds_klops (implemented with various libraries) to gdsp y 
    
    Cell made with cut polygons from the z profile 
    'cs'       = contours FROM matplotlif contourf function
    'cellname' = string cellname, e.g. 'TOP' 
    Returns:[0] gdspy library, [1] cell with polygons  
    ##By default the levels start at 0 and go to the number of levels in the contourplot 
    """

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
    
################################################################################################################################3

class GDSMask():
    """
    Class GDSMask:
        Creates a class to integrate the GDS library and layout corresponding to a mask aperture.
        Receives an aperture and provides methods to calculate the corresponding GDS layout
    
    Args:
        mask: aperture object 
    
    Methods:
        mask
        levels: number of unique discretized levels in mask 
        layers: list of layers in layout
        cells: list of cells in layout
        total_polygons: total number of polygons in layout in all layers
        total_vertices:  total number of vertices in all polygons in all layers

        create_layout(): Creates layout (raster or vector or edges etc todo)
        merge(layer): merges polygons in layer
        plot(): plots mask
        viewer(): Opens LayoutViewer of gdspy
        save_gds(filename): saves layout to gds file

    """
    def __init__(self, mask, units=1e-6, precision=1e-9, verbose=True):
        assert type(mask) is Aperture, "aperture must be of type Aperture"
        self.mask = mask
        self.gdslib = None
        self.units = units
        self.precision = precision
        self.verbose=verbose
        self._gdslib_init_error = "Error: gdsmask not created yet. Run GDSMask.create_layout()"
        self.layers = None
    
    @property
    def levels(self):
        levels = self.mask.levels
        assert levels is not None, "Cannot access GDSMask.levels as aperture is not yet discretized"
        return self.mask.levels
    
    @property
    def aperture(self):
        aperturevar = self.mask.aperture_discretized
        assert aperturevar is not None, "Cannot access GDSMask.aperture as aperture is not yet discretized"
        return self.mask.aperture_discretized
        
    
    @property
    def total_polygons(self):
        assert self.gdslib is not None, "%s"%self._gdslib_init_error
        
    @property
    def total_vertices(self):
        assert self.gdslib is not None, "%s"%self._gdslib_init_error
    

    
    def _init_layout(self):
        """ Initializes the gdslib"""
        
        self.gdslib = gdspy.GdsLibrary()
        
        
    def view_layout(self):
        """ Runs gdspy.LayoutViewer"""
        
        assert self.gdslib is not None, "Must create_layout() from discretized aperture first"
        gdspy.LayoutViewer(self.gdslib)
    
    def write_gds(self, filename, cells=None, timestamp=None, binary_cells=None):
        """ Writes layout to gds file using gdspy library"""
        self.gdslib.write_gds(filename, cells, timestamp, binary_cells)
        print("Saved %s"%(filename))
        
        
        
#<<<<<<< dev-refactor-compatibility
#    def create_layout(self, mode="raster", cellname='TOP', merge=True, break_vertices=250):
#=======
    def create_layout(self, mode="raster", cellname='TOP', merge=False, break_vertices=250):
#>>>>>>> dev-refactor
        """
        Creates GDS layout of the discretized aperture
        
        Args:
            mode: default Raster. TO DO
            cellname: name of the topcell to include all the merged polygons
            merge: default False. If True, will merge the individual pixel polygons 
            break_vertices: threshold value to speed up the merging of polygons
        
        Returns:
            gdslib: library with topcell will update the class internal gdslib
        """

        if self.gdslib is None:
            self._init_layout()
        assert self.aperture is not None, "Cannot create_layout() as aperture is not yet discretized"
        
        
        if mode is "raster":
            return self._create_layout_raster(cellname=cellname, merge=merge, break_vertices=break_vertices)
        #else:
        #    raise ValueError('Supported modes are "raster"')
        elif mode is "contour": 
            return self._create_layout_contour(cellname = cellname)
        else: 
            raise ValueError("Unsuported option!")
        
        
        
    def _create_layout_raster(self, cellname='mask', merge=True, break_vertices=250):
        """
        Creates the gds layout using raster mode where each data point is a pixel rectangle to
        be defined in the layout
        
        """
        self.gdslib = gdspy.GdsLibrary()
        
        cell = gdspy.Cell(cellname,exclude_from_current=True)
        
        self.layers = np.arange(len(self.levels))
        total_layers = len(self.layers)
        total_points = self.mask.shape[0]*self.mask.shape[1]
        
        size_x, size_y = self.mask.shape
        XX = self.mask.XX
        YY = self.mask.YY

        half_pixel_x = self.mask.pixel_x/2
        half_pixel_y = self.mask.pixel_y/2

        # normalize to units:
        XX = XX/self.units
        YY = YY/self.units
        half_pixel_x = half_pixel_x/self.units
        half_pixel_y = half_pixel_y/self.units
        datatype = 0
        
        # Creates top cell that will reference the remaining cells
        topcell = gdspy.Cell(cellname, exclude_from_current=True)
    
        # Make evaluation of computational effort
        if self.verbose:
            print("Mask has %d number of points distributed in %d layers"%(total_points, total_layers))       

        with Timer("Total time converting to GDS"):
            # Creates a list of cells where each cell will correspond to a layer, to enable merging afterwards.
            list_cells = []
            for layer in self.layers:
                cellname = "tempcelllayer%d"%(layer)

                cell = gdspy.Cell(cellname,exclude_from_current=True)
                list_cells.append(cell)

            # Iterates to create a layout in each layer corresponding to one discretized level
            if self.verbose:
                print("Creating individual pixel polygons")
            with Timer("Create Polygons"):
                for i in range(size_x):
                    for j in range(size_y):

                        current_point = i*size_x+j
                        aperture_value = self.aperture[i,j]

                        current_layer = aperture_value
                        x = XX[i,j]
                        y = YY[i,j]
                        
                        # Creates corners of rectangle with center at the data value
                        rectangle_first_corner = (x-half_pixel_x,y-half_pixel_y)
                        rectangle_second_corner = (x+half_pixel_x,y+half_pixel_y)
                        
                        #Creates rectangle and adds it to the cell corresponding to the current layer
                        rect = gdspy.Rectangle(rectangle_first_corner, rectangle_second_corner, current_layer, datatype)
                        cell = list_cells[current_layer]
                        
                        # Saves each rectangle to a separate cell so we can merge more easily afterwards
                        cell.add(rect)
                        
                    if self.verbose:
                        progress_bar(current_point/total_points)
                        
                if self.verbose: 
                    progress_bar(1)

            list_cell_refs = []
            list_merged_polygons = []
            for i, cell in enumerate(list_cells):
                
                # Run merging of polygons 
                if merge:
                    total_polygons = len(cell.polygons)
                    print("Merging layer %d of %d with %d polygons:"%(i,len(list_cells)-1, total_polygons))
                    with Timer():
                        merged = merge_polygons(cell.polygons, break_vertices=break_vertices, verbose=self.verbose)
                else:
                    merged = cell.polygons
                for m in merged:
                    list_merged_polygons.append(m)
                # remove temporary cells from memory
                self.gdslib.remove(cell)
                
            # add all polygons to topcell
            topcell.add(list_merged_polygons)
            # add topcell to library
            self.gdslib.add(topcell)

            return self.gdslib
            
            
    def _create_layout_contour(self, cellname='TOP'):
        """
        Creates the gds layout using contour mode via matplotlib library 
        """
        self.gdslib = gdspy.GdsLibrary()
        
        cell = gdspy.Cell(cellname,exclude_from_current=True)
        
        self.layers = np.arange(len(self.levels))
        total_layers = len(self.layers)
        total_points = self.mask.shape[0]*self.mask.shape[1]
        
        size_x, size_y = self.mask.shape
        XX = self.mask.XX
        YY = self.mask.YY

        half_pixel_x = self.mask.pixel_x/2
        half_pixel_y = self.mask.pixel_y/2

        # normalize to units:
        XX = XX/self.units
        YY = YY/self.units
        half_pixel_x = half_pixel_x/self.units
        half_pixel_y = half_pixel_y/self.units
        datatype = 0
        
        # Creates top cell that will reference the remaining cells
        topcell = gdspy.Cell(cellname, exclude_from_current=True)
    

        with Timer("Total time converting to GDS"):

            # Iterates to create a layout in each layer corresponding to one discretized level
            if self.verbose:
                print("Creating contours ")
            with Timer("Create Contours"):
            ## consider to make the discretization as for the raster, to have fixed number of levels 
            ###TODO: change for the actual position of in zlevs 
                cs = plt.contourf(XX,YY,self.mask.aperture_discretized, len(self.levels))

            self.gdslib, cell1 = cell_wpol_gdspy(cs, 'TOP', prec = self.precision, mpoints=1e9)
           

            return self.gdslib
