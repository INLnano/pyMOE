"""
GDS converter module

"""





import gdspy
import numpy as np


# import pyMOE as moe
from pyMOE.aperture import Aperture
from pyMOE.utils import progress_bar, Timer


def count_vertices(pols):
    """ Counts vertices of polygons
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



class GDSMask():
    """
    Class GDSMask:
        Creates a class to integrate the GDS library and layout corresponding to a mask aperture.
        Receives an aperture and provides methods to calculate the corresponding GDS layout
    
    Args:
        mask:
    
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
        
        
        
    def create_layout(self, mode="raster", cellname='mask', merge=True, break_vertices=250):
        """
        Creates GDS layout of the discretized aperture
        
        Args:
            mode: default Raster. TO DO
            cellname: name of the topcell to include all the merged polygons
            merge: default True. If True, will merge the individual pixel polygons 
            break_vertices: threshold value to speed up the merging of polygons
        
        Returns:
            gdslib: library with topcell will update the class internal gdslib
        """

        if self.gdslib is None:
            self._init_layout()
        assert self.aperture is not None, "Cannot create_layout() as aperture is not yet discretized"
        
        
        if mode is "raster":
            return self._create_layout_raster(cellname=cellname, merge=merge, break_vertices=break_vertices)
        else:
            raise ValueError('Supported modes are "raster"')
        
        
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
