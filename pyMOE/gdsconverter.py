"""
gdsconverter.py
GDS converter module

"""
import numpy as np

from pyMOE.aperture import Aperture
from pyMOE.utils import progress_bar, Timer

import matplotlib.pyplot as plt 
from scipy.constants import milli, nano, micro

from scipy.spatial.distance import cdist

import pya


def count_vertices(pols):
    """ Counts vertices of polygons
    
    Args: 
        :pols: polygon or polygon set 
    
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
        :polygons:              must be a list of polygons, polygonset, rectangles etc that the boolean operation accepts
        :layer:                 default 0 and ignored. The merged list of merged polygons will have the same layer as the input polygons (assumed the same for all)
        :assume_non_overlap:    default True, if True the function will break the merged list of polygons when reaching break_vertices
        :break_vertices:        approximate number of vertices to consider in the boolean operation before breaking the list.
        :verbose:               default True. Prints the progress bar.
        
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

def cell_wpol_gdspy(cs, cellname, prec=1e-6, mpoints=1e9):
    """
    ###Initially coming from gds_klops 
    ###TODO: migrate functions from the gds_klops (implemented with various libraries) to gdsp y 
    
    Cell made with cut polygons from the z profile 
    
    Args:
        :cs:        contours FROM matplotlif contourf function
        :cellname:  string cellname, e.g. 'TOP' 
    
    Returns:
        :[0]:       gdspy library, 
        :[1]:       cell with polygons  
    
    ##By default the levels start at 0 and go to the number of levels in the contourplot 
    """
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
            if polset.layers!=[]: 
                layes = polset.layers[0]
            if polset.datatypes!=[]:
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
        :mask: aperture object 
    
    Methods:
        :levels:            number of unique discretized levels in mask 
        :layers:            list of layers in layout
        :cells:             list of cells in layout
        :total_polygons:    total number of polygons in layout in all layers
        :total_vertices:    total number of vertices in all polygons in all layers

        :create_layout():   Creates layout (raster or vector or edges etc todo)
        :merge(layer):      merges polygons in layer
        :plot():            plots mask
        :viewer():          Opens LayoutViewer of gdspy
        :save_gds(filename): saves layout to gds file

    """
    def __init__(self, mask, units=1e-6, precision=1e-9, verbose=True):
        print(type(mask ))
        # assert type(mask) is Aperture, "aperture must be of type Aperture"
        self.mask = mask
        self.gdslib = None
        self.units = units
        self.precision = precision
        self.verbose=verbose
        self._gdslib_init_error = "Error: gdsmask not created yet. Run GDSMask.create_layout()"
        self.layers = None
        self.grayvalues = None
    
    @property
    def levels(self):
        levels = self.mask.levels
        assert levels is not None, "Cannot access GDSMask.levels as aperture is not yet discretized. Please run aperture.discretize(n) before, where n is the number of levels."
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
        """ Writes layout to gds file using gdspy library
        
        DEPRECATED """
    #     self.gdslib.write_gds(filename, cells, timestamp, binary_cells)
    #     print("Saved %s"%(filename))
    
    # def write_gds(self, filename, cells=None, timestamp=None, binary_cells=None):
        self.write_layout(filename)        


    def write_layout(self, filename):
        """ Writes layout togds file using klayout pya library"""
        print("Saving file to %s"%(filename))
        with Timer("Saving GDS file"):
            self.layout.write(filename,)
        print("Saved %s"%(filename))
        
        

    def create_layout(self, mode="raster", cellname='TOP', merge=False, break_vertices=250):
        """
        Creates GDS layout of the discretized aperture
        
        Args:
            :mode:          default Raster. (can also accept contour)
            :cellname:      name of the topcell to include all the merged polygons
            :merge:         default False. If True, will merge the individual pixel polygons 
            :break_vertices: threshold value to speed up the merging of polygons
        
        Returns:
            # :gdslib: l      ibrary with topcell will update the class internal gdslib
        """

        # if self.gdslib is None:
        #     self._init_layout()
        assert self.aperture is not None, "Cannot create_layout() as aperture is not yet discretized"
        
        
        if mode == "raster":
            return self._create_layout_raster_klayout(cellname=cellname, merge=merge, break_vertices=break_vertices)
        elif mode == "contour": 
            return self._create_layout_contour(cellname = cellname)
        else: 
            raise ValueError("Unsuported option!")
        
        
        
    def _create_layout_raster(self, cellname='mask', merge=True, break_vertices=250):
        """
        Creates the gds layout using raster mode where each data point is a pixel rectangle to
        be defined in the layout
        
        """
        import gdspy

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
                        
                        if not np.isnan(aperture_value):
                            current_layer = int(np.rint(aperture_value))
                        
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
        import gdspy

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
        #topcell = gdspy.Cell(cellname, exclude_from_current=True)
    

        with Timer("Total time converting to GDS"):

            # Iterates to create a layout in each layer corresponding to one discretized level
            if self.verbose:
                print("Creating contours ")
            with Timer("Create Contours"):
            ## consider to make the discretization as for the raster, to have fixed number of levels 
            ###TODO: change for the actual position of in zlevs 
                plt.ioff()
                cs = plt.contourf(XX,YY,self.mask.aperture_discretized, len(self.levels))
                plt.close()
                plt.ion()

            self.gdslib, cell1 = cell_wpol_gdspy(cs, cellname, prec = self.precision, mpoints=1e9)
           

            return self.gdslib
        












    
    def _create_layout_raster_klayout(self, cellname='top', merge=True, break_vertices=250, layer_name_prefix="Level", 
        level_scale=micro):
        """
        Creates the gds layout using raster mode where each data point is a pixel rectangle to
        be defined in the layout
        Args:
            :cellname:          name of the topcell to include all the merged polygons
            :merge:             default False. If True, will merge the individual pixel polygons 
            :break_vertices:    threshold value to speed up the merging of polygons
            :layer_name_prefix: prefix of the layer name
            :level_scale:       scaling factor to apply to the level value  (default 1e-6 for micro)
        
        """

        import pya
        #initialize
        layout = pya.Layout()
        top = layout.create_cell(cellname)
        self.layout = layout
        
        # cell = gdspy.Cell(cellname,exclude_from_current=True)
        
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
        # topcell = gdspy.Cell(cellname, exclude_from_current=True)
        maskcell = layout.create_cell('pixelmask')

    
        # Make evaluation of computational effort
        if self.verbose:
            print("Mask has %d number of points distributed in %d layers"%(total_points, total_layers))       

        with Timer("Total time converting to GDS"):
            # Creates a list of cells where each cell will correspond to a layer, to enable merging afterwards.
            list_cells = []
            # for layer in self.layers:
            for layer,level in zip(self.layers, self.levels):



                # layer_name = "%s%03d"%(layer_name_prefix,layer)
                layer_name = "%s%03d_%0.3f"%(layer_name_prefix,layer, level/level_scale)
                # layer_i = layout.layer(int(layer), datatype, layer_name)
                layer_i = layout.layer(layer_name)


            # Iterates to create a layout in each layer corresponding to one discretized level
            if self.verbose:
                print("Creating individual pixel polygons")
            with Timer("Create Polygons"):
                for i in range(size_x):
                    for j in range(size_y):

                        current_point = i*size_x+j
                        aperture_value = self.aperture[i,j]
                        
                        if not np.isnan(aperture_value):
                            current_layer = int(np.rint(aperture_value))
                            # print('current layer', current_layer)
                            x = XX[i,j]
                            y = YY[i,j]
                            
                            # Creates corners of rectangle with center at the data value
                            rectangle_first_corner = (x-half_pixel_x,y-half_pixel_y)
                            rectangle_second_corner = (x+half_pixel_x,y+half_pixel_y)
                            
                            # Creates rectangle and adds it to the cell corresponding to the current layer
                            
                            maskcell.shapes(current_layer).insert(pya.DBox(rectangle_first_corner[0],rectangle_first_corner[1],
                                                                    rectangle_second_corner[0],rectangle_second_corner[1]))
                        
                    if self.verbose:
                        progress_bar(current_point/total_points)
                        
                if self.verbose: 
                    progress_bar(1)

            if merge:
                if self.verbose:
                    print("Merging polygons inside layers")

                with Timer("Merging polygons"):

                    # creates a merged cell to add the regions into
                    mergedcell = layout.create_cell('mask')


                    total_layers = len(layout.layer_infos())
                    for layer_i, layer in enumerate(layout.layer_infos()):

                      
                        if self.verbose: 
                            progress_bar(layer_i/total_layers)

                        # Creates a region of all polygons within a layer
                        region = pya.Region(maskcell.begin_shapes_rec(layer_i))
                        # merges polygons in region
                        region.merge()
                        mergedcell.shapes(layer_i).insert(region)
                    
                    # deletes original cell from 
                    maskcell.clear()
                    maskcell.delete()
                    if self.verbose: 
                        progress_bar(1)
                    trans = pya.Trans(pya.Point(0,0))
                    
                    new_instance = pya.DCellInstArray(mergedcell.cell_index(), trans, pya.Vector(50, 0 ), pya.Vector(0, 50), 1,1)
            else:
                new_instance = pya.DCellInstArray(maskcell.cell_index(), trans, pya.Vector(50, 0 ), pya.Vector(0, 50), 1,1)

            top.insert(new_instance)



            return self.layout 
        

def find_closest_indices(x, y):
    from scipy.spatial.distance import cdist

    # Reshape y to a column vector
    y = y.reshape(-1, 1)

    # Calculate the distances between each value of y and x
    distances = cdist(y, x.reshape(-1, 1))

    # Find the index of the closest value in x for each value of y
    closest_indices = np.argmin(distances, axis=1)

    return closest_indices


def levels2grayvalue(levels, calibration_height, calibration_grayvalue, ):
    idx = find_closest_indices(calibration_height, levels)

    grayvalue = calibration_grayvalue[idx]
    
    return grayvalue


def load_grayscale_contrast(filename):
    grayvalue, height = np.loadtxt(filename, delimiter=',', unpack=True)
    grayvalue = np.array(grayvalue).astype(int)
    return grayvalue, height






def create_poly_dicing_corner(length, width):
    poly = pya.DPolygon([ 
    pya.DPoint(0, 0), pya.DPoint(0, length),  pya.DPoint(width, length),
    pya.DPoint(width, width),  pya.DPoint(length, width), pya.DPoint(length, width),
    pya.DPoint(length, 0)
    ])
    return poly

def create_corners_cell(layout, field_width, field_height, corner_length, corner_width, layer="layer127"):
    cell_corner_single = layout.create_cell('corner_instance')

    cell_corner = layout.create_cell('corners')
    layer_corners = layout.layer(layer)
    corner_polygon = create_poly_dicing_corner(corner_length, corner_width)

    cell_corner_single.shapes(layer_corners).insert(corner_polygon)


    instance = pya.DCellInstArray(cell_corner_single.cell_index(),
                                pya.DTrans(pya.DTrans.R0, pya.DPoint(-field_width/2, -field_height/2))
                                )
    cell_corner.insert(instance)

    instance = pya.DCellInstArray(cell_corner_single.cell_index(),
                                pya.DTrans(pya.DTrans.R270, pya.DPoint(-field_width/2, field_height/2))
                                )
    cell_corner.insert(instance)


    instance = pya.DCellInstArray(cell_corner_single.cell_index(),
                                pya.DTrans(pya.DTrans.R180, pya.DPoint(field_width/2, field_height/2))
                                )
    cell_corner.insert(instance)

    instance = pya.DCellInstArray(cell_corner_single.cell_index(),
                                pya.DTrans(pya.DTrans.R90, pya.DPoint(field_width/2, -field_height/2))
                                )
    cell_corner.insert(instance)


    cell_corner.flatten(-1, True)

    return cell_corner



def create_label_cell(layout, text, position=(0,0), rotate=True, mag=1000, layer="layer127"):
    

    cell_label = layout.create_cell('label')
    layer_label = layout.layer(layer)


    gen = pya.TextGenerator.default_generator()
    
    region = gen.text(text, layout.dbu, mag)
    if rotate:
        cell_label.shapes(layer_label).insert(region, pya.DTrans(pya.DTrans.R90,pya.DVector(position[0], position[1])))
    else:
        cell_label.shapes(layer_label).insert(region, pya.DTrans(pya.DVector(position[0], position[1])))
    return cell_label


class GrayscaleCalibration():
    """
    Class GrayscaleCalibration:
        Creates a class to integrate the GDS library and layout corresponding to a mask aperture.
        Receives an aperture and provides methods to calculate the corresponding GDS layout
    
    Args:
        :mask: aperture object
        """
    def __init__(self, mask_height_scale=micro, verbose=True):
        self.calibration_grayvalue = None
        self.calibration_height = None
        self.calibration_height_range = None
        self.layout = None
        self.verbose = verbose

        self.mask_heights = None
        self.level_heights = None
        self.level_height_offset = None
        self.mask_height_normalization = None
        self.mask_height_scale = mask_height_scale
        self.mask_grayvalues = None
        self.mask_height_range = None


        self.layer_names_calibrated = False
        
    def load_calibration(self, filename):
        """
        Loads the calibration file from a csv file with two columns, the first column with the grayscale value and the second column with the height value
        
        Args:
            :filename: string with the path to the calibration file
        """
        self.calibration_grayvalue, self.calibration_height = load_grayscale_contrast(filename)

        self.calibration_height_range = np.ptp(self.calibration_height)
        if self.verbose == True:
            print("Grayscale calibration range is %0.3f um"%(self.calibration_height_range))



    def load_gdsfile(self, gdsfile):
        """
        Loads the mask file
        
        Args:
            :gdsfile: string with the path to the gds file
        """
        self.layout = pya.Layout()
        self.layout.read(gdsfile)
        

    def get_levels_from_layer_names(self):
        """
        Extracts the levels from the layer names in the mask
        Assumes that the layer name is "LevelDDD_xxx" where xxx is the height in um
        """

        assert self.layout is not None, "No layout loaded. Run load_gdsfile() first."


        # extracting the level and height from mask
        total_layers = len(self.layout.layer_infos())

        heights = []
        for layer_i, layer in enumerate(self.layout.layer_infos()):
            layername = str(layer.name).replace('\'\'\n', '')
            # print(layername)
            # print(layer.name)
            if "Level" not in layername:
                continue
            _, height = layername.split("_")
            height = float(height)
            heights.append(height)
            # print(layer_i, layer, layername, height)

        assert len(heights) >0, "No levels found in the mask. Check if the mask is correctly generated."

        self.mask_heights = np.array(heights)

        if self.verbose:
            print("Loaded mask with %d levels"%len(self.mask_heights))
            
        self.mask_height_range = np.ptp(self.mask_heights)
        if self.verbose == True:
            print("Mask range is %0.3f um"%(self.mask_height_range))
            if self.calibration_height_range is not None:
                if self.mask_height_range>self.calibration_height_range:
                    print("Warning! Calibration range is smaller than mask height range.")
        
        if self.calibration_grayvalue is not None:
            # Adjusts the grayvalues of the mask to the calibration values
            self.adjust_grayvalues( self.level_height_offset)



    def get_levels_from_height_range(self, mask_height_range, endpoint=True):
        """
        Extracts the levels from the given height range and number of layers on the mask

        """

        assert self.layout is not None, "No layout loaded. Run load_gdsfile() first."


        # extracting the level and height from mask
        total_layers = len(self.layout.layer_infos())

        heights = np.linspace(0, -mask_height_range, total_layers, endpoint=endpoint)
        self.mask_heights = np.array(heights)

        if self.verbose:
            print("Loaded mask with %d levels"%len(self.mask_heights))
            
        self.mask_height_range = np.ptp(self.mask_heights)
        if self.verbose == True:
            print("Mask range is %0.3f um"%(self.mask_height_range))
            if self.calibration_height_range is not None:
                if self.mask_height_range>self.calibration_height_range:
                    print("Warning! Calibration range is smaller than mask height range.")
        
        if self.calibration_grayvalue is not None:
            # Adjusts the grayvalues of the mask to the calibration values
            self.adjust_grayvalues( self.level_height_offset)



    def plot_calibration(self):
        assert self.calibration_grayvalue is not None, "No calibration data loaded."

        plt.figure()
        plt.plot(self.calibration_grayvalue, self.calibration_height)

        if self.level_heights is not None:
            for level in self.level_heights:
                plt.axhline(level, color='Gray', linestyle='-', lw=1)
        plt.scatter(self.mask_grayvalues, self.level_heights, s=4)
        plt.xlabel("Grayvalue")
        plt.ylabel("Height [um]")

    def adjust_grayvalues(self, height_offset=None):
        """
        Adjusts the grayvalues of the mask to the calibration values
        """
        self.mask_height_normalization = np.max(self.mask_heights)-np.max(self.calibration_height) 
        if height_offset is None:
            height_offset = -self.mask_height_normalization
        else:
            height_offset = -self.mask_height_normalization+height_offset
        self.level_height_offset = height_offset
        self.level_heights = self.mask_heights+self.level_height_offset

        self.mask_grayvalues = levels2grayvalue(self.level_heights, 
                                                self.calibration_height, self.calibration_grayvalue)





    def calibrate_layer_names(self, force_naming=False, datatype=0):
        """
        Changes the layer names to the calibrated grayvalues
        """

        assert self.layer_names_calibrated == False, "Layer names already changed."
        self.layer_names_calibrated = True

        for layer_i, layer in enumerate(self.layout.layer_infos()):
            layername = str(layer)
            # print(layername)
            if not force_naming:
                if "Level" not in layername:
                    continue
                prefix, height = layername.split("_")
                _,idx = prefix.split("Level")
                idx = int(idx)
                ly = self.layout.layer(layername)

            else:
                idx = layer_i
                ly = self.layout.layer(idx, datatype )

            info = self.layout.get_info(ly)
            info.name = "layer%03d"%self.mask_grayvalues[idx]
            newinfo = pya.LayerInfo(info.name)

            # print(info)
            self.layout.set_info(ly, newinfo)
            # print("Layer handle " + str(ly) + " refers to " + str(self.layout.get_info(ly)))


    
    def save_calibrated_gdsfile(self, filename, force_save=False):

        if not force_save:
            assert self.layer_names_calibrated == True, "Must calibrate layer names first"
        

        with Timer("Saving calibrated GDS file"):
            self.layout.write(filename)
        print("Saved %s"%(filename))

    def add_corners(self, field_width=10000, field_height=10000, corner_length=500, corner_width=100, layer="layer127"):
        """ Adds corners around the mask
        """
        trans = pya.Trans(pya.Point(0,0))

        top_cell = self.layout.top_cell()

        cell_corner = create_corners_cell(self.layout, field_width, field_height, corner_length, corner_width, layer)
        top_cell.insert(pya.DCellInstArray(cell_corner.cell_index(), trans))

    def add_label(self, label, position=(-4000, -4000),rotate=True, mag=100, layer="layer127"):
        """ Adds text label to the mask file
        """
        top_cell = self.layout.top_cell()
        trans = pya.Trans(pya.Point(0,0))

        cell_label = create_label_cell(self.layout, label, position,rotate, mag, layer)
        top_cell.insert(pya.DCellInstArray(cell_label.cell_index(), trans))





# gdsfile = "testmask4.dxf"   #name of gds file 

# filename = "Pyra2_E80_F-20_greylevels.csv"
# calibration_grayvalue, calibration_height = load_grayscale_contrast(filename)
# plt.plot(calibration_grayvalue, calibration_height)

# import klayout
# import pya


# # Loads the mask file
# layout = pya.Layout()
# layout.read(gdsfile)



# # extracting the level and height from mask
# total_layers = len(layout.layer_infos())

# heights = []
# for layer_i, layer in enumerate(layout.layer_infos()):
#     layername = str(layer)
#     if "Level" not in layername:
#         continue
#     _, height = layername.split("_")
#     height = float(height)
#     heights.append(height)
#     # print(layer_i, layer, layername, height)

# assert len(heights) >0, "No levels found in the mask. Check if the mask is correctly generated."

# heights = np.array(heights)
    
# plt.plot(calibration_grayvalue, calibration_height)
# height_offset = -1.4
# grayvalues = levels2grayvalue(heights, calibration_height, calibration_grayvalue, height_offset=height_offset)

# for level in heights:
#     plt.axhline(level+height_offset, color='red', linestyle='-', lw=1)
# # plt.plot(mask1.levels, grayvalue)
# grayvalues, heights



# for layer_i, layer in enumerate(layout.layer_infos()):
#     layername = str(layer)
    
#     if "Level" not in layername:
#         continue
#     prefix, height = layername.split("_")
#     _,idx = prefix.split("Level")
#     idx = int(idx)
#     print(idx)

#     ly = layout.layer(layername)
#     info = layout.get_info(ly)
#     print(info)
#     info.name = "layer%03d"%grayvalues[idx]
#     # print(info)
#     layout.set_info(ly, info)
#     print("Layer handle " + str(ly) + " refers to " + str(layout.get_info(ly)))

# outfile = gdsfile.replace(".dxf", "_grayscale.dxf")
# layout.write(outfile)