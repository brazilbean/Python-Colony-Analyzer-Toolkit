## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.threshold Module

# Imports
import numpy as np

# Bean's bag-o-tricks
import bean # https://github.com/brazilbean/bean-python-toolkit

# pycat imports
import grid as pygrid
#import pycat.io as io
#import pycat.threshold as pythresh

## Methods
def find_colony_borders( box, bbox=None ):
    """ Finds the rectangular bounding box around the center colony. """
    
    midr = box.shape[0]/2
    midc = box.shape[1]/2
    
    w = midr / 3
    
    # North
    mi = np.argmin(box[:midr, midc-w:midc+w], 0)
    rmin = np.median(mi)
    if bbox is not None: # get closer to the colony
        tmp = np.any(bbox[midr:rmin:-1, midc-w:midc+w], 1)
        rmin = midr - (bean.find(tmp,0) or len(tmp))
    
    # South
    mi = np.argmin(box[midr:, midc-w:midc+w], 0)
    rmax = np.median(mi)+midr
    if bbox is not None:
        tmp = np.any(bbox[midr:rmax, midc-w:midc+w], 1)
        rmax = midr + (bean.find(tmp,0) or len(tmp))
        
    # West
    mi = np.argmin(box[midr-w:midr+w, :midc],1)
    cmin = np.median(mi)
    if bbox is not None:
        tmp = np.any(bbox[midr-w:midr+w, midc:cmin:-1], 0)
        cmin = midc - (bean.find(tmp,0) or len(tmp))
        
    # East
    mi = np.argmin(box[midr-w:midr+w, midc:],1)
    cmax = np.median(mi) + midc
    if bbox is not None:
        tmp = np.any(bbox[midr-w:midr+w, midc:cmax], 0)
        cmax = midc + (bean.find(tmp,0) or len(tmp))
        
    return (rmin, rmax, cmin, cmax)
    
def threshold_bounded( plate, grid, r, c ):
    """ Quatifies the colony size using a bouding box on the original image
        and the pixels defined by the binary plate image. """
    
    box = pygrid.get_box( plate, grid.r[r,c], grid.c[r,c], grid.win )
    bbox = pygrid.get_box( grid.thresh, grid.r[r,c], grid.c[r,c], grid.win )
    rmin, rmax, cmin, cmax = find_colony_borders( box, bbox )
    
    bbox[:rmin,:] = False
    bbox[rmax:,:] = False
    bbox[:,:cmin] = False
    bbox[:,cmax:] = False
    
    return np.sum(bbox)
    
def clear_adjacent_colonies(box, bbox):
    # Find colony borders
    rmin, rmax, cmin, cmax = find_colony_borders( box, bbox )
    
    # Block out everything but colony in binary image
    bbox[:rmin,:] = False
    bbox[rmax:,:] = False
    bbox[:,:cmin] = False
    bbox[:,cmax:] = False
    
    # Block out everything but colony in pixel image
    box[bbox] = np.nan
    
    return box, bbox
    
## Classes
class ColonyQuantifier(object):
    label = None
    
    def __init__(self, label = 'size'):
        self.label = label
    
    def __call__(self, *args, **kwargs):
        return self.quantify(*args, **kwargs)
        
    def parse_box(self, *args):
        if len(args) == 4:
            # Plate, grid, and row/column indexes passed
            plate, grid, row, col = args
            box = pygrid.gridtools.get_box( 
                plate, grid.r[row,col], grid.c[row,col], grid.win )
            
            # Determine binary image
            bbox = pygrid.gridtools.get_box( 
               grid.thresh, grid.r[row,col], grid.c[row,col], grid.win )
           
        elif len(args) == 2:
            # Pixel and binary colony boxes passed
            box, bbox = args
            
        elif len(args) == 1:
            # Binary colony box passed
            box = []
            bbox = args[0]
            
        else:
            raise Exception('Incorrect number of arguments: %i'%(len(args)))
            
        return box, bbox
            
class ColonyArea(ColonyQuantifier):
    def __init__(self):
        super(ColonyArea, self).__init__('area')
    
    def quantify(self, *args):
        box, bbox = self.parse_box(*args)
        if len(box) == 0:
            raise Exception('This method requires pixel information.')
        
        
        # Remove adjacent colonies
        _, bbox = clear_adjacent_colonies(box, bbox)
        
        # Sum area
        return np.sum(bbox)
        
class ColonySumIntensity(ColonyQuantifier):
    def __init__(self):
        super(ColonySumIntensity, self).__init__('sumintensity')
        
    def quantify(self, *args):
        box, bbox = self.parse_box(*args)
        if len(box) == 0:
            raise Exception('This method requires pixel information.')
        
        
        # Remove adjacent colonies
        box, bbox = clear_adjacent_colonies(box, bbox)
        
        # Sum pixel intensity
        return np.nansum(box)
            
            
            
    
