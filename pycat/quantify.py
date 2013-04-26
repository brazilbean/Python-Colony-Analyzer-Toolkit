## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.threshold Module

# Imports
import numpy as np

# Bean's bag-o-tricks
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

# pycat imports
import pycat.grid as pygrid
import pycat.io as io
import pycat.threshold as pythresh

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
    
def measure_colony_sizes( plate, **params ):
    bean.default_param( params, \
        manualgrid = False, \
        thresholdobject = pythresh.local_fitted(), \
        localthreshold = True, \
        sizefunction = threshold_bounded )

    # TODO
    # - check for plate type and dimensions
    # - manual grid
    
    # Load plate
    if isinstance(plate, basestring):
        # "plate" is the file name
        plate = io.load_plate(plate, **params)
    
    # Determine grid
    if params['manualgrid']:
        print "Warning, manual grid not yet implemented\n"
        return None, None
    else:
        if 'grid' not in params:
            params['grid'] = pygrid.determine_colony_grid( plate, **params )
    
    grid = params['grid']
    
    # Intensity thresholds
    if not hasattr(grid, 'thresh'):
        grid.thresh = pythresh.local_gaussian( plate, grid )
        
        #if params['localthreshold']:
        #    grid.thresh = \
        #        pythresh.compute_local_thresholds( plate, grid, **params )
        #else:
        #    grid.thresh = \
        #        pythresh.compute_global_threshold( plate, grid, **params )
    
    # Measure colony size
    size_fun = params['sizefunction']
    #sizes = bean.nans(grid.dims)
    #grid.threshed = np.zeros(plate.shape)==1
    #for r in range(0, grid.dims[0]):
    #    for c in range(0, grid.dims[1]):
    #        sizes[r,c] = size_fun( plate, grid, r, c )
    #        #pygrid.set_box(grid.threshed, tmp, grid.r[r,c], grid.c[r,c])
    #        
    rrr, ccc = range(0, grid.dims[0]), range(0, grid.dims[1])
    sizes = [ [ size_fun( plate, grid, r, c ) for c in ccc] for r in rrr];
    return np.array(sizes), grid
    
 