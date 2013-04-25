## Python Colony Analyzer Toolkit
# Colony size module
# Gordon Bean, April 2013

# Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import stats as pystat
import sys

sys.path.append('/cellar/users/gbean/toolkits/python_colony_analyzer/lib')
from load_and_crop import *
from gridtools import *
from intensitythreshold import *

sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
from bean import *

## Methods
def find_colony_borders( box, it=None ):
    midr = box.shape[0]/2
    midc = box.shape[1]/2
    
    w = midr / 3
    
    # North
    mi = np.argmin(box[:midr, midc-w:midc+w], 0)
    rmin = np.median(mi)
    if it is not None: # get closer to the colony
        tmp = np.max(box[midr:rmin:-1, midc-w:midc+w], 1) < it
        rmin = midr - (find(tmp,0) or len(tmp))
    
    # South
    mi = np.argmin(box[midr:, midc-w:midc+w], 0)
    rmax = np.median(mi)+midr
    if it is not None:
        tmp = np.max(box[midr:rmax, midc-w:midc+w], 1) < it
        rmax = midr + (find(tmp,0) or len(tmp))
        
    # West
    mi = np.argmin(box[midr-w:midr+w, :midc],1)
    cmin = np.median(mi)
    if it is not None:
        tmp = np.max(box[midr-w:midr+w, midc:cmin:-1], 0)
        cmin = midc - (find(tmp,0) or len(tmp))
        
    # East
    mi = np.argmin(box[midr-w:midr+w, midc:],1)
    cmax = np.median(mi) + midc
    if it is not None:
        tmp = np.max(box[midr-w:midr+w, midc:cmax], 0) < it
        cmax = midc + (find(tmp,0) or len(tmp))
        
    return (rmin, rmax, cmin, cmax)
    
def threshold_bounded( plate, grid, r, c ):
    box = get_box( plate, grid.r[r,c], grid.c[r,c], grid.win )
    rmin, rmax, cmin, cmax = find_colony_borders( box, grid.thresh[r,c] )
    
    box = box > grid.thresh[r,c]
    box[:rmin,:] = False
    box[rmax:,:] = False
    box[:,:cmin] = False
    box[:,cmax:] = False
    
    return np.sum(box), box
    
def measure_colony_sizes( plate, **params ):
    default_param( params, \
        manualgrid = False, \
        thresholdobject = local_fitted(), \
        localthreshold = True, \
        sizefunction = threshold_bounded )

    # TODO
    # - check for plate type and dimensions
    # - manual grid
    
    # Load plate
    if isinstance(plate, basestring):
        # "plate" is the file name
        plate = load_plate(plate, **params)
    
    # Determine grid
    if params['manualgrid']:
        print "Warning, manual grid not yet implemented\n"
        return None, None
    else:
        if 'grid' not in params:
            params['grid'] = determine_colony_grid( plate, **params )
    
    grid = params['grid']
    
    # Intensity thresholds
    if not hasattr(grid, 'thresh'):
        if params['localthreshold']:
            grid.thresh = compute_local_thresholds( plate, grid, **params )
        else:
            grid.thresh = compute_global_threshold( plate, grid, **params )
    
    # Measure colony size
    size_fun = params['sizefunction']
    sizes = bean.nans(grid.dims)
    grid.threshed = np.zeros(plate.shape)==1
    for r in range(0, grid.dims[0]):
        for c in range(0, grid.dims[1]):
            sizes[r,c], tmp = size_fun( plate, grid, r, c )
            set_box(grid.threshed, tmp, grid.r[r,c], grid.c[r,c])
    
    return sizes, grid
    
    
    
