## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.analysis Module

# Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as si

# Bean's bag-o-tricks
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

## Functions
def pseudoplate( colsizes, **params ):
    bean.default_param( params, 
        colorbar = True, 
        emap = False, 
        style = 'circles', 
        notemarker = '\leftarrow' )
    
    # Find dimensions
    n = bean.numel(colsizes)
    dims = (np.array([8, 12]) * np.sqrt( n / 96 )).astype('int')
    
    if 'max_circ' not in params:
        rad = np.min(np.diff( plt.gca().bbox.get_points(), axis=0) / 1.5 / dims)
        params['max_circ'] = np.pi * rad**2
    
    # Reorder plate, if necessary
    # TODO - port this functionality
    
    # Draw figure
    if params['style'] == 'imagesc':
        bean.imagesc( np.reshape(colsizes, dims) )
        
    elif params['style'] == 'circles':
        xx, yy = np.meshgrid( range(0,dims[1]), range(0,dims[0]) )
        
        max_circ = params['max_circ']
        tmp = np.reshape( colsizes, dims )
        
        # Calculate size of spots
        if 'min_size' in params:
            min_size = params['min_size']
        else:
            min_size = np.min(tmp)
        
        tmp_size = tmp - min_size + 0.0
        tmp_size[tmp_size < 0] = 0.0
        
        if 'max_size' in params:
            max_size = params['max_size'] - min_size
        else:
            max_size = np.max(tmp_size)
        
        if max_size != 0:
            tmp_size = tmp_size / (max_size + 0.0)
            tmp_size[tmp_size > 1] = 1
            tmp_size = (max_circ - 1) * tmp_size
        else:
            tmp_size = np.ones( dims ) * (max_circ - 1)
                
        plt.scatter( xx, yy, s=tmp_size + 1, c=tmp, edgecolors='none' )
        plt.xlim(-1, dims[1])
        plt.ylim(-1, dims[0])
        plt.gca().invert_yaxis()
        plt.draw()
        
    else:
        pass

##---------------------##
##    Colony Filters   ##
##---------------------##
def h2_colony_filter( plate ):
    plate = np.array(plate).astype(float)
    shape = plate.shape / np.array([8, 12])
    rows = np.array([7, 8]) * shape
    cols = np.array([1, 2]) * shape
    plate[rows[0]:rows[1], cols[0]:cols[1]] = np.nan
    return plate
    
def empty_spot_filter( plate, emptythreshold=50 ):
    plate = np.array(plate).astype(float)
    plate[plate < emptythreshold] = np.nan
    return plate
    
##---------------------##
## Spatial Corrections ##
##---------------------##
def empty_neighbor( plate, **params ):
    bean.default_param( params, 
        emptyspots=np.isnan(plate) )
    
    footprint = np.ones((3,3))
    footprint[0,0], footprint[0,-1], footprint[-1,0], footprint[-1,-1] = 0,0,0,0
    ens = si.maximum_filter( params['emptyspots'], footprint=footprint )
    nnans = ~np.isnan(plate)
    ensm = bean.pmode( plate[np.logical_and(ens, nnans)] )
    nensm = bean.pmode( plate[np.logical_and(~ens, nnans)] )
    
    pmed = np.median(plate)
    filt = np.ones(plate.shape)
    filt[ens] = ensm / nensm
    
    return np.array(plate / filt), filt * pmed
    
def simple_border( plate, **params ):
    bean.default_param( params, 
        depth=2 )
    
    # Compute row/column medians
    d = params['depth']
    pmed = np.median(plate)
    plate_ = bean.fil( plate, np.isnan, pmed )
    filt = np.array([[max(rm,cm) for cm in np.median(plate_,0)] \
        for rm in np.median(plate_,1)])
    filt[d:-d,d:-d] = pmed
    
    # Apply correction
    return np.array( plate / filt * pmed ), filt

def median_spatial( plate, **params ):
    bean.default_param( params, 
        size=(7,7), 
        mode='reflect' )
    
    # Compute filter
    pmed = np.median(plate)
    plate_ = bean.fil( plate, np.isnan, pmed )
    filt = si.median_filter( plate_, size=params['size'], mode=params['mode'] )
    
    # Scale
    return np.array( plate / filt * pmed ), filt
    
def apply_spatial_corrections( plate, **params ):
    """ Perform a spatial correction on the colony size data. """
    bean.default_param( params, 
        corrections = [median_spatial, simple_border], 
        internalx4 = False, 
        returninfo = False )
        
    ## Reshape plate
    shape = plate.shape
    n = bean.numel(plate)
    dims = (np.array([8, 12]) * np.sqrt( n / 96 )).astype('int')
    plate = np.reshape( plate, dims )
    
    ## Apply corrections
    outs = []
    for correction in params['corrections']:
        plate, out = correction( plate, **params )
        outs.append(np.reshape(out,shape))
    
    # Return results
    if params['returninfo']:
        return np.reshape(plate,shape), outs
    else:
        return np.reshape(plate,shape)
    
## Normalization ##
def normalize_colony_sizes( plate, **params ):
    plate = np.array(plate)
    return plate / bean.pmode(bean.ind(plate, bean.nisnan))
    
    
    
    
    
    
    


