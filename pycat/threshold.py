## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.threshold Module

# Imports
import numpy as np
import scipy
from scipy import stats as pystat

# Bean's bag-o-tricks
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

# pycat imports
import pycat.grid as pygrid

## Classes
class threshold_object:
    def get_colony_box( self, plate, grid, rr, cc ):
        return None
    def determine_threshold( self, plate, grid, rr, cc ):
        return None
        
class max_min_mean( threshold_object ):
    def get_colony_box(self, plate, grid, row, col):
        return pygrid.get_box( plate, \
            grid.r[row,col], grid.c[row,col], grid.win )
    def determine_threshold(self, plate, grid=None, row=None, col=None):
        if grid is None:
            box = plate
        else:
            box = self.get_colony_box(plate, grid, row, col)
        return (np.max(box) + np.min(box)) / 2
        
class local_fitted( threshold_object ):
    pthresh = 4.5;
    
    def __init__(self, pthresh=4.5):
        self.pthresh = pthresh
        
    def _get_pm_std( self, box ):
        mb = np.min(box)
        it = (np.max(box) + mb) / 2
        #pm = bean.pmode(box)
        pm = pygrid.pixelmode(box)
        c = 5
        while c > 0 and pm > it:
            it = (it + mb) / 2
            #pm = bean.pmode(box[box<it])
            pm = pygrid.pixelmode(box[box<it])
            c -= 1
        tmp = box[box<pm] - pm
        st = np.std(np.vstack((tmp,-tmp)))
        return pm, st
        
    def _get_pvals( self, b, bpm, st ):
        return -np.log10( 1 - pystat.norm.cdf( b-bpm, 0, st ) )
        
    def get_colony_box( self, plate, grid, rr, cc):
        # Determine window size
        win = grid.win
        if np.prod(grid.dims) == 6144:
            win = win * 2
        elif np.prod(grid.dims) == 24576:
            win = win * 4
        
        # Get box
        return pygrid.get_box( plate, grid.r[rr,cc], grid.c[rr,cc], win )
        
    def determine_threshold( self, plate, grid=None, row=None, col=None ):
        if grid is None:
            box = plate
        else:
            box = self.get_colony_box(plate, grid, row, col)
        
        # Fit background pixels
        bpm, st = self._get_pm_std( box )
        
        # Compute p-values
        pvals = self._get_pvals( box, bpm, st )
        
        # Threshold
        bb = pvals > self.pthresh
        return np.min( box[bb] ) - 1
    
class _local_gaussian(threshold_object):
    gplate = None
    mode = None
    sigma = None

    def __init__(self):
        self.gplate = None
        self.mode = 'reflect'
        self.sigma = None
        
    def _set_gplate(self, plate, grid):
        self.sigma = grid.win/2
        self.gplate = np.zeros(plate.shape, 'double')
        scipy.ndimage.gaussian_filter(plate, self.sigma, 
            output=self.gplate, mode=self.mode)
            
    def get_colony_box( self, plate, grid, row, col ):
        return pygrid.get_box( plate, grid, row, col, grid.win )
        
    def determine_threshold( self, plate, grid=None, row=None, col=None ):
        if self.gplate is None:
            self._set_gplate(plate, grid)
        
        box = self.get_colony_box(plate, grid, row, col)
        gbox = self.get_colony_box(self.gplate, grid, row, col)

        return np.min(box[box>gbox])
    
## Methods
def local_gaussian(plate, grid, sigma=None, mode='reflect', offset=0):
    if sigma is None:
        sigma = grid.win/2
        
    gplate = np.zeros(plate.shape, 'double')
    scipy.ndimage.gaussian_filter(plate, sigma, output=gplate, mode=mode)
    return plate > (gplate - offset)
    
def compute_global_threshold( plate, grid_, **params ):
    grid = grid_
    bean.default_param( params, thresholdobject=local_fitted() )
    
    if not isinstance(params['thresholdobject'], threshold_object):
        raise Exception('thresholdobject must inherit from threshold_object');
    
    # Get global box
    box = plate[np.round(np.max(grid.r[0,:])) : np.round(np.min(grid.r[-1,:])),
        np.round(np.max(grid.c[:,0])) : np.round(np.min(grid.c[:,-1])) ]
    
    # Get threshold
    filt = np.random.rand(*box.shape) < 0.05
    return np.zeros(grid.dims) \
        + params['thresholdobject'].determine_threshold(box[filt])
    
def compute_local_thresholds( plate, grid_, **params ):
    grid = grid_
    bean.default_param( params, \
        thresholdobject = local_fitted(), \
        smoothing = True, \
        parallel = False )
    
    if not isinstance(params['thresholdobject'], threshold_object):
        raise Exception('thresholdobject must inherit from threshold_object');
    
    # Iterate over grid positions
    # - get the box
    # - estimate the intensity
    thobj = params['thresholdobject'];
    if params['parallel']:
        # Parallel processing using IPython DirectView
        raise Exception('Parallel mode not yet supported')
        
    else:
        det_thresh = lambda r, c: thobj.determine_threshold( \
            thobj.get_colony_box(plate, grid, r, c) )
        rrr, ccc = range(0, grid.dims[0]), range(0, grid.dims[1])
        its = np.array([ [ det_thresh(r,c) for c in ccc ] for r in rrr ])
        
        #its = bean.nans(grid.dims)
        #for rr in range(0, grid.dims[0]):
        #    for cc in range(0, grid.dims[1]):
        #        box = thobj.get_colony_box(plate, grid, rr, cc)
        #        its[rr,cc] = thobj.determine_threshold(box)
            
    # Smoothing
    if params['smoothing']:
        # TODO - write the spatial correction function
        print "Smoothing not yet implemented...\n"
    
    return its