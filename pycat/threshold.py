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
        
class max_min_mean_threshold( threshold_object ):
    def get_colony_box(self, plate, grid, row, col):
        return pygrid.get_box( plate, \
            grid.r[row,col], grid.c[row,col], grid.win )
    def determine_threshold(self, plate, grid=None, row=None, col=None):
        if grid is None:
            box = plate
        else:
            box = self.get_colony_box(plate, grid, row, col)
        return (np.max(box) + np.min(box)) / 2
        
class local_fitted_threshold( threshold_object ):
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
    
class local_gaussian_threshold(threshold_object):
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
def local_gaussian(plate, grid=None, sigma=None, mode='reflect', offset=0):
    if sigma is None:
        if grid is not None:
            sigma = grid.win/2
        else:
            sigma = np.min(plate.shape)/4
        
    gplate = np.zeros(plate.shape, 'double')
    scipy.ndimage.gaussian_filter(plate, sigma, output=gplate, mode=mode)
    return plate > (gplate - offset)
    
def local_fitted( plate, grid ):
    bplate = np.zeros(plate.shape) == 1
    locfit = local_fitted_threshold()
    for r in range(0, grid.dims[0]):
        for c in range(0, grid.dims[1]):
            box = locfit.get_colony_box( plate, grid, r, c )
            it = locfit.determine_threshold( box )
            pygrid.set_box( bplate, box > it, r, c )
    return bplate
            
