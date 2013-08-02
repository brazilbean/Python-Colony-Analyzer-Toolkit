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
class IntensityThreshold( object ):
    def get_colony_box( self, plate, grid, row, col ):
        return pygrid.get_box( plate, \
            grid.r[row,col], grid.c[row,col], grid.win )
            
    def determine_threshold( self, plate, grid, rr, cc ):
        return None

    def apply_threshold( self, plate, grid ):
        # Iterate through each position in grid
        # - Get the box
        # - Estimate the threshold
        # - save the thresholded box
        thrplate = np.zeros(plate.shape) == 1 
        for r in range(0, grid.dims[0]):
            for c in range(0, grid.dims[1]):
                box = self.get_colony_box(plate, grid, r, c)
                it = self.determine_threshold( box )
                pygrid.set_box(thrplate, box>it, grid.r[r,c], grid.c[r,c])
        return thrplate
                
class MaxMinMeanThreshold( IntensityThreshold ):
    def get_colony_box(self, plate, grid, row, col):
        return pygrid.get_box( plate, \
            grid.r[row,col], grid.c[row,col], grid.win )
    def determine_threshold(self, plate, grid=None, row=None, col=None):
        if grid is None:
            box = plate
        else:
            box = self.get_colony_box(plate, grid, row, col)
        return (np.max(box) + np.min(box)) / 2
        
class LocalFittedThreshold( IntensityThreshold ):
    bins = None
    fdr = None
    fast = None
    fullplate = None
    num_background_iters = None
    upper_threshold_function = None
    pthresh = 4.5;
    
    def __init__(self, fdr = 0.01, fast = True, full = False, \
                 pthresh=4.5, num_background_iters = 5, \
                 upper_threshold_function = None):
        self.pthresh = pthresh
        self.fdr = fdr
        self.fast = fast
        self.full = full
        self.num_background_iters = num_background_iters
        self.upper_threshold_function = upper_threshold_function
        
        # Full plate settings
        if self.full and self.upper_threshold_function is not None:
            self.upper_threshold_function = lambda box: \
                (np.min(box) + np.max(box))/2
        
    def _get_pm_std( self, box ):
        if self.fast:
            pmfun = pygrid.pixelmode
        else:
            pmfun = bean.pmode
        pm = pmfun( box )
        
        if self.upper_threshold_function is not None:
            it = self.upper_threshold_function( box )
            c = self.num_background_iters
            while c > 0 and pm > it:
                it = (it + np.min(box))/2
                pm = pmfun(box[box < it])
                c = c - 1
        
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
    
class LocalGaussianThreshold( IntensityThreshold ):
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
        
    def apply_threshold( self, plate, grid ):
        if self.gplate is None:
            self._set_gplate(plate, grid)

        return self.gplate
    
class BackgroundOffset( IntensityThreshold ):
    offset = None
    fullplate = None
    background_max = None
    
    def __init__(self, offset = 1.25, fullplate = False, background_max = None):
        self.offset = offset
        self.fullplate = fullplate
        self.background_max = background_max
        
    def calibrate(self, plate, grid):
        mid = pygrid.get_box \
            (plate, np.mean(grid.r), np.mean(grid.c), grid.win*5)
        self.background_max = (np.min(mid) + np.max(mid)) / 2
        
    def determine_threshold(self, box):
        if self.fullplate:
            bg = pygrid.pixelmode(box[box < self.background_max])
        else:
            bg = pygrid.pixelmode(box)
            
        it = bg * self.offset
        
    def apply_threshold(self, plate, grid):
        if self.fullplate and self.background_max is None:
            self.calibrate(plate, grid)
        return super(BackgroundOffset, self).apply_threshold(plate, grid)
        
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
            
