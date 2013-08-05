## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.grid Module

# Imports
import numpy as np
import sys
import copy

# Bean's bag-o-tricks
import bean # https://github.com/brazilbean/bean-python-toolkit

# pycat imports
import pycat.grid.gridtools as gtools
import pycat.io as io

## Classes
class Grid:
    # Properties
    info = None
    dims = None
    win = None
    r = None
    c = None
    thresh = None
    threshed = None
    
    def __init__ (self):
        self.info = {}
        self.dims = None
        self.win = None
        self.r = None
        self.c = None
        self.thresh = None
        self.threshed = None
    
    def __repr__(self):
        foo = vars(self)
        s = '';
        for k,v in foo.iteritems():
            if isinstance(v, np.ndarray) and bean.numel(v)>5:
                s += '{}\t{},{}\n'.format(k, type(v), v.shape)
            else:
                s += '{}:\t{}\n'.format(k, v)
        return s

    def copy(self):
        return copy.deepcopy(self)
        
        #other = Grid()
        #for k,v in vars(self).iteritems():
        #    if isinstance(v, np.ndarray):
        #        setattr(other, k, np.array(v))
        #    else:
        #        setattr(other, k, v)
        #return other
        
## Grid Initialization Methods
class GridMethod(object):
    gridspacing = None
    dimensions = None
    
    def __init__(self, gridspacing = None, dimensions = None, **kwargs):
        self.gridspacing = gridspacing
        self.dimensions = dimensions
        
    def __call__(self, *args, **kwargs):
        return self.fit_grid(*args, **kwargs)
        
    def fit_grid(self, plate):
        pass
        
    def initialize_grid(self, plate):
        grid = Grid()
        
        # Grid spacing
        if self.gridspacing is not None:
            grid.win = self.gridspacing
        else:
            grid.win = gtools.estimate_grid_spacing(plate)
            
        # Grid dimensions
        if self.dimensions is not None:
            grid.dims = self.dimensions
        else:
            grid.dims = gtools.estimate_dimensions(plate, grid.win)
        
        # Other elements
        grid.r = np.zeros(grid.dims) * np.nan
        grid.c = np.zeros(grid.dims) * np.nan
        
        return grid
        

class GridOffsetMethod(GridMethod):
    _size_standard = None
    _init_grid_dims = None
    
    def __init__(self, **kwargs):
        GridMethod.__init__(self)
        self._size_standard = np.array([1853.0, 2765.0])
        self._init_grid_dims = np.array([8, 8])
        
    def fit_grid(self, plate):
        # Initialize grid 
        grid = self.initialize_grid(plate)
        
        # Identify grid orientation
        grid = self.estimate_grid_orientation(plate, grid)
        
        # Compute initial placement
        grid = self.compute_initial_placement(plate, grid)
        
        # Perform initial adjustment
        grid = self.perform_initial_adjustment(plate, grid)
        
        # Determine grid-colony overlap
        overlap = self.determine_overlap(plate, grid)
        
        # Correct grid offset
        grid = self.correct_offset_grid(grid, overlap)
        
        # Make final adjustments
        grid = self.make_final_adjustments(plate, grid)
        
        # Sign the package
        grid.info['grid'] = self

    def estimate_grid_orientation(self, plate, grid):
        tang = self._size_standard[0] / self._size_standard[1]
        ratiofun = lambda xp, yp: np.arctan( -(yp-xp*tang)/(yp*tang-xp) )
        yp, xp = plate.shape
        
        theta = ratiofun( xp, yp )
        if (np.mean(plate[0,np.floor(xp/2):]) \
            > np.mean(plate[0,0:np.floor(xp/2)])):
            # The left side of the top border is brighter than the right side
            #  so the plate is rotated right instead of left
            theta = -theta
            
        return theta

    def compute_initial_placement(self, plate, grid):
        # Drop into plate - get initial position
        rpos, cpos = np.array(plate.shape)/2
        r0, c0 = gtools.adjust_spot( plate, rpos, cpos, grid.win )
        
        # Determine the initial grid positions
        colpositions = range(0, self._init_grid_dims[1]-1) * grid.win
        rowpositions = range(0, self._init_grid_dims[0]-1) * grid.win
        cc0, rr0 = np.meshgrid(colpositions, rowpositions)
        
        # Define the initial grid coordinates (top-left corner of grid)
        ri = range(0, self._init_grid_dims[0]-1)
        ci = range(0, self._init_grid_dims[1]-1)
        grid.r[ri,ci] = rr0
        grid.c[ri,ci] = cc0
        
        # Rotate the grid according to the orientation estimate
        theta = self.estimate_grid_orientation(plate, grid)
        rotmat = np.array([[np.cos(theta), -np.sin(theta)], \
                           [np.sin(theta),  np.cos(theta)]])
        val = ~np.isnan(grid.r)
        tmp = np.dot(rotmat, np.vstack((grid.c[val], grid.r[val])))
        
        # Assign position coordinates
        grid.r[val] = r0 + tmp[1,:]
        grid.c[val] = c0 + tmp[0,:]
        
        return grid

    def perform_initial_adjustment( self, plate, grid ):
        # Extrapolate grid
        rrr, ccc = np.nonzero(~np.isnan(grid.r))
        grid = gtools.adjust_grid_polar( plate, grid, rrr, ccc )
        
        # Fit diagonal of grid (to get a more robust fit)
        cc, rr = np.meshgrid(range(0,grid.dims[1]), range(0,grid.dims[0]))
        foo = np.round(grid.dims[0]/grid.dims[1] * cc) == rr
        rrr, ccc = np.nonzero(foo)
        iii = range(0, grid.dims[1])
        grid = gtools.adjust_grid_polar( plate, grid, rrr[iii], ccc[iii] )
        
        return grid

# Estimate initial grid
def estimate_initial_grid_offset( plate, **params ):
    grid = initialize_grid(plate, **params)
    
    # Drop into plate - get initial position
    rpos, cpos = np.array(plate.shape)/2
    rpos, cpos = adjust_spot( plate, rpos, cpos, grid.win )
    
    box = get_box(plate, rpos, cpos, grid.win*5)
    
    # Get initial placement
    it = estimate_intensity_threshold(plate)
    cents = bean.centroids(box>it)
    cents = cents[ ~np.any(np.isnan(cents),1), :]
    
    tmp = np.round(cents/grid.win)
    nix = np.any(np.logical_or(tmp < 1,tmp==np.max(tmp,0)),1)
    rrr, ccc = tmp[~nix,0].astype(int), tmp[~nix,1].astype(int)
    
    grid.r, grid.c = np.zeros(grid.dims)+np.nan, np.zeros(grid.dims)+np.nan
    grid.r[rrr,ccc] = cents[~nix,0] + rpos - grid.win*5
    grid.c[rrr,ccc] = cents[~nix,1] + cpos - grid.win*5
    
    grid = adjust_subgrid( plate, grid, rrr, ccc )
    
    # Do major offset
    val = np.logical_and( grid.r < plate.shape[0], grid.c < plate.shape[1] )
    grr = np.floor(grid.r[val]).astype(int)
    gcr = np.floor(grid.c[val]).astype(int)
    pdata = np.zeros(grid.dims)
    for ii,jj in enumerate(bean.find(bean.ind(val))):
        r,c = bean.ind2sub(jj, grid.dims)
        pdata[r,c] = plate[grr[ii],gcr[ii]]
    
    radj = grid.dims[0] - \
        bean.find(sum(pdata>it,1) < plate.shape[1]/grid.win-grid.dims[1],0)
    cadj = grid.dims[1] - \
        bean.find(sum(pdata>it,0) < plate.shape[0]/grid.win-grid.dims[0],0)
    
    grid.r -= radj * grid.win + grid.win/2
    grid.c -= cadj * grid.win + grid.win/2
    
    return grid
    
def estimate_initial_grid_aspect( plate, **params ):
    if 'sizestandard' not in params:
        params['sizestandard'] = [1853, 2765]
    
    # Compute grid spacing and dimensions
    grid = initialize_grid(plate, **params)
    dims = grid.dims
    win = grid.win
    
    # Identify the grid orientation
    tang = params['sizestandard'][0] / (params['sizestandard'][1] + 0.0)
    ratiofun = lambda xp, yp: np.arctan( -(yp-xp*tang)/(yp*tang-xp) )
    yp, xp = plate.shape
    
    theta = ratiofun( xp, yp );
    if (np.mean(plate[0,np.floor(xp/2):]) > np.mean(plate[0,0:np.floor(xp/2)])):
        theta = -theta
        
    rotmat = np.array([[np.cos(theta), -np.sin(theta)],\
        [np.sin(theta), np.cos(theta)]])
    y, x = dims * win
    
    tmp = np.array([[-x/2, -y/2],[x/2, -y/2],[x/2, y/2],[-x/2, y/2]])
    coords = np.transpose(np.dot( rotmat, np.transpose(tmp) ))
    coords = coords + np.array([plate.shape[1], plate.shape[0]])/2
    
    return determine_grid_from_corners( coords, grid )

def estimate_initial_grid( plate, **params ):
    bean.default_param( params, 
        mode='aspectratio' )
    
    if params['mode'] == 'aspectratio':
        return estimate_initial_grid_aspect( plate, **params )
    elif params['mode'] == 'offset':
        return estimate_initial_grid_offset( plate, **params )
    else:
        raise Exception('Unsupported option for MODE: %s' % params['mode'])
        
# Determine colony grid
def determine_colony_grid( plate, **params ):
    # Compute Grid
    if 'initialgrid' in params:
        grid = params['initialgrid'];
    else:
        grid = estimate_initial_grid( plate, **params );
         
    # Adjust Grid
    return adjust_grid( plate, grid, **params );
       
def manual_grid( plate, **params ):
    
    # Get corner points
    fig = plt.figure(figsize=(12,8))
    fig.suptitle('Manual Grid Alignment')
    ax = fig.add_subplot(111)
    ax.set_title('Please select the four corners of the colony grid')
    ax.imshow(plate, cmap='gray')
    def add_num(x, y, numstr):
        ax.text(x, y, numstr, transform=ax.transAxes, color='red', 
            fontsize=16, fontweight='bold')
    add_num(0.05, 0.95, '1')
    add_num(0.95, 0.95, '2')
    add_num(0.95, 0.05, '3')
    add_num(0.05, 0.05, '4')
    fig.canvas.draw()
    corners = bean.get_points(4, color='r')
    plt.close(fig)
    
    # Determine grid from corners
    print "Computing grid position..."
    sys.stdout.flush()
    grid = initialize_grid(plate)
    grid = determine_grid_from_corners( corners, grid )

    # Adjust Grid
    grid = adjust_grid( plate, grid, **params )
    
    # Verify
    fig = plt.figure(figsize=(12,8))
    fig.suptitle('')
    ax = fig.add_subplot(111)
    ax.set_title('Is the colony grid correct? (Respond in terminal)')
    ax.imshow(plate, cmap='gray')
    xl, yl = plt.xlim(), plt.ylim()
    ax.scatter(grid.c, grid.r, s=5)
    plt.xlim(xl), plt.ylim(yl)
    fig.canvas.draw()
    plt.pause(0.01)
    
    resp = raw_input('Is the colony grid correct? Y/n/cancel: ')
    plt.close(fig)
    
    if resp.lower()[0] == 'y':
        return grid
    elif resp.lower()[0] == 'n':
        return manual_grid(plate, **params)
    else:
        return None
    

