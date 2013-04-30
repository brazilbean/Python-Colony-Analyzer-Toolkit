## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.grid Module

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Bean's bag-o-tricks
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

# pycat imports
import pycat.io as io

## Classes
class Grid:
    
    def __init__ (self, plate=None, manual=False, **params ):
        if plate is None:
            return
            
        if manual:
            pass
        else:
            self = estimate_initial_grid( plate, **params )
            return self
    
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
        other = Grid()
        for k,v in vars(self).iteritems():
            if isinstance(v, np.ndarray):
                setattr(other, k, np.array(v))
            else:
                setattr(other, k, v)
        return other
        
## Methods
# Get box
def get_box(plate, rpos, cpos, win):
    """ Returns a copy of the data in plate centered at (rpos,cpos) 
         within a window of 2*win """
    ir = lambda x: int(round(x))
    return np.array( 
        plate[max(0,ir(rpos-win)) : min(plate.shape[0],ir(rpos+win)),
        max(0,ir(cpos-win)) : min(plate.shape[1],ir(cpos+win))] )
    
def set_box(plate, box, rpos, cpos):
    win = (box.shape[0])/2 
    
    rpos = np.round(rpos)
    cpos = np.round(cpos)
    
    if box.dtype is np.dtype(bool):
        plate[rpos-win : rpos+win, cpos-win : cpos+win] |= box
    else:
        plate[rpos-win : rpos+win, cpos-win : cpos+win] = box
        
# Pixel mode
def pixelmode(box):
    boxs = np.sort(bean.ind(box))
    n, xx = np.histogram(boxs, range(int(np.min(boxs)), int(np.max(boxs))))

    start = 0
    pos = 0
    while pos < len(boxs):
        if boxs[pos] == boxs[start]:
            pass
        else:
            boxs[start:pos] = np.linspace(boxs[start], boxs[pos], pos-start)
            start = pos
        pos += 1

    w = int( max(n) * 2 )
    #tmp = [ np.sum(boxs[ii-w:ii+w]-boxs[ii-w]) for ii in range(w,len(boxs)-w) ]
    tmp = [ boxs[ii+w]-boxs[ii-w] for ii in range(w,len(boxs)-w) ]
    return boxs[np.argmin(tmp) + w]

# Estimate intensity threshold
def estimate_intensity_threshold( plate ):
    # Find middle box
    w = np.floor(min(plate.shape)/10);
    tmp = np.floor(np.array(plate.shape)/2);
    rmid = tmp[0];
    cmid = tmp[1];
    box = plate[rmid-w:rmid+w, cmid-w:cmid+w];
    
    thresh = (min(bean.ind(box)) + max(bean.ind(box))) / 2;
    return thresh;

# Estimate grid spacing
def estimate_grid_spacing( plate ):
    # Estimate intensity threshold
    thresh = estimate_intensity_threshold( plate );
    
    # Get middle of plate
    w = np.floor(min(plate.shape)/10);
    tmp = np.floor(np.array(plate.shape)/2);
    rmid = tmp[0];
    cmid = tmp[1];
    box = plate[rmid-w:rmid+w, cmid-w:cmid+w];
        
    # Identify spot centers
    cents = bean.centroids(box>thresh)
    areas = bean.areas(box>thresh)
    
    cents = cents[ areas > 10, :];
    cents = cents[:,:,np.newaxis];
    
    dists = np.sqrt( np.sum( \
    (np.transpose(cents, (0, 2, 1)) - np.transpose(cents, (2, 0, 1)))**2, -1) )
    tmp = np.max(dists);
    for ii in range(0, dists.shape[0]): dists[ii,ii] = tmp;
    
    return int(np.round( bean.pmode( np.min( dists, 0 ) ) ));

# Estimate grid dimensions
def estimate_dimensions( image, win ):
    tmp = np.log( np.array(image.shape) / np.array([8, 12]) / win )/np.log(2)
    tmp = np.array([8, 12]) * 2 ** np.floor( tmp )
    return np.array([int(ii) for ii in tmp])
    
# Determine grid from corners
def determine_grid_from_corners( corners, grid ):
    dims = grid.dims
    cc, rr = np.meshgrid( range(0, dims[1]), range(0, dims[0]) )
    
    # Grid factors
    rrr = np.transpose(np.array([[0, 0, dims[0]-1, dims[0]-1]]))
    ccc = np.transpose(np.array([[0, dims[1]-1, dims[1]-1, 0]]))
    
    lq = np.linalg.lstsq
    rfact = lq( np.hstack((np.ones((4,1)), rrr, ccc)), corners[:,1])[0]
    cfact = lq( np.hstack((np.ones((4,1)), rrr, ccc)), corners[:,0])[0]
    
    # Compute grid position
    n = np.prod(rr.shape)
    tmp = np.hstack((np.ones((n,1)), \
        bean.ind(rr,nx1=True), bean.ind(cc,nx1=True)))
    grid.r = np.reshape( np.dot( tmp, rfact ), dims )
    grid.c = np.reshape( np.dot( tmp, cfact ), dims )
    
    return grid
    
# Estimate initial grid
def estimate_initial_grid( plate, **params ):
    if 'sizestandard' not in params:
        params['sizestandard'] = [1853, 2765];
    
    # Compute grid spacing and dimensions
    if 'gridspacing' not in params:
        params['gridspacing'] = estimate_grid_spacing( plate );
    win = params['gridspacing'];
    
    if 'dimensions' not in params:
        params['dimensions'] = estimate_dimensions( plate, win );
    dims = params['dimensions'];
    
    # Initialize grid
    grid = Grid()
    grid.win = win
    grid.dims = dims
    
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

# Adjust grid
def adjust_grid( plate, grid, **params ):
    bean.default_param( params, \
        convergethresh = 3, \
        adjustmentwindow = int(np.round( grid.dims[0]/8.0 )), \
        finaladjust = True, \
        fitfunction = lambda r, c: \
            np.hstack((np.ones((np.prod(r.shape),1)), \
            bean.ind(r,nx1=True), bean.ind(c,nx1=True))) )
        
    # Subroutines
    def measure_offset( box ):
        win = (box.shape[0]-1)/2
        
        # Get centroid
        thresh = (np.min(box) + np.max(box))/2
        cents = bean.centroids( box > thresh )
        cents = cents - (win+1)
        tmp = np.all(abs(cents) < win/2,1)
        if sum(tmp)>0:
            return cents[bean.find(tmp,0),:]
        else:
            return np.nan
        
    def adjust_spot( plate, rpos, cpos, win ):
        box = get_box( plate, rpos, cpos, win )
        off = np.round(measure_offset( box ))
        
        if np.any(off > win/2):
            raise Exception('Offset too large - decrease the adjustment window')
        if np.any(np.isnan(off)):
            return np.nan, np.nan
        else:
            return (rpos + off[0], cpos + off[1])
            
    def minor_adjust_grid( plate, grid, rrr, ccc ):
        r0 = grid.r[0,0]
        dims = grid.dims
        win = grid.win
        
        rtmp, ctmp = bean.nans(grid.r.shape), bean.nans(grid.r.shape)
        
        for rr in rrr:
            for cc in ccc:
                rtmp[rr,cc], ctmp[rr,cc] = adjust_spot \
                    (plate, grid.r[rr,cc], grid.c[rr,cc], win)
        
        # Estimate grid parameters
        iii = np.logical_and(~np.isnan(rtmp), ~np.isnan(ctmp))
        cc, rr = np.meshgrid( range(0, dims[1]), range(0, dims[0]) )
        Afun = params['fitfunction']
        
        lq = np.linalg.lstsq
        rfact = lq( Afun(rr[iii], cc[iii]), rtmp[iii] )[0]
        cfact = lq( Afun(rr[iii], cc[iii]), ctmp[iii] )[0]
        
        # TODO - include grid.factors?
        
        # Compute grid position
        grid.r = np.reshape(np.dot(Afun(rr,cc), rfact), grid.r.shape)
        grid.c = np.reshape(np.dot(Afun(rr,cc), cfact), grid.c.shape)
        
        # Estimate convergence
        fitfact = np.abs( grid.r[0,0] - r0 )
        
        return grid, fitfact

    # End of minor_adjust_grid    

    fitfact = grid.win
    aw = params['adjustmentwindow']
    
    # Initial adjustment
    #while fitfact > params['convergethresh']:
    for ii in [1]:
        # Adjust internal spots
        rrr = grid.dims[0]/2 + range(-aw, aw+1)
        ccc = grid.dims[1]/2 + range(-aw, aw+1)
        
        grid, fitfact = minor_adjust_grid( plate, grid, rrr, ccc )
        
    # Final adjustment
    if params['finaladjust']:
        rrr = np.round( np.linspace(0, grid.dims[0]-1, 2*aw) )
        ccc = np.round( np.linspace(0, grid.dims[1]-1, 2*aw) )
        
        grid, fitfact = minor_adjust_grid( plate, grid, rrr, ccc )
        
    # Extra stuff
    return grid
        
# Determine colony grid
def determine_colony_grid( plate, **params ):
    # Compute Grid
    if 'initialgrid' in params:
        grid = params['initialgrid'];
    else:
        grid = estimate_initial_grid( plate, **params );
         
    # Adjust Grid
    return adjust_grid( plate, grid, **params );
       
# View plate image
def view_plate_image( filename, **params ):
    # Default parameters
    bean.default_param( params, \
        showimage = True, \
        interactive = False, \
        showgrid = False, \
        shownotes = True, \
        showaxes = False, \
        notes = [], \
        returnnotes = False, \
        applythreshold = False, \
        maskthreshold = False )
    bean.default_param( params, newfigure = params['interactive'] )
    bean.default_param( params, \
        gridspecs = {'s': 10, 'c':'#0000FF', 'marker':'o'}, \
        notespecs = {'s': 20, 'c':'#FF0000', 'marker':'o','alpha':0.5} )
    
    # Load image and image info
    if isinstance(filename, basestring):
        plate = io.load_plate( filename )
    else:
        # Assume filename is the plate
        plate = filename
    
    # Show image
    if params['newfigure']:
        fig = plt.figure()
    else:
        fig = plt.gcf()
    
    def show_plate():
        plt.clf()
        plt.imshow(plate, aspect='auto', interpolation='none', cmap='gray' )
    show_plate()
    xl = plt.xlim();
    yl = plt.ylim();
    ax = fig.gca();
        
    # Show grid
    if 'grid' not in params:
        grid = determine_colony_grid( plate )
    else:
        grid = params['grid']
    
    
    if params['showgrid']:
        gax = fig.add_axes(ax.get_position().bounds, frameon=False)
        plt.xlim(xl)
        plt.ylim(yl)
        def show_grid():
            plt.sca(gax)
            plt.cla()
            plt.scatter(grid.c, grid.r, **params['gridspecs']) 
            plt.xlim(xl)
            plt.ylim(yl)
            #fig.canvas.draw()
        show_grid()

    # Show notes
    notes = np.zeros(grid.dims) == 1
    for note in params['notes']:
        r, c = bean.ind2sub(note, grid.dims)
        notes[r,c] = True
    if params['shownotes']:
        nax = fig.add_axes(ax.get_position().bounds, frameon=False)
        plt.xlim(xl)
        plt.ylim(yl)
        def show_notes(notes):
            plt.sca(nax)
            plt.cla()
            plt.scatter(grid.c[notes], grid.r[notes], **params['notespecs'])
            plt.xlim(xl)
            plt.ylim(yl)
            fig.canvas.draw()
        if np.any(notes):
            show_notes(notes)
    
    # Interactive figure
    if params['interactive']:
        class clicked:
            click = False
            def __init__(self):
                self.click = False
            def on(self):
                self.click = True
            def off(self):
                self.click = False
            def isclicked(self):
                return self.click
        click = clicked()
                
        def get_coords(event):
            r, c = event.ydata, event.xdata
            foo = np.abs(grid.r - r) + abs(grid.c-c)
            r, c = bean.ind2sub(np.argmin(foo), grid.dims)
            return r, c
            
        def toggle_note(r, c):
            notes[r,c] = ~notes[r,c]
            pass
            
        def onclick(event):
            if click.isclicked():
                return
            else:
                click.on()
                # Get grid coordinate
                r, c = get_coords(event)
                #print r, c
                
                # Toggle note
                toggle_note(r, c)
                #print notes[r,c]
                
                # Refresh image
                if params['shownotes'] and np.any(notes):
                    show_notes(notes)
                click.off()
            
        fig.canvas.mpl_connect('button_press_event', onclick)
        
    plt.draw()
    plt.show()
    if params['returnnotes']:
        return notes


