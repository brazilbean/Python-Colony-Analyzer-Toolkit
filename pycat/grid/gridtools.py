## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.grid.gridtools Module

import numpy as np

# Bean's bag-o-tricks
import bean # https://github.com/brazilbean/bean-python-toolkit

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
    
## Adjust Grid Functions ##
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

# Default coefficient function
def default_coefficient_function( r, c ):
    return np.hstack(( 
        np.ones((np.prod(r.shape),1)), 
        bean.ind(r,nx1=True), 
        bean.ind(c,nx1=True) ))
        
# adjust_subgrid
def adjust_grid_linear( plate, grid, rrr, ccc, 
    coeffunction = default_coefficient_function ):
        
    dims = grid.dims
    win = grid.win
    
    rtmp, ctmp = bean.nans(grid.r.shape), bean.nans(grid.r.shape)
    
    # Find true colony locations
    for rr, cc in zip(bean.ind(rrr),bean.ind(ccc)):
        rtmp[rr,cc], ctmp[rr,cc] = adjust_spot \
            (plate, grid.r[rr,cc], grid.c[rr,cc], win)
    
    # Estimate grid parameters
    iii = np.logical_and(~np.isnan(rtmp), ~np.isnan(ctmp))
    cc, rr = np.meshgrid( range(0, dims[1]), range(0, dims[0]) )
    
    lq = np.linalg.lstsq
    rfact = lq( coeffunction(rr[iii], cc[iii]), rtmp[iii] )[0]
    cfact = lq( coeffunction(rr[iii], cc[iii]), ctmp[iii] )[0]
    
    # TODO - include grid.factors?
    
    # Compute grid position
    grid.r = np.reshape(np.dot(coeffunction(rr,cc), rfact), grid.r.shape)
    grid.c = np.reshape(np.dot(coeffunction(rr,cc), cfact), grid.c.shape)
    
    return grid

def adjust_grid_polar( plate, grid, rrr, ccc ):
    '''Adjust the grid using polar coordinates.
    This method's strength is that it combines spacing and positions information
    from both the rows and the columns. It tends to be biased in accuracy
    towards the top-left corner.
    Note that rrr and ccc are paired subscripts, as returned from numpy.nonzero.
    '''
    
    dims = grid.dims
    win = grid.win
    rtmp, ctmp = bean.nans(grid.r.shape), bean.nans(grid.r.shape)
    
    # Find true colony locations
    for rr, cc in zip(bean.ind(rrr),bean.ind(ccc)):
        rtmp[rr,cc], ctmp[rr,cc] = adjust_spot \
            (plate, grid.r[rr,cc], grid.c[rr,cc], win)
    
    ## Convert to polar coordinates
    # Set top-left coordinate as reference
    ii = bean.find(~np.isnan(rtmp),1)
    r0, c0 = rtmp[ii], ctmp[ii]
    
    # Set all positions relative to reference
    rpos = rtmp - r0
    cpos = ctmp - c0
    
    # Compute rho (radius)
    rho = np.sqrt(rpos**2 + cpos**2)
    
    # Compute theta
    theta = np.atan2(-rpos, cpos)
    
    ## Compute expected positions (in polar)
    cc, rr = np.meshgrid( range(0, dims[1]), range(0, dims[0]) )
    r0i, c0i = np.ind2sub(dims, ii)
    rr = rr - r0i
    cc = cc - c0i
    
    rho_exp = np.sqrt(rr**2 + cc**2)
    theta_exp = np.atan2(-rr,cc)
    
    # Get theta factor
    theta_fact = np.median(theta - theta_exp)
    
    # Get rho factor
    rho_fact = np.median(rho / rho_exp)
    
    ## Return cartesian, updated coordinates
    grid.r = -rho_fact * rho_exp * np.sin(theta_exp + theta_fact) + r0
    grid.c = rho_fact * rho_exp * np.cos(theta_exp + theta_fact) + c0
    
    grid.info['theta'] = theta_fact
    
    if np.all(np.isnan(grid.r)):
        # TODO - throw an error - NaN grid.
        pass
    
    return grid
    
    
# Adjust grid
#def adjust_grid( plate, grid, **params ):
#    bean.default_param( params, \
#        convergethresh = 3, \
#        adjustmentwindow = int(np.round( grid.dims[0]/8.0 )), \
#        finaladjust = True )
#    
#    aw = params['adjustmentwindow']
#    
#    # Initial adjustment
#    #while fitfact > params['convergethresh']:
#    for ii in [1]:
#        # Adjust internal spots
#        rrr = grid.dims[0]/2 + np.array(range(-aw, aw+1))
#        ccc = grid.dims[1]/2 + np.array(range(-aw, aw+1))
#        ccc, rrr = np.meshgrid( ccc, rrr )
#        
#        grid = adjust_subgrid( plate, grid, rrr, ccc, **params )
#        
#    # Final adjustment
#    if params['finaladjust']:
#        rrr = np.round( np.linspace(0, grid.dims[0]-1, 2*aw) )
#        ccc = np.round( np.linspace(0, grid.dims[1]-1, 2*aw) )
#        ccc, rrr = np.meshgrid(ccc, rrr)
#        
#        grid = adjust_subgrid( plate, grid, rrr, ccc, **params )
#        
#    # Extra stuff
#    return grid
#        