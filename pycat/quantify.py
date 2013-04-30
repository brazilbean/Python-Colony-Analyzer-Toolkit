## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.threshold Module

# Imports
import numpy as np
import pickle
import os

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
    bean.default_param( params, 
        manualgrid = False, 
        thresholdfunction = pythresh.local_gaussian, 
        sizefunction = threshold_bounded )

    # TODO
    # - check for plate type and dimensions
    # - manual grid
    
    # Load plate
    if isinstance(plate, basestring):
        # "plate" is the file name
        plate = io.load_plate(plate)
    
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
        grid.thresh = params['thresholdfunction']( plate, grid )
    
    # Measure colony size
    size_fun = params['sizefunction']
    rrr, ccc = range(0, grid.dims[0]), range(0, grid.dims[1])
    sizes = [ [ size_fun( plate, grid, r, c ) for c in ccc] for r in rrr];
  
    return np.array(sizes), grid
    
def analyze_image( filename, **params ):
    """ Quantify the colony sizes in the image and save the output to file. """

    output_extension = '.cs.txt'
    
    # Measure colony sizes
    cs, grid = measure_colony_sizes( filename, **params )
    if cs is None:
        # The analysis failed, or was canceled
        return
    
    # Print the .txt file
    txtfile = filename + output_extension
    cc, rr = np.meshgrid(range(0, grid.dims[1]), range(0, grid.dims[0]))
    with open(txtfile, 'w') as fid:
        fid.write("row\tcolumn\tsize\n")
        for r, c, s in zip(bean.ind(rr), bean.ind(cc), bean.ind(cs)):
            fid.write("%i\t%i\t%i\n" % (r+1, c+1, s))
        
    # Save grid data
    gridfile = filename + '.info.pkl'
    with open(gridfile, 'w') as fid:
        pickle.dump(grid, fid, pickle.HIGHEST_PROTOCOL)
    
def analyze_directory_of_images( imagedir, **params ):
    bean.default_param( params, 
        filepattern = '\.JPG$', 
        verbose = False, 
        parallel = False )
    
    # Get image files
    files = bean.dirfiles( imagedir, params['filepattern'] )
    
    # Scan the files
    verb = params['verbose']
    if params['parallel']:
        #print('Parallel mode not yet supported.')
        #pass
        # Number of threads
        if type(params['parallel']) is int:
            nthreads = params['parallel']
        else:
            nthreads = 4

        bean.verbose( verb, 
            ' Analyzing %s with %i threads' % (imagedir, nthreads) )
            
        fun = lambda f: analyze_image(f, **params)
        map( fun, files )
        #bean.multimap( fun, files, None, nthreads )
        
    else:
        for ff in files:
            try:
                bean.verbose( verb, ' Analyzing: %s' % ff )
                analyze_image( ff, **params )
            except Exception:
                print "Image failed: \n %s \n" % ff
            
            
    
    
