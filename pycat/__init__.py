## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu

## Modules
import io, grid, threshold, quantify, analysis

# __all__ definintion
__all__ = ['io','grid','threshold','quantify', 'analysis']

# More Imports
import numpy as np
import pickle

# Bean's bag-o-tricks
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

# pycat imports
import pycat.grid as pygrid
import pycat.threshold as pythresh
import pycat.quantify as pyquant

def measure_colony_sizes( plate, **params ):
    bean.default_param( params, 
        manualgrid = False, 
        thresholdfunction = pythresh.local_gaussian, 
        sizefunction = pyquant.threshold_bounded )

    # Load plate
    if isinstance(plate, basestring):
        # "plate" is the file name
        plate = io.load_plate(plate)
    
    # Average plate across RGB
    if len(plate.shape) > 2:
        plate = np.mean(plate,2)
        
    # Determine grid
    if params['manualgrid']:
        params['grid'] = pygrid.manual_grid(plate, **params)
        if params['grid'] is None:
            # The user canceled the alignment
            return None, None
    else:
        if 'grid' not in params:
            params['grid'] = pygrid.determine_colony_grid( plate, **params )
    
    pgrid = params['grid']
    
    # Intensity thresholds
    if not hasattr(pgrid, 'thresh'):
        pgrid.thresh = params['thresholdfunction']( plate, pgrid )
    
    # Measure colony size
    size_fun = params['sizefunction']
    rrr, ccc = range(0, pgrid.dims[0]), range(0, pgrid.dims[1])
    sizes = [ [ size_fun( plate, pgrid, r, c ) for c in ccc] for r in rrr];
  
    return np.array(sizes), pgrid
    
        
def analyze_image( filename, **params ):
    """ Quantify the colony sizes in the image and save the output to file. """

    output_extension = '.cs.txt'
    
    # Measure colony sizes
    cs, pgrid = measure_colony_sizes( filename, **params )
    if cs is None:
        # The analysis failed, or was canceled
        return
    
    # Print the .txt file
    txtfile = filename + output_extension
    cc, rr = np.meshgrid(range(0, pgrid.dims[1]), range(0, pgrid.dims[0]))
    with open(txtfile, 'w') as fid:
        fid.write("row\tcolumn\tsize\n")
        for r, c, s in zip(bean.ind(rr), bean.ind(cc), bean.ind(cs)):
            fid.write("%i\t%i\t%i\n" % (r+1, c+1, s))
        
    # Save grid data
    gridfile = filename + '.info.pkl'
    with open(gridfile, 'w') as fid:
        pickle.dump(pgrid, fid, pickle.HIGHEST_PROTOCOL)
    
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
        print('Parallel mode not yet supported.')
        pass
#        # Number of threads
#        if type(params['parallel']) is int:
#            nthreads = params['parallel']
#        else:
#            nthreads = 4
#
#        bean.verbose( verb, 
#            ' Analyzing %s with %i threads' % (imagedir, nthreads) )
#            
#        fun = lambda f: analyze_image(f, **params)
#        map( fun, files )
#        #bean.multimap( fun, files, None, nthreads )
        
    else:
        for ff in files:
            try:
                bean.verbose( verb, ' Analyzing: %s' % ff )
                analyze_image( ff, **params )
            except Exception:
                print "Image failed: \n %s \n" % ff
            
    