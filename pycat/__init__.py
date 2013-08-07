## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu

## TODO
# Load grid coordinate option
# Multiple metrics option
# Fix view_plate_image

## Modules
from . import io as _io
from . import grid as _grid
from . import threshold as _threshold
from . import quantify as _quantify

# More Imports
import numpy as np
import matplotlib.pyplot as plt
import pickle

# Bean's bag-o-tricks
import bean # https://github.com/brazilbean/bean-python-toolkit

def measure_colony_sizes( plate, 
    imgloader = _io.PlateLoader(),
    grid = _grid.OffsetGridMethod(),
    threshold = _threshold.BackgroundOffset(),
    metric = _quantify.ColonyArea() ):

    ## Load plate
    record_loader = False
    if isinstance(plate, basestring):
        # "plate" is the file name
        filename = plate
        record_loader = True
        plate = imgloader(plate)
        
    else:
        # "plate" is the actual image - assume it is pre-processed
        pass
        
    ## Determine grid
    if isinstance(grid, _grid.Grid):
        pgrid = grid
    else:
        pgrid = grid(plate)
    if record_loader:
        pgrid.info['imgloader'] = imgloader
        pgrid.info['file'] = filename
        
    ## Intensity thresholds
    if pgrid.thresh is None:
        pgrid.thresh = threshold( plate, pgrid )
        pgrid.info['threshold'] = threshold
    
    # Measure colony size
    pgrid.info['metric'] = metric
    rrr, ccc = range(0, pgrid.dims[0]), range(0, pgrid.dims[1])
    sizes = [ [ metric( plate, pgrid, r, c ) for c in ccc] for r in rrr];
  
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
        plate = _io.load_plate( filename )
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
        grid = _grid.determine_colony_grid( plate )
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
             
## Finish up __init__

import io, grid, threshold, quantify, analysis
__all__ = ['io','grid','threshold','quantify', 'analysis']

