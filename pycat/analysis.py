## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.analysis Module

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Bean's bag-o-tricks
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

# pycat imports
import pycat.grid as pygrid
import pycat.io as io
import pycat.threshold as pythresh

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
        rad = np.min( np.diff( plt.gca().bbox.get_points(), axis=0) / 2 / dims )
        params['max_circ'] = np.pi * rad**2
    
    # Reorder plate, if necessary
    # TODO - port this functionality
    
    # Draw figure
    if params['style'] == 'imagesc':
        pass
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
        
        tmp_size = tmp_size / (max_size + 0.0)
        tmp_size[tmp_size > 1] = 1
        
        tmp_size = (max_circ - 1) * tmp_size
        
        plt.scatter( xx, yy, s=tmp_size + 1, c=tmp )
        plt.xlim(-1, dims[1])
        plt.ylim(-1, dims[0])
        plt.gca().invert_yaxis()
        plt.draw()
        
    else:
        pass
    




