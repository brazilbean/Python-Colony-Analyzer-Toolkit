## PyCAT - Python Colony Analyzer Toolkit
# Gordon Bean, gbean@ucsd.edu
# April 2013

## pycat.io Module

# Imports
import sys
sys.path.append('/cellar/users/gbean/Dropbox/pyfiles')
import bean

import numpy as np
import matplotlib.image as mpimg

# Functions
def rgb2gray(rgb):
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray
    
def load_image(filename, method='grayscale'):
    """ Loads the image and averages across the RGB channels"""
    
    # Load image
    img = mpimg.imread(filename);
    
    # Average across the three channels
    if method == 'mean':
        return np.mean(img,2).astype('uint8')

    elif method == 'original':
        return np.mean(img,2)

    elif method == 'grayscale':
        return rgb2gray(img)

    else:
        raise Exception('Method %s not supported.'%method)
        
def crop_image(img, offset=0, background_thresh=0.9):
    """ Crops the image to the plate"""
    
    # Estimate background intensity
    foo = np.mean(img,0);
    fmid = np.floor( len(foo)/2 );
    w = np.floor(fmid/6);
    bthresh = min( foo[fmid-2*w : fmid+2*w] );
    
    # Crop
    crop = np.array([0]*4);
    tmpb = np.mean(img < bthresh,1) < background_thresh;
    crop[0] = max(1, bean.find(tmpb,0) - offset);
    crop[1] = min(len(tmpb), bean.find(tmpb, -1) + offset);
    
    tmpb = np.mean(img < bthresh,0) < background_thresh;
    crop[2] = max(1, bean.find(tmpb,0) - offset);
    crop[3] = min(len(tmpb), bean.find(tmpb,-1) + offset);
    
    return img[crop[0]:crop[1], crop[2]:crop[3]];

def load_plate(filename):
    """ Loads and crops the plate image file"""
    img = load_image(filename);
    return crop_image(img);
    
