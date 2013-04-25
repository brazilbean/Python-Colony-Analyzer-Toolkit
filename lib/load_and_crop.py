## Python Colony Analyzer Toolkit
# Load and Crop Module
# Gordon Bean, April 2013

# Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import imp
bean = imp.load_source('bean','/cellar/users/gbean/Dropbox/pyfiles/bean.py');

# Functions
def rgb2gray(rgb):
    """ Borrowed from stackoverflow and wikipedia: 
        http://stackoverflow.com/questions/12201577/
            convert-rgb-image-to-grayscale-in-python
        http://en.wikipedia.org/wiki/Grayscale#Converting_color_to_grayscale""" 
    
    r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b

    return gray
    
def load_image(filename, method='grayscale'):
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
    img = load_image(filename);
    return crop_image(img);
    

