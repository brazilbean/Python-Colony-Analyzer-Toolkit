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
def load_image(filename):
    # Load image
    img = mpimg.imread(filename);
    
    # Average across the three channels
    return np.mean(img,2);
    
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
    

