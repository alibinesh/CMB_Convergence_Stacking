#import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

from flipper import liteMap


'''All functions to filter the galaxies and map to make them ready to stack'''

def mask(ra, dec, m, r):
    '''
    masks out data within a disk, this does the same as the
    flipper mask function
    '''
    x,y = m.skyToPix(ra,dec)
    dat = m.data
    array_size = dat.shape
    mask = np.zeros(array_size)
    num_pix = 0 #number of non-zero pixels
    for i in np.arange(0,len(mask[0])):
        for j in np.arange(0, len(mask)):
            dist = np.sqrt((i-x)**2 + (j-y)**2)
            if dist < r:
                mask[j][i] = 1
                num_pix +=1
    return mask*dat, num_pix


def generate_filtsmap(smap):
    '''
    Returns a list of filtered submap stamps around the void locations
    from smap_list()
    '''
    el = np.arange(10000)
    Fel = (1. - np.exp(-el**2/(2*10.**2)))*np.exp(-el**2/(2*300.**2))
    #Fel = 1/(1 + np.exp(-(el-1000))) # I just changed this, previously this was 1000
    filteredMap = smap.filterFromList([el,Fel])
    return filteredMap

def delt_T_filt(void_RA, void_Dec, r, m):
    '''
    returns a list of a list of delta_T filters meant to be applied to
    mode filtered maps
    '''
    #change r to r_pix
    #input r in degrees
    r_arcmin = 60*r #degrees to arcmins
    pix_scale = (m.pixScaleX)*((180/np.pi) *(60))
    r_pix = r_arcmin / pix_scale 
    mapout = m.mask(void_RA,void_Dec,r_pix,mask_lo=15, mask_hi=25)
    mapin = m.mask(void_RA,void_Dec,np.sqrt(2)*r_pix,mask_lo=15, mask_hi=25)
    outer_circle = m.copy()
    outer_circle.data = m.data - mapin.data
    inner_circle = m.copy()
    inner_circle.data = m.data - mapout.data
    diffmap = inner_circle.copy()
    diffmap.data[:]= outer_circle.data[:]-inner_circle.data[:]
    return diffmap, inner_circle


def mod_filt(filtmap, delTmap):
    '''
    takes in the mode filtered maps, and the delta T filter and subtracts them
    to get a modified T and filtered submap ready for stacking
    '''
    modmap = filtmap.copy()
    #fix the annulus mean bit!
    num_pix_index = np.where(delTmap.data != 0.)
    num_pix = len(num_pix_index[0])
    annulus_mean = np.sum(delTmap.data)/num_pix
    modmap.data = filtmap.data - annulus_mean
    return modmap

