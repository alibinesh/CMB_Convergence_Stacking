#import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords
import scipy.ndimage.interpolation as sni

from flipper import liteMap



def min_dim_of_list_of_arrays(list_of_arrays):
    '''
    Will return list (y,x) of the minimum dimensions of
    all of the submaps
    '''
    #print "List of arrays: ", list_of_arrays
    dim_list_x = []
    dim_list_y = []
    for z in np.arange(0, len(list_of_arrays)):
    	#print list_of_arrays[z]
        dat = list_of_arrays[z]
        #print dat
        dim_list_x.append(dat.shape[1])
        dim_list_y.append(dat.shape[0])
        
        #get min of x and y dim of these
    min_dim_x = min(dim_list_x)
    min_dim_y = min(dim_list_y)

    return (min_dim_y, min_dim_x)

#def Stack_on_Positions(map,N,RA,Dec,N_objects,Radius):
#    stack = np.zeros([Radius*2,Radius*2])
#    counter = 0
#    i = 0
#    while (i < N_objects):
#        xc = RA[i]
#        yc = Dec[i]
#        if ((xc > Radius) and (xc < N-Radius)):
#            if ((yc > Radius) and (yc < N-Radius)):
#                stack += map[xc-Radius:xc+Radius,yc-Radius:yc+Radius]
#                counter +=1
#        i = i + 1
#    return(stack/counter)

def crop_arraynew(array1, array2):
    '''
    returns the same array back with it cropped to ensure that it has the same dimensions
    as what the desired input is
    '''
    if array1.Nx == array2.Nx and array1.Ny == array2.Ny:
        return array1
    if array1.Nx > array2.Nx or array1.Ny > array2.Ny:
#        print array1
        smap2 = array2.copy()
        smap2.data[:] = 0
#        print array1.data[:].shape, "array1 shape"
#        print array2.data[:].shape, "array2 shape"
#        print smap2.data[:].shape, "Smap2 shape"
        smap2.data[:] = array1.data[0:smap2.Ny , 0:smap2.Nx] #+ array2.data[0:smap2.Ny, 0:smap2.Nx]
        array1 = smap2.copy()
#        print array1.data[:].shape, "new array"
        
        return array1
    if array1.Nx < array2.Nx or array1.Ny < array2.Ny:
        smap2 = array1.copy()
        smap2.data[:] = 0
        smap2.data[:] = array1.data[0:smap2.Ny , 0:smap2.Nx] #+ array2.data[0:smap2.Ny, 0:smap2.Nx]
        array1 = smap2.copy()
        
        return array1
def crop_array(array, dim_list):
    '''
    returns the same array back with it cropped to ensure that it has the same dimensions
    as what the desired input is
    '''
    if array.shape == dim_list:
        #print "original array", array.shape
        return array
    if array.shape != dim_list:
        new_array = array[:dim_list[0], :dim_list[1]]
        
        
        return new_array

def stack_rescale(stacked_map, R_filt, R_filt_ref):
    
    rescale_param = R_filt_ref/R_filt
    #rescale_param = R_filt/R_filt_ref #flipped it around, making smaller maps larger
    
    print "The rescaling parameter from the filters is: ", rescale_param

    #now we have to make sure it is a multiple of the pixel scale
    stacked_dim_x = stacked_map.shape[0] #I just picked x, but it should realistically be the same
    pix_frac = 1.0/stacked_dim_x #rescale_param should be some multiple of this

    rescale_frac = float(rescale_param)/float(pix_frac) # how evenly does this go into how many pixels there are?

    #now we want to get to the closest multiple of this, so we would round it to the nearest whole integer value
    #rescale_factor should be a decimal value between 0 and 1 which will resize your matrix
    rescale_factor = round(rescale_frac)*pix_frac
    
    #print "The rescaling factor from the pixels is: ", rescale_factor

    rescaled_stacked_map = sni.zoom(stacked_map, rescale_factor)

    return rescaled_stacked_map

def stack_all_rescaled_maps(rescale_stack_array):
    
    #get the dimensions of smallest stacked map bin
    #this is because sometimes they can be off by one or two pixel rows
    
    dim_list = min_dim_of_list_of_arrays(rescale_stack_array)
    
    total_stack = np.zeros(dim_list)
    
    for i in np.arange(0, len(rescale_stack_array)):
        stack_map_crop = crop_array(rescale_stack_array[i], dim_list)
        total_stack += stack_map_crop

    total_stack = total_stack/len(rescale_stack_array)
        
    return total_stack

def degrees_to_pix(r_deg,ref_map):
    r_arcmin = 60*r_deg #degrees to arcmins
    #pix_scale_x = 0.00014532290052643554 # for act_map_deep5_2.pixScaleX
    pix_scale_x = ref_map.pixScaleX
    pix_scale = (pix_scale_x)*((180/np.pi) *(60)) #I chose a random specific map to use
    r_pix = r_arcmin / pix_scale 
    return r_pix

def array_mask(data_array, R_filt, ref_map):
    '''
    masks out data outside a disk, this does the same as the
    flipper mask function. but then we take the average temp in that disk - this represents 
    the \Delta T of the map - since this would be used after already subtracting the outer circle
    annulus from the individual maps
    r is in arc minutes - we change it to pixels in there
    '''
    #calculate the middle pixel of the map
    cent_x = data_array.shape[0]/2.0
    cent_y = data_array.shape[1]/2.0

    #changing r to pixel scales
    r_pix = degrees_to_pix(R_filt, ref_map)
    
    #print "yo i changed this"
    
    y_size, x_size = data_array.shape

    x = np.arange(0,x_size)
    y = np.arange(0,y_size)

    xv, yv = np.meshgrid(x,y)

    distances = np.sqrt((cent_x - xv)**2. + (cent_y - yv)**2.)

    circ_y,circ_x = np.where(distances < r_pix)

    circ_ind = np.where(distances < r_pix)

    circ_vals = data_array[circ_ind]
    
    delta_T = np.sum(circ_vals)/float(len(circ_vals))
    
    '''
    array_size = data_array.shape
    mask = np.zeros(array_size)
    num_pix = 0 #number of non-zero pixels
    for i in np.arange(0,len(mask[0])):
        for j in np.arange(0, len(mask)):
            dist = np.sqrt((i-cent_x)**2 + (j-cent_y)**2)
            if dist < r_pix:
                mask[j][i] = 1
                num_pix +=1

    inner_circ = mask*data_array
    delta_T = np.sum(inner_circ)/float(num_pix)
    '''
    return delta_T
    



