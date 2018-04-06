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

from flipper import *


'''
Any functions which relate to actual stacking of the voids, pre analysis
'''

def min_dim_of_list(list_of_smaps):
    '''
    Will return list (y,x) of the minimum dimensions of
    all of the submaps
    '''
    dim_list_x = []
    dim_list_y = []
    for z in np.arange(0, len(list_of_smaps)):
        dat = list_of_smaps[z].data
        dim_list_x.append(dat.shape[1])
        dim_list_y.append(dat.shape[0])
        
        #get min of x and y dim of these
    min_dim_x = min(dim_list_x)
    min_dim_y = min(dim_list_y)

    return (min_dim_y, min_dim_x)

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


def crop_array(array, dim_list):
    '''
    returns the same array back with it cropped to ensure that it has the same dimensions
    as what the desired input is
    '''
    if array.shape == dim_list:
        return array
    if array.shape != dim_list:
        new_array = array[:dim_list[0], :dim_list[1]]
        return new_array

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

def stacked_voids(filt_array):
    '''
    Generates a stacked images of all input void submaps
    and outputs that image
    '''
    
    if filt_array == []: #empty
        return
        
    else:
        #get data from the litemaps
        filt_data_array = []
        for z in np.arange(0, len(filt_array)):
            filt_data_array.append(filt_array[z].data)
        
        #get the smallest submap dimensions
        dim_list = min_dim_of_list(filt_array)
        filt_data_sum = np.zeros(dim_list)
        counter = 0.
        for k in np.arange(0, len(filt_array)):
            #croppping all subamps to make sure they all have the same dimensions
            filt_dat_crop = crop_array(filt_data_array[k], dim_list)
            filt_data_sum = filt_data_sum + filt_dat_crop
            print "stacked info for void # ", str(k), np.mean(filt_data_sum), np.std(filt_data_sum)
            counter += 1.
            
        #filt_data_avg = filt_data_sum/len(filt_array)
        print "the counter we're dividing by: ", counter
        filt_data_avg = filt_data_sum/float(counter)
        print "Final dimension of stacked map: ", filt_data_avg.shape
        
        '''
        # we dont care to use litemaps anymore
        #new liteMap for stacked data
        filt_map = filt_array[0] #first filtered liteMap in input list
        stackLiteMap = filt_map.copy() #copying dimensions of liteMap 
        stackLiteMap.data = filt_data_avg #inputting stacked data into new liteMap
        '''
        
        return filt_data_avg
    
'''
def stack_rescale(stacked_map, R_filt, R_filt_ref):
    
    rescale_param = R_filt_ref/R_filt
    #we have R_filt_ref/R_filt instead of the other way around
    #because we want to stack smaller maps on a larger scale
    #because of the way the zoom function works
    rescaled_stacked_map = sni.zoom(stacked_map, rescale_param)

    return rescaled_stacked_map
'''

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
    
def void_temperature_profile(delt_T_for_bin):
    #creating a binned version of Cai et al. (2016)
    #this didnt work out so well - come back to this
    #void temperature profile
    
    err_delT_bin_upper = [] #will be std of delt_T for that bin
    err_delT_bin_lower = []
    avg_delT_bin = []
    avg_r_eff_bin = [] #not necessary?
    
    for z in np.arange(0, len(delt_T_for_bin)):
        std_delT = np.std(delt_T_for_bin[z])
        avg_delT = np.mean(delt_T_for_bin[z])
        avg_delT_bin.append(avg_delT)
        err_delT_bin_upper.append(std_delT + avg_delT)
        err_delT_bin_lower.append(avg_delT - std_delT)
        
	#print len(avg_delT_bin)
	#print avg_delT_bin
	'''
	print "average delta_T in bins: ", avg_delT_bin
	print "upper bound error on average delta_T in bins", err_delT_bin_upper
	print "lower bound error on average delta_T in bins", err_delT_bin_lower
	'''
	
    return avg_delT_bin, err_delT_bin_upper, err_delT_bin_lower
    
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
    
    print "yo i changed this"
    
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
    
    
############################################################################
#S/N stuff

def SN_std(mod_smaps_list, condition, R):
	#using the whole map
	
	if condition == "whole map":
		std_temp_list = np.zeros(len(mod_smaps_list))
		for k in np.arange(0,len(mod_smaps_list)):
			mod_dat = mod_smaps_list[k].data
			std_temp_list[k] = np.std(mod_dat)
		
	#using only the temp outside the inner circle
	if condition == "outside temp":
		mod_map = mod_smaps_list[k]
		std_temp_list = np.zeros(len(mod_smaps_list))
		mapout = mod_map.mask(void_RA,void_Dec,r_pix,mask_lo=15, mask_hi=25)
		std_mapout = np.std(mapout.data)
		std_temp_list[k].append(std_mapout)
		
		#num_pix_index = np.where(delTmap.data != 0.)
		
		
	return std_temp_list

def SN_significance(delT_list, std_temp_list):
	#calculates the significance of the deltaT measure in each stamp
	#by dividing that value by the standard deviation of the temperature in that stamp
	sig = np.absolute(np.array(delT_list)/np.array(std_temp_list))
	return sig


