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


def neg_gal_RA(gal_RA):
    '''takes any negative RA values and converts them to positive degrees
    input is a list, returns the same list but with positive values only'''
    for i in np.arange(0,len(gal_RA)):
       if gal_RA[i] < 0.0: #if negative
            gal_RA[i] = 360.0 + gal_RA[i]
    return gal_RA
    
def chunks(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

def gen_rand_bounds(m, N, bound_x, bound_y):
    ''' generates radom numbers taking into account galaxy stamp size so that patches
    wont go over map bounds'''
    
    if m.x0 > m.x1:
        a = m.x1 + bound_x
        b = m.x0 - bound_x
        c = m.y0 + bound_y
        d = m.y1 - bound_y
        random_RA = []
        random_Dec = []
        for k in np.arange(0,N):
            random_RA.append(random.uniform(a, b))
            random_Dec.append(random.uniform(c, d))
    
    if m.x0 < m.x1:
        a = m.x0 + bound_x
        b = m.x1 - bound_x
        c = m.y0 + bound_y
        d = m.y1 - bound_y
        random_RA = []
        random_Dec = []
        for k in np.arange(0,N):
            random_RA.append(random.uniform(a,b))
            random_Dec.append(random.uniform(c,d))
    
        
    return random_RA, random_Dec

def del_galaxies_close_to_submap_bounds(RA_list, Dec_list, bound_x, bound_y, m):
    #def del_galaxys_close_to_submap_bounds(RA_list, Dec_list, bound_x, bound_y, m):
    ''' Deletes galaxy locations which are too close to the submap bounds
this gets rid of assertion errors'''
    index_del = []
    
    if m.x0 > m.x1:
        for i in np.arange(0,len(RA_list)):
            if RA_list[i] - bound_x < m.x1:
                index_del.append(i)
            if RA_list[i] + bound_x > m.x0:
                index_del.append(i)
            if Dec_list[i] + bound_y > m.y1:
                index_del.append(i)
            if Dec_list[i] - bound_y < m.y0:
                index_del.append(i)
    
    if m.x0 < m.x1:
        for i in np.arange(0,len(RA_list)):
            if RA_list[i] + bound_x > m.x1:
                index_del.append(i)
            if RA_list[i] - bound_x < m.x0:
                index_del.append(i)
            if Dec_list[i] + bound_y > m.y1:
                index_del.append(i)
            if Dec_list[i] - bound_y < m.y0:
                index_del.append(i)
    
    RA_list_new = np.delete(RA_list, index_del)
    Dec_list_new = np.delete(Dec_list, index_del)
    return RA_list_new, Dec_list_new

#def del_galaxies_close_to_submap_bounds(RA_list, Dec_list, r_ang_list, r_eff_list, bound_x, bound_y, m):
#    #def del_galaxys_close_to_submap_bounds(RA_list, Dec_list, bound_x, bound_y, m):
#    ''' Deletes galaxy locations which are too close to the submap bounds
#this gets rid of assertion errors'''
#    index_del = []
#    
#    if m.x0 > m.x1:
#        for i in np.arange(0,len(RA_list)):
#            if RA_list[i] - bound_x < m.x1:
#                index_del.append(i)
#            if RA_list[i] + bound_x > m.x0:
#                index_del.append(i)
#            if Dec_list[i] + bound_y > m.y1:
#                index_del.append(i)
#            if Dec_list[i] - bound_y < m.y0:
#                index_del.append(i)
#    
#    if m.x0 < m.x1:
#        for i in np.arange(0,len(RA_list)):
#            if RA_list[i] + bound_x > m.x1:
#                index_del.append(i)
#            if RA_list[i] - bound_x < m.x0:
#                index_del.append(i)
#            if Dec_list[i] + bound_y > m.y1:
#                index_del.append(i)
#            if Dec_list[i] - bound_y < m.y0:
#                index_del.append(i)
#    
#    RA_list_new = np.delete(RA_list, index_del)
#    Dec_list_new = np.delete(Dec_list, index_del)
#    r_eff_new = np.delete(r_eff_list, index_del)
#    r_ang_new = np.delete(r_ang_list, index_del)
#    return RA_list_new, Dec_list_new, r_ang_new, r_eff_new

def radial_bin(condition, gal_r, gal_RA, gal_Dec, gal_r_eff, n):
    '''note: the input gal_RA, gal_Dec, and gal_r are LISTS'''
    #if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
	#if the argument is 'angular radius', we'll bin by angular radius
	
    if condition == 'angular radius':
    	#meaning we're binning by the angular radius

		#set up a dictionary which related gal_index to other gal properties
		#dictionary of gal_RA, gal_Dec, gal_r
		gal_dict = dict((z[0],list(z[1:])) for z in zip(gal_r,gal_RA,gal_Dec,gal_r_eff))
		gal_r_sort = sorted(gal_r, key=float)
		gal_r_split = chunks(gal_r_sort, n)
		#now len(gal_r_split) should be n
		#empty pre beinned gal_RA and Dec lists, based off of the input bin number
		gal_RA_split = [[] for _ in range(n)]
		gal_Dec_split = [[] for _ in range(n)]
		gal_r_eff_split = [[] for _ in range(n)]
		gal_r_avg = []

	
		for i in np.arange(0,len(gal_r_split)):
			#compute average gal_r for each bin
			#this will be used for the R_filt = 0.7*gal_r_avg
			gal_r_avg.append(np.mean(gal_r_split[i]))
			for k in np.arange(0, len(gal_r_split[i])):
				#splitting up gal_RA based on bin
				gal_RA_split[i].append(gal_dict[gal_r_split[i][k]][0])
				gal_Dec_split[i].append(gal_dict[gal_r_split[i][k]][1])
				gal_r_eff_split[i].append(gal_dict[gal_r_split[i][k]][2])
			
		gal_r_eff_avg = [np.mean(gal_r_eff_split[i]) for i in range(0, len(gal_r_split))]

		return gal_r_split, gal_RA_split, gal_Dec_split, gal_r_avg, gal_r_eff_split, gal_r_eff_avg
	
    elif condition == 'effective radius':
		#set up a dictionary which related gal_index to other gal properties
		#dictionary of gal_RA, gal_Dec, gal_r
		gal_dict = dict((z[0],list(z[1:])) for z in zip(gal_r_eff,gal_RA,gal_Dec,gal_r))
		gal_reff_sort = sorted(gal_r_eff, key=float)
		gal_reff_split = chunks(gal_reff_sort, n)
		#now len(gal_r_split) should be n
		#empty pre beinned gal_RA and Dec lists, based off of the input bin number
		gal_RA_split = [[] for _ in range(n)]
		gal_Dec_split = [[] for _ in range(n)]
		gal_r_split = [[] for _ in range(n)]
		gal_r_eff_avg = []
		
		for i in np.arange(0,len(gal_reff_split)):
			#compute average gal_r for each bin
			#this will be used for the R_filt = 0.7*gal_r_avg
			gal_r_eff_avg.append(np.mean(gal_reff_split[i]))
			for k in np.arange(0, len(gal_reff_split[i])):
				#splitting up gal_RA based on bin
				gal_RA_split[i].append(gal_dict[gal_reff_split[i][k]][0])
				gal_Dec_split[i].append(gal_dict[gal_reff_split[i][k]][1])
				gal_r_split[i].append(gal_dict[gal_reff_split[i][k]][2])
				
		gal_r_avg = [np.mean(gal_r_split[i]) for i in range(0, len(gal_reff_split))]
			
		return gal_r_split, gal_RA_split, gal_Dec_split, gal_r_avg, gal_reff_split, gal_r_eff_avg


def degrees_to_pix(r_deg):
    r_arcmin = 60.0*r_deg #degrees to arcmins
    pix_scale_x = 0.00014532290052643554 # for act_map_deep5_2.pixScaleX
    pix_scale = (pix_scale_x)*((180/np.pi) *(60)) #I chose a random specific map to use
    r_pix = r_arcmin / pix_scale 
    return r_pix

def arcmin_to_degrees(r_arcmin):
	'''inputs a numpy array of radii in arcminutes and outputs that list in degrees'''
	r_degrees = r_arcmin/60.0
	return r_degrees
def arcsec_to_arcmin(r_arcsec):
    '''input the list of radii in arc seconds and outputs that list in arc minutes'''
    r_arcmin = r_arcsec/60.0
    return r_arcmin
