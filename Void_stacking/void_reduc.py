#import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

from flipper import *

'''
These functions contain manipulations to the map and voids with no flipper dependence
it's mainly for cutting out voids to be ready to use for flipper analysis
'''

def void_split(RA_list, Dec_list, r_eff_list, map_1, map_2):
    '''
    Will split up the input void locations into whether they are in which cut of the act map
    This function ONLY works for the original ACT map given by Renee, must change if we want
    more general
    '''
    RA_map1 = []
    Dec_map1 = []
    r_eff_list_1 = []
    RA_map2 = []
    Dec_map2 = []
    r_eff_list_2 = []
    for i in np.arange(0,len(RA_list)):
        #print RA_list[i]
        if RA_list[i] < map_1.x0 and RA_list[i] > map_1.x1:
            if Dec_list[i] < map_1.y1 and Dec_list[i] > map_1.y0:
                #print i
                RA_map1.append(RA_list[i])
                Dec_map1.append(Dec_list[i])
                r_eff_list_1.append(r_eff_list[i])
        if RA_list[i] < map_2.x0 and RA_list[i] > map_2.x1:
            if Dec_list[i] < map_2.y1 and Dec_list[i] > map_2.y0:
                RA_map2.append(RA_list[i])
                Dec_map2.append(Dec_list[i])
                r_eff_list_2.append(r_eff_list[i])
    return RA_map2, Dec_map2, r_eff_list_2, RA_map1, Dec_map1, r_eff_list_1

def filter_voids(m, ra_list, dec_list):
    '''
    filters out voids which do not reside within the desired CMB map
    '''
    x1 = m.x1 #RA upper bound
    x0 = m.x0 #RA lower bound
    y1 = m.y1 #Dec upper bound
    y0 = m.y0 #Dec lower bound
    print(x0,x1,y0,y1)
    new_RA_list = []
    new_Dec_list = []
    for j in np.arange(0, len(ra_list)):
        #if void_RA_cmass1[j] > x0 and void_RA_cmass1[j] < x1:
        if x1 > ra_list[j] > x0:
       	    print j
            #if void_dec_cmass1[j] > y0 and void_dec_cmass1[j] < y1:
            if y1 > dec_list[j] > y0:
            	print j
                new_RA_list.append(ra_list[j])
                new_Dec_list.append(dec_list[j])
    return new_RA_list, new_Dec_list

def gen_rand_points(m, N):
    ''' Generates random points withtin the bounds of the map (m) desired,
    Will generate N random points'''
    a = m.x1
    b = 360
    c = m.y0
    d = m.y1

    #print a, b, c, d
    
    random_RA = []
    random_Dec = []
    for k in np.arange(0,N):
        random_RA.append(random.uniform(a, b))
        random_Dec.append(random.uniform(c, d))
    return random_RA, random_Dec
    
def gen_rand_bounds(m, N, bound_x, bound_y):
    ''' generates radom numbers taking into account void stamp size so that patches
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

def del_voids_close_to_submap_bounds(RA_list, Dec_list, r_ang_list, r_eff_list, bound_x, bound_y, m):
    #def del_voids_close_to_submap_bounds(RA_list, Dec_list, bound_x, bound_y, m):
    ''' Deletes void locations which are too close to the submap bounds
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
    r_eff_new = np.delete(r_eff_list, index_del)
    r_ang_new = np.delete(r_ang_list, index_del)
    return RA_list_new, Dec_list_new, r_ang_new, r_eff_new

def rand_del_voids_close_to_submap_bounds(RA_list, Dec_list, bound_x, bound_y, m):
    #def del_voids_close_to_submap_bounds(RA_list, Dec_list, bound_x, bound_y, m):
    ''' Deletes void locations which are too close to the submap bounds
this gets rid of assertion errors, for random values'''
    index_del = []
    for i in np.arange(0,len(RA_list)):
        if RA_list[i] + bound_x <= m.x1:
            index_del.append(i)
        if RA_list[i] - bound_x >= m.x0:
            index_del.append(i)
        if Dec_list[i] + bound_x >= m.y1:
            index_del.append(i)
        if Dec_list[i] - bound_x <= m.y0:
            index_del.append(i)
    RA_list_new = np.delete(RA_list, index_del)
    Dec_list_new = np.delete(Dec_list, index_del)
    #r_eff_new = np.delete(r_eff_list, index_del)
    return RA_list_new, Dec_list_new #,r_eff_new
    
########################## 
#New functions updated as of September 13th
##########################

'''
The aim of this code is to generate a function which will take a desired bin input
value and split of the voids based on their void radius.  To do this we will
take a histogram of voids (or order them by value) and split up that list
evenly by the desired bin number size.
The output of this function would n lists (n being the input bin size) of void
RA and Dec split up based on their r_v, it would also then compute mean(r_v) of
each of the bins. This mean value would be used as the annulus size for each of the
bins.
'''

def chunks(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

def radial_bin(condition, void_r, void_RA, void_Dec, void_r_eff, n):
    '''note: the input void_RA, void_Dec, and void_r are LISTS'''
    #if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
	#if the argument is 'angular radius', we'll bin by angular radius
	
    if condition == 'angular radius':
    	#meaning we're binning by the angular radius

		#set up a dictionary which related void_index to other void properties
		#dictionary of void_RA, void_Dec, void_r
		void_dict = dict((z[0],list(z[1:])) for z in zip(void_r,void_RA,void_Dec,void_r_eff))
		void_r_sort = sorted(void_r, key=float)
		void_r_split = chunks(void_r_sort, n)
		#now len(void_r_split) should be n
		#empty pre beinned void_RA and Dec lists, based off of the input bin number
		void_RA_split = [[] for _ in range(n)]
		void_Dec_split = [[] for _ in range(n)]
		void_r_eff_split = [[] for _ in range(n)]
		void_r_avg = []

	
		for i in np.arange(0,len(void_r_split)):
			#compute average void_r for each bin
			#this will be used for the R_filt = 0.7*void_r_avg
			void_r_avg.append(np.mean(void_r_split[i]))
			for k in np.arange(0, len(void_r_split[i])):
				#splitting up void_RA based on bin
				void_RA_split[i].append(void_dict[void_r_split[i][k]][0])
				void_Dec_split[i].append(void_dict[void_r_split[i][k]][1])
				void_r_eff_split[i].append(void_dict[void_r_split[i][k]][2])
			
		void_r_eff_avg = [np.mean(void_r_eff_split[i]) for i in range(0, len(void_r_split))]

		return void_r_split, void_RA_split, void_Dec_split, void_r_avg, void_r_eff_split, void_r_eff_avg
	
    elif condition == 'effective radius':
		#set up a dictionary which related void_index to other void properties
		#dictionary of void_RA, void_Dec, void_r
		void_dict = dict((z[0],list(z[1:])) for z in zip(void_r_eff,void_RA,void_Dec,void_r))
		void_reff_sort = sorted(void_r_eff, key=float)
		void_reff_split = chunks(void_reff_sort, n)
		#now len(void_r_split) should be n
		#empty pre beinned void_RA and Dec lists, based off of the input bin number
		void_RA_split = [[] for _ in range(n)]
		void_Dec_split = [[] for _ in range(n)]
		void_r_split = [[] for _ in range(n)]
		void_r_eff_avg = []
		
		for i in np.arange(0,len(void_reff_split)):
			#compute average void_r for each bin
			#this will be used for the R_filt = 0.7*void_r_avg
			void_r_eff_avg.append(np.mean(void_reff_split[i]))
			for k in np.arange(0, len(void_reff_split[i])):
				#splitting up void_RA based on bin
				void_RA_split[i].append(void_dict[void_reff_split[i][k]][0])
				void_Dec_split[i].append(void_dict[void_reff_split[i][k]][1])
				void_r_split[i].append(void_dict[void_reff_split[i][k]][2])
				
		void_r_avg = [np.mean(void_r_split[i]) for i in range(0, len(void_reff_split))]
			
		return void_r_split, void_RA_split, void_Dec_split, void_r_avg, void_reff_split, void_r_eff_avg
		
		
##########################
def neg_void_RA(void_RA):
    '''takes any negative RA values and converts them to positive degrees
    input is a list, returns the same list but with positive values only'''
    for i in np.arange(0,len(void_RA)):
       if void_RA[i] < 0.0: #if negative
            void_RA[i] = 360.0 + void_RA[i]
    return void_RA
    
    
#########################
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
    







            




