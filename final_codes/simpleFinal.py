#Stacking yo
import pyfits
import numpy as np
import os

#import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

from flipper import liteMap

# to time how long each step takes
import time


#importing functions from other files
from gal_cat_load import *
from map_filt import *
from map_upload import *
from plots import *
from reductions import *
from stack import *



maps_pwd = '/Users/alibinesh/SURP/Ali/Data/actpolS2Lensing/actPolS2LensProducts/'

code_pwd = '/Users/alibinesh/SURP/Ali/final_codes/'

#Where the plots will be saved 
map_dump = '/Users/alibinesh/SURP/Ali/final_plots/'

#UPLOADING ACT CMB MAPS


act_map_d5, act_map_d6 = ACT_upload()


#splitting up Deep5 and Deep56 map using submaps, because flipper has trouble when the RA coordinates wrap around

act_map_d5_1, act_map_d5_2 = ACT_split(act_map_d5)

#Define the maps we are using
maps = [act_map_d5_2, act_map_d6, act_map_d5_1]
map_names = ["act_map_d5_2", "act_map_d6", "act_map_d5_1"]

map_dict = dict(zip(maps,map_names))

#pixel scale
delta_xi = 0.499532/60 #degrees
delta_yi = 0.498333/60 #degrees

#load in galaxy info from fits files
gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
#ang_R = arcsec_to_arcmin(ang_R)
#eff_R = arcsec_to_arcmin(eff_R)
ang_R = arcmin_to_degrees(ang_R)
eff_R = arcmin_to_degrees(eff_R)
#print ang_R
#print eff_R

#change all gal_RA to positive ones
gal_RA = neg_gal_RA(gal_RA)

print "Total length of the catalogue: ", len(gal_RA)

print "The galaxy effective radius range of the catalogue is: ", np.amin(eff_R),"-", np.amax(eff_R)

#number of radial bins

n = 3

#if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
#if the argument is 'angular radius', we'll bin by angular radius
#this splits up the galaxies into n lists, index containing a list of the galaxy details for that radial bin!

ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('angular radius', ang_R, gal_RA, gal_dec, eff_R, n)


print "galaxy effective radius average before loop: ", eff_R_avg

#keeping track of the delta T lists of all the galaxies in each bin
delt_T_for_bin = [[] for _ in range(n)]

#all stacked maps
#this is so that I have them all in an array ready to rescale and restack into a single graph

all_stacked_maps = []

all_stacked_maps_mod = []

galaxy_bin_temp_rescale = []

#a list of this is also necessary for rescaling stacked galaxies later on in the code

R_filt_list = []

total_galaxies_num = 0

ref_map = maps[0]

#looping through each radial bin
for j in np.arange(0,len(gal_RA)):
    print "                    "
    print "############# Stacking for Bin ", j+1
    print "                    "
    
    num_in_stack = 0 #starting the count
    
    #make sure ang_R and eff_R are in degrees
    
    R = 0.7 * ang_R_avg[j]
    
    R_filt_list.append(R)
    
    #R = 0.5 # choosing constant radius right now
    
    print "The filter radius is ", R
    
    bound_x = 2.0 * R 
    bound_y = 2.0 * R
    
    #empty lists which contains the galaxy stamp maps at different stages
    smaps = []
    filt_smaps = []
    annulus_list = []
    inner_circ_list = []
    mod_smaps = []
    
    
    gal_r_eff_in_bin =[]    
    

    #so we can stack all of the maps together
    #will loop through all of the maps individually
    
    t1=time.time()
    
    stacked_mod_smap = 'hello'
    stacked_smap = 'hello'
    
    #now looping through each of the CMB maps so we can stack
    
    for orig_map in maps:
        print "..... Map ", map_dict[orig_map]

        #filter the whole map
        #check to make sure that the filter is what you want!
        map = generate_filtsmap(orig_map) #Why are we filtering it?
        
        #get rid of galaxy locations who are located less than the sub map bound
        #distance from the corner
        
        print("Total Number of Galaxies ", len(gal_RA[j]))
        gal_RA_new, gal_dec_new, ang_R_new, eff_R_new  = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], ang_R[j], eff_R[j], bound_x, bound_y, map)
        print("Galaxies within map bounds ",len(gal_RA_new))
        
        for a in np.arange(0,len(eff_R_new)):
            gal_r_eff_in_bin.append(eff_R_new[a])
        
        #looping through each individual galaxy which falls in the map boundary
        
        if len(gal_RA_new) !=0:
        	if num_in_stack == 0:
        		print "first map"
        		smap = map.selectSubMap(gal_RA_new[0]-bound_x, gal_RA_new[0] + bound_x, gal_dec_new[0] - bound_y, gal_dec_new[0] + bound_y)
        		smap_data = smap.data[:]
        		stacked_mod_smap = np.zeros(smap_data.shape) #will need to set shape I think
        		stacked_smap = np.zeros(smap_data.shape) #will neep to set shape I think
        		galaxy_bin_temp = 0.
        	
        	for k in np.arange(0,len(gal_RA_new)):
				num_in_stack +=1
				total_galaxies_num += 1
				smap = map.selectSubMap(gal_RA_new[k]-bound_x, gal_RA_new[k] + bound_x, gal_dec_new[k] - bound_y, gal_dec_new[k] + bound_y)
				annulus, inner_circ = delt_T_filt(gal_RA_new[k], gal_dec_new[k], R, smap)
				mod_smap = mod_filt(smap, annulus)
				smap_data = smap.data[:]
				mod_smap_data = mod_smap.data[:]
				arrays = [smap_data,stacked_smap]
				dim_list = min_dim_of_list_of_arrays(arrays)
				#cropping the stacked map
				stacked_smap = crop_array(stacked_smap, dim_list)
				stacked_mod_smap = crop_array(stacked_mod_smap, dim_list)
				#cropping the smaps
				smap_data = crop_array(smap_data, dim_list)
				mod_smap_data = crop_array(mod_smap_data, dim_list)
				stacked_mod_smap += mod_smap_data
				stacked_smap += smap_data
				inner_circ_data = inner_circ.data[:]
				annulus_data = annulus.data[:]
				galaxy_bin_temp += (np.mean(inner_circ_data[inner_circ_data != 0.]) - np.mean(annulus_data[annulus_data != 0.]))
        else:
            continue
        
    if num_in_stack != 0:
        #dividing by counter
        print 'stacking all maps for bin: ', j+1
        stacked_mod_smap *= 1./float(num_in_stack)
        stacked_smap *= 1./float(num_in_stack)
        galaxy_bin_temp *= 1./float(num_in_stack)
        galaxy_bin_temp_rescale.append(galaxy_bin_temp)
        
        all_stacked_maps.append(stacked_smap)
        all_stacked_maps_mod.append(stacked_mod_smap)
        
        
        #plotting the stacked maps per bin!
        
        gal_cent_x = stacked_smap.shape[0]/2.0
    	gal_cent_y = stacked_smap.shape[1]/2.0
   		 
    	fig,ax = plt.subplots(1, figsize=(10,10))
    	ax.set_aspect('equal')
    	ax.imshow(stacked_smap, origin = 'lower')
    	circle1 = Circle((gal_cent_x,gal_cent_y), radius=degrees_to_pix(R, ref_map), color='k', linewidth = 2, fill=False)
    	circle2 = Circle((gal_cent_x,gal_cent_y), radius=np.sqrt(2)*degrees_to_pix(R, ref_map), color='k', linewidth = 2, fill=False)
    	ax.text(stacked_smap.shape[0]*(1./10.), stacked_smap.shape[1]*(1./10.), r'$\Delta$T = ' +str(galaxy_bin_temp) + r' $\mu$K', style='italic',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    	
    	ax.add_patch(circle1)
    	ax.add_patch(circle2)
    	cax = ax.imshow(stacked_smap, origin = 'lower')
    	cbar = fig.colorbar(cax)
    	fig.canvas.draw()
    	image_size_x = stacked_smap.shape[1]
    	image_size_y = stacked_smap.shape[0]
    	num_ticks = 7
        xlabels, ylabels = pix_to_deg_ticks(stacked_smap,ref_map)
        ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
        ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
        xlabels = xlabels[ticks_num_x]
        ylabels = ylabels[ticks_num_y]
        xlabels = np.around(xlabels, decimals=5)
        ylabels = np.around(ylabels, decimals=5)
        plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
        plt.yticks(ticks_num_y, ylabels, rotation='vertical')
    	plt.title("Stacked Galaxy Locations for Bin " + str(j+1)+ " for " + str(num_in_stack) + " galaxies")
    	plt.xlabel("[Arcminutes]")
    	plt.ylabel("[Arcminutes]")
    	plt.savefig(map_dump+'stacked_bin_'+str(j+1))
           
           
    else:
        continue
    
    
    #total number of galaxies in the bin stack
    print 'The total number of galaxiess stacked in bin '+ str(j+1) + ' is ' + str(num_in_stack)
    t2=time.time()
    print "Stacking for this bin took " + str((t2-t1)/60.) + " minutes."


#Now it's time to look at the list with all of the stacked maps
#and rescale them so that we can stack them all in one plot!

print "                            "

print "----------------------------"
print "Rescaling all binned stacked maps......"
print "----------------------------"

rescaled_stacked_maps = []

#corresponds to the largest galaxy bin
#because we want to re scale all smaller galaxies larger, that's just the way the function works

 
R_filt_ref = R_filt_list[0]

#I'm doing the rescale only for the *NON* modified (annulus subtracted) bins
# If you want to change that just change "all_stacked_maps" to "all_stacked_maps_mod" below

for z in np.arange(0, len(all_stacked_maps_mod)):
	#rescaling for individual radial bin 'z'
	#we're rescaling all maps to the smallest sized map
    rescaled_map = stack_rescale(all_stacked_maps_mod[z], R_filt_list[z], R_filt_ref)
    rescaled_stacked_maps.append(rescaled_map)
    #plotting each of the individual rescaled stacked maps:
    rescaled_stacked_map_per_bin_plot(all_stacked_maps_mod[z],rescaled_map,z,R_filt_list[z], R_filt_ref,map_dump, ref_map)

#now we will put together all of the maps into one rescaled stack
total_rescaled_stack =  stack_all_rescaled_maps(rescaled_stacked_maps) 

#getting delta_T for the new rescaled map

#computing the "outer circle" of the rescaled galaxy data array
outer_circ = array_mask(total_rescaled_stack, np.sqrt(2)*R_filt_ref, ref_map)

#computing the "inner circle" of the rescaled galaxy data array
inner_circ = array_mask(total_rescaled_stack, R_filt_ref, ref_map)

#computing the annulus by subtracting the two (to get everything between the inner and the outer circle
annulus =  outer_circ - inner_circ

#annulus subtracting the inner circle
galaxy_bin_temp_rescale = inner_circ - annulus


print "Rescaled delta T: ", galaxy_bin_temp_rescale


total_rescaled_stacked_map_plot2(total_rescaled_stack, R_filt_ref, galaxy_bin_temp_rescale,map_dump, total_galaxies_num, ref_map)


print "DONE CHECK PLOTS FOLDER"
