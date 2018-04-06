#Stacking yo
import pyfits
import numpy as np
import os
from matplotlib import cm
import sys
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
ell_dd, C_phi = np.loadtxt("/Users/alibinesh/SURP/CMBAnalysis_SummerSchool/camb_49411559_scalcls.dat", usecols=(0, 4), unpack=True) 


act_map_d5, act_map_d6 = ACT_upload()
act_map_d6 = act_map_d6.selectSubMap(32.02, 38.0, -1.0, -7.0)

#splitting up Deep5 and Deep56 map using submaps, because flipper has trouble when the RA coordinates wrap around

act_map_d5_1, act_map_d5_2 = ACT_split(act_map_d5)
act_map_d5_1 = act_map_d5_1.selectSubMap(0.1, 2.4, -1.3, 1.0)
act_map_d5_2 = act_map_d5_2.selectSubMap(354.0, 359.96, -2.98, 2.981)
act_map_d5_2.plot()
act_map_d6.plot()
act_map_d5_1.plot()

#Define the maps we are using
maps = [act_map_d5_2, act_map_d6, act_map_d5_1]
map_names = ["act_map_d5_2", "act_map_d6", "act_map_d5_1"]

map_dict = dict(zip(maps,map_names))

#pixel scale
delta_xi = 0.501097/60 #degrees
delta_yi = 0.498237/60 #degrees

#load in galaxy info from fits files
gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
ang_R = arcsec_to_arcmin(ang_R)
eff_R = arcsec_to_arcmin(eff_R)
#ang_R = arcmin_to_degrees(ang_R)
#eff_R = arcmin_to_degrees(eff_R)
#print ang_R
#print eff_R

#change all gal_RA to positive ones
gal_RA = neg_gal_RA(gal_RA)

print "Total length of the catalogue: ", len(gal_RA)

print "The galaxy effective radius range of the catalogue is: ", np.amin(eff_R),"-", np.amax(eff_R)

#number of radial bins

n = 1

#if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
#if the argument is 'angular radius', we'll bin by angular radius
#this splits up the galaxies into n lists, index containing a list of the galaxy details for that radial bin!

ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('effective radius', ang_R, gal_RA, gal_dec, eff_R, n)


print "galaxy effective radius average before loop: ", eff_R_avg



#all stacked maps
#this is so that I have them all in an array ready to rescale and restack into a single graph

all_stacked_smaps = []
all_stacked_rmaps = []

galaxy_bin_temp_rescale = []

#a list of this is also necessary for rescaling stacked galaxies later on in the code

R_filt_list = []

total_galaxies_num = 0

ref_map = act_map_d6.selectSubMap(35.0 - .155, 35.1666667, -4.1666667, -4.0 + .155)
ref_map.data[:] = 0.0

def Stack_on_Positions(map,RA, Dec,N_objects,Radius):
    stack = ref_map
    #print stack.shape
    counter = 0
    i = 0
    while (i < N_objects):
        sys.stdout.write("\r galaxy number : %d of %d" % (i,N_objects))
        xc = RA[i]
        yc = Dec[i]
        smap = map.selectSubMap(xc-Radius, xc+Radius, yc-Radius, yc+Radius)
        smap = crop_arraynew(smap, stack)
        stack.data[:] += smap.data[:]
        counter +=1
        i = i + 1
        
    print counter
    stack.data[:] = stack.data[:]/counter
    
    return stack, counter


#looping through each radial bin
for j in np.arange(0,len(gal_RA)):
    print "                    "
    print "############# Stacking for Bin ", j+1
    print "                    "
    
    num_in_stack = 0 #starting the count
    
    
    R = eff_R_avg[j]
    
    R_filt_list.append(R)
    
    #R = 0.5 
    
    #print "The filter radius is ", R
    
    bound_x = 0.1666667 
    bound_y = 0.1666667
    print "bound_x", bound_x
    print "bound_y", bound_y
    
    

    #will loop through all of the maps individually
    
    t1=time.time()
    
    #now looping through each of the CMB maps so we can stack
    all_stacked_smaps = np.zeros([40,40])
    
    
    for orig_map in maps:
        print "Map ", map_dict[orig_map]

        #filter the whole map
        #check to make sure that the filter is what you want!
        #map = generate_filtsmap(orig_map) 
        filt_min = 2000.0
        filt_max = 5782.0
        pix_size = 0.5
        apodizing_window = cosine_window(orig_map)
        orig_map.data[:] = apodizing_window * orig_map.data[:]
        filt = make_2d_filter(orig_map,pix_size,ell_dd,filt_min, filt_max)
        map = generate_filt_map(orig_map, filt)
        map.plot()
        #map = orig_map
        
        #get rid of galaxy locations who are located less than the sub map bound
        #distance from the corner
        
        print("Total Number of Galaxies ", len(gal_RA[j]))
        gal_RA_new, gal_dec_new, ang_R_new, eff_R_new  = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], ang_R[j], eff_R[j], bound_x, bound_y, map)
        print("Galaxies within map bounds ",len(gal_RA_new))
        
        
        #looping through each individual galaxy which falls in the map boundary
        N_objects = len(gal_RA_new)
        Radius = 0.1666667 
        smap, counter = Stack_on_Positions(map,gal_RA_new,gal_dec_new,N_objects,Radius)
        ref_map.data[:] += smap.data[:]
        
            
    #print "Number of random maps after crop" , len(rand_RA)    



        
        
        #plotting the stacked maps per bin!
        
    	#tempmap.plot(useImagePlot=False, colBarOrient='vertical', colBarShrink=0.7)

   		 
#    	fig,ax = plt.subplots(1)
#    	ax.set_aspect('equal')
#    	ax.imshow(tempmap.data[:], origin = 'lower')
#
#    	
#    	cax = ax.imshow(tempmap.data[:], origin = 'lower', cmap=cm.jet)
#    	#Max = plt.scatter([15], [17])
#    	cbar = fig.colorbar(cax)
#    	fig.canvas.draw()
#    	image_size_x = tempmap.data[:].shape[1]
#    	image_size_y = tempmap.data[:].shape[0]
#    	num_ticks = 5
#        xlabels, ylabels = pix_to_deg_ticks(tempmap.data[:],ref_map)
#        ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
#        ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
#        xlabels = xlabels[ticks_num_x] * 60.0
#        ylabels = ylabels[ticks_num_y] * 60.0
#        xlabels = np.around(xlabels, decimals=1)
#        ylabels = np.around(ylabels, decimals=1)
#        plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
#        plt.yticks(ticks_num_y, ylabels, rotation='vertical')
#    	plt.title("Stacked Locations for " + str(num_in_stack) + " galaxies")
#    	plt.xlabel("[Arcminutes]")
#    	plt.ylabel("[Arcminutes]")
    	#plt.savefig(map_dump+'stacked_bin_'+str(j+1))
        
    	from mpl_toolkits.axes_grid1 import make_axes_locatable
    	print("map mean:",np.mean(ref_map.data[:]),"map rms:",np.std(ref_map.data[:]))
    	plt.figure(figsize=[7,7])
    	im = plt.imshow(ref_map.data[:], origin='lower',cmap=cm.jet, vmin = -0.05, vmax = 0.05)
    	ax=plt.gca()
    	divider = make_axes_locatable(ax)
#    	cax = divider.append_axes("right", size="5%", pad=0.05)
    	cbar = plt.colorbar(im)
    	num_ticks = 5
        xlabels, ylabels = pix_to_deg_ticks(ref_map.data[:],ref_map)
        ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
        ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
        xlabels = xlabels[ticks_num_x] * 60.0
        ylabels = ylabels[ticks_num_y] * 60.0
        xlabels = np.around(xlabels, decimals=1)
        ylabels = np.around(ylabels, decimals=1)
        plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
        plt.yticks(ticks_num_y, ylabels, rotation='vertical')
    	plt.title("Stacked Locations for " + str(counter) + " galaxies")
    	plt.xlabel("[Arcminutes]")
    	plt.ylabel("[Arcminutes]")

    	cbar.set_label('lensing deflection', rotation=270)
    	plt.show()



    else:
        continue
    
    
    #total number of galaxies in the bin stack
    print 'The total number of galaxiess stacked in bin '+ str(j+1) + ' is ' + str(counter)
    t2=time.time()
    print "Stacking for this bin took " + str((t2-t1)/60.) + " minutes."
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print("map mean:",np.mean(ref_map.data[:]),"map rms:",np.std(ref_map.data[:]))
    plt.figure(figsize=[7,7])
    im = plt.imshow(ref_map.data[:], origin='lower',cmap=cm.jet, vmin = -0.05, vmax = 0.05)
    ax=plt.gca()
    divider = make_axes_locatable(ax)
#    	cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im)
    num_ticks = 5
    xlabels, ylabels = pix_to_deg_ticks(ref_map.data[:],ref_map)
    ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
    ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
    xlabels = xlabels[ticks_num_x] * 60.0
    ylabels = ylabels[ticks_num_y] * 60.0
    xlabels = np.around(xlabels, decimals=1)
    ylabels = np.around(ylabels, decimals=1)
    plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
    plt.yticks(ticks_num_y, ylabels, rotation='vertical')
    plt.title("Stacked Locations for " + str(counter) + " galaxies")
    plt.xlabel("[Arcminutes]")
    plt.ylabel("[Arcminutes]")

    cbar.set_label('lensing deflection', rotation=270)
    plt.show()

