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

#splitting up Deep5 and Deep56 map using submaps, because flipper has trouble when the RA coordinates wrap around

act_map_d5_1, act_map_d5_2 = ACT_split(act_map_d5)



############Forming square maps
#act_map_d6 = act_map_d6.selectSubMap(32.02, 38.0, -1.0, -7.0)
#act_map_d5_1 = act_map_d5_1.selectSubMap(0.1, 2.4, -1.3, 1.0)
#act_map_d5_2 = act_map_d5_2.selectSubMap(354.0, 359.96, -2.98, 2.981)

#Define the maps we are using
maps = [act_map_d5_2, act_map_d6, act_map_d5_1]
#maps = [act_map_d6]
map_names = ["act_map_d5_2", "act_map_d6", "act_map_d5_1"]
#map_names = ["act_map_d6"]

#act_map_d6 = act_map_d6.selectSubMap(32.02, 38.0, -1.0, -7.0)




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


#####Filtering Maps
 
filt_min = 2000.0
filt_max = 5782.0
pix_size = 0.5

D6expanded = np.zeros([1200,1200])
D51expanded = np.zeros([719,719])
D52expanded = np.zeros([839,839])

if len(D6expanded) < len(act_map_d6.data[:]):
    D6expanded = act_map_d6.data[:].copy()
    D6expanded[:len(D6expanded)] += D6expanded
else:
    D6expanded = D6expanded.copy()
    D6expanded[:len(act_map_d6.data[:])] += act_map_d6.data[:]
    
#if len(D51expanded) == len(act_map_d5_1.data[:]):
#    D51expanded = act_map_d5_1.data[:].copy()
#    D51expanded[:len(D51expanded)] += D51expanded
#else:
#    D51expanded = D51expanded.copy()
#    D51expanded[:len(act_map_d5_1.data[:])] += act_map_d5_1.data[:]
    
if len(D52expanded) < len(act_map_d5_2.data[:]):
    D52expanded = act_map_d5_2.data[:].copy()
    D52expanded[:len(D52expanded)] += D52expanded
else:
    D52expanded = D52expanded.copy()
    D52expanded[:len(act_map_d5_2.data[:])] += act_map_d5_2.data[:]

def Hann_window(smap):
    "makes a cosine window for apodizing to avoid edges effects in the 2d FFT" 
    # make a 2d coordinate system
    Nx = smap.shape[0]
    Ny = smap.shape[1]
    inds_x  = (np.arange(Nx)+.5 - Nx/2.)/Nx *np.pi ## eg runs from -pi/2 to pi/2
    #X1 = X[0].reshape(N,1)
    #Y = np.concatenate((X,X1),axis=1)
    inds_y  = (np.arange(Ny)+.5 - Ny/2.)/Ny *np.pi ## eg runs from -pi/2 to pi/2
  
    # make a window map
    window_map = np.outer(np.cos(inds_y), np.cos(inds_x))
    #window_map = np.reshape(window_map, (smap_data.shape[0],smap_data.shape[1] ))
   
    # return the window map
    return(window_map)

D6apodizing_window = Hann_window(D6expanded)
D51apodizing_window = Hann_window(D51expanded)
D52apodizing_window = Hann_window(D52expanded)

D6expanded = D6apodizing_window * D6expanded
#act_map_d5_1.data[:] = D51apodizing_window * act_map_d5_1.data[:]
D52expanded = D52apodizing_window * D52expanded

D6filt = make_2d_filter(D6expanded,pix_size,ell_dd,filt_min, filt_max)
#D51filt = make_2d_filter(act_map_d5_1,pix_size,ell_dd,filt_min, filt_max)
D52filt = make_2d_filter(D52expanded,pix_size,ell_dd,filt_min, filt_max)

D6expanded = generate_filt_map(D6expanded, D6filt)
#act_map_d5_1 = generate_filt_map(act_map_d5_1, D51filt)
D52expanded = generate_filt_map(D52expanded, D52filt)

#act_map_d6.plot()
#act_map_d5_1.plot()
#act_map_d5_2.plot()

#####Filtering Maps




#all stacked maps
#this is so that I have them all in an array ready to rescale and restack into a single graph

all_stacked_smaps = []
all_stacked_rmaps = []

galaxy_bin_temp_rescale = []

#a list of this is also necessary for rescaling stacked galaxies later on in the code

R_filt_list = []

total_galaxies_num = 0

ref_map = maps[0]
Radius = 20
def Stack_on_Positions(map,Nx, Ny, cat,N_objects,Radius):
    stack = np.zeros([2*Radius,2*Radius])
    print map.shape
#    Nx = map.Nx
#    #print Nx
#    Ny = map.Ny
#    #print Ny
    counter = 0
    i = 0
    while (i < N_objects):
        sys.stdout.write("\r galaxy number : %d of %d" % (i,N_objects))
        xc = int(cat[0,i])
        yc = int(cat[1,i])
#        print xc
#        print yc
        
        if ((xc > Radius) and (xc < Nx-Radius)):
            if ((yc > Radius) and (yc < Ny-Radius)):
                stack = stack + map[xc-Radius:xc+ Radius,yc-Radius:yc+Radius]
                counter +=1
                sys.stdout.write("\r galaxy number : %d of %d" % (counter,N_objects))

        i = i + 1
    
    return(stack) , counter


###########################For D6
bound_x = 0.1666667 
bound_y = 0.1666667

Nx = D6expanded.shape[0]
#print Nx
Ny = D6expanded.shape[1]
#print Ny

gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('effective radius', ang_R, gal_RA, gal_dec, eff_R, n)

D6 = D6expanded
#D6 = D6.transpose()

for j in np.arange(0,len(gal_RA)):
    gal_RA_new, gal_dec_new = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], bound_x, bound_y, act_map_d6)
print("Galaxies within map bounds ",len(gal_RA_new))

D6_rand_RA, D6_rand_Dec = gen_rand_bounds(act_map_d6, 5307, bound_x, bound_y)
#5307
D6_rand_RA = np.asarray(D6_rand_RA)
D6_rand_Dec = np.asarray(D6_rand_Dec)

D6RA = -1*(D6_rand_RA - 40.000000) * 60.0 * (1.0/0.501095)
D6Dec = (D6_rand_Dec + 7.799062) * 60.0 * (1.0/0.498239)

catD6 = np.zeros([2,len(D6RA)])
i = 0
while i < len(D6RA):
    catD6[0, i] = D6RA[i]
    catD6[1, i] = D6Dec[i]
    i += 1

D6_objects = len(D6RA)

D6_smap, D6_Counter = Stack_on_Positions(D6, Nx, Ny,catD6,D6_objects,Radius)
#D6_smap = D6_smap.transpose()

plt.figure(1)
from mpl_toolkits.axes_grid1 import make_axes_locatable
print("map mean:",np.mean(D6_smap),"map rms:",np.std(D6_smap))
plt.figure(figsize=[7,7])
im = plt.imshow(D6_smap, origin='lower',cmap=cm.jet)
ax=plt.gca()
divider = make_axes_locatable(ax)
#    	cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im)
num_ticks = 5
xlabels, ylabels = pix_to_deg_ticks(D6_smap,ref_map)
ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
xlabels = xlabels[ticks_num_x] * 60.0
ylabels = ylabels[ticks_num_y] * 60.0
xlabels = np.around(xlabels, decimals=1)
ylabels = np.around(ylabels, decimals=1)
plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
plt.yticks(ticks_num_y, ylabels, rotation='vertical')
plt.title("Stacked Locations for " + str(D6_Counter) + " galaxies")
plt.xlabel("[Arcminutes]")
plt.ylabel("[Arcminutes]")

cbar.set_label('lensing deflection', rotation=270)
plt.show()


#################### D_51
#bound_x = 0.1666667 
#bound_y = 0.1666667
#
#Nx = act_map_d5_1.Nx
##print Nx
#Ny = act_map_d5_1.Ny
##print Ny
#
#gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
#ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('effective radius', ang_R, gal_RA, gal_dec, eff_R, n)
#
#D5_1 = act_map_d5_1.data[:]
#D5_1 = D5_1.transpose()
#
#for j in np.arange(0,len(gal_RA)):
#    gal_RA_new, gal_dec_new = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], bound_x, bound_y, act_map_d5_1)
#print("Galaxies within map bounds ",len(gal_RA_new))
#
#D51_RA = -1*(gal_RA_new - 2.500000) * 60.0 * (1.0/0.498333)
#D51_Dec = (gal_dec_new + 3.001372) * 60.0 * (1.0/0.499532)
#
#catD51 = np.zeros([2,len(D51_RA)])
#i = 0
#while i < len(D51_RA):
#    catD51[0, i] = D51_RA[i]
#    catD51[1, i] = D51_Dec[i]
#    i += 1
#
#D51_objects = len(D51_RA)
#
#D51_smap, D51_Counter = Stack_on_Positions(D5_1, Nx, Ny,catD51,D51_objects,Radius)
##D51_smap = D51_smap.transpose()
#
#plt.figure(2)
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#print("map mean:",np.mean(D51_smap),"map rms:",np.std(D51_smap))
#plt.figure(figsize=[7,7])
#im = plt.imshow(D51_smap, origin='lower',cmap=cm.jet)
#ax=plt.gca()
#divider = make_axes_locatable(ax)
##    	cax = divider.append_axes("right", size="5%", pad=0.05)
#cbar = plt.colorbar(im)
#num_ticks = 5
#xlabels, ylabels = pix_to_deg_ticks(D51_smap,ref_map)
#ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
#ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
#xlabels = xlabels[ticks_num_x] * 60.0
#ylabels = ylabels[ticks_num_y] * 60.0
#xlabels = np.around(xlabels, decimals=1)
#ylabels = np.around(ylabels, decimals=1)
#plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
#plt.yticks(ticks_num_y, ylabels, rotation='vertical')
#plt.title("Stacked Locations for " + str(D51_Counter) + " galaxies")
#plt.xlabel("[Arcminutes]")
#plt.ylabel("[Arcminutes]")
#
#cbar.set_label('lensing deflection', rotation=270)
#plt.show()
#

################### D_52
bound_x = 0.1666667 
bound_y = 0.1666667

Nx = D52expanded.shape[0]
#print Nx
Ny = D52expanded.shape[1]
#print Ny

gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('effective radius', ang_R, gal_RA, gal_dec, eff_R, n)

D5_2 = D52expanded
#D5_2 = D5_2.transpose()

for j in np.arange(0,len(gal_RA)):
    gal_RA_new, gal_dec_new = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], bound_x, bound_y, act_map_d5_2)
print("Galaxies within map bounds ",len(gal_RA_new))

D52_rand_RA, D52_rand_Dec = gen_rand_bounds(act_map_d5_2, 3352, bound_x, bound_y)

D52_rand_RA = np.asarray(D52_rand_RA)
D52_rand_Dec = np.asarray(D52_rand_Dec)


D52_RA = -1*(D52_rand_RA - 360.000000) * 60.0 * (1.0/0.499404)
D52_Dec = (D52_rand_Dec + 3.001372) * 60.0 * (1.0/0.499532)

catD52 = np.zeros([2,len(D52_RA)])
i = 0
while i < len(D52_RA):
    catD52[0, i] = D52_RA[i]
    catD52[1, i] = D52_Dec[i]
    i += 1

D52_objects = len(D52_RA)



D52_smap, D52_Counter = Stack_on_Positions(D5_2, Nx,Ny, catD52, D52_objects,Radius)
#D52_smap = D52_smap.transpose()

plt.figure(3)
from mpl_toolkits.axes_grid1 import make_axes_locatable
print("map mean:",np.mean(D52_smap),"map rms:",np.std(D52_smap))
plt.figure(figsize=[7,7])
im = plt.imshow(D52_smap, origin='lower',cmap=cm.jet)
ax=plt.gca()
divider = make_axes_locatable(ax)
#    	cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im)
num_ticks = 5
xlabels, ylabels = pix_to_deg_ticks(D52_smap,ref_map)
ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
xlabels = xlabels[ticks_num_x] * 60.0
ylabels = ylabels[ticks_num_y] * 60.0
xlabels = np.around(xlabels, decimals=1)
ylabels = np.around(ylabels, decimals=1)
plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
plt.yticks(ticks_num_y, ylabels, rotation='vertical')
plt.title("Stacked Locations for " + str(D52_Counter) + " galaxies")
plt.xlabel("[Arcminutes]")
plt.ylabel("[Arcminutes]")

cbar.set_label('lensing deflection', rotation=270)
plt.show()




###### combining all 3

stacked_map = (D6_smap + D52_smap)
total_counter = D6_Counter  + D52_Counter
stacked_map = -1.0*(stacked_map/total_counter)



levels = 2


plt.figure(4)
from mpl_toolkits.axes_grid1 import make_axes_locatable
print("map mean:",np.mean(stacked_map),"map rms:",np.std(stacked_map))
plt.figure(figsize=[7,7])
im = plt.imshow(stacked_map, origin='lower',cmap=cm.jet, vmin = -0.05, vmax = 0.05)
#plt.contour(stacked_map, levels, colors='k', origin='lower')
ax=plt.gca()
divider = make_axes_locatable(ax)
#    	cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im)
num_ticks = 9
xlabels, ylabels = pix_to_deg_ticks(stacked_map,ref_map)
ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
xlabels = xlabels[ticks_num_x] * 60.0
ylabels = ylabels[ticks_num_y] * 60.0
xlabels = np.around(xlabels, decimals=1)
ylabels = np.around(ylabels, decimals=1)
plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
plt.yticks(ticks_num_y, ylabels, rotation='vertical')
#plt.title("Stacked Locations for " + str(total_counter) + " random coordinates")
plt.xlabel("[Arcminutes]")
plt.ylabel("[Arcminutes]")

cbar.set_label('lensing deflection', rotation=270)
plt.show()