#PSMap
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
#act_map_d5_1 = act_map_d5_1.selectSubMap(0.1, 2.4, -1.3, 1.0)
#act_map_d5_2 = act_map_d5_2.selectSubMap(354.0, 359.96, -2.98, 2.981)
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

Number_of_SZ_Clusters  = 100.
Mean_Amplitude_of_SZ_Clusters = 50.
SZ_beta = 2.0
SZ_Theta_core = 3.0

def SZ_source_component(N,pix_size,Number_of_SZ_Clusters,Mean_Amplitude_of_SZ_Clusters,SZ_beta,SZ_Theta_core,do_plots):
    "makes a realization of a naive SZ map"
    SZMap = np.zeros([N,N])
    SZcat = np.zeros([3,Number_of_SZ_Clusters]) ## catalogue of SZ sources, X, Y, amplitude
    # make a distribution of point sources with varying amplitude
    i = 0.
    while (i < Number_of_SZ_Clusters):
        pix_x = N*np.random.rand() 
        pix_y = N*np.random.rand() 
        pix_amplitude = 10.0
        SZMap[pix_x,pix_y] += pix_amplitude
#        SZcat[0,i] = pix_x
#        SZcat[1,i] = pix_y
        SZcat[0,i] = (pix_x * 0.500754 * (1.0/60.0)) + 32.025000
        SZcat[1,i] = -1.008385 - (pix_y * 0.498085 * (1.0/60.0))
        SZcat[2,i] = pix_amplitude
        
        i = i + 1
    
    # make a beta function
    beta = beta_function(N,pix_size,SZ_beta,SZ_Theta_core)
    
    # convovle the beta function with the point source amplitude to get the SZ map
    FT_beta = np.fft.fft2(np.fft.fftshift(beta))
    FT_SZMap = np.fft.fft2(np.fft.fftshift(SZMap))
    SZMap = np.fft.fftshift(np.real(np.fft.ifft2(FT_beta*FT_SZMap)))
    
    # return the SZ map
    return(SZMap,SZcat)    
  ############################### 

def beta_function(N,pix_size,SZ_beta,SZ_Theta_core):
  # make a beta function
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) * pix_size
    X = np.outer(ones,inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)
    
    beta = (1 + (R/SZ_Theta_core)**2.)**((1-3.*SZ_beta)/2.)

    # return the beta function map
    return(beta)
  ############################### 

#def Exponential_source_component(N,pix_size,Number_of_Sources_EX,Amplitude_of_Sources_EX):
#    "makes a realization of a naive point source map"
#    PSMap = np.zeros([N,N])
#    PScat = np.zeros([3,Number_of_Sources_EX])
#    i = 0.
#    while (i < Number_of_Sources_EX):
#        pix_x = np.random.uniform(250, 500)
#        pix_y = N*np.random.rand(250, 500) 
#        pix_amplitude = 50.0
#        PSMap[pix_x,pix_y] += pix_amplitude
#        PScat[0,i] = (pix_x * 0.500754 * (1.0/60.0)) + 32.025000
#        PScat[1,i] = -1.008385 - (pix_y * 0.498085 * (1.0/60.0))
#        PScat[2,i] = pix_amplitude
#        
#        
#        i = i + 1
#
#    return(PSMap, PScat)



#load in galaxy info from fits files
gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
print gal_RA
print gal_dec
print "number of galaxies before = ", gal_RA.shape
N = 718
pix_size = 0.5
Number_of_Sources_EX = 1200.
Amplitude_of_Sources_EX = 1000.
PSMap, PScat = SZ_source_component(N,pix_size,Number_of_SZ_Clusters,Mean_Amplitude_of_SZ_Clusters,SZ_beta,SZ_Theta_core,True)
#PSMap, PScat = Exponential_source_component(N,pix_size,Number_of_Sources_EX,Amplitude_of_Sources_EX)

gal_RA[:] = 0.0
gal_dec[:] = 0.0
gal_RA = np.append(PScat[0, :], gal_RA)
gal_dec = np.append( PScat[1, :], gal_dec)
print gal_RA
print gal_dec
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
print "number of galaxies before = ",  gal_RA.shape
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

ref_map = act_map_d6.selectSubMap(35.0 - .15555, 35.1666667, -4.1666667, -4.0 + .15555)
ref_map.data[:] = 0.0
ref_map1 = ref_map.copy()
#print ref_map1.data[:]
ref_map2 = ref_map.copy()
#print ref_map2.data[:]
ref_map3 = ref_map.copy()
#print ref_map3.data[:]




num_in_stack = 0

def Plot_CMB_Map(Map_to_Plot,X_width,Y_width):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print("map mean:",np.mean(Map_to_Plot),"map rms:",np.std(Map_to_Plot))
    plt.figure(figsize=(7,7))
    im = plt.imshow(np.real(Map_to_Plot), interpolation='bilinear', origin='lower',cmap=cm.RdBu_r)
    #im.set_clim(c_min,c_max)
    ax=plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(im, cax=cax)
    #cbar = plt.colorbar()
    im.set_extent([0,X_width,0,Y_width])
    plt.ylabel('angle $[^\circ]$')
    plt.xlabel('angle $[^\circ]$')
    cbar.set_label('lensing deflection', rotation=270)
    plt.show()
    return(0)

def Stack_on_Positions(map,N,cat,N_objects,Radius):
    stack = np.zeros([2*Radius,2*Radius])
    counter = 0
    i = 0
    while (i < N_objects):
        sys.stdout.write("\r galaxy number : %d of %d" % (i,N_objects))
        xc = cat[0,i]
        yc = cat[1,i]
        if ((xc > Radius) and (xc < N-Radius)):
            if ((yc > Radius) and (yc < N-Radius)):
                stack += map.data[xc-Radius:xc+Radius,yc-Radius:yc+Radius]
                counter +=1
        i = i + 1
    return(stack/counter) , counter


#def Stack_on_Positions(map,ref_map, RA, Dec,N_objects,Radius):
#    
#    #print stack.shape
#    counter = 0
#    i = 0
#    while (i < N_objects):
#        sys.stdout.write("\r galaxy number : %d of %d" % (i,N_objects))
#        xc = RA[i]
#        yc = Dec[i]
#        smap = map.selectSubMap(xc-Radius, xc+Radius, yc-Radius, yc+Radius)
#        #smap.plot()
#        smap = crop_arraynew(smap, ref_map)
#        ref_map.data[:] += smap.data[:]
#        counter +=1
#        i = i + 1
#        
#    print counter
#    ref_map.data[:] = ref_map.data[:]/counter
#    
#    return ref_map, counter




#    
#    PSMap = np.zeros([N,N])
#    PScat = np.zeros([3,Number_of_Sources_EX])
#    i = 0.
#    while (i < Number_of_Sources_EX):
##        pix_x = np.random.uniform(34.0, 36.0)
##        pix_y = np.random.uniform(-5.0, -3.0)
#        pix_x = N*np.random.rand() 
#        pix_y = N*np.random.rand() 
#        pix_amplitude = 10.0
#        PScat[0,i] = pix_x
#        PScat[1,i] = pix_y
#        PScat[2,i] = pix_amplitude
#        PSMap[pix_x,pix_y] += pix_amplitude
#        
#        i = i + 1
#
#    return(PSMap, PScat)



#pix_size = 0.5
#for i in PScat[0:i]:
#   i = N*pix_size/60
#
#for i in PScat[1:i]:
#    i = N*pix_size/60



bound_x = 0.1667
bound_y = 0.1667
    
filt_min = 2000.0
filt_max = 70000.0
pix_size = 0.5
X_width = N*pix_size/60
Y_width = N*pix_size/60
#apodizing_window_D51 = cosine_window(act_map_d5_1)
#act_map_d5_1.data[:] = apodizing_window_D51 * act_map_d5_1.data[:]
#filt_d51 = make_2d_filter(act_map_d5_1,pix_size,ell_dd,filt_min, filt_max)
#act_map_d5_1 = generate_filt_map(act_map_d5_1, filt_d51)   
#act_map_d5_1.plot()
#
#apodizing_window_D52 = cosine_window(act_map_d5_2)
#act_map_d5_2.data[:] = apodizing_window_D52 * act_map_d5_2.data[:]
#filt_d52 = make_2d_filter(act_map_d5_2,pix_size,ell_dd,filt_min, filt_max)
#act_map_d5_2 = generate_filt_map(act_map_d5_2, filt_d52)
#act_map_d5_2.plot()

D6 = act_map_d6.copy()
D6.data[:]=0.0
D6.data[:] = PSMap
#act_map_d6.data[:] = act_map_d6.data[:]
#apodizing_window_D6 = cosine_window(act_map_d6)
#act_map_d6.data[:] = apodizing_window_D6 * act_map_d6.data[:]
#filt_d6 = make_2d_filter(act_map_d6,pix_size,ell_dd,filt_min, filt_max)
#act_map_d6 = generate_filt_map(act_map_d6, filt_d6)
D6.plot()
#act_map_d6.plot()
#act_map_d6.data[:] = act_map_d6.data[:] + D6.data[:]
#act_map_d6.plot()



#gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
#ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('effective radius', ang_R, gal_RA, gal_dec, eff_R, n)
#       
#for j in np.arange(0,len(gal_RA)):
#    gal_RA_new, gal_dec_new, ang_R_new, eff_R_new  = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], ang_R[j], eff_R[j], bound_x, bound_y, act_map_d5_1)
#
#N_objects = len(gal_RA_new)
#Radius = 0.1667 
#smap_d51, counter = Stack_on_Positions(act_map_d5_1,ref_map1,gal_RA_new,gal_dec_new,N_objects,Radius)
#ref_map1 = smap_d51.copy()
#print ref_map1.data[:]
#num_in_stack += counter

#gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()
#ang_R, gal_RA, gal_dec, ang_R_avg, eff_R, eff_R_avg = radial_bin('effective radius', ang_R, gal_RA, gal_dec, eff_R, n)
#for j in np.arange(0,len(gal_RA)):
#    gal_RA_new, gal_dec_new, ang_R_new, eff_R_new  = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], ang_R[j], eff_R[j], bound_x, bound_y, act_map_d5_2)
#
#N_objects = len(gal_RA_new)
#Radius = 0.1667 
#smap_d52, counter = Stack_on_Positions(act_map_d5_2,ref_map2,gal_RA_new,gal_dec_new,N_objects,Radius)
#ref_map2 = smap_d52.copy()
#print ref_map2.data[:]
#num_in_stack += counter


for j in np.arange(0,len(gal_RA)):
    gal_RA_new, gal_dec_new  = del_galaxies_close_to_submap_bounds(gal_RA[j], gal_dec[j], bound_x, bound_y, act_map_d6)

RA = (gal_RA_new - 32.025000) * 60.0 * (1.0/0.500754)

Dec = ((-1.0 * gal_dec_new) - 1.008385) * 60.0 * (1.0/0.498085)

cat = np.zeros([2,len(RA)])
i = 0
while i < len(RA):
    cat[0, i] = RA[i]
    cat[1, i] = Dec[i]
    i += 1

print "gal_RA_new = ", len(RA)
N_objects = len(RA)
Radius =  20.0
N = D6.Nx
smap_d6, counter = Stack_on_Positions(D6, N,cat,N_objects,Radius)
#ref_map3 = smap_d6.copy()
num_in_stack += counter
#print ref_map3.data[:]
        
#ref_map.data[:] = ref_map1.data[:] + ref_map2.data[:] + ref_map3.data[:]
#ref_map = ref_map3.copy()
        
        #looping through each individual galaxy which falls in the map boundary

        
            
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
print("map mean:",np.mean(smap_d6),"map rms:",np.std(smap_d6))
plt.figure(figsize=[7,7])
im = plt.imshow(smap_d6, origin='lower',cmap=cm.jet)
ax=plt.gca()
divider = make_axes_locatable(ax)
#    	cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im)
num_ticks = 5
xlabels, ylabels = pix_to_deg_ticks(smap_d6,smap_d6)
ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
xlabels = xlabels[ticks_num_x] * 60.0
ylabels = ylabels[ticks_num_y] * 60.0
xlabels = np.around(xlabels, decimals=1)
ylabels = np.around(ylabels, decimals=1)
plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
plt.yticks(ticks_num_y, ylabels, rotation='vertical')
plt.title("Stacked Locations for " + str(num_in_stack) + " galaxies")
plt.xlabel("[Arcminutes]")
plt.ylabel("[Arcminutes]")

cbar.set_label('lensing deflection', rotation=270)
plt.show()




    
    
    #total number of galaxies in the bin stack
