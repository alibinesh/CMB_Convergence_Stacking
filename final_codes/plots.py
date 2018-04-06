import pyfits
import numpy as np
import os
from matplotlib import cm

#import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords
from matplotlib.patches import Circle

from flipper import *



import matplotlib.ticker as ticker


def Plot_CMB_Map(Map_to_Plot,X_width,Y_width):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print("map mean:",np.mean(Map_to_Plot),"map rms:",np.std(Map_to_Plot))
    plt.figure(figsize=(7,7))
    im = plt.imshow(np.real(Map_to_Plot), interpolation='bilinear', origin='lower',cmap=cm.jet)
    #im.set_clim(c_min,c_max)
    ax=plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(im, cax=cax)
    #cbar = plt.colorbar()
    im.set_extent([-X_width,X_width,-Y_width,Y_width])
    plt.ylabel('angle $[^\circ]$')
    plt.xlabel('angle $[^\circ]$')
    cbar.set_label('lensing deflection', rotation=270)
    plt.show()
    return(0)

def calculate_2d_spectrum(Map1,Map2,N):
    "calcualtes the power spectrum of a 2d map by FFTing, squaring, and azimuthally averaging"
    
    # make a 2d ell coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
    kX = np.outer(ones,inds) / (0.50/60. * np.pi/180.)
    kY = np.transpose(kX)
    K = np.sqrt(kX**2. + kY**2.)
    ell_scale_factor = 2. * np.pi 
    ell2d = K * ell_scale_factor
    
    # make an array to hold the power spectrum results
    N_bins = int(10000.0/50.0)
    ell_array = np.arange(N_bins)
    CL_array = np.zeros(N_bins)
    
    # get the 2d fourier transform of the map
    FMap1 = np.fft.ifft2(np.fft.fftshift(Map1))
    FMap2 = np.fft.ifft2(np.fft.fftshift(Map2))
    PSMap = np.fft.fftshift(np.real(np.conj(FMap1) * FMap2))
    # fill out the spectra
    i = 0
    while (i < N_bins):
        ell_array[i] = (i + 0.5) * 50.0
        inds_in_bin = ((ell2d >= (i* 5.0)) * (ell2d < ((i+1)* 5.0))).nonzero()
        CL_array[i] = np.mean(PSMap[inds_in_bin])
        #print i, ell_array[i], inds_in_bin, CL_array[i]
        i = i + 1
 
    # return the power spectrum and ell bins
    return(ell_array,CL_array*np.sqrt(pix_size /60.* np.pi/180.)*2.)


def degrees_to_pix(r_deg,map):
    r_arcmin = 60*r_deg #degrees to arcmins
    #pix_scale_x = 0.00014532290052643554 # for act_map_deep5_2.pixScaleX, change this to the maps you're using map.pixScaleX
    pix_scale_x = map.pixScaleX
    pix_scale = (pix_scale_x)*((180/np.pi) *(60)) 
    r_pix = r_arcmin / pix_scale 
    return r_pix

def pix_to_deg(r_pix,map):
	#pix_scale_x = 0.00014532290052643554
	pix_scale_x = map.pixScaleX
	pix_scale = (pix_scale_x)*((180/np.pi) *(60)) 
	r_arcmin = r_pix*pix_scale
	r_deg = r_arcmin/60.
	return r_deg
    

def radial_pix_ticks(stacked_image):
	y_shape, x_shape = stacked_image.shape
	pix_ticks_x = np.arange(0,x_shape)
	pix_ticks_y = np.arange(0,y_shape)
	
	#find the x any y middle pixel
	
	middle_pix_x = np.round(x_shape/2.)
	middle_pix_y = np.round(y_shape/2.)
	
	#dealing with x only right now
	
	second_half_x = pix_ticks_x[np.where(pix_ticks_x < middle_pix_x)]
	
	new_second_half_x = np.arange(0, len(second_half_x))
	
	reversed_new_first_half_x = np.arange(1, len(second_half_x))
	
	new_first_half_x = reversed_new_first_half_x[::-1]
	
	#putting it all together
	
	new_pix_ticks_x = np.concatenate((new_first_half_x,new_second_half_x),axis = 0)

	#dealing with y now
	
	second_half_y = pix_ticks_x[np.where(pix_ticks_y < middle_pix_y)]
	
	new_second_half_y = np.arange(0, len(second_half_y))
	
	reversed_new_first_half_y = np.arange(1, len(second_half_y))
	
	new_first_half_y = reversed_new_first_half_y[::-1]
	
	#putting it all together
	
	new_pix_ticks_y = np.concatenate((new_first_half_y,new_second_half_y),axis = 0)
	
	return new_pix_ticks_x, new_pix_ticks_y

def pix_to_deg_ticks(stacked_image, ref_map):
	new_pix_ticks_x, new_pix_ticks_y = radial_pix_ticks(stacked_image)
	
	deg_ticks_x = np.zeros(len(new_pix_ticks_x))
	deg_ticks_y = np.zeros(len(new_pix_ticks_y))
	#dealing with x:
	for i in np.arange(0,len(new_pix_ticks_x)):
		deg_ticks_x[i] = pix_to_deg(new_pix_ticks_x[i],ref_map)
	for j in np.arange(0,len(new_pix_ticks_y)):
		deg_ticks_y[j] = pix_to_deg(new_pix_ticks_y[j], ref_map)
	
	return deg_ticks_x, deg_ticks_y

def rescaled_stacked_map_per_bin_plot(stacked_map_bin,rescaled_map,z,R_filt,R_filt_ref,map_dump, ref_map):
    '''
    Plots both the rescaled map and the original map to see the pixel number difference
    all_stacked_maps[z] = stack_map_bin
    '''
    
    x_original = stacked_map_bin.shape[0]
    y_original = stacked_map_bin.shape[1]
    x_rescaled = rescaled_map.shape[0]
    y_rescaled = rescaled_map.shape[1]
    
    rescale_param = R_filt_ref/R_filt
    
    galaxy_cent_x = stacked_map_bin.shape[0]/2.
    galaxy_cent_y = stacked_map_bin.shape[1]/2.
    
    rescaled_galaxy_cent_x = rescaled_map.shape[0]/2.
    rescaled_galaxy_cent_y = rescaled_map.shape[1]/2.
    
    fig = plt.figure()
    fig.subplots_adjust(hspace=.5)
    ax = fig.add_subplot(211)
    
    #original stacked map
    
    ax.imshow(stacked_map_bin, origin = 'lower', cmap=cm.jet)
    #circ= Circle((1,3), 1,color='k', linewidth = 2, fill=False)
    #circle1 = Circle((galaxy_cent_x,galaxy_cent_y), radius=degrees_to_pix(R_filt, ref_map), color='k', linewidth = 2, fill=False)
    #circle2 = Circle((galaxy_cent_x,galaxy_cent_y), radius=np.sqrt(2)*degrees_to_pix(R_filt, ref_map), color='k', linewidth = 2, fill=False)
    #ax.add_patch(circle1)
    #ax.add_patch(circle2)
    #ax.plot(galaxy_cent_x,galaxy_cent_y,'wo')
    plt.xlim(0,x_original)
    plt.ylim(0,y_original)
    plt.title("Original galaxy Stacked Map, Bin "+str(z))
    print "                   "
    print "                   "
    print "                   "
    
    #rescaled stacked map
    
    ax1 = fig.add_subplot(212)
    ax1.imshow(rescaled_map, origin = 'lower', cmap=cm.jet)
    #circ= Circle((2,3), 1,color='k', linewidth = 2, fill=False)
    #circle1 = Circle((rescaled_galaxy_cent_x,rescaled_galaxy_cent_y), radius=rescale_param*degrees_to_pix(R_filt, ref_map), color='k', linewidth = 2, fill=False)
    #circle2 = Circle((rescaled_galaxy_cent_x,rescaled_galaxy_cent_y), radius=rescale_param*np.sqrt(2)*degrees_to_pix(R_filt, ref_map), color='k', linewidth = 2, fill=False)
    #ax1.add_patch(circle1)
    #ax1.add_patch(circle2)
    #ax1.plot(rescaled_galaxy_cent_x,rescaled_galaxy_cent_y,'wo')
    plt.xlim(0,x_rescaled)
    plt.ylim(0,y_rescaled)
    plt.title("Rescaled galaxy Stacked Map, Bin "+str(z))
    fig.savefig(map_dump + "rescale_orig_bin_"+str(z)+"_galaxy.png")
    return

def total_rescaled_stacked_map_plot2(total_rescaled_stack,map_dump, total_galaxys_num,ref_map):
    galaxy_cent_x = total_rescaled_stack.shape[0]/2.0
    galaxy_cent_y = total_rescaled_stack.shape[1]/2.0
    fig,ax = plt.subplots(1, figsize=(10,10))
    ax.set_aspect('equal')
    ax.imshow(total_rescaled_stack, origin = 'lower')
    #circle1 = Circle((galaxy_cent_x,galaxy_cent_y), radius=degrees_to_pix(R,ref_map), color='k', linewidth = 2, fill=False)
    #circle2 = Circle((galaxy_cent_x,galaxy_cent_y), radius=np.sqrt(2)*degrees_to_pix(R,ref_map), color='k', linewidth = 2, fill=False)
    #ax.text(total_rescaled_stack.shape[0]*(1./10.), total_rescaled_stack.shape[1]*(1./10.), r'$\Delta$T = ' +str(rescaled_dT) + r' $\mu$K', style='italic',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    #ax.add_patch(circle1)
    #ax.add_patch(circle2)
    cax = ax.imshow(total_rescaled_stack, origin = 'lower', cmap=cm.jet)
    cbar = fig.colorbar(cax)
    fig.canvas.draw()
    image_size_x = total_rescaled_stack.shape[1]
    image_size_y = total_rescaled_stack.shape[0]
    num_ticks = 5
    xlabels, ylabels = pix_to_deg_ticks(total_rescaled_stack,ref_map)
    ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
    ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
    xlabels = xlabels[ticks_num_x] * 60.0
    ylabels = ylabels[ticks_num_y] * 60.0
    xlabels = np.around(xlabels, decimals=1)
    ylabels = np.around(ylabels, decimals=1)
    plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
    plt.yticks(ticks_num_y, ylabels, rotation='vertical')
    #plt.title("Rescaled Stacked Map for " + str(total_galaxys_num) + " galaxies")
    plt.xlabel("[Arcminutes]")
    plt.ylabel("[Arcminutes]")
    plt.savefig(map_dump+'rescaled_stack_galaxys.png')

def total_rescaled_stacked_random_map_plot2(total_rescaled_stack,map_dump, total_galaxys_num,ref_map):
    galaxy_cent_x = total_rescaled_stack.shape[0]/2.0
    galaxy_cent_y = total_rescaled_stack.shape[1]/2.0
    fig,ax = plt.subplots(1, figsize=(10,10))
    ax.set_aspect('equal')
    ax.imshow(total_rescaled_stack, origin = 'lower')
    #circle1 = Circle((galaxy_cent_x,galaxy_cent_y), radius=degrees_to_pix(R,ref_map), color='k', linewidth = 2, fill=False)
    #circle2 = Circle((galaxy_cent_x,galaxy_cent_y), radius=np.sqrt(2)*degrees_to_pix(R,ref_map), color='k', linewidth = 2, fill=False)
    #ax.text(total_rescaled_stack.shape[0]*(1./10.), total_rescaled_stack.shape[1]*(1./10.), r'$\Delta$T = ' +str(rescaled_dT) + r' $\mu$K', style='italic',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    #ax.add_patch(circle1)
    #ax.add_patch(circle2)
    cax = ax.imshow(total_rescaled_stack, origin = 'lower', cmap=cm.jet)
    cbar = fig.colorbar(cax)
    fig.canvas.draw()
    image_size_x = total_rescaled_stack.shape[1]
    image_size_y = total_rescaled_stack.shape[0]
    num_ticks = 5
    xlabels, ylabels = pix_to_deg_ticks(total_rescaled_stack,ref_map)
    ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
    ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
    xlabels = xlabels[ticks_num_x] * 60.0
    ylabels = ylabels[ticks_num_y] * 60.0
    xlabels = np.around(xlabels, decimals=1)
    ylabels = np.around(ylabels, decimals=1)
    plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
    plt.yticks(ticks_num_y, ylabels, rotation='vertical')
    #plt.title("Rescaled Stacked Map for " + str(total_galaxys_num) + " galaxies")
    plt.xlabel("[Arcminutes]")
    plt.ylabel("[Arcminutes]")
    plt.savefig(map_dump+'test_map.png')






