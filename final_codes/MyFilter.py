
import numpy as np
import matplotlib
import sys
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import flipper.liteMap as lm

ell_dd, C_phi = np.loadtxt("/Users/alibinesh/SURP/CMBAnalysis_SummerSchool/camb_49411559_scalcls.dat", usecols=(0, 4), unpack=True) 


## variables to set up the size of the map
N = 718.  # this is the number of pixels in a linear dimension
            ## since we are using lots of FFTs this should be a factor of 2^N
pix_size  = 0.5 # size of a pixel in arcminutes

## variables to set up the map plots
c_min = -600.0  # minimum for color bar
c_max = 600.0   # maximum for color bar
X_width = N*pix_size/60.  # horizontal map width in degrees
Y_width = N*pix_size/60.  # vertical map width in degrees

delta_ell = 500.
ell_max = 10000.


def make_CMB_dd_map(N,pix_size,ell_dd,C_dd):
    "makes a realization of a simulated CMB sky map"

    # convert Dl to Cl
    ClTT = C_dd * 2 * np.pi / ((ell_dd*(ell_dd+1.))**2)
    ClTT[0] = 0.
    ClTT[1] = 0.

    # make a 2d coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
    X = np.outer(ones,inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)
    
    # now make a 2d CMB power spectrum
    ell_scale_factor = 2. * np.pi / (pix_size/60. * np.pi/180.)
    ell2d = R * ell_scale_factor
    ClTT_expanded = np.zeros(ell2d.max()+1)
    ClTT_expanded[0:(ClTT.size)] = ClTT
    CLTT2d = ClTT_expanded[ell2d.astype(int)]
    ## make a plot of the 2d cmb power spectrum, note the x and y axis labels need to be fixed
    #Plot_CMB_Map(CLTT2d,ell2d.max(),ell2d.max())  ###
 
    # now make a realization of the CMB with the given power spectrum in fourier space
    ramdomn_array_for_T = np.fft.fft2(np.random.normal(0,1,(N,N)))    
    FT_2d = np.sqrt(CLTT2d) * ramdomn_array_for_T
    ## make a plot of the 2d cmb simulated map in fourier space, note the x and y axis labels need to be fixed
    #Plot_CMB_Map(np.real(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)**2/2/np.pi),0,np.max(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)**2/2/np.pi),ell2d.max(),ell2d.max())  ###
    CMB_dd = np.fft.ifft2(np.fft.fftshift(FT_2d)) /(pix_size /60.* np.pi/180.)
    CMB_dd = np.real(CMB_dd)

    ## return the map
    return(CMB_dd)
  ###############################

def Plot_CMB_Map(Map_to_Plot,X_width,Y_width):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print("map mean:",np.mean(Map_to_Plot),"map rms:",np.std(Map_to_Plot))
    plt.figure(figsize=(5,5))
    im = plt.imshow(np.real(Map_to_Plot), interpolation='bilinear', origin='lower',cmap=cm.jet)
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
  ###############################

## make a CMB T map


def make_2d_filter(N,pix_size,ell_dd,filt_min, filt_max):
    "makes a realization of a simulated CMB sky map"

    # convert Dl to Cl
    Fel = ((1. - np.exp(-ell_dd**2/(2*np.int(filt_min)**2)))*np.exp(-ell_dd**2/(2*np.int(filt_max)**2)))
    Fel = Fel/Fel.max()
    plt.plot(ell_dd, Fel)
    #plt.xlabel('$\ell$')
    

    # make a 2d coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
    X = np.outer(ones,inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)
    
    # now make a 2d CMB power spectrum
    ell_scale_factor = 2. * np.pi / (pix_size/60. * np.pi/180.)
    ell2d = R * ell_scale_factor
    Fel_expanded = np.zeros(ell2d.max()+1)
    print Fel_expanded.size
    print Fel.size
    Fel_expanded[0:(Fel.size)] = Fel
    Fel2d = Fel_expanded[ell2d.astype(int)]
    ## make a plot of the 2d cmb power spectrum, note the x and y axis labels need to be fixed
    #Plot_CMB_Map(Fel2d,ell2d.max(),ell2d.max())  ###
 
    # now make a realization of the CMB with the given power spectrum in fourier space
    FT_2d = Fel2d
    ## make a plot of the 2d cmb simulated map in fourier space, note the x and y axis labels need to be fixed
    #Plot_CMB_Map(np.real(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)**2/2/np.pi),ell2d.max(),ell2d.max())  ###
    filt = np.fft.ifft2(np.fft.fftshift(FT_2d))
    filt = np.real(filt)

    ## return the map
    return(filt)

def Plot_Filtered_Map(Map_to_Plot,X_width,Y_width):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    print("map mean:",np.mean(Map_to_Plot),"map rms:",np.std(Map_to_Plot))
    plt.figure(figsize=[5,5])
    im = plt.imshow(np.real(Map_to_Plot), interpolation='bilinear', origin='lower',cmap=cm.jet)
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

def cosine_window(N):
    "makes a cosine window for apodizing to avoid edges effects in the 2d FFT" 
    # make a 2d coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.)/N *np.pi ## eg runs from -pi/2 to pi/2
    X = np.outer(ones,inds)
    Y = np.transpose(X)
  
    # make a window map
    window_map = np.cos(X) * np.cos(Y)
   
    # return the window map
    return(window_map)


def calculate_2d_spectrum(Map1,Map2,delta_ell,ell_max,pix_size,N):
    "calcualtes the power spectrum of a 2d map by FFTing, squaring, and azimuthally averaging"
    
    # make a 2d ell coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
    kX = np.outer(ones,inds) / (pix_size/60. * np.pi/180.)
    kY = np.transpose(kX)
    K = np.sqrt(kX**2. + kY**2.)
    ell_scale_factor = 2. * np.pi 
    ell2d = K * ell_scale_factor
    print ell2d.shape
    
    # make an array to hold the power spectrum results
    N_bins = int(ell_max/delta_ell)
    ell_array = np.arange(N_bins)
    print ell_array.shape
    #print ell_array
    CL_array = np.zeros(N_bins)
    print CL_array.shape
    #print CL_array
    
    # get the 2d fourier transform of the map
    FMap1 = np.fft.ifft2(np.fft.fftshift(Map1))
    FMap2 = np.fft.ifft2(np.fft.fftshift(Map2))
    PSMap = np.fft.fftshift(np.real(np.conj(FMap1) * FMap2))
    #Plot_CMB_Map(PSMap,ell2d.max(),ell2d.max())
    # fill out the spectra
    i = 0
    while (i < N_bins):
        ell_array[i] = (i + 0.5) * delta_ell
        #print "ell2d", ell2d
        #print "i*delta_ell", i* delta_ell
        inds_in_bin = ((ell2d >= (i* delta_ell)) * (ell2d < ((i+1)* delta_ell))).nonzero()
        CL_array[i] = np.mean(PSMap[inds_in_bin])
        #print CL_array.shape
        i = i + 1
        
    # return the power spectrum and ell bins
    return(ell_array,CL_array*np.sqrt(pix_size /60.* np.pi/180.)*2.)


Deep6 = lm.liteMapFromFits('/Users/alibinesh/SURP/Ali/Data/actpolS2Lensing/actPolS2LensProducts/realKappaCoadd_6.fits')
D6 = Deep6.selectSubMap(36.98,31.0, -7.0, -1.0)
D6.plot()

smap = D6.data[:]
smap_plot = Plot_CMB_Map(smap, X_width, Y_width)

#ell_smap, spectrum_smap = calculate_2d_spectrum(smap, smap, delta_ell, ell_max, pix_size, N)
#
#plt.plot(ell_smap,spectrum_smap, color = 'r')
##plt.plot(ell_dd, C_kk, color = 'b')
##plt.plot(ell_dd, Fel, color = 'y')
#plt.ylabel('ACTPol Patch $C_{kk}$ [$\mu$K$^2$]')
#plt.xlabel('$\ell$')
#plt.show()


filt = make_2d_filter(N, pix_size, ell_dd, 2000.0, 5000.0)

window = (cosine_window(N))
    
window_map = window * smap


filtered_map = np.fft.ifft2(np.fft.fftshift(np.fft.fft2((smap))) * np.fft.fftshift(np.fft.fft2((filt))))

filt_map = (np.fft.ifft2(np.fft.fft2((window_map)) * ((np.fft.fft2((filt))))))

unchanged = np.fft.ifft2(np.fft.fft2((smap)))

#Fel = ((1. - np.exp(-ell_dd**2/(2*6000.0**2)))*np.exp(-ell_dd**2/(2*8000.0**2)))


#plot_filtered_map = Plot_Filtered_Map(filtered_map,X_width,Y_width)
#plot_filtered_map = Plot_Filtered_Map(unchanged,X_width,Y_width)

plot_filt_map = Plot_Filtered_Map(filt_map,X_width,Y_width)

#ell_from_filtmap, spectrum_from_filtmap = calculate_2d_spectrum(filtered_map,filtered_map,delta_ell,ell_max,pix_size,N)
#
#plt.plot(ell_from_filtmap,spectrum_from_filtmap, color = 'r')
#plt.plot(ell_smap, spectrum_smap, color = 'b')
##plt.plot(ell_dd, Fel, color = 'y')
#plt.ylabel('filtered $C_{kk}$ [$\mu$K$^2$]')
#plt.xlabel('$\ell$')
#plt.show()

#ones = np.ones(N)
#inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
#kX = np.outer(ones,inds) / (pix_size/60. * np.pi/180.)
#kY = np.transpose(kX)
#K = np.sqrt(kX**2. + kY**2.)
#ell_scale_factor = 2. * np.pi 
#ell2d = K * ell_scale_factor
#    
## make an array to hold the power spectrum results
#N_bins = int(ell_max/delta_ell)
#ell_array = np.arange(N_bins)
#CL_array = np.zeros(N_bins)
#    
## get the 2d fourier transform of the map
#FMap1 = np.fft.ifft2(np.fft.fftshift(smap))
#PSMap = np.fft.fftshift(np.real(FMap1))
#
#smap_spectrum = Plot_CMB_Map(PSMap, X_width, Y_width)



