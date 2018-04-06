from astropy.io import fits
import numpy as np

#code to open fits file and sort through columns for use in code
#fits file is a new galaxy catalogue

#D6Cat = '/Users/alibinesh/SURP/Ali/Catalogues/D6.fits'
#D5Cat = '/Users/alibinesh/SURP/Ali/Catalogues/D5.fits'




def get_gal_info_from_fits():
	'''
	Get galaxy info from the fits file, instead of text file.
	'''

	gal_catalogue = fits.open('/Users/alibinesh/SURP/Ali/Catalogues/combined_full.fits')
    
	gal_header = gal_catalogue[1].header

	gal_data = gal_catalogue[1].data


	gal_RA = gal_data['ra']
	gal_dec = gal_data['dec']
	gal_z =  gal_data['z']
	eff_R = gal_data['expRad_i']
	ang_R = gal_data['petroRad_i']

	#switching R and theta because they seem to be swapped
	#I can take this away if I find out it's actually incorrect
	#gal_theta, gal_R = np.copy(gal_R), np.copy(gal_theta)

	return gal_RA, gal_dec, gal_z, eff_R, ang_R

#gal_RA, gal_dec, gal_z, eff_R, ang_R = get_gal_info_from_fits()

#	eff_R = gal_data['deVRad_u']
#	ang_R = gal_data['expRad_u']