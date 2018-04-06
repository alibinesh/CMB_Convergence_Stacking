from astropy.io import fits
import numpy as np

#code to open fits file and sort through columns for use in code
#fits file is a new void catalogue

def get_void_info_from_fits():
	'''
	Get void info from the fits file, instead of text file.
	'''

	void_catalogue = fits.open('voids_los-size10-50_nside256_vth21_nlev23.fit')

	void_header = void_catalogue[1].header

	void_data = void_catalogue[1].data

	void_RA = void_data['ra      ']
	void_dec = void_data['dec     ']
	void_z =  void_data['z       ']
	void_rlos = void_data['rlos    ']	
	void_R = void_data['R       ']
	void_theta = void_data['theta   ']
	void_los_size = void_data['los_size']
	void_dens = void_data['dens    ']
	void_f_vol = void_data['f_vol   ']


	#switching R and theta because they seem to be swapped
	#I can take this away if I find out it's actually incorrect
	#void_theta, void_R = np.copy(void_R), np.copy(void_theta)

	return void_RA, void_dec, void_z, void_rlos, void_R, void_theta, void_los_size, void_dens, void_f_vol






