import pyfits
import numpy as np
import os

#import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

from flipper import liteMap


maps_pwd = '/Users/alibinesh/SURP/Ali/Data/actpolS2Lensing/actPolS2LensProducts/'

def ACT_upload():

    #uploading ACT cmb maps

	act_map_d5 = liteMap.liteMapFromFits(maps_pwd+'realKappaCoadd_5.fits')
	act_map_d6 = liteMap.liteMapFromFits(maps_pwd+'realKappaCoadd_6.fits')

	return act_map_d5, act_map_d6


def ACT_split(act_map_d5):

	#splitting up d5 and Deep6 map using submaps

	act_map_d5_1 = act_map_d5.selectSubMap(0.000005, act_map_d5.x0 -0.000005,act_map_d5.y0 + 0.000005,act_map_d5.y1 -0.000005) #RA = 0 to x0
	act_map_d5_2 = act_map_d5.selectSubMap(act_map_d5.x1+ 0.000005,360- 0.000005,act_map_d5.y0 + 0.000005,act_map_d5.y1 - 0.000005) #RA = x1 to 360


	return act_map_d5_1, act_map_d5_2