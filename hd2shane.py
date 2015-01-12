import numpy as np
import matplotlib as m
import matplotlib.pyplot as mp
import simshane as ss
import pyfits as pf

from astropy.io import fits
m.interactive(True)

# read in spectrum
hddata, hdhdr = fits.getdata('spectra/stars/hd165459_stisnic_002.fits', header=True)

# column 1 is wavelengths in Angstroms, assign to array and convert to nm
lambdas = hddata.field(0)/10
# column 2 is F_lam in ergs s-1 cm-2 A-1, convert to ergs s-1 m-2 nm-1
hdflux  = hddata.field(1) * 1.0e5

#K-band strehl
kstrehl = 0.8
kss     = str(kstrehl)

# Send star through telescope
hd = np.column_stack((lambdas, hdflux)) # create 2 column array
# 16 x 16 subaps
# J
hdout  = ss.simshane(hd,skymod='cp',df=0.01,odel=0.01,aoT=23, subaps=8, aperrad=35.0)
# H
hdout  = ss.simshane(hd,skymod='cp',df=0.01,odel=0.01,aoT=23, subaps=8, aperrad=20.0)
# Ks
hdout  = ss.simshane(hd,skymod='cp',df=0.01,odel=0.01,aoT=23, subaps=8, aperrad=15.0)
#hdout8 = ss.simshane(hd,skymod='cp',df=0.03,odel=0.01,subaps=8,aoT=18.0)
