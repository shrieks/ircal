import numpy as np
import matplotlib.pyplot as mp
import pyfits as pf
#import simIrcal as si
import simshane as ss
from astropy.io import fits
#from itertools import cycle

# Read in calspec star flux file
# The unit of flux in all files is erg s-1 cm-2 A-1.
# wavelengths in Angstroms
#bdfits = pf.getdata('spectra/stars/bd26d2606_stis_002.fits',0)
bdfits, bdhdr = fits.getdata('spectra/stars/bd26d2606_stis_002.fits', header=True)
#bdfits, bdhdr = fits.getdata('spectra/stars/bd21d0607_stis_003.fits', header=True)
#bd.names
# ['WAVELENGTH', 'FLUX', 'STATERROR', 'SYSERROR', 'FWHM', 'DATAQUAL', 'TOTEXP']
waves = bdfits.field(0) / 10.0 # convert to nm
# convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
flux = bdfits.field(1) * 1.0e5

bd = np.column_stack([waves, flux])

#vega = np.loadtxt('spectra/alpha_lyr_stis_005.txt')
# convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
#vegaflux = vega[:,1]*1.0e5
#vega[:,1] = vegaflux

#vegalamb = vega[:,0]/10.0   # table has wavelength in angstroms
#vega[:,0] = vegalamb

# K band strehl ratios
#olstr = 0.6
#ols = str(olstr)
#shastr = 0.8
#shas = str(shastr)
# Send bd & vega through the telescope
temp = 10.6
dufr = 0.01
bd60 = ss.simshane(bd,skymod='cp',df=dufr,odel=0.01,aoT=temp,aperrad=60.0)
bd80 = ss.simshane(bd,skymod='cp',df=dufr,odel=0.01,aoT=temp,aperrad=80.0)
bd10 = ss.simshane(bd,skymod='cp',df=dufr,odel=0.01,aoT=temp,aperrad=100.0)
# Angie's darks taken at 20C
#bd20 = ss.simshane(bd,skymod='cp',df=0.02,odel=0.01,aoT=20.0)
