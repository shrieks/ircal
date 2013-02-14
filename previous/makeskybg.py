import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline # spline interpolator


# Zero point fluxes for various filters (i.e. flux of mag 0 object)
# from ngao_bkg_filters.dat
Jzp  = 3.015e+09 # photons/(s m^2 arcsec^2 nm)? 
Hzp  = 2.785e+09
Kszp = 1.506e+09 # Ks
Kzp  = 1.505e+09 # K
Kpzp = 1.661E+09 # Kprime

# From Connie notes
Jzpc  = 3.031e+09 # photons/(s m^2)? 
Hzpc  = 2.817e+09
Kszpc = 1.5066e+09 # Ks
Kpzpc = 1.658e+09  # Kprime
Kzpc  = 1.5062e+09 # K

# zeropoint array is currently J, H, Kprime
zeropoint = np.array([Jzp, Hzp, Kpzp])

# read in mauna kea nearir sky background in 900-5600 nm range
# 1.6mm H2O vapor column, 1.5 airmass
# mksky = np.loadtxt('skybgfiles/nearIR_skybg_16_15_r5.dat')
# 1.6 mm H2O vapor, 1.0 airmass
mksky = np.loadtxt('skybgfiles/mk_skybg_zm_16_10_ph.dat')
lambdas = mksky[:,0] # sampling of 0.02nm and a resolution of 0.04nm
skyflux = mksky[:,1] # ph/sec/arcsec^2/nm/m^2

fileroot='./emissdatafiles'
filters = ['Jband', 'Hband', 'Kprime']
filterflux = np.zeros(len(filters), dtype='float64')
filtermags = np.zeros(len(filters), dtype='float64')

# Get flux in the filter by multiplying sky flux by spline interpolated
# filter transmission vector, but interpolation can only be done for 
# wavelengths less than maximum in *band.txt file
for i, filt in enumerate(filters):
    # read in transmission vector for filter
    filename = fileroot + '/' + filt + '.txt'
    print filename
    vec = np.loadtxt(filename)
    wavevec = vec[:,0]            # first column - wavelengths in um
    wavevec = wavevec * 1000.0    # convert to nm

    transvec = vec[:,1]           # second column - transmissivity in percent
    transvec = transvec/100.0     # convert to fraction

    fspline = InterpolatedUnivariateSpline(wavevec, transvec, k=3)

    # choose interpolation max limit for filter
    print max(wavevec)
    filterlambdas = lambdas[np.where( lambdas < max(wavevec) )]
    filterskyflux = skyflux[np.where( lambdas < max(wavevec) )]
    # interpolate to get transmissivity at each wavelength for skybg
    filtertrans = fspline(filterlambdas)
    # multiply at each wavelength to get flux past filter
    filteroutput  = filterskyflux * filtertrans
    filterflux[i] = filteroutput.sum()

# m1 - m2 = -2.5 log(F1/F2), m2 = 0 here
filtermags = -2.5 * np.log10(filterflux/zeropoint)

print filterflux
print filtermags
