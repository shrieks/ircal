import pyfits
import math
import numpy as N
from scipy import interpolate
import os.path

def mtham_sky( wave, iphase,tfil= 'lick_sky_d55_2011aug29.fits',tdir = '.'):


#     Sky brightness at Mt Hamilton
#..............................................................................
#

    ## Using an empirical 'dark' sky measurement
    mtham_hdu = pyfits.open(os.path.join(tdir,tfil))
    mtham_swv = mtham_hdu[0].data
    mtham_sfx = mtham_hdu[1].data

    ## Interpolate
    flam = interpolate.splev(wave, interpolate.splrep(mtham_sfx, mtham_swv))
    a = N.where(wave < mtham_swv)
    flam[a] = mtham_sfx[0]

    fnu = flam / (3e10) * wave * (wave * 1e-8)
    msky = -2.5 * alog10(fnu)  - 48.6
    return msky

def maunakea_sky( wave, iphase,  tdir='.', empiri=True,  flg_sky=0):

#    nwv = len(wave)
    msky = zeros_like(wave)

    minwave = wave.min()
    maxwave = wave.max()

    if empiri:
        if flg_sky:
            tfil = 'mkea_sky_newmoon_DEIMOS_1200_2011oct.fits'
        else:
            tfil = 'mkea_sky_newmoon_DEIMOS_600_2011oct.fits'

        skyhdu = pyfits.open(os.path.join(tdir,tfil))
        mkea_swv = skyhdu[0].data
        mkea_sfx = skyhdu[1].data

        emp_wave = N.where(N.logical_and(wave > minwave,wave < maxwave))

        flam = interpolate.splev(wave, interpolate.splrep(mkea_swv, mkea_sfx))
        fnu = flam / (3e10) * wave[emp_wave] * (wave[emp_wave] * 1e-8)

        msky = -2.5 * alog10(fnu)  - 48.6
    else:

        xwave = N.array([3500, 4200, 5500, 6700, 7800, 22000.])
        xphase = N.array([0., 3., 7., 10., 14.])
        xsky = N.array([ [22.4, 23.0, 21.9, 21.2, 19.9, 12.0], 
                         [21.5, 22.4, 21.7, 20.8, 19.9, 12.0],  
                         [19.9, 21.6, 21.4, 20.6, 19.7, 12.0], 
                         [18.5, 20.7, 20.7, 20.3, 19.5, 12.0], 
                         [17.0, 19.5, 20.0, 19.9, 19.2, 12.0] ]) ## Vega values

        xsky[:,0] = xsky[:,0] + 0.71
        xsky[:,1] = xsky[:,1] - 0.11
        xsky[:,3] = xsky[:,3] + 0.199
        xsky[:,4] = xsky[:,4] + 0.454 # AB offsets
        xsky = N.ravel(xsky)
        bxwave = N.tile(xwave,len(xphase))
        bxphase = N.repeat(xphase,len(xwave))
        
        emp_wave =  N.where(N.logical_and(wave > minwave,wave < maxwave))

        inphase = N.repeat(iphase,len(emp_wave))
        bxp=N.column_stack((bxwave,bxphase))
        msky= interpolate.griddata(bxp,xsky,(emp_wave, inphase), method='cubic')


    return(msky)
