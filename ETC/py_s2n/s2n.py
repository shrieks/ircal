from math import sqrt
import numpy as N
import pyfits
from scipy import interpolate
import os.path
import instrument
import telescope
import observation

def gaussslit(w, h, xo, yo):
      
      psf  = [1.000,   .995,   .985,   .971,   .954,        
              .933,   .911,   .886,   .860,   .833, 
              .804,   .774,   .743,   .713,   .682, 
              .651,   .620,   .594,   .559,   .529, 
              .500,   .471,   .443,   .417,   .391, 
              .366,   .342,   .319,   .297,   .276, 
              .256,   .237,   .218,   .202,   .187, 
              .172,   .158,   .145,   .132,   .122, 
              .113,   .104,   .097,   .089,   .082, 
              .077,   .072,   .065,   .059,   .057, 
              .052,   .049,   .046,   .042,   .039, 
              .037,   .034,   .032,   .029,   .027, 
              .026,   .024,   .023,   .021,   .019, 
              .018,   .017,   .017,   .016,   .016, 
              .015,   .014,   .013,   .012,   .011,    
              .010,   .010,   .009,   .009,   .008, 
              .008,   .007,   .007,   .006,   .006,    
              .005,   .005,   .005,   .004,   .004, 
              .004,   .004,   .003,   .003,   .003, 
              .003,   .003,   .002,   .002,   .002]

      width = 20*w
      height = 20*h 
      xoff = 40*xo
      yoff = 40*yo

      xin=0.
      xout=0.
      for i in range(1,200) :
            y=float(100-i)
            dy=abs(y-yoff)
            for j in range(1,200) :
                  x=float(j-100)
                  dx=abs(x-xoff)
                  radius=sqrt(x*x+y*y)
                  if radius >= 99. :
                        flux = 0.
                  else :
                        irad=long(radius)
                        drad=radius-long(radius)
                        flux=(1.-drad)*psf[irad]+drad*psf[irad+1]

                  if dy < height and  dx < width :
                        xin=xin+flux
                  if dy > height or  dx > width :
                        xout=xout+flux
                  if (dy == height and  dx <= width) or (dx == width and  dy <= height) :
                        xin=xin+.5*flux
                        xout=xout+.5*flux
          
      

      slit=xin/(xin+xout)
      return slit



      

def getsens(grating,dichroic,sens_dir = '.'):
    if grating == 'G2':
        if instr.dichroic == 'd46':
            sens_fil = 'sens_Kastb600_4310_d55.fits.gz'
        elif instr.dichroic == 'd55':
            sens_fil = 'sens_Kastb600_4310_d55.fits.gz'
    elif grating == 'G3':
        if instr.dichroic == 'd46':
            sens_fil = 'sens_Kastb830_3460_d46.fits.gz'
        elif instr.dichroic == 'd55':
            sens_fil = 'sens_Kastb830_3460_d46.fits.gz'
    elif grating == "600/7500":
        if instr.dichroic == 'd46':
            sens_fil = 'sens_Kastr600_7500_d55.fits.gz'
        elif instr.dichroic == 'd55':
            sens_fil = 'sens_Kastr600_7500_d55.fits.gz'

    senshdu = pyfits.open(os.path.join(sens_dir,sens_fil))

    sens = senshdu[2].data()
    return(sens)

            
def kast_thruput(wave,instr,sens_dir = '.'):
    thru = N.zeros_like(wave,dtype=N.float64)
    bidx = -1
    ridx = -1

    if instr[0].dichroic == 'd46':
        bidx = N.where(wave < 4600)
        ridx = N.where(wave >= 4600)
    elif instr[0].dichroic == 'd55':
        bidx = N.where(wave < 5500)
        ridx = N.where(wave >= 5500)


    sens = getsens(instr[0].grating,instr[0].dichroic,sens_dir)
    thru[bidx] = interpolate.splev(wave[bidx],interpolate.splrep(sens.WAV[0,:],sens.EFF[0,:]))

    sens = getsens(instr[1].grating,instr[1].dichroic,sens_dir)    
    thru[ridx] = interpolate.splev(wave[ridx],interpolate.splrep(sens.WAV[0,:],sens.EFF[0,:]))

    thru[N.where(thru <= 1e-5)] = 1e-5

    return(thru)


def fluxjohnson(wave):

    xwave = [3062.,3312,3562,3812,4062,4212,4462,4712,4962,5425,5925,6425,6925,7750,8350,11050]
    xflux = [569.,586,590,1326,1770,1707,1530,1356,1257,1054,886,749,641,502,435,263]
    wflux = interpolate.splev(wave, interpolate.splrep(xwave, xflux)   )
    return(wflux)


def spec_calcs2n(wave,thru,tel,instr,obs):

    """
    This is an old, old, old piece of code.
      REVISION HISTORY:
     
      ??/??/?? Written by G. Donald Penrod 
      11/13/89 Modified by S. Vogt       - made inputs less confusing
      06/08/92 Modified by M. Keane      - changed from Hamilton to HIRES 
      02/07/96 Modified by C. Churchill  - structured queries by function
                                         - set defaults for Decker C1
                                         - added comments
      20-Oct-2005 Ported to IDL by JXP
      23-Mar-2011 Generatlized for all Keck spectrometers
      23-Aug-2011 Generatlized for all spectrometers

      12/1/2011 - 22+ years later Brad Holden ports this to Python, requires numpy
    """

    height = instr.sheight
    width = instr.swidth
  
    dark = instr.dark
  
    binc = instr.bind
    binr = instr.bins
  
    seeing = obs.seeing
    phase = obs.mphase
    air = obs.airmass
    time = obs.exptime
  
    ## SED
    if obs.mstar < 1:
          mstar = 17.
    else:
          mstar = obs.mstar
    if obs.mtype == 0:
          mtype = 1
    else:
          mtype = obs.mtype
  
    if mtype == 1 :
        nj0     = fluxjohnson(wave)
        n0 = nj0
        
    elif mtype == 2: 
#          nab0    = 3.54e-9 * wave *1e-8 / c.c / c.h  ## photons/s/cm^2/Ang
        nab0    = 10.0**(-0.4*48.6) / 6.626e-27 / wave ## photons/s/cm^2/Ang
        n0 = nab0
        
    else: 
        return(-1)
      
    ## Spreading out the light
    slit0   = x_gaussslit(width/seeing, height/seeing, 0., 0.)
    rows    = long(3*seeing/instr.SCALE_PERP+0.999)
    columns = 2. > ( (width < 3.*seeing)/instr.SCALE_PARA )
    nsky    = float((long(height/instr.SCALE_PERP+0.999) - rows))/rows     
  
    ## Readno
    read = instr.readno * sqrt(rows / float(binr))
  
    #####
    ## S/N
    pixel    = binc*(wave/instr.R)  ## Ang
    slarea   = 3*seeing*width 
  
    ## Grab extinction
    if tel.name == 'KeckI':
          extinct = maunakea_trans(wave)
    elif tel.name == 'KeckII':
          extinct  = maunakea_trans(wave)
    elif tel.name == 'Lick-3m':
          extinct  = mtham_trans(wave)

    slit1 = slit0*10**(-0.4*extinct*air)
    magsky = 0
    ## Sky
    if tel.name == 'KeckI':
          magsky = maunakea_sky(wave, phase, empiri=False)
    elif tel.name =='KeckII': 
          if instr.name == 'DEIMOS': 
              phase = 0 ## Only New Moon so far
              if instr.grating == '1200':
                    magsky = maunakea_sky(wave, phase, empiri=True, flg_sky=1)
              else:
                    magsky = maunakea_sky(wave, phase)
          elif  instr.name =='ESI': 
                phase = 0 ## Only New Moon so far
                magsky = maunakea_sky(wave, phase)
        
          else:
                magsky   = maunakea_sky(wave, phase,  empiri=False)
    elif 'Lick-3m':
          magsky   = mtham_sky(wave, phase)

    ## Object
    star     = n0*(10**(-0.4*mstar))*tel.area*thru*slit1*pixel*time
    dstar    = sqrt(star)
    noise    = read
    #  noise    = read/sqrt(float(binr*binc))
    projslit = columns*(wave/instr.R)
  
    ##  compute the sky counts and noise
    sky   = n0*(10**(-0.4*magsky))*slarea*tel.area*thru*pixel*time
    dsky  = sqrt(sky)
    #
    #  compute the full signal to noise
    ndark  = binc*dark*rows*time/3600.
    ddark  = sqrt(ndark)
    tnoise = sqrt(star+(1.+1./nsky)*(noise*noise+sky+ndark))
    sn     = star/tnoise


  
    fstrct = { 
          sn: sn, 
          star: star, 
          tnoise: tnoise, 
          extinct: extinct, 
          noise: noise, 
          ndark: ndark, 
          slit0: slit0, 
          thru: thru, 
          sky: sky 
          }

    return fstrct

