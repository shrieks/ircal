import numpy as np
import matplotlib.pyplot as mp
from background import *

def mintime(srcflux, df=0.03, adel=0.0265, odel=0.01,
            source='pts', filts='irc', aoT=18.0):

    print "Running Shane New AO minimum time estimator:"
    if source == 'csb':
        print 'Source type: constant surface brightness patch 1 arcsec^2'
    elif source == 'ext':
        print 'Source type: extended with surface brightness profile'
    else:
        print 'Source type: point source'
    print 'AO temp = {0:.1f}, dust frac = {1:.3f}'.format(aoT,df)
    print 'Al delta = {0:.4f}, Other delta   = {1:.4f}'.format(adel,odel)

    hc = 1.98645e-16      # in cgs - erg cm

    # Set filter curves based on MOSFIRE or Mt Ham filters
    if filts == 'mos':
        filters =  { 'J': MJ, 'H': MH, 'K': MK, 'Ks': MKs}
        filcurve = 'mosfire_JHKKs'
    else:
        filters =  { 'J': J, 'H': H, 'K': K, 'Ks': Ks}
        filcurve = 'JHKpK'

    # Normal usage mode is 16x16 subapertures
    # reference filter for Strehl
    sfil = K
    sref  = 0.8
    print 'Filter set = {0}, Strehl ref = {1}'.format(filts, sref)

    strehl =  { 'J': np.exp(-(sfil['cwvl']/J['cwvl'])**2 * np.log(1/sref)),  
                'H': np.exp(-(sfil['cwvl']/H['cwvl'])**2 * np.log(1/sref)),
                'K': np.exp(-(sfil['cwvl']/K['cwvl'])**2 * np.log(1/sref)),
                'Ks': np.exp(-(sfil['cwvl']/Ks['cwvl'])**2 * np.log(1/sref)),
                'K\'': np.exp(-(sfil['cwvl']/Kp['cwvl'])**2 * np.log(1/sref))}

    # H2RG detector parameters
    ircam = { 'invgain': 2.0,    # electrons per ADU or DN - TBD
              'pscale' : 0.035,  # arcsec pixel-1
              'pixsize': 18.0,   # micron pixel pitch
              'fov'    : 20.0,   # arcsec square
              'fovdia' : 28.0,   # arcsecs - circular aperture equiv on square
              'pixels' : 1024.0, # pixels squared
              'rdnoise': 14.6,   # electrons per CDS read
              'cdstime': 10.6,   # CDS exposure time in seconds
              'rnfowl' : 5.25,   # fowler read noise - 32 samples - e- pixel-1
              'rnfowl16': 6.0,   # fowler 16 read noise
              'fullwell': 1.0e5  # electrons
              }

    # Surfaces in paths/NewShaneAO.txt
    shaneao = ['Aluminum', 'Aluminum',
               'FM3_7', 'FM3_7', 'ALPAO', 'FM3_7', 'TT_DichBS1-refl', 
               'FM3_7', 'MEMS_window', 'MEMS_window', 'BAU', 'MEMS_window', 
               'MEMS_window', 'WFS_DichBS4-trans', 'FM3_7', 'FM3_7', 'CaF', 
               'CaF', 'BAU', 'BAU', filcurve, 'ColdStop', 'BAU', 'BAU', 
               'H2RG_QE']


    shaneaoT= [aoT, aoT, 
               aoT, aoT, aoT, aoT, aoT, aoT, aoT, aoT, 
               aoT, aoT, aoT, aoT, aoT, aoT, aoT, 
               -196.0, -196.0, -196.0, -196.0, -196.0, -196.0, -196.0, -196.0]

    shaneaodel = [adel, adel, 
                  odel, odel, odel, odel, odel, odel, odel, odel,
                  odel, odel, odel, odel, odel, odel, odel, 
                  odel, odel, odel, odel, odel, odel, odel, odel]

    # build telescope data structure to send to telescope sim
    shane =  {'prdia'  : 3.048, # m
              'secdia' : 0.99,  # m - 33 mean dia on Mt Ham site
              'flength': 53.41, # m at Cassegrain focus
              'fratioc': 17.5,  # focal length/diameter at Cass focus
              'fratioa': 28.5,  # at aperture wheel in dewar where slit sits
              'fratiod': 35.0,  # at detector
              'slitwd' : 100.0e-6 # 100 microns 
              }

    # Cerro Pachon sky at 10 mm H2O column and airmass 1.0, 9000' up
    sky  = np.loadtxt('skybgfiles/cp_skybg_zm_100_10_ph.dat')
    skylambdas = sky[:,0]  
    inir = np.where(skylambdas < 2500)
    lambdas = skylambdas[inir]
    dlam   = round(lambdas[1]-lambdas[0],2)

    skyflux = sky[:,1][inir] 
    trans = np.loadtxt('skybgfiles/cptrans_zm_100_10.dat')

    totel,tstel,tetel,thtel = emissivity(lambdas, shaneao, skyflux, shaneaoT,
                                         shaneaodel, dustfrac=df)

    # column 0 of trans is wavelength in microns but it corresponds to
    # lambdas (nm) above - so it can be ignored
    skytrans = trans[:,1][inir]

    # linear interpolate to sky lambdas and conver to photons s-1 m-2 nm-1
    srclint = np.interp(lambdas, srcflux[:,0], srcflux[:,1])*lambdas*100.0/(1.0e9*hc) 
    # mag0 flux at the telescope after propagating through the sky
    # mag0attel = mag0flux * skytrans
    srcattel =  srclint * skytrans

    # Source at detector is source flux at primary multiplied by system
    # transmission vector
    # mag0atout = mag0attel * thtel
    srcatout = srcattel * thtel

    # collecting area = primary area reduced by secondary area
    collarea = np.pi * ((shane['prdia']**2-shane['secdia']**2)/4.0)

    # slit width in arcseconds at detector
    slitangle = (206265*shane['slitwd']*(shane['fratioc']/shane['fratioa'])/
                 shane['flength'])
    # slit angular radius in arcsecs
    slitarad  = slitangle/2.0
    # slit radius in pixels
    slitprad  = slitarad/ircam['pscale']

    # Compute time at which sky background equals read noise in a pixel
    # Sky background per pixel. totel in units ph s-1 m-2 nm-1 arcsec-2
    # pixskybg has units ph s-1 pixel-1
    pixskybg = totel * collarea * dlam * ircam['pscale']**2
    tmin = ircam['rnfowl']**2/pixskybg

    mp.clf()
    mp.plot(lambdas, tmin)
