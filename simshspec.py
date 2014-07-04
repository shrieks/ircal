import numpy as np
from background import *
from resamplesky import * 

def simspec(srcflux, speclam, skymod='cp', df=0.03, adel=0.0265, odel=0.01, 
            subaps=16, source='pts', filts='irc', aoT=18.0, R=500.0, resample=0):
    """
    Simulate Shane AO 3m AO system + instrument emissivity and througput
    Input: 
          srcflux- 2-d array with wavelengths in nm in column 0
                   and flux in units of photons s-1 m-2 nm-1 in col 1
          speclam- wavelengths of interest for spectroscopy SNR
          skymod - which sky model to use
                   cp = Gemini Cerro Pachon 10 mm H2O, airmass 1.0
                   mh = modified Cerro Pachon (adds 10C blackbody)
                   mk = 2 x Mauna Kea flux
          df     - dust fraction (usually 0.0 or 0.008 (0.8%)
          adel   - aluminum reflectivity delta to reduce Al reflectivity
          odel   - delta by which to reduce other surface reflectivity
                   mimics degradation for both
          subaps - number of subapertures. 8 = 8x8, 16 = 16x16
          source - flag to signal input source:
                   pts = point source
                   csb = constant surface brightness - over 1 arcsec^2 area
                         so surface brightness equation reduces to mu = m
                   ext = extended source 
                   for csb and ext flux is in units of
                   photons s-1 m-2 nm-1 arcsec-2
          filts  - 'irc' - for current IRCAL filter set
                   'mos' - for MOSFIRE filter set
          aoT    - AO system temperature. 10C is nominal. 18C is max
                   summer temp
          R      - spectroscopy resolution: lambda/dlambda
    """

    print "--------------"
    print "Running Shane New AO simulation model with parameters:"
    if source == 'csb':
        print 'Source type: constant surface brightness patch 1 arcsec^2'
    elif source == 'ext':
        print 'Source type: extended with surface brightness profile'
    else:
        print 'Source type: point source'
    print 'skymodel = {0}, AO temp = {1:.1f}, dust frac = {2:.3f}'.format(skymod,aoT,df)
    print 'Al delta = {0:.4f}, Other delta   = {1:.4f}, R = {2}'.format(adel,odel, R)

    # ST mag system m = 0 = -2.5log10 F_lam - 21.10
    # zeroptflux = 3.631e-4 # ergs s-1 nm-1 m-2
    hc = 1.98645e-16      # in cgs - erg cm

    # Set filter curves based on MOSFIRE or Mt Ham filters
    if filts == 'mos':
        filters =  { 'J': MJ, 'H': MH, 'Ks': MKs}
        filcurve = 'mosfire_JHKKs'
    else:
        filters =  { 'J': J, 'H': H, 'Ks': Ks}
        filcurve = 'JHKpK'

    # Normal usage mode is 16x16 subapertures
    # reference filter for Strehl
    sfil = K
    if subaps == 8:
        sref = 0.6
    else:
        sref  = 0.8
    print 'Filter set = {0}, Strehl ref = {1}, subapertures={2}'.format(filts, sref, subaps)

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
              'rnfowl' : 5.25,   # fowler read noise - 32 samples
              'rnfowl16': 6.0,   # fowler 16 read noise
              'fullwell': 1.0e5  # electrons
              }

    # Surface list for Shane AO
    # Surfaces : Primary, Secondary,
    #            Fold Mirror 1, OAP 1, DM1, OAP 2	          
    #            TT dichroic/Fold Mirror 2 (LGS/NGS mode), OAP 3, 
    #            MEMs window, MEMs window, MEMs DM, MEMs window, MEMs window, 
    #            NaDichroicSciencePath, OAP 4, Fold mirror 3, IRCAL window
    #          These surfaces at -196C (in Dewar):
    #          Window, Turn Mirror, OAP 1, Filter, ColdStop, OAP 2, 
    #          Turn Mirror, Detector
    #          flat 95% for Cold stop, 100 mm dia extra around aperture

    # Surfaces in paths/NewShaneAO.txt
    shaneao = ['Aluminum', 'Aluminum',
               'FM3_7', 'FM3_7', 'ALPAO', 'FM3_7', 'TT_DichBS1-refl', 
               'FM3_7', 'MEMS_window', 'MEMS_window', 'BAU', 'MEMS_window', 
               'MEMS_window', 'WFS_DichBS4-trans', 'FM3_7', 'FM3_7', 'CaF', 
               'CaF', 'BAU', 'BAU', filcurve, 'Grism', 'ColdStop', 'BAU', 
               'BAU', 'H2RG_QE']


    shaneaoT= [aoT, aoT, 
               aoT, aoT, aoT, aoT, aoT, aoT, aoT, aoT, 
               aoT, aoT, aoT, aoT, aoT, aoT, aoT, 
               -196.0, -196.0, -196.0, -196.0, -196.0, 
               -196.0, -196.0, -196.0, -196.0]
    shaneaodel = [adel, adel, 
                  odel, odel, odel, odel, odel, odel, odel, odel,
                  odel, odel, odel, odel, odel, odel, odel, 
                  odel, odel, odel, odel, odel, odel, odel, odel, odel]

    # build telescope data structure to send to telescope sim
    shane =  {'prdia'  : 3.048, # m
              'secdia' : 0.99,  # m - mean dia on Mt Ham site
              'flength': 53.41, # m at Cassegrain focus
              'fratioc': 17.5,  # focal length/diameter at Cass focus
              'fratioa': 28.5,  # at aperture wheel in dewar where slit sits
              'fratiod': 35.0,  # at detector
              'slitwd' : 100.0e-6 # 100 microns 
              }

    # Sky background data
    if skymod == 'cp':
        # Cerro Pachon sky at 10 mm H2O column and airmass 1.0, 9000' up
        sky  = np.loadtxt('skybgfiles/cp_skybg_zm_100_10_ph.dat')
        #dlam = 0.02 # d lambda = 0.02 nm for above data set
        #gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints/ir-background-spectra
    # The following is derived from the CP sky. Adds 10C Blackbody to spectrum
    # to try and reproduce Mt. Ham K values, J&H are unchanged
    elif skymod == 'mh':
        sky  = np.loadtxt('skybgfiles/lick_sky_zenith.txt')
        #dlam = 0.02 # d lambda = 0.02 nm for above data set
    elif skymod == 'mk':
        # ShaneAO spreadsheet uses MK data x 2 absent other data - worth a try
        # sky  = np.loadtxt('skybgfiles/nearIR_skybg_16_15_r5.dat')
        sky  = np.loadtxt('skybgfiles/mk_skybg_zm_16_15_ph.dat')
        # dlam = 0.1
        #dlam = 0.02

    skylambdas = sky[:,0]   # sampling of 0.02nm and a resolution of 0.4nm
                            # or 0.04 nm - Gemini website above has typo
    inir = np.where(skylambdas <= 2600)
    lambdas = skylambdas[inir]
    #dlam   = round(lambdas[1]-lambdas[0],2)
    
    # Sky transmission data
    if skymod == 'mk':
        # Multiply mauna kea flux by 2 and divide by airmass to get airmass=1.0
        # estimate
        skyflux = 2*sky[:,1][inir]/1.5 # ph s-1 m-2 nm-1 arcsec-2
        trans = np.loadtxt('skybgfiles/mktrans_zm_16_10.dat')
    else:
        skyflux = sky[:,1][inir] 
        trans = np.loadtxt('skybgfiles/cptrans_zm_100_10.dat')

    totel,tstel,tetel,thtel = emissivity(lambdas, shaneao, skyflux, shaneaoT,
                                         shaneaodel, dustfrac=df)

    # column 0 of trans is wavelength in microns but it corresponds to
    # lambdas (nm) above - so it can be ignored
    skytrans = trans[:,1][inir]

    # convert lambdas to cm and then create flux table of mag 0 object
    # in units of ph/s/nm/m^2 by taking zeroptflux/(hc/lambda)
    # this is the flux at the top of the atmosphere
    # mag0flux = zeroptflux*lambdas*100.0/(1.0e9*hc)

    # linear interpolate to sky lambdas and conver to photons s-1 m-2 nm-1
    srclint = np.interp(lambdas, srcflux[:,0], srcflux[:,1])*lambdas*100.0/(1.0e9*hc) 
    # mag0 flux at the telescope after propagating through the sky
    # mag0attel = mag0flux * skytrans
    srcattel =  srclint * skytrans

    # Source at detector is source flux at primary multiplied by system
    # transmission vector
    # mag0atout = mag0attel * thtel
    srcatout = srcattel * thtel

    if resample:
        fluxvector = np.zeros([len(lambdas), 2], np.float)
        fluxvector[:,0] = lambdas
        fluxvector[:,1] = srcatout
        reflux = resamplesky(fluxvector, R=R)

    # telescope characteristics
    # collecting area = primary area reduced by secondary area
    collarea = np.pi * ((shane['prdia']**2-shane['secdia']**2)/4.0)

    # slit width in arcseconds at detector
    slitangle = (206265*shane['slitwd']*(shane['fratioc']/shane['fratioa'])/
                 shane['flength'])
    # slit angular radius in arcsecs
    slitarad  = slitangle/2.0
    # slit radius in pixels
    slitprad  = slitarad/ircam['pscale']

    output = np.zeros((len(speclam),), dtype=[('trans', np.float),
                                              ('Rsky', np.float),
                                              ('Rsrc', np.float),
                                              ('RN2n', np.float),
                                              ('Rstar', np.float),
                                              ('npix', np.float),
                                              ('limmag',np.float)])


    for i, lamb in enumerate(speclam):
        # dlambda for given wavelength -> R = lambda/sdlam
        sdlam = lamb/R

        # TBD - numerical integration of the 2D Airy function
        #       to figure out what the slit losses are
        # Airy core radius in arcsecs and pixels
        airycorearf = (3600*(180.0/np.pi)*1.21966*lamb/(shane['prdia']*1.0e9))
        airycorearp = airycorearf / ircam['pscale']
        # check if airy disk radius is greater or less than slit angular radius
        # print airycorearf, slitarad
        if airycorearf > slitarad:
            # if greater i.e. not-diffraction limited, chop off side segments
            # angle subtended by segment in radians
            thetaa = 2 * np.arccos(slitarad/airycorearf)
            thetap = 2 * np.arccos(slitprad/airycorearp)
            # segment area
            segaa  = airycorearf**2 * (thetaa - np.sin(thetaa))/2.0
            segap  = airycorearp**2 * (thetap - np.sin(thetap))/2.0
            # area in airy disk within slit in arcsec^2 and whole pixels
            aresel = np.pi*airycorearf**2 - 2*segaa
            nresel = np.round(np.pi*airycorearp**2 - 2*segap)
            # Fraction of airy disk within slit
            areafrac = aresel/(np.pi*airycorearf**2)
        else:
            # otherwise just calculate area of airy disk in arcsec^2 and pixels
            aresel = np.pi * airycorearf**2
            nresel = np.round(np.pi * airycorearp**2)
            areafrac = 1.0

        output['npix'][i] = nresel
        # figure out which filter the wavelength is in and set strehl
        # accordingly
        # If it's not in the filters, the sky clobbers any transmission
        sratio = 1.0
        for filt in filters.keys():
            if (lamb > filters[filt]['start']) & (lamb < filters[filt]['end']):
                # sdlam = filters[filt]['end'] - filters[filt]['start']
                sratio = strehl[filt]
        # print sratio

        # fraction of flux in airy disk        
        corefrac = 0.86 * sratio

        # indices in the lambdas array in the interval sdlambda around
        # wavelength of interest
        indices = (np.where((lambdas > (lamb - sdlam/2.0)) &
                            (lambdas < (lamb + sdlam/2.0))))

        # input flux from sky
        siflux = dlam*skyflux[indices].sum()

        # mean transmission through the sdlam window
        output['trans'][i] = np.mean(thtel[indices])
        trans    = thtel[indices]/len(indices[0])

        # photons s-1 m-2 <arcsec-2>
        # Flux from sky only at detector
        soflux = dlam*tstel[indices].sum()
        # Flux from emissivity only at detector
        eflux  = dlam*tetel[indices].sum()
        # Source + emissivity at detector
        oflux  = dlam*totel[indices].sum()

        # units photons s-1        
        # Flux from source at detector
        sflux = dlam*srcatout[indices].sum()

        # Flux in resel from sky & emissivity (R_sky) - photons s-1
        output['Rsky'][i] = collarea * aresel * oflux
        # Flux in airy disk - scaled by area within slit - photons s-1
        output['Rsrc'][i] = areafrac * corefrac * collarea * sflux

        # Read noise term - RN^2 * npix in airy core - fowler sampling for now
        # In photons or DN [n pixels * (e-/pixel)^2]
        output['RN2n'][i] = nresel * ircam['rnfowl']**2 

        # For t=300s and SNR of 5, what limiting magnitude w.r.t. vega
        # do we reach
        expt = 300.0 # s
        snr = 5.0
        A = expt**2
        B = -expt * snr**2
        C = -snr**2 * (output['Rsky'][i] * expt + output['RN2n'][i])
        Rstar = (-B + np.sqrt(B**2 - 4*A*C))/(2*A)
        output['Rstar'][i] = Rstar
        if output['Rsrc'][i] > 0:
            output['limmag'][i] = -2.5*np.log10(Rstar/(output['Rsrc'][i]))
        else:
            output['limmag'][i] = 0.0

    # Round printout to 3 decimal places
    np.set_printoptions(precision=3)
    # Suppress exponential notation for small numbers
    np.set_printoptions(suppress=True)

#    print 'thruput  : ', output['trans']
#    print 'Rsrc     : ', output['Rsrc']
#    print 'Rsky     : ', output['Rsky'] 
#    print 'RN2n     : ', output['RN2n']
#    print 'Fowler limiting mag:', output['limmag']

    return output


