import numpy as np
from background import *

def simshane(srcflux, skymod='cp', df=0.03, adel=0.0265, odel=0.01,
             subaps=16, source='pts', filts='irc', aoT=18.0, aperrad=40.0):
    """
    Simulate Shane AO 3m AO system + instrument emissivity and througput
    Input: 
          srcflux- 2-d array with wavelengths in nm in column 0
                   and flux in units of ergs s-1 m-2 nm-1 in col 1 
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
    print 'Al delta = {0:.4f}, Other delta   = {1:.4f}, aperture radius = {2}'.format(adel,odel,aperrad)


    # ST mag system m = 0 = -2.5log10 F_lam - 21.10
    # zeroptflux = 3.631e-4 # ergs s-1 nm-1 m-2
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
    ircam = { 'invgain': 2.25,   # electrons per ADU or DN 
              'pscale' : 0.035,  # arcsec pixel-1
              'pixsize': 18.0,   # micron pixel pitch
              'fov'    : 20.0,   # arcsec square
              'fovdia' : 28.0,   # arcsecs - circular aperture equiv on square
              'pixels' : 2048.0, # pixels squared
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
    # cass secondary mirror is 33" dia, but light baffle is 39" dia
    shane =  {'prdia': 3.048, # m
              'secdia': 0.99, # m - mean dia on Mt Ham site
              'surfaces': shaneao,
              'temps': shaneaoT,
              'deltas': shaneaodel
              }

    # Sky background data
    if skymod == 'cp':
        # Cerro Pachon sky at 10 mm H2O column and airmass 1.0, 9000' up
        sky  = np.loadtxt('skybgfiles/cp_skybg_zm_100_10_ph.dat')
        # wavelength in nm, flux in ph s-1 m-2 nm-1 arcsec-2
    # The following is derived from the CP sky. Adds 10C Blackbody to spectrum
    # to try and reproduce Mt. Ham K values, J&H are unchanged
    elif skymod == 'mh':
        sky  = np.loadtxt('skybgfiles/lick_sky_zenith.txt')
    elif skymod == 'mk':
        # ShaneAO spreadsheet uses MK data x 2 absent other data - worth a try
        # sky  = np.loadtxt('skybgfiles/nearIR_skybg_16_15_r5.dat')
        sky  = np.loadtxt('skybgfiles/mk_skybg_zm_16_15_ph.dat')


    skylambdas = sky[:,0]  
    inir = np.where(skylambdas < 2500)
    lambdas = skylambdas[inir]
    dlam   = round(lambdas[1]-lambdas[0],2)
    
    # Sky transmission data
    if skymod == 'mk':
        # Multiply mauna kea flux by 2 and divide by airmass to get airmass=1.0
        # estimate
        skyflux = 2*sky[:,1][inir]/1.5 # ph s-1 m-2 nm-1 arcsec-2
        trans = np.loadtxt('skybgfiles/mktrans_zm_16_10.dat')
    else:
        skyflux = sky[:,1][inir] 
        trans = np.loadtxt('skybgfiles/cptrans_zm_100_10.dat')
    
    #total output (sky+emiss), sky output, emissivity and telescope transmission
    totel,tstel,tetel,thtel = emissivity(lambdas, shaneao, skyflux, shaneaoT,
                                         shaneaodel, dustfrac=df)

    # column 0 of trans is wavelength in microns but it corresponds to
    # lambdas (nm) above - so it can be ignored
    skytrans = trans[:,1][inir]

    # convert lambdas to cm and then create flux table of mag 0 object
    # in units of ph/s/nm/m^2 by taking zeroptflux/(hc/lambda)
    # this is the flux at the top of the atmosphere
    # mag0flux = zeroptflux*lambdas*100.0/(1.0e9*hc)

    # linear interpolate to sky lambdas and convert to photons s-1 m-2 nm-1
    # hc in ergs cm-1, lambdas converted to cm
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

    output = np.zeros((len(filters),), dtype=[('filter', np.str_,2),
                                              ('siflux', np.float),
                                              ('imag', np.float),
                                              ('trans', np.float),
                                              ('soflux', np.float),
                                              ('eflux', np.float),
                                              ('oflux', np.float),
                                              ('omag', np.float),
                                              ('Rsky', np.float),
                                              ('Remi', np.float),
                                              #('mag0', np.float),
                                              ('sflux',np.float),
                                              ('Rsrc', np.float),
                                              ('RN2n', np.float),
                                              ('limmag',np.float)])

    for i, filt in enumerate(filters.keys()):
        # indices in the lambdas array corresponding to the filter
        indices = (np.where((lambdas > filters[filt]['start']) &
                            (lambdas < filters[filt]['end'])))
        # set filter name
        output['filter'][i] = filt
        # input flux magnitude (needed for sky magnitudes)
        output['siflux'][i] = dlam*skyflux[indices].sum()
        output['imag'][i]  = (-2.5*np.log10(output['siflux'][i]/
                                            filters[filt]['zpc']))
        # mean transmission through the filter
        output['trans'][i] = np.mean(thtel[indices])

        trans    = thtel[indices]/len(indices[0])
        # photons s-1 m-2 <arcsec-2>
        # Flux from sky only at detector
        output['soflux'][i]= dlam*tstel[indices].sum()
        # Flux from emissivity only at detector
        output['eflux'][i] = dlam*tetel[indices].sum()
        # Sky + emissivity at detector
        output['oflux'][i] = dlam*totel[indices].sum()
        # source at detector - photons s-1 m-2
        output['sflux'][i] = dlam*srcatout[indices].sum()

        # Divide flux by average trans to get photons incident
        # at top of atmosphere, then divide by zero point to get mag
        output['omag'][i]  = (-2.5 * np.log10(output['oflux'][i]/
                                              (trans.sum()*
                                               filters[filt]['zpc'])))
        
        # Airy core radius - now user defined based on photometry radius in real images
        #airycorearf = (3600*(180.0/np.pi)*1.21966*filters[filt]['cwvl']/(shane['prdia'])) # arcsecs
        airycorearp = aperrad #airycorearf / ircam['pscale'] # pixels
        airycorearf = airycorearp * ircam['pscale']
        # angular area of airy core
        airycoreaa = np.pi*airycorearf**2 # in arcsec^2
        airycoreap = np.ceil(np.pi*airycorearp**2) # in whole pixels

        corefrac =  1.0 # 0.86 * strehl[filt]
        #print 'core fraction = ', corefrac
        # Flux in airy disk from sky & emissivity (R_sky) - photons s-1
        if source == 'pts':
            output['Rsky'][i] = collarea*airycoreaa*output['oflux'][i]
            output['Remi'][i] = collarea*airycoreaa*output['eflux'][i]
            # units photons s-1        
            output['Rsrc'][i] = corefrac*collarea*output['sflux'][i]

        # Read noise term - RN^2 * npix in airy core - fowler sampling for now
        # In photons or DN [n pixels * (e-/pixel)^2]
            output['RN2n'][i] = airycoreap * ircam['rnfowl']**2 
        elif source == 'csb':
            npix = np.round(1.0/(ircam['pscale']**2))
            # rates in photons or e- / s-1
            # area assumed to be 1 arcsec^2 for constant surface brightness obj
            output['Rsky'][i] = collarea*output['oflux'][i]
            output['Remi'][i] = collarea*output['eflux'][i]
            output['Rsrc'][i] = strehl[filt]*collarea*output['sflux'][i]
            output['RN2n'][i] = npix * ircam['rnfowl']**2 
        # elif source == 'ext':
            

        # For t=300s and SNR of 5, what limiting magnitude w.r.t. vega
        # do we reach
        expt = 300.0 # s
        snr = 5.0
        A = expt**2
        B = -expt * snr**2
        C = -snr**2 * (output['Rsky'][i] * expt + output['RN2n'][i])
        Rstar = (-B + np.sqrt(B**2 - 4*A*C))/(2*A)
        output['limmag'][i] = -2.5*np.log10(Rstar/(output['Rsrc'][i]))


    # Round printout to 3 decimal places
    np.set_printoptions(precision=3)
    # Suppress exponential notation for small numbers
    np.set_printoptions(suppress=True)

    print '              {}'.format(output['filter'])
#    print 'sky iflux: ', output['siflux']
    print 'sky  imag: ', output['imag'] 
    print 'thruput  : ', output['trans']
#    print 'sky oflux: ', output['soflux']
#    print 'sky eflux: ', output['eflux']
    print 'totalflux: ', output['oflux'] 
    print 'sky  omag: ', output['omag'] 
    print 'src  flux: ', output['sflux']
    print 'Rsrc     : ', output['Rsrc']
    print 'Rsky+emi : ', output['Rsky']
    print 'Remi     : ', output['Remi'] 
    print 'RN2n     : ', output['RN2n']
    print 'Tot flux : ', output['Rsrc'] + output['Rsky']
    print 'Fowler limiting mag:', output['limmag']

    return output


