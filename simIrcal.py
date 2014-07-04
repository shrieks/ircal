import numpy as np
from background import *

def simIrcal(srcflux, df=0.03, adel=0.0265, odel=0.01, skymod='cp', aoT=18.0,
             sref=0.65):
    """
    Simulate Shane AO 3m AO system + IRCAL emissivity and througput
    Input: 
          srcflux- 2-d array with wavelengths in nm in column 0
                   and flux in units of photons s-1 m-2 nm-1 in col 1 
          df     - dust fraction (usually 0.0 or 0.008 (0.8%)
          adel   - aluminum reflectivity delta to reduce Al reflectivity
          odel   - delta by which to reduce other surface reflectivity
                   mimics degradation for both
          skymod - *cp* - unmodified Cerro Pachon sky, am 1.0, 10 mm H2O
                   *mh* - Connie's modified sky based on CP
          aoT    - AO system temperature in C
          sref   - reference Strehl
    """

    print "Running Shane Ircal simulation model with parameters:"
    print 'skymodel = {0}, dust fraction = {1:.3f}'.format(skymod,df)
    print 'Al delta = {0:.4f}, Other delta   = {1:.4f}'.format(adel,odel)
    print 'AO sys temp = {0:.1f}, K strehl = {1:.2f}'.format(aoT, sref)

    # ST mag system m = 0 = -2.5log10 F_lam - 21.10
    # zeroptflux = 3.631e-4 # ergs s-1 nm-1 m-2
    hc = 1.98645e-16      # in cgs - erg cm

    filters =  { 'J': J, 'H': H, 'K': K, 'Ks': Ks }

    # K band is Nyquist sampled by IRCAL
    # Changed for K filters 0.65 at 2 um, 0.65 with NGS - Olivier/Gavel paper
    # 8x8 subapertures
    # J & H is determined by formula: S = exp[-(lambda_0/lambda)^2 log(1/S_0)]
    # IRCAL published numbers: astro.berkeley.edu/~jrg/ircal/spie/ircal.html
    # JRGraham suggests Strehl is around 0.35 in K band.
    # collecting area is a degree of freedom (whether to subtract 
    # secondary or not) (minor), temp is the other - 18 C in Mt Ham summer
    # reference filter and strehl
    sfil = K
    # sref = 0.65

#    strehl8 =  { 'J': 0.25,    'H': 0.45,    'K': 0.42,     'Ks': 0.42, 
#                 'K\'': 0.42} 
    # Incorporates Don's guesses from spreadsheet for J, H, K
#    strehl8 =  { 'J': 0.1,    'H': 0.2,    'K': 0.3,     'Ks': 0.3, 
#                 'K\'': 0.3} 

    strehl8 =  { 'J': np.exp(-(sfil['cwvl']/J['cwvl'])**2 * np.log(1/sref)),    
                 'H': np.exp(-(sfil['cwvl']/H['cwvl'])**2 * np.log(1/sref)),
                 'K': np.exp(-(sfil['cwvl']/K['cwvl'])**2 * np.log(1/sref)),
                 'Ks': np.exp(-(sfil['cwvl']/Ks['cwvl'])**2 * np.log(1/sref)),
                 'K\'': np.exp(-(sfil['cwvl']/Kp['cwvl'])**2 * np.log(1/sref))} 
    # print strehl8

    # 16x16 subapertures - not used on IRCAL

    # PICNIC detector parameters
    ircal = { 'invgain': 10.0,   # electrons per ADU or DN
              'pscale' : 0.0756, # arcsec pixel-1
              'pixsize': 40.0,   # micron pixels
              'fov'    : 19.4,   # arcsec square
              'pixels' : 256.0,  # pixels squared
              'rdnoise': 30.0,   # electrons per CDS read
              'exptime': 57.0,   # one full read time in ms (3 us per pixel)
              'rnfowl' : 12.0,   # fowler read noise in e-/pixel
              'fowtime': 912.0   # 16 reads in ms ~1s here
              }


    # Primary, Secondary, Turn 1, Tip/tilt, OAP1, DM, OAP2, Na Dichroic,
    # IR mirror1, IR mirror2, Dewar window
#    oldshaneao = ['Aluminum', 'Aluminum', 'oldao_1turn', 'oldao_tiptilt', 
#                  'oldao_oap1', 'oldao_dm', 'oldao_oap2', 'NaR_Splitter', 
#                  'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap']

    oldshaneao = ['Aluminum', 'Aluminum', 'oldao_1turn', 'oldao_tiptilt', 
                  'oldao_oap1', 'oldao_dm', 'oldao_oap2', 'NaR_Splitter', 
                  'X1_Silver_Extrap', 'X1_Silver_Extrap', 'CaF']

    oldaotlist = [aoT, aoT, aoT, aoT, aoT, aoT, aoT, aoT, aoT, aoT, aoT]

    oldaodlist = [adel, adel, 
                  odel, odel, odel, odel, odel, odel, odel, odel, odel]

    # Window, Turn mirror, OAP1, Filter, OAP2, Turn mirror, cold stop, detector
#    ircal =   ['X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 
#               'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 
#               'ColdStop', 'PICNIC_QE']

    ircalslist =   ['CaF', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 'JHKpK', 
                    'X1_Silver_Extrap', 'X1_Silver_Extrap', 'ColdStop', 
                    'PICNIC_QE']

    ircaltlist = [-196.0, -196.0, -196.0, -196.0, -196.0, -196.0, -196.0, 
                   -196.0]

    ircaldlist = [odel, odel, odel, odel, odel, odel, odel, odel]


    # build telescope data structure to send to telescope sim
    shane =  {'prdia': 3.048, # m
              'secdia': 0.99, #m - 29.92" average from lick diag
              'surfaces': oldshaneao,
              'temps': oldaotlist,
              'deltas': oldaodlist
              }

    # Sky background
    if skymod == 'mh':
        sky = np.loadtxt('skybgfiles/lick_sky_zenith.txt')
    else:
        sky = np.loadtxt('skybgfiles/cp_skybg_zm_100_10_ph.dat')

    skylambdas = sky[:,0]   # sampling of 0.02nm and a resolution of 0.04nm
    inir = np.where(skylambdas < 2500)
    lambdas = skylambdas[inir]
    dlam   = round(lambdas[1]-lambdas[0],2)

    skyflux = sky[:,1][inir]

    # sky transmission. Wavelengths same as for sky background
    trans = np.loadtxt('skybgfiles/cptrans_zm_100_10.dat')
    skytrans = trans[:,1][inir]

    # ouput of emissivity program: total output photons from telescope,
    # photons from source, photons from emissivity, effective throughput
    # of system
    # telescope + AO
    taotot,taosrc,taoem,taoth = emissivity(lambdas, oldshaneao, skyflux, 
                                           oldaotlist, oldaodlist,
                                           dustfrac=df)
    # estimate effect of oversized cold stop - allows in emissivity from
    # rest of system outside dewar
    cstot, cssrc, csem, csth = emissivity(lambdas, oldshaneao, skyflux*0.0, 
                                          oldaotlist, oldaodlist,
                                          dustfrac=df, allemiss=1.0)
    # ircal input equal skyflux + 0.181*extra emissivity let in by cold
    # stop
    dewarin = taotot + 0.181*cstot
    totel,tstel,tetel,ircth = emissivity(lambdas, ircalslist, dewarin, 
                                         ircaltlist, ircaldlist,
                                         dustfrac=df)
    # telescope throughput is tel+ao thruput * dewar contents thruput
    thtel = taoth * ircth
        

    # convert lambdas to cm and then create flux table of mag 0 object
    # in units of ph/s/nm/m^2 by taking zeroptflux/(hc/lambda)
    # this is the flux at the top of the atmosphere
    # mag0flux = zeroptflux*lambdas*100.0/(1.0e9*hc)

    # linear interpolate to sky lambdas and convert to photons s-1 m-2 nm-1
    srclint = np.interp(lambdas, srcflux[:,0], srcflux[:,1])*lambdas*100.0/(1.0e9*hc) 
    # src flux at the telescope after propagating through the sky
    # mag0attel = mag0flux * skytrans
    srcattel =  srclint * skytrans

    # Source at detector is source flux at primary multiplied by system
    # transmission vector - photons s-1 m-2 nm-1
    # mag0atout = mag0attel * thtel
    srcatout = srcattel * thtel

    # collecting area (area of primary reduced by mean secondary area
    collarea = np.pi * ((shane['prdia']**2-shane['secdia']**2)/4.0)
    # collarea = np.pi * (shane['prdia']**2/4.0)

    output = np.zeros((len(filters),), dtype=[('filter', np.str_,2),
                                              ('siflux', np.float),
                                              ('imag', np.float),
                                              ('trans', np.float),
                                              ('soflux', np.float),
                                              ('oflux', np.float),
                                              ('eflux', np.float),
                                              ('omag', np.float),
                                              ('Rsky', np.float),
                                              # ('mag0', np.float),
                                              ('Rsrc', np.float),
                                              ('RN2n', np.float),
                                              ('limmagv', np.float)])

    for i, filt in enumerate(filters.keys()):
        # indices in the lambdas array corresponding to the filter
        indices = (np.where((lambdas > filters[filt]['start']) &
                            (lambdas < filters[filt]['end'])))

        # set filter name
        output['filter'][i] = filt
        # input flux magnitude (needed for sky magnitudes) in filter
        output['siflux'][i] = dlam*skyflux[indices].sum()
        output['imag'][i]  = (-2.5*np.log10(output['siflux'][i]/
                                            filters[filt]['zpc']))
        # mean transmission through the filter
        output['trans'][i] = np.mean(thtel[indices])

        trans    = thtel[indices]/len(indices[0])
        irctr    = np.mean(ircth[indices])
        # print filt, irctr
        # photons s-1 m-2 <arcsec-2>
        # Flux from sky only at detector  - photons s-1 m-2 arcsec-2
        output['soflux'][i]= dlam*tstel[indices].sum()
        # Flux from emissivity only at detector - photons s-1 m-2 arcsec-2
        output['eflux'][i] = dlam*tetel[indices].sum()
        # Source + emissivity at detector - photons s-1 m-2 arcsec-2
        output['oflux'][i] = dlam*totel[indices].sum()

        # Divide flux by average trans to get photons incident
        # at top of atmosphere, then divide by zero point to get sky
        # magnitude as measured by telescope
        output['omag'][i]  = (-2.5 * np.log10(output['oflux'][i]/
                                              (trans.sum()*
                                               filters[filt]['zpc'])))
        
        # Airy core radius in arcsecs using center wavelength of filter
        # then converted to pixels using instrument plate scale
        airycorearf = (3600*(180.0/np.pi)*1.21966*filters[filt]['cwvl']/
                       (shane['prdia']))    # arcsecs

        airycorearp = airycorearf / ircal['pscale'] # pixels
        # angular area of airy core
        airycoreaa = np.pi*airycorearf**2 # in arcsec^2
        airycoreap = np.ceil(np.pi*airycorearp**2) # in whole pixels

        corefrac = 0.86 * strehl8[filt]

        # Flux in airy core from sky & emissivity (R_sky) - photons s-1
        output['Rsky'][i] = collarea*airycoreaa*output['oflux'][i]

        # units photons s-1        
        output['Rsrc'][i] = corefrac*collarea*dlam*srcatout[indices].sum()

        # Read noise term - RN^2 * npix in airy core - for CDS read
        output['RN2n'][i] = airycoreap*(ircal['rnfowl'])**2

        # For t=300s and SNR of 5, what limiting magnitude w.r.t. vega
        # do we reach
        expt = 300.0 # s
        snr = 5.0
        A = expt**2
        B = -expt * snr**2
        C = -snr**2 * (output['Rsky'][i] * expt + output['RN2n'][i])
        Rstar = (-B + np.sqrt(B**2 - 4*A*C))/(2*A)
        output['limmagv'][i] = -2.5*np.log10(Rstar/(output['Rsrc'][i]))


    # Round printout to 3 decimal places
    np.set_printoptions(precision=3)
    # Suppress exponential notation for small numbers
    np.set_printoptions(suppress=True)

    print '              {}'.format(output['filter'])
#    print 'sky iflux: ', output['siflux']
    print 'sky imag : ', output['imag'] 
    print 'sky trans: ', output['trans']
#    print 'sky oflux: ', output['soflux']
#    print 'sky eflux: ', output['eflux']
#    print 'totalflux: ', output['oflux'] 
    print 'sky omag : ', output['omag'] 
    print ''
    print 'Rsky : ', output['Rsky'] 
    print 'Rsrc : ', output['Rsrc'] 
    print 'RN2n : ', output['RN2n'] 
    print 'lim magv:', output['limmagv']

    return output

