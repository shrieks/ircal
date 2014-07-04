import numpy as np
# change to different background file with Keck filters or move Keck
# filters in here. 
from background import *
# from telesim import *

def simkeck(am=1.0, df=0.0, adel=0.0265, odel=0.0, strehl=1.0, fscale=1.0):

    print "Running Keck NGAO simulation model with parameters:"
    print 'airmass  = {0:.1f}, dust fraction = {1:.3f}'.format(am,df)
    print 'Al delta = {0:.4f}, Other delta   = {1:.4f}'.format(adel,odel)

    # ST mag system m = 0 = -2.5log10 F_lam - 21.10
    zeroptflux = 3.631e-4 # ergs s-1 nm-1 m-2
    hc = 1.98645e-16      # in cgs - erg cm

    filters =  { 'J': Jparams, 'H': Hparams, 'K': Kparams}

    # Surface list for keck AO
    keckao = ['Aluminum', 'Aluminum', 'Aluminum', 
              'Aluminum', 'X1_Silver_Extrap', 'Aluminum', 'X1_Silver_Extrap', 
              'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 
              'TT_Dichroic_Refl', 'CAF2_BULK', 'NIR_AR', 'NIR_AR', 'CAF2_BULK']
    # temperature for keck AO surfaces
    keckaoT = np.array([2.6, 2.6, 2.6, 5,5,5,5,5,5,5,5,5,5,5,5])
    # delta by which surface transmissivity is modified
    keckaodel = np.array([adel, adel, adel, adel, 
                          odel, adel, odel, 
                          odel, odel, odel, 
                          odel, odel,odel,odel,odel])
    # build telescope data structure to send to telescope sim
    # to test SNR - change dia to 0.25 m to get 1 arcsec diffraction limit in J
    keck = { 'dia': 10, # m
             'surfaces': keckao,
             'temps': keckaoT,
             'deltas': keckaodel
             }
    
    # read in mauna kea nearir sky background in 900-5600 nm range
    if am==1.5:
        # 1.6mm H2O vapor column, 1.5 airmass
        mksky = np.loadtxt('skybgfiles/nearIR_skybg_16_15_r5.dat')
        # dlam = 0.1
    else:
        # 1.6 mm H2O vapor, 1.0 airmass
        mksky = np.loadtxt('skybgfiles/mk_skybg_zm_16_10_ph.dat')
        # dlam = 0.02 # d lambda = 0.02 nm for the above data set

    lambdas = mksky[:,0] # sampling of 0.02nm and a resolution of 0.04nm
    dlam   = round(lambdas[1]-lambdas[0],2)

    skyflux = mksky[:,1]/am # ph/s/arcsec^2/nm/m^2, scaled to get flux
                            # at airmass=1.0

    totel,tstel,tetel,thtel = emissivity(lambdas, keckao, skyflux, keckaoT,
                                         keckaodel, dustfrac=df)

    if am == 1.0:
        # Build sources - readin sky transmission first
        mktrans = np.loadtxt('skybgfiles/mktrans_zm_16_10.dat')
        # column 0 of mktrans is wavelength in microns but it corresponds to
        # lambdas (nm) above - so it can be ignored
        skytrans = mktrans[:,1]

        # convert lambdas to cm and then create flux table of mag 0 object
        # in units of ph s-1 nm-1 m-2 by taking zeroptflux/(hc/lambda)
        # erg s-1 nm-1 m-2 * cm / (erg cm) = # s-1 nm-1 m-2
        # this is the flux at the top of the atmosphere
        mag0flux = zeroptflux*lambdas*100.0/(1.0e9*hc*fscale)
        vega = np.loadtxt('spectra/alpha_lyr_stis_005.txt')
        # convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
        vegaflux = vega[:,1]*1.0e5/fscale
        vegalamb = vega[:,0]/10.0   # table has wavelength in angstroms
        # linear interpolate to sky lambdas and conver to photons s-1 m-2 nm-1
        vegalint = np.interp(lambdas, vegalamb, vegaflux)*lambdas*100.0/(1.0e9*hc) 
        # mag0/vega flux at the telescope primary after propagating through the sky
        srcattel = mag0flux * skytrans
        vegattel = vegalint * skytrans

        # Source at detector is source flux at primary multiplied by system
        # transmission vector
        srcatout = srcattel * thtel
        vegatout = vegattel * thtel

    # primary area (actually should reduce by secondary area)
    prarea = np.pi * (keck['dia']/2.0)**2
    # angular radius of airy disk core for each wavelength, make sure lambda
    # and dia have same units
    airycorearl = (3600*(180.0/np.pi)*1.21966*(lambdas*1.0e-9)/keck['dia'])
    # airycorear = 1.5 # arcsecs
    # Two ways of dealing with area - angular radius from filter center lambda
    # or at each radius. Using center lambda for now
    corefrac  = 0.86 * strehl

    output = np.zeros((3,), dtype=[('filter', np.str_,2),
                                   ('siflux', np.float),
                                   ('imag', np.float),
                                   ('trans', np.float),
                                   ('soflux', np.float),
                                   ('oflux', np.float),
                                   ('eflux', np.float),
                                   ('omag', np.float),
                                   ('aflux', np.float),
                                   ('src', np.float),
                                   ('vega', np.float)])


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
        # Source + emissivity at detector
        output['oflux'][i] = dlam*totel[indices].sum()

        # Divide flux by average trans to get photons incident
        # at top of atmosphere, then divide by zero point to get mag
        output['omag'][i]  = (-2.5 * np.log10(output['oflux'][i]/
                                              (trans.sum()*
                                               filters[filt]['zpc'])))
        if am == 1.0:
            # Flux in airy core from sky and emissivity
            airycorearf = (3600*(180.0/np.pi)*1.21966*filters[filt]['cwvl']/
                           (keck['dia']))
            # angular area of airy core
            airycoreaa = np.pi*airycorearf**2

            output['aflux'][i] = corefrac * airycoreaa * output['oflux'][i]
            # units photons m-2 s-1        

            output['src'][i] = corefrac * dlam*srcatout[indices].sum()
            #        output['snr'][i] = (np.sqrt(prarea)*output['src'][i]/
            #                            np.sqrt(output['src'][i]+output['aflux'][i]))
            output['vega'][i]= dlam*vegatout[indices].sum()


    # Round printout to 3 decimal places
    np.set_printoptions(precision=3)
    # Suppress exponential notation for small numbers
    np.set_printoptions(suppress=True)

    print '              {}'.format(output['filter'])
    print 'sky iflux: ', output['siflux']
    print 'sky imag : ', output['imag'] 
    print 'sky trans: ', output['trans']
    print 'sky sflux: ', output['soflux'] 
    print 'sky omag : ', output['omag'] 
