import numpy as np
from background import *

# airmass, am = 1.0 or 1.5; df = 0.0 or 0.008; aldelta = 0.0265 or 0.0165
# otherdel = 0.0 or 0.005
# At some point, clean this code up to put filters in a structure with
# name, central lambda, bandwidth, zero point etc. 
# iterate rather than having sequential code - makes it easier to modify
def mkskybg(am=1.5, df=0.0, aldelta=0.0265, otherdel=0.0, amscale=True): 

    # Collecting area of keck telescope in m^2
    dia = 10.0 # m
    prarea = np.pi*(dia/4)**2

    # read in mauna kea nearir sky background in 900-5600 nm range
    if am==1.5:
        # 1.6mm H2O vapor column, 1.5 airmass
        mksky = np.loadtxt('skybgfiles/nearIR_skybg_16_15_r5.dat')
        dlam = 0.1
    else:
        # 1.6 mm H2O vapor, 1.0 airmass
        mksky = np.loadtxt('skybgfiles/mk_skybg_zm_16_10_ph.dat')
        dlam = 0.02 # d lambda = 0.02 nm for the above data set

    lambdas = mksky[:,0] # sampling of 0.02nm and a resolution of 0.04nm

    if amscale == False:
        am=1.0

    skyflux = mksky[:,1]/am # ph/sec/arcsec^2/nm/m^2, scaled to get flux
                            # at airmass=1.0

    # Filter central wavelengths from keckngao/ngao_bkg_filters.dat in nm
    Jcwvl = 1.250E-06/1e-9
    Hcwvl = 1.635E-06/1e-9
    Kpcwvl= 2.120E-06/1e-9 
    Kscwvl= 2.150E-06/1e-9
    Kcwvl = 2.200E-06/1e-9

    # Zero point fluxes for various filters (i.e. flux of mag 0 object)
    # from ngao_bkg_filters.dat
    Jzp  = 3.015e+09 # photons/(s m^2)? 
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


    # jidx is a tuple here whose only element is an array
    jidx  = np.where((lambdas > 1170) & (lambdas < 1330))
    hidx  = np.where((lambdas > 1490) & (lambdas < 1780))
    ksidx = np.where((lambdas > 1990) & (lambdas < 2310))
    kidx  = np.where((lambdas > 2030) & (lambdas < 2370))
    kpidx = np.where((lambdas > 1950) & (lambdas < 2290))

    jsky = dlam*skyflux[jidx].sum()    # photons s-1 m-2 arcsec-2
    jskymag = -2.5*np.log10(jsky/Jzpc) # mag arcsec-2

    hsky = dlam*skyflux[hidx].sum()
    hskymag = -2.5*np.log10(hsky/Hzpc)

    ksky = dlam*skyflux[kidx].sum()
    kskymag = -2.5*np.log10(ksky/Kzpc)

    kssky = dlam*skyflux[ksidx].sum()
    ksskymag = -2.5*np.log10(kssky/Kszpc)

    kpsky = dlam*skyflux[kpidx].sum()
    kpskymag = -2.5*np.log10(kpsky/Kpzpc)

    # Keck surfaces at various temperatures
    # Primary, secondary etc. at 2.6 C

# one silvered surface for the rotator; 
# 4 silvered surfaces for the AO system:
# 2 OAPs, Tip-Tilt mirror, DM; dichroic first surface, dichroic second surface,
# instrument window. Could add another surface for the bulk but that probably 
# over-counts it.

# three aluminum mirrors for the telescope
# aluminum-silver-aluminum for the rotator
# 4 silvered surfaces for the AO system: 2 OAPs, Tip-Tilt mirror, DM
# dichroci first surface, dichroic second surface, instrument window. 
# assume the second surface of the window is cold, so no emissivity,
# and nothing for the bulk because that is trivial. look later to check. 
    keck = ['Aluminum', 'Aluminum', 'Aluminum', 
            'Aluminum', 'X1_Silver_Extrap', 'Aluminum', 
            'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 
            'X1_Silver_Extrap', 'TT_Dichroic_Refl', 
            'X1_Silver_Extrap', 'X1_Silver_Extrap']

    keckTarr = np.array([2.6, 2.6, 2.6, 
                         5.0,5.0,5.0,
                         5.0,5.0,5.0,5.0,5.0,5.0,5.0])

    deltas =   np.array([aldelta, aldelta, aldelta, 
                         aldelta, otherdel, aldelta, 
                         otherdel,otherdel,otherdel,otherdel,
                         otherdel,otherdel,otherdel])

    keckao = ['Aluminum', 'Aluminum', 'Aluminum', 
              'Aluminum', 'X1_Silver_Extrap', 'Aluminum', 'X1_Silver_Extrap', 
              'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 
#              'TT_Dichroic_Refl', 'CaF', 'ARIR', 'ARIR', 'CaF']
              'TT_Dichroic_Refl', 'CAF2_BULK', 'NIR_AR', 'NIR_AR', 'CAF2_BULK']

    keckaoT = np.array([2.6, 2.6, 2.6, 5,5,5,5,5,5,5,5,5,5,5,5])

    keckaodel = np.array([aldelta, aldelta, aldelta, aldelta, 
                          otherdel, aldelta, otherdel, 
                          otherdel, otherdel, otherdel, 
                          otherdel, otherdel,otherdel,otherdel,otherdel])


    # Total emissions of different systems
    tekeck, thkeck = emissivity(lambdas, keck, skyflux, keckTarr, deltas, 
                                dustfrac=df)

    tekeckao, thkeckao = emissivity(lambdas, keckao, skyflux, keckaoT, 
                                    keckaodel, dustfrac=df)


    # Calculate average transmissivity in each filter
    javgtrans = np.mean(thkeckao[jidx])
    havgtrans = np.mean(thkeckao[hidx])
    kavgtrans = np.mean(thkeckao[kidx])
    # print 'AO trans: {0:.3f}'.format(thkeckao)

    # See Connie email/note for expalanation of jtran and why needed
    jtran = thkeck[jidx]/len(jidx[0])
    jflux = dlam*tekeck[jidx].sum()/jtran.sum() # photons s-1 m-2 arcsec-2
    jmag  = -2.5 * np.log10(jflux/Jzpc)         # mag arcsec-2

    aojtran = thkeckao[jidx]/len(jidx[0])
    aojflux = dlam*tekeckao[jidx].sum()/aojtran.sum()
    aojmag  = -2.5 * np.log10(aojflux/Jzpc) 
    # Angular area of airy disk core for filter, in arcsec^2
    # 3600 * (deg/radian) * 1.22 * cwvl/D
    # need to do fancier stuff later
    jangarea = 3600*(180.0/np.pi)*1.21966*Jcwvl/dia
    # Flux (photons/s) in Airy disk core at detector (100% Strehl)
    aojairy = aojflux * prarea * jangarea

    htran = thkeck[hidx]/len(hidx[0])
    hflux = dlam*tekeck[hidx].sum()/htran.sum()
    hmag  = -2.5 * np.log10(hflux/Hzpc) 

    aohtran = thkeckao[hidx]/len(hidx[0])
    aohflux = dlam*tekeckao[hidx].sum()/aohtran.sum()
    aohmag  = -2.5 * np.log10(aohflux/Hzpc) 
    hangarea = 3600*(180.0/np.pi)*1.21966*Hcwvl/dia
    aohairy = aohflux * prarea * hangarea

    ktran = thkeck[kidx]/len(kidx[0])
    kflux = dlam*tekeck[kidx].sum()/ktran.sum()
    kmag  = -2.5 * np.log10(kflux/Kzpc) 

    aoktran = thkeckao[kidx]/len(kidx[0])
    aokflux = dlam*tekeckao[kidx].sum()/aoktran.sum()
    aokmag  = -2.5 * np.log10(aokflux/Kzpc) 
    kangarea = 3600*(180.0/np.pi)*1.21966*Kcwvl/dia
    aokairy = aokflux * prarea * kangarea
    
    kstran = thkeck[ksidx]/len(ksidx[0])
    ksflux = dlam*tekeck[ksidx].sum()/kstran.sum()
    ksmag  = -2.5 * np.log10(ksflux/Kszpc) 

    kptran = thkeck[kpidx]/len(kpidx[0])
    kpflux = dlam*tekeck[kpidx].sum()/kptran.sum()
    kpmag  = -2.5 * np.log10(kpflux/Kpzpc) 

    #print jflux, hflux, kflux, ksflux, kpflux
    print '              J       H       K'
    print 'sky  mag : {0:.3f} {1:.3f} {2:.3f} '.format(jskymag, hskymag, kskymag)
    print 'upto dich: {0:.3f} {1:.3f} {2:.3f} '.format(jmag, hmag, kmag)
    print 'incl bulk: {0:.3f} {1:.3f} {2:.3f}\n'.format(aojmag, aohmag, aokmag)

    return kmag, aokmag


ams = 'On'
# def mkskybg(model='mk', df=0.0, aldelta=0.0265, otherdel=0.0):
print "MK, am=1.0, df=0.0, aldelta = 0.0265, otherdelta=0.0:   "
km5,aokm5 = mkskybg(am=1.0)
print "MK, am=1.0, df=0.008, aldelta = 0.0265, otherdelta=0.0:   "
km6,aokm6  = mkskybg(am=1.0,df=0.008)

print "Airmass scaling is ", ams
if ams == 'On':
    print "NG, am=1.5, df=0.0, aldelta = 0.0265, otherdelta=0.0:   "
    km1,aokm1  = mkskybg(am=1.5)
    print "NG, am=1.5, df=0.008, aldelta = 0.0265, otherdelta=0.0:   "
    km2,aokm2  = mkskybg(am=1.5,df=0.008)

    print "NG, am=1.5, df=0.0, aldelta = 0.0165, otherdelta=0.005:   "
    km3,aokm3  = mkskybg(am=1.5,aldelta=0.0165,otherdel=0.005)
    print "NG, am=1.5, df=0.008, aldelta = 0.0165, otherdelta=0.005:   "
    km4,aokm4  = mkskybg(am=1.5,df=0.008,aldelta=0.0165,otherdel=0.005)

    print "NG, am=1.5, df=0.0, aldelta = 0.0165, otherdelta=0.0: "
    km3,aokm3  = mkskybg(am=1.5,aldelta=0.0165)
    print "NG, am=1.5, df=0.008, aldelta = 0.0165, otherdelta=0.0:"
    km4,aokm4  = mkskybg(am=1.5,df=0.008,aldelta=0.0165)
else:
    km1,aokm1  = mkskybg(am=1.5, amscale=False)
    km2,aokm2  = mkskybg(am=1.5,df=0.008, amscale=False)

    km3,aokm3  = mkskybg(am=1.5,aldelta=0.0165,otherdel=0.005,amscale=False)
    km4,aokm4  = mkskybg(am=1.5,df=0.008,aldelta=0.0165,otherdel=0.005,
                         amscale=False)


