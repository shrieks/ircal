import numpy as np

adel = 0.0265

odel = 0.0

oldshaneao = ['Aluminum', 'Aluminum', 'Aluminum', 'oldao_1turn',
              'oldao_tiptilt', 'oldao_oap1', 'oldao_dm', 'oldao_oap2',
              'NaR_Splitter', 'X1_Silver_Extrap', 'X1_Silver_Extrap',
              'X1_Silver_Extrap']

oldaotlist = [10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.]

oldaodlist = [adel, adel, odel, odel, odel, odel, odel, odel, odel,
              odel, odel, odel]

sky = np.loadtxt('skybgfiles/lick_sky_zenith.txt')

skyflux = sky[:,1]

lambdas = sky[:,0] 

dlam   = round(lambdas[1]-lambdas[0],2)

trans = np.loadtxt('skybgfiles/cptrans_zm_100_10.dat')

skytrans = trans[:,1]

from background import *
taotot,taosrc,taoem,taoth = emissivity(lambdas, oldshaneao, skyflux,oldaotlist, oldaodlist,dustfrac=0.0)

cstot, cssrc, csem, csth = emissivity(lambdas, oldshaneao, skyflux*0.0,oldaotlist, oldaodlist,dustfrac=0.0, allemiss=1.0)

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

javgtrans = np.mean(taoth[jidx])
havgtrans = np.mean(taoth[hidx])
kavgtrans = np.mean(taoth[kidx])
ksavgtrans = np.mean(taoth[ksidx])
kpavgtrans = np.mean(taoth[kpidx])

# See Connie email/note for expalanation of jtran and why needed
jtran = taoth[jidx]/len(jidx[0])
jflux = dlam*taotot[jidx].sum()/jtran.sum() # photons s-1 m-2 arcsec-2
jmag  = -2.5 * np.log10(jflux/Jzpc)         # mag arcsec-2

htran = taoth[hidx]/len(hidx[0])
hflux = dlam*taotot[hidx].sum()/htran.sum()
hmag  = -2.5 * np.log10(hflux/Hzpc) 

ktran = taoth[kidx]/len(kidx[0])
kflux = dlam*taotot[kidx].sum()/ktran.sum()
kmag  = -2.5 * np.log10(kflux/Kzpc) 

kstran = taoth[ksidx]/len(ksidx[0])
ksflux = dlam*taotot[ksidx].sum()/kstran.sum()
ksmag  = -2.5 * np.log10(ksflux/Kszpc) 

kptran = taoth[kpidx]/len(kpidx[0])
kpflux = dlam*taotot[kpidx].sum()/kptran.sum()
kpmag  = -2.5 * np.log10(kpflux/Kpzpc) 

print '              J       H       K      Ks       Kp'
print 'sky  mag : {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f}'.format(jmag, hmag, kmag, ksmag, kpmag)
print 'sky  flux : {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f}'.format(jflux, hflux, kflux, ksflux, kpflux)

dewarin = taotot + 0.181*cstot

jflux = dlam*dewarin[jidx].sum()/jtran.sum() # photons s-1 m-2 arcsec-2
jmag  = -2.5 * np.log10(jflux/Jzpc)         # mag arcsec-2

hflux = dlam*dewarin[hidx].sum()/htran.sum()
hmag  = -2.5 * np.log10(hflux/Hzpc) 

kflux = dlam*dewarin[kidx].sum()/ktran.sum()
kmag  = -2.5 * np.log10(kflux/Kzpc) 

ksflux = dlam*dewarin[ksidx].sum()/kstran.sum()
ksmag  = -2.5 * np.log10(ksflux/Kszpc) 

kpflux = dlam*dewarin[kpidx].sum()/kptran.sum()
kpmag  = -2.5 * np.log10(kpflux/Kpzpc) 

print '              J       H       K      Ks       Kp'
print 'sky  mag : {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f}'.format(jmag, hmag, kmag, ksmag, kpmag)
print 'sky  flux : {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f}'.format(jflux, hflux, kflux, ksflux, kpflux)
