%load_ext autoreload

%autoreload 2

from emissivity import *

namelist = ['MEMS_window', 'PICNIC_QE', 'IR_transmission', 'Kprime', 
            'Hband', 'Jband', 'ALPAO', 'Aluminum', 'AluminumPlusDust', 
            'ARIR', 'CaF', 'FSG98', 'H2RG_QE', 'LGSWFS_Dich', 'NaHG', 
            'NaR_Splitter', 'TT_Dichroic_Refl','X1_Silver_Extrap']

lambdas = np.arange(600, 2350, 50)

t1 = 10

out, theff = emissivity(lambdas, namelist, t1)

output = np.column_stack((lambdas.flatten(),out.flatten(),theff.flatten()))

np.savetxt('emissivity-py.out',output)

clf()

grid(True)

yscale('log')

pyplot.xlim(800.0,2300.0)

pyplot.ylim(1e-15,1e5)

plot(lambdas, out, '-')

pyplot.xlabel('Wavelength (nm)')

pyplot.ylabel(r'#/(s m$^2$ arcsec$^2$ nm)')

#-------------------------------------------------------
# vegaflux
lamveg = vega[:,0]/10.0  # convert Angstroms to nm

clf()
plot(lamveg, vegaflux, '-')
vegain = vegaflux*lamveg*100/(1.0e9*hc)

clf()
plot(lamveg, vegain, '-')
plot(lambdas, vegalint, 'r-')
plot(lambdas, mag0flux, 'r-')
plot(lambdas, vegalint*skytrans, 'g-')
plot(lambdas, mag0flux*skytrans, 'b-')
xlim(900,2600)
ylim(0,5e8)

plot([2370,2370],[0,0.7e9],'y--')
plot([2030,2030],[0,1.6e9],'y--')
plot([1780,1780],[0,1.6e9],'y--')
plot([1490,1490],[0,1.6e9],'y--')
plot([1330,1330],[0,1.6e9],'y--')
plot([1170,1170],[0,1.6e9],'y--')
from matplotlib import rc
xlabel(r'Wavelength ($nm$)')
ylabel(r'$\gamma\ s^{-1} m^{-2}$')
text(1250, 4.5e8, 'J')
text(1635, 4.5e8, 'H')
text(2200, 4.5e8, 'K')
#-------------------------------------------------------------

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

#--------------------------------------

# fluxes in Ks band from sky+em and vega at detector
oflux = 10238.297152247636 # ph s-1 m-2 arcsec-2
sflux = 188741250.57290003 # ph s-1 m-2

corefrac = 0.86 * 0.6
collarea = np.pi*(3.048/2)**2 - np.pi*(0.76/2)**2
airycoreaa = np.pi*(206265*1.22*2.15e-6/3.048)**2
airycoreap = np.ceil(np.pi*(206265*1.22*2.15e-6/(3.048/0.0756))**2)

Rsky = collarea*airycoreaa*oflux
Rsrc = corefrac*collarea*sflux
Rn2n = 18*12**2

A = 300.0**2
B = -5**2 * 300
C = -5**2 * (Rsky * 300 + Rn2n)

Rstar = (-B + np.sqrt(B**2 - 4*A*C))/(2*A)
lmag = -2.5*np.log10(Rstar/Rsrc)
print lmag

#-----------------------------------------
# MOSFIRE fillters - creating a single envelope curve
mks = np.loadtxt('emissdatafiles/mosfire_Ks.txt')
mk = np.loadtxt('emissdatafiles/mosfire_K.txt')
mh = np.loadtxt('emissdatafiles/mosfire_H.txt')
mj = np.loadtxt('emissdatafiles/mosfire_J.txt')
# reverse array
mh = mh[::-1]
# concatenate
mjh = np.vstack((mj,mh))
# check for dupes in column 0 
[item for item, count in collections.Counter(mjh[:,0]).iteritems() if count > 1]
# there are none
# sort by column 0
a = mjh
mjhs = a[a[:,0].argsort()]
plot(mjhs[:,0], mjhs[:,1])
# cut mj above 1.4 and 1.4 < mh < 1.87
mj14 = mj[np.where(mj[:,0] < 1.4)]
mhc = mhc = mh[np.where((mh[:,0] > 1.4) & (mh[:,0] < 1.87))]
mkc = mk[np.where(mk[:,0] > 1.87)]

# k and ks
mk = mk[::-1]
mks = mks[::-1]
mkks = np.vstack((mk, mks))
a = mkks
mkks = a[a[:,0].argsort()]
[item for item, count in collections.Counter(mkks[:,0]).iteritems() if count > 1]
# lots of duplicates
for i, w in enumerate(mk[:,0]):
    mkks[i,0] = w
    if w in mks[:,0]:
        row = np.where(mks[:,0] == w)
        mkks[i,1] = maximum(mk[i,1], mks[row, 1])
    else:    
        mkks[i,1] = mk[i,1]

plot(mk[:,0], mk[:,1], 'r-.')
plot(mks[:,0], mks[:,1], 'g--')
plot(mkks[:,0], mkks[:,1], 'b-')

#-----------------------------------------------------------------
# Transmission curves

from matplotlib import rc
clf()
xlim(900,2500)
yscale('log')
xlabel('Wavelength (nm)')
ylabel(r'photons s$^{-1}$ m$^{-2}$ nm$^{-1}$ arcsec$^{-2}$')
#em, = plot(lambdas, taoem, 'r-')
#sk, = plot(lambdas, taosrc, 'g-')
#cs, = plot(lambdas, 0.181*cstot, 'm-')
#de, = plot(lambdas, dewarin, 'b-')

em, = plot(lambdas, stetel, 'r-')
sk, = plot(lambdas, ststel, 'g-')
to, = plot(lambdas, stotel, 'b-')
legend([sk, em, to],['Sky','Telescope+AO','Total'],loc='lower right')
#legend([sk, em, cs, de],['Sky','Telescope+AO','Cold Stop excess', 'Total'],loc='lower right')
title('Background and Emissivity - ShARCS')

clf()
xlim(900,2500)
ylim(0,1)
xlabel('Wavelength (nm)')
ylabel('Transmission')
sk, = plot(lambdas, skytrans, 'g-')
ao, = plot(lambdas, taoth, 'r-')
ir, = plot(lambdas, ircth, 'm-')
ai, = plot(lambdas, thtel, 'c-.')
to, = plot(lambdas, skytrans*thtel, 'b-')
title('System Transmission Curve - IRCAL')


de, = plot(lambdas, dewarin, 'r-')
to, = plot(lambdas, stotel, 'b-')
legend([de, to],['IRCAL','ShARCS'],loc='lower right')
plot([2310,2310],[0,0.7e9],'y--')
ylim(1e-5,1e4)
plot([2310,2310],[1e-5,1e4],'y--')
plot([1990,1990],[1e-5,1e4],'y--')
text(2150, 1e-3, 'Ks')

#--------------------------------------------------
# AGNs

agnprops = at.read('spectra/AGN/properties-edit.txt')
fileroot = 'spectra/AGN/'

for i, gal in enumerate(agnprops.Name):
    agns = np.loadtxt(fileroot+gal+'.txt')
    agnflux = agns[:,1]*1.0e5
    # table has wavelength in microns
    # shift lambda to rest frame
    agnlamb = agns[:,0]*1000.0/(1+agnprops.z[i])   
    
    plot(agnlamb, agnflux)
    grid(true)
    xlim(1100, 2400)
    xlabel(r'Rest Frame Wavelength ($nm$)')
    ylabel(r'Flux - ergs s$^{-1}$ m$^{-2}$ nm$^{-1}$')

#---------------------------------------------------
# Dwarfs

vlmb = at.read('tables/VLM_binaries.tsv')
ltyd = at.read('tables/LTYDwarfs.csv') #1280 all from 2MASS

ltyj = ltyd.jmag.astype(np.float)
ltyjlim = ltyj[np.where((ltyj < 22.28) & (ltyj > 0))]
ltyjlim.size # 1192 J, 1153-1157 H, 1020-1029 K - above lim mag in sample
             # ~80% of sample
             # visible in northern hemi and less than airmass 2 (60 deg)
             # means dec > -22.49 -> 983 J, 955 H, 840 K ~65%
# Worst case:
ltyjm = ltyd.jmag.astype(np.float) + ltyd.jmag_error.astype(np.float)
ltyjlim = ltyjm[np.where((ltyjm < 22.28) & (ltyjm > 0)& (ltyd.decl > -10.68))]
             # For airmass 1.5 892 J, 873 H, 752 K
             # ~60%

for vi, name in enumerate(vlmb.col1):
    if name in ltyd.designation:
        li = np.where(ltyd.designation == name)
        print name, ltyd.jmag[li], ltyd.hmag[li], ltyd.kmag[li]

#VLM Binaries & 2MASS intersection
#                            J         H         K
#2MASS J03202839-0446358 [ 13.259] [ 12.535] [ 12.134]
#2MASS J14044941-3159329 [ 15.577] [ 14.955] [ 14.538]
#2MASS J22551861-5713056 [ 14.083] [ 13.189] [ 12.579]
#------------------------------------------------------
import matplotlib.patches as mpa
pasb = 1281.5
clf()
plot(lambdas, np.median(specarr, axis=0), '-')
xlim(1200, 1350)
grid(True)
xlabel(r'Wavelength ($nm$)')
ylabel(r'Flux - ergs s$^{-1}$ m$^{-2}$ nm$^{-1}$')
title('Generic Type I AGN spectrum')

plot([fe2fb,fe2fb],[1e-11,8e-11],'g--')
plot([pasb,pasb],[1e-11,8e-11],'r--')

mp.text(1600, 1.5e-11, 'Continuum J band magnitude = 15.0')

mp.annotate('[FeII]', xy=(fe2fb, 7.5e-11),  xycoords='data',
            xytext=(-50, 0), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"))
mp.annotate(r'Pas-$\beta$', xy=(pasb, 7.5e-11),  xycoords='data',
            xytext=(25, 0), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"))

dlamf = fe2fb/R
dlamp = pasb/R

xyf = fe2fb-dlamf/2, 1e-11
widthf, htf = dlamf, 3e-11
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.7', edgecolor='0.7')
mp.gca().add_patch(pf)

xyp = pasb-dlamp/2, 1e-11
widthp, htp = dlamp, 4.5e-11
pp = mpa.Rectangle(xyp, widthp, htp, facecolor='0.7', edgecolor='0.7')
mp.gca().add_patch(pp)
mp.draw

#-------------------------------------------------------------
# AGN exp time and SNR plots

# Paschen beta plot
ymin = 1
ymax = 300

fig = mp.figure()
ax1 = fig.add_subplot(211)
ax1.grid(True)
ax2 = ax1.twiny()

ax1.set_xlim(xmin, xmax)
lo_ticks = np.array([1400., 1600., 1800., 2000., 2200.])
new_tick_locations = (lo_ticks-xmin)/(xmax-xmin)

ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(lo_ticks.astype(np.int))
ax2.set_xlabel(r"Wavelength $(nm)$")

#ax1.set_title(r'Paschen-$\beta$ (Top), [FeII] (Bot) SNR v. $z$')

def pb_tick_function(X):
    V = ((xmin + X * (xmax-xmin))/ pasb)-1.0
    return ["%.2f" % z for z in V]

lo_ticks = np.array([1400., 1600., 1800., 2000., 2200.])
new_tick_locations = (lo_ticks-xmin)/(xmax-xmin)
ax1.set_xticks(lo_ticks)
ax1.set_xticklabels(pb_tick_function(new_tick_locations))

#ax1.set_xticklabels([])
ax1.set_ylim(ymin, ymax)
ax1.set_yscale('log')
ax1.set_ylabel('SNR')
# reduce number of y-axis ticks
#ax1.set_yticks([1e-2,1,1e2])

# Plot Paschen-beta SNR vs. z for each mag
curves = []
for k,m in enumerate(mags):
    p, = ax1.plot(pasb*(1+zrange), snrarr[:,0,k]+1e-9, lines[k])
    curves.append(p)

ax1.text(xmin+10, 1.4, r'Paschen-$\beta$', size='large', weight='bold')
ax1.plot([xmin,xmax],[5,5],'k-',linewidth=2)
ax1.text(xmin+100, 5.5, 'SNR=5')
#ax1.plot([H['start'],H['start']],[ymin,ymax],'--',color='0.75',linewidth=1)
ax1.plot([H['start'],H['start']],[ymin,ymax],'b--',linewidth=1)
ax1.plot([H['end'],H['end']],[ymin,ymax],'b--',linewidth=1)
ax1.text(H['cwvl']*1e9, 1.4, 'H', color='b')
ax1.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--',linewidth=1)
ax1.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--',linewidth=1)
ax1.text(Ks['cwvl']*1e9, 1.4, 'Ks', color='r')

xy1 = xmin, ymin
width1, ht1 = H['start']-xmin, ymax
p1 = mpa.Rectangle(xy1, width1, ht1, facecolor='0.8', edgecolor='0.8', zorder=2)
ax1.add_patch(p1)

xy2 = H['end'], ymin
width2, ht2 = Ks['start']-H['end'], ymax
p2 = mpa.Rectangle(xy2, width2, ht2, facecolor='0.8', edgecolor='0.8', zorder=2)
ax1.add_patch(p2)


#  [FeII] plot
ax3 = fig.add_subplot(212)
ax3.grid(True)
#ax4 = ax3.twiny()

#ax3.set_title(r'[FeII] SNR')

ax3.set_xlim(xmin, xmax)
def fe_tick_function(X):
    V = ((xmin + X * (xmax-xmin))/ fe2fb)-1.0
    return ["%.2f" % z for z in V]

ax3.set_xticks(lo_ticks)
ax3.set_xticklabels(fe_tick_function(new_tick_locations))
ax3.set_xlabel(r'Redshift ($z$)')

ax3.set_ylim(ymin, ymax)
ax3.set_yscale('log')
ax3.set_ylabel('SNR')
#ax3.set_yticks([1e-2,1,1e2])

# Plot [FeII] SNR vs. z for each mag
curves = []
for k,m in enumerate(mags):
    p, = ax3.plot(fe2fb*(1+zrange), snrarr[:,1,k]+1e-9, lines[k])
    curves.append(p)

ax3.text(xmin+10, 1.4, r'[FeII]', size='large', weight='bold')
ax3.plot([xmin,xmax],[5,5],'k-',linewidth=2)
ax3.text(xmax-150, 5.5, 'SNR=5')
ax3.plot([H['start'],H['start']],[ymin,ymax],'b--',linewidth=1)
ax3.plot([H['end'],H['end']],[ymin,ymax],'b--',linewidth=1)
ax3.text(H['cwvl']*1e9, 1.4, 'H', color='b')
ax3.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--',linewidth=1)
ax3.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--',linewidth=1)
ax3.text(Ks['cwvl']*1e9, 1.4, 'Ks', color='r')

xyf = xmin, ymin
widthf, htf = H['start']-xmin, ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=2)
ax3.add_patch(pf)

xyf = H['end'], ymin
widthf, htf = Ks['start']-H['end'], ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=2)
ax3.add_patch(pf)

mp.legend(curves, mags, loc='upper left', title='J-band mags',
          prop={'size':12})

mp.show()

#-------------------
# Exposure time plot
fig = mp.figure()

ymin = 1
ymax = 2e3
# Paschen beta plot

ax1 = fig.add_subplot(211)
ax1.grid(True)
ax2 = ax1.twiny()
lo_ticks = np.array([1400., 1600., 1800., 2000., 2200.])
new_tick_locations = (lo_ticks-xmin)/(xmax-xmin)

ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(lo_ticks.astype(np.int))
ax2.set_xlabel(r"Wavelength $(nm)$")

def pb_tick_function(X):
    V = ((xmin + X * (xmax-xmin))/ pasb)-1.0
    return ["%.2f" % z for z in V]

ax1.set_xlim(xmin, xmax)
lo_ticks = np.array([1400., 1600., 1800., 2000., 2200.])
new_tick_locations = (lo_ticks-xmin)/(xmax-xmin)
ax1.set_xticks(lo_ticks)
ax1.set_xticklabels(pb_tick_function(new_tick_locations))

#ax1.set_xticklabels([])
ax1.set_ylim(ymin, ymax)
ax1.set_yscale('log')
ax1.set_ylabel(r'Exposure time ($s$)')
# reduce number of y-axis ticks
#ax1.set_yticks([1e-2,1,1e2])

# Plot Paschen-beta SNR vs. z for each mag
curves = []
for k,m in enumerate(mags):
    p, = ax1.plot(pasb*(1+zrange), exptarr[:,0,k]+1e-9, lines[k])
    curves.append(p)

ax1.text(xmin+10, 1.5, r'Paschen-$\beta$', size='large', weight='bold')
ax1.plot([xmin,xmax],[600, 600],'k-',linewidth=2)
ax1.text(xmin+10, 700, '10min Fowler-32')
#ax1.plot([H['start'],H['start']],[ymin,ymax],'--',color='0.75',linewidth=1)
ax1.plot([H['start'],H['start']],[ymin,ymax],'b--',linewidth=1)
ax1.plot([H['end'],H['end']],[ymin,ymax],'b--',linewidth=1)
ax1.text(H['cwvl']*1e9, 1.1, 'H', color='b')
ax1.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--',linewidth=1)
ax1.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--',linewidth=1)
ax1.text(Ks['cwvl']*1e9, 1.1, 'Ks', color='r')

xy1 = xmin, ymin
width1, ht1 = H['start']-xmin, ymax
p1 = mpa.Rectangle(xy1, width1, ht1, facecolor='0.8', edgecolor='0.8', zorder=2)
ax1.add_patch(p1)

xy2 = H['end'], ymin
width2, ht2 = Ks['start']-H['end'], ymax
p2 = mpa.Rectangle(xy2, width2, ht2, facecolor='0.8', edgecolor='0.8', zorder=2)
ax1.add_patch(p2)


#  [FeII] plot
ax3 = fig.add_subplot(212)
ax3.grid(True)
#ax4 = ax3.twiny()

def fe_tick_function(X):
    V = ((xmin + X * (xmax-xmin))/ fe2fb)-1.0
    return ["%.2f" % z for z in V]

ax3.set_xlim(xmin, xmax)
ax3.set_xticks(lo_ticks)
ax3.set_xticklabels(fe_tick_function(new_tick_locations))
ax3.set_xlabel(r'Redshift ($z$)')

ax3.set_ylim(ymin, ymax)
ax3.set_yscale('log')
ax3.set_ylabel(r'Exposure time ($s$)')
#ax3.set_yticks([1e-2,1,1e2])

# Plot [FeII] SNR vs. z for each mag
curves = []
for k,m in enumerate(mags):
    p, = ax3.plot(fe2fb*(1+zrange), exptarr[:,1,k]+1e-9, lines[k])
    curves.append(p)

ax3.text(xmin+10, 1.5, r'[FeII]', size='large', weight='bold')
ax3.plot([xmin,xmax],[600, 600],'k-',linewidth=2)
ax3.text(xmax-280, 700, '10min Fowler-32')

ax3.plot([H['start'],H['start']],[ymin,ymax],'b--',linewidth=1)
ax3.plot([H['end'],H['end']],[ymin,ymax],'b--',linewidth=1)
ax3.text(H['cwvl']*1e9, 1.1, 'H', color='b')
ax3.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--',linewidth=1)
ax3.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--',linewidth=1)
ax3.text(Ks['cwvl']*1e9, 1.1, 'Ks', color='r')

xyf = xmin, ymin
widthf, htf = H['start']-xmin, ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=2)
ax3.add_patch(pf)

xyf = H['end'], ymin
widthf, htf = Ks['start']-H['end'], ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=2)
ax3.add_patch(pf)

mp.legend(curves, mags, loc='upper left', title='J-band mags',
          prop={'size':12})

mp.show()

#-------------------------------------------------------------------
# Minimum exposure time for spec

import matplotlib.patches as mpa

limlam = np.arange(1100, 2314, 1.0)

vega = np.loadtxt('spectra/alpha_lyr_stis_005.txt')
# convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
vegaflux = vega[:,1]*1.0e5
vega[:,1] = vegaflux

vegalamb = vega[:,0]/10.0   # table has wavelength in angstroms
vega[:,0] = vegalamb

limit = simspec (vega, limlam, subaps=8)
rnf = 5.25 # e- pixel-1 for fowler-32 read
dc = 0.01  # e- s-1 pixel-1
bglimt = rnf/(dc + limit['Rsky']/limit['npix'])
srclimt = rnf/(limit['Rsky']/limit['npix'])

clf()
ymin = 1
ymax = 2e4 + 5e3
grid(True)
ylim(ymin, ymax)
yscale('log')
xlim(1100,2300)

pb, = plot(limlam, bglimt, 'k-')
ps, = plot(limlam, srclimt, 'r-')
legend([pb, ps],['Background limited', 'Source limited'], loc='upper right')

plot([J['end'],J['end']],[ymin,ymax],'b--')
plot([H['start'],H['start']],[ymin,ymax],'g--')
plot([H['end'],H['end']],[ymin,ymax],'g--')
plot([K['start'],K['start']],[ymin,ymax],'r--')

mp.text(J['cwvl']*1e9, 2, 'J', color='b')
mp.text(H['cwvl']*1e9, 2, 'H', color='g')
mp.text(Ks['cwvl']*1e9, 2, 'Ks', color='r')

#mp.plot([J['start'],K['end']],[300,300],'-', linewidth=2)
#mp.text(2100, 350, '300s Fowler')
#mp.text(2100, 200, 'Exposure')

xlabel(r'Wavelength ($nm$)')
ylabel('Time ($s$)')
title('Minimum Exposure Time for Spectroscopy-ShARCS 8') 

xyf = J['end'], ymin
widthf, htf = H['start']-J['end'], ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=3)
gca().add_patch(pf)

xyf = H['end'], ymin
widthf, htf = Ks['start']-H['end'], ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=3)
gca().add_patch(pf)

#-------------------------------------
# SKy spectrum with airmass/water column

#sky10 = np.loadtxt('skybgfiles/cp_skybg_zm_100_10_ph.dat')
sky10 = np.loadtxt('skybgfiles/cptrans_zm_100_10.dat')
#sky20 = np.loadtxt('skybgfiles/cp_skybg_zm_100_20_ph.dat')
sky20 = np.loadtxt('skybgfiles/cptrans_zm_100_20.dat')

clf()
p10, = plot(sky10[:,0]*1000, sky10[:,1], 'k-') 
#yscale('log')

xlim(900,2500)

p20, = plot(sky20[:,0]*1000, sky20[:,1], 'r-')

#ymin = 1e-5
#ymax = 1e5
ymin = 0
ymax = 1.4
xmin = 900
xmax = 2500
ylim(ymin,ymax)

plot([J['start'],J['start']],[ymin,ymax],'b--')
plot([J['end'],J['end']],[ymin,ymax],'b--')
plot([H['start'],H['start']],[ymin,ymax],'g--')
plot([H['end'],H['end']],[ymin,ymax],'g--')
plot([K['start'],K['start']],[ymin,ymax],'r--')
plot([Ks['end'],Ks['end']],[ymin,ymax],'r--')

#ftext = 2e4
ftext = 1.05
mp.text(J['cwvl']*1e9, ftext, 'J', color='b')
mp.text(H['cwvl']*1e9, ftext, 'H', color='g')
mp.text(Ks['cwvl']*1e9, ftext, 'Ks', color='r')

xlabel(r'Wavelength ($nm$)')
#ylabel(r'photons s$^{-1}$ m$^{-2}$ nm$^{-1}$ arcsec$^{-2}$')
ylabel('Transmission')
title('Cerro Pachon Sky Transmission')
mp.legend([p10, p20], ['10 mm H2O', '20 mm H2O'], loc='upper right', 
          title = 'Water column')

#------------------------------------------------
# Water column

import asciitable as at
import matplotlib.pyplot as mp

def movingaverage(interval, window_size):
    window= numpy.ones(int(window_size))/float(window_size)
    return numpy.convolve(interval, window, 'same')

lutz = at.read('skybgfiles/LUTZ-pwv.txt')

import matplotlib.dates as dates
import matplotlib.ticker as ticker

clf()
fig = figure()
ax = fig.add_subplot(111)
pd, = ax.plot(range(lutz.col1.size), lutz.col2, 'r-')
pw, = ax.plot(range(lutz.col1.size), movingaverage(lutz.col2,350), 'b-', 
              linewidth=2)
pm, = ax.plot(range(lutz.col1.size), movingaverage(lutz.col2,1300), 'k-',
              linewidth=2)

mticks = np.array([1300,2600,3900, 5200,6500])
months = np.array(['<Jan', '<Feb', '<Mar', '<Apr', '<May'])
hmonth = np.array([650, 1950, 3250, 4550, 5850])

ax.set_xticks(mticks)
ax.set_xticklabels(months, ha='right')
ax.set_xlim(0,6500)
ax.set_ylim(0, 35)
ax.set_ylabel('Precipitable Water Vapor (mm)')
ax.set_title('PWV for San Jose, 2013')

mp.legend([pd,pw,pm],['Daily', 'Weekly moving avg.', 'Monthly moving avg.'],
          loc = 'upper left')
ax.text(3900, 32, 'Mean, median = 11 mm')
