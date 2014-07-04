import numpy as np
import matplotlib.pyplot as mp
import matplotlib.patches as mpa
import simshspec as ss
from background import *
from itertools import cycle

# Read in vega flux file
vega = np.loadtxt('spectra/alpha_lyr_stis_005.txt')
# convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
vegaflux = vega[:,1]*1.0e5
vega[:,1] = vegaflux

vegalamb = vega[:,0]/10.0   # table has wavelength in angstroms
vega[:,0] = vegalamb

# K band strehl ratios
shastr = 0.8
shas = str(shastr)
olstr = 0.6
ols = str(olstr)

#Paschen-beta line
pasb = 1280 # nm

# Generate array of wavelengths for spectroscopy
R = 500.0
sn = 5.0
ft = 300.0 # fowler time in secs
filters = {'J': J, 'H': H, 'Ks': Ks}
# filters = ['J': MJ, 'H': MH, 'Ks': MKs] # mosfire

#for i, f in enumerate(filters.keys()):
    # space = filters[f]['cwvl']*1.0e9/(2*R) # in nm
space = 2.0 # nm
lamJ = np.arange(J['start'], J['end']+space, space)
lamH = np.arange(H['start'], H['end']+space, space)
lamKs= np.arange(Ks['start'], Ks['end']+space, space)
lj = len(lamJ)
lh = len(lamH)
lks= len(lamKs)
lams = np.hstack((lamJ, lamH, lamKs))
# mosfire
#mlamJ = np.arange(MJ['start'], MJ['end']+space, space)
#mlamH = np.arange(MH['start'], MH['end']+space, space)
#mlamKs= np.arange(MKs['start'], MKs['end']+space, space)
#mlams = np.hstack((mlamJ, mlamH, mlamKs))
#lmj = len(mlamJ)
#lmh = len(mlamH)
#lmks= len(mlamKs)

# Send through the telescope
vegaspec = ss.simspec(vega, lams)
#vegamos  = ss.simspec(vega, mlams, filts='mos')
vegamos  = ss.simspec(vega, lams, subaps=8)

# Plot of spectroscopy magnitude limits
ymin = 15
ymax = 19.5
mp.clf()
mp.grid(True)
mp.ylim(ymin, ymax)
mp.xlim(J['start'],Ks['end'])

pj, = mp.plot(lamJ, vegaspec['limmag'][0:lj], 'r-')
pmj,= mp.plot(lamJ,  vegamos['limmag'][0:lj], 'b-')
#pmj,= mp.plot(mlamJ, vegamos['limmag'][0:lmj], 'b-')

ph, = mp.plot(lamH, vegaspec['limmag'][lj:lj+lh], 'r-')
pmh,= mp.plot(lamH,  vegamos['limmag'][lj:lj+lh], 'b-')
#pmh,= mp.plot(mlamH, vegamos['limmag'][lmj:lmj+lmh], 'b-')

pks,= mp.plot(lamKs, vegaspec['limmag'][lj+lh:lj+lh+lks], 'r-')
pks,= mp.plot(lamKs,  vegamos['limmag'][lj+lh:lj+lh+lks], 'b-')
#pmks,=mp.plot(mlamKs, vegamos['limmag'][lmj+lmh:lmj+lmh+lmks], 'b-')
mp.xlabel('Wavelength (nm)')
mp.ylabel('Magnitudes')
mp.title('Spectroscopy magnitude limits, R=500')
#leg1 = mp.legend([pmj,pmh,pmks],['J','H','Ks'], loc='lower center', 
#                 title='ShARCS 8')
#mp.gca().add_artist(leg1)
#mp.legend([pj, ph, pks],['J', 'H','Ks'], loc='lower left', title='ShARCS 16')
mp.legend([pj, pmj],['ShARCS 16', 'ShARCS 8'], loc='lower left')

mp.plot([J['start'],J['start']],[ymin,ymax],'b--')
mp.plot([J['end'],J['end']],[ymin,ymax],'b--')
mp.plot([H['start'],H['start']],[ymin,ymax],'g--')
mp.plot([H['end'],H['end']],[ymin,ymax],'g--')
mp.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--')
mp.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--')

mp.text(J['cwvl']*1e9, ymax-0.5, 'J', color='b')
mp.text(H['cwvl']*1e9, ymax-0.5, 'H', color='g')
mp.text(Ks['cwvl']*1e9,ymax-0.5, 'Ks', color='r')

xyf = J['end'], ymin
widthf, htf = H['start']-J['end'], ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=3)
mp.gca().add_patch(pf)

xyf = H['end'], ymin
widthf, htf = Ks['start']-H['end'], ymax
pf = mpa.Rectangle(xyf, widthf, htf, facecolor='0.8', edgecolor='0.8', zorder=3)
mp.gca().add_patch(pf)

mp.savefig('figs/web/SpecMagLims.png')
#mp.show()
# End plot
#-------------------------------------

fi = open('tables/spec-ircFiltsR500.txt', 'w')
for i, lam in enumerate(lams):
    fi.write('{0:.1f} & {1:.2f} & {2:.3e} & {3:.3e} & {4:.3e}\\\\\n'.format(lam, vegaspec['limmag'][i], vegaspec['Rstar'][i]*ft, vegaspec['Rsky'][i]*ft, vegaspec['RN2n'][i]))
    fi.write('\hline\n')
    
fi.close()

#fm = open('tables/spec-mosFiltsR500.txt', 'w')
#for i, lam in enumerate(mlams):
#    fm.write('{0:.1f} & {1:.2f} & {2:.3e} & {3:.3e} & {4:.3e}\\\\\n'.format(lam#, vegamos['limmag'][i], vegamos['Rstar'][i]*ft, vegamos['Rsky'][i]*ft, vegamos[#'RN2n'][i]))
#    fm.write('\hline\n')
#    
#fm.close()


# magnitude range
mags = np.arange(15.0, 20.5, 1.0)

# exposure time array mags (rows) x filters (columns)
exptarr = np.zeros((len(lams), mags.size), np.float)
snrarr  = np.zeros((len(lams), mags.size), np.float)

ief = open('tables/spec-iFR500expt.txt', 'w')
isf = open('tables/spec-iFR500snr.txt', 'w')

for i, lam in enumerate(lams):
    slam = str(lam)

    # Sky photons s-1 in diffraction limited core for filter
    Rsky  = vegaspec['Rsky'][i]
    Rn2n  = vegaspec['RN2n'][i]

    for j, m in enumerate(mags):
        # m - 0 = -2.5 log (Rstar/Rvega)
        # Star photons s-1 in slit given magnitude
        Rstar = vegaspec['Rsrc'][i] * 10**(-m/2.5)

        # Solve the SNR quadratic equation for t
        exptarr[i,j] = ((sn**2*(Rstar + Rsky) + 
                         np.sqrt(sn**4*(Rstar + Rsky)**2 +
                                 4*sn**2*Rstar**2*Rn2n))
                        /(2*Rstar**2))

        ief.write('{0:.1f} & {1:.1f} & {2:.3e} & {3:.3e} & {4:.3e} & {5:.3e}\\\\\n'.format(lam, m, exptarr[i,j], Rstar*exptarr[i,j], Rsky*exptarr[i,j], Rn2n))
        ief.write('\hline\n')

        snrarr[i,j] = ((Rstar*ft)/
                       np.sqrt(Rstar*ft + Rsky*ft + Rn2n))
        isf.write('{0:.1f} & {1:.1f} & {2:.2f} & {3:.3e} & {4:.3e} & {5:.3e}\\\\\n'.format(lam, m, snrarr[i,j], Rstar*ft, Rsky*ft, Rn2n ))
        isf.write('\hline\n')

    ief.write('\hline\n')
    isf.write('\hline\n')

ief.close()
isf.close()

# MOSFIRE filters
mexptarr = np.zeros((len(mlams), mags.size), np.float)
msnrarr  = np.zeros((len(mlams), mags.size), np.float)

mef = open('tables/spec-mFR500expt.txt', 'w')
msf = open('tables/spec-mFR500snr.txt', 'w')

for i, lam in enumerate(mlams):
    slam = str(lam)

    # Sky photons s-1 in diffraction limited core for filter
    Rsky  = vegamos['Rsky'][i]
    Rn2n  = vegamos['RN2n'][i]

    for j, m in enumerate(mags):
        # m - 0 = -2.5 log (Rstar/Rvega)
        # Star photons s-1 in slit given magnitude
        Rstar = vegamos['Rsrc'][i] * 10**(-m/2.5)

        # Solve the SNR quadratic equation for t
        mexptarr[i,j] = ((sn**2*(Rstar + Rsky) + 
                          np.sqrt(sn**4*(Rstar + Rsky)**2 +
                                  4*sn**2*Rstar**2*Rn2n))
                         /(2*Rstar**2))

        mef.write('{0:.1f} & {1:.1f} & {2:.3e} & {3:.3e} & {4:.3e} & {5:.3e}\\\\\n'.format(lam, m, mexptarr[i,j], Rstar*mexptarr[i,j], Rsky*mexptarr[i,j], Rn2n))
        mef.write('\hline\n')

        msnrarr[i,j] = ((Rstar*ft)/
                        np.sqrt(Rstar*ft + Rsky*ft + Rn2n))
        msf.write('{0:.1f} & {1:.1f} & {2:.2f} & {3:.3e} & {4:.3e} & {5:.3e}\\\\\n'.format(lam, m, msnrarr[i,j], Rstar*ft, Rsky*ft, Rn2n ))
        msf.write('\hline\n')

    mef.write('\hline\n')
    msf.write('\hline\n')

mef.close()
msf.close()


lines = ["-","--","-.",":"]
linecycler = cycle(lines)

lam2 = [1152.5, 1252.5, 1352.5, 1508.0, 1608.0, 1708.0, 1808.0, 
        2000.0, 2100.0, 2200.0, 2300.0]
mlam2= [1153.0, 1253.0, 1353.0, 1506.5, 1606.5, 1706.5, 1806.5, 
        2000.0, 2100.0, 2200.0, 2300.0]

mp.clf()
mp.yscale('log')
mp.grid(True)
mp.title('Spectroscopy - Exposure time vs. Magnitude for SNR=5')
mp.xlabel('Magnitude (Vega system)')
mp.ylabel('Exposure time (s)')

curves = []
for lam in lam2:
    pts = np.where(lams == lam)
    p, = mp.plot(mags, exptarr[pts[0],:][0], next(linecycler))
    curves.append(p)

mp.legend(curves, lam2, loc='upper left', ncol=2, title='ShARCS 16')
mp.plot([min(mags),max(mags)],[300,300],'-', linewidth=2)
mp.text(np.min(mags)+0.1, 350, '300s Fowler')
mp.text(np.min(mags)+0.1, 175, 'Exposure')

mp.clf()
mp.yscale('log')
mp.grid(True)
mp.title('Spectroscopy - Exposure time vs. Magnitude for SNR=5')
mp.xlabel('Magnitude (Vega system)')
mp.ylabel('Exposure time (s)')
curves = []
for lam in mlam2:
    pts = np.where(mlams == lam)
    p, = mp.plot(mags, mexptarr[pts[0],:][0], next(linecycler))
    curves.append(p)
mp.legend(curves, mlam2, loc='upper left', ncol=2, title='ShARCS 8')
mp.plot([min(mags),max(mags)],[300,300],'-', linewidth=2)
mp.text(np.min(mags)+0.1, 350, '300s Fowler')
mp.text(np.min(mags)+0.1, 175, 'Exposure')


mp.clf()
mp.yscale('log')
mp.grid(True)
mp.title('Spectroscopy - SNR vs. Magnitude for 300s Fowler-32')
mp.xlabel('Magnitude (Vega system)')
mp.ylabel('SNR')
curves = []
for lam in lam2:
    pts = np.where(lams == lam)
    p, = mp.plot(mags, snrarr[pts[0],:][0], next(linecycler))
    curves.append(p)
mp.legend(curves, lam2, loc='lower left', ncol=2, title='IRCAL Filters')
mp.plot([min(mags),max(mags)],[5,5],'-',linewidth=2)
mp.text(np.max(mags)-0.5, 5.5, 'SNR=5')

mp.clf()
mp.yscale('log')
mp.grid(True)
mp.title('Spectroscopy - SNR vs. Magnitude for 300s Fowler-32')
mp.xlabel('Magnitude (Vega system)')
mp.ylabel('SNR')
curves = []
for lam in mlam2:
    pts = np.where(mlams == lam)
    p, = mp.plot(mags, msnrarr[pts[0],:][0], next(linecycler))
    curves.append(p)
mp.legend(curves, mlam2, loc='lower left', ncol=2, title='MOSFIRE Filters')
mp.plot([min(mags),max(mags)],[5,5],'-',linewidth=2)
mp.text(np.max(mags)-0.5, 5.5, 'SNR=5')


