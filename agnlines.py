import numpy as np
import matplotlib.pyplot as mp
import asciitable as at
import simshspec as sp
from background import *
from matplotlib import rc
from itertools import cycle

def agnlines():

    R  = 500.0 # lambda/delta lambda
    sn = 5.0   # 5 sigma detection
    ft = 600.0 # fowler time in secs

    lines = ["-","--","-.",":","-","--"]
    linecycler = cycle(lines)

    # Wavelength scale from start of J to end of K filters
    xmin = 1100 # nm
    xmax = 2400

    # Rest frame wavelengths of interest
    pasb  = 1282.18072 # nm - Paschen beta line
    fe2fb = 1256.68    # [FeII]
    h2    = 2121.83    # H2
    brg   = 2166.12    # Br-gamma
    #lams = np.array((pasb, fe2fb, h2, brg))
    #lname= ['Paschen-$\beta$', '[FeII]', 'H$_2$', 'Br-$\gamma$']
    lams = np.array((pasb, fe2fb))
    lname= ['Paschen-$\beta$', '[FeII]']

    # Read in AGN flux file
    # From 2008 Landt et al. Mrk110 is a narrow-line Seyfert 1 galaxy
    # z = 0.035 (from NED); J mag = 13.9 (from 2MASS); V mag = 16.4 (Veron-Cetty)
    #mkn110 = np.loadtxt('spectra/AGN/mkn110.txt')
    #agnz = 0.035  # mrk110 has z=0.035
    #agnJ = 13.872   # J mag in Johnson system
    #agnH = 13.159   # H mag
    # convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
    #agnflux = mkn110[:,1]*1.0e5
    #agnlamb = mkn110[:,0]*1000.0/(1+agnz)   # table has wavelength in microns
                                            # shift lambda to rest frame

    # J-band magnitude range
    mags = np.arange(15.0, 21.0, 1.0)
    zrange = np.arange(0.01, 0.1, 0.01)

    agns = at.read('spectra/AGN/properties-edit.txt')
    fileroot = 'spectra/AGN/'
    plotroot = 'figs/agns/'

    # exposure time array redshift (rows) x wavelengths (columns) x mags
    exptarr = np.zeros((len(agns.z), len(zrange), lams.size, mags.size), 
                       np.float)
    snrarr  = np.zeros((len(agns.z), len(zrange), lams.size, mags.size), 
                       np.float)


    for n, gal in enumerate(agns.Name):
        #if n == 1:
        #    break
        print 'Reading AGN '+gal
        agnspec = np.loadtxt(fileroot+gal+'.txt')
        # convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
        agnflux = agnspec[:,1]*1.0e5
        # table has wavelength in microns, convert to nm & 
        # shift lambda to rest frame
        agnlamb = agnspec[:,0]*1000.0/(1+agns.z[n])   

        #exptarr = np.zeros((len(zrange), lams.size, mags.size), np.float)
        #snrarr  = np.zeros(len(zrange), lams.size, mags.size), np.float)

        for i, z in enumerate(zrange):

            # Shift spectrum by redshift
            agn = np.column_stack((agnlamb*(1+z), agnflux))
            # Feed through telescope with wavelengths of interest shifted too
            agnspec = sp.simspec(agn, lams*(1+z))

            for j, l in enumerate(lams):
                # Sky photons s-1 in diffraction limited core for filter
                Rsky  = agnspec['Rsky'][j]
                Rn2n  = agnspec['RN2n'][j]
                if Rsky < 0:
                    continue

                for k, m in enumerate(mags):
                    # m - Jmag = -2.5 log (Ragn/Rsrc)
                    # Source photons s-1 in slit given magnitude
                    Rsrc = agnspec['Rsrc'][j] 
                    if  Rsrc > 0:
                        Ragn = Rsrc * 10**(-(m-agns.J[n])/2.5)
                    else:
                        continue

                    # Solve the SNR quadratic equation for t
                    exptarr[n,i,j,k] = ((sn**2*(Ragn + Rsky) + 
                                         np.sqrt(sn**4*(Ragn + Rsky)**2 +
                                                 4*sn**2*Ragn**2*Rn2n))
                                        /(2*Ragn**2))

                    snrarr[n,i,j,k] = ((Ragn*ft)/
                                       np.sqrt(Ragn*ft + Rsky*ft + Rn2n))

#    return exptarr, snrarr
        fig = mp.figure()

        ax = fig.add_subplot(311)
        ax.set_xlim(xmin, xmax)
        ax.set_xticklabels([])
        ax.grid(True)
        ax.plot(agnlamb, agnflux)
        ax.set_ylabel(r'Flux - ergs s$^{-1}$ m$^{-2}$ nm$^{-1}$')
        ax.text(0.5, 0.9,'Spectrum and SNR for AGN '+gal, ha='center', 
                va='center', transform=ax.transAxes)

        ax2 = ax.twiny()
        lo_ticks = np.array([1400., 1600., 1800., 2000., 2200.])
        new_tick_locations = (lo_ticks-xmin)/(xmax-xmin)

        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(lo_ticks)
        ax2.set_xlabel(r"Wavelength $(nm)$")

        # Paschen beta plot
        ax1 = fig.add_subplot(312)
        ymin = 1e-2
        ymax = 1e3
        ax1.grid(True)

        ax1.set_title(r'Paschen-$\beta$ (Top), [FeII] (Bot) SNR v. $z$')
        ax1.set_xlim(xmin, xmax)
        #ax1.set_xlabel(r'Wavelength ($nm$)')
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
        #ax1.locator_params(axis='y',tight=True, nbins=4)
        ax1.set_yticks([1e-2,1,1e2])

        # Plot Paschen-beta SNR vs. z for each mag
        curves = []
        for k,m in enumerate(mags):
            p, = ax1.plot(pasb*(1+zrange), snrarr[n,:,0,k]+1e-9, lines[k])
            curves.append(p)

        ax1.plot([xmin,xmax],[5,5],'k-',linewidth=2)
        ax1.text(xmax-160, 5.5, 'SNR=5')
        #ax1.plot([H['start'],H['start']],[ymin,ymax],'--',color='0.75',linewidth=1)
        ax1.plot([H['start'],H['start']],[ymin,ymax],'b--',linewidth=1)
        ax1.plot([H['end'],H['end']],[ymin,ymax],'b--',linewidth=1)
        ax1.text(H['cwvl']*1e9, 0.1, 'H', color='b')
        ax1.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--',linewidth=1)
        ax1.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--',linewidth=1)
        ax1.text(Ks['cwvl']*1e9, 0.1, 'Ks', color='r')

        mp.legend(curves, mags, loc='upper left', title='J-band mags',
                  prop={'size':10})

        #  [FeII] plot
        ax3 = fig.add_subplot(313)
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
        ax3.set_yticks([1e-2,1,1e2])

        # Plot [FeII] SNR vs. z for each mag
        curves = []
        for k,m in enumerate(mags):
            p, = ax3.plot(fe2fb*(1+zrange), snrarr[n,:,1,k]+1e-9, lines[k])
            curves.append(p)

        ax3.plot([xmin,xmax],[5,5],'k-',linewidth=2)
        ax3.text(xmin+100, 5.5, 'SNR=5')
        ax3.plot([H['start'],H['start']],[ymin,ymax],'b--',linewidth=1)
        ax3.plot([H['end'],H['end']],[ymin,ymax],'b--',linewidth=1)
        ax3.text(H['cwvl']*1e9, 0.1, 'H', color='b')
        ax3.plot([Ks['start'],Ks['start']],[ymin,ymax],'r--',linewidth=1)
        ax3.plot([Ks['end'],Ks['end']],[ymin,ymax],'r--',linewidth=1)
        ax3.text(Ks['cwvl']*1e9, 0.1, 'Ks', color='r')

        #ax4.set_xticks(new_tick_locations)
        #ax4.set_xticklabels(tick_function(new_tick_locations))

        #mp.tight_layout()
        #mp.show()
        mp.savefig(plotroot+gal+'.png')

    return exptarr, snrarr
