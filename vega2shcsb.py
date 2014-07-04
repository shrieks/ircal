import numpy as np
import matplotlib.pyplot as mp
import simshane as ss
#from itertools import cycle

# Read in vega flux file
vega = np.loadtxt('spectra/alpha_lyr_stis_005.txt')
# convert from ergs s-1 cm-2 A-1 to ergs s-1 m-2 nm-1
vegaflux = vega[:,1]*1.0e5
vega[:,1] = vegaflux

vegalamb = vega[:,0]/10.0   # table has wavelength in angstroms
vega[:,0] = vegalamb

vegacsb = ss.simshane(vega, source='csb')
vegamcsb = ss.simshane(vega, source='csb', subaps=8)
#vegamcsb = ss.simshane(vega, source='csb', filts='mos')

mags = np.arange(14.0, 24.5, 0.5)

exptarr = np.zeros((mags.size, vegacsb['filter'].size), np.float)
snrarr = np.zeros((mags.size, vegacsb['filter'].size), np.float)

mexptarr = np.zeros((mags.size, vegamcsb['filter'].size), np.float)
msnrarr = np.zeros((mags.size, vegamcsb['filter'].size), np.float)

sn = 5.0
ft = 300.0 # fowler time in secs

for j, f in enumerate(vegacsb['filter']):
    ef = open('tables/sharcexpt-csb'+f+'.txt', 'w')
    sf = open('tables/sharcsnr-csb' +f+'.txt', 'w')
    mef = open('tables/mosfexpt-csb'+f+'.txt', 'w')
    msf = open('tables/mossnr-csb' +f+'.txt', 'w')

    for i, m in enumerate(mags):
        # m - 0 = -2.5 log (Rstar/Rvega)
        # Star photons s-1 in psf coer given magnitude
        Rstar = vegacsb['Rsrc'][j] * 10**(-m/2.5)
        mRstar = vegamcsb['Rsrc'][j] * 10**(-m/2.5)
        # Sky photons s-1 in diffraction limited core for filter
        Rsky  = vegacsb['Rsky'][j]
        mRsky  = vegamcsb['Rsky'][j]
        # Solve the SNR quadratic equation for t
        exptarr[i,j] = ((sn**2*(Rstar + Rsky) + 
                         np.sqrt(sn**4*(Rstar + Rsky)**2 +
                                 4*sn**2*Rstar**2*vegacsb['RN2n'][j]))
                        /(2*Rstar**2))
        ef.write('{0:.1f} & {1:.3e} & {2:.3e} & {3:.3e} & {4:.3e}\\\\\n'.format(m, exptarr[i,j], Rstar*exptarr[i,j], Rsky*exptarr[i,j], vegacsb['RN2n'][j]))
        ef.write('\hline\n')

        mexptarr[i,j]= ((sn**2*(mRstar + mRsky) + 
                         np.sqrt(sn**4*(mRstar + mRsky)**2 +
                                 4*sn**2*mRstar**2*vegamcsb['RN2n'][j]))
                        /(2*mRstar**2))
        mef.write('{0:.1f} & {1:.3e} & {2:.3e} & {3:.3e} & {4:.3e}\\\\\n'.format(m, mexptarr[i,j], mRstar*mexptarr[i,j], mRsky*mexptarr[i,j], vegamcsb['RN2n'][j]))
        mef.write('\hline\n')


        snrarr[i,j] = ((Rstar*ft)/
                       np.sqrt(Rstar*ft + Rsky*ft + vegacsb['RN2n'][j]))
        sf.write('{0:.1f} & {1:.2f} & {2:.3e} & {3:.3e} & {4:.3e}\\\\\n'.format(m, snrarr[i,j], Rstar*ft, Rsky*ft,vegacsb['RN2n'][j] ))
        sf.write('\hline\n')

        msnrarr[i,j]= ((mRstar*ft)/
                       np.sqrt(mRstar*ft + mRsky*ft + vegamcsb['RN2n'][j]))
        msf.write('{0:.1f} & {1:.2f} & {2:.3e} & {3:.3e} & {4:.3e}\\\\\n'.format(m, msnrarr[i,j], mRstar*ft, mRsky*ft,vegamcsb['RN2n'][j] ))
        msf.write('\hline\n')

    ef.close()
    sf.close()
    mef.close()
    msf.close()

#lines = ["-","--","-.",":"]
lstyles = ['r-', 'b-', 'g-', 'c-']
lstylei = ['r--', 'b--', 'g--', 'c--']
lstylem = ['r-.', 'b-.', 'g-.', 'c-.']
#linecycler = cycle(lines)

mp.clf()
mp.yscale('log')
mp.grid(True)
mp.title('Exposure time vs. Magnitude - Const. Surf. Brightness')
mp.xlabel('Magnitude (Vega system)')
mp.ylabel('Exposure time (s)')
sj, = mp.plot(mags, exptarr[:,3], lstyles[3])
sh, = mp.plot(mags, exptarr[:,1], lstyles[1])
sk, = mp.plot(mags, exptarr[:,2], lstyles[2])
sks, = mp.plot(mags, exptarr[:,0], lstyles[0])

mj, = mp.plot(mags, mexptarr[:,3], lstylem[3])
mh, = mp.plot(mags, mexptarr[:,1], lstylem[1])
mk, = mp.plot(mags, mexptarr[:,2], lstylem[2])
mks, = mp.plot(mags, mexptarr[:,0], lstylem[0])

mp.plot([np.min(mags),np.max(mags)],[300,300],'--')
mp.text(np.min(mags)+0.1, 350, '300s Fowler')
mp.text(np.min(mags)+0.1, 150, 'Exposure')
leg1 = mp.legend([sj,sh,sk,sks],['J','H','K','Ks'], loc='upper left', 
                 title='ShARCS 16')
mp.legend([mj,mh,mk,mks],['J','H','K','Ks'], loc='upper center', 
          title='ShARCS 8')
mp.gca().add_artist(leg1)
#mp.savefig('figs/exptVmag+mos-csb.png')
mp.savefig('figs/web/exptVmag-csb.png')
#mp.show()

mp.clf()
mp.yscale('log')
mp.grid(True)
mp.title('SNR vs. Magnitude - Const. Surf. Brightness')
mp.xlabel('Magnitude (Vega system)')
mp.ylabel('SNR')

sj, = mp.plot(mags, snrarr[:,3], lstyles[3])
sh, = mp.plot(mags, snrarr[:,1], lstyles[1])
sk, = mp.plot(mags, snrarr[:,2], lstyles[2])
sks, = mp.plot(mags, snrarr[:,0], lstyles[0])

mj, = mp.plot(mags, msnrarr[:,3], lstylem[3])
mh, = mp.plot(mags, msnrarr[:,1], lstylem[1])
mk, = mp.plot(mags, msnrarr[:,2], lstylem[2])
mks, = mp.plot(mags, msnrarr[:,0], lstylem[0])

mp.plot([np.min(mags), np.max(mags)],[5,5],'--')
mp.text(np.min(mags)+0.1, 6, 'SNR=5')

leg1 = mp.legend([sj,sh,sk,sks],['J','H','K','Ks'], loc='lower left', 
                 title='ShARCS 16')
mp.legend([mj,mh,mk,mks],['J','H','K','Ks'], loc='lower center', 
          title='ShARCS 8')
mp.gca().add_artist(leg1)
#mp.savefig('figs/snrVmag+mos-csb.png')
mp.savefig('figs/web/snrVmag-csb.png')
#mp.show()
