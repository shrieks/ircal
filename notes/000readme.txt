mk_skybg_zm_16_10_ph.dat - sky background model file

agnavg.py - Takes wavelengths of interest for AGNs, list of AGNs in
	    spectra/AGN/properties-edit.txt. For each AGN in list, scales 
	    continuum to J=15th mag, interpolates to list of wavelengths 
	    (lambdas). Takes median flux of all AGNs as representative 
	    spectrum for AGNs. This rep spectrum is red-shifted through a 
	    range of values and fed through shane spectroscopy simulator
	    (simshspec). Solves for exposure time (to get to SNR=5) and SNR
	    for lines of interest and plots them. 

agnlines.py - Reads in spectra for each AGN in list above and makes plots 
              for each individually rather than taking shifted median spectrum

coldstop.py - calculations to figure the effect of the oversized cold stop
	      and the extra flux allowed in by it

cosmocalc.py - cosmology calculator from web cosmology calculator

makebbsky.pro - adds (scalefac)*(1 - transmissionvec) * Planck function     
	      	of temp T to existing skyvec
		wavevec is in nm
		transmissionvec and skyvec use wavevec too
		Uses convgauss function - gaussian smooths data
		Just there as an option, no particular reason. 
		Connie says this should not be used - attempt to add
		blackbody characteristics to sky background

mintime.py - minimum exposure time needed for sky background to overcome
             read noise and for source photons to equal sky background (?)

mkskybg.py - Sky background for Mauna Kea sky models in various filters

peakdetect.py - 3rd party routine to find local peaks in 1-D array

simIrcal.py - emissivity and throughput for Shane+IRCAL 

simkeck.py - Keck emissivity, throughput and sky background

simshane.py - ShARCS simulation for imaging point source and constant surface 
	      brightness patch (add in brightness profile)

simshspec.py - ShARCS simulation for spectroscopy

vega2shane.py - Send point source through ShARCS 
vega2spec     - Send pt source through ShARCS spectro
vega2shcsb - CSB patch through ShARCS
vega2spcsb - CSB patch through ShARCS spectro

getcoatingdata.pro/py -  Takes a list of surfaces, searches for data
		         files corresponding to those surfaces,
		         converts all wavelengths to nm, all
		         coefficients to fractions, then fits a spline
		         to the coating data vs. lambda curve and
		         returns the wavelength vectors, coating
		         vectors and spline function
			 At some point test other fits like linear,
		         poly etc. 

throughput.pro/py - Takes a list of coatings and wavelengths
		    gets coating data and finds coating vec
		    subtracts delta value from final coating value
		    Multiplies coat value at each lambda if given a
		    list
		    If given a single name, gives transmissivity for
		    that layer
		    Output starts at 1 for all lambdas, each coating
		    reduces amount. FInal answer should be less than 1
		    Valid in 800-2300 nm range (for coatings now)
		    Needs to work from 1.1 um to 2.45 um for full JHK coverage

		    What does a negative value mean?
		    Negative values are bad spline fit points. Ignore

		    delta value is a parameter to reduce Al trans because
		    not perfect or pristine

emissiter.py - takes a list of optics and finds throughput & emissivity
	       for they system by finding tranmissivity at each layer
	       (throughput/getcoatingdata) and adding in
	       blackbody*emissivity. 
	       Given photons in, gives photons out for the system
	       Computes overall transfer function or transmissivity of
	       the system 

emissivity.py seems to give the right output, but methodology is wrong
emissivity2.py is mod of above and output looks right, spline is off by a bit

emissiter seems to make the most sense from rewriting (getcoatingdata
and throughput combined into one program - getcoatdata that just takes
one layer and wavelength list, computes 3rd order spline fit to
coating data and then returns transmissivity at each wavelength from
spline fit)

makeskybg.py - own program. Must use lick sky or Cerro Tololo or Mauna Kea 
	       (mk) sky. Zeropoints in ngao_bkg_filters - photon flux is
	       ZP - reference flux for 0 mag object. Filter curves in
	       same keckngao directory - apply filter and narrow graph
	       See gemini.edu page under Sciops for sky files for
	       different water vapor and airmass. Sky bg has
	       emissivity

background.py - is the shizzle. Does getcoatingdata (adds in capacity to
	        read from Keck material file and apply flat IRtrans),
		throughput and emissivity

* Add in source+sky photons at input of telescope. 
* Zero pts. for J,H,K in photons/s-m^2-nm-arcsec^2 in ngao_bkg_filters.dat

* Take system thruput at f(lambda); background as f(lambda)
* Signal per lambda bin
* Bg per lambda bin 

+++++++++++++++++++++++++++++++++++
For my MK sky model, airmass 1.0:
K' sky bg is 14.77 mag/arcsec^2
Ks sky bg is 14.80 mag/arcsec^2
K sky bg is 14.87 mag/arcsec^2
H sky bg is 13.80 mag/arcsec^2
J sky bg is 16.08 mag/arcsec^2
+++++++++++++++++++++++++++++++++++

Using nearIR_skybg_16_15_r5.dat:
KAON 501	zero-points /am		mk_skybg_zm_16_10_ph.dat:
Table 1
J 16.01		3.031e9	    16.0017	16.076
H 13.78		2.817e9	    13.7964	13.805
K 14.91		1.505e9     14.8642	14.875

mksky = np.loadtxt('skybgfiles/mk_skybg_zm_16_10_ph.dat')
dlam = 0.02
lambdas = mksky[:,0]
skyflux = mksky[:,1]
jidx  = np.where((lambdas > 1170) & (lambdas < 1330))
hidx  = np.where((lambdas > 1490) & (lambdas < 1780))
kidx  = np.where((lambdas > 2030) & (lambdas < 2370))
jsky = dlam*skyflux[jidx].sum()
-2.5*log10(jsky/Jzpc)

Connie's notes with my output interleaved
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
readcol, 'mk_skybg_zm_16_10_ph.dat', /silent, format='f,f', skylamall, skyfluxall

keckao = ['Aluminum', 'Aluminum', 'Aluminum', 'Aluminum', 'X1_Silver_Extrap', 'Aluminum', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap',  'X1_Silver_Extrap', 'TT_Dichroic_Refl', 'X1_Silver_Extrap',  'X1_Silver_Extrap']

keckaotlist = [2.6, 2.6, 2.6, 5,5,5,5,5,5,5,5,5,5]

deltalist = [0.0265, 0.0265, 0.0265, 0.0265, 0.0, 0.0265, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

tkeck = background(lambda, keckao, keckaotlist, skyflux, deltalist, dustfac=0.0, edust=1.0, total_trans=ttrans)

avtran = total(ttrans[kidx])/(n_elements(kidx))
kflux = (total(tkeck[kidx])*0.02)/avtran
print, kflux
print, -2.5*alog10(kflux/kzp)

For dustfrac = 0, K mag of background at airmass 1.0 is 13.09
for dustrac = 0.008, the ngao kludge value, BG is 12.84
(I think this is leaving out the last two surfaces, check)

=======================================
MK, am=1.0, df=0.0, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
upto dich: 16.074 13.808 12.997 

MK, am=1.0, df=0.008, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
upto dich: 16.074 13.807 12.729 
=======================================

NGAO used a different sky model, at 1.5 airmasses and different sampling (get
that right in the integration of the flux below if you recheck this! 0.02 ->
0.1), nearIR_skybg_16_15_r5.dat

My model, without anything after the TT dichroic because the bulk contribution
is trivial and the Keck NGAO ar coatings are nearly perfect, with the sky at 1.5
airmasses gives: 
K bg = 12.98 for dustfrac=0
K bg = 12.75 for dustfac = 0.008

=======================================
NG, am=1.5, df=0.0, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
upto dich: 15.999 13.799 12.995 

NG, am=1.5, df=0.008, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
upto dich: 15.999 13.798 12.728 
=======================================

the ngao study used higher aluminum reflectivity by 0.01 in K but slightly
lower protected silver by 0.005 So:
deltalist = [0.0165, 0.0165, 0.0165, 0.0165, 0.005, 0.0165, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]

If I use those deltas I get:
K bg of 13.03 with dustfrac=0.0 
K bg 12.80 with dustfrac=0.008

===============================================
NG, am=1.5, df=0.0, aldelta = 0.0165, otherdelta=0.005:   
              J       H       K
upto dich: 15.999 13.799 13.015 

NG, am=1.5, df=0.008, aldelta = 0.0165, otherdelta=0.005:   
              J       H       K
upto dich: 15.999 13.798 12.746 
===============================================

so even less in line with their results. 

If I add two more silver surfaces I get the NGAO BG values, 12.9 for
dustfrac=0, 12.65 for dustfrac=0.008  But I think that's unfair to the Keck
NGAO model.  

A tweak to any of the coatings would make the same kind of change, but should
check the real ngao model too.

keckao = ['Aluminum', 'Aluminum', 'Aluminum', 'Aluminum', 'X1_Silver_Extrap', 'Aluminum', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 'X1_Silver_Extrap', 'TT_Dichroic_Refl', 'CAF2_BULK', 'NIR_AR', 'NIR_AR', 'CAF2_BULK' ]

keckaotlist = [2.6, 2.6, 2.6, 5,5,5,5,5,5,5,5,5,5,5,5]

deltalist = [0.0265, 0.0265, 0.0265, 0.0265, 0.0, 0.0265, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0]

tkeck = background(lambda, keckao, keckaotlist, skyflux, deltalist, dustfrac=0.0, fileroot=IRCALROOT, total_trans=ttrans)

for airmass 1.0, 
K bg is 13.06 mag/arcsec^2 for dustfrac 0.0
K bg is 12.73 mag/arcsec^2 with dustfrac 0.008 

==================================================
MK, am=1.0, df=0.0, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
incl bulk: 16.073 13.807 13.016

MK, am=1.0, df=0.008, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
incl bulk: 16.073 13.806 12.706
==================================================

Try with the sky model at airmass 1.5:
K bg is 12.95 mag/arcsec^2 with dustfrac=0.0
K bg is 12.65 mag/arcsec^2 with dustfrac=0.008

================================================
NG, am=1.5, df=0.0, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
incl bulk: 15.999 13.798 13.014

NG, am=1.5, df=0.008, aldelta = 0.0265, otherdelta=0.0:   
              J       H       K
incl bulk: 15.999 13.797 12.705
================================================


* ngao_bkg.pro has mat struct that reads in materials and calculates
  nir trans and emiss - check and see how they do it
* Check throughput with Keck numbers and verify

* Shane materials list - make up from Reni's spreadsheet
- Hamilton sky - Elli
  * Ellie sent Gemini sky mags
  x Will check for Gemini sky measurements in March
  * Connie's modified Cerro Pachon sky gives closest to Gemini skybg numbers
    below
* Pt. source through (any mag) - flux calibrated spectrum eventually
* Make up a spectrum, flat in F_nu? Flux goes as 1/lambda
  All photons from source collected by primary (over area pi*(D/4)^2)
  end up focused in Airy rings. Core has 0.86*Strehl*Nph(source)
  Profile & SNR for pt. source
  - Take Airy disk sized area at focus of telescope (instrument for example)
    Calculate sky photons in that area per wavelength bin by integrating 
    flux over it.
  - Take a point source of some magnitude above atmosphere - this will
    give flux - photons s-1 m-2 nm-2. Propagate through atmosphere  
    (get atmosphere transmission) and then figure out how many photons
    at focus in Airy disk
  - Ratio of source photons/Sqrt(Source+sky photons) = SNR
    Poisson stats in noise limit of photons
    Add in system noise or read noise

* SNR too high! * Check your flux conversion.
      	  	- Check if SNR scales as expected with PSF (expand Airy
		  disk to 1.5 arcsec), telescope dia and src magnitude (yes)
  - Confusion cleared up:
    Source flux * avg trans is source flux at instrument
    avg trans is needed to calculate magnitude above atmosphere wrt mag0
    source - so divide output flux by avg trans to get photons incident to
    telescope
    This mag is separate from SNR at the detector - Detector only cares 
    about actual photons received there. 

- Randomness that came from changing SNR to include emissivity in denom
  and checking to make sure emissivity was same for all sources
  - If sky total output is included - it includes emissivity
    Emissivity units are photons s-1 m-2 nm-1 arcsec-2
    Has to be included only with sky (the multiply by airy area * 0.86 * strehl)
    Have to find a way to exclude emissivity from source going through
    * Can simplify telescope - run sky through to get sky + emissivity
      and average transmission
    * Pt source flux at detector is pt source influx * average transmission
    * SNR = (Nphotons source at det)/(Nphotons src + Nphotons sky&emis at det)
    * background is fine. 
      mod simkeck to send sky through background. Put filter info where?
      In background? telesim not needed
    - can actually recreate telesim with generic part of simkeck & shane
    
Tools to plot thruput curve
Exp time calculator from Brad's code

Shane sky background:
--------------------
From Gemini spectrograph manual 
(http://mtham.ucolick.org/techdocs/instruments/gemini/gemini_summary_table.html)
Gemini is at 2-ary cassegrain focus, so only looks at 1-ary and 2-ary mirrors.
Good check to find transmission of AO system - reliable sky bg measurement.
Check zero points for Gemini and IRCAL?

Gemini Backgrounds:		IRCAL backgrounds:
J = 16.0 mag per square arcsec	J 16 mag/arcsec^2   
H = 14.0 mag per square arcsec	H 14.4 mag/arcsec^2 
K = 13.0 mag per square arcsec	K 9.3 mag/arcsec^2  
				Ks 10.3 mag/arcsec^2
				 
Cerro Pachon (am 1.0):
0.008,0.0165,0.005 0.008,0.0265,0.0 0.0,0.0265,0.0 0.0,0.0165,0.005
16.103		   16.103	    16.103	   16.103
13.823		   13.823	    13.825	   13.824
12.093             12.194	    12.555	   12.435

Connie's mh:	   
16.103	 	   16.103	    16.103	   16.103
13.821		   13.821	    13.822	   13.822
11.606		   11.669	    11.878	   11.812

MK (am 1.0)
15.370		   15.370	    15.370	   15.370
13.065		   13.065	    13.065	   13.065
12.063		   12.161	    12.509	   12.394

Zeropoints (1 DN/s):	
J = 21.4
H = 20.9
K = 20.4

Connie's results:
----------------
My final, adjusted sky model predicts the Mt. Hamilton sky backgrounds:
J sky bg is 16.11 mag/arcsec^2
H skyi bg is 13.82 mag/arcsec^2
K' sky bg is 13.19 mag/arcsec^2
Ks sky bg is  12.97 mag/arcsec^2
K sky bg is 12.48 mag/arcsec^2

With older materials list:
 newshaneao = ['Aluminum', 'Aluminum', 'NaHG', 'NaHG', 'ALPAO', 'NaHG', 'TT_Dichroic_Refl', 'NaHG','MEMS_window', 'LGSWFS_Dich', 'FSG98', 'FSG98', 'ARIR', 'ARIR', 'FSG98', 'FSG98', 'FSG98', 'FSG98', 'H2RG_QE']

df = 0.0:
background mags:	my results:  throughput:		my results: 
K mag 11.93		11.967	     K 0.56			0.556	    
Ks mag 12.45		12.49	     Ks 0.56			0.559	    
Kp mag 12.70			     Kp 0.56				    
H mag 13.76		13.821	     H 0.52			0.522	    
J mag 16.04		16.104	     J 0.48			0.484       

df = 0.05:
background mags:	my results:
K mag 11.02		11.053
Ks mag 11.56		11.596
Kp mag 11.83		
H mag 13.70		13.815
J mag 15.98		16.104

newshaneao = ['Aluminum', 'Aluminum', 'LickCoating', 'LickCoating', 'ALPAO', 'LickCoating', 'TTLGS_R', 'LickCoating','MEMS_window', 'MEMS_window', 'BAU', 'MEMS_window', 'MEMS_window','WFS_T', 'FSG98', 'FSG98', 'ARIR', 'ARIR', 'BAU', 'BAU', 'ColdStop', 'FSG98', 'FSG98', 'H2RG_QE']

df = 0.0
Background:	C surfs:		my surfs:
K mag 11.92	11.956		       12.555  
Ks mag 12.45	12.488		       13.08
Kp mag 12.70				
H mag 13.76	13.821		       13.825
J mag 16.04	16.102		       16.103

Throughput:  my results:
J 0.46	     0.464
H 0.51	     0.508
Ks 0.53      0.533

df = 0.05	C surfs	My surfs	df = 0.06  my surfs:
K mag 10.76	10.789	10.84		10.57	   10.575
Ks mag 11.32	11.342	11.396		11.126	   11.133
Kp mag 11.59
H mag 13.70	13.812	13.814		13.809	   13.81
J mag 15.98	16.102	16.103		16.102	   16.103

The IRCAL manual gives the following measured background values:
J 16 mag/arcsec^2
H 14.4 mag/arcsec^2
Ks 10.3 mag/arcsec^2
K 9.3 mag/arcsec^2 

From Shane spreadsheet:
Note in IRCAL user manual says K band sky is 13'th mag/arcsec^2 and 1-2 mag fainter in J and H			
			
			
wvl	mag/as^2	ph/s/band/as^2/m^2	
	reported	reported	model
1200-1300	15	8900		2672.63554180365
1500-1800	14	7786.8479376797	27074.5638392369
2000-2300	13	8833.4028227227	2398.96755720328
			
mJ=0 star	0	8.90E+09	
mH=0 star	0	3.10E+09	
mK=0 star	0	1.40E+09	K-short	       

----------------------------------------------------------
Including excess emissivity through cold stop:

IDL code:
Old AO only
J        1160.0393       16.042783
H        8813.3453       13.761615
K        28457.091       11.809230
Ks        17720.601       12.323799
Kp        15586.771       12.567071

Old AO + bg
J        1071.9688       16.128509
H        8440.1054       13.808597
K        93549.716       10.517100
Ks        56811.833       11.058898
Kp        48713.354       11.329841

Python:
In [37]: run coldstop.py
I              J       H       K      Ks       Kp
sky  mag : 16.105 13.819 11.809 12.325 12.570	     # not icluding coldst
sky  flux : 1095.101 8358.689 28450.434 17693.462 15542.351
              J       H       K      Ks       Kp     # incl coldstop extra
sky  mag : 16.105 13.805 10.448 10.990 11.261
sky  flux : 1095.113 8465.918 99654.626 60520.836 51908.046

In [38]: run vega2shane # includes ircal
              ['Ks'     'H'      'K'     'J']        # not icluding coldstop
sky omag :  [ 12.324  13.822  11.809  16.107]

df = 0.0, sky is mh
              ['Ks'     'H'      'K'     'J']        # incl coldstop extra
sky omag :  [ 10.988  13.808  10.448  16.107]
df = 0.0, sky is cp
sky omag :  [ 11.134  13.81   10.596  16.107]

df = 0.008, sky is mh
sky omag :  [ 10.905  13.806  10.365  16.107]
df = 0.008, sky is cp
sky omag:   [ 11.04   13.809  10.501  16.107]

df = 0.07, sky is mh
sky omag :  [ 10.221  13.787   9.676  16.107]
df = 0.07, sky is cp
sky omag :  [ 10.291  13.79    9.746  16.107]

df = 0.07, sky is mh, al = 0.0165, odel = 0.005
sky omag :  [ 10.2    13.786   9.654  16.107]
df = 0.07, sky is cp, al = 0.0165, odel = 0.005
sky omag :  [ 10.268  13.789   9.723  16.107]

Find:
* Combo of dust and delta modification that will approximate ircal numbers
(Aluminum and other deltas)
* Throughput for ircal
* Sky mag at top of telescope with CP and MH skies? - Find numbers
- Should sky photons be propagated through sky transmission vector?
  - Ask Gemini Sci Ops!


* Keck & Lick
  ===========
  List surfaces
  refl & trans curves
  sky model & transmission

  Sky mag (at top of tele)
  background & thruput/sensitivity (NG Keck AO, Shane+IRCAL, Shane+new IR)

Connie suggests (Apr 12, 2013):
Stick to matching Ks band mag. According to Mt Hamilton LIRC II manual, K' flux
can vary by 2x (0.75 mag) with 9 C variation in temp. 

Reprod 10.8 Ks mag and thruput with df=0.03-0.04 and Al delta=0.0265
No df thruput is > measured for IRCAL. df=0.05, delta=0.02 thruput is 
50% < measured.
 

Table: SNR, R*, Rsky*npix, RN^2*npix/Nreads - differing reads 2,4,8,16,32 reads
(1.45 s for 2 reads). Fastest exposure is 1.45s (2 read times per pixel)
Pick number reads where src or sky dominate over RN term


