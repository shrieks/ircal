import os
import numpy as np
import asciitable as at
from scipy.interpolate import InterpolatedUnivariateSpline # spline interpolator
from scipy.interpolate import interp1d

# Define filter parameters - Mt Ham filters
J =  { 'name' : 'J', 
       'cwvl' : 1.238E-06,      # center wavelength in m
       'start': 1102.5,         # band start in nm
       'end'  : 1373.5,         # band end in nm
       'zpc'  : 3.031e+09       # Connie's zero point photons s-1 m-2
       }

H =  { 'name' : 'H',
       'cwvl' : 1.656E-06,      # center wavelength in m
       'start': 1508.0,         # band start in nm
       'end'  : 1804.0,         # band end in nm
       'zpc'  : 2.817e+09       # Connie's zero point photons s-1 m-2
       }

K =  { 'name' : 'K',
       'cwvl' : 2.195E-06,      # center wavelength in m
       'start': 1989.5,         # band start in nm
       'end'  : 2400.5,         # band end in nm
       'zpc'  : 1.5062e+09      # Connie's zero point photons s-1 m-2
       }

# Ks
Ks = { 'name' : 'Ks',
       'cwvl' : 2.150E-06,      # center wavelength in m
       'start': 1990.0,         # band start in nm
       'end'  : 2310.0,         # band end in nm
       'zpc'  : 1.5066e+09      # Connie's zero point photons s-1 m-2
       }

# K'
Kp =  { 'name' : 'K\'',
        'cwvl' : 2.115E-06,      # center wavelength in m
        'start': 1909.5,         # band start in nm
        'end'  : 2320.5,         # band end in nm
        'zpc'  : 1.658e+09       # Connie's zero point photons s-1 m-2
        }


# MOSFIRE filter set
MJ = { 'name' : 'MJ', 
       'cwvl' : 1.253E-06,      # center wavelength in m
       'start': 1153,           # band start in nm
       'end'  : 1353,           # band end in nm
       'zpc'  : 3.031e+09       # Connie's zero point photons s-1 m-2
       }

MH = { 'name' : 'MH',
       'cwvl' : 1.637E-06,      # center wavelength in m
       'start': 1466.5,         # band start in nm
       'end'  : 1807.5,         # band end in nm
       'zpc'  : 2.817e+09       # Connie's zero point photons s-1 m-2
       }

MK = { 'name' : 'MK',
       'cwvl' : 2.162E-06,      # center wavelength in m
       'start': 1920.5,         # band start in nm
       'end'  : 2403.5,         # band end in nm
       'zpc'  : 1.5062e+09      # Connie's zero point photons s-1 m-2
       }

# Ks
MKs ={ 'name' : 'MKs',
       'cwvl' : 2.147E-06,      # center wavelength in m
       'start': 1990.0,         # band start in nm
       'end'  : 2304.0,         # band end in nm
       'zpc'  : 1.5066e+09      # Connie's zero point photons s-1 m-2
       }

# Original Keck filters
# Define filter parameters
Jparams =  { 'cwvl' : 1.250E-06,      # center wavelength in m
             'start': 1170.0,         # band start in nm
             'end'  : 1330.0,         # band end in nm
             'zpc'  : 3.031e+09       # Connie's zero point photons s-1 m-2
             }

Hparams =  { 'cwvl' : 1.635E-06,      # center wavelength in m
             'start': 1490.0,         # band start in nm
             'end'  : 1780.0,         # band end in nm
             'zpc'  : 2.817e+09       # Connie's zero point photons s-1 m-2
             }

Kparams =  { 'cwvl' : 2.200E-06,      # center wavelength in m
             'start': 2030.0,         # band start in nm
             'end'  : 2370.0,         # band end in nm
             'zpc'  : 1.5062e+09      # Connie's zero point photons s-1 m-2
             }

# Ks
Ksparams=  { 'cwvl' : 2.150E-06,      # center wavelength in m
             'start': 1990.0,         # band start in nm
             'end'  : 2310.0,         # band end in nm
             'zpc'  : 1.5066e+09      # Connie's zero point photons s-1 m-2
             }

# K'
Kpparams=  { 'cwvl' : 2.120E-06,      # center wavelength in m
             'start': 1950.0,         # band start in nm
             'end'  : 2290.0,         # band end in nm
             'zpc'  : 1.658e+09       # Connie's zero point photons s-1 m-2
             }


def getcoatingdata(namelist, fileroot='./emissdatafiles', irdust=0.8):
    """
    read in coating data. 
    INPUT: 
    namelist: list of coatings. 
    fileroot: directory containing coating data files
    irdust: transmission coefficient (dafault = 0.8)

    RETURNS 
    splinelist: cubic spline fitting function of transmissivity vs lambda
                for each coating in namelist
    """
    
    # coating data files with wavelengths in nm
    nmfiles = ['H2RG_QE', 'PICNIC_QE', 'lick_ir', 'FM3_7', 'WFS_DichBS4-trans',
               'TT_DichBS1-refl', 'JHKpK']
    # coating data files with refl/trans in decimal 
    fracfiles = ['IR_transmission', 'JHKpK', 'mosfire_JHKKs']
    # Find number of names in list
    nnames = len(namelist)

    splinelist = []

    for i, s in enumerate(namelist):

        if 'PlusDust' in s:
            # Drop the 'PlusDust' suffix in the filename and set the scale
            # accordingly
            filename = fileroot + '/' + s[0:-len('PlusDust')] + '.txt'
            scalefac = irdust
        elif 'BULK' in s:
            # Drop the 'BULK' suffix in the filename
            s = s[0:-len('_BULK')]
            filename = fileroot + '/' + s + '.txt'            
        else:
            filename = fileroot + '/' + s + '.txt'
            scalefac = 1.0

        # Construct alternate wavevector for coatings with no files
        altwavevec = np.arange(1.0, 2.350, 0.05)
        # Reader flag names arrays with header column names as strings
        othercoats = at.read('keckngao/ngao_bkg_materials.txt',
                             Reader=at.CommentedHeader)
        fromtable = False

        # Check if file exists
        if os.path.isfile(filename):
            # These are 2-column files
            vec = np.loadtxt(filename)
            wavevec = vec[:,0]            # first column
            coatvec = scalefac * vec[:,1] # second column scaled
        elif s in othercoats['ID']:
            wavevec = altwavevec
            # Find the constant IR transmissivity value in file
            othertau = 1-othercoats['IRTRAN'][np.where(othercoats['ID'] == s)]
            coatvec = np.ones(altwavevec.size)*othertau[0]
            fromtable = True
        else:
            # Find a better way to gracefully handle this
            # Right now file missing causes program to barf and stop
            # catastrophically
            print 'File ' + s + '.txt does not exist. Skipping...'

        # Convert to nm from microns for all files except thos in nmfiles
        if s not in nmfiles:
            #print "Converting to nm"
            wavevec = wavevec * 1000.0

        # All files except IR_transmission and IR trans read from tables
        # have coating percentages. Convert them to fractions
        # Change this to just check if any element in coatvec > 1
        # and divide by 100 if it is
        if s not in fracfiles and fromtable==False:
            #print "converting to fraction"
            coatvec = coatvec/100.0

        # cubic spline fit for now - check other fits later
        spline = InterpolatedUnivariateSpline(wavevec, coatvec, k=3)
        # k = 3rd degree smoothing spline
        # check if the spline blew up & switch to linear fit if it did
        if np.isnan(np.sum(spline(wavevec))):
            #spline = InterpolatedUnivariateSpline(wavevec, coatvec, k=1)
            spline = interp1d(wavevec, coatvec, kind='linear')

        splinelist.append(spline)

    return splinelist


def throughput(lambdas, namelist, delta=0.0, fileroot='./emissdatafiles'):
    """
    uses getcoatingdata, and is ONLY meant to work at wavelengths in the
    range 800 nm - 2300 nm
    namelist is a list of coatings to be multipled to give a throughput 
    lambdas is in nm
    delta is subtracted from the reflectivity - usually applies to Al only
    """

    # populate list of arrays with wavelengths, coating and spline func
    # corresponding to each name in namelist
    sl = getcoatingdata(namelist)

    nlambs = len(lambdas)
    nnames = len(namelist)

    # set initial value to 1 at each lambda - perfect transmission
    tauvec = np.ones(nlambs, dtype='float64') 

    for i, s in enumerate(namelist):
        # Find the amount let through by each coating or layer for each lambda
        coatvalue = sl[i](lambdas) - delta
        # Multiply by the output vector
        tauvec = tauvec * coatvalue

    # return the total transmissivity of the system
    return tauvec



def emissivity(lambdas, names, skyemis, temparrinC, deltaarr, dustfrac=0.0, 
               edust = 1.0, allemiss = 0.0, fileroot='.'):
    """
    this relies on throughput.pro which ONLY works for 1 um and redder
    assumes everyting in <names> is at the same temperature
    lambda:   is wavelength array in nm
    names :   list of surface names for the telescope/AO/instrument
    skyemis:  sky emission in the lambdas wavelengths
    temparr:  array of all surface temperatures in C
    delta:    is subtracted from the reflectivity of the coating before it
              is returned by the throughput function          
    dustfrac: is a fraction of the surface with emissivity=1, possibly
              from dust or other contamination
    edust:    emissivity of dust
    set allemis to 1 to add background radiation of emissivity 1.0 from each
    surface without seting the throughput to zero. This is to account for, e.g.,
    the oversized cold stop in ircal

    """

    nlambs = len(lambdas)  # Number of wavelengths
    nsurfs = len(names)    # Number of surfaces

    # define array for output flux from optics and system transmissivity
    output = np.zeros(nlambs, dtype='float64')
    tau_sys = np.zeros(nlambs, dtype='float64')   # not optional any more
                                                  # return it always
    # convert lambda from nm to cm
    lambdacm = lambdas*100.0/1.0e9
    # pre-compute various constants 
    # 2hc^2 in W/(m^2 ster cm^(-4)))
    hc2x2 = 1.191044e-8     # e-16 in MKS
    # hc/k [K cm] - k is k_B (Boltzmann's constant)
    hcoverk = 1.438769      # e-2 in MKS
    # hc [W-cm-s]
    hc = 1.98645e-23        # e-25 in MKS

    for i,s in enumerate(names):

        # Find throughput after all surfaces, dimensionless
        thvec = throughput(lambdas, [s], delta=deltaarr[i])

        # emissivity is 1-thruput at each lambda
        epsilon = 1.0 - thvec

        # deal with places where epsilon is less than 0.2?
        # is this necessary?
        # epsilon[np.where( epsilon < 0.2 )] = 0.2

        # account for particle contamination - assume particles have 
        # emiss = edust
        # if the temperature of a sruface is very low, assume it's in
        # a dewar or similar enclosure and is not contaminated by dust
        if temparrinC[i] < -20.0:
            epsilon_eff = epsilon
        else:
            epsilon_eff = (1.0 - dustfrac)*epsilon + dustfrac*edust

        tau_eff     = 1.0 - epsilon_eff   # new effective transmissivity

        if allemiss == 1.0:
            epsilon_eff = epsilon*0.0 + 1.0 # set all to 1 in vector

        tempinK = temparrinC[i] + 273.0
        # Planck's blackbody formula - units W/(m^2 cm sr)
        intensity = ((hc2x2/lambdacm**5) * 
                     (1.0/(np.exp(hcoverk/(lambdacm*tempinK)) - 1.0)))
        # Divide by energy per photon = hc/lamba i.e. Joules or W-s
        # photons per second per m^2-steradian-wavelength(cm)
        phpersmsr = intensity/(hc/lambdacm)    # photons s-1 m-2 cm-1 sr-1
        # 1 sr = 1 rad^2 = 206265 arcsec^2
        # photons s-1 m-2 nm-1 arcsec-2
        phpersmas = phpersmsr*100.0/(206265.0**2 * 1.0e9) 


        # apply effective emissivity
        onesurface = phpersmas * epsilon_eff

        # output of optic[i] = (output of prior optic)*(throughput of optic[i]
        #                      + emissivity of optic[i]

        if i==0:
            # keep track of source photons and emissivity separately
            srcout = skyemis*tau_eff
            emisout= onesurface
            # add in emissivity, photons s-1 m-2 nm-1 arcsec-2
            output = srcout + onesurface
            tau_sys = tau_eff
        else:
            srcout = srcout*tau_eff
            emisout= emisout*tau_eff + onesurface
            output = output*tau_eff + onesurface
            tau_sys = tau_sys * tau_eff

        # inefficient, srcout is just skyemis * tau_sys
   
    # return output, tau_sys
    return output, srcout, emisout, tau_sys
