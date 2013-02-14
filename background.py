import os
import numpy as np
import asciitable as at
from scipy.interpolate import InterpolatedUnivariateSpline # spline interpolator


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

        # Convert to nm from microns for all files except:
        if 'H2RG_QE' not in s and 'PICNIC_QE' not in s and 'lick_ir' not in s:
            wavevec = wavevec * 1000.0

        # All files except IR_transmission and IR trans read from tables
        # have coating percentages. Convert them to fractions
        # Change this to just check if any element in coatvec > 1
        # and divide by 100 if it is
        if 'IR_transmission' not in s and fromtable==False:
            coatvec = coatvec/100.0

        # cubic spline fit for now - check other fits later
        spline = InterpolatedUnivariateSpline(wavevec, coatvec, k=3)
        # k = 3rd degree smoothing spline

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
        if allemiss == 0:
            epsilon_eff = (1.0 - dustfrac)*epsilon + dustfrac*edust
        else:
            epsilon_eff = 1.0
        tau_eff     = 1.0 - epsilon_eff   # new effective transmissivity

        tempinK = temparrinC[i] + 273.0
        # Planck's blackbody formula - units W/(m^2 cm sr)
        intensity = ((hc2x2/lambdacm**5) * 
                     (1.0/(np.exp(hcoverk/(lambdacm*tempinK)) - 1.0)))
        # Divide by energy per photon = hc/lamba i.e. Joules or W-s
        # photons per second per m^2-steradian-wavelength(cm)
        phpersmsr = intensity/(hc/lambdacm)    # photons/(s m^2 cm sr)
        # 1 sr = 1 rad^2 = 206265 arcsec^2
        # photons per second m^2-nm-arcsec^2
        phpersmas = phpersmsr*100.0/(206265.0**2 * 1.0e9) 
                                               # photons/(s m^2 nm arcsec^2)


        # apply effective emissivity
        onesurface = phpersmas * epsilon_eff

        # output of optic[i] = (output of prior optic)*(throughput of optic[i]
        #                      + emissivity of optic[i]

        if i==0:
            output = skyemis*tau_eff + onesurface
            tau_sys = tau_eff
        else:
            output = output*tau_eff + onesurface
            tau_sys = tau_sys * tau_eff
    
    return output, tau_sys
