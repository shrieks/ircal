import os
import numpy as np
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
            filename = fileroot + '/' + s[0:len('PlusDust')] + '.txt'
            scalefac = irdust
        else:
            filename = fileroot + '/' + s + '.txt'
            scalefac = 1.0

        # Check if file exists
        if os.path.isfile(filename):
            # These are 2-column files
            vec = np.loadtxt(filename)
            wavevec = vec[:,0]            # first column
            coatvec = scalefac * vec[:,1] # second column scaled
        else:
            # Find a better way to gracefully handle this
            # Right now file missing causes program to barf and stop
            # catastrophically
            print 'File ' + s + '.txt does not exist. Skipping...'

        # Convert to nm from microns for all files except:
        if 'H2RG_QE' not in s and 'PICNIC_QE' not in s and 'lick_ir' not in s:
            wavevec = wavevec * 1000.0

        # All files except IR_transmission have coating percentages
        # convert them to fractions
        # convert this to just check if any element in coatvec > 1
        # and divide by 100 if it is
        if 'IR_transmission' not in s:
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



def emissivity(lambdas, names, tempinC, dustfrac=0.0, dustemis = 1.0, 
               delta=0.0, fileroot='.'):
    """
    this relies on throughput.pro which ONLY works for 1 um and redder
    assumes everyting in <names> is at the same temperature
    input lambda is in nm
    temperature in C
    delta is subtracted from the reflectivity of the coating before it
    is returned by the throughput function
    dustfrac is a fraction of the surface with emissivity=1, possibly
    from dust or other contamination

    call this function multiple times and add the results for different temps
    """

    nlambs = len(lambdas)  # Number of wavelengths
    nsurfs = len(names)    # Number of surfaces

    # define array for output flux from optics
    output = np.zeros(nlambs, dtype='float64')
    thruput_eff = np.zeros(nlambs, dtype='float64') # not optional any more
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

    tempinK = tempinC + 273.0
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

    for i,s in enumerate(names):

        # Might need a more lax condition here to account for lower
        # case or 'withdust' etc. - substring Al?
        if 'luminum' in s:
        # if s=='Aluminum':
            thvec = throughput(lambdas, [s], delta=delta)
        else:
            thvec = throughput(lambdas, [s])

        # Find throughput after all surfaces, dimensionless
        # thvec = throughput(lambdas, names, delta=delta_i)
        # emissivity is 1-thruput at each lambda
        epsilon = 1.0 - thvec

        # deal with places where epsilon is less than 0.2?
        # is this necessary?
        # epsilon[np.where( epsilon < 0.2 )] = 0.2

        # account for particle contamination - assume particles have 
        # emiss = dustemis
        epsilon_eff = (1.0 - dustfrac)*epsilon + dustfrac*dustemis
        tau_eff     = 1.0 - epsilon_eff   # new effective transmissivity

        # apply effective emissivity
        onesurface = phpersmas * epsilon_eff

        # output of optic[i] = (output of prior optic)*(throughput of optic[i]
        #                      + emissivity of optic[i]
        output = output*tau_eff + onesurface

        if i==0:
            thruput_eff = tau_eff
        else:
            thruput_eff = thruput_eff * tau_eff

    
    return output, thruput_eff



