import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline # spline interpolator


def getcoatdata(lambdas, name, delta=0.0, irdust=0.8, 
                   fileroot='./emissdatafiles'):
    """
    read in coating data. 
    INPUT: 
    name    : coating 
    fileroot: directory containing coating data files
    irdust: transmission coefficient (dafault = 0.8)

    RETURNS 
    spline: cubic spline fitting function

    ONLY meant to work at wavelengths in the
    range 800 nm - 2300 nm
    """

    if 'PlusDust' in name:
        # Drop the 'PlusDust' suffix in the filename and set the scale
        # accordingly
        filename = fileroot + '/' + name[0:len('PlusDust')] + '.txt'
        scalefac = irdust
    else:
        filename = fileroot + '/' + name + '.txt'
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
        print 'File ' + name + '.txt does not exist. Skipping...'

    # Convert to nm from microns for all files except:
    if 'H2RG_QE' not in name and 'PICNIC_QE' not in name and 'lick_ir' not in name:
        wavevec = wavevec * 1000.0

    # All files except IR_transmission have coating percentages
    # convert them to fractions
    # convert this to just check if any element in coatvec > 1
    # and divide by 100 if it is
    if 'IR_transmission' not in name:
        coatvec = coatvec/100.0

    # spline = interp1d(wavevec, coatvec, kind='cubic')
    spline = InterpolatedUnivariateSpline(wavevec, coatvec, k=3)
    # k = 3rd degree smoothing spline

    transvec = spline(lambdas) - delta

    # return the transmissivity vector of the layer
    return transvec


def emissivity(lambdas, names, temp1, temp2='unset',n_at_temp1='unset', 
               dustfrac=0.0, dustemis = 1.0, delta=0.0, fileroot='.'):
    """
    this works for 1 um and redder
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

    # change output to set to sky background i.e. input
    output = np.zeros(nlambs, dtype='float64')
    thruput_eff = np.array(nlambs, dtype='float64') # not optional any more
                                                    # return it always

    if n_at_temp1 == 'unset':
        if temp2 != 'unset':
            n_at_temp1 = 3
            print "Number of surfaces at temp1 not specified"
            print "Assuming Keck system of 3 surfaces at temp1"
        else:
            n_at_temp1 = nsurfs
            print "Applying temp1 to all surfaces"

    if n_at_temp1 > nsurfs or n_at_temp1 < 0:
        n_at_temp1 = nsurfs
        print "Number of surfaces at temp1 not realistic"
        print "Setting equal to number of surfaces (size of name list)"

    # convert lambda from nm to cm
    lambdacm = lambdas*100.0/1.0e9
    # pre-compute various constants 
    # 2hc^2 in W/(m^2 ster cm^(-4)))
    hc2x2 = 1.191044e-8     # e-16 in MKS
    # hc/k [K cm] - k is k_B (Boltzmann's constant)
    hcoverk = 1.438769      # e-2 in MKS
    # hc [W-cm-s]
    hc = 1.98645e-23        # e-25 in MKS

    temp1K = temp1 + 273.0
    # Planck's blackbody formula - units W/(m^2 cm sr)
    intensity1 = ((hc2x2/lambdacm**5) * 
                  (1.0/(np.exp(hcoverk/(lambdacm*temp1K)) - 1.0)))
    # Divide by energy per photon = hc/lamba i.e. Joules or W-s
    # photons per second per m^2-steradian-wavelength(cm)
    phpersmsr1 = intensity1/(hc/lambdacm)    # photons/(s m^2 cm sr)
    # 1 sr = 1 rad^2 = 206265 arcsec^2
    # photons per second m^2-nm-arcsec^2
    phpersmas1 = phpersmsr1*100.0/(206265.0**2 * 1.0e9) 
                                             # photons/(s m^2 nm arcsec^2)

    if temp2 != 'unset':
        temp2K = temp2 + 273.0
        intensity2 = ((hc2x2/lambdacm**5) * 
                      (1.0/(np.exp(hcoverk/(lambdacm*temp2K)) - 1.0)))
        phpersmsr2 = intensity2/(hc/lambdacm)
        phpersmas2 = rphpermsr2*100.0/(206265.0**2 * 1.0e9)

    for i,s in enumerate(names):

        # Might need a more lax condition here to account for lower
        # case or 'withdust' etc. - substring Al?
        # if 'luminum' in s:
        if s=='Aluminum':
            trans = getcoatdata(lambdas, s, delta)
        else:
            trans = getcoatdata(lambdas, s)

        # emissivity is 1-thruput at each lambda
        epsilon = 1.0 - trans

        # deal with places where epsilon is less than 0.2?
        # is this necessary?
        # epsilon[np.where( epsilon < 0.2 )] = 0.2

        # account for particle contamination - assume particles have 
        # emiss = dustemis
        epsilon_eff = (1.0 - dustfrac)*epsilon + dustfrac*dustemis
        tau_eff     = 1.0 - epsilon_eff   # new effective transmissivity

        # set appropriate photons per second per m^2 per arcsec^2 per lambda
        if i > n_at_temp1 and temp2 != 'unset':
            phpersmas = phpersmas2     # surfaces at temp2
        else:
            phpersmas = phpersmas1     # surfaces at temp1

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
