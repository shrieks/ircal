import numpy as np

def resamplesky(sky, R=500.0, endlambda=2500.0):
    """
    Resamples a given sky spectrum from constant delta-lambda to constant R.
    Inputs:
          sky    - 2 column vector with wavelength (col 0) and flux (col 1)                   
          R      - spectroscopy resolution: lambda/dlambda
    Output:
          resky  - resampled sky spectrum
                   2 column array. Column 1 is wavelength in nm
                                   Column 2 is flux in photons s-1 m-2 nm-1 arcsec-2
    
    Resampling is done by computing dlambda at the starting wavelength (lambda0) 
    and adding all photons between lambda0 and lambda0 + dlambda. Process repeats
    with lambda0 + dlambda as a starting point. Continues until the end point
    (currently 2500 nm)
    """

    startlambda = sky[0,0]
    lambdas = []
    flux    = []
    while startlambda < endlambda:
        lambdas.append(startlambda)
        dlambda = startlambda / R
        # find end point - startlambda + dlambda
        nextl, = where(sky[:,0] <= startlambda + dlambda)
        # count all photons between startlambda & startlambda + dlambda
        fluxsum = np.sum(sky[nextl[0]:nextl[-1], 1])
        flux.append(fluxsum)
        startlambda = sky[next[-1]+1,0]

    # Create an array of the right length with 2 columns
    resky = np.zeros([len(lambdas), 2], np.float)
    np.shape(resky)
    resky[:,0] = lambdas
    resky[:,1] = flux
    return resky
