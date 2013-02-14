import numpy as np
from getcoatdata import * 

def throughput(lambdas, namelist, delta=0.0, fileroot='./emissdatafiles'):
    """
    uses getcoatingdata, and is ONLY meant to work at wavelengths in the
    range 800 nm - 2300 nm
    namelist is a list of coatings to be multipled to give a throughput 
    lambdas is in nm
    delta is subtracted from the reflectivity 
    """

    # populate list of arrays with wavelengths, coating and spline func
    # corresponding to each name in namelist
    wl, cl, sl = getcoatingdata(namelist)

    nlambs = len(lambdas)
    nnames = len(namelist)

    outvec = np.zeros(nlambs, dtype='float64')
    outvec = outvec + 1.0     # set initial value to 1 at each lambda

    for i, s in enumerate(namelist):
        # Find the amount let through by each coating or layer for each lambda
        coatvalue = sl[i](lambdas) - delta
        # Multiply by the output vector
        outvec = outvec * coatvalue

    return outvec

        
