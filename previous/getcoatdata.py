# import math
import numpy as np
import os
from scipy.interpolate import interp1d # spline interpolator
from scipy.interpolate import splrep # spline interpolator
from scipy.interpolate import InterpolatedUnivariateSpline # spline interpolator

def getcoatingdata(namelist, fileroot='emissdatafiles', irdust=0.8):
    """
    read in coating data. 
    INPUT: 
    namelist: list of coatings. 
    fileroot: directory containing coating data files
    irdust: transmission coefficient (dafault = 0.8)

    RETURNS 
    wavelist: wavelength vectors in nm 
    coatlist: coating data, in the order given in <namelist>
    splinelist: cubic spline fitting function
    """

    # Find number of names in list
    nnames = len(namelist)

    wavelist = []
    coatlist = []
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

        # spline = interp1d(wavevec, coatvec, kind='cubic')
        spline = InterpolatedUnivariateSpline(wavevec, coatvec, k=3)
        # k = 3rd degree smoothing spline

        wavelist.append(wavevec)
        coatlist.append(coatvec)
        splinelist.append(spline)

    return wavelist, coatlist, splinelist


#namelist = ['MEMS_window', 'PICNIC_QE', 'IR_transmission', 'Kprime', 
#            'Hband', 'Jband', 'ALPAO', 'Aluminum', 'AluminumPlusDust', 
#            'ARIR', 'CaF', 'FSG98', 'H2RG_QE', 'LGSWFS_Dich', 'NaHG', 
#            'NaR_Splitter', 'TT_Dichroic_Refl','X1_Silver_Extrap']
#            'lick_ir_background', 'MEMS_window', 'TT_Dichroic_Trans']
#            'deimos_qe_data', 'DEIMOS_Camera', 'Grating_900_550']

#getcoatingdata(namelist)
