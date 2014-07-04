class Observation:

    def __init__(self):
        self.seeing = 1.0
        self.mphase = 0.0
        self.airmass = 1.0
        self.exptime = 1.0
        self.mstar = 17.0
        self.mtype = 2.0 # 1.0 - Vega mag 2.0 AB mag - convention in code this is cloned from
        self.filter = '' # the filter for normalizing the spectrum
                         # a blank value means that the spectrum will be assumed to be of Vega or
                         # an flat f_nu spectrum
        self.spectrum = '' # a blank spectrum means a flat fnu spectrum
        self.redshift = 0.0 # the redshift of the spectrum
        
    def __repr__(self):
         return '<Observation  %f %f %f %f %f %f %s %s %f>' % (self.mstar, self.airmass, self.seeing, self.mphase, self.exptime, self.mtype, self.filter, self.spectrum, self.redshift)
