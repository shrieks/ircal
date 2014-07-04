class Telescope(object):

    def __init__(self):
        self.name = '' # KeckI, KeckII, Lick-3m 
        self.area= 0.
        self.plate_scale= 0. 

    def __repr__(self):
         return '<Telescope "%s" %f %f>' % (self.name, self.area, self.plate_scale)
    
    def lick(self):
        self.area  = 63617. # cm^2 -- 3m telescope with 10% central obscuration
        self.plate_scale = 1.379 # "/mm  [Should confirm; secondary dependent]
        self.name = 'Lick-3m'
        return(self)

    def keck(self):
        self.area  = 723674. # cm^2 -- 10m telescope with 7.9% central obscuration
        self.plate_scale = 1.379 # "/mm
        return(self)

    def keckone(self):
        self = Telescope.keck()
        self.name = 'KeckI'
        return(self)

    def kecktwo(self):
        self = Telescope.keck()
        self.name = 'KeckII'
        return(self)
