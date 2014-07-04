import telescope

class instrument:

    def __init__(self):
        self.name = ""
        self.sheight = 0.0 # Slit height
        self.swidth = 0.0 # Slit width
        self.dark = 0.0 # dark current
        self.readnoise = 0.0 # read noise
        self.bind = 1 # dispersion binning
        self.bins = 1 # spatial binning
        self.scale_para = 1.0 # pixel scale
        self.scale_perp = 1.0 # pixel scale
        self.mag_para = 1.0 # 
        self.mag_perp = 1.0 #
        self.R = 1.0 # lambda/delta-lambda
        self.pixel_size= 1.
        self.mlambda= 0.0
        self.dichroic= ''
        self.grating= ''
        self.cwave= 0. # Central wavelength
        self.wvmnx= [] #
        self.dely= 0.0 


def kast(tel,slit=1.5,grism="G2",grating="600/7500",):
    blue = instrument()
    red = instrument()

    blue.mag_perp = 20.9
    blue.mag_para = 20.9
    red.mag_perp = 20.9
    red.mag_para = 20.9
    blue.pixel_size = 15.0
    red.pixel_size = 27.0

    blue.scale_perp = tel.plate_scale* blue.mag_perp*blue.pixel_size/1000
    red.scale_perp = tel.plate_scale* red.mag_perp*red.pixel_size/1000
    blue.scale_para = tel.plate_scale* blue.mag_para*blue.pixel_size/1000
    red.scale_para = tel.plate_scale* red.mag_para*red.pixel_size/1000

    blue.grating = grism
    if grism == "G1":
        blue.R = 2344
    elif grism == "G2":
        blue.R = 4254
    elif grism == "G3":
        blue.R = 5492
    red.grating = grating
    if grating == "600/7500":
        red.R = 3164


    blue.readno = 3.7
    blue.dark = 0.001

    red.readno = 12.5
    red.dark = 0.001

    blue.bind = 1
    blue.bins = 1
    blue.wvmnx = [3000., 6000.]

    red.bind = 1
    red.bins = 1
    red.wvmnx = [5000., 10000.]

    red.swidth = slit
    blue.swidth = slit
    red.sheight = 120.0
    blue.sheight = 120.0

    return [blue,red]

def deimos(tel,slit=1.0,grating="600"):
    dei = instrument()
    dei.MAG_PERP = 8.03  # Modified to give observed resolution 0.75" maps to 4.5 pixels
    dei.MAG_PARA = 8.03  # Same
    dei.PIXEL_SIZE = 15.0 # in microns

    dei.SCALE_PERP = tel.PLATE_SCALE*dei.MAG_PERP*(dei.PIXEL_SIZE/1000.) # Arcsec
    dei.SCALE_PARA = tel.PLATE_SCALE*dei.MAG_PARA*(dei.PIXEL_SIZE/1000.) # Arcsec
          
    dei.grating = grating
    if grating == "600":
        dei.R     = 11538.5    # 1 pixel (native) dispersion
    elif grating == "1200":
        dei.R    = 22727.3    # 1 pixel (native) dispersion
        

    ## Detector
    dei.readno = 2.6
    dei.dark = 4.  # electrons/pix/hr
    dei.bind = 1
    dei.bins = 1

    ## Wavelength range
    dei.wvmnx = [4000., 10000.]

    ## Slit
    dei.swidth = slit 
    dei.sheight = 10.0  # arcsec

    return dei
