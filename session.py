%load_ext autoreload

%autoreload 2

from emissivity import *

namelist = ['MEMS_window', 'PICNIC_QE', 'IR_transmission', 'Kprime', 
            'Hband', 'Jband', 'ALPAO', 'Aluminum', 'AluminumPlusDust', 
            'ARIR', 'CaF', 'FSG98', 'H2RG_QE', 'LGSWFS_Dich', 'NaHG', 
            'NaR_Splitter', 'TT_Dichroic_Refl','X1_Silver_Extrap']

lambdas = np.arange(600, 2350, 50)

t1 = 10

out, theff = emissivity(lambdas, namelist, t1)

output = np.column_stack((lambdas.flatten(),out.flatten(),theff.flatten()))

np.savetxt('emissivity-py.out',output)

clf()

grid(True)

yscale('log')

pyplot.xlim(800.0,2300.0)

pyplot.ylim(1e-15,1e5)

plot(lambdas, out, '-')

pyplot.xlabel('Wavelength (nm)')

pyplot.ylabel(r'#/(s m$^2$ arcsec$^2$ nm)')

