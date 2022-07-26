import numpy as np
from scipy.stats import poisson
from poppars import *
import h5py
from astropy.io import fits as pyfits


np.random.seed(10)
# places sources in a circle enclosing the radial caustic

circ_caust_rat = 1.2 # ratio between circle radius and caustic radius

modelname = 'fiducial_100sqdeg'
pop = h5py.File('%s_lenses.hdf5'%modelname, 'r')

nsamp = pop.attrs['nsamp']

sourcecat = pyfits.open('/Users/alessandro/catalogs/skills_sourceonly_zcut.fits')[1].data

nsource_tot = len(sourcecat)

sourceind = np.arange(nsource_tot)
# shuffles source catalog (it's originally ranked by redshift)
np.random.shuffle(sourceind)

outlines = []
outlines.append('# lens_id rcirc nsources source_indices(list)\n')

sourcecount = 0
for i in range(nsamp):
    rcirc = circ_caust_rat*pop['tcaust'][i]/pop['q'][i]**0.5

    area = np.pi*rcirc**2

    lam = area * nbkg / 3600. # expectation value of the number of sources in the circle

    nsources = poisson.rvs(lam)

    line = '%d %f %d'%(i, rcirc, nsources)
    for n in range(nsources):
        line += ' %d'%sourceind[sourcecount]
        sourcecount += 1
    line += '\n'

    outlines.append(line)

print('Total number of sources: %d'%sourcecount)
f = open('%s_sources.cat'%modelname, 'w')
f.writelines(outlines)
f.close()

