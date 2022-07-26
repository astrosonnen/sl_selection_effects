import numpy as np
import os
import glafic
import h5py
from astropy.io import fits as pyfits


modelname = 'fiducial_100sqdeg'
sourcecat = pyfits.open('/Users/alessandro/catalogs/skills_sourceonly_zcut.fits')[1].data
pop = h5py.File('%s_lenses.hdf5'%modelname, 'r')

f = open('%s_sources.cat'%modelname, 'r')
sourcelines = f.readlines()
f.close()

np.random.seed(20)

# primary parameters
omegaM = 0.3
omegaL = 0.7
weos = -1.
hubble = 0.7
prefix = 'out'
xmin = -3.
ymin = -3.
xmax = 3.
ymax = 3.
pix_ext = 0.1
pix_poi = 0.1
maxlev = 5

glafic.init(omegaM, omegaL, weos, hubble, prefix, xmin, ymin, xmax, ymax, pix_ext, pix_poi, maxlev, verb = 0)
glafic.set_secondary('flag_halodensity 2')
glafic.set_secondary('nfw_users 1')
glafic.set_secondary('halodensity 200')

# sets parameters of the power-law lens
rein = 1.
gamma0 = 2.
q = 1.
e = 1. - q
zd = 0.3
zs = 1.5

glafic.startup_setnum(2, 1, 0)
glafic.set_lens(1, 'gnfw', zd, 1e13, 0.0, 0.0, e, 90.0, 10., 1.5)
glafic.set_lens(2, 'sers', zd, 1e11, 0.0, 0.0, e, 90.0, 1., 4.)
glafic.set_extend(1, 'sersic', 1.5, 0.5, 0.3, 0., 0., 0., 0.1, 1.)

#glafic.set_primary('xmax 5.')
#glafic.set_primary('xmin -5.')
glafic.set_primary(omegaM, omegaL, weos, hubble, prefix, xmin, ymin, xmax, ymax, pix_ext, pix_poi, maxlev, verb = 0)

# model_init needs to be done again whenever model parameters are changed
glafic.model_init(verb = 0)

a = glafic.writeimage()
print(np.array(a).shape)

pyfits.PrimaryHDU(a).writeto('test.fits', overwrite=True)


