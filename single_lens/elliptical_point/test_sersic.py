import numpy as np
import os
import glafic
import h5py


np.random.seed(0)

# primary parameters
omegaM = 0.3
omegaL = 0.7
weos = -1.
hubble = 0.7
prefix = 'sersic'
xmin = -3.
ymin = -3.
xmax = 3.
ymax = 3.
pix_ext = 0.2
pix_poi = 0.1
maxlev = 5

glafic.init(omegaM, omegaL, weos, hubble, prefix, xmin, ymin, xmax, ymax, pix_ext, pix_poi, maxlev, verb = 0)

# sets parameters of the Sersic profile
reff = 1.
nser = 4.
q = 1.
e = 1. - q
zd = 0.3
zs = 1.5
Mtot = 1e11 # this should be in units of M_Sun/h

glafic.startup_setnum(1, 0, 1)
glafic.set_lens(1, 'sers', zd, Mtot, 0.0, 0.0, e, 0.0, reff, nser)
glafic.set_point(1, 1.5, 0.5, 0.)

# model_init needs to be done again whenever model parameters are changed
glafic.model_init(verb = 0)

# writes critical curves and caustics
#glafic.writecrit(zs)

# calculates the Einstein radius
glafic.calcein2(zs, 0., 0., 1)

#os.system('mv out_crit.dat sersic_crit.dat')

glafic.quit()

