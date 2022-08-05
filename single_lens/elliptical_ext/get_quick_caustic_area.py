import numpy as np
from sl_profiles import sersic
from scipy.special import gamma as gfunc
from scipy.interpolate import splrep, splev
import os
import glafic
import h5py
from scipy.signal import convolve2d
from astropy.io import fits as pyfits
from wl_profiles import gnfw, deVaucouleurs as deV
from wl_cosmology import Mpc, c, G, M_Sun
import wl_cosmology
from scipy.optimize import brentq


nrein = 11
ltein_grid = np.linspace(-1., 0., nrein)

zs_ref = 1.5
pix = 0.05

# primary parameters
omegaM = 0.3
omegaL = 0.7
weos = -1.
hubble = 0.7
prefix = 'tmp'
xmin = -2.
ymin = -2.
xmax = 2.
ymax = 2.
pix_ext = pix
pix_poi = 0.1
maxlev = 5

glafic.init(omegaM, omegaL, weos, hubble, prefix, xmin, ymin, xmax, ymax, pix_ext, pix_poi, maxlev, verb = 0)
glafic.set_secondary('flag_hodensity 2')
glafic.set_secondary('nfw_users 1')
glafic.set_secondary('hodensity 200')

# WARNING: setting ellipticity to zero, for the purpose of computing
# the Einstein radius. Need to bring it back up to 0.3 later.
glafic.startup_setnum(2, 1, 0)
glafic.set_lens(1, 'gnfw', 0.3, 2.021e12, 0.0, 0.0, 0., 90.0, 10., 1.5)
glafic.set_lens(2, 'sers', 0.3, 1.087e11, 0.0, 0.0, 0., 90.0, 1., 4.)
glafic.set_extend(1, 'sersic', 1.5, 1., 0.3, 0., 0., 0., 0.1, 1.)

glafic.model_init(verb=0)


# calculates the Einstein radius on a grid of source redshifts
nz = 21
zs_grid = np.linspace(0.32, 1.8, nz)
tein_grid = 0.*zs_grid

for i in range(nz):
    tein_grid[i] = glafic.calcein2(zs_grid[i], 0., 0.)

zs_spline = splrep(tein_grid, zs_grid)

output_file = h5py.File('quick_caustic_area.hdf5', 'w')

cs_grid = 0.*ltein_grid

for i in range(nrein):

    tein = 10.**ltein_grid[i]

    zs = splev(tein, zs_spline)

    glafic.writecrit(zs)

    f = open('tmp_crit.dat', 'r')
    table = np.loadtxt(f)
    f.close()

    xs1 = table[:, 2]
    ys1 = table[:, 3]
    xs2 = table[:, 6]
    ys2 = table[:, 7]

    a = max(abs(xs1).max(), abs(xs2).max())
    b = max(abs(ys1).max(), abs(ys2).max())

    area = np.pi*a*b
    cs_grid[i] = area

output_file.create_dataset('ltein_grid', data=ltein_grid)
output_file.create_dataset('cs_grid', data=cs_grid)
