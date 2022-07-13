import numpy as np
import pylab
from astropy.io import fits as pyfits
from matplotlib import rc
rc('text', usetex=True)


suffix = ['x0.0_y0.0', 'x0.4_y0.0', 'x0.122_y0.105', 'x0.6_y0.0']

level = 1.8

leftm = 0.15
fig, ax = pylab.subplots(1, 4, figsize=(12, 4))
pylab.subplots_adjust(left=leftm, right=1.00, bottom=0.15, top=0.99, wspace=0.)

for i in range(4):
    img = pyfits.open('composite_%s_image.fits'%suffix[i])[0].data

    ax[i].contourf(img,[level,1e8], colors='b')
    ax[i].set_aspect(1.)

pylab.show()


