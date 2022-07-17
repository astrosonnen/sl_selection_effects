import numpy as np
import pylab
import sys
from astropy.io import fits as pyfits
from skimage import measure
from matplotlib import rc
rc('text', usetex=True)


filename = sys.argv[1]
sky_rms = float(sys.argv[2])
thresh = 2.*sky_rms

x0 = 39.5
y0 = 39.5

fsize = 18

fig, ax = pylab.subplots(1, 1, figsize=(6, 6))
pylab.subplots_adjust(left=0., right=1.00, bottom=0., top=1., wspace=0.)

img = pyfits.open(filename)[0].data

footprint = img > thresh
labels = measure.label(footprint)
nreg = labels.max()
print('%d regions'%nreg)
for n in range(nreg):
    npix_here = (labels==n+1).sum()
    signal = img[labels==n+1].sum()
    noise = npix_here**0.5 * sky_rms
    img_sn = signal/noise
    if img_sn < 10.:
        img[labels==n+1] = 0.
    print(n, img_sn)

ax.contourf(img,[thresh,1e8], colors='b')
ax.set_aspect(1.)
ax.scatter(x0, y0, color='k')

pylab.show()

