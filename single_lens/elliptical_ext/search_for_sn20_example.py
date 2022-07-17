import numpy as np
from astropy.io import fits as pyfits
from skimage import measure
import h5py


ftot = 200.

nsim = 10000

logre = -0.8
lfrat = 1.
sky_rms = 10.**(1. - lfrat) 
nsigma = 2.

# frat is the ratio between the source intrinsic flux and the sky rms over an area equal to the square of the Einstein radius.
# In my sims, the Einstein radius is 1 arcsec. Pixels have a 0.05" size.
# Then, the sky rms over a 1 square arcsec is 20 times that over a single pixel
# All simulations have a source intrinsic flux of 200. Then an frat=10 
# corresponds to a sky rms per pixel of 1.

x0 = 39.5
y0 = 39.5

X, Y = np.meshgrid(np.arange(80), np.arange(80))

preamble = 'ftot200_logre%2.1f'%logre

for i in range(nsim):
    img = pyfits.open('mockdir/'+preamble+'/'+preamble+'_%04d_image.fits'%i)[0].data
    footprint = img > nsigma * sky_rms

    labels = measure.label(footprint)
    nreg = labels.max()
    nimg = 0
    sn_close_to_20 = False
    for n in range(nreg):
        npix_here = (labels==n+1).sum()
        signal = img[labels==n+1].sum()
        noise = npix_here**0.5 * sky_rms
        img_sn = signal/noise
        if img_sn < 20.:
            img[labels==n+1] = 0.
        else:
            nimg += 1
            if img_sn < 21.:
                sn_close_to_20 = True

    if nimg == 2 and sn_close_to_20:
        print(i)

