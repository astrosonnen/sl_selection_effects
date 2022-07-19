import numpy as np
from astropy.io import fits as pyfits
from skimage import measure
import h5py
import pylab


np.random.seed(0)

ftot = 200.
nre = 11
logre_grid = np.linspace(-1., 0., nre)

nsim = 10000

# first re-generates source positions to match those of the sims
rsim = 1.5
source_area = np.pi * rsim**2

source_x = []
source_y = []

r = np.random.rand(nsim)**0.5 * rsim
phi = 2.*np.pi*np.random.rand(nsim)
source_x = r*np.cos(phi)
source_y = r*np.sin(phi)

lfrat = 2.
# frat is the ratio between the source intrinsic flux and the sky rms over an area equal to the square of the Einstein radius.
# In my sims, the Einstein radius is 1 arcsec. Pixels have a 0.05" size.
# Then, the sky rms over a 1 square arcsec is 20 times that over a single pixel
# All simulations have a source intrinsic flux of 200. Then an frat=10 
# corresponds to a sky rms per pixel of 1.

nsigma = 2.

x0 = 39.5
y0 = 39.5

X, Y = np.meshgrid(np.arange(80), np.arange(80))

sky_rms = 10.**(1. - lfrat) 

re_inds = [0]

for m in re_inds:
    logre = logre_grid[m]

    preamble = 'ftot200_logre%2.1f'%logre
    nlens = 0
    islens_sim = np.zeros(nsim, dtype=bool)
    nimages_sim = np.zeros(nsim, dtype=int)
    aperture_sim = np.zeros(nsim)

    for i in range(nsim):
        img = pyfits.open('mockdir/'+preamble+'/'+preamble+'_%04d_image.fits'%i)[0].data
        footprint = img > nsigma * sky_rms

        labels = measure.label(footprint)
        nreg = labels.max()
        nimg = 0
        for n in range(nreg):
            npix_here = (labels==n+1).sum()
            signal = img[labels==n+1].sum()
            noise = npix_here**0.5 * sky_rms
            img_sn = signal/noise
            if img_sn < 10.:
                img[labels==n+1] = 0.
            else:
                nimg += 1

        # checks subtended angle

        new_footprint = img > nsigma * sky_rms

        ypix = Y[new_footprint]
        xpix = X[new_footprint]
        rpix = ((xpix - x0)**2 + (ypix - y0)**2)**0.5
        cospix = (xpix - x0)/rpix
        sinpix = (ypix - y0)/rpix

        npix = len(xpix)

        max_aperture = 0.
        for j in range(npix):
            cosdiff = cospix[j]*cospix + sinpix[j]*sinpix
            aperture = 180.*np.arccos(cosdiff).max()/np.pi
            if aperture > max_aperture:
                max_aperture = aperture

        aperture_sim[i] = max_aperture

        #print('%d, %d pixels, %d images. Subtended angle: %d degrees'%(i, footprint.sum(), nimg, max_aperture))

        if nimg > 1:
            nlens += 1
            islens_sim[i] = True
            nimages_sim[i] = nimg
        elif nimg == 1:
            nimages_sim[i] = 1
            if max_aperture > 90.:
                islens_sim[i] = True

    print(m, islens_sim.sum())
    pylab.scatter(source_x[islens_sim], source_y[islens_sim])

pylab.show()

