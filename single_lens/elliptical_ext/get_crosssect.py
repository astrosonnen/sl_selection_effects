import numpy as np
from astropy.io import fits as pyfits
from skimage import measure
import h5py


ftot_list = [200.]
re_list = [0.2]

nsim = 1000

rsim = 1.
source_area = np.pi * rsim**2

sky_rms = 1.
nsigma = 2.

x0 = 39.5
y0 = 39.5

X, Y = np.meshgrid(np.arange(80), np.arange(80))

for ftot in ftot_list:
    for re in re_list:
        preamble = 'ftot%d_re%2.1f'%(ftot, re)
        nlens = 0
        islens_sim = np.zeros(nsim, dtype=bool)
        nimages_sim = np.zeros(nsim, dtype=int)
        aperture_sim = np.zeros(nsim)

        for i in range(nsim):
            img = pyfits.open('mockdir/'+preamble+'_%04d_image.fits'%i)[0].data
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

            print('%d, %d pixels, %d images. Subtended angle: %d degrees'%(i, footprint.sum(), nimg, max_aperture))

            if nimg > 1:
                nlens += 1
                islens_sim[i] = True
                nimages_sim[i] = nimg
            elif nimg == 1:
                nimages_sim[i] = 1
                if max_aperture > 90.:
                    islens_sim[i] = True

output = h5py.File('%s_finding_results.hdf5'%preamble, 'w')

output.create_dataset('islens', data=islens_sim)
output.create_dataset('nimg', data=nimages_sim)
output.create_dataset('aperture', data=aperture_sim)

