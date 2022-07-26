import numpy as np
from simpars import *
from skimage import measure


def detect_lens(img):

    ny, nx = img.shape
    x0 = nx/2. - 0.5
    y0 = ny/2. - 0.5

    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))

    footprint = img > nsigma_pixdet * sky_rms

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

    new_footprint = img > nsigma_pixdet * sky_rms

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

    islens = False

    if nimg > 1:
        islens = True
    elif nimg == 1:
        if max_aperture > min_angle:
            islens = True

    return islens, nimg

