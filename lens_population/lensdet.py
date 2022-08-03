import numpy as np
from simpars import *
from skimage import measure


def detect_lens(img):

    ny, nx = img.shape
    x0 = nx/2. - 0.5
    y0 = ny/2. - 0.5

    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))

    footprint = img > nsigma_pixdet * sky_rms

    img_detected = img.copy()

    labels = measure.label(footprint)
    nreg = labels.max()
    nimg = 0
    for n in range(nreg):
        npix_here = (labels==n+1).sum()
        signal = img[labels==n+1].sum()
        noise = npix_here**0.5 * sky_rms
        img_sn = signal/noise
        if img_sn < 10. or npix_here < npix_min:
            img_detected[labels==n+1] = 0.
        else:
            nimg += 1

    # checks subtended angle

    new_footprint = img_detected > nsigma_pixdet * sky_rms

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

    # tries to maximise the number of detected images, by varying the threshold
    nimg_max = nimg

    sb_max = img.max()
    sb_min = nsigma_pixdet * sky_rms

    sorted_img = img.flatten().copy()
    sorted_img.sort()
    above_threshold = sorted_img[sorted_img > sb_min]
    nup = len(above_threshold)

    i = 0
    nimg_here = nimg
    sb_maxlim = sb_min

    while nimg_here >= nimg_max and i < nup:
        sb_lim = above_threshold[i]

        img_detected = img.copy()
        footprint_here = img > sb_lim

        labels = measure.label(footprint_here)
        nreg = labels.max()
        nimg_here = 0
        for n in range(nreg):
            npix_here = (labels==n+1).sum()
            signal = img[labels==n+1].sum()
            noise = npix_here**0.5 * sky_rms
            img_sn = signal/noise
            if img_sn >= 10. and npix_here >= npix_min:
                nimg_here += 1

        if nimg_here > nimg_max:
            nimg_max = nimg_here
            sb_maxlim = sb_lim

        i += 1

    if nimg_max > nimg:
        nmax_footprint = img_detected > sb_maxlim
    else:
        nmax_footprint = new_footprint

    return islens, nimg, nimg_max, new_footprint, nmax_footprint

