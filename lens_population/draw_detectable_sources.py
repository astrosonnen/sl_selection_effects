import numpy as np
from scipy.stats import poisson
from simpars import *
from sc_profiles import sersic
import h5py
import os
from scipy.signal import convolve2d
from astropy.io import fits as pyfits
from skimage import measure
import sys


np.random.seed(20)

ndraw = 100000
sourcecat = pyfits.open('/Users/alessandro/catalogs/skills_sourceonly_zcut.fits')[1].data

psf = pyfits.open('psf.fits')[0].data
pix_arcsec = 0.1
sb_min = nsigma_pixdet * sky_rms

nsource_tot = len(sourcecat)

sourceind = np.arange(nsource_tot)
# shuffles source catalog (it's originally ranked by redshift)
np.random.shuffle(sourceind)

sreff_draw = sourcecat['Re_arcsec_CM'][:ndraw]
nser_draw = sourcecat['sersic_n_CM'][:ndraw]
sq_draw = sourcecat['axis_ratio_CM'][:ndraw]
smag_draw = sourcecat['g_SDSS_apparent_corr'][:ndraw]
zs_draw = sourcecat['zobs'][:ndraw]
pa_draw = sourcecat['PA_random'][:ndraw]

f = open('preamble.input', 'r')
prelines = f.readlines()
f.close()

outlines = prelines.copy()

for i in range(ndraw):
    ind = sourceind[i]
    outlines.append('reset_par simdir/source_%07d\n'%ind)
    outlines.append('reset_extend 1 1 %f\n'%zs_draw[i])

    ftot = 10.**(-2./5.*(smag_draw[i] - zeropoint))
    I0 = ftot/(2.*np.pi*(sreff/pix_arcsec)**2*nser/sersic.b(nser)**(2*nser)*gfunc(2.*nser))

    outlines.append('reset_extend 1 2 %f\n'%I0)
    outlines.append('reset_extend 1 5 %f\n'%(1. - sq_draw[i]))
    outlines.append('reset_extend 1 6 %f\n'%pa_draw[i])
    outlines.append('reset_extend 1 7 %f\n'%sreff_draw[i])
    outlines.append('reset_extend 1 8 %f\n'%nser_draw[i])

    outlines.append('writeimage_ori\n')
    outlines.append('\n')

outlines.append('quit\n')

f = open('unlensed.input', 'w')
f.writelines(outlines)
f.close()

os.system('glafic unlensed.input')

detected = np.zeros(ndraw, dtype=bool)

for i in range(ndraw):
    ind = sourceind[i]

    img = pyfits.open('simdir/source_%07d_source.fits'%ind)[0].data
    img_wseeing = convolve2d(img, psf, mode='same')

    footprint = img_wseeing > sb_min

    labels = measure.label(footprint)
    nreg = labels.max()
    npix_tmp = (labels==1).sum()
    signal = img[labels==1].sum()
    noise = npix_tmp**0.5 * sky_rms
    img_sn = signal/noise
    if img_sn >= 10. and npix_tmp >= npix_min:
        detected[i] = True

ndet = detected.sum()
print('%d detected sources'%ndet)

output_file = h5py.File('detectable_sources.hdf5', 'w')
output_file.create_dataset('sreff', data=sreff_draw[detected])
output_file.create_dataset('smag', data=smag_draw[detected])
output_file.create_dataset('zs', data=zs_draw[detected])
output_file.create_dataset('sq', data=sq_draw[detected])
output_file.create_dataset('spa', data=spa_draw[detected])
output_file.create_dataset('nser', data=nser_draw[detected])
output_file.create_dataset('index', data=sourceind[detected])

