import numpy as np
from astropy.io import fits as pyfits
from skimage import measure
from scipy.special import gamma as gfunc
from scipy.optimize import brentq
from sl_profiles import sersic
from lensdet import detect_lens
import h5py


ftot = 200.
nre = 21
logre_grid = np.linspace(-1., 1., nre)

nsim = 10000

rsim = 1.5
source_area = np.pi * rsim**2

pix = 0.05
nser = 1.

nfrat = 5
lfrat_grid = np.linspace(1., 3., nfrat)
# frat is the ratio between the source intrinsic flux and the sky rms over an area equal to the square of the Einstein radius.
# In my sims, the Einstein radius is 1 arcsec. Pixels have a 0.05" size.
# Then, the sky rms over a 1 square arcsec is 20 times that over a single pixel
# All simulations have a source intrinsic flux of 200. Then an frat=10 
# corresponds to a sky rms per pixel of 1.

nsigma = 2.

minmag = 3.

x0 = 39.5
y0 = 39.5

X, Y = np.meshgrid(np.arange(80), np.arange(80))

cs_grid = np.zeros((nfrat, nre))

for l in range(nfrat):
    lfrat = lfrat_grid[l]
    sky_rms = 10.**(1. - lfrat) 
    sb_min = nsigma * sky_rms

    for m in range(nre):
        print(l, m)
        logre = logre_grid[m]

        preamble = 'ftot200_logre%2.1f'%logre
        nlens = 0
        islens_sim = np.zeros(nsim, dtype=bool)
        nimages_sim = np.zeros(nsim, dtype=int)

        img_file = h5py.File('mockdir/logre%2.1f_images.hdf5'%logre, 'r')

        for i in range(nsim):
        
            img = img_file['lens_%04d'%i][()]

            islens, nimg_std, nimg_max, std_footprint, best_footprint = detect_lens(img, sky_rms)
            islens_sim[i] = islens

        cs_grid[l, m] = islens_sim.sum()/float(nsim) * source_area

output = h5py.File('crosssect.hdf5', 'w')

output.create_dataset('lfrat_grid', data=lfrat_grid)
output.create_dataset('logre_grid', data=logre_grid)
output.create_dataset('cs_grid', data=cs_grid)

