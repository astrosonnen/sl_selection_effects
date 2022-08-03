import numpy as np
from astropy.io import fits as pyfits
from skimage import measure
from scipy.special import gamma as gfunc
from scipy.signal import convolve2d
from scipy.optimize import brentq
from sl_profiles import sersic
from lensdet import detect_lens
import h5py


ntein = 5
ltein_grid = np.linspace(-0.9, -0.1, ntein)

nsim = 10000

nser = 1.

sky_rms = 1.
sb_min = 2. * sky_rms

cs_grid = np.zeros(ntein)

#for m in range(ntein):
for m in range(3, 4):

    print(m)
    ltein = ltein_grid[m]

    nlens = 0
    islens_sim = np.zeros(nsim, dtype=bool)
    nimages_sim = np.zeros(nsim, dtype=int)

    img_file = h5py.File('mockdir/ltein%2.1f_images.hdf5'%ltein, 'r')

    rmax = img_file.attrs['rmax']
    source_area = np.pi*rmax**2

    for i in range(nsim):

        img = img_file['lens_%04d_wseeing'%i][()]

        res = detect_lens(img, sky_rms, npix_min=1)
        islens_sim[i] = res[0]
        print(i, res[0])

    cs_grid[m] = islens_sim.sum()/float(nsim) * source_area

output = h5py.File('largesource_crosssect.hdf5', 'w')

output.create_dataset('ltein_grid', data=ltein_grid)
output.create_dataset('cs_grid', data=cs_grid)

