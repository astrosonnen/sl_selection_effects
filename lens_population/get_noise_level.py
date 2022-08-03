import numpy as np
from astropy.io import fits as pyfits
from scipy.interpolate import splrep, splev


nrms = 41
rms_grid = np.linspace(0.001, 0.005, nrms)

sci = pyfits.open('stdsource_source.fits')[0].data

sn_grid = 0.*rms_grid

nsigma = 2.

for i in range(nrms):
    footprint = sci > nsigma * rms_grid[i]
    
    npix_here = footprint.sum()

    signal = sci[footprint].sum()
    noise = npix_here**0.5 * rms_grid[i]

    sn_grid[i] = signal/noise

rms_spline = splrep(np.flipud(sn_grid), np.flipud(rms_grid))

sky_rms = splev(10., rms_spline)
print(sky_rms)

"""

import pylab
pylab.plot(rms_grid, sn_grid)
pylab.axvline(sky_rms)
pylab.show()
"""


