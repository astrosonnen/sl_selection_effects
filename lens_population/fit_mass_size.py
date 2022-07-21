import numpy as np
from scipy.optimize import leastsq, least_squares
from astropy.io import fits as pyfits
from sl_cosmology import Dang
from scipy.interpolate import splrep, splev


samp_table = pyfits.open('sdss_zbin_latetype_red.fits')[1].data

cut = (samp_table['logmstar_dev'] > 11.) & (samp_table['logmstar_dev'] < 12.2)

lmstar = samp_table['logmstar_dev'][cut]
lmstar_piv = 11.4

reff_arcsec = samp_table['dev_reff'][cut]

zmin = 1.9
zmax = 2.1
zgrid = np.linspace(zmin, zmax, 5)
dd_grid = 0.*zgrid
for i in range(5):
    dd_grid[i] = Dang(zgrid[i])

dd_spline = splrep(zgrid, dd_grid)

reff_kpc = reff_arcsec * np.deg2rad(1./3600.) * splev(samp_table['z'][cut], dd_spline) * 1000.

lreff = np.log10(reff_kpc)
lreff_arr = lreff.reshape((1, len(lreff)))
lmstar_arr = lmstar.reshape((1, len(lmstar)))

def fitfunc(p, x):
    return p[0] + p[1]*(x - lmstar_piv)

def errfunc(p, x, y):
    #return lreff - fitfunc(p)
    return y - fitfunc(p, x)

r = leastsq(errfunc, [1., 1.], (lmstar[:10], lreff[:10]))
print(r)
print(errfunc([1., 1.], lmstar[:10], lreff[:10]))


