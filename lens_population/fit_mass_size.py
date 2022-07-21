import numpy as np
from scipy.optimize import leastsq, least_squares, minimize
from astropy.io import fits as pyfits
from sl_cosmology import Dang
from scipy.interpolate import splrep, splev
import pylab


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

# I have no idea why I had to do this, but leastsq wouldn't work otherwise.
lmstar0 = np.zeros(len(lmstar))
lmstar0 += lmstar
lmstar = lmstar0.copy()

def fitfunc(p):
    return p[0] + p[1]*(lmstar - lmstar_piv)

def errfunc(p):
    return fitfunc(p) - lreff

p0 = [1., 1.]
s = leastsq(errfunc, p0)
print(s[0])

model_lreff = s[0][0] + s[0][1] * (lmstar - lmstar_piv)
scatter = (((model_lreff - lreff)**2).sum()/len(lreff))**0.5
print(scatter)

pylab.scatter(lmstar, lreff)
xlim = pylab.xlim()
xs = np.linspace(xlim[0], xlim[1])
pylab.plot(xs, s[0][0] + s[0][1]*(xs - lmstar_piv), linestyle='--', color='k')
pylab.show()

