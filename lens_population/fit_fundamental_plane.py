import numpy as np
from scipy.optimize import leastsq, least_squares, minimize
from astropy.io import fits as pyfits
from sl_cosmology import Dang
from scipy.interpolate import splrep, splev
import pylab


samp_table = pyfits.open('sdss_zbin_latetype_red.fits')[1].data

cut = (samp_table['logmstar_dev'] > 11.) & (samp_table['logmstar_dev'] < 12.2) & (samp_table['veldisp'] > 0.) & (samp_table['veldisperr'] < samp_table['veldisp'])

lmstar = samp_table['logmstar_dev'][cut]
lmstar_piv = 11.4
lreff_piv = 1.2

veldisp = samp_table['veldisp'][cut]
veldisp_err = samp_table['veldisperr'][cut]

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

lsigma = np.zeros(len(lmstar))
lsigma_err = np.zeros(len(lmstar))
lsigma += np.log10(veldisp)
lsigma_err += 0.5*(np.log10(veldisp+veldisp_err) - np.log10(veldisp-veldisp_err))

def nlogpfunc(p):
    lsigma_model = p[0] + p[1]*(lmstar - lmstar_piv) + p[2]*(lreff - lreff_piv)
    logp_samp = -0.5*(lsigma - lsigma_model)**2/(lsigma_err**2 + p[3]**2) - 0.5*np.log((p[3]**2 + lsigma_err**2))
    return -logp_samp.sum()

p0 = [2.3, 0.2, -0.2, 0.1]
bounds = [(2., 3.), (-1., 1.), (-1., 1.), (0., 1.)]

print(nlogpfunc(p0))

res = minimize(nlogpfunc, p0, method='L-BFGS-B', bounds=bounds)
print(res.x)

