import numpy as np
import emcee
from astropy.io import fits as pyfits
from sl_cosmology import Dang
from scipy.interpolate import splrep, splev
from scipy.stats import truncnorm
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

mu_par = {'name': 'mu', 'lower': -3., 'upper': 3., 'guess': 1., 'spread': 0.1}
beta_par = {'name': 'beta', 'lower': -3., 'upper': 3., 'guess': 1., 'spread': 0.1}
sig_par = {'name': 'sig', 'lower': 0., 'upper': 1., 'guess': 0.2, 'spread': 0.1}

pars = [mu_par, beta_par, sig_par]
npars = len(pars)

bounds = []
for par in pars:
    bounds.append((par['lower'], par['upper']))

def logprior(p):
    for i in range(npars):
        if p[i] < bounds[i][0] or p[i] > bounds[i][1]:
            return -1e300
    return 0.

def logpfunc(p):

    lprior = logprior(p)
    if lprior < 0.:
        return -1e300

    mu, beta, sig = p

    lreff_model = mu + beta * (lmstar - lmstar_piv)
    
    sumlogp = (-0.5*(lreff_model - lreff)**2/sig**2 - np.log(sig)).sum()

    return sumlogp

nwalkers = 100

sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc)

start = []
for i in range(nwalkers):
    tmp = np.zeros(npars)
    for j in range(npars):
        a, b = (bounds[j][0] - pars[j]['guess'])/pars[j]['spread'], (bounds[j][1] - pars[j]['guess'])/pars[j]['spread']
        p0 = truncnorm.rvs(a, b, size=1)*pars[j]['spread'] + pars[j]['guess']
        tmp[j] = p0
    start.append(tmp)

print("Sampling")

sampler.run_mcmc(start, 1000)

blobchain = sampler.blobs

for n in range(npars):
    print('%s: %4.3f'%(pars[n]['name'], np.median(sampler.chain[:, 300:, n])))

