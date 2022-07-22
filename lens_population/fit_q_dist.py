import numpy as np
from scipy.optimize import leastsq, least_squares, minimize
from scipy.stats import beta
from astropy.io import fits as pyfits
from sl_cosmology import Dang
from scipy.interpolate import splrep, splev
import pylab
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
rc('text', usetex=True)


# fits a Beta distribution to the q distribution of the sample

samp_table = pyfits.open('sdss_zbin_latetype_red.fits')[1].data

cut = (samp_table['logmstar_dev'] > 11.) & (samp_table['logmstar_dev'] < 12.2) & (samp_table['dev_ba'] > 0.) & (samp_table['dev_ba'] < 1.)

"""
# optional: plots q for three mass bins, to see if there's a trend
# with stellar mass. There isn't.
qbins = np.linspace(0., 1., 11)

mbins = np.linspace(11., 12., 4)

for i in range(3):
    mbin = (samp_table['logmstar_dev'] > mbins[i]) & (samp_table['logmstar_dev'] < mbins[i+1])
    ngal = mbin.sum()
    weights = np.ones(ngal)/ngal
    pylab.hist(samp_table['dev_ba'][mbin], weights=weights, bins=qbins, histtype='step')
pylab.show()
"""

q_samp = samp_table['dev_ba'][cut]

def logpfunc(p):
    a, b = p
    return -beta.logpdf(q_samp, a, b).sum()

p0 = [2., 2.]
print(logpfunc(p0))

s = minimize(logpfunc, p0)
print(s)
a, b = s.x

qbins = np.linspace(0., 1., 11)

weights = np.ones(len(q_samp))/len(q_samp)/(qbins[1] - qbins[0])

fsize = 18

fig, ax = pylab.subplots()

pylab.subplots_adjust(left=0.14, right=1., bottom=0.14, top=0.98)

ax.hist(q_samp, bins=qbins, weights=weights, histtype='stepfilled', label='Observed')
q_arr = np.linspace(0., 1., 101)

ax.plot(q_arr, beta.pdf(q_arr, a, b), label='Model', linewidth=2)

ax.legend(loc='upper left', fontsize=fsize)

ax.tick_params(axis='both', which='both', top=True, right=True, labelsize=fsize, direction='in')#, width=1)

ax.set_xlabel('Axis ratio $q$', fontsize=fsize)
ax.set_ylabel('P$(q)$', fontsize=fsize)

ax.xaxis.set_major_locator(MultipleLocator(0.2))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))

ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))

pylab.savefig('../paper/q_dist.eps')
pylab.show()

