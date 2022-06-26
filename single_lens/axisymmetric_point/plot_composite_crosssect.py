import numpy as np
import pylab
import h5py
from scipy.optimize import brentq
from plotters import probcontour
from labellines import labelLine, labelLines
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
rc('text', usetex=True)


fsize = 18

logscale = True

if logscale:
    leftm = 0.14
else:
    leftm = 0.11

cs_file = h5py.File('composite_rein1.0reff_grid.hdf5', 'r')

dms_grid = cs_file['dms_grid'][()]
ndms = len(dms_grid)

gammadm_grid = cs_file['gammadm_grid'][()]
ngamma = len(gammadm_grid)

fdm_grid = cs_file['fdm_grid'][()]
nfdm = len(fdm_grid)

cs_gammadm = cs_file['cs_vs_gammadm'][()]
cs_fdm = cs_file['cs_vs_fdm'][()]

fig, ax = pylab.subplots(1, 2, figsize=(8, 4))
pylab.subplots_adjust(left=leftm, right=0.99, bottom=0.13, top=0.99, wspace=0.)

colseq = pylab.rcParams['axes.prop_cycle'].by_key()['color']

for i in range(ndms):
    ax[0].plot(gammadm_grid, cs_gammadm[:, i], label='$\Delta m_s=%2.1f$'%dms_grid[i], linewidth=2, color=colseq[i])
    ax[1].plot(fdm_grid, cs_fdm[:, i], label='$\Delta m_s=%2.1f$'%dms_grid[i], linewidth=2, color=colseq[i])

ax[0].xaxis.set_major_locator(MultipleLocator(0.5))
ax[0].xaxis.set_minor_locator(MultipleLocator(0.1))

ax[0].tick_params(axis='both', which='both', top=True, right=True, labelsize=fsize, direction='in')

ax[0].set_xlabel('$\gamma_{\mathrm{DM}}$', fontsize=fsize)
ax[0].set_ylabel('$\sigma_{\mathrm{SL}}$', fontsize=fsize)

ax[1].xaxis.set_major_locator(MultipleLocator(0.2))
ax[1].xaxis.set_minor_locator(MultipleLocator(0.05))

ax[1].tick_params(axis='both', which='both', top=True, right=True, labelsize=fsize, direction='in', labelleft=False)

ax[1].set_xlabel('$f_{\mathrm{DM}}$', fontsize=fsize)

"""
lines = pylab.gca().get_lines()
gamma_vals = 2.3 - 0.05*np.arange(6)
labelLines(lines, xvals=gamma_vals, fontsize=fsize, backgroundcolor='white')

ax.plot(gamma_grid[gamma_grid<=2.], np.pi*cs_file['beta_caust_grid'][()][gamma_grid<=2.]**2, linestyle='--', color='k', label='Multiple images')
ax.plot([2., 2.], [np.pi*cs_file['beta_caust_grid'][()][gamma_grid<=2.][-1]**2, 100.], linestyle='--', color='k')

ax.xaxis.set_major_locator(MultipleLocator(0.2))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))

if logscale:
    ax.set_yscale('log')
    ax.set_ylim(0., 13.)
else:
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

ax.tick_params(axis='both', which='both', top=True, labelsize=fsize, direction='in')

ax.set_xlabel('$\gamma$', fontsize=fsize)
ax.set_ylabel('$\sigma_{\mathrm{SL}}$ (arcsec$^2$)', fontsize=fsize)

#ax.set_xlim(gamma_grid[0], gamma_grid[-1])
ax.set_xlim(1.5, 2.5)

#ax.legend(loc='upper left', fontsize=fsize)

"""
pylab.show()


