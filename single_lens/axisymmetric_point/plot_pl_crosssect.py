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

cs_file = h5py.File('pl_point_crosssect_grid.hdf5', 'r')

dms_grid = cs_file['dms_grid'][()]
nms = len(dms_grid)

nind = 6
dms_indices = 20*np.arange(nind)

gamma_grid = cs_file['gamma_grid'][()]
ngamma = len(gamma_grid)

cs_grid = cs_file['crosssect_grid'][()]

fig = pylab.figure()
pylab.subplots_adjust(left=leftm, right=0.99, bottom=0.13, top=0.99, wspace=0., hspace=0.)

ax = fig.add_subplot(1, 1, 1)

for i in range(nind):
    ax.plot(gamma_grid, cs_grid[:, dms_indices[i]], label='$\Delta m_s=%2.1f$'%dms_grid[dms_indices[i]], linewidth=2)

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

pylab.savefig('../../paper/axisymm_pl_crosssect.eps')
pylab.show()


