import numpy as np
import pylab
import h5py
from scipy.interpolate import splrep, splev
from plotters import probcontour
from labellines import labelLine, labelLines
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
rc('text', usetex=True)


fsize = 18

ylim = (0.1, 1.5)

logscale = True

if logscale:
    leftm = 0.14
else:
    leftm = 0.08

cs_file = h5py.File('crosssect.hdf5', 'r')

lfrat_grid = cs_file['lfrat_grid'][()]
nfrat = len(lfrat_grid)

logre_grid = cs_file['logre_grid'][()]
nre = len(logre_grid)

cs_grid = cs_file['cs_grid'][()]

fig, ax = pylab.subplots(1, 1)#, figsize=(8, 4))
pylab.subplots_adjust(left=leftm, right=0.98, bottom=0.14, top=0.97, wspace=0.)

colseq = pylab.rcParams['axes.prop_cycle'].by_key()['color']

for i in range(nfrat):
    ax.loglog(10.**logre_grid, cs_grid[i, :]/np.pi, linewidth=2, color=colseq[i])
    cs_spline = splrep(10.**logre_grid, cs_grid[i, :]/np.pi, k=1)

nfr2 = nfrat
cs_fixedfr2 = np.zeros((nfr2, nfrat))
re_fixedfr2 = np.zeros((nfr2, nfrat))

for i in range(nfr2):
    for j in range(nfrat):
        if i+j>1:
            re_fixedfr2[i, j] = 10.**logre_grid[i+j-2]
            cs_fixedfr2[i, j] = cs_grid[j, i+j-2]

for i in range(nfr2):
    if i<2:
        istart = 2-i
        iend = 4
    else:
        istart = 0
        iend = 4
    ax.loglog(re_fixedfr2[i, istart:iend], cs_fixedfr2[i, istart:iend]/np.pi, color='grey', linewidth=0.5)

#lines = labelLines(ax.get_lines()[2:], xvals = [0.65, 0.9, 0.58], fontsize=fsize, backgroundcolor='white')

#ax.xaxis.set_major_locator(MultipleLocator(0.5))
#ax.xaxis.set_minor_locator(MultipleLocator(0.1))

ax.tick_params(axis='both', which='both', top=True, right=True, labelsize=fsize, direction='in')

ax.set_xlabel('$\\theta_{\mathrm{e,s}}/\\theta_{\mathrm{Ein}}$', fontsize=fsize)
ax.set_ylabel('$\sigma_{\mathrm{SL}}/(\pi\\theta_{\mathrm{Ein}}^2)$', fontsize=fsize)

ax.set_ylim(ylim[0], ylim[1])
#ax.set_xlim(0.5, 1.)

#pylab.savefig('../../paper/ell_ext_cs.eps')
pylab.show()


