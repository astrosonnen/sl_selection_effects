import numpy as np
import pylab
from astropy.io import fits as pyfits
import h5py
from simpars import *
from scipy.optimize import leastsq
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib import rc
rc('text', usetex=True)


fsize = 20

sims = ['fiducial_1000sqdeg', 'highscatter_1000sqdeg']
labels = ['Fiducial', 'High scatter']
nsims = len(sims)

colseq = pylab.rcParams['axes.prop_cycle'].by_key()['color']

fig, ax = pylab.subplots(4, 1, figsize=(6, 12))

pylab.subplots_adjust(left=0.2, right=1.00, bottom=0.05, top=1., wspace=0., hspace=0.)

ntein = 21
tein_arr = np.linspace(0., 2., ntein)

for n in range(nsims):

    galpop = h5py.File('%s_galaxies.hdf5'%sims[n], 'r')

    laspsmed_gal = np.median(galpop['lasps'][()])
    laspsmed_arr = np.zeros(ntein)
    laspserr_arr = np.zeros(ntein)

    lm200med_gal = mu_h
    lm200med_arr = np.zeros(ntein)
    lm200err_arr = np.zeros(ntein)

    for i in range(ntein):
        lenscut = galpop['islens'][()] & (galpop['tein_zs'][()] > tein_arr[i])
        nlens = lenscut.sum()
        print('Theta_Ein > %2.1f. %d lenses'%(tein_arr[i], lenscut.sum()))
        laspsmed_arr[i] = np.median(galpop['lasps'][lenscut])
        laspserr_arr[i] = np.std(galpop['lasps'][lenscut])/float(nlens)**0.5

        # fits for the stellar-halo mass relation

        lm200_here = galpop['lm200'][lenscut]
        lmstar_here = galpop['lmstar'][lenscut]

        def fitfunc(p):
            return p[0] + p[1] * (lmstar_here - lmstar_piv)

        def errfunc(p):
            return fitfunc(p) - lm200_here

        pfit = leastsq(errfunc, [13., 1.])

        mu_h_here = pfit[0][0]
        lm200_scat = np.std(fitfunc(pfit[0]) - lm200_here)

        lm200med_arr[i] = mu_h_here
        lm200err_arr[i] = lm200_scat/float(nlens)**0.5

    ax[0].errorbar(tein_arr, laspsmed_arr, yerr=laspserr_arr, color=colseq[n], label=labels[n])
    ax[1].errorbar(tein_arr, lm200med_arr, yerr=lm200err_arr, color=colseq[n])
    #ax[2].errorbar(tein_arr, lreffmed_arr, yerr=lrefferr_arr, color=colseq[n])
    #ax[3].errorbar(tein_arr, qmed_arr, yerr=qerr_arr, color=colseq[n])

ax[0].set_ylabel('Median $\log{\\alpha_{\mathrm{SPS}}}$', fontsize=fsize)
ax[0].axhline(laspsmed_gal, color='k', linestyle='--', label='Parent population')

ax[0].yaxis.set_major_locator(MultipleLocator(0.02))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.005))

ax[0].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[0].legend(loc='lower left', fontsize=fsize)

ax[1].axhline(mu_h, color='k', linestyle='--')
ax[1].set_ylabel('$\mu_{\mathrm{h}}$', fontsize=fsize)

ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.02))

#ax[2].axhline(lreffmed_gal, color='k', linestyle='--')
#ax[2].set_ylabel('$\mu_{\mathrm{R}}$', fontsize=fsize)

ax[2].yaxis.set_major_locator(MultipleLocator(0.05))
ax[2].yaxis.set_minor_locator(MultipleLocator(0.01))

#ax[3].axhline(qmed_gal, color='k', linestyle='--')
#ax[3].set_ylabel('Median $q$', fontsize=fsize)
ax[3].set_xlabel('Minimum $\\theta_{\mathrm{Ein}}$', fontsize=fsize)

ax[3].yaxis.set_major_locator(MultipleLocator(0.02))
ax[3].yaxis.set_minor_locator(MultipleLocator(0.005))

ax[1].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[2].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[3].tick_params(axis='both', which='both', direction='in', labelsize=fsize, right=True, top=True)

for j in range(4):
    ax[j].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[j].xaxis.set_minor_locator(MultipleLocator(0.1))


pylab.savefig('../paper/lens_mass_bias.eps')
pylab.show()


