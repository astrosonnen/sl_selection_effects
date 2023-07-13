import numpy as np
import pylab
import h5py
from simpars import lmobs_piv, lreff_piv
from scipy.optimize import leastsq
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib import rc
rc('text', usetex=True)


fsize = 20

smag_min = 24.8
smag_max = 25.2

lpopnames = ['fiducial_1000sqdeg_lenses.hdf5', 'fiducial_1000sqdeg_minmulenses.hdf5']
labels = ['Fiducial', '$\mu_{\mathrm{tot}} > 3$']
lensdefs = ['islens', 'minmulens']
teindefs = ['tein_zs', 'minmutein_zs']
nsims = len(lpopnames)

colseq = pylab.rcParams['axes.prop_cycle'].by_key()['color']
colors = [colseq[0], colseq[5]]
linestyles = ['-', '--']

fig, ax = pylab.subplots(3, 1, figsize=(6, 10))

pylab.subplots_adjust(left=0.22, right=1.00, bottom=0.07, top=1., wspace=0., hspace=0.)

ntein = 21
tein_arr = np.linspace(0., 2., ntein)

detectable = h5py.File('detectable_sources.hdf5', 'r')

zmed_det = np.median(detectable['zs'][()])
smagmed_det = np.median(detectable['smag'][()])

smagbin_det = (detectable['smag'][()] > smag_min) & (detectable['smag'][()] < smag_max)

nsermed_det = np.median(detectable['nser'][smagbin_det])
sreffmed_det = np.median(detectable['sreff'][smagbin_det])

qmed_det = np.median(detectable['sq'][()])

galpop = h5py.File('fiducial_1000sqdeg_galaxies.hdf5', 'r')
laspsmed_gal = np.median(galpop['lasps'][()])

for n in range(nsims):

    laspsmed_arr = np.zeros(ntein)
    laspserr_arr = np.zeros(ntein)

    lmdm5med_arr = np.zeros(ntein)
    lmdm5err_arr = np.zeros(ntein)

    def mdm5fitfunc(p):
        return p[0] + p[1] * (galpop['lmobs'][()] - lmobs_piv) + p[2] * (galpop['lreff'][()] - lreff_piv)

    def mdm5errfunc(p):
        return mdm5fitfunc(p) - galpop['lmdm5'][()]

    pmdm5fit = leastsq(mdm5errfunc, [11., 0., 0.])

    lmdm5med_gal = pmdm5fit[0][0]

    for i in range(ntein):
        lenscut = galpop[lensdefs[n]][()] & (galpop[teindefs[n]][()] > tein_arr[i])

        nlens = lenscut.sum()
        print('Theta_Ein > %2.1f. %d lenses'%(tein_arr[i], lenscut.sum()))

    lenspop = h5py.File(lpopnames[n], 'r')

    zmed_arr = np.zeros(ntein)
    zerr_arr = np.zeros(ntein)

    smagmed_arr = np.zeros(ntein)
    smagerr_arr = np.zeros(ntein)

    sreffmed_arr = np.zeros(ntein)
    srefferr_arr = np.zeros(ntein)

    for i in range(ntein):
        teincut = (lenspop['tein_zs'][()] > tein_arr[i])
        nlens = teincut.sum()

        zmed_arr[i] = np.median(lenspop['zs'][teincut])
        zerr_arr[i] = np.std(lenspop['zs'][teincut])/float(nlens)**0.5

        smagmed_arr[i] = np.median(lenspop['smag'][teincut])
        smagerr_arr[i] = np.std(lenspop['smag'][teincut])/float(nlens)**0.5

        smagbin = teincut & (lenspop['smag'][()] > smag_min) & (lenspop['smag'][()] < smag_max)
        nbin = smagbin.sum()

        sreffmed_arr[i] = np.median(lenspop['sreff'][smagbin])
        srefferr_arr[i] = np.std(lenspop['sreff'][smagbin])/float(nbin)**0.5

    ax[0].errorbar(tein_arr, zmed_arr, yerr=zerr_arr, color=colors[n], label=labels[n])
    ax[1].errorbar(tein_arr, smagmed_arr, yerr=smagerr_arr, color=colors[n], label=labels[n])
    ax[2].errorbar(tein_arr, sreffmed_arr, yerr=srefferr_arr, color=colors[n], label=labels[n], linestyle=linestyles[n])

# NEED TO REMOVE INDENTATION
ax[0].axhline(zmed_det, color='k', linestyle='--', label='Detectable population')
ax[0].set_ylabel('Median $z_{\mathrm{s}}$', fontsize=fsize)

ax[0].yaxis.set_major_locator(MultipleLocator(0.1))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.02))

ax[1].axhline(smagmed_det, color='k', linestyle='--', label='Detectable pop.')
ax[1].set_ylabel('Median $m_{\mathrm{s}}$', fontsize=fsize)

ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.02))

ax[2].axhline(sreffmed_det, color='k', linestyle='--')
ax[2].set_ylabel("Median $\\theta_{\mathrm{s}}$\n at $m_{\mathrm{s}}=25\, ('')$", fontsize=fsize)
ax[2].set_xlabel('Minimum $\\theta_{\mathrm{Ein}}$', fontsize=fsize)

ax[0].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[1].legend(loc=(0.04, 0.4), fontsize=fsize)

ax[1].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[2].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)

for j in range(3):
    ax[j].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[j].xaxis.set_minor_locator(MultipleLocator(0.1))

pylab.savefig('../paper/minmu_source_bias.eps')
pylab.show()


