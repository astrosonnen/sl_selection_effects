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

fig, ax = pylab.subplots(4, 1, figsize=(6, 13))

pylab.subplots_adjust(left=0.2, right=1.00, bottom=0.05, top=1., wspace=0., hspace=0.)

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
        laspsmed_arr[i] = np.median(galpop['lasps'][lenscut])
        laspserr_arr[i] = np.std(galpop['lasps'][lenscut])/float(nlens)**0.5

        # fits for the stellar-halo mass relation

        lmobs_here = galpop['lmobs'][lenscut]
        lreff_here = galpop['lreff'][lenscut]
        lmdm5_here = galpop['lmdm5'][lenscut]

        def mdm5fitfunc(p):
            return p[0] + p[1] * (lmobs_here - lmobs_piv) + p[2] * (lreff_here - lreff_piv)

        def mdm5errfunc(p):
            return mdm5fitfunc(p) - lmdm5_here

        pmdm5fit = leastsq(mdm5errfunc, [11., 0., 0.])

        mu_mdm5_here = pmdm5fit[0][0]
        lmdm5_scat = np.std(mdm5fitfunc(pmdm5fit[0]) - lmdm5_here)

        lmdm5med_arr[i] = mu_mdm5_here
        lmdm5err_arr[i] = lmdm5_scat/float(nlens)**0.5

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

    ax[0].errorbar(tein_arr, laspsmed_arr, yerr=laspserr_arr, color=colseq[n], label=labels[n])
    ax[1].errorbar(tein_arr, lmdm5med_arr, yerr=lmdm5err_arr, color=colseq[n], label=labels[n])
    ax[2].errorbar(tein_arr, smagmed_arr, yerr=smagerr_arr, color=colors[n], label=labels[n])
    ax[3].errorbar(tein_arr, sreffmed_arr, yerr=srefferr_arr, color=colors[n], label=labels[n], linestyle=linestyles[n])
    #ax[3].errorbar(tein_arr, nsermed_arr, yerr=nsererr_arr, color=colseq[n], label=labels[n])
    #ax[4].errorbar(tein_arr, qmed_arr, yerr=qerr_arr, color=colseq[n], label=labels[n])

# NEED TO REMOVE INDENTATION
ax[0].axhline(laspsmed_gal, color='k', linestyle='--')#, label='Parent population')
ax[0].set_ylabel('Median $\log{\\alpha_{\mathrm{sps}}}$', fontsize=fsize)

ax[0].yaxis.set_major_locator(MultipleLocator(0.05))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.01))

ax[1].axhline(lmdm5med_gal, color='k', linestyle='--', label='Detectable pop.')
ax[1].set_ylabel('Median $m_{\mathrm{s}}$', fontsize=fsize)

ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.02))

ax[2].axhline(smagmed_det, color='k', linestyle='--')
ax[2].set_ylabel('Median $m_{\mathrm{s}}$', fontsize=fsize)

ax[2].yaxis.set_major_locator(MultipleLocator(0.2))
ax[2].yaxis.set_minor_locator(MultipleLocator(0.05))

ax[3].axhline(sreffmed_det, color='k', linestyle='--')
ax[3].set_ylabel("Median $\\theta_{\mathrm{s}}$\n at $m_{\mathrm{s}}=25\, ('')$", fontsize=fsize)
ax[3].set_xlabel('Minimum $\\theta_{\mathrm{Ein}}$', fontsize=fsize)

#ax[3].yaxis.set_major_locator(MultipleLocator(0.02))
#ax[3].yaxis.set_minor_locator(MultipleLocator(0.005))

ax[0].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
#ax[1].legend(loc='upper left', fontsize=fsize)
ax[2].legend(loc=(0.02, 0.8), fontsize=fsize, framealpha=1.)

ax[1].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[2].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[3].tick_params(axis='both', which='both', direction='in', labelsize=fsize, right=True, top=True)

#for j in range(4):
#    ax[j].xaxis.set_major_locator(MultipleLocator(0.5))
#    ax[j].xaxis.set_minor_locator(MultipleLocator(0.1))

#pylab.savefig('../paper/source_bias.eps')
pylab.show()


