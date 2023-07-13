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

fig, ax = pylab.subplots(5, 1, figsize=(6, 14))

pylab.subplots_adjust(left=0.23, right=1.00, bottom=0.05, top=1., wspace=0., hspace=0.)

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
qmed_gal = np.median(galpop['q'][()])

for n in range(nsims):

    lenspop = h5py.File(lpopnames[n], 'r')

    laspsmed_arr = np.zeros(ntein)
    laspserr_arr = np.zeros(ntein)

    lmdm5med_arr = np.zeros(ntein)
    lmdm5err_arr = np.zeros(ntein)

    qmed_arr = np.zeros(ntein)
    qerr_arr = np.zeros(ntein)

    def mdm5fitfunc(p):
        return p[0] + p[1] * (galpop['lmobs'][()] - lmobs_piv) + p[2] * (galpop['lreff'][()] - lreff_piv)

    def mdm5errfunc(p):
        return mdm5fitfunc(p) - galpop['lmdm5'][()]

    pmdm5fit = leastsq(mdm5errfunc, [11., 0., 0.])

    lmdm5med_gal = pmdm5fit[0][0]

    gammamed_arr = np.zeros(ntein)
    gammaerr_arr = np.zeros(ntein)

    lmstmed_arr = np.zeros(ntein)
    lmsterr_arr = np.zeros(ntein)

    def gammafitfunc(p):
        return p[0] + p[1] * (galpop['lmobs'][()] - lmobs_piv) + p[2] * (galpop['lreff'][()] - lreff_piv)

    def gammaerrfunc(p):
        return gammafitfunc(p) - galpop['gammadm'][()]

    pgammafit = leastsq(gammaerrfunc, [1.5, 0., 0.])

    gammamed_gal = pgammafit[0][0]

    for i in range(ntein):
        lenscut = galpop[lensdefs[n]][()] & (galpop[teindefs[n]][()] > tein_arr[i])
        lenspopcut = (lenspop['tein_zs'][()] > tein_arr[i])

        nlens = lenscut.sum()
        print('Theta_Ein > %2.1f. %d lenses'%(tein_arr[i], lenscut.sum()))
        laspsmed_arr[i] = np.median(galpop['lasps'][lenscut])
        laspserr_arr[i] = np.std(galpop['lasps'][lenscut])/float(nlens)**0.5

        qmed_arr[i] = np.median(galpop['q'][lenscut])
        qerr_arr[i] = np.std(galpop['q'][lenscut])/float(nlens)**0.5

        # fits for the stellar-halo mass relation

        lmobs_here = galpop['lmobs'][lenscut]
        lreff_here = galpop['lreff'][lenscut]
        lmdm5_here = galpop['lmdm5'][lenscut]
        lmst_here = lenspop['lmst'][lenspopcut].flatten()
        gammadm_here = galpop['gammadm'][lenscut]

        def mdm5fitfunc(p):
            return p[0] + p[1] * (lmobs_here - lmobs_piv) + p[2] * (lreff_here - lreff_piv)

        def mdm5errfunc(p):
            return mdm5fitfunc(p) - lmdm5_here

        pmdm5fit = leastsq(mdm5errfunc, [11., 0., 0.])

        mu_mdm5_here = pmdm5fit[0][0]
        lmdm5_scat = np.std(mdm5fitfunc(pmdm5fit[0]) - lmdm5_here)

        lmdm5med_arr[i] = mu_mdm5_here
        lmdm5err_arr[i] = lmdm5_scat/float(nlens)**0.5

        def gammafitfunc(p):
            return p[0] + p[1] * (lmobs_here - lmobs_piv) + p[2] * (lreff_here - lreff_piv)

        def gammaerrfunc(p):
            return gammafitfunc(p) - gammadm_here

        pgammafit = leastsq(gammaerrfunc, [11., 0., 0.])

        mu_gamma_here = pgammafit[0][0]
        gamma_scat = np.std(gammafitfunc(pgammafit[0]) - gammadm_here)

        gammamed_arr[i] = mu_gamma_here
        gammaerr_arr[i] = gamma_scat/float(nlens)**0.5

        lmstmed_arr[i] = np.median(lmst_here)
        lmsterr_arr[i] = np.std(lmst_here)/float(nlens)**0.5

    ax[0].errorbar(tein_arr, laspsmed_arr, yerr=laspserr_arr, color=colors[n], label=labels[n])
    ax[1].errorbar(tein_arr, lmdm5med_arr, yerr=lmdm5err_arr, color=colors[n], label=labels[n])
    ax[2].errorbar(tein_arr, gammamed_arr, yerr=gammaerr_arr, color=colors[n])#, label=labels[n])
    ax[3].errorbar(tein_arr, qmed_arr, yerr=qerr_arr, color=colors[n], label=labels[n], linestyle=linestyles[n])
    ax[4].errorbar(tein_arr, lmstmed_arr, yerr=lmsterr_arr, color=colors[n], label=labels[n])

    ax[2].axhline(gammamed_gal, color=colors[n], linestyle='--')

ax[0].axhline(laspsmed_gal, color='k', linestyle='--')#, label='Parent population')
ax[0].set_ylabel('Median $\log{\\alpha_{\mathrm{sps}}}$', fontsize=fsize)

ax[0].yaxis.set_major_locator(MultipleLocator(0.05))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.01))

ax[1].axhline(lmdm5med_gal, color='k', linestyle='--', label='General pop.')

ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.02))
ax[1].set_ylabel('$\mu_{\mathrm{DM},0}$ (Mean $\log{M_{\mathrm{DM},5}}$ \n at fixed $M_*^{(\mathrm{obs})}$, $R_{\mathrm{e}}$)', fontsize=fsize)

ax[2].set_ylabel('$\mu_{\gamma,0}$ (Mean $\gamma_{\mathrm{DM},5}$\n at fixed $M_*^{(\mathrm{obs})}$, $R_{\mathrm{e}}$)', fontsize=fsize)

ax[2].yaxis.set_major_locator(MultipleLocator(0.01))
ax[2].yaxis.set_minor_locator(MultipleLocator(0.002))
#ax[2].set_ylim(10.98, 11.36)

ax[3].axhline(qmed_gal, color='k', linestyle='--')
ax[3].set_ylabel('Median $q$', fontsize=fsize)
#ax[3].set_xlabel('Minimum $\\theta_{\mathrm{Ein}}$', fontsize=fsize)

ax[3].yaxis.set_major_locator(MultipleLocator(0.02))
ax[3].yaxis.set_minor_locator(MultipleLocator(0.005))

ax[4].set_ylabel('Median $\lambda_{\mathrm{int}}$', fontsize=fsize)
ax[4].set_xlabel('Minimum $\\theta_{\mathrm{Ein}}$', fontsize=fsize)

ax[4].yaxis.set_major_locator(MultipleLocator(0.02))
ax[4].yaxis.set_minor_locator(MultipleLocator(0.005))

ax[1].legend(loc='upper left', fontsize=fsize)

ax[0].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[1].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[2].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[3].tick_params(axis='both', which='both', direction='in', labelbottom=False, labelsize=fsize, right=True, top=True)
ax[4].tick_params(axis='both', which='both', direction='in', labelsize=fsize, right=True, top=True)

for j in range(5):
    ax[j].xaxis.set_major_locator(MultipleLocator(0.5))
    ax[j].xaxis.set_minor_locator(MultipleLocator(0.1))

pylab.savefig('../paper/minmu_lens_bias.eps')
pylab.show()


