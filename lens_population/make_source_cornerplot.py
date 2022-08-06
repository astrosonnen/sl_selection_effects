import pylab
import numpy as np
import h5py
from plotters import probcontour
from astropy.io import fits as pyfits
import os
import sys
import h5py
from matplotlib.ticker import MultipleLocator
from matplotlib import rc
rc('text', usetex=True)


fsize = 10
source_color = (1., 0.2, 1.)
lensed_color = 'b'

modelname = 'fiducial_1000sqdeg'

lenspop = h5py.File('%s_lenses.hdf5'%modelname, 'r')
sourcecat = pyfits.open('/Users/alessandro/catalogs/skills_sourceonly_zcut.fits')[1].data

catpars = ['zobs', 'g_SDSS_apparent_corr', 'Re_arcsec_CM', 'sersic_n_CM', 'axis_ratio_CM']
pars = ['zs', 'smag', 'sreff', 'nser', 'sq']
npars = len(pars)

nlens = len(lenspop['z'][()])
nsource = len(sourcecat)

source_samp = {}
lensed_samp = {}
for i in range(npars):
    source_samp[pars[i]] = sourcecat[catpars[i]].copy()
    lensed_samp[pars[i]] = lenspop['%s'%pars[i]][()].copy()

nbins = 20

labels = ['$z_{\mathrm{s}}$', '$m_{\mathrm{s}}$', '$\\theta_{\mathrm{s,e}}$', '$n_{\mathrm{s}}$', '$q_{\mathrm{s}}$']

lims = [(0.8, 2.5), (21., 28.5), (0., 1.2), (0., 5.), (0., 1.)]

major_step = [1, 2, 0.5, 2, 0.5]
minor_step = [0.2, 0.5, 0.1, 0.5, 0.1]

fig = pylab.figure()
pylab.subplots_adjust(left=0.08, right=0.99, bottom=0.12, top=0.99, hspace=0.1, wspace=0.1)
#pylab.figtext(0.45, 0.95, 'Extended model', fontsize=fsize+3, backgroundcolor=(0., 1., 0.))

for i in range(npars):

    ax = fig.add_subplot(npars, npars, (npars+1)*i + 1)

    bins = np.linspace(lims[i][0], lims[i][1], nbins+1)

    sweights = np.ones(nsource)/float(nsource)/(bins[1] - bins[0])
    pylab.hist(source_samp[pars[i]], bins=bins, color=source_color, histtype='stepfilled', weights=sweights, linewidth=2)#, label='General population')

    lweights = np.ones(nlens)/float(nlens)/(bins[1] - bins[0])
    pylab.hist(lensed_samp[pars[i]], bins=bins, color=lensed_color, histtype='step', weights=lweights, linewidth=2)#, label='General population')

    if i==0:
        ylim = pylab.ylim()
        pylab.ylim(ylim[0], ylim[1])

        box = ax.get_position()
        ax.legend(loc='upper right', bbox_to_anchor=(7., 1.), fontsize=fsize, scatterpoints=3)

    ax.set_xlim((lims[i][0], lims[i][1]))
    ax.tick_params(which='both', direction='in', labelrotation=45)
    ax.xaxis.set_major_locator(MultipleLocator(major_step[i]))
    ax.xaxis.set_minor_locator(MultipleLocator(minor_step[i]))

    ax.set_yticks(())
    if i == npars-1:
        ax.set_xlabel(labels[i], fontsize=fsize)
    else:
        ax.tick_params(axis='x', labelbottom=False)

donelabel = False
for j in range(1, npars): # loops over rows
    if j == npars-1:
        xvisible = True
    else:
        xvisible = False

    for i in range(j): # loops over columns
        ax = pylab.subplot(npars, npars, npars*j+i+1)

        probcontour(source_samp[pars[i]], source_samp[pars[j]], color=source_color, style='filled', linewidths=2)
        probcontour(lensed_samp[pars[i]], lensed_samp[pars[j]], color=lensed_color, style='lines', linewidths=2, smooth=5)
        #probcontour(lens_samp[lenspars[i]], lens_samp[lenspars[j]], color=color, style='solid')
        #ax.scatter(lensed_samp[pars[i]], lensed_samp[pars[j]], color=color, s=3, linewidth=0)

        ax.set_xlim(lims[i])
        ax.set_ylim(lims[j])
        ax.tick_params(which='both', direction='in', labelsize=fsize, labelrotation=45)

        ax.xaxis.set_major_locator(MultipleLocator(major_step[i]))
        ax.xaxis.set_minor_locator(MultipleLocator(minor_step[i]))

        ax.yaxis.set_major_locator(MultipleLocator(major_step[j]))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_step[j]))

        if i == 0:
            yvisible = True
            ax.set_ylabel(labels[j], fontsize=fsize)
            #if j == 1:
            #    box = ax.get_position()
            #    ax.legend(loc='upper right', bbox_to_anchor=(5., 2.), fontsize=fsize)

        else:
            ax.tick_params(axis='y', labelleft=False)

        if xvisible:
            ax.set_xlabel(labels[i], fontsize=fsize)
        else:
            ax.tick_params(axis='x', labelbottom=False)

pylab.savefig('../paper/source_cornerplot.eps')
pylab.show()

