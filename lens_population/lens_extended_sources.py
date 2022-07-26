import numpy as np
import os
import glafic
import h5py
from simpars import *
from astropy.io import fits as pyfits
from scipy.interpolate import splev
from sl_profiles import sersic
from scipy.special import gamma as gfunc
from lensdet import detect_lens


modelname = 'fiducial_100sqdeg'
sourcecat = pyfits.open('/Users/alessandro/catalogs/skills_sourceonly_zcut.fits')[1].data
pop = h5py.File('%s_lenses.hdf5'%modelname, 'r+')

modeldir = 'extended_sims/%s/'%modelname
if not os.path.isdir(modeldir):
    os.system('mkdir %s'%modeldir)

nsamp = pop.attrs['nsamp']
islens_samp = np.zeros(nsamp, dtype=bool)

f = open('%s_sources.cat'%modelname, 'r')
sourcelines = f.readlines()[1:]
f.close()

# primary parameters
omegaM = 0.3
omegaL = 0.7
weos = -1.
hubble = 0.7
prefix = 'out'
xmin = -5.
ymin = -5.
xmax = 5.
ymax = 5.
pix_ext = pix_arcsec
pix_poi = 0.1
maxlev = 5

glafic.init(omegaM, omegaL, weos, hubble, prefix, xmin, ymin, xmax, ymax, pix_ext, pix_poi, maxlev, verb = 0)
glafic.set_secondary('flag_halodensity 2')
glafic.set_secondary('nfw_users 1')
glafic.set_secondary('halodensity 200')

glafic.startup_setnum(2, 1, 0)
glafic.set_lens(1, 'gnfw', 0.3, 1e13, 0.0, 0.0, 0., 90.0, 10., 1.5)
glafic.set_lens(2, 'sers', 0.3, 1e11, 0.0, 0.0, 0., 90.0, 1., 4.)
glafic.set_extend(1, 'sersic', zs_ref, 0.5, 0.3, 0., 0., 0., 0.1, 1.)

for i in range(nsamp):
    line = sourcelines[i].split()
    nsource = int(line[2])
    rmax = float(line[1])
    islens = False
    if nsource > 0:
        arcsec2kpc = np.deg2rad(1./3600.) * splev(pop['z'][i], dd_spline) * 1000.
        reff_kpc = 10.**pop['lreff'][i]
        reff_arcsec = reff_kpc/arcsec2kpc
        rs_arcsec = pop['rs'][i]/arcsec2kpc

        glafic.set_lens(1, 'gnfw', pop['z'][i], 10.**pop['lm200'][i]*hubble, 0., 0., 1. - pop['q'][i], 90., rs_arcsec, pop['gammadm'][i])
        glafic.set_lens(2, 'sers', pop['z'][i], 10.**pop['lmstar'][i]*hubble, 0., 0., 1. - pop['q'][i], 90., reff_arcsec, 4.)

        n = 0
        while not islens and n < nsource:
            sourcestrs = line[3+n].split(',')
            sourceind = int(sourcestrs[0])
            xpos = float(sourcestrs[1])
            ypos = float(sourcestrs[2])

            nser = sourcecat['sersic_n_CM'][sourceind]
            sreff = sourcecat['Re_arcsec_CM'][sourceind]
            sq = sourcecat['axis_ratio_CM'][sourceind]
            spa = sourcecat['PA_random'][sourceind]
            zs = sourcecat['zobs'][sourceind]
            smag = sourcecat['g_SDSS_apparent_corr'][sourceind]

            ftot = 10.**(-2./5.*(smag - zeropoint))
            I0 = ftot/(2.*np.pi*(sreff/pix_arcsec)**2*nser/sersic.b(nser)**(2*nser)*gfunc(2.*nser))

            glafic.set_extend(1, 'sersic', zs, I0, xpos, ypos, 1.-sq, spa, sreff, nser)

            # model_init needs to be done again whenever model parameters are changed
            glafic.model_init(verb = 0)

            img = np.array(glafic.writeimage())

            detection, nimg = detect_lens(img)

            if detection:
                islens = True

                print('%d is a lens'%i)

                """
                # creates an hdf5 file for the lens
                lens_file = h5py.File(modeldir+'lens_%06d.hdf5'%i, 'w')
                lens_file.attrs['z'] = pop['z'][i]
                lens_file.attrs['tein'] = pop['tein'][i]
                lens_file.attrs['reff_kpc'] = 10.**pop['lreff'][i]
                lens_file.attrs['reff_arcsec'] = 10.**pop['lreff'][i]
                lens_file.attrs['lm200'] = pop['lm200'][i]
                lens_file.attrs['lm200'] = pop['lm200'][i]
                lens_file.attrs['lmstar'] = pop['lmstar'][i]
                lens_file.attrs['zs'] = zs
                lens_file.attrs['source_xpos'] = xpos
                lens_file.attrs['source_ypos'] = ypos
                lens_file.attrs['source_nser'] = nser
                lens_file.attrs['source_q'] = sq
                lens_file.attrs['source_pa'] = spa
                lens_file.attrs['source_mag'] = smag
                lens_file.attrs['source_re'] = sreff
                lens_file.attrs['source_index'] = sourceind

                lens_file.create_dataset('img', data=img)

                lens_file.close()
                """

                # creates a noisy version of the image
                img_wnoise = img + np.random.normal(0., sky_rms, img.shape)

                hdr = pyfits.Header()

                # creates an hdf5 file for the lens
                hdr['galno'] = i
                hdr['zlens'] = pop['z'][i]
                hdr['tein'] = pop['tein'][i]
                hdr['reff_ang'] = reff_arcsec
                hdr['reff_kpc'] = reff_kpc
                hdr['lm200'] = pop['lm200'][i]
                hdr['lmstar'] = pop['lmstar'][i]
                hdr['lens_q'] = pop['q'][i]
                hdr['zs'] = zs
                hdr['src_x'] = xpos
                hdr['src_y'] = ypos
                hdr['src_nser'] = nser
                hdr['src_q'] = sq
                hdr['src_pa'] = spa
                hdr['src_mag'] = smag
                hdr['src_re'] = sreff
                hdr['src_ind'] = sourceind
                hdr['nimg'] = nimg

                # calculates the total magnification
                obsftot = img.sum()
                avg_mu = obsftot/ftot

                hdr['avg_mu'] = avg_mu

                phdu = pyfits.PrimaryHDU(header=hdr)

                ihdu = pyfits.ImageHDU(header=hdr, data=img)
                nhdu = pyfits.ImageHDU(header=hdr, data=img_wnoise)

                hdulist = pyfits.HDUList([phdu, ihdu, nhdu])

                hdulist.writeto(modeldir+'/lens_%06d.fits'%i, overwrite=True)

            n += 1

    islens_samp[i] = islens

if 'islens' in pop:
    data = pop['islens']
    data[()] = islens_samp

else:
    pop.create_dataset('islens', data=islens_samp)

pop.close()

