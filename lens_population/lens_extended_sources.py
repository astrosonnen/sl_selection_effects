import numpy as np
import os
import glafic
import h5py
from simpars import *
from astropy.io import fits as pyfits
from scipy.interpolate import splev
from scipy.optimize import brentq
from sl_profiles import sersic
from scipy.special import gamma as gfunc
from lensdet import detect_lens


modelname = 'fiducial_100sqdeg'
sourcecat = pyfits.open('/Users/alessandro/catalogs/skills_sourceonly_zcut.fits')[1].data
pop = h5py.File('%s_galaxies.hdf5'%modelname, 'r+')

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
prefix = 'tmp'
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

zs_list = []
xpos_list = []
ypos_list = []
nser_list = []
sreff_list = []
sq_list = []
spa_list = []
smag_list = []
avg_mu_list = []
nimg_list = []
nmax_list = []
tein_zs_list = []

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

            # measures detectable source unlensed flux
            def zerofunc(R):
                return ftot * sersic.Sigma(R, nser, sreff/pix_arcsec) - nsigma_pixdet * sky_rms

            if zerofunc(0.) < 0.:
                fdet = 0.
            elif zerofunc(10.*sreff/pix_arcsec) > 0.:
                fdet = ftot
            else:
                Rmax = brentq(zerofunc, 0., 10.*sreff/pix_arcsec)
                fdet = ftot * sersic.M2d(Rmax, nser, sreff/pix_arcsec)

            detection, nimg, nmax, footprint, nmax_footprint = detect_lens(img)

            if detection:
                islens = True

                #glafic.writelens(zs)

                # reads in lens properties
                #lensprop = pyfits.open('tmp_lens.fits')[0].data
                #mu_grid = lensprop[6, :, :]**(-1)

                tein_zs = glafic.calcein2(zs, 0., 0.)
                tein_zs_list.append(tein_zs)

                print('%d is a lens'%i)

                zs_list.append(zs)
                xpos_list.append(xpos)
                ypos_list.append(ypos)
                nser_list.append(nser)
                sreff_list.append(sreff)
                sq_list.append(sq)
                spa_list.append(spa)
                smag_list.append(smag)
                nimg_list.append(nimg)
                nmax_list.append(nmax)

                # creates a noisy version of the image
                img_wnoise = img + np.random.normal(0., sky_rms, img.shape)

                hdr = pyfits.Header()

                # creates an fits file for the lens
                hdr['galno'] = i
                hdr['zlens'] = pop['z'][i]
                hdr['tein_zrf'] = pop['tein'][i]
                hdr['tein_zs'] = tein_zs
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

                # calculates the average magnification over the footprint
                #footprint = img > nsigma_pixdet * sky_rms

                avg_mu = abs(img[footprint]).sum()/fdet
                avg_mu_list.append(avg_mu)

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

# makes file with lenses

lens_file = h5py.File('%s_lenses.hdf5'%modelname, 'w')

lens_file.create_dataset('z', data=pop['z'][islens_samp])
lens_file.create_dataset('index', data=np.arange(nsamp)[islens_samp])
lens_file.create_dataset('lmobs', data=pop['lmobs'][islens_samp])
lens_file.create_dataset('lmstar', data=pop['lmstar'][islens_samp])
lens_file.create_dataset('lasps', data=pop['lasps'][islens_samp])
lens_file.create_dataset('lm200', data=pop['lm200'][islens_samp])
lens_file.create_dataset('lmdm5', data=pop['lmdm5'][islens_samp])
lens_file.create_dataset('r200', data=pop['r200'][islens_samp])
lens_file.create_dataset('lreff', data=pop['lreff'][islens_samp])
lens_file.create_dataset('tein_zref', data=pop['tein'][islens_samp])
lens_file.create_dataset('tein_zs', data=np.array(tein_zs_list))
lens_file.create_dataset('tcaust', data=pop['tcaust'][islens_samp])
lens_file.create_dataset('q', data=pop['q'][islens_samp])
lens_file.create_dataset('rs', data=pop['rs'][islens_samp])
lens_file.create_dataset('gammadm', data=pop['gammadm'][islens_samp])

lens_file.create_dataset('zs', data=np.array(zs_list))
lens_file.create_dataset('xpos', data=np.array(xpos_list))
lens_file.create_dataset('ypos', data=np.array(ypos_list))
lens_file.create_dataset('nser', data=np.array(nser_list))
lens_file.create_dataset('sreff', data=np.array(sreff_list))
lens_file.create_dataset('sq', data=np.array(sq_list))
lens_file.create_dataset('spa', data=np.array(spa_list))
lens_file.create_dataset('smag', data=np.array(smag_list))
lens_file.create_dataset('avg_mu', data=np.array(avg_mu_list))
lens_file.create_dataset('nimg', data=np.array(nimg_list))
lens_file.create_dataset('nmax', data=np.array(nmax_list))

pop.close()
lens_file.close()


