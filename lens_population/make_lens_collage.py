import numpy as np
from PIL import Image, ImageDraw
from astropy.io import fits as pyfits
import pyplz_rgbtools
import h5py
from simpars import *
from lensdet import detect_lens

modelname = 'fiducial_100sqdeg'
lens_file = h5py.File('%s_lenses.hdf5'%modelname, 'r')

nshow = 100

npix = 100

vmax = 5.* sky_rms
vmin = -1.*sky_rms
vdet = nsigma_pixdet*sky_rms

fullim = Image.new('L', (4*npix, nshow*npix), 'black')

for i in range(nshow):

    hdulist = pyfits.open('extended_sims/%s/lens_%06d.fits'%(modelname, lens_file['index'][i]))
    #hdulist = pyfits.open('extended_sims/%s/lens_%06d.fits'%(modelname, 30759))
    img = hdulist[1].data
    img_wnoise = hdulist[2].data

    islens, nimg, nimg_max, footprint, nmax_footprint = detect_lens(img)

    fp = 0.*img
    nm = 0.*img

    fp[footprint] = 255.
    nm[nmax_footprint] = 255.

    sim = Image.new('L', (npix, npix), 'black')
    nim = Image.new('L', (npix, npix), 'black')
    fim = Image.new('L', (npix, npix), 'black')
    nmim = Image.new('L', (npix, npix), 'black')

    img[img < vdet] = vmin
    img[img > vmax] = vmax
    img += vmin
    img *= 255/(vmax - vmin)
    
    img_wnoise[img_wnoise < vmin] = vmin
    img_wnoise[img_wnoise > vmax] = vmax
    img_wnoise += vmin
    img_wnoise *= 255/(vmax - vmin)
     
    sim.putdata(img.flatten())
    nim.putdata(img_wnoise.flatten()) 
    fim.putdata(fp.flatten()) 
    nmim.putdata(nm.flatten()) 

    draw = ImageDraw.Draw(sim)
    #draw.text((5, 10), 'zs=%3.2f'%lens_file['zs'][i], 'white')
    #draw.text((5, 30), 'zl=%3.2f'%lens_file['z'][i], 'white')
    draw.text((5, 10), 'mu=%2.1f'%lens_file['avg_mu'][i], 'white')
    draw.text((50, 10), 'tein=%2.1f'%lens_file['tein_zs'][i], 'white')
    draw.text((5, 30), 'nmax=%d'%lens_file['nmax'][i], 'white')
    draw.text((50, 30), 'nimg=%d'%lens_file['nimg'][i], 'white')
    draw.text((5, 80), '%06d'%lens_file['index'][i], 'white')

    fullim.paste(sim, (0, i*npix))
    fullim.paste(nim, (npix, i*npix))
    fullim.paste(fim, (2*npix, i*npix))
    fullim.paste(nmim, (3*npix, i*npix))

fullim.save('%s_lens_collage.png'%modelname)

