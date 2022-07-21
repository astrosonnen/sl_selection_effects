import numpy as np
from astropy.io import fits as pyfits
import os


# defines a sample of galaxies with the following criteria:
# - in SDSS DR7 
# - in the redshift range 0.19 < z < 0.21
# - above some threshold in rest-frame color (to be defined)
# - with best-fit Sersic index > 2.5

# Stores the following quantities:
# - RA, dec, specz
# - r-band de Vaucouleurs magnitude from Meert et al. (2015)
# - r-band de Vaucouleurs half-light radius from Meert et al. (2015)
# - r-band de Vaucouleurs axis ratio from Meert et al. (2015)
# - Stellar mass from Mendel et al. (2014), renormalized to de Vaucouleurs photometry of Meert et al. (2015)

# sdss table from Meert et al. (2015)

gdir = os.environ.get('GDRIVE')

meert_catalog = pyfits.open(gdir+'/catalogs/sdss_upenn_photometry_mstar/meert_et_al_data_tables_v2/UPenn_PhotDec_CAST.fits')[1]
meert_table = meert_catalog.data

meert_dev = pyfits.open(gdir+'/catalogs/sdss_upenn_photometry_mstar/meert_et_al_data_tables_v2/UPenn_PhotDec_Models_rband.fits')[2].data # de Vaucouleurs fit table
meert_sersic = pyfits.open(gdir+'/catalogs/sdss_upenn_photometry_mstar/meert_et_al_data_tables_v2/UPenn_PhotDec_Models_rband.fits')[3].data # Sersic fit table
meert_serexp = pyfits.open(gdir+'/catalogs/sdss_upenn_photometry_mstar/meert_et_al_data_tables_v2/UPenn_PhotDec_Models_rband.fits')[5].data # Sersic+Exponential fit table
f = open(gdir+'/catalogs/sdss_upenn_photometry_mstar/UPenn_PhotDec_Mstar_mlMendel14.dat', 'r')
meert_mstar_serexp = np.loadtxt(f, usecols=(2, ))

meert_mstar_dev = meert_mstar_serexp - 1./2.5 * (meert_dev['m_bulge'] - meert_serexp['m_tot'])

nmeert = len(meert_table)
meert_indices = np.arange(nmeert)

zmin = 0.19
zmax = 0.21

ramin = 125.
ramax = 228.
decmin = -3.
decmax = 6.

meert_cut = (meert_table['z'] > zmin) & (meert_table['z'] < zmax) & ((meert_dev['finalflag'] & 2**20) == 0) & ((meert_sersic['finalflag'] & 2**20) == 0) & ((meert_serexp['finalflag'] & 2**20) == 0) & (meert_mstar_serexp > 0.) & (meert_sersic['n_bulge'] > 2.5)

meert_ra = meert_table['ra'][meert_cut]
meert_dec = meert_table['dec'][meert_cut]
meert_z = meert_table['z'][meert_cut]
meert_cutindices = meert_indices[meert_cut]

ncut = meert_cut.sum()

print(meert_cut.sum())

# matches with SDSS to get model magnitudes and velocity dispersion
sdss_table = pyfits.open(gdir+'/catalogs/sdss_dr16_legacy_galaxy_modelmags.fits')[1].data

sdss_cut = (sdss_table['z'] > zmin) & (sdss_table['z'] < zmax)
sdss_ra = sdss_table['ra'][sdss_cut]
sdss_dec = sdss_table['dec'][sdss_cut]
sdss_modelgmr = sdss_table['modelMag_g'][sdss_cut] - sdss_table['modelMag_r'][sdss_cut]

keep_ind = []
modelgmr_list = []

for i in range(ncut):
    dists = ((meert_ra[i] - sdss_ra)**2 * np.cos(np.deg2rad(meert_dec[i]))**2 + (meert_dec[i] - sdss_dec)**2)**0.5 * 3600.
    closest = dists.argmin()
    if dists[closest] < 1.: # found match in SDSS catalog
        if sdss_modelgmr[closest] > 1.2:
            keep_ind.append(meert_cutindices[i])
            modelgmr_list.append(sdss_modelgmr[closest])

print(len(keep_ind))

# writes catalog with galaxies that were kept
new_table = meert_table[keep_ind]
newhdu = pyfits.BinTableHDU(data=new_table)

mstar_column = pyfits.Column(name='logMstar_Serexp', format='E', array=meert_mstar_serexp[keep_ind])
mstar_column = pyfits.Column(name='logMstar_deV', format='E', array=meert_mstar_dev[keep_ind])
gmr_column = pyfits.Column(name='model_g-i', format='5E', array=np.array(modelgmr_list))
r_dev = pyfits.Column(name='dev_reff', format='E', unit='arcsec', array=meert_dev['r_bulge'][keep_ind])
m_dev = pyfits.Column(name='dev_rmag', format='E', unit='mag', array=meert_dev['m_bulge'][keep_ind])
ba_dev = pyfits.Column(name='dev_ba', format='E', array=meert_dev['ba_bulge'][keep_ind])
pa_dev = pyfits.Column(name='dev_pa', format='E', array=meert_dev['pa_bulge'][keep_ind])
flag_dev = pyfits.Column(name='dev_finalflag', format='J', array=meert_dev['finalflag'][keep_ind])
r_serexp = pyfits.Column(name='serexp_reff', format='E', unit='arcsec', array=meert_serexp['r_tot'][keep_ind])
m_serexp = pyfits.Column(name='serexp_rmag', format='E', unit='mag', array=meert_serexp['m_tot'][keep_ind])

extra_columns = pyfits.ColDefs([mstar_column, gmr_column, r_dev, m_dev, ba_dev, pa_dev, flag_dev, r_serexp, m_serexp])

new_columns = newhdu.columns + extra_columns

hdu = pyfits.BinTableHDU.from_columns(new_columns)
hdu.writeto('sdss_zbin_latetype_red.fits', overwrite=True)

