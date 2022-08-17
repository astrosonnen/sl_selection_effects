import numpy as np
from scipy.interpolate import splrep, splev, splint
import sl_cosmology
from sl_cosmology import Mpc, M_Sun, c, G
from astropy.io import fits as pyfits


# number density of background sources
nbkg = 70. # per square arcminute

# photometric zeropoint
zeropoint = 25.

# pixel size
pix_arcsec = 0.1

# lens detection pars
nsigma_pixdet = 2.
sky_rms = 0.00433
min_angle = 90.
minmag = 3.
npix_min = 10

# stellar mass function
phi_muz = 1.009e-3
alpha_muz = -0.92
lmstar_muz = 11.21

# lens redshift range
zmin = 0.1
zmax = 0.7

# mass-size relation
mu_R = 1.20
beta_R = 0.63
sigma_R = 0.14
lreff_piv = 1.2

# halo mass distribution
mu_h = 13.
beta_h = 1.

# initial concentration before contraction
c200_0 = 5.

# stellar population synthesis mismatch parameter (log of the mean)
mu_sps = 0.1

lmobs_piv = 11.4
lmstar_piv = 11.5
lmobs_min = 11.
lmobs_max = 12.5

# axis ratio distribution
alpha_q = 6.28
beta_q = 2.05

fourpi_volume = 4*np.pi/3. * (sl_cosmology.comovd(zmax)**3 - sl_cosmology.comovd(zmin)**3)

# observed stellar mass function
def lmobsfunc(lmobs):
    return phi_muz * 10.**((lmobs - lmstar_muz) * (alpha_muz + 1)) * np.exp(-10.**(lmobs - lmstar_muz))

nmobs = 101
lmobs_grid = np.linspace(lmobs_min, lmobs_max, nmobs)

lmobsfunc_spline = splrep(lmobs_grid, lmobsfunc(lmobs_grid))
# number of galaxies within a cubic comoving Mpc
ngal_1Mpc3 = splint(lmobs_min, lmobs_max, lmobsfunc_spline)

cumfunc_lmobs_grid = 0.*lmobs_grid
for i in range(nmobs):
    cumfunc_lmobs_grid[i] = splint(lmobs_min, lmobs_grid[i], lmobsfunc_spline) / ngal_1Mpc3

invcum_lmobs_spline = splrep(cumfunc_lmobs_grid, lmobs_grid)

# defines the redshift distribution function
nz = 71
zgrid = np.linspace(zmin, zmax, nz)
zfunc_grid = 0.*zgrid
for i in range(nz):
    zfunc_grid[i] = sl_cosmology.comovd(zgrid[i])**2 * sl_cosmology.dcomovdz(zgrid[i])

zfunc_spline = splrep(zgrid, zfunc_grid)
zfunc_norm = splint(zmin, zmax, zfunc_spline)

cumfunc_z_grid = 0.*zgrid
for i in range(nz):
    cumfunc_z_grid[i] = splint(zmin, zgrid[i], zfunc_spline) / zfunc_norm

invcum_z_spline = splrep(cumfunc_z_grid, zgrid)

# makes a spline for the angular diameter distance
dd_grid = 0.*zgrid
for i in range(nz):
    dd_grid[i] = sl_cosmology.Dang(zgrid[i])

dd_spline = splrep(zgrid, dd_grid)

# and one for the critical surface mass density at the highest source redshift
zs_ref = 2.5

ds_ref = sl_cosmology.Dang(zs_ref)

dds_grid = 0.*zgrid
for i in range(nz):
    dds_grid[i] = sl_cosmology.Dang(zgrid[i], zs_ref)

kpc = Mpc/1000.

# S_cr at the reference source redshift, as a function of lens redshift
s_cr_grid = c**2/(4.*np.pi*G)*ds_ref/dds_grid/dd_grid/Mpc/M_Sun*kpc**2 # units: M_Sun/kpc**2

s_cr_spline = splrep(zgrid, s_cr_grid)

# Quasar luminosity function (from Manti et al. 2017)

def lphistar_func(z):
    return -6.0991 + 0.0209*z + 0.0171*z**2

def Mstar_func(z):
    return -22.5216 + 1.6510*z + 0.2869*z**2

alpha_Q = -1.35
beta_Q = -3.23

qso_hdulist = pyfits.open('optical_nir_qso_sed_001.fits')

qso_data = qso_hdulist[1].data
qso_wav = qso_data['wavelength']
qso_nu = 1./qso_wav
qso_flambda = qso_data['flux']
qso_restfnu = qso_flambda * qso_wav**2 # I don't care about normalization here

ref_lam = 1450. # Wavelength at which the Manti et al. (2017) luminosity function is defined
restuv_window = (qso_wav > ref_lam-50.) & (qso_wav < ref_lam+50.)

qso_uvnorm = np.median(qso_restfnu[restuv_window])

absmag_0 = -2.5*np.log10(qso_uvnorm/(1e-5)**2)

# loads the SDSS i-band filter
f = open('i_SDSS.res', 'r')
iband_wav, iband_t = np.loadtxt(f, unpack=True)
f.close()

iband_nu = 1./iband_wav
iband_spline = splrep(np.flipud(iband_nu), np.flipud(iband_t))

# at each source redshift, I calculate the transformation from rest-frame UV absolute magnitude to
# observed frame i-band magnitude

nzs = 18
zs_grid = np.linspace(0.8, zs_ref, nzs)
deltamag = 0. * zs_grid
for i in range(nzs):
    dlum = sl_cosmology.Dang(zs_grid[i]) * (1. + zs_grid[i])**2

    scale = 10.**(



