import numpy as np
from scipy.interpolate import splrep, splev, splint
import sl_cosmology


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


