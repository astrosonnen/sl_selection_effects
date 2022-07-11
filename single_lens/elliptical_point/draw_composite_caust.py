import numpy as np
import os
#import glafic
import h5py
from wl_profiles import gnfw, deVaucouleurs as deV
from wl_cosmology import Mpc, c, G, M_Sun
import wl_cosmology
from scipy.optimize import brentq


np.random.seed(0)

# primary parameters
omegaM = 0.3
omegaL = 0.7
weos = -1.
hubble = 0.7
prefix = 'out'
xmin = -3.
ymin = -3.
xmax = 3.
ymax = 3.
pix_ext = 0.2
pix_poi = 0.1
maxlev = 5

#glafic.init(omegaM, omegaL, weos, hubble, prefix, xmin, ymin, xmax, ymax, pix_ext, pix_poi, maxlev, verb = 0, nfw_users=1, flag_hodensity=2, hodensity=200.)

# defines the lens parameters
zd = 0.3
zs = 1.5

dd = wl_cosmology.Dang(zd)
ds = wl_cosmology.Dang(zs)
dds = wl_cosmology.Dang(zs, zd)

kpc = Mpc/1000.
arcsec2rad = np.deg2rad(1./3600.)
arcsec2kpc = arcsec2rad * dd * 1000.

s_cr = c**2/(4.*np.pi*G)*ds/dds/dd/Mpc/M_Sun*kpc**2 # critical surface mass density, in M_Sun/kpc**2

# I want a lens with a 1" Einstein radius
rein_arcsec = 1.
rein_phys = rein_arcsec * arcsec2kpc
mein_phys = np.pi * rein_phys**2 * s_cr

# and with Reff = rein
reff_arcsec = rein_arcsec
reff_phys = reff_arcsec * arcsec2kpc

# and with f_dm = 0.5
f_dm = 0.5
mstar_phys = (1. - f_dm) * mein_phys / deV.M2d(reff_phys, reff_phys)

# and with dark matter slope of 1.5
gammadm = 1.5

# Finds the virial mass of the halo
rs_phys = 100.
rs_arcsec = rs_phys / arcsec2kpc

mein_dm = f_dm * mein_phys
gnfw_norm = mein_dm / gnfw.M2d(rein_phys, rs_phys, gammadm)

rhoc = wl_cosmology.rhoc(zd)

def r200_zerofunc(r200):
    m3d_here = gnfw_norm * gnfw.M3d(r200, rs_phys, gammadm)
    volume = 4./3.*np.pi * (r200/1000.)**3
    avg_rho = m3d_here / volume
    return avg_rho - 200.*rhoc

r200_phys = brentq(r200_zerofunc, 10., 1000.)
m200_phys = gnfw_norm * gnfw.M3d(r200_phys, rs_phys, gammadm)

m200_glafic = m200_phys * hubble
mstar_glafic = mstar_phys * hubble

print('GNFW component parameters. Mvir=%4.3e, rs=%4.3f arcsec'%(m200_glafic, rs_arcsec))
print('Sersic component parameters. Mtot=%4.3e, reff=%4.3f arcsec'%(mstar_glafic, reff_arcsec))

"""

reff_kpc = 7.
lmstar_phys = 11.5
mstar_phys = 10.**lmstar_phys
rs_phys = 100.
m200_phys = 1e13
gammadm = 1.5

rs_arcsec = rs_phys / arcsec2kpc
reff_arcsec = reff_phys / arcsec2kpc

rs_phys = rs_arcsec * arcsec2kpc


M200_glafic = 1e13
h = 0.7
m200_phys = M200_glafic / h
lm200 = np.log10(m200_phys)

rhoc = wl_cosmology.rhoc(zd)

# defines lensing-related functions
def alpha_dm(x, gnfw_norm, rs, gammadm, s_cr):
    # deflection angle (in kpc)
    return gnfw_norm * gnfw.fast_M2d(abs(x), rs, gammadm) / np.pi/x/s_cr

def alpha_star(x, mstar, reff, s_cr): 
    # deflection angle (in kpc)
    return mstar * deV.M2d(abs(x), reff) / np.pi/x/s_cr

def alpha(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr):
    return alpha_dm(x, gnfw_norm, rs, gammadm, s_cr) + alpha_star(x, mstar, reff, s_cr)

def kappa(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr): 
    # dimensionless surface mass density
    return (mstar * deV.Sigma(abs(x), reff) + gnfw_norm * gnfw.fast_Sigma(abs(x), rs, gammadm))/s_cr
   
def mu_r(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr):
    # radial magnification
    return (1. + alpha(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr)/x - 2.*kappa(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr))**(-1)

def mu_t(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr):
    # tangential magnification
    return (1. - alpha(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr)/x)**(-1)

dx = 0.0001
dx_search = 0.001

Rfrac_min = gnfw.R_grid[0]
Rfrac_max = gnfw.R_grid[-1]

def get_reff_phys(mstar, reff, m200, rs, gammadm, s_cr):

    xmin = max(deV.rgrid_min*reff, Rfrac_min*rs)

    r200 = (m200*3./200./(4.*np.pi)/rhoc)**(1./3.) * 1000.

    gnfw_norm = m200 / gnfw.M3d(r200, rs, gammadm)

    def zerofunc(x):
        return alpha(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr) - x

    reff_phys = brentq(zerofunc, 0.1, 100.)
    return reff_phys

# sets parameters of the power-law lens
rein = 1.
gamma0 = 2.
q = 0.8
e = 1. - q
zd = 0.3
zs = 1.5

nms = 5
dms_grid = np.linspace(-2., 2., nms)

ngamma = 9
gamma_grid = np.linspace(1.6, 2.4, ngamma)

cs_grid = np.zeros((nms, ngamma))
quadcs_grid = np.zeros((nms, ngamma))

glafic.startup_setnum(1, 0, 1)
glafic.set_lens(1, 'pow', zd, zs, 0.0, 0.0, e, 0.0, rein, gamma0)
glafic.set_point(1, 1.5, 0.5, 0.)

# model_init needs to be done again whenever model parameters are changed
glafic.model_init(verb = 0)

# sets parameters of the source position
npoint = 10000
beta_max = 2.
t = np.random.rand(npoint)
phi = 2.*np.pi*np.random.rand(npoint)
beta = beta_max * t**0.5
beta_x = beta * np.cos(phi)
beta_y = beta * np.sin(phi)

source_plane_area = np.pi * beta_max**2

ndet_grid = np.zeros((nms, ngamma, npoint))

for g in range(ngamma):
    print(g)
    glafic.set_lens(1, 'pow', zd, zs, 0.0, 0.0, e, 0.0, rein, gamma_grid[g])
    glafic.model_init(verb = 0)

    nlens_det = np.zeros(nms, dtype=int)
    nquad_det = np.zeros(nms, dtype=int)

    for n in range(npoint):
    
        a = glafic.point_solve(1.5, beta_x[n], beta_y[n], verb = 0)
        nimg = len(a)
    
        # loop over source magnitudes
        for m in range(nms):
            nimg_det = 0
            for imgno in range(nimg):
                absmu = abs(a[imgno][2])
                if -2.5*np.log10(absmu) < -dms_grid[m]:
                    nimg_det += 1
            if nimg_det > 1:
                nlens_det[m] += 1
                if nimg_det > 3:
                    nquad_det[m] += 1

            ndet_grid[m, g, n] = nimg_det

    cs_grid[:, g] = source_plane_area * nlens_det / npoint
    quadcs_grid[:, g] = source_plane_area * nquad_det / npoint

glafic.quit()

cs_output = h5py.File('ell%2.1f_powerlaw_point_crosssect.hdf5'%e, 'w')
cs_output.attrs['e'] = e

cs_output.create_dataset('gamma_grid', data=gamma_grid)
cs_output.create_dataset('dms_grid', data=dms_grid)
cs_output.create_dataset('cs_grid', data=cs_grid)
cs_output.create_dataset('quadcs_grid', data=quadcs_grid)
cs_output.create_dataset('ndet_grid', data=ndet_grid)
"""

