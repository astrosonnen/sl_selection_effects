import numpy as np
from wl_profiles import gnfw, deVaucouleurs as deV
from wl_cosmology import Mpc, c, G, M_Sun
import wl_cosmology
from scipy.interpolate import splrep, splev, splint
from scipy.optimize import brentq
from scipy.stats import truncnorm
from scipy.special import erf
import emcee
import h5py
import sys


zd = 0.3
zs = 1.5
rs_arcsec = 5.
gammadm = 1.5

dd = wl_cosmology.Dang(zd)
ds = wl_cosmology.Dang(zs)
dds = wl_cosmology.Dang(zs, zd)

kpc = Mpc/1000.
arcsec2rad = np.deg2rad(1./3600.)
arcsec2kpc = arcsec2rad * dd * 1000.

rs_phys = rs_arcsec * arcsec2kpc

s_cr = c**2/(4.*np.pi*G)*ds/dds/dd/Mpc/M_Sun*kpc**2 # critical surface mass density, in M_Sun/kpc**2

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

def get_rein_kpc(mstar, reff, m200, rs, gammadm, s_cr):

    xmin = max(deV.rgrid_min*reff, Rfrac_min*rs)

    r200 = (m200*3./200./(4.*np.pi)/rhoc)**(1./3.) * 1000.

    gnfw_norm = m200 / gnfw.M3d(r200, rs, gammadm)

    def zerofunc(x):
        return alpha(x, gnfw_norm, rs, gammadm, mstar, reff, s_cr) - x

    rein_kpc = brentq(zerofunc, 0.1, 100.)
    return rein_kpc

rein_kpc = get_rein_kpc(0., 1., m200_phys, rs_phys, gammadm, s_cr)
print(rein_kpc / arcsec2kpc)

