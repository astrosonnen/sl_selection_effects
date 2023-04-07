import h5py
import numpy as np
import sys
from sl_profiles import gnfw, deVaucouleurs as deV
from sl_cosmology import *
from simpars import *


modelname = sys.argv[1]

pop = h5py.File(modelname+'_lenses.hdf5', 'a')

nlens = len(pop['z'][()])

psi1 = np.zeros(nlens)
psi2_pl = np.zeros(nlens)
psi2 = np.zeros(nlens)
psi3 = np.zeros(nlens)
lmst = np.zeros(nlens)

# defines lensing-related functions (for computation of Einstein radius)
def alpha_dm(x, gnfw_norm, rs, gammadm, s_cr):
    # deflection angle (in kpc)
    return gnfw_norm * gnfw.fast_M2d(abs(x), rs, gammadm) / np.pi/x/s_cr

def alpha_star(x, mstar, reff, s_cr): 
    # deflection angle (in kpc)
    return mstar * deV.M2d(abs(x), reff) / np.pi/x/s_cr

kpc = Mpc/1000.
dx = 0.01

for i in range(nlens):

    ds = sl_cosmology.Dang(pop['zs'][i])
    dds = sl_cosmology.Dang(pop['z'][i], pop['zs'][i])
    dd = splev(pop['z'][i], dd_spline)
    s_cr = c**2/(4.*np.pi*G)*ds/dds/dd/Mpc/M_Sun*kpc**2
    arcsec2kpc = np.deg2rad(1./3600.) * dd * 1000.

    rein = pop['tein_zs'][i] * arcsec2kpc

    gnfw_norm = 10.**pop['lm200'][i] / gnfw.fast_M3d(pop['r200'][i], pop['rs'][i], pop['gammadm'][i])

    def alpha(x):
        return alpha_dm(x, gnfw_norm, pop['rs'][i], pop['gammadm'][i], s_cr) + alpha_star(x, 10.**pop['lmstar'][i], 10.**pop['lreff'][i], s_cr)

    psi1[i] = rein
    psi2[i] = (alpha(rein + dx) - alpha(rein - dx))/(2.*dx)
    psi3[i] = (alpha(rein + dx) + alpha(rein - dx) - 2.*alpha(rein))/dx**2

psi2_pl = -psi3/(1. - psi2) * psi1
lmst = (1. - psi2_pl)/(1. - psi2)

pop.create_dataset('psi2', data=psi2)
pop.create_dataset('psi3', data=psi3)
pop.create_dataset('lmst', data=lmst)

