import numpy as np
from qsopars import *
from scipy.interpolate import splrep, splev
import h5py


nqso = 100000

zqso_samp = splev(np.random.rand(nqso), invcum_zqso_spline)
t_samp = np.random.rand(nqso)

ms_samp = np.zeros(nqso)

for i in range(nqso):
    ind = ztoind(zqso_samp[i])
    ms_here = splev(t_samp[i], invcum_phiqso_splines[ind])
    ms_samp[i] = ms_here

output_file = h5py.File('parent_qsos.hdf5', 'w')
output_file.create_dataset('qsomag', data=ms_samp)
output_file.create_dataset('zqso', data=zqso_samp)

