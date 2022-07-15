import numpy as np
from sl_profiles import sersic
from scipy.special import gamma as gfunc


ftot = 20.
re_arr = [0.1, 0.2, 0.5, 1.]
pix = 0.05
nser = 1.

for re in re_arr:
    I0 = ftot/(2.*np.pi*(re/pix)**2*nser/sersic.b(nser)**(2*nser)*gfunc(2.*nser))
    print(I0)


