import numpy as np
from sl_profiles import sersic
from scipy.special import gamma as gfunc
from astropy.io import fits as pyfits


ftot = 1.
re = 0.5
pix = 0.1
nser = 1.

I0 = ftot/(2.*np.pi*(re/pix)**2*nser/sersic.b(nser)**(2*nser)*gfunc(2.*nser))
print(I0)

# double-check

sci = pyfits.open('stdsource_source.fits')[0].data

print(sci.sum())

