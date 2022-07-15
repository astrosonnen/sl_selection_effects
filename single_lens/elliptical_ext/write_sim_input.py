import numpy as np
from sl_profiles import sersic
from scipy.special import gamma as gfunc


np.random.seed(0)

# prepares a Glafic input file

ftot = 200.
re = 0.2
pix = 0.05
nser = 1.

prefix = 'ftot%d_re%2.1f'%(ftot, re)

f = open('preamble.input', 'r')
lines = f.readlines()
f.close()

I0 = ftot/(2.*np.pi*(re/pix)**2*nser/sersic.b(nser)**(2*nser)*gfunc(2.*nser))

lines.append('\n')
lines.append('start_command\n')
lines.append('\n')
lines.append('reset_extend 1 2 %f\n'%I0)
lines.append('\n')

# generates N sources within a circle of radius 1''
nsource = 1000

rmax = 1.
r = np.random.rand(nsource)**0.5 * rmax
phi = 2.*np.pi*np.random.rand(nsource)
x = r * np.cos(phi)
y = r * np.sin(phi)

for n in range(nsource):
    lines.append('reset_par prefix %s_%04d\n'%(prefix, n))
    lines.append('reset_extend 1 3 %f\n'%(x[n]))
    lines.append('reset_extend 1 4 %f\n'%(y[n]))
    lines.append('writeimage\n')
    lines.append('reset_par prefix %s_%04d_wnoise\n'%(prefix, n))
    lines.append('writeimage 0. 1.\n')
    lines.append('\n')

f = open('mockdir/%s.input'%prefix, 'w')
f.writelines(lines)
f.close()

