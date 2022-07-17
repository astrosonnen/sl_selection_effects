import numpy as np
from sl_profiles import sersic
from scipy.special import gamma as gfunc
import os


np.random.seed(0)

# prepares a Glafic input file

nre = 11
logre_grid = np.linspace(-1., 0., nre)

ftot = 200.
pix = 0.05
nser = 1.

f = open('preamble.input', 'r')
prelines = f.readlines()
f.close()

nsource = 10000
    
rmax = 1.

   
prefix = 'ftot%d_sourceonly'%(ftot)

lines = prelines.copy()
    
lines.append('\n')
lines.append('start_command\n')
lines.append('\n')

dirname = 'mockdir/%s/'%prefix
    
if not os.path.isdir(dirname):
    os.system('mkdir %s'%dirname)
 
for i in range(nre):

    re = 10.**logre_grid[i]
     
    I0 = ftot/(2.*np.pi*(re/pix)**2*nser/sersic.b(nser)**(2*nser)*gfunc(2.*nser))
    
    lines.append('reset_par prefix %s\n'%(prefix + '_logre%2.1f'%logre_grid[i]))
    lines.append('reset_extend 1 2 %f\n'%I0)
    lines.append('reset_extend 1 7 %f\n'%10.**logre_grid[i])
    lines.append('\n')
    
    lines.append('writeimage_ori\n')
    lines.append('\n')
    
f = open(dirname+'/%s.input'%prefix, 'w')
f.writelines(lines)
f.close()

