import numpy as np
import pylab
from sl_profiles import deVaucouleurs as deV, gnfw
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
from scipy.optimize import brentq
rc('text', usetex=True)


fsize = 18

e_list = [0., 0.1, 0.2, 0.3]
ne = len(e_list)

colseq = pylab.rcParams['axes.prop_cycle'].by_key()['color']

fig, ax = pylab.subplots(figsize=(5.5, 5.5))

pylab.subplots_adjust(left=0., right=1., bottom=0., top=1.)

for i in range(ne):

    e_here = e_list[i]
    q_here = 1. - e_here

    f = open('composite_e%2.1f_crit.dat'%e_here)
    table = np.loadtxt(f)
    f.close()

    xs1 = table[:, 2]
    ys1 = table[:, 3]
    xs2 = table[:, 6]
    ys2 = table[:, 7]

    nseg = len(xs1)

    for j in range(nseg):
        if j==0:
            pylab.plot([xs1[j], xs2[j]], [ys1[j], ys2[j]], color=colseq[i], label='$q=%2.1f$'%q_here)
        else:
            pylab.plot([xs1[j], xs2[j]], [ys1[j], ys2[j]], color=colseq[i])

#ax.set_aspect(1.)

pylab.legend(loc='lower right', fontsize=fsize)

#pylab.axis('off')
pylab.show()

