from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

from Eadyinfo_oc_muvar import *

print 'U_0 is:', U_0
print 'Lambda is:', Lambda

print '\nH is', H
print 'L is', L
print 'T is', T

###############################################################################################

#Set resolution
npts = 50

#Set coordinates
xmin = -5
xmax = 5

ymin = -1
ymax = 1

zmin = 0
zmax = 1

mu_min = 1.7
mu_max = 2.5

timestep = 10

xvalues = np.linspace(xmin, xmax, npts)
yvalues = np.linspace(ymin, ymax, npts)
zvalues = np.linspace(zmin, zmax, npts)
muvalues = np.linspace(mu_min, mu_max, npts)


xlength = len(xvalues)
ylength = len(yvalues)
zlength = len(zvalues)
mulength = len(muvalues)

###############################################################################################

APEdensity_t = np.empty((zlength,ylength,xlength,mulength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        for m in range(0,xlength,1):
            for n in range(0,mulength,1):
                APEdensity_t[i,j,m,n] = APEconversion(x = xvalues[m], y = yvalues[j], z = zvalues[i], t = timestep, muvar = muvalues[n])

APE_t = []

for n in range(0,mulength,1):
    APE_t.append(np.sum(APEdensity_t[i,j,m,n] for i in range(0,zlength,1) for j in range(0,ylength,1) for m in range(0,xlength,1)))

APE_t_max = max(APE_t)
mupos = APE_t.index(APE_t_max)
mumax = muvalues[mupos]

print '\nmaximum occurs at: (%s, %s)' % (mumax, APE_t_max)

###############################################################################################

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$\mu$', fontsize=16)
ax.set_ylabel(r'$\frac{d}{dt}APE \;$ (arbitrary units)', fontsize = 16)
ax.plot(muvalues, APE_t)
ax.scatter(mumax, APE_t_max, s=15, color='k')

plt.show()
