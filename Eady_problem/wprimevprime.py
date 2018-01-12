from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

from Eadyinfo_oc import *

###############################################################################################

#Set resolution
npts = 100

tmin = 0
tmax = 50

tvalues = np.linspace(tmin,tmax,npts)
tlength = len(tvalues)

ypos = 0
zpos = 0.5
timestep = 10

zmin = 0
zmax = 1

zvalues = np.linspace(zmin,zmax,npts)
zlength = len(zvalues)

prediction = (-f0*Lambda)/(N*N) #From scale analysis of w' and v' at steering level

###############################################################################################

vprimes = []
wprimes = []

#Evaluating vprime at the origin (and at z=0.5) moving with the wave, changing time
for i in range(0, tlength, 1):
    vprimes.append((L/T)*vprime(x=c.real*tvalues[i]*(T/L), y=ypos, z=zpos, t=tvalues[i]))

for i in range(0, tlength, 1):
    wprimes.append((H/T)*wprime(x=c.real*tvalues[i]*(T/L), y=ypos, z=zpos, t=tvalues[i]))


vprimes2 = []
wprimes2 = []

#Evaluating vprime at the origin moving with the wave, changing z
for i in range(0, zlength, 1):
    vprimes2.append((L/T)*vprime(x=c.real*tvalues[i]*(T/L), y=ypos, z=zvalues[i], t=tvalues[i]))

for i in range(0, zlength, 1):
    wprimes2.append((H/T)*wprime(x=c.real*tvalues[i]*(T/L), y=ypos, z=zvalues[i], t=tvalues[i]))


wprimeovervprime = []
for i in range(0, tlength, 1):
    wprimeovervprime.append(wprimes[i]/vprimes[i])

wprimeovervprime2 = []
for i in range(0, zlength, 1):
    wprimeovervprime2.append(wprimes2[i]/vprimes2[i])

'''
#This is just a check for vprimes2 and wprimes2 method:
vprimes3 = []
wprimes3 = []

#Evaluating vprime at the origin moving with the wave, changing z
#Using WForigin functions
for i in range(0, zlength, 1):
    vprimes3.append((L/T)*WForiginvprime(y=ypos, z=zvalues[i], t=tvalues[i]))

for i in range(0, zlength, 1):
    wprimes3.append((H/T)*WForiginwprime(y=ypos, z=zvalues[i], t=tvalues[i]))

wprimeovervprime3 = []
for i in range(0, zlength, 1):
    wprimeovervprime3.append(np.abs(wprimes3[i]/vprimes3[i]))
'''

###############################################################################################

'''As expected, the ratio of wprime to vprime in the frame of the wave is constant with time. However, it is
not equal to (f0*Lambda)/(N^2) as predicted from theory at the steering level. Why is it smaller?'''
ratiofig = plt.figure()
ratioax = ratiofig.add_subplot(111)
ratioax.set_title(r"w'/v' at the steering level, z=H/2")
ratioax.set_xlabel(r'$time \; (T)$', fontsize=16)
ratioax.set_ylabel(r"$|w'/v'|$", fontsize=16)
ratioax.plot(tvalues, np.abs(wprimeovervprime), color='black', label='Calculated') #Take absolute value
ratioax.axhline(y=prediction, ls='dashed', color='red', label='Predicted (slope of isentropes)')
ratioax.set_ylim([0,0.005])
ratioax.legend()

'''Peaks at z=0 for some reason, and decreases monotonically as z increases.
Nothing special about the steering level here. This doesn't make much sense given the boundary conditions
for the vertical velocity. Note that plot should be independent of y and t.'''
ratiofig2 = plt.figure()
ratioax2 = ratiofig2.add_subplot(111)
ratioax2.set_xlabel(r'$z \; (H)$', fontsize=16)
ratioax2.set_ylabel(r"$|w'/v'|$", fontsize=16)
ratioax2.plot(zvalues, np.abs(wprimeovervprime2), color='black') #Take absolute value
ratioax2.axhline(y=prediction, ls='dashed', color='red', label='Slope of isentropes)')
ratioax2.axhline(y=prediction/2., ls='dashed', color='green', label='Half-slope of isentropes')
ratioax2.axvline(x=0.5, ls='dashed', color='black')
ratioax2.set_ylim([0,0.005])
ratioax2.legend()

plt.show()
