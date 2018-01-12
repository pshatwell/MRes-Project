from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

from Eadyinfo_oc import *

print 'U_0 is:', U_0
print 'Lambda is:', Lambda

print '\nH is', H
print 'L is', L
print 'T is', T

print '\nk is', k
print 'l is', l
print 'mu is', mu

print '\nc is:', c
print 'c^2 is:', np.abs(c)**2
print 'Real part of c is:', c.real
print 'sigma is:', sigma

###############################################################################################

#Set resolution
npts = 200

#Set coordinates
xmin = -np.pi/(k*L) #Extent in x covers one full wavelength
xmax = np.pi/(k*L)

zmin = 0
zmax = 1

tmin = 0
tmax = 51

kmin = 0
kmax = 1e-4

xpos = 0
zpos = 0.5
timestep = 0

xvalues = np.linspace(xmin, xmax, npts)
zvalues = np.linspace(zmin, zmax, npts)
tvalues = np.arange(tmin, tmax, 1)
kvalues = np.linspace(kmin, kmax, npts)

xlength = len(xvalues)
zlength = len(zvalues)
tlength = len(tvalues)
klength = len(kvalues)

###############################################################################################

#This is slopefunction2 but for x = c.real*(T/L)*t at the wave frame origin
def slopefunction(z):
    return -(k*(c.real - Lambda*H*z)*dphi_idz(z) + sigma*dphi_rdz(z))/(Lambda*k*phi_i(z)) - 1

def slopefunction2(x,z,t):
    return -((k*(c.real - Lambda*H*z)*dphi_rdz(z) - sigma*dphi_idz(z))*np.sin(k*(L*x - c.real*T*t)) + (k*(c.real - Lambda*H*z)*dphi_idz(z) + sigma*dphi_rdz(z))*np.cos(k*(L*x - c.real*T*t)))/(Lambda*k*(phi_r(z)*np.sin(k*(L*x - c.real*T*t)) + phi_i(z)*np.cos(k*(L*x - c.real*T*t)))) - 1

slopefunctionvalues = []

slopefunctionvalues2 = np.empty((zlength,xlength))

for i in range(0,zlength,1):
    slopefunctionvalues.append(slopefunction(zvalues[i]))

for i in range(0,zlength,1):
    for j in range(0,xlength,1):
        if np.abs(slopefunction2(x=xvalues[j], z=zvalues[i], t=timestep)) > 1.: #We filter out any values where the parcel slope is steeper than the isentropic slope
            slopefunctionvalues2[i,j] = 0
        else:
            slopefunctionvalues2[i,j] = np.abs(slopefunction2(x=xvalues[j], z=zvalues[i], t=timestep)) #Take absolute value

###############################################################################################

def chi_function(x,z,t):
    return phi_r(z)*np.sin(k*(L*x - c.real*T*t)) + phi_i(z)*np.cos(k*(L*x - c.real*T*t))

chivalues = np.empty((zlength,xlength))

for i in range(0,zlength,1):
    for j in range(0,xlength,1):
        chivalues[i,j] = chi_function(x=xvalues[j], z=zvalues[i], t=timestep)

###############################################################################################

'''
#The value of 0.67 at the steering level seems consistent with wprimevprime script, but, again, it is not symmetric
#about the steering level. For small z, it is NOT TRUE that the parcel slope is close to the isentropic slope.
#This is only a section at x=ct of the contour plot in figure 2.
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$z \; (H)$', fontsize=16)
ax.set_ylabel(r'$|slope^{parcel}/slope^{\bar{\Theta}}|$', fontsize=16)
ax.plot(zvalues, np.abs(slopefunctionvalues)) #Taking absolute value
ax.axhline(y=0.5, ls='dashed', color='black')
ax.axvline(x=0.5, ls='dashed', color='black')
'''

'''This is a contour plot of the parcel slope/isentrope slope in the x-z plane (can't do y-z plane as functions
have no y-dependence). At least it peaks near the steering level now!
Changing timestep only translates this figure zonally.'''
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xlabel(r'$x \; (L)$', fontsize=16)
ax2.set_ylabel(r'$z \; (H)$', fontsize=16)
#contour = ax2.contourf(slopefunctionvalues2, extent=[xmin,xmax,zmin,zmax], origin='lower', aspect='auto', cmap='magma')
contour = ax2.matshow(slopefunctionvalues2, extent=[xmin,xmax,zmin,zmax], origin='lower', aspect='auto', cmap='magma')
plt.colorbar(contour)
ax2.axhline(y=0.5, ls='dashed', color='black')

'''
#x-z structure of vprime. The parcel slope contour in figure 2 blows up where
#this plot equals zero. How do you calculate where this is zero analytically?
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.set_xlabel(r'$x \; (L)$', fontsize=16)
ax3.set_ylabel(r'$z \; (H)$', fontsize=16)
chi_contour = ax3.contourf(chivalues, extent=[xmin,xmax,zmin,zmax], origin='lower', aspect='auto', cmap='magma')
plt.colorbar(chi_contour)
ax3.axhline(y=0, ls='dashed', color='black')
'''


plt.show()
