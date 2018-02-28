from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

from Eadyinfo_oc_muvar import *

print 'U_0 (m/s) is:', U_0
print 'Lambda (m/s per m) is:', Lambda

print '\nH (m) is', H
print 'L (m) is', L
print 'T (s) is', T

print '\nMax Eady growth rate (s^(-1)) is approximately:', (-0.31*Lambda*f0)/N #Vallis 9.89
print 'this corresponds to a timescale (s) of:', -N/(0.31*Lambda*f0)
print 'which, in days, is:', -N/(0.31*Lambda*f0*24*60*60)
print 'in units of T this is:', -N/(0.31*Lambda*f0*T)

###############################################################################################

#Set resolution
npts = 50

#Set coordinates
zmin = 0
zmax = 1

mu_min = 1.58
mu_max = 2.5

timestep = 10

zvalues = np.linspace(zmin, zmax, npts)
muvalues = np.linspace(mu_min, mu_max, npts)

zlength = len(zvalues)
mulength = len(muvalues)

###############################################################################################

def E_function(z, muvar):
    return k_function(muvar)*(c_function(muvar).real - Lambda*H*z)*dphi_rdz_function(z, muvar) - sigma_function(muvar)*dphi_idz_function(z, muvar) + k_function(muvar)*Lambda*phi_r_function(z,muvar)

def F_function(z, muvar):
    return k_function(muvar)*(c_function(muvar).real - Lambda*H*z)*dphi_idz_function(z, muvar) + sigma_function(muvar)*dphi_rdz_function(z, muvar) + k_function(muvar)*Lambda*phi_i_function(z,muvar)

def APEconversion_function(z, muvar):
    return -0.25*(f0**2/N**2)*np.exp(2*sigma_function(muvar)*T*timestep)*(F_function(z, muvar)*dphi_rdz_function(z, muvar) - E_function(z, muvar)*dphi_idz_function(z, muvar))
    #return -(F_function(z, muvar)*dphi_rdz_function(z, muvar) - E_function(z, muvar)*dphi_idz_function(z, muvar))

APE_rate_xy = np.empty((zlength, mulength))

#Calculate zonal and meridional mean rate of APE conversion (see project book)
for i in range(0,zlength,1):
    for j in range(0,mulength,1):
        APE_rate_xy[i,j] = APEconversion_function(z=zvalues[i], muvar=muvalues[j])

APE_rate_xyz = []

scaling = 5*10**5

#Calculate vertical mean of the APE conversion rate
for j in range(0,mulength,1):
    APE_rate_xyz.append(scaling*(1./H)*np.sum(APE_rate_xy[i,j] for i in range(0,zlength,1)))


APE_max = max(APE_rate_xyz)
mupos = APE_rate_xyz.index(APE_max)
mumax = muvalues[mupos]

print '\n maximum APE conversion occurs at: (%s, %s)' % (mumax, APE_max)
print '\n length scale (m) of max APE conversion is:', (2*np.pi/mumax)*L

###############################################################################################

def PHT_function(z, muvar):
    return 0.25*((k_function(muvar)*f0)/(alpha*g))*np.exp(2*sigma_function(muvar)*T*timestep)*(phi_i_function(z, muvar)*dphi_rdz_function(z, muvar) - phi_r_function(z, muvar)*dphi_idz_function(z, muvar))
    #return -k_function(muvar)*(phi_i_function(z, muvar)*dphi_rdz_function(z, muvar) - phi_r_function(z, muvar)*dphi_idz_function(z, muvar))

PHT_xy = np.empty((zlength, mulength))

#Calculate zonal and meridional mean poleward heat transport (see project book)
for i in range(0,zlength,1):
    for j in range(0,mulength,1):
        PHT_xy[i,j] = PHT_function(z=zvalues[i], muvar=muvalues[j])

PHT_xyz = []

#Calculate vertical mean of the poleward heat transport
for j in range(0,mulength,1):
    PHT_xyz.append((1./H)*np.sum(PHT_xy[i,j] for i in range(0,zlength,1)))


PHT_max = max(PHT_xyz)
mupos2 = PHT_xyz.index(PHT_max)
mumax2 = muvalues[mupos2]

print '\n maximum PHT occurs at: (%s, %s)' % (mumax2, PHT_max)
print '\n length scale (m) of max PHT is:', (2*np.pi/mumax2)*L

###############################################################################################

sigma_values = []

scaling2 = 10**(-3)

for i in range(0,mulength,1):
    sigma_values.append(scaling2*sigma_function(muvalues[i]))

sigma_max = max(sigma_values)
mupos3 = sigma_values.index(sigma_max)
mumax3 = muvalues[mupos3]

###############################################################################################

c_i_values = []

scaling3 = 5*10**(-8)

for i in range(0,mulength,1):
    c_i_values.append(scaling3*muvalues[i]*c_function(muvalues[i]).imag)

c_i_max = max(c_i_values)
mupos4 = c_i_values.index(c_i_max)
mumax4 = muvalues[mupos4]

###############################################################################################

#Note that APE_rate_xyz, sigma and c_i are all scaled arbitrarily!!

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$\mu$', fontsize=16)
#ax.set_ylabel(r'$\frac{d}{dt}APE \;$ (arbitrary units)', fontsize = 16)
ax.axhline(y=0, ls='dashed', color='k')
ax.axvline(x=1.61, ls='dashed', color='k')
ax.plot(muvalues, APE_rate_xyz, label='APE conversion rate')
ax.plot(muvalues, PHT_xyz, label='Poleward heat transport')
ax.plot(muvalues, sigma_values, label='Growth rate')
ax.plot(muvalues, c_i_values, label='Imaginary phase speed')
ax.scatter(mumax, APE_max, s=15, color='k')
ax.scatter(mumax2, PHT_max, s=15, color='k')
ax.scatter(mumax3, sigma_max, s=15, color='k')
ax.scatter(mumax4, c_i_max, s=15, color='k')
ax.axvline(x=mumax, ls='dotted', color='k')
ax.axvline(x=mumax2, ls='dotted', color='k')
ax.axvline(x=mumax3, ls='dotted', color='k')
#ax.axvline(x=mumax4, ls='dotted', color='k')
ax.legend()

plt.show()
