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
npts = 500

mu_min = 0
mu_max = 2.5

muvalues = np.linspace(mu_min, mu_max, npts)

###############################################################################################

def c_i_function(muvar):
    return 100*U_0*np.sqrt(((1./np.tanh(muvar/2.)) - muvar/2.)*(muvar/2. - np.tanh(muvar/2.)))

def k_mu_function(muvar):
    return 20000*np.sqrt((muvar/L)**2 - l**2)

def sigma_mu_function(muvar):
    return k_mu_function(muvar)*c_i_function(muvar)

###############################################################################################

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(r'$\mu$', fontsize=16)
ax.plot(muvalues, c_i_function(muvalues), label= r'$c_i \mu$')
ax.plot(muvalues, k_mu_function(muvalues), label = r'$k_{\mu}$')
ax.plot(muvalues, sigma_mu_function(muvalues), label = r'$\sigma_{\mu} = k_{\mu} c_i \mu$')
ax.legend(loc='upper left')

plt.show()
