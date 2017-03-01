import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sympy import *
from scipy.integrate import odeint

from sympyEadyinfo import *

#TEST SCRIPT USING SYMPY TO CALCULATE STREAMFUNCTION DERIVATIVES 

############################## SET UP COORDS FOR PLOTTING ################################################################

xvalues = np.linspace(0,10,50)
zvalues = np.linspace(-1,1,50)

xlength = len(xvalues)
zlength = len(zvalues)

time = 0

print 'L is (m):', L
print 'H is (m):', H
print 'H_R is (m):', H_R
print 'H/H_R is:', Hratio
print 'T is (s):', T
print 'k is (m^-1):', k
print 'sigma_max is (s^-1):', sigma_max
print 'N is (s^-1):', N
print 'U is (m*s^-1):', U
print 'velocity shear is (s^-1):', shear
print 'velocity scaling (to dimensionalise) is:', L*f0

############################## CREATE MATRIX OF PHIPRIME VALUES ################################################################

#phiprime is the streamfunction perturbation

#Initialise matrix of zeros for x-z plane values for phiprime
phiprime_matrix = np.zeros((zlength,xlength))

#Evaluate phiprime matrix values
for i in range(0, zlength, 1):
    for j in range(0, xlength, 1):
        phiprime_matrix[i,j] = phiprime.subs(x,xvalues[j]).subs(z,zvalues[i]).subs(t,0).evalf()

########################### FIND QUANTITIES FROM STREAMFUNCTION ###################################################################

vprime = (1./(L*(f0**2)))*phiprime_x #meridional perturbation velocity, Gill 12.9.3
wprime = (1./(H*f0))*(1./N**2)*(shear*phiprime_x - phiprime_zt - shear*z*phiprime_zx) #vertical perturbation velocity, Gill 12.9.6
thetaprime = (1./(alpha*g))*phiprime_z #potential temperature perturbation, Gill 13.2.5

print 'wprime is:', wprime
print 'phiprime_x is:', phiprime_x
print 'phiprime_zt is:', phiprime_zt
print 'phiprime_zx is:', phiprime_zx

print 'wprime at z = -H is:', wprime.subs(z,-1).subs(t,1).subs(x, xlength/2.).evalf()

#print 'vprime is:', vprime
#print 'wprime is:', wprime

#Initialise matrix of zeros for vprime values in x-z plane
vprime_matrix = np.zeros((zlength,xlength))

#Evaluate values for vprime matrix
for i in range(0, zlength, 1):
    for j in range(0, xlength, 1):
        vprime_matrix[i,j] = vprime.subs(x,xvalues[j]).subs(z,zvalues[i]).subs(t,0).evalf()


#Initialise matrix of zeros for wprime values in x-z plane
wprime_matrix = np.zeros((zlength,xlength))

#Evaluate values for wprime matrix
for i in range(0, zlength, 1):
    for j in range(0, xlength, 1):
        wprime_matrix[i,j] = wprime.subs(x,xvalues[j]).subs(z,zvalues[i]).subs(t,0).evalf()


#Initialise matrix of zeros for wprime values in x-z plane
thetaprime_matrix = np.zeros((zlength,xlength))

#Evaluate values for wprime matrix
for i in range(0, zlength, 1):
    for j in range(0, xlength, 1):
        thetaprime_matrix[i,j] = thetaprime.subs(x,xvalues[j]).subs(z,zvalues[i]).subs(t,0).evalf()

##############################################################################################

fig = plt.figure()
plt.set_cmap('inferno')
ax1 = fig.add_subplot(221)
ax1.set_title('Streamfunction')
ax1.set_xlabel('x (L)')
ax1.set_ylabel('z (H)')
phicontour = ax1.contourf(phiprime_matrix, origin='lower', extent=[0,10,-1,1], aspect='auto')
plt.colorbar(phicontour)

ax2 = fig.add_subplot(222)
ax2.set_title('v')
ax2.set_xlabel('x (L)')
ax2.set_ylabel('z (H)')
vcontour = ax2.contourf(vprime_matrix, origin='lower', extent=[0,10,-1,1], aspect='auto')
plt.colorbar(vcontour)

ax3 = fig.add_subplot(223)
ax3.set_title('w')
ax3.set_xlabel('x (L)')
ax3.set_ylabel('z (H)')
wcontour = ax3.contourf(wprime_matrix, origin='lower', extent=[0,10,-1,1], aspect='auto')
plt.colorbar(wcontour)

ax4 = fig.add_subplot(224)
ax4.set_title('theta')
ax4.set_xlabel('x (L)')
ax4.set_ylabel('z (H)')
thetacontour = ax4.contourf(thetaprime_matrix, origin='lower', extent=[0,10,-1,1], aspect='auto')
plt.colorbar(thetacontour)

plt.show()
