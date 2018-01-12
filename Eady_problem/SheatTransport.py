from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
from matplotlib import rc
rc('text', usetex=True)

from Eadyinfo_oc import *

###############################################################################################

'''DOUBLE CHECK ALL THIS. STILL DOESN'T LOOK AS EXPECTED.'''

#Set resolution
npts = 20

#Define coordinates
xmin = -8
xmax = 8

ymin = -1
ymax = 1

zmin = 0
zmax = 1

tmin = 20
tmax = 31

xvalues = np.linspace(xmin, xmax, npts)
yvalues = np.linspace(ymin, ymax, npts)
zvalues = np.linspace(zmin, zmax, npts)

tvalues = np.arange(tmin, tmax, 1)

xlength = len(xvalues)
ylength = len(yvalues)
zlength = len(zvalues)
tlength = len(tvalues)

###############################################################################################

#Slope (angle) of mean isentropic surfaces
gamma = np.arctan((dthetady*L)/(dthetadz*H))

print 'gamma is:', gamma

#Angle values for inclination of slanted sheet
alpha_values = np.linspace(0.001,gamma,npts)
print 'alpha_values are:', alpha_values

#Gradient of the mean isentropic surfaces
grad_theta = np.array((0,dthetady*L,dthetadz*H))

print 'grad_theta is:', grad_theta

#Vector along mean isentropic surface
#Note grad_theta dotted with s0 is zero, as expected
s0 = np.array((0,-dthetadz*H,dthetady*L))

s = np.zeros((npts,3))
s_hats = np.zeros((npts,3))

s_y = 1. #y-component of vector along slanted sheet, held constant

#Create array of vectors for set of slanted sheets angled at each value of alpha
for i in range(npts):
    s[i] = np.array((0, -s_y, s_y*np.tan(alpha_values[i])))
    s_hats[i] = s[i]/la.norm(s[i]) #vectors are normalised

print 's_hats are:', s_hats

normal = np.zeros((npts,3))
normal_hats = np.zeros((npts,3))

#Create array of normals for the slanted sheets
for i in range(npts):
    normal[i] = np.array((0,s[i,2],-s[i,1]))
    normal_hats[i] = normal[i]/la.norm(normal[i])

print 'normal_hats are:', normal_hats

###############################################################################################

#Create array for velocity fields in slanted sheet
velocityfield_p = np.zeros((len(normal_hats), ylength, xlength, tlength, 3))

velocityfield_p2 = np.zeros((len(normal_hats), ylength, xlength, tlength, 2))

#Create projected velocity fields for Green motion for sheets defined by
#normals from the set of normal_hats
for sheet in range(len(normal_hats)):

    #Equation for plane
    planezvalues = 0.5 - yvalues*(normal_hats[sheet,1]/normal_hats[sheet,2])
    planezlength = len(planezvalues)

    planeyvalues = (0.5 - planezvalues)*(normal_hats[sheet,2]/normal_hats[sheet,1])
    planeylength = len(planeyvalues)

    #Evaluate 3d velocity values on the sheet
    for i in range(0, planeylength, 1):
        for j in range(0, xlength, 1):
            for m in range(0, tlength, 1):
                velocityfield_p[sheet,i,j,m,0] = uprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m]) + umeanflow(z=planezvalues[i]) - c.real*(T/L) #Subtracting phase speed, so we're in wave frame
                velocityfield_p[sheet,i,j,m,1] = vprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m])
                velocityfield_p[sheet,i,j,m,2] = wprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m])

    #Project the 3d velocity vectors from velocityfield_p onto 2d plane

    #These define the vectors that span the plane
    xtilde = np.array((1,0,0))
    ytilde = np.cross(normal_hats[sheet],xtilde)

    for i in range(0, planeylength, 1):
        for j in range(0, xlength, 1):
            for m in range(0, tlength, 1):
                velocityfield_p2[sheet,i,j,m,0] = np.dot(velocityfield_p[sheet,i,j,m], xtilde)
                velocityfield_p2[sheet,i,j,m,1] = np.dot(velocityfield_p[sheet,i,j,m], ytilde)

'''velocityfield_p2 contains every 2d projected velocity field for each slanted sheet at each timestep'''

###############################################################################################

v_max_values = np.zeros((len(normal_hats),tlength))

for sheet in range(len(normal_hats)):
    for m in range(0, tlength, 1):
        v_max_values[sheet,m] = (L/T)*np.amax(velocityfield_p2[sheet,:,:,m,1]) #Note (L/T) factor to dimensionalise (m/s)

'''v_max_values contains the max meridional velocity on the sheet,
for each slanted sheet and each timestep'''

###############################################################################################

theta_gradients = np.zeros((len(s_hats), tlength))

for sheet in range(len(s_hats)):
    for m in range(0, tlength, 1):
        theta_gradients[sheet,m] = np.dot(grad_theta, s_hats[sheet])

'''theta_gradients contains rate of change of theta along meridional
direction of sheet, for each slanted sheet'''

###############################################################################################

temp_heat_transport = rho_0*c_p*H*L*L*np.multiply(v_max_values, theta_gradients)

heat_transport = np.zeros((len(s_hats), tlength))

print 'cosine of alpha_values are:', np.cos(alpha_values)

for i in range(len(alpha_values)):
    heat_transport[i,:] = temp_heat_transport[i,:]*np.cos(alpha_values[i])

print 'shape of heat_transport is:', heat_transport.shape

###############################################################################################

fig1 = plt.figure()
vmaxax = fig1.add_subplot(121)
vmaxax.set_title(r'$v_{max}^{sheet} \; (cm/s)$', fontsize=18)
vmaxax.set_xlabel(r'$Time \; (T)$', fontsize=16)
vmaxax.set_ylabel(r'$\alpha \; (rad)$', fontsize=16)
vmaxax.axhline(y=gamma/2., color='black', ls='dashed')
#Note v_max_values multiplied by 100 to get units in cm/s
vmaximage = vmaxax.imshow(100*v_max_values, extent=[tmin, tmax, 0, gamma], origin='lower', aspect='auto', cmap='viridis')
plt.colorbar(vmaximage)

thetagradax = fig1.add_subplot(122)
thetagradax.set_title(r'$\frac{\partial \Theta}{\partial \tilde{y}}|^{sheet} \; (K/m)$', fontsize=18)
thetagradax.set_xlabel(r'$Time \; (T)$', fontsize=16)
thetagradax.set_ylabel(r'$\alpha \; (rad)$', fontsize=16)
thetagradax.axhline(y=gamma/2., color='black', ls='dashed')
thetagradimage = thetagradax.imshow(theta_gradients, extent=[tmin, tmax, 0, gamma], origin='lower', aspect='auto', cmap='magma')
plt.colorbar(thetagradimage)


fig2 = plt.figure()
heattransportax = fig2.add_subplot(111)
heattransportax.set_title(r'$\dot{Q}_{y}^{sheet} \; (W)$', fontsize=18)
heattransportax.set_xlabel(r'$Time \; (T)$', fontsize=16)
heattransportax.set_ylabel(r'$\alpha \; (rad)$', fontsize=16)
heattransportax.axhline(y=gamma/2., color='black', ls='dashed')
heattransportimage = heattransportax.contourf(heat_transport, extent=[tmin, tmax, 0, gamma], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(heattransportimage)

plt.show()
