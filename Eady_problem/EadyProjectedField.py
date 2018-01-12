
'''Script to project Eady Model velocity field within a sheet slanted upwards
 towards the pole. A slope equal to half that of the mean isentropes is meant to
 represent the main axis of the instability.

 An attempt is also made to plot the heat transport within the sheet.'''

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
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

print '\ninitial uprime amplitude is:', l
print 'initial vprime amplitude is:', k
print 'initial wprime amplitude is:', -(f0/(N**2))*Lambda*k

print '\ntime (s) when vprime amplitude is equal to U_0 is:', (1./sigma)*np.log(U_0/(k))
print 'in periods of T, that is:', ((1./sigma)*np.log(U_0/(k)))/T

###############################################################################################

#Define coordinates

npts = 30

xmin = -(2*np.pi)/(k*L)
xmax = (2*np.pi)/(k*L)

ymin = -1
ymax = 1

zmin = 0
zmax = 1

tmin = 0
tmax = 101

xvalues = np.linspace(xmin, xmax, npts)
yvalues = np.linspace(ymin, ymax, npts)
zvalues = np.linspace(zmin, zmax, npts)

tvalues = np.arange(tmin, tmax, 1)

xlength = len(xvalues)
ylength = len(yvalues)
zlength = len(zvalues)
tlength = len(tvalues)

###############################################################################################

#Create the projected velocity field for Green motion

#Normal to surface at half the slope of theta surface
#Maximum energy release for baroclinic instability typically within this half-slope surface
normal_i = np.array((0, dthetady*L, 2*dthetadz*H)) #i for instability
#normal_i = np.array((0, dthetady*L, dthetadz*H))
normal_ihat = normal_i/(la.norm(normal_i))

#Slope (angle) of mean isentropic surfaces
gamma = np.arctan((dthetady*L)/(dthetadz*H))

timestep = 20 #Choose from 0-100 (in units of T)

print 'gamma is:', gamma
print 'gamma/2 is:', gamma/2.

zconst = 0.5 #height of plane

planezvalues = zconst - yvalues*(normal_i[1]/normal_i[2]) #Equation for plane
planezlength = len(planezvalues)

planeyvalues = (zconst - planezvalues)*(normal_i[2]/normal_i[1])
planeylength = len(planeyvalues)


#Create array for velocity field in half-isentrope plane
velocityfield_p = np.zeros((planeylength, xlength, tlength, 3))

#Evaluate 3d velocity values on the plane
for i in range(0, planeylength, 1):
    for j in range(0, xlength, 1):
        for m in range(0, tlength, 1):
            #velocityfield_p[i,j,m,0] = uprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m]) + umeanflow(z=planezvalues[i])
            velocityfield_p[i,j,m,0] = (L/T)*uprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m]) + (L/T)*umeanflow(z=planezvalues[i]) - c.real #Subtracting phase speed, so we're in wave frame
            velocityfield_p[i,j,m,1] = (L/T)*vprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m])
            velocityfield_p[i,j,m,2] = (H/T)*wprime(x=xvalues[j], y = planeyvalues[i], z=planezvalues[i], t=tvalues[m])


#Project the 3d velocity vectors from velocityfield_p onto 2d plane
velocityfield_p2 = np.zeros((planeylength, xlength, tlength, 2))

#These define the vectors that span the plane
xtilde = np.array((1,0,0))
ytilde = np.cross(normal_ihat,xtilde)

for i in range(0, planeylength, 1):
    for j in range(0, xlength, 1):
        for m in range(0, tlength, 1):
            velocityfield_p2[i,j,m,0] = np.dot(velocityfield_p[i,j,m], xtilde)
            velocityfield_p2[i,j,m,1] = np.dot(velocityfield_p[i,j,m], ytilde)

###############################################################################################

#Create projected mean theta field
thetafield_p2 = np.zeros((planeylength, xlength))

for i in range(0, planeylength, 1):
    for j in range(0, xlength, 1):
        thetafield_p2[i,j] = theta(y = planeyvalues[i], z = planezvalues[i])

#Create projected theta perturbation field
thetaprimefield_p2 = np.zeros((planeylength, xlength, tlength))

for i in range(0, planeylength, 1):
    for j in range(0, xlength, 1):
        for m in range(0, tlength, 1):
            thetaprimefield_p2[i,j,m] = thetaprime(x=xvalues[j], y=planeyvalues[i], z=planezvalues[i], t=tvalues[m])


#Calculate heat transport within the sheet
vtildetransport = rho_0*c_p*np.multiply(velocityfield_p2[:,:,timestep,1], thetaprimefield_p2[:,:,timestep])

###############################################################################################

#Plot the vector velocity field

#Define index slices for plotting
planexslice = int(round(xlength/2.))
planeyslice = int(round(planeylength/2.))

#Define grid for quiver plot
Xp, Yp = np.meshgrid(xvalues, planeyvalues)

#Velocity magnitudes
vels = np.hypot(100*velocityfield_p2[:,:,timestep,0], 100*velocityfield_p2[:,:,timestep,1]) #in cm/s


#Plot the projected velocity field
velfieldfig = plt.figure()
velfieldax = velfieldfig.add_subplot(111)
velfieldax.set_title(r'$\alpha = \gamma/2, \quad T = %s, \quad z = %s$' % (str(timestep), str(zconst)), fontsize='16')
velfieldax.set_xlabel(r'$\tilde{x} \; (L)$', fontsize='16')
velfieldax.set_ylabel(r'$\tilde{y} \; (L)$', fontsize='16')
velfieldax.set_ylim([-1.1,1.1])
Q = velfieldax.quiver(Xp, Yp, 100*velocityfield_p2[:,:,timestep,0], 100*velocityfield_p2[:,:,timestep,1], vels, units='xy', scale=1, cmap='magma') #Note values are multiplied by 100 (cm/s)
velfieldax.quiverkey(Q, 0.8, 0.93, 0.4, '0.4 cm/s', labelpos='E', coordinates = 'figure')
plt.colorbar(Q)

#Plot the projected velocity field with background mean isentropes
velfieldfig2 = plt.figure()
velfieldax2 = velfieldfig2.add_subplot(111)
velfieldax2.set_title(r'$\alpha = \gamma/2, \quad T = %s, \quad z = %s$' % (str(timestep), str(zconst)), fontsize='16')
velfieldax2.set_xlabel(r'$\tilde{x} \; (L)$', fontsize='16')
velfieldax2.set_ylabel(r'$\tilde{y} \; (L)$', fontsize='16')
velfieldax2.set_ylim([-1.1,1.1])
thetacontour = velfieldax2.contourf(thetafield_p2, extent=[xmin,xmax,ymin,ymax], origin='lower', aspect='auto', cmap='viridis')
plt.colorbar(thetacontour)
Q2 = velfieldax2.quiver(Xp, Yp, 100*velocityfield_p2[:,:,timestep,0], 100*velocityfield_p2[:,:,timestep,1], units='xy', scale=1, cmap='magma') #Note values are multiplied by 100 (cm/s)
velfieldax2.quiverkey(Q2, 0.64, 0.93, 0.4, '0.4 cm/s', labelpos='E', coordinates = 'figure')

###############################################################################################

#Plot the heat transport within the sheet

heattransportfig = plt.figure()
heattransportax = heattransportfig.add_subplot(111)
heattransportax.set_title(r"$\tilde{v'}\tilde{\theta'} \quad (Wm^{-2}) \qquad (\alpha = \gamma/2, \quad T = %s, \quad z = %s)$" % (str(timestep), str(zconst)), fontsize='16')
heattransportax.set_xlabel(r'$\tilde{x} \; (L)$', fontsize='16')
heattransportax.set_ylabel(r'$\tilde{y} \; (L)$', fontsize='16')
heattransportcontour = heattransportax.contourf(vtildetransport, extent=[xmin,xmax,ymin,ymax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(heattransportcontour)

###############################################################################################
'''
vprimetilde_list = []

#Note velocity values are multiplied by (L/T)

#Calculating zonal mean of absolute value of vprime along the centre of the plane
for i in range(0, tlength, 1):
    vprimetilde = 100*(L/T)*(1./xlength)*np.sum(np.abs(velocityfield_p2[planeyslice,j,i,1]) for j in range(xlength)) #in cm/s
    vprimetilde_list.append(vprimetilde)

vprime_amp_values = 100*k*np.exp(sigma*T*tvalues) #in cm/s

#Plotting growth of vprime with time
vprimetildefig = plt.figure()
vprimetildeax = vprimetildefig.add_subplot(111)
vprimetildeax.set_xlabel(r'Time (T)', fontsize='16')
vprimetildeax.set_ylabel( r'(cm/s)', fontsize='16')
vprimetildeax.plot(tvalues, vprime_amp_values, label=r"$|v'|$ amplitude")
vprimetildeax.plot(tvalues, vprimetilde_list, label=r"$\bar{|\tilde{v'}|}$ zonal mean")
vprimetildeax.axhline(y=5, color='black', ls='dashed')
vprimetildeax.axhline(y=10, color='black', ls='dashed')
vprimetildeax.legend(loc='upper left')
vprimetildeax.grid()
'''

plt.show()
