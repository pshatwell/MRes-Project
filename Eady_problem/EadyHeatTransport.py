
'''Script to plot structure of meridional and vertical heat transports
in the Eady Model.'''

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
from matplotlib import rc
rc('text', usetex=True)

from Eadyinfo_oc import *

###############################################################################################

#Set resolution
npts = 20

#Define coordinates
xmin = -(np.pi)/(k*L)
xmax = (np.pi)/(k*L)
xslice = 10

ymin = -1.
ymax = 1.

zmin = 0.
zmax = 1.

xvalues = np.linspace(xmin, xmax, npts)
yvalues = np.linspace(ymin, ymax, npts)
zvalues = np.linspace(zmin, zmax, npts)

xlength = len(xvalues)
ylength = len(yvalues)
zlength = len(zvalues)

###############################################################################################

'''HOW DOES THE ZONAL MEAN STRUCTURE DEPEND ON THE TIMESTEP?'''

timestep = 20 #Evaluate Eady fields at this time

#Calculate the poleward and vertical heat transports for Eady model

vprime_matrix = np.zeros((zlength,ylength,xlength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        for m in range(0,xlength,1):
            vprime_matrix[i,j,m] = (L/T)*vprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=timestep)


wprime_matrix = np.zeros((zlength,ylength,xlength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        for m in range(0,xlength,1):
            wprime_matrix[i,j,m] = (H/T)*wprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=timestep)


thetaprime_matrix = np.zeros((zlength,ylength,xlength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        for m in range(0,xlength,1):
            thetaprime_matrix[i,j,m] = thetaprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=timestep)


poleward_heat_transport = rho_0*c_p*np.multiply(vprime_matrix,thetaprime_matrix)
vertical_heat_transport = rho_0*c_p*np.multiply(wprime_matrix,thetaprime_matrix)

###############################################################################################

#Calculate zonal mean of poleward heat transport
zonal_mean_transport_p = np.zeros((zlength,ylength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zonal_mean_transport_p[i,j] = (1./xlength)*np.sum(poleward_heat_transport[i,j,m] for m in range(0,xlength,1))


#Calculate zonal mean of vertical heat transport
zonal_mean_transport_v = np.zeros((zlength,ylength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zonal_mean_transport_v[i,j] = (1./xlength)*np.sum(vertical_heat_transport[i,j,m] for m in range(0,xlength,1))


#Create zonal mean heat transport vector in meridional plane
zmean_transport_vector = np.zeros((zlength,ylength,2))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zmean_transport_vector[i,j,0] = zonal_mean_transport_p[i,j]
        zmean_transport_vector[i,j,1] = zonal_mean_transport_v[i,j]

transport_vector = np.zeros((zlength,ylength,2))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        transport_vector[i,j,0] = poleward_heat_transport[i,j,xslice]
        transport_vector[i,j,1] = vertical_heat_transport[i,j,xslice]


halftheta_matrix = np.zeros((zlength,ylength))

#Evaluate values for half-slope surfaces
for i in range(0, zlength, 1):
    for j in range(0, ylength, 1):
        halftheta_matrix[i,j] = halftheta(y=yvalues[j], z=zvalues[i])

###############################################################################################

#Calculate zonal means of perturbation fields

#zonal mean of vprime
zmean_vprime = np.zeros((zlength,ylength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zmean_vprime[i,j] = (1./xlength)*np.sum(vprime_matrix[i,j,m] for m in range(0,xlength,1))

#zonal mean of wprime
zmean_wprime = np.zeros((zlength,ylength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zmean_wprime[i,j] = (1./xlength)*np.sum(wprime_matrix[i,j,m] for m in range(0,xlength,1))

zmean_vector = np.zeros((zlength,ylength,2))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zmean_vector[i,j,0] = zmean_vprime[i,j]
        zmean_vector[i,j,1] = zmean_wprime[i,j]

#zonal mean of thetaprime
zmean_thetaprime = np.zeros((zlength,ylength))

for i in range(0,zlength,1):
    for j in range(0,ylength,1):
        zmean_thetaprime[i,j] = (1./xlength)*np.sum(thetaprime_matrix[i,j,m] for m in range(0,xlength,1))

###############################################################################################

#Slices for plotting give midpoint indices of x,y,z values
xslice = int(round(xlength/2.))
yslice = int(round(ylength/2.))
zslice = int(round(zlength/2.))


#Plot the poleward heat transport
fig = plt.figure()
fig.suptitle(r"Poleward Heat Transport $v'\theta' \; (Wm^{-2})$", fontsize=18)
pheatax_xy = fig.add_subplot(311)
pheatax_xy.set_xlabel(r'$x \; (L)$', fontsize=16)
pheatax_xy.set_ylabel(r'$y \; (L)$', fontsize=16)
pheatax_xy_contour = pheatax_xy.contourf(poleward_heat_transport[zslice,:,:], extent=[xmin,xmax,ymin,ymax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(pheatax_xy_contour)

pheatax_xz = fig.add_subplot(312)
pheatax_xz.set_xlabel(r'$x \; (L)$', fontsize=16)
pheatax_xz.set_ylabel(r'$z \; (H)$', fontsize=16)
pheatax_xz_contour = pheatax_xz.contourf(poleward_heat_transport[:,yslice,:], extent=[xmin,xmax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(pheatax_xz_contour)

pheatax_yz = fig.add_subplot(313)
pheatax_yz.set_xlabel(r'$y \; (L)$', fontsize=16)
pheatax_yz.set_ylabel(r'$z \; (H)$', fontsize=16)
pheatax_yz_contour = pheatax_yz.contourf(poleward_heat_transport[:,:,xslice], extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(pheatax_yz_contour)


#Plot the vertical heat transport
fig2 = plt.figure()
fig2.suptitle(r"Vertical Heat Transport $w'\theta' \; (Wm^{-2})$", fontsize=18)
vheatax_xy = fig2.add_subplot(311)
vheatax_xy.set_xlabel(r'$x \; (L)$', fontsize=16)
vheatax_xy.set_ylabel(r'$y \; (L)$', fontsize=16)
vheatax_xy_contour = vheatax_xy.contourf(vertical_heat_transport[zslice,:,:], extent=[xmin,xmax,ymin,ymax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(vheatax_xy_contour)

vheatax_xz = fig2.add_subplot(312)
vheatax_xz.set_xlabel(r'$x \; (L)$', fontsize=16)
vheatax_xz.set_ylabel(r'$z \; (H)$', fontsize=16)
vheatax_xz_contour = vheatax_xz.contourf(vertical_heat_transport[:,yslice,:], extent=[xmin,xmax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(vheatax_xz_contour)

vheatax_yz = fig2.add_subplot(313)
vheatax_yz.set_xlabel(r'$y \; (L)$', fontsize=16)
vheatax_yz.set_ylabel(r'$z \; (H)$', fontsize=16)
vheatax_yz_contour = vheatax_yz.contourf(vertical_heat_transport[:,:,xslice], extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(vheatax_yz_contour)

###############################################################################################

#Plot the zonal mean poleward heat transport
#Independent of height and peaks in middle of channel as expected
fig3 = plt.figure()
fig3.suptitle(r'Zonal-mean heat transport $(Wm^{-2})$', fontsize=18)
zmeanax = fig3.add_subplot(211)
zmeanax.set_title(r"$\bar{v'\theta'}$", fontsize=16)
#zmeanax.set_xlabel(r'$y \; (L)$', fontsize=16)
zmeanax.set_ylabel(r'$z \; (H)$', fontsize=16)
zmeancontour = zmeanax.contourf(zonal_mean_transport_p, extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(zmeancontour)

zmeanax2 = fig3.add_subplot(212)
zmeanax2.set_title(r"$\bar{w'\theta'}$", fontsize=16)
zmeanax2.set_xlabel(r'$y \; (L)$', fontsize=16)
zmeanax2.set_ylabel(r'$z \; (H)$', fontsize=16)
zmeancontour2 = zmeanax2.contourf(zonal_mean_transport_v, extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(zmeancontour2)

'''
#Plot the zonal mean perturbation fields
fig4 = plt.figure()
fig4.suptitle(r'Zonal mean perturbation fields', fontsize=18)
vprimeax = fig4.add_subplot(311)
vprimeax.set_title(r"$\bar{v'}$")
#vprimeax.set_xlabel(r'$y \; (L)$', fontsize=16)
vprimeax.set_ylabel(r'$z \; (H)$', fontsize=16)
vprimecontour = vprimeax.contourf(zmean_vprime, extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(vprimecontour)

wprimeax = fig4.add_subplot(312)
wprimeax.set_title(r"$\bar{w'}$")
#wprimeax.set_xlabel(r'$y \; (L)$', fontsize=16)
wprimeax.set_ylabel(r'$z \; (H)$', fontsize=16)
wprimecontour = wprimeax.contourf(zmean_wprime, extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(wprimecontour)

thetaprimeax = fig4.add_subplot(313)
thetaprimeax.set_title(r"$\bar{\theta'}$")
thetaprimeax.set_xlabel(r'$y \; (L)$', fontsize=16)
thetaprimeax.set_ylabel(r'$z \; (H)$', fontsize=16)
thetaprimecontour = thetaprimeax.contourf(zmean_thetaprime, extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='inferno')
plt.colorbar(thetaprimecontour)
'''

#Define grid for quiver plot
Yp, Zp = np.meshgrid(yvalues, zvalues)

heatfig = plt.figure()
heatax = heatfig.add_subplot(111)
heatax.set_title(r'Zonal-mean heat transport',fontsize=18)
heatax.set_xlabel(r'$y \; (L)$', fontsize=16)
heatax.set_ylabel(r'$z \; (H)$', fontsize=16)
meantheta = heatax.contourf(halftheta_matrix, extent=[ymin,ymax,zmin,zmax], origin='lower', aspect='auto', cmap='viridis')
Q = heatax.quiver(Yp, Zp, zmean_transport_vector[:,:,0], zmean_transport_vector[:,:,1], units='xy')
#Q = heatax.quiver(Yp, Zp, transport_vector[:,:,0], transport_vector[:,:,1], units='xy')
heatax.quiverkey(Q, 0.8, 0.92, -2, r'$2 \; Wm^{-2}$', labelpos='N', coordinates = 'figure')

plt.show()
