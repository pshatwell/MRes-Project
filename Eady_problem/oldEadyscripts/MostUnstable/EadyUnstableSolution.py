import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from EadyUnstableinfoTEST import *

###############################################################################################

#Set up coordinates for plotting

xmin = 0
xmax = 6

ymin = -1
ymax = 1

zmin = -1
zmax = 1

xvalues = np.linspace(xmin, xmax, 50)
yvalues = np.linspace(ymin, ymax, 50)
zvalues = np.linspace(zmin, zmax, 50)

time = 0

xlength = len(xvalues)
ylength = len(yvalues)
zlength = len(zvalues)

###############################################################################################

#Create empty matrix for streamfunction
phiprime_matrix = np.zeros((zlength,ylength,xlength))

#Evaluate streamfunction matrix values
for i in range(0, zlength, 1):
    for j in range(0, ylength, 1):
        for m in range(0, xlength, 1):
            phiprime_matrix[i,j,m] = phiprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=time)


#Create empty matrix for meridional velocity perturbation
vprime_matrix = np.zeros((zlength,ylength,xlength))

#Evaluate vprime matrix values
for i in range(0, zlength, 1):
    for j in range(0, ylength, 1):
        for m in range(0, xlength, 1):
            vprime_matrix[i,j,m] = vprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=time)


#Create empty matrix for vertical velocity perturbation
wprime_matrix = np.zeros((zlength,ylength,xlength))

#Evaluate wprime matrix values
for i in range(0, zlength, 1):
    for j in range(0, ylength, 1):
        for m in range(0, xlength, 1):
            wprime_matrix[i,j,m] = wprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=time)

#Check boundary conditions (should be zero at vertical boundaries)
#print 'wprime at z=-H is:', wprime_matrix[0,10]
#print 'wprime at z=H is:', wprime_matrix[-1,10]


#Create empty matrix for potential temperature perturbation
thetaprime_matrix = np.zeros((zlength,ylength,xlength))

#Evaluate thetaprime matrix values:
for i in range(0, zlength, 1):
    for j in range(0, ylength, 1):
        for m in range(0, xlength, 1):
            thetaprime_matrix[i,j,m] = thetaprime(x=xvalues[m], y=yvalues[j], z=zvalues[i], t=time)


#Create matrix for background potential temperature distribution
theta_matrix = np.zeros((zlength,ylength))

thetayvalues = np.linspace(ymin,ymax,ylength)
thetazvalues = np.linspace(zmin,zmax,zlength)

for i in range(0, zlength, 1):
    for j in range(0, ylength, 1):
        theta_matrix[i,j] = theta(y=thetayvalues[j], z=thetazvalues[i])

###############################################################################################

#Plot the solution in x-z plane

fig = plt.figure()
fig.suptitle('x-z- plane')
plt.set_cmap('inferno')
ax1 = fig.add_subplot(221)
ax1.set_title('Streamfunction')
ax1.set_xlabel('x (L)')
ax1.set_ylabel('z (H)')
phicontour = ax1.contourf(phiprime_matrix[:,0,:], origin='lower', extent=[xmin, xmax, zmin, zmax], aspect='auto')
plt.colorbar(phicontour)

ax2 = fig.add_subplot(222)
ax2.set_title('v')
ax2.set_xlabel('x (L)')
ax2.set_ylabel('z (H)')
vcontour = ax2.contourf(vprime_matrix[:,0,:], origin='lower', extent=[xmin, xmax, zmin, zmax], aspect='auto')
plt.colorbar(vcontour)

ax3 = fig.add_subplot(223)
ax3.set_title('w')
ax3.set_xlabel('x (L)')
ax3.set_ylabel('z (H)')
wcontour = ax3.contourf(wprime_matrix[:,0,:], origin='lower', extent=[xmin, xmax, zmin, zmax], aspect='auto')
plt.colorbar(wcontour)

ax4 = fig.add_subplot(224)
ax4.set_title('theta')
ax4.set_xlabel('x (L)')
ax4.set_ylabel('z (H)')
thetaprimecontour = ax4.contourf(thetaprime_matrix[:,0,:], origin='lower', extent=[xmin, xmax, zmin, zmax], aspect='auto')
plt.colorbar(thetaprimecontour)

###############################################################################################

#Plot the solution in x-y plane

fig2 = plt.figure()
fig2.suptitle('x-y plane')
plt.set_cmap('inferno')
ax5 = fig2.add_subplot(221)
ax5.set_title('Streamfunction')
ax5.set_xlabel('x (L)')
ax5.set_ylabel('y (L)')
phicontour = ax5.contourf(phiprime_matrix[0,:,:], origin='lower', extent=[xmin, xmax, ymin, ymax], aspect='auto')
plt.colorbar(phicontour)

ax6 = fig2.add_subplot(222)
ax6.set_title('v')
ax6.set_xlabel('x (L)')
ax6.set_ylabel('y (L)')
vcontour = ax6.contourf(vprime_matrix[0,:,:], origin='lower', extent=[xmin, xmax, ymin, ymax], aspect='auto')
plt.colorbar(vcontour)

ax7 = fig2.add_subplot(223)
ax7.set_title('w')
ax7.set_xlabel('x (L)')
ax7.set_ylabel('y (L)')
wcontour = ax7.contourf(wprime_matrix[0,:,:], origin='lower', extent=[xmin, xmax, ymin, ymax], aspect='auto')
plt.colorbar(wcontour)

ax8 = fig2.add_subplot(224)
ax8.set_title('theta')
ax8.set_xlabel('x (L)')
ax8.set_ylabel('y (L)')
thetaprimecontour = ax8.contourf(thetaprime_matrix[0,:,:], origin='lower', extent=[xmin, xmax, ymin, ymax], aspect='auto')
plt.colorbar(thetaprimecontour)

plt.show()
#plt.close()

#Save figure as pdf file
#plt.savefig('figures/EadyUnstableSolution.pdf')
