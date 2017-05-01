import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as la

from Eadyinfo import *

print 'H is', H
print 'L is', L
print 'T is', T
#print 'e-folding time (T) is 1/(sigma*T):', 1./(sigma*T)

###############################################################################################

def main(t):

    xmin = -3
    xmax = 3

    ymin = -1
    ymax = 1

    zmin = -1
    zmax = 1

    xvalues = np.linspace(xmin, xmax, 50)
    yvalues = np.linspace(ymin, ymax, 50)
    zvalues = np.linspace(zmin, zmax, 50)

    time = t

    xlength = len(xvalues)
    ylength = len(yvalues)
    zlength = len(zvalues)

###############################################################################################

    #Create empty matrix for streamfunction
    phi_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                phi_matrix[i,j,m] = streamfunction(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    print 'type of phi_matrix is:', type(phi_matrix)
    print 'shape of phi_matrix is:', phi_matrix.shape

    #Create empty matrix for wprime
    wprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                wprime_matrix[i,j,m] = wprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    print 'type of wprime_matrix is:', type(wprime_matrix)
    print 'shape of wprime_matrix is:', wprime_matrix.shape

###############################################################################################
    zheight = 0.5

    def velocity(s,t):
        x,y,z = s
        dsdt = [umeanflow(z=z) + uprime(x=x,y=y,z=z,t=t), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt

    tmin = 0
    tmax = 20

    t = np.linspace(tmin, tmax, 200)

    s0 = np.array((0,0,0))

    sol = odeint(velocity, s0, t)

###############################################################################################

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Streamfunction')
    ax1 = fig1.add_subplot(311)
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    #Plotting in x-y plane at z=0
    xy_contour = ax1.contourf(phi_matrix[24,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour)

    ax2 = fig1.add_subplot(312)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    #Plotting in x-z plane at y=0
    xz_contour = ax2.contourf(phi_matrix[:,24,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour)

    ax3 = fig1.add_subplot(313)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    yz_contour = ax3.contourf(phi_matrix[:,:,24], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour)

    fig2 = plt.figure()
    ax3 = fig2.add_subplot(111, projection='3d')
    ax3.plot(sol[:,0], sol[:,1], sol[:,2])

    plt.show()

main(t=0)
