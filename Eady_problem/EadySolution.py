
'''Script to plot solution for Eady model.'''


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as la

from Eadyinfo_oc import *


print 'H is', H
print 'L is', L
print 'T is', T

print 'k is', k
print 'l is', l
print 'mu is', mu

print 'c is:', c
print 'sigma is:', sigma

###############################################################################################

def main(t):

    xmin = -3
    xmax = 3

    ymin = -1
    ymax = 1

    zmin = 0
    zmax = 1

    xvalues = np.linspace(xmin, xmax, 50)
    yvalues = np.linspace(ymin, ymax, 50)
    zvalues = np.linspace(zmin, zmax, 50)

    time = t

    xlength = len(xvalues)
    ylength = len(yvalues)
    zlength = len(zvalues)

###############################################################################################

    #Create arrays to display solution graphically

    #Create empty matrix for streamfunction perturbation
    psiprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                psiprime_matrix[i,j,m] = psiprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for uprime
    uprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                uprime_matrix[i,j,m] = uprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for vprime
    vprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                vprime_matrix[i,j,m] = vprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for wprime
    wprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                wprime_matrix[i,j,m] = wprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for thetaprime
    thetaprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate thetaprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                thetaprime_matrix[i,j,m] = thetaprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

###############################################################################################

    #Plot the streamfunction and vertical velocity perturbations

    #Plot streamfunction perturbation

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('psiprime')
    ax1 = fig1.add_subplot(311)
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    #Plotting in x-y plane at z=H/2
    xy_contour = ax1.contourf(psiprime_matrix[24,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour)

    ax2 = fig1.add_subplot(312)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour = ax2.contourf(psiprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour)

    ax3 = fig1.add_subplot(313)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour = ax3.contourf(psiprime_matrix[:,:,0], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour)

###############################################################################################

    #Plot zonal velocity perturbation

    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('uprime')
    ax4 = fig2.add_subplot(311)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    #Plotting in x-y plane at z=H/2
    xy_contour_u = ax4.contourf(uprime_matrix[24,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_u)

    ax5 = fig2.add_subplot(312)
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_u = ax5.contourf(uprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_u)

    ax6 = fig2.add_subplot(313)
    ax6.set_xlabel('y (L)')
    ax6.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_u = ax6.contourf(uprime_matrix[:,:,0], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour_u)

###############################################################################################

    #Plot meridional velocity perturbation

    fig3 = plt.figure()
    plt.set_cmap('inferno')
    fig3.suptitle('vprime')
    ax7 = fig3.add_subplot(311)
    ax7.set_xlabel('x (L)')
    ax7.set_ylabel('y (L)')
    #Plotting in x-y plane at z=H/2
    xy_contour_v = ax7.contourf(vprime_matrix[24,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_v)

    ax8 = fig3.add_subplot(312)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_v = ax8.contourf(vprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_v)

    ax9 = fig3.add_subplot(313)
    ax9.set_xlabel('y (L)')
    ax9.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_v = ax9.contourf(vprime_matrix[:,:,0], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour_v)

###############################################################################################

    #Plot vertical velocity perturbation

    fig4 = plt.figure()
    plt.set_cmap('inferno')
    fig4.suptitle('wprime')
    ax10 = fig4.add_subplot(311)
    ax10.set_xlabel('x (L)')
    ax10.set_ylabel('y (L)')
    #Plotting in x-y plane at z=H/2
    xy_contour_w = ax10.contourf(wprime_matrix[24,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_w)

    ax11 = fig4.add_subplot(312)
    ax11.set_xlabel('x (L)')
    ax11.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_w = ax11.contourf(wprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_w)

    ax12 = fig4.add_subplot(313)
    ax12.set_xlabel('y (L)')
    ax12.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_w = ax12.contourf(wprime_matrix[:,:,0], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour_w)

###############################################################################################

    #Plot potential temperature perturbation

    fig5 = plt.figure()
    plt.set_cmap('inferno')
    fig5.suptitle('thetaprime')
    ax13 = fig5.add_subplot(311)
    ax13.set_xlabel('x (L)')
    ax13.set_ylabel('y (L)')
    #Plotting in x-y plane at z=H/2
    xy_contour_theta = ax13.contourf(thetaprime_matrix[24,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_theta)

    ax14 = fig5.add_subplot(312)
    ax14.set_xlabel('x (L)')
    ax14.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_theta = ax14.contourf(thetaprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_theta)

    ax15 = fig5.add_subplot(313)
    ax15.set_xlabel('y (L)')
    ax15.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_theta = ax15.contourf(thetaprime_matrix[:,:,0], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour_theta)

###############################################################################################

    #Compare vertical structure of different fields
    fig6 = plt.figure()
    fig6.suptitle('Eady solutions', fontsize=15)
    plt.set_cmap('inferno')
    ax16 = fig6.add_subplot(221)
    ax16.set_title('Streamfunction', fontsize=14)
    #ax16.set_xlabel('x (L)')
    ax16.set_ylabel('z (H)')
    ax16.contourf(psiprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

    ax17 = fig6.add_subplot(222)
    ax17.set_title('v perturbation', fontsize=14)
    #ax17.set_xlabel('x (L)')
    ax17.set_ylabel('z (H)')
    ax17.contourf(vprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

    ax18 = fig6.add_subplot(223)
    ax18.set_title('w perturbation', fontsize=14)
    ax18.set_xlabel('x (L)')
    ax18.set_ylabel('z (H)')
    ax18.contourf(wprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

    ax19 = fig6.add_subplot(224)
    ax19.set_title('theta perturbation', fontsize=14)
    ax19.set_xlabel('x (L)')
    ax19.set_ylabel('z (H)')
    ax19.contourf(thetaprime_matrix[:,0,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])


    plt.show()

###############################################################################################


main(t=0)
