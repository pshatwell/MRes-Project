
'''Script to plot solution for Eady model.'''

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la

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
print 'sigma is:', sigma

print '\ninitial uprime amplitude is:', 100*l
print 'initial vprime amplitude is:', 100*k
print 'initial wprime amplitude is:', -100*(f0/(N**2))*Lambda*k

###############################################################################################

def main(t):

    npts = 35

    xmin = -8
    xmax = 8

    ymin = -1
    ymax = 1

    zmin = 0
    zmax = 1

    xvalues = np.linspace(xmin, xmax, npts)
    yvalues = np.linspace(ymin, ymax, npts)
    zvalues = np.linspace(zmin, zmax, npts)

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

    print 'type of psiprime_matrix is:', type(psiprime_matrix)
    print 'shape of psiprime_matrix is:', psiprime_matrix.shape

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


    #Make array for full perturbation velocity field
    velocityfield = np.zeros((zlength, ylength, xlength, 3))

    print 'shape of velocity field is:', velocityfield.shape

    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                velocityfield[i,j,m,0] = uprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time) + umeanflow(z=zvalues[i]) - c.real*(T/L) #Subtracting phase speed, so we're in wave frame
                velocityfield[i,j,m,1] = vprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)
                velocityfield[i,j,m,2] = wprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for thetaprime
    thetaprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate thetaprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                thetaprime_matrix[i,j,m] = thetaprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

###############################################################################################

    #PLOTTING

    #Define index slices for plotting

    zslice = int(round(zlength/2.))
    yslice = int(round(ylength/2.))
    xslice = int(round(xlength/2.))


    #Plot streamfunction perturbation

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('psiprime')
    ax1 = fig1.add_subplot(311)
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    #Plotting in x-y plane at z=H/2
    xy_contour = ax1.contourf(psiprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour)

    ax2 = fig1.add_subplot(312)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour = ax2.contourf(psiprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour)

    ax3 = fig1.add_subplot(313)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour = ax3.contourf(psiprime_matrix[:,:,xslice], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
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
    xy_contour_u = ax4.contourf(uprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_u)

    ax5 = fig2.add_subplot(312)
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_u = ax5.contourf(uprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_u)

    ax6 = fig2.add_subplot(313)
    ax6.set_xlabel('y (L)')
    ax6.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_u = ax6.contourf(uprime_matrix[:,:,xslice], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
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
    xy_contour_v = ax7.contourf(vprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_v)

    ax8 = fig3.add_subplot(312)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_v = ax8.contourf(vprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_v)

    ax9 = fig3.add_subplot(313)
    ax9.set_xlabel('y (L)')
    ax9.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_v = ax9.contourf(vprime_matrix[:,:,xslice], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
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
    xy_contour_w = ax10.contourf(wprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_w)

    ax11 = fig4.add_subplot(312)
    ax11.set_xlabel('x (L)')
    ax11.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_w = ax11.contourf(wprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_w)

    ax12 = fig4.add_subplot(313)
    ax12.set_xlabel('y (L)')
    ax12.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_w = ax12.contourf(wprime_matrix[:,:,xslice], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
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
    xy_contour_theta = ax13.contourf(thetaprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_theta)

    ax14 = fig5.add_subplot(312)
    ax14.set_xlabel('x (L)')
    ax14.set_ylabel('z (H)')
    #Plotting in x-z plane
    xz_contour_theta = ax14.contourf(thetaprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_theta)

    ax15 = fig5.add_subplot(313)
    ax15.set_xlabel('y (L)')
    ax15.set_ylabel('z (H)')
    #Plotting in y-z plane
    yz_contour_theta = ax15.contourf(thetaprime_matrix[:,:,xslice], origin='lower', aspect='auto', extent = [ymin,ymax,zmin,zmax])
    plt.colorbar(yz_contour_theta)

###############################################################################################

    #Compare vertical structure of different fields
    fig6 = plt.figure()
    fig6.suptitle('Eady solutions', fontsize=17)
    plt.set_cmap('inferno')
    ax16 = fig6.add_subplot(221)
    ax16.set_title('Streamfunction', fontsize=16)
    #ax16.set_xlabel('x (L)')
    ax16.set_ylabel('z (H)', fontsize=16)
    ax16.contourf(psiprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

    ax17 = fig6.add_subplot(222)
    ax17.set_title('v perturbation', fontsize=16)
    #ax17.set_xlabel('x (L)')
    ax17.set_ylabel('z (H)', fontsize=16)
    ax17.contourf(vprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

    ax18 = fig6.add_subplot(223)
    ax18.set_title('w perturbation', fontsize=16)
    ax18.set_xlabel('x (L)', fontsize=16)
    ax18.set_ylabel('z (H)', fontsize=16)
    ax18.contourf(wprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

    ax19 = fig6.add_subplot(224)
    ax19.set_title('theta perturbation', fontsize=16)
    ax19.set_xlabel('x (L)', fontsize=16)
    ax19.set_ylabel('z (H)', fontsize=16)
    ax19.contourf(thetaprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])

###############################################################################################

    #Plot the vector velocity field

    #Define grids for quiver plots
    X, Y = np.meshgrid(xvalues, yvalues)
    X1, Z1 = np.meshgrid(xvalues, zvalues)
    Y2, Z2 = np.meshgrid(yvalues, zvalues)

    #Plot velocity field in each plane
    velfieldfig1 = plt.figure()
    velfieldax1 = velfieldfig1.add_subplot(111)
    velfieldax1.set_xlabel('x (L)', fontsize='15')
    velfieldax1.set_ylabel('y (L)', fontsize='15')
    Q1 = velfieldax1.quiver(X,Y, velocityfield[zslice,:,:,0], velocityfield[zslice,:,:,1], headlength=3, headaxislength=3)

    velfieldfig2 = plt.figure()
    velfieldax2 = velfieldfig2.add_subplot(111)
    velfieldax2.set_xlabel('x (L)', fontsize='15')
    velfieldax2.set_ylabel('z (H)', fontsize='15')
    Q2 = velfieldax2.quiver(X1,Z1, velocityfield[:,yslice,:,0], velocityfield[:,yslice,:,2], headlength=3, headaxislength=3)

    #Note w is multiplied by 10 here in order to be seen
    velfieldfig3 = plt.figure()
    velfieldax3 = velfieldfig3.add_subplot(111)
    velfieldax3.set_xlabel('y (L)', fontsize='15')
    velfieldax3.set_ylabel('z (H)', fontsize='15')
    Q3 = velfieldax3.quiver(Y2,Z2, velocityfield[:,:,xslice,1], 10*velocityfield[:,:,xslice,2], headlength=3, headaxislength=3)

###############################################################################################

    #Show all the plots

    plt.show()

###############################################################################################

#Run the programme, choosing time (in T) to evaluate fields

main(t=2)
