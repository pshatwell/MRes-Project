
'''Script to demonstrate Green's idea of isentropic relative flow
using unstable solutions of the Eady model.'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import numpy.linalg as la

from Eadyinfo import *

print 'H is (m)', H
print 'L is (m)', L
print 'T is (s)', T

print 'k is (m^-1)', k
print 'l is (m^-1)', l
print 'mu is', mu

print 'c is (ms^-1)', c
print 'sigma is (s^-1)', sigma
if sigma != 0:
    print 'e-folding time (T) is:', 1./(sigma*T)

###############################################################################################

def main(start, stop):

    #Define velocity function for 3d parcel trajectories

    def velocity3d(s,t):
        x,y,z = s
        dsdt = [uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        #dsdt = [uprime(x=x,y=y,z=z,t=t), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt


###############################################################################################

    #Define timesteps for integration

    tmin = start #Should this necessarily be zero?
    tmax = stop

    t = np.linspace(tmin, tmax, 500)

###############################################################################################

    #Define initial positions of parcels

    #setting both 'cold' and 'warm' parcels off at the same height but different y
    #seems to give closest reproduction of Green picture
    zpos = 0.5
    xpos = 0
    xshift = np.pi/(k*L) #nondimensional shift of half a wavelength
    ypos = 0.85

    #Set of 5 'warm' parcels (positive y)
    s0_a = np.array((xpos, ypos, zpos ))
    s0_b = np.array((xpos+0.5, ypos, zpos))
    s0_c = np.array((xpos+1, ypos, zpos))
    s0_d = np.array((xpos+1.5, ypos, zpos))
    s0_e = np.array((xpos+2, ypos, zpos))

    #Set of 5 'cold' parcels (negative y)
    s0_f = np.array((xpos+xshift, -ypos, zpos))
    s0_g = np.array((xpos+xshift+0.5, -ypos, zpos))
    s0_h = np.array((xpos+xshift+1, -ypos, zpos))
    s0_i = np.array((xpos+xshift+1.5, -ypos, zpos))
    s0_j = np.array((xpos+xshift+2, -ypos, zpos))

###############################################################################################

    #Solve for parcel trajectories

    sol_a = odeint(velocity3d, s0_a, t)
    sol_b = odeint(velocity3d, s0_b, t)
    sol_c = odeint(velocity3d, s0_c, t)
    sol_d = odeint(velocity3d, s0_d, t)
    sol_e = odeint(velocity3d, s0_e, t)

    sol_f = odeint(velocity3d, s0_f, t)
    sol_g = odeint(velocity3d, s0_g, t)
    sol_h = odeint(velocity3d, s0_h, t)
    sol_i = odeint(velocity3d, s0_i, t)
    sol_j = odeint(velocity3d, s0_j, t)

    EFsolutions = [sol_a, sol_b, sol_c, sol_d, sol_e, sol_f, sol_g, sol_h, sol_i, sol_j]

###############################################################################################

    #Distance calculations between parcels a and b

    #Calculate initial parcel separation
    d_i = la.norm(s0_b - s0_a)
    print 'initial separation is:', d_i

    #Calculate final parcel separation
    d_f = la.norm(sol_b[-1] - sol_a[-1])
    print 'final separation is:', d_f

    distances = []
    for i in range(len(t)):
        distances.append(la.norm(sol_b[i] - sol_a[i]))

###############################################################################################

    #Transform to wave frame by shifting x-coordinates
    shift = np.zeros_like(sol_a)
    shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

    #Define new x-shifted trajectories
    rel_sol_a = sol_a + shift
    rel_sol_b = sol_b + shift
    rel_sol_c = sol_c + shift
    rel_sol_d = sol_d + shift
    rel_sol_e = sol_e + shift

    rel_sol_f = sol_f + shift
    rel_sol_g = sol_g + shift
    rel_sol_h = sol_h + shift
    rel_sol_i = sol_i + shift
    rel_sol_j = sol_j + shift

    WFsolutions = [rel_sol_a, rel_sol_b, rel_sol_c, rel_sol_d, rel_sol_e, rel_sol_f, rel_sol_g, rel_sol_h, rel_sol_i, rel_sol_j]
    warmWFsolutions = [rel_sol_a, rel_sol_b, rel_sol_c, rel_sol_d, rel_sol_e]
    coldWFsolutions = [rel_sol_f, rel_sol_g, rel_sol_h, rel_sol_i, rel_sol_j]

###############################################################################################

    #Create matrix for background potential temperature distribution
    theta_matrix = np.zeros((50,50))

    ymin = min(sol_a[:,1])
    ymax = max(sol_a[:,1])

    zmin = min(sol_a[:,2])
    zmax = max(sol_a[:,2])

    thetayvalues = np.linspace(ymin,ymax,50)
    thetazvalues = np.linspace(zmin,zmax,50)

    for i in range(0, 50, 1):
        for j in range(0, 50, 1):
            theta_matrix[i,j] = theta(y=thetayvalues[j], z=thetazvalues[i])

    print 'dthetady is:', dthetady*L
    print 'dthetadz is:', dthetadz*H

    gradient = (dthetadz*H)/(dthetady*L)
    print 'gradient is:', gradient

###############################################################################################

    #Define normal to theta surfaces

    normal = np.array((0, dthetady*L, dthetadz*H))
    normalhat = normal/(la.norm(normal))

    print 'normal to theta surfaces is:', normalhat

    sol_a_p = np.zeros_like(rel_sol_a)
    sol_b_p = np.zeros_like(rel_sol_a)
    sol_c_p = np.zeros_like(rel_sol_a)
    sol_d_p = np.zeros_like(rel_sol_a)
    sol_e_p = np.zeros_like(rel_sol_a)
    sol_f_p = np.zeros_like(rel_sol_a)
    sol_g_p = np.zeros_like(rel_sol_a)
    sol_h_p = np.zeros_like(rel_sol_a)
    sol_i_p = np.zeros_like(rel_sol_a)
    sol_j_p = np.zeros_like(rel_sol_a)

    #Project WF solutions onto theta surface
    for i in range(len(t)):
        sol_a_p[i] = rel_sol_a[i] - np.dot(rel_sol_a[i], normalhat)*normalhat
        sol_b_p[i] = rel_sol_b[i] - np.dot(rel_sol_b[i], normalhat)*normalhat
        sol_c_p[i] = rel_sol_c[i] - np.dot(rel_sol_c[i], normalhat)*normalhat
        sol_d_p[i] = rel_sol_d[i] - np.dot(rel_sol_d[i], normalhat)*normalhat
        sol_e_p[i] = rel_sol_e[i] - np.dot(rel_sol_e[i], normalhat)*normalhat
        sol_f_p[i] = rel_sol_f[i] - np.dot(rel_sol_f[i], normalhat)*normalhat
        sol_g_p[i] = rel_sol_g[i] - np.dot(rel_sol_g[i], normalhat)*normalhat
        sol_h_p[i] = rel_sol_h[i] - np.dot(rel_sol_h[i], normalhat)*normalhat
        sol_i_p[i] = rel_sol_i[i] - np.dot(rel_sol_i[i], normalhat)*normalhat
        sol_j_p[i] = rel_sol_j[i] - np.dot(rel_sol_j[i], normalhat)*normalhat

    projected_solutions = [sol_a_p, sol_b_p, sol_c_p, sol_d_p, sol_e_p, sol_f_p, sol_g_p, sol_h_p, sol_i_p, sol_j_p]

###############################################################################################

    #Plot the full 3d trajectories in the Earth frame
    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Earth frame')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.set_zlabel('z (H)')
    for i in EFsolutions:
        ax1.plot(i[:,0], i[:,1], i[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    for i in EFsolutions:
        ax2.plot(i[:,0], i[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    for i in EFsolutions:
        ax3.plot(i[:,1], i[:,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    for i in EFsolutions:
        ax4.plot(i[:,0], i[:,1])

###############################################################################################

    #Plot the full 3d trajectories in the Wave frame
    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('Wave frame')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('y (L)')
    ax5.set_zlabel('z (H)')
    for i in WFsolutions:
        ax5.plot(i[:,0], i[:,1], i[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    for i in WFsolutions:
        ax6.plot(i[:,0], i[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    for i in WFsolutions:
        ax7.plot(i[:,1], i[:,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    for i in WFsolutions:
        ax8.plot(i[:,0], i[:,1])

###############################################################################################
    '''
    #Plot parcels a and b separation with time

    fig3 = plt.figure()
    fig3.suptitle('Parcel separation')
    ax9 = fig3.add_subplot(111)
    ax9.plot(t,distances, color='black')
    ax9.set_xlabel('time (T)')
    ax9.set_ylabel('separation')
    ax9.axhline(y=d_i,color='black',ls='dotted')
    '''

###############################################################################################
    '''
    #Plot background potential temperature distribution

    fig4 =plt.figure()
    plt.set_cmap('inferno')
    fig4.suptitle('Basic state potential temperature')
    ax10 = fig4.add_subplot(111)
    ax10.set_xlabel('y (L)')
    ax10.set_ylabel('z (H)')
    thetacontour = ax10.contourf(theta_matrix, origin='lower', aspect='auto', extent=[ymin,ymax,zmin,zmax])
    plt.colorbar(thetacontour)

    ax10.plot(thetayvalues, gradient*thetayvalues)
    '''
###############################################################################################

    #Earth frame/Wave frame comparison figure

    fig5 = plt.figure()
    fig5.suptitle('Earth frame, left; and Wave frame, right')
    ax11 = fig5.add_subplot(221, projection = '3d')
    ax11.set_xlabel('x (L)')
    ax11.set_ylabel('y (L)')
    ax11.set_zlabel('z (H)')
    for i in EFsolutions:
        ax11.plot(i[:,0], i[:,1], i[:,2])

    ax12 = fig5.add_subplot(222, projection = '3d')
    ax12.set_xlabel('x (L)')
    ax12.set_ylabel('y (L)')
    ax12.set_zlabel('z (H)')
    for i in WFsolutions:
        ax12.plot(i[:,0], i[:,1], i[:,2])

    ax13 = fig5.add_subplot(223)
    ax13.set_xlabel('x (L)')
    ax13.set_ylabel('y (L)')
    for i in EFsolutions:
        ax13.plot(i[:,0], i[:,1])

    ax14 = fig5.add_subplot(224)
    ax14.set_xlabel('x (L)')
    ax14.set_ylabel('y (L)')
    for i in WFsolutions:
        ax14.plot(i[:,0], i[:,1])

###############################################################################################

    #WF isentropic surface projection figure, THIS IS THE GREEN PICTURE WE WANT

    fig6 = plt.figure()
    fig6.suptitle('Relative motion within isentropic surface')
    ax15 = fig6.add_subplot(111, projection = '3d')
    ax15.set_xlabel('x (L)')
    ax15.set_ylabel('y (L)')
    ax15.set_zlabel('z (H)')
    #for i in WFsolutions:
    #    ax15.plot(i[:,0], i[:,1], i[:,2])
    for i in projected_solutions:
        ax15.plot(i[:,0], i[:,1], i[:,2])

    times = np.arange(0,500,20)

    '''
    #Add a dot along trajectories every 20 timesteps to indicate time evolution
    for i in range(len(times)):
        for j in projected_solutions:
            ax15.scatter(j[times[i],0], j[times[i],1], j[times[i],2], marker='o', c='black', s=8)
    '''
###############################################################################################

    #Plot of meridional oscillation extent with time
    fig7 = plt.figure()
    fig7.suptitle('Meridional extent with time')
    ax16 = fig7.add_subplot(111)
    meanyseparation = np.zeros_like(t)

    for i in range(len(t)):
        meanyseparation[i]=(1./len(WFsolutions))*np.sum((j[i,1] - j[0,1])**2 for j in WFsolutions)

    for i in WFsolutions:
        ax16.plot(t, (i[:,1] - i[0,1])**2)
        ax16.plot(t, meanyseparation, linewidth=2.5, color='black')


    plt.show()

###############################################################################################

#Run the programme

main(start=0, stop=70) #60 time periods is ~28 days
