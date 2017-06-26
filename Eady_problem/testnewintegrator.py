
'''Script to demonstrate Green's idea of isentropic relative flow
using unstable solutions of the Eady model.'''

'''THETA IS NOT MATERIALLY CONSERVED?'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
from scipy.integrate import ode
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

def main(start, stop, zpos):

    #Define velocity function for 3d parcel trajectories

    def velocity3d(t,s):
        x,y,z = s
        dsdt = [uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = start #Should this necessarily be zero?
    tmax = stop

    t = np.linspace(tmin, tmax, 500)

###############################################################################################

    #Define initial positions of parcels

    #Setting both 'cold' and 'warm' parcels off at the same height but different y
    #seems to give closest reproduction of Green picture

    #zpos = 0.01 #zpos=0.5 determines steering level, i.e. where c.real matches mean flow speed
    xpos = 0
    xshift = np.pi/(k*L) #nondimensional shift of half a wavelength
    ypos = 0.8

    #Set of 5 'warm' parcels (positive y)
    s0_a = np.array((xpos, ypos, zpos))

###############################################################################################

    #Solve for parcel trajectories
    solver1 = ode(velocity3d)
    solver1.set_integrator('dop853')
    solver1.set_initial_value(s0_a, tmin)

    solver2 = ode(velocity3d)
    solver2.set_integrator('lsoda')
    solver2.set_initial_value(s0_a, tmin)

    sol_a = np.empty((len(t), 3))
    sol_a[0] = s0_a

    sol_b = np.empty((len(t), 3))
    sol_b[0] = s0_a

    step1 = 1
    while solver1.successful() and solver1.t < tmax:
        solver1.integrate(t[step1])
        sol_a[step1] = solver1.y
        step1 += 1

    step2 = 1
    while solver2.successful() and solver2.t < tmax:
        solver2.integrate(t[step2])
        sol_b[step2] = solver2.y
        step2 += 1

###############################################################################################

    #Transform to wave frame by shifting x-coordinates
    shift = np.zeros_like(sol_a)
    shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

    #Define new x-shifted trajectories
    rel_sol_a = sol_a + shift #rel for relative motion
    rel_sol_b = sol_b + shift

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

    #Plot the full 3d trajectories in the Earth frame
    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Earth frame')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.set_zlabel('z (H)')
    ax1.plot(sol_a[:,0], sol_a[:,1], sol_a[:,2])
    ax1.plot(sol_b[:,0], sol_b[:,1], sol_b[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_title('x-z plane', fontsize=10)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    ax2.plot(sol_a[:,0], sol_a[:,2])
    ax2.plot(sol_b[:,0], sol_b[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_title('y-z plane', fontsize=10)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    ax3.plot(sol_a[:,1], sol_a[:,2])
    ax3.plot(sol_b[:,1], sol_b[:,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_title('x-y plane', fontsize=10)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    ax4.plot(sol_a[:,0], sol_a[:,1])
    ax4.plot(sol_b[:,0], sol_b[:,1])

###############################################################################################

    #Plot the full 3d trajectories in the Wave frame
    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('Wave frame')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('y (L)')
    ax5.set_zlabel('z (H)')
    ax5.plot(rel_sol_a[:,0], rel_sol_a[:,1], rel_sol_a[:,2])
    ax5.plot(rel_sol_b[:,0], rel_sol_b[:,1], rel_sol_b[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_title('x-z plane', fontsize=10)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    ax6.plot(rel_sol_a[:,0], rel_sol_a[:,2])
    ax6.plot(rel_sol_b[:,0], rel_sol_b[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_title('y-z plane', fontsize=10)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    ax7.plot(rel_sol_a[:,1], rel_sol_a[:,2])
    ax7.plot(rel_sol_b[:,1], rel_sol_b[:,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_title('x-y plane', fontsize=10)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    ax8.plot(rel_sol_a[:,0], rel_sol_a[:,1])
    ax8.plot(rel_sol_b[:,0], rel_sol_b[:,1])

    plt.show()

###############################################################################################

#Run the programme

main(start=0, stop=40, zpos=0.5) #60 time periods is ~28 days
