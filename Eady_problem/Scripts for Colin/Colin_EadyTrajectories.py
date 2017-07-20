from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

from Colin_Eadyinfo import *

###############################################################################################

#Programme to calculate and plot fluid parcel trajectories within an unstable Eady wave
#Trajectories are caluculated using both the implicit midpoint rule and the lsoda integrator

#Note that trajectories are confined to a meridional channel with -1 < y < + 1
#Confined by vertical boundaries with 0 < z < H
#Periodic boundary conditions in the x direction

def main(start, stop, zpos=0.5, ypos=0.5):

    def velocity(t,s):
        x,y,z = s
        dsdt = np.array((uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)))
        return dsdt

    def velocity2(s,t): #Need this for lsoda integrator; it takes s and t in opposite order
        x,y,z = s
        dsdt = [uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = start
    tmax = stop
    nt = 500

    t = np.linspace(tmin, tmax, nt)

    dt = (tmax - tmin)/nt

###############################################################################################

    #Define initial positions of parcels

    #Note: zpos=0.5 determines steering level, i.e. where c.real matches mean flow speed

    xpos = 0

    s0_a = np.array((xpos, -ypos, zpos))
    s0_b = np.array((xpos, 0, zpos))
    s0_c = np.array((xpos, ypos, zpos))


    #Define empty arrays for solutions and set their initial positions

    empty_a = np.empty((len(t), 3))
    empty_b = np.empty((len(t), 3))
    empty_c = np.empty((len(t), 3))

    empty_a[0] = s0_a
    empty_b[0] = s0_b
    empty_c[0] = s0_c

###############################################################################################

    #Solve for parcel trajectories

    #Integrate numerically using implicit midpoint rule

    def implicitmidpoint(solution):
        step = 1
        while t[step] < tmax:
            solution[step] = solution[step-1] + dt*velocity(t=(t[step]+t[step-1])/2., s=(solution[step]+solution[step-1])/2.) #integrate using rule
            step += 1 #increase the step value by 1
        return solution

    sol_a = implicitmidpoint(empty_a)
    sol_b = implicitmidpoint(empty_b)
    sol_c = implicitmidpoint(empty_c)

    #Create list of Earth frame trajectories
    EFsolutions = np.array((sol_a, sol_b, sol_c))

###############################################################################################

    #Transform to wave frame by shifting x-coordinates

    shift = np.zeros_like(sol_a)
    shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

    #Define new x-shifted trajectories
    rel_sol_a = sol_a + shift #rel for relative motion
    rel_sol_b = sol_b + shift
    rel_sol_c = sol_c + shift

    #Create list of Wave frame trajectories
    WFsolutions = np.array((rel_sol_a, rel_sol_b, rel_sol_c))

###############################################################################################

    #Solve for parcel trajectories AGAIN but using Scipy's odeint - 'lsoda' integrator

    sol_a2 = odeint(velocity2, s0_a, t)
    sol_b2 = odeint(velocity2, s0_b, t)
    sol_c2 = odeint(velocity2, s0_c, t)

    EFsolutions2 = [sol_a2, sol_b2, sol_c2]

    #Define new x-shifted trajectories
    rel_sol_a2 = sol_a2 + shift #rel for relative motion
    rel_sol_b2 = sol_b2 + shift
    rel_sol_c2 = sol_c2 + shift

    #Create list of Wave frame trajectories
    WFsolutions2 = [rel_sol_a2, rel_sol_b2, rel_sol_c2]

###############################################################################################

    #Plot the trajectories
    #Note, we plot the solution sol[:-1, . ] i.e. up to the last element in the array
    #This is because how the integration scheme is written, the last element will be zero

    #Plot in the Earth frame

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Earth frame')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.set_zlabel('z (H)')
    for i in EFsolutions:
    #for i in EFsolutions2:
        ax1.plot(i[:-1,0], i[:-1,1], i[:-1,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_title('x-z plane', fontsize=10)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    for i in EFsolutions:
    #for i in EFsolutions2:
        ax2.plot(i[:-1,0], i[:-1,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_title('y-z plane', fontsize=10)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    for i in EFsolutions:
    #for i in EFsolutions2:
        ax3.plot(i[:-1,1], i[:-1,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_title('x-y plane', fontsize=10)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    for i in EFsolutions:
    #for i in EFsolutions2:
        ax4.plot(i[:-1,0], i[:-1,1])

###############################################################################################

    #Plot in the Wave frame

    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('Wave frame')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('y (L)')
    ax5.set_zlabel('z (H)')
    for i in WFsolutions:
    #for i in WFsolutions2:
        ax5.plot(i[:-1,0], i[:-1,1], i[:-1,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_title('x-z plane', fontsize=10)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    for i in WFsolutions:
    #for i in WFsolutions2:
        ax6.plot(i[:-1,0], i[:-1,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_title('y-z plane', fontsize=10)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    for i in WFsolutions:
    #for i in WFsolutions2:
        ax7.plot(i[:-1,1], i[:-1,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_title('x-y plane', fontsize=10)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    for i in WFsolutions:
    #for i in WFsolutions2:
        ax8.plot(i[:-1,0], i[:-1,1])


    plt.show()

###############################################################################################

#Run the programme

main(start=0, stop=30, zpos=0.7, ypos=0.5)
