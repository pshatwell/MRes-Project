
'''Script to plot fluid parcel trajectories for Eady model.'''


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

    tmin = start #WHAT SHOULD THIS BE? DO WE WANT THERE TO ALREADY BE A LARGE INSTABILITY?
    tmax = stop

    t = np.linspace(tmin, tmax, 500)

###############################################################################################

    #Define initial positions of parcels

    s0_a_3d = np.array((0.1, 0.5, 0.7))
    s0_b_3d = np.array((0.2, 0.5, 0.7))
    s0_c_3d = np.array((0.3, 0.5, 0.7))
    s0_d_3d = np.array((3.1, 0.5, 0.8))
    s0_e_3d = np.array((3.2, 0.5, 0.8))

    #Solve for parcel trajectories
    sol_a_3d = odeint(velocity3d, s0_a_3d, t)
    sol_b_3d = odeint(velocity3d, s0_b_3d, t)
    sol_c_3d = odeint(velocity3d, s0_c_3d, t)
    sol_d_3d = odeint(velocity3d, s0_d_3d, t)
    sol_e_3d = odeint(velocity3d, s0_e_3d, t)

###############################################################################################

    #Distance calculations

    #Calculate initial parcel separation
    d_i = la.norm(s0_b_3d - s0_a_3d)
    print 'initial separation is:', d_i

    #Calculate final parcel separation
    d_f = la.norm(sol_b_3d[-1] - sol_a_3d[-1])
    print 'final separation is:', d_f

    distances = []
    for i in range(len(t)):
        distances.append(la.norm(sol_b_3d[i] - sol_a_3d[i]))

###############################################################################################


    #Transform to wave frame by shifting x-coordinates
    shift = np.zeros_like(sol_a_3d)
    shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

    #Define new x-shifted trajectories
    rel_sol_a_3d = sol_a_3d + shift
    rel_sol_b_3d = sol_b_3d + shift
    rel_sol_c_3d = sol_c_3d + shift
    rel_sol_d_3d = sol_d_3d + shift
    rel_sol_e_3d = sol_e_3d + shift


###############################################################################################

    #Create matrix for background potential temperature distribution
    theta_matrix = np.zeros((50,50))

    ymin = sol_a_3d[0,1]
    ymax = sol_a_3d[-1,1]

    zmin = sol_a_3d[0,2]
    zmax = sol_a_3d[-1,2]

    thetayvalues = np.linspace(ymin,ymax,50)
    thetazvalues = np.linspace(zmin,zmax,50)

    for i in range(0, 50, 1):
        for j in range(0, 50, 1):
            theta_matrix[i,j] = theta(y=thetayvalues[j], z=thetazvalues[i])

###############################################################################################

    #Plot the full 3d trajectories in the Earth frame
    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Earth frame')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.set_zlabel('z (H)')
    ax1.plot(sol_a_3d[:,0], sol_a_3d[:,1], sol_a_3d[:,2])
    ax1.plot(sol_b_3d[:,0], sol_b_3d[:,1], sol_b_3d[:,2])
    ax1.plot(sol_c_3d[:,0], sol_c_3d[:,1], sol_c_3d[:,2])
    ax1.plot(sol_d_3d[:,0], sol_d_3d[:,1], sol_d_3d[:,2])
    ax1.plot(sol_e_3d[:,0], sol_e_3d[:,1], sol_e_3d[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    ax2.plot(sol_a_3d[:,0], sol_a_3d[:,2])
    ax2.plot(sol_b_3d[:,0], sol_b_3d[:,2])
    ax2.plot(sol_c_3d[:,0], sol_c_3d[:,2])
    ax2.plot(sol_d_3d[:,0], sol_d_3d[:,2])
    ax2.plot(sol_e_3d[:,0], sol_e_3d[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    ax3.plot(sol_a_3d[:,1], sol_a_3d[:,2])
    ax3.plot(sol_b_3d[:,1], sol_b_3d[:,2])
    ax3.plot(sol_c_3d[:,1], sol_c_3d[:,2])
    ax3.plot(sol_d_3d[:,1], sol_d_3d[:,2])
    ax3.plot(sol_e_3d[:,1], sol_e_3d[:,2])
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    ax4.plot(sol_a_3d[:,0], sol_a_3d[:,1])
    ax4.plot(sol_b_3d[:,0], sol_b_3d[:,1])
    ax4.plot(sol_c_3d[:,0], sol_c_3d[:,1])
    ax4.plot(sol_d_3d[:,0], sol_d_3d[:,1])
    ax4.plot(sol_e_3d[:,0], sol_e_3d[:,1])

###############################################################################################

    #Plot the full 3d trajectories in the Wave frame
    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('Wave frame')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('y (L)')
    ax5.set_zlabel('z (H)')
    ax5.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,1], rel_sol_a_3d[:,2])
    ax5.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,1], rel_sol_b_3d[:,2])
    ax5.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,1], rel_sol_c_3d[:,2])
    ax5.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,1], rel_sol_d_3d[:,2])
    ax5.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,1], rel_sol_e_3d[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    ax6.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,2])
    ax6.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,2])
    ax6.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,2])
    ax6.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,2])
    ax6.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    ax7.plot(rel_sol_a_3d[:,1], rel_sol_a_3d[:,2])
    ax7.plot(rel_sol_b_3d[:,1], rel_sol_b_3d[:,2])
    ax7.plot(rel_sol_c_3d[:,1], rel_sol_c_3d[:,2])
    ax7.plot(rel_sol_d_3d[:,1], rel_sol_d_3d[:,2])
    ax7.plot(rel_sol_e_3d[:,1], rel_sol_e_3d[:,2])
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    ax8.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,1])
    ax8.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,1])
    ax8.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,1])
    ax8.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,1])
    ax8.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,1])

###############################################################################################

    #Plot parcel separation with time

    fig3 = plt.figure()
    fig3.suptitle('Parcel separation')
    ax9 = fig3.add_subplot(111)
    ax9.plot(t,distances, color='black')
    ax9.set_xlabel('time (T)')
    ax9.set_ylabel('separation')
    ax9.axhline(y=d_i,color='black',ls='dotted')


###############################################################################################

    #Plot background potential temperature distribution
    fig4 =plt.figure()
    plt.set_cmap('inferno')
    fig4.suptitle('Basic state potential temperature')
    ax10 = fig4.add_subplot(111)
    ax10.set_xlabel('y (L)')
    ax10.set_ylabel('z (H)')
    thetacontour = ax10.contourf(theta_matrix, origin='lower', aspect='auto', extent=[ymin,ymax,zmin,zmax])
    plt.colorbar(thetacontour)

    plt.show()

###############################################################################################

#Run the programme

main(start=15, stop=35)
