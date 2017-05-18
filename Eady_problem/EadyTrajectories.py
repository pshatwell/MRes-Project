
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

    top = 0.53
    bottom = 0.47

    s0_a = np.array((0, -0.5, top))
    s0_b = np.array((0.5, -0.5, top))
    s0_c = np.array((1, -0.5, top))
    s0_d = np.array((1.5, -0.5, top))
    s0_e = np.array((2, -0.5, top))

    s0_f = np.array((3, 0.5, bottom))
    s0_g = np.array((3.5, 0.5, bottom))
    s0_h = np.array((4, 0.5, bottom))
    s0_i = np.array((4.5, 0.5, bottom))
    s0_j = np.array((5, 0.5, bottom))


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

###############################################################################################

    #Distance calculations

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


###############################################################################################

    #Create matrix for background potential temperature distribution
    theta_matrix = np.zeros((50,50))

    ymin = sol_a[0,1]
    ymax = sol_a[-1,1]

    zmin = sol_a[0,2]
    zmax = sol_a[-1,2]

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
    ax1.plot(sol_a[:,0], sol_a[:,1], sol_a[:,2])
    ax1.plot(sol_b[:,0], sol_b[:,1], sol_b[:,2])
    ax1.plot(sol_c[:,0], sol_c[:,1], sol_c[:,2])
    ax1.plot(sol_d[:,0], sol_d[:,1], sol_d[:,2])
    ax1.plot(sol_e[:,0], sol_e[:,1], sol_e[:,2])

    ax1.plot(sol_f[:,0], sol_f[:,1], sol_f[:,2])
    ax1.plot(sol_g[:,0], sol_g[:,1], sol_g[:,2])
    ax1.plot(sol_h[:,0], sol_h[:,1], sol_h[:,2])
    ax1.plot(sol_i[:,0], sol_i[:,1], sol_i[:,2])
    ax1.plot(sol_j[:,0], sol_j[:,1], sol_j[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    ax2.plot(sol_a[:,0], sol_a[:,2])
    ax2.plot(sol_b[:,0], sol_b[:,2])
    ax2.plot(sol_c[:,0], sol_c[:,2])
    ax2.plot(sol_d[:,0], sol_d[:,2])
    ax2.plot(sol_e[:,0], sol_e[:,2])

    ax2.plot(sol_f[:,0], sol_f[:,2])
    ax2.plot(sol_g[:,0], sol_g[:,2])
    ax2.plot(sol_h[:,0], sol_h[:,2])
    ax2.plot(sol_i[:,0], sol_i[:,2])
    ax2.plot(sol_j[:,0], sol_j[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    ax3.plot(sol_a[:,1], sol_a[:,2])
    ax3.plot(sol_b[:,1], sol_b[:,2])
    ax3.plot(sol_c[:,1], sol_c[:,2])
    ax3.plot(sol_d[:,1], sol_d[:,2])
    ax3.plot(sol_e[:,1], sol_e[:,2])

    ax3.plot(sol_f[:,1], sol_f[:,2])
    ax3.plot(sol_g[:,1], sol_g[:,2])
    ax3.plot(sol_h[:,1], sol_h[:,2])
    ax3.plot(sol_i[:,1], sol_i[:,2])
    ax3.plot(sol_j[:,1], sol_j[:,2])


    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    ax4.plot(sol_a[:,0], sol_a[:,1])
    ax4.plot(sol_b[:,0], sol_b[:,1])
    ax4.plot(sol_c[:,0], sol_c[:,1])
    ax4.plot(sol_d[:,0], sol_d[:,1])
    ax4.plot(sol_e[:,0], sol_e[:,1])

    ax4.plot(sol_f[:,0], sol_f[:,1])
    ax4.plot(sol_g[:,0], sol_g[:,1])
    ax4.plot(sol_h[:,0], sol_h[:,1])
    ax4.plot(sol_i[:,0], sol_i[:,1])
    ax4.plot(sol_j[:,0], sol_j[:,1])

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
    ax5.plot(rel_sol_c[:,0], rel_sol_c[:,1], rel_sol_c[:,2])
    ax5.plot(rel_sol_d[:,0], rel_sol_d[:,1], rel_sol_d[:,2])
    ax5.plot(rel_sol_e[:,0], rel_sol_e[:,1], rel_sol_e[:,2])

    ax5.plot(rel_sol_f[:,0], rel_sol_f[:,1], rel_sol_f[:,2])
    ax5.plot(rel_sol_g[:,0], rel_sol_g[:,1], rel_sol_g[:,2])
    ax5.plot(rel_sol_h[:,0], rel_sol_h[:,1], rel_sol_h[:,2])
    ax5.plot(rel_sol_i[:,0], rel_sol_i[:,1], rel_sol_i[:,2])
    ax5.plot(rel_sol_j[:,0], rel_sol_j[:,1], rel_sol_j[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    ax6.plot(rel_sol_a[:,0], rel_sol_a[:,2])
    ax6.plot(rel_sol_b[:,0], rel_sol_b[:,2])
    ax6.plot(rel_sol_c[:,0], rel_sol_c[:,2])
    ax6.plot(rel_sol_d[:,0], rel_sol_d[:,2])
    ax6.plot(rel_sol_e[:,0], rel_sol_e[:,2])

    ax6.plot(rel_sol_f[:,0], rel_sol_f[:,2])
    ax6.plot(rel_sol_g[:,0], rel_sol_g[:,2])
    ax6.plot(rel_sol_h[:,0], rel_sol_h[:,2])
    ax6.plot(rel_sol_i[:,0], rel_sol_i[:,2])
    ax6.plot(rel_sol_j[:,0], rel_sol_j[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    ax7.plot(rel_sol_a[:,1], rel_sol_a[:,2])
    ax7.plot(rel_sol_b[:,1], rel_sol_b[:,2])
    ax7.plot(rel_sol_c[:,1], rel_sol_c[:,2])
    ax7.plot(rel_sol_d[:,1], rel_sol_d[:,2])
    ax7.plot(rel_sol_e[:,1], rel_sol_e[:,2])

    ax7.plot(rel_sol_f[:,1], rel_sol_f[:,2])
    ax7.plot(rel_sol_g[:,1], rel_sol_g[:,2])
    ax7.plot(rel_sol_h[:,1], rel_sol_h[:,2])
    ax7.plot(rel_sol_i[:,1], rel_sol_i[:,2])
    ax7.plot(rel_sol_j[:,1], rel_sol_j[:,2])


    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    ax8.plot(rel_sol_a[:,0], rel_sol_a[:,1])
    ax8.plot(rel_sol_b[:,0], rel_sol_b[:,1])
    ax8.plot(rel_sol_c[:,0], rel_sol_c[:,1])
    ax8.plot(rel_sol_d[:,0], rel_sol_d[:,1])
    ax8.plot(rel_sol_e[:,0], rel_sol_e[:,1])

    ax8.plot(rel_sol_f[:,0], rel_sol_f[:,1])
    ax8.plot(rel_sol_g[:,0], rel_sol_g[:,1])
    ax8.plot(rel_sol_h[:,0], rel_sol_h[:,1])
    ax8.plot(rel_sol_i[:,0], rel_sol_i[:,1])
    ax8.plot(rel_sol_j[:,0], rel_sol_j[:,1])

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

main(start=0, stop=50)
