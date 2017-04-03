import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

from EadyUnstableinfo import *
from EadyUnstableSolution import *

#MAKE DOUBLY SURE THE TRANSFORMATION IS CORRECT.
#AT THE MOMENT WE SIMPLY SHIFT THE X-COORDS USING THE PHASE SPEED.

#READ EADY MODEL SOLUTION IN VALLIS. FIND AND READ EADY TRAJECTORY PAPERS.

#WHY ARE THE ISOTHERMALS SLANTING THE WRONG WAY??

###############################################################################################

#Print out physical quantities

print 'Horizontal lengthscale, L is (m):', L
print 'Height of vertical boundaries, H is (m):', H
print 'Rossby height, H_R is (m):', H_R
print 'H/H_R is:', Hratio
print 'Time period, T is (s):', T
print 'Horizontal wavenumber, k is (m^-1):', k
print 'Maximum growth rate, sigma_max is (s^-1):', sigma_max
print 'e-folding time (T) is (sigma_max*T)^(-1):', 1./(sigma_max*T)
print 'Buoyancy frequency, N is (s^-1):', N
print 'Maximum zonal velocity, U is (m*s^-1):', U
print 'Velocity shear is (s^-1):', shear
print 'Velocity scaling (to dimensionalise) is:', L*np.abs(f0)
print 'Phase speed c is (m*s^-1):', c
print 'Angular frequency, k*c is (s^-1):', k*c

###############################################################################################

#Use main function to avoid global variables
#Parameter p defines position of initial cluster, time (in T) defines period over which to
#integrate for parcel trajectories.

def main(p, time):

###############################################################################################

    #Define velocity function for 3d parcel trajectories

    def velocity3d(s,t):
        x,y,z = s
        dsdt = [umeanflow(z=z), vprime(x=x,z=z,t=t), wprime(x=x,z=z,t=t)]
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = 0
    tmax = time

    t = np.linspace(tmin, tmax, 200)

###############################################################################################

    #Define initial positions of parcels
    s0_a_3d = [p, p, p/5.]
    s0_b_3d = [p, p-0.1, (p-0.1)/5.]
    s0_c_3d = [p, p+0.1, (p+0.1)/5.]
    s0_d_3d = [p, p-0.2, (p-0.2)/5.]
    s0_e_3d = [p, p+0.2, (p+0.2)/5.]
    s0_f_3d = [p, p-0.3, (p-0.3)/5.]
    s0_g_3d = [p, p+0.3, (p+0.3)/5.]
    s0_h_3d = [p, p-0.4, (p-0.4)/5.]
    s0_i_3d = [p, p+0.4, (p+0.4)/5.]

    #Solve for parcel trajectories
    sol_a_3d = odeint(velocity3d, s0_a_3d, t)
    sol_b_3d = odeint(velocity3d, s0_b_3d, t)
    sol_c_3d = odeint(velocity3d, s0_c_3d, t)
    sol_d_3d = odeint(velocity3d, s0_d_3d, t)
    sol_e_3d = odeint(velocity3d, s0_e_3d, t)
    sol_f_3d = odeint(velocity3d, s0_f_3d, t)
    sol_g_3d = odeint(velocity3d, s0_g_3d, t)
    sol_h_3d = odeint(velocity3d, s0_h_3d, t)
    sol_i_3d = odeint(velocity3d, s0_i_3d, t)

###############################################################################################

    #Transform to wave frame by shifting x-coordinates
    shift = np.zeros_like(sol_a_3d)
    shift[:,0] = -c*T*t

    #Define new x-shifted trajectories
    rel_sol_a_3d = sol_a_3d + shift
    rel_sol_b_3d = sol_b_3d + shift
    rel_sol_c_3d = sol_c_3d + shift
    rel_sol_d_3d = sol_d_3d + shift
    rel_sol_e_3d = sol_e_3d + shift
    rel_sol_f_3d = sol_f_3d + shift
    rel_sol_g_3d = sol_g_3d + shift
    rel_sol_h_3d = sol_h_3d + shift
    rel_sol_i_3d = sol_i_3d + shift

###############################################################################################

    #Plot the full 3d trajectories in the Earth frame
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.set_zlabel('z (H)')
    ax1.plot(sol_a_3d[:,0], sol_a_3d[:,1], sol_a_3d[:,2])
    ax1.plot(sol_b_3d[:,0], sol_b_3d[:,1], sol_b_3d[:,2])
    ax1.plot(sol_c_3d[:,0], sol_c_3d[:,1], sol_c_3d[:,2])
    ax1.plot(sol_d_3d[:,0], sol_d_3d[:,1], sol_d_3d[:,2])
    ax1.plot(sol_e_3d[:,0], sol_e_3d[:,1], sol_e_3d[:,2])
    ax1.plot(sol_f_3d[:,0], sol_f_3d[:,1], sol_f_3d[:,2])
    ax1.plot(sol_g_3d[:,0], sol_g_3d[:,1], sol_g_3d[:,2])
    ax1.plot(sol_h_3d[:,0], sol_h_3d[:,1], sol_h_3d[:,2])
    ax1.plot(sol_i_3d[:,0], sol_i_3d[:,1], sol_i_3d[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    ax2.plot(sol_a_3d[:,0], sol_a_3d[:,2])
    ax2.plot(sol_b_3d[:,0], sol_b_3d[:,2])
    ax2.plot(sol_c_3d[:,0], sol_c_3d[:,2])
    ax2.plot(sol_d_3d[:,0], sol_d_3d[:,2])
    ax2.plot(sol_e_3d[:,0], sol_e_3d[:,2])
    ax2.plot(sol_f_3d[:,0], sol_f_3d[:,2])
    ax2.plot(sol_g_3d[:,0], sol_g_3d[:,2])
    ax2.plot(sol_h_3d[:,0], sol_h_3d[:,2])
    ax2.plot(sol_i_3d[:,0], sol_i_3d[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    ax3.plot(sol_a_3d[:,1], sol_a_3d[:,2])
    ax3.plot(sol_b_3d[:,1], sol_b_3d[:,2])
    ax3.plot(sol_c_3d[:,1], sol_c_3d[:,2])
    ax3.plot(sol_d_3d[:,1], sol_d_3d[:,2])
    ax3.plot(sol_e_3d[:,1], sol_e_3d[:,2])
    ax3.plot(sol_f_3d[:,1], sol_f_3d[:,2])
    ax3.plot(sol_g_3d[:,1], sol_g_3d[:,2])
    ax3.plot(sol_h_3d[:,1], sol_h_3d[:,2])
    ax3.plot(sol_i_3d[:,1], sol_i_3d[:,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    ax4.plot(sol_a_3d[:,0], sol_a_3d[:,1])
    ax4.plot(sol_b_3d[:,0], sol_b_3d[:,1])
    ax4.plot(sol_c_3d[:,0], sol_c_3d[:,1])
    ax4.plot(sol_d_3d[:,0], sol_d_3d[:,1])
    ax4.plot(sol_e_3d[:,0], sol_e_3d[:,1])
    ax4.plot(sol_f_3d[:,0], sol_f_3d[:,1])
    ax4.plot(sol_g_3d[:,0], sol_g_3d[:,1])
    ax4.plot(sol_h_3d[:,0], sol_h_3d[:,1])
    ax4.plot(sol_i_3d[:,0], sol_i_3d[:,1])

    plt.savefig('figures/EarthFrame/c_%s_p_%s_t_%s.pdf' % (str(c), str(p), str(tmax)))
###############################################################################################

    #Plot the full 3d trajectories in the Wave frame
    fig2 = plt.figure()
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('y (L)')
    ax5.set_zlabel('z (H)')
    ax5.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,1], rel_sol_a_3d[:,2])
    ax5.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,1], rel_sol_b_3d[:,2])
    ax5.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,1], rel_sol_c_3d[:,2])
    ax5.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,1], rel_sol_d_3d[:,2])
    ax5.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,1], rel_sol_e_3d[:,2])
    ax5.plot(rel_sol_f_3d[:,0], rel_sol_f_3d[:,1], rel_sol_f_3d[:,2])
    ax5.plot(rel_sol_g_3d[:,0], rel_sol_g_3d[:,1], rel_sol_g_3d[:,2])
    ax5.plot(rel_sol_h_3d[:,0], rel_sol_h_3d[:,1], rel_sol_h_3d[:,2])
    ax5.plot(rel_sol_i_3d[:,0], rel_sol_i_3d[:,1], rel_sol_i_3d[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    ax6.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,2])
    ax6.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,2])
    ax6.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,2])
    ax6.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,2])
    ax6.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,2])
    ax6.plot(rel_sol_f_3d[:,0], rel_sol_f_3d[:,2])
    ax6.plot(rel_sol_g_3d[:,0], rel_sol_g_3d[:,2])
    ax6.plot(rel_sol_h_3d[:,0], rel_sol_h_3d[:,2])
    ax6.plot(rel_sol_i_3d[:,0], rel_sol_i_3d[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    ax7.plot(rel_sol_a_3d[:,1], rel_sol_a_3d[:,2])
    ax7.plot(rel_sol_b_3d[:,1], rel_sol_b_3d[:,2])
    ax7.plot(rel_sol_c_3d[:,1], rel_sol_c_3d[:,2])
    ax7.plot(rel_sol_d_3d[:,1], rel_sol_d_3d[:,2])
    ax7.plot(rel_sol_e_3d[:,1], rel_sol_e_3d[:,2])
    ax7.plot(rel_sol_f_3d[:,1], rel_sol_f_3d[:,2])
    ax7.plot(rel_sol_g_3d[:,1], rel_sol_g_3d[:,2])
    ax7.plot(rel_sol_h_3d[:,1], rel_sol_h_3d[:,2])
    ax7.plot(rel_sol_i_3d[:,1], rel_sol_i_3d[:,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    ax8.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,1])
    ax8.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,1])
    ax8.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,1])
    ax8.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,1])
    ax8.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,1])
    ax8.plot(rel_sol_f_3d[:,0], rel_sol_f_3d[:,1])
    ax8.plot(rel_sol_g_3d[:,0], rel_sol_g_3d[:,1])
    ax8.plot(rel_sol_h_3d[:,0], rel_sol_h_3d[:,1])
    ax8.plot(rel_sol_i_3d[:,0], rel_sol_i_3d[:,1])

    #plt.savefig('figures/WaveFrame/c_%s_p_%s_t_%s.pdf' % (str(c), str(p), str(tmax)))
    #plt.close()

    plt.show()

###############################################################################################

#Run the programme

#pvalues = [i for i in np.arange(-4,4,0.5)]

#for i in pvalues:
    #main(p=i, time=5)

main(p=0, time=5)
