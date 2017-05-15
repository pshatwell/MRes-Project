import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import numpy.linalg as la

from EadyUnstableinfoTEST import *
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
print 'Zonal wavenumber, k is (m^-1):', k
print 'Meridional wavenumber, l is (m^-1)', l
#print 'Maximum growth rate, sigma_max is (s^-1):', sigma_max
print 'Growth rate, sigma is (s^-1):', sigma
print 'e-folding time (T) is (sigma*T)^(-1):', 1./(sigma*T)
print 'Buoyancy frequency, N is (s^-1):', N
print 'Maximum zonal velocity, U0 is (m*s^-1):', U0
print 'Velocity shear is (s^-1):', shear
print 'Velocity scaling (to dimensionalise) is:', L*np.abs(f0)
print 'Phase speed c is (m*s^-1):', c*L
print 'Angular frequency, k*c is (s^-1):', k*c.real*L

###############################################################################################

#Use main function to avoid global variables
#Parameter p defines position of initial cluster, time (in T) defines period over which to
#integrate for parcel trajectories.

def main(p, time):

###############################################################################################

    #Define velocity function for 3d parcel trajectories

    def velocity3d(s,t):
        x,y,z = s
        dsdt = [uprime(x=x,y=y,z=z,t=t)+umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = 0
    tmax = time

    t = np.linspace(tmin, tmax, 300)

###############################################################################################

    #Define initial positions of parcels
    '''
    s0_a_3d = [p, p, p/5.]
    s0_b_3d = [p, p-0.1, (p-0.1)/5.]
    s0_c_3d = [p, p+0.1, (p+0.1)/5.]
    s0_d_3d = [p, p-0.2, (p-0.2)/5.]
    s0_e_3d = [p, p+0.2, (p+0.2)/5.]
    s0_f_3d = [p, p-0.3, (p-0.3)/5.]
    s0_g_3d = [p, p+0.3, (p+0.3)/5.]
    s0_h_3d = [p, p-0.4, (p-0.4)/5.]
    s0_i_3d = [p, p+0.4, (p+0.4)/5.]
    '''
    s0_a_3d = np.array((3, 0.5, 0.1))
    s0_b_3d = np.array((4, 0.5, -0.1))

    #Solve for parcel trajectories
    sol_a_3d = odeint(velocity3d, s0_a_3d, t)
    sol_b_3d = odeint(velocity3d, s0_b_3d, t)

    '''
    sol_c_3d = odeint(velocity3d, s0_c_3d, t)
    sol_d_3d = odeint(velocity3d, s0_d_3d, t)
    sol_e_3d = odeint(velocity3d, s0_e_3d, t)
    sol_f_3d = odeint(velocity3d, s0_f_3d, t)
    sol_g_3d = odeint(velocity3d, s0_g_3d, t)
    sol_h_3d = odeint(velocity3d, s0_h_3d, t)
    sol_i_3d = odeint(velocity3d, s0_i_3d, t)
    '''

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
    shift[:,0] = -c.real*T*t

    #Define new x-shifted trajectories
    rel_sol_a_3d = sol_a_3d + shift
    rel_sol_b_3d = sol_b_3d + shift

    '''
    rel_sol_c_3d = sol_c_3d + shift
    rel_sol_d_3d = sol_d_3d + shift
    rel_sol_e_3d = sol_e_3d + shift
    rel_sol_f_3d = sol_f_3d + shift
    rel_sol_g_3d = sol_g_3d + shift
    rel_sol_h_3d = sol_h_3d + shift
    rel_sol_i_3d = sol_i_3d + shift
    '''
###############################################################################################

    #Plot the full 3d trajectories in the Earth frame
    fig1 = plt.figure()
    fig1.suptitle('Earth frame')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.set_zlabel('z (H)')
    ax1.plot(sol_a_3d[:,0], sol_a_3d[:,1], sol_a_3d[:,2])
    ax1.plot(sol_b_3d[:,0], sol_b_3d[:,1], sol_b_3d[:,2])
    '''
    ax1.plot(sol_c_3d[:,0], sol_c_3d[:,1], sol_c_3d[:,2])
    ax1.plot(sol_d_3d[:,0], sol_d_3d[:,1], sol_d_3d[:,2])
    ax1.plot(sol_e_3d[:,0], sol_e_3d[:,1], sol_e_3d[:,2])
    ax1.plot(sol_f_3d[:,0], sol_f_3d[:,1], sol_f_3d[:,2])
    ax1.plot(sol_g_3d[:,0], sol_g_3d[:,1], sol_g_3d[:,2])
    ax1.plot(sol_h_3d[:,0], sol_h_3d[:,1], sol_h_3d[:,2])
    ax1.plot(sol_i_3d[:,0], sol_i_3d[:,1], sol_i_3d[:,2])
    '''
    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    ax2.plot(sol_a_3d[:,0], sol_a_3d[:,2])
    ax2.plot(sol_b_3d[:,0], sol_b_3d[:,2])
    '''
    ax2.plot(sol_c_3d[:,0], sol_c_3d[:,2])
    ax2.plot(sol_d_3d[:,0], sol_d_3d[:,2])
    ax2.plot(sol_e_3d[:,0], sol_e_3d[:,2])
    ax2.plot(sol_f_3d[:,0], sol_f_3d[:,2])
    ax2.plot(sol_g_3d[:,0], sol_g_3d[:,2])
    ax2.plot(sol_h_3d[:,0], sol_h_3d[:,2])
    ax2.plot(sol_i_3d[:,0], sol_i_3d[:,2])
    '''


    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    ax3.plot(sol_a_3d[:,1], sol_a_3d[:,2])
    ax3.plot(sol_b_3d[:,1], sol_b_3d[:,2])
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    '''
    ax3.plot(sol_c_3d[:,1], sol_c_3d[:,2])
    ax3.plot(sol_d_3d[:,1], sol_d_3d[:,2])
    ax3.plot(sol_e_3d[:,1], sol_e_3d[:,2])
    ax3.plot(sol_f_3d[:,1], sol_f_3d[:,2])
    ax3.plot(sol_g_3d[:,1], sol_g_3d[:,2])
    ax3.plot(sol_h_3d[:,1], sol_h_3d[:,2])
    ax3.plot(sol_i_3d[:,1], sol_i_3d[:,2])
    '''
    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    ax4.plot(sol_a_3d[:,0], sol_a_3d[:,1])
    ax4.plot(sol_b_3d[:,0], sol_b_3d[:,1])
    '''
    ax4.plot(sol_c_3d[:,0], sol_c_3d[:,1])
    ax4.plot(sol_d_3d[:,0], sol_d_3d[:,1])
    ax4.plot(sol_e_3d[:,0], sol_e_3d[:,1])
    ax4.plot(sol_f_3d[:,0], sol_f_3d[:,1])
    ax4.plot(sol_g_3d[:,0], sol_g_3d[:,1])
    ax4.plot(sol_h_3d[:,0], sol_h_3d[:,1])
    ax4.plot(sol_i_3d[:,0], sol_i_3d[:,1])
    '''

    #plt.savefig('figures/EarthFrame/c_%s_p_%s_t_%s.pdf' % (str(c), str(p), str(tmax)))
###############################################################################################

    #Plot the full 3d trajectories in the Wave frame
    fig2 = plt.figure()
    fig2.suptitle('Wave frame')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)')
    ax5.set_ylabel('y (L)')
    ax5.set_zlabel('z (H)')
    ax5.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,1], rel_sol_a_3d[:,2])
    ax5.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,1], rel_sol_b_3d[:,2])
    '''
    ax5.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,1], rel_sol_c_3d[:,2])
    ax5.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,1], rel_sol_d_3d[:,2])
    ax5.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,1], rel_sol_e_3d[:,2])
    ax5.plot(rel_sol_f_3d[:,0], rel_sol_f_3d[:,1], rel_sol_f_3d[:,2])
    ax5.plot(rel_sol_g_3d[:,0], rel_sol_g_3d[:,1], rel_sol_g_3d[:,2])
    ax5.plot(rel_sol_h_3d[:,0], rel_sol_h_3d[:,1], rel_sol_h_3d[:,2])
    ax5.plot(rel_sol_i_3d[:,0], rel_sol_i_3d[:,1], rel_sol_i_3d[:,2])
    '''
    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    ax6.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,2])
    ax6.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,2])
    '''
    ax6.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,2])
    ax6.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,2])
    ax6.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,2])
    ax6.plot(rel_sol_f_3d[:,0], rel_sol_f_3d[:,2])
    ax6.plot(rel_sol_g_3d[:,0], rel_sol_g_3d[:,2])
    ax6.plot(rel_sol_h_3d[:,0], rel_sol_h_3d[:,2])
    ax6.plot(rel_sol_i_3d[:,0], rel_sol_i_3d[:,2])
    '''
    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    ax7.plot(rel_sol_a_3d[:,1], rel_sol_a_3d[:,2])
    ax7.plot(rel_sol_b_3d[:,1], rel_sol_b_3d[:,2])
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    '''
    ax7.plot(rel_sol_c_3d[:,1], rel_sol_c_3d[:,2])
    ax7.plot(rel_sol_d_3d[:,1], rel_sol_d_3d[:,2])
    ax7.plot(rel_sol_e_3d[:,1], rel_sol_e_3d[:,2])
    ax7.plot(rel_sol_f_3d[:,1], rel_sol_f_3d[:,2])
    ax7.plot(rel_sol_g_3d[:,1], rel_sol_g_3d[:,2])
    ax7.plot(rel_sol_h_3d[:,1], rel_sol_h_3d[:,2])
    ax7.plot(rel_sol_i_3d[:,1], rel_sol_i_3d[:,2])
    '''
    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    ax8.plot(rel_sol_a_3d[:,0], rel_sol_a_3d[:,1])
    ax8.plot(rel_sol_b_3d[:,0], rel_sol_b_3d[:,1])
    '''
    ax8.plot(rel_sol_c_3d[:,0], rel_sol_c_3d[:,1])
    ax8.plot(rel_sol_d_3d[:,0], rel_sol_d_3d[:,1])
    ax8.plot(rel_sol_e_3d[:,0], rel_sol_e_3d[:,1])
    ax8.plot(rel_sol_f_3d[:,0], rel_sol_f_3d[:,1])
    ax8.plot(rel_sol_g_3d[:,0], rel_sol_g_3d[:,1])
    ax8.plot(rel_sol_h_3d[:,0], rel_sol_h_3d[:,1])
    ax8.plot(rel_sol_i_3d[:,0], rel_sol_i_3d[:,1])
    '''

    #plt.savefig('figures/WaveFrame/c_%s_p_%s_t_%s.pdf' % (str(c), str(p), str(tmax)))
    #plt.close()

    fig3 = plt.figure()
    fig3.suptitle('Parcel separation')
    ax9 = fig3.add_subplot(111)
    ax9.plot(t,distances, color='black')
    ax9.set_xlabel('time (T)')
    ax9.set_ylabel('separation')
    ax9.axhline(y=d_i,color='black',ls='dotted')


    plt.show()

###############################################################################################

#Run the programme

#pvalues = [i for i in np.arange(-4,4,0.5)]

#for i in pvalues:
    #main(p=i, time=5)

main(p=0, time=3)
