from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import numpy.linalg as la

from Eadyinfo import *

'''WHY DOES MERIDIONAL EXTENT OF CHANNEL APPEAR TO BE -2 TO +2??'''

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

    def velocity(t,s):
        x,y,z = s
        dsdt = np.array((uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)))
        return dsdt

    def velocity2(s,t):
        x,y,z = s
        dsdt = [uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = start
    tmax = stop
    nt = 1000

    t = np.linspace(tmin, tmax, nt)

    dt = (tmax - tmin)/nt

    print 'length of t is', len(t)
    print 'dt is', dt

###############################################################################################

    #Define initial positions of parcels

    #Note: zpos=0.5 determines steering level, i.e. where c.real matches mean flow speed

    xpos = 0
    xshift = 2*np.pi/(k*L) #nondimensional shift of a wavelength
    ypos = 1.5

    #Set of 5 'warm' parcels (positive y, near top)
    s0_a = np.array((xpos, ypos, zpos))
    s0_b = np.array((xpos+0.5, ypos, zpos))
    s0_c = np.array((xpos+1, ypos, zpos))
    s0_d = np.array((xpos+1.5, ypos, zpos))
    s0_e = np.array((xpos+2, ypos, zpos))

    #Set of 5 'cold' parcels (negative y, near bottom)
    s0_f = np.array((xpos+xshift, -ypos, zpos))
    s0_g = np.array((xpos+xshift+0.5, -ypos, zpos))
    s0_h = np.array((xpos+xshift+1, -ypos, zpos))
    s0_i = np.array((xpos+xshift+1.5, -ypos, zpos))
    s0_j = np.array((xpos+xshift+2, -ypos, zpos))

    '''
    #Set of 5 'warm' parcels (positive y, near top)
    s0_a = np.array((xpos, -1.7, 0.9))
    s0_b = np.array((xpos+0.5, -1.7, 0.9))
    s0_c = np.array((xpos+1, -1.7, 0.9))
    s0_d = np.array((xpos+1.5, -1.7, 0.9))
    s0_e = np.array((xpos+2, -1.7, 0.9))

    #Set of 5 'cold' parcels (negative y, near bottom)
    s0_f = np.array((xpos+xshift, 0.5, 0.5))
    s0_g = np.array((xpos+xshift+0.5, 0.5, 0.5))
    s0_h = np.array((xpos+xshift+1, 0.5, 0.5))
    s0_i = np.array((xpos+xshift+1.5, 0.5, 0.5))
    s0_j = np.array((xpos+xshift+2, 0.5, 0.5))
    '''

    #Define empty arrays for solutions and set their initial positions

    empty_a = np.empty((len(t), 3))
    empty_b = np.empty((len(t), 3))
    empty_c = np.empty((len(t), 3))
    empty_d = np.empty((len(t), 3))
    empty_e = np.empty((len(t), 3))

    empty_f = np.empty((len(t), 3))
    empty_g = np.empty((len(t), 3))
    empty_h = np.empty((len(t), 3))
    empty_i = np.empty((len(t), 3))
    empty_j = np.empty((len(t), 3))


    empty_a[0] = s0_a
    empty_b[0] = s0_b
    empty_c[0] = s0_c
    empty_d[0] = s0_d
    empty_e[0] = s0_e

    empty_f[0] = s0_f
    empty_g[0] = s0_g
    empty_h[0] = s0_h
    empty_i[0] = s0_i
    empty_j[0] = s0_j

    print 'solution shape is', empty_a.shape
    print 'velocity type is', type(velocity(t=0, s=(0,0,0)))
    print 'velocity shape is', velocity(t=t,s=(0,0,0)).shape

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
    sol_d = implicitmidpoint(empty_d)
    sol_e = implicitmidpoint(empty_e)

    sol_f = implicitmidpoint(empty_f)
    sol_g = implicitmidpoint(empty_g)
    sol_h = implicitmidpoint(empty_h)
    sol_i = implicitmidpoint(empty_i)
    sol_j = implicitmidpoint(empty_j)

    #Create list of Earth frame trajectories
    EFsolutions = [sol_a, sol_b, sol_c, sol_d, sol_e, sol_f, sol_g, sol_h, sol_i, sol_j]

    print 'last element in sol_a is', sol_a[-1]
    print sol_a

###############################################################################################

    #Transform to wave frame by shifting x-coordinates

    shift = np.zeros_like(sol_a)
    shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

    #Define new x-shifted trajectories
    rel_sol_a = sol_a + shift #rel for relative motion
    rel_sol_b = sol_b + shift
    rel_sol_c = sol_c + shift
    rel_sol_d = sol_d + shift
    rel_sol_e = sol_e + shift

    rel_sol_f = sol_f + shift
    rel_sol_g = sol_g + shift
    rel_sol_h = sol_h + shift
    rel_sol_i = sol_i + shift
    rel_sol_j = sol_j + shift

    #Create list of Wave frame trajectories
    WFsolutions = [rel_sol_a, rel_sol_b, rel_sol_c, rel_sol_d, rel_sol_e, rel_sol_f, rel_sol_g, rel_sol_h, rel_sol_i, rel_sol_j]


###############################################################################################

    #Solve for parcel trajectories AGAIN but using Scipy's odeint - 'lsoda' integrator

    sol_a2 = odeint(velocity2, s0_a, t)
    sol_b2 = odeint(velocity2, s0_b, t)
    sol_c2 = odeint(velocity2, s0_c, t)
    sol_d2 = odeint(velocity2, s0_d, t)
    sol_e2 = odeint(velocity2, s0_e, t)

    sol_f2 = odeint(velocity2, s0_f, t)
    sol_g2 = odeint(velocity2, s0_g, t)
    sol_h2 = odeint(velocity2, s0_h, t)
    sol_i2 = odeint(velocity2, s0_i, t)
    sol_j2 = odeint(velocity2, s0_j, t)

    EFsolutions2 = [sol_a2, sol_b2, sol_c2, sol_d2, sol_e2, sol_f2, sol_g2, sol_h2, sol_i2, sol_j2]

    #Define new x-shifted trajectories
    rel_sol_a2 = sol_a2 + shift #rel for relative motion
    rel_sol_b2 = sol_b2 + shift
    rel_sol_c2 = sol_c2 + shift
    rel_sol_d2 = sol_d2 + shift
    rel_sol_e2 = sol_e2 + shift

    rel_sol_f2 = sol_f2 + shift
    rel_sol_g2 = sol_g2 + shift
    rel_sol_h2 = sol_h2 + shift
    rel_sol_i2 = sol_i2 + shift
    rel_sol_j2 = sol_j2 + shift

    #Create list of Wave frame trajectories
    WFsolutions2 = [rel_sol_a2, rel_sol_b2, rel_sol_c2, rel_sol_d2, rel_sol_e2, rel_sol_f2, rel_sol_g2, rel_sol_h2, rel_sol_i2, rel_sol_j2]

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

    #Define normals to project trajectories onto surfaces

    #Normal of theta surface
    normal = np.array((0, dthetady*L, dthetadz*H))
    normalhat = normal/(la.norm(normal))

    #Normal to surface at half the slope of theta surface
    #Maximum energy release for baroclinic instability typically within this half-slope surface
    normal_i = np.array((0, dthetady*L, 2*dthetadz*H)) #i for instability
    normal_ihat = normal_i/(la.norm(normal_i))

    print 'normal to theta surfaces is:', normalhat
    print 'normal to half-slope surface is:', normal_ihat

    sol_a_p = np.zeros_like(rel_sol_a) #p for projected
    sol_b_p = np.zeros_like(rel_sol_a)
    sol_c_p = np.zeros_like(rel_sol_a)
    sol_d_p = np.zeros_like(rel_sol_a)
    sol_e_p = np.zeros_like(rel_sol_a)
    sol_f_p = np.zeros_like(rel_sol_a)
    sol_g_p = np.zeros_like(rel_sol_a)
    sol_h_p = np.zeros_like(rel_sol_a)
    sol_i_p = np.zeros_like(rel_sol_a)
    sol_j_p = np.zeros_like(rel_sol_a)

    sol_a_pi = np.zeros_like(rel_sol_a) #i for instability
    sol_b_pi = np.zeros_like(rel_sol_a)
    sol_c_pi = np.zeros_like(rel_sol_a)
    sol_d_pi = np.zeros_like(rel_sol_a)
    sol_e_pi = np.zeros_like(rel_sol_a)
    sol_f_pi = np.zeros_like(rel_sol_a)
    sol_g_pi = np.zeros_like(rel_sol_a)
    sol_h_pi = np.zeros_like(rel_sol_a)
    sol_i_pi = np.zeros_like(rel_sol_a)
    sol_j_pi = np.zeros_like(rel_sol_a)

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

    #Project WF solutions onto sheet at half-slope of theta surface
    for i in range(len(t)):
        sol_a_pi[i] = rel_sol_a[i] - np.dot(rel_sol_a[i], normal_ihat)*normal_ihat
        sol_b_pi[i] = rel_sol_b[i] - np.dot(rel_sol_b[i], normal_ihat)*normal_ihat
        sol_c_pi[i] = rel_sol_c[i] - np.dot(rel_sol_c[i], normal_ihat)*normal_ihat
        sol_d_pi[i] = rel_sol_d[i] - np.dot(rel_sol_d[i], normal_ihat)*normal_ihat
        sol_e_pi[i] = rel_sol_e[i] - np.dot(rel_sol_e[i], normal_ihat)*normal_ihat
        sol_f_pi[i] = rel_sol_f[i] - np.dot(rel_sol_f[i], normal_ihat)*normal_ihat
        sol_g_pi[i] = rel_sol_g[i] - np.dot(rel_sol_g[i], normal_ihat)*normal_ihat
        sol_h_pi[i] = rel_sol_h[i] - np.dot(rel_sol_h[i], normal_ihat)*normal_ihat
        sol_i_pi[i] = rel_sol_i[i] - np.dot(rel_sol_i[i], normal_ihat)*normal_ihat
        sol_j_pi[i] = rel_sol_j[i] - np.dot(rel_sol_j[i], normal_ihat)*normal_ihat

    projected_solutions = [sol_a_p, sol_b_p, sol_c_p, sol_d_p, sol_e_p, sol_f_p, sol_g_p, sol_h_p, sol_i_p, sol_j_p]
    projected_solutions_i = [sol_a_pi, sol_b_pi, sol_c_pi, sol_d_pi, sol_e_pi, sol_f_pi, sol_g_pi, sol_h_pi, sol_i_pi, sol_j_pi]

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
        ax1.plot(i[:-1,0], i[:-1,1], i[:-1,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    ax2.set_title('x-z plane', fontsize=10)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    for i in EFsolutions:
        ax2.plot(i[:-1,0], i[:-1,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    ax3.set_title('y-z plane', fontsize=10)
    ax3.set_xlabel('y (L)')
    ax3.set_ylabel('z (H)')
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    for i in EFsolutions:
        ax3.plot(i[:-1,1], i[:-1,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    ax4.set_title('x-y plane', fontsize=10)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')
    for i in EFsolutions:
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
        ax5.plot(i[:-1,0], i[:-1,1], i[:-1,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    ax6.set_title('x-z plane', fontsize=10)
    ax6.set_xlabel('x (L)')
    ax6.set_ylabel('z (H)')
    for i in WFsolutions:
        ax6.plot(i[:-1,0], i[:-1,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    ax7.set_title('y-z plane', fontsize=10)
    ax7.set_xlabel('y (L)')
    ax7.set_ylabel('z (H)')
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    for i in WFsolutions:
        ax7.plot(i[:-1,1], i[:-1,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    ax8.set_title('x-y plane', fontsize=10)
    ax8.set_xlabel('x (L)')
    ax8.set_ylabel('y (L)')
    for i in WFsolutions:
        ax8.plot(i[:-1,0], i[:-1,1])

###############################################################################################

    #Plot of mean parcel displacement with time

    meanxseparation = np.zeros_like(t)
    meanyseparation = np.zeros_like(t)
    meanzseparation = np.zeros_like(t)

    for i in range(len(t)):
        meanxseparation[i]=(1./len(WFsolutions))*np.sum((j[i,0] - j[0,0])**2 for j in WFsolutions)
        meanyseparation[i]=(1./len(WFsolutions))*np.sum((j[i,1] - j[0,1])**2 for j in WFsolutions)
        meanzseparation[i]=(1./len(WFsolutions))*np.sum((j[i,2] - j[0,2])**2 for j in WFsolutions)

    displacement = np.zeros_like(t)
    for i in range(len(t)):
        displacement[i] = np.sqrt(meanxseparation[i] + meanyseparation[i] + meanzseparation[i]) #Note lack of 'squares' as separations are already distances squared


    #Do the same for 'lsoda' integrated trajectories
    meanxseparation2 = np.zeros_like(t)
    meanyseparation2 = np.zeros_like(t)
    meanzseparation2 = np.zeros_like(t)

    for i in range(len(t)):
        meanxseparation2[i]=(1./len(WFsolutions2))*np.sum((j[i,0] - j[0,0])**2 for j in WFsolutions2)
        meanyseparation2[i]=(1./len(WFsolutions2))*np.sum((j[i,1] - j[0,1])**2 for j in WFsolutions2)
        meanzseparation2[i]=(1./len(WFsolutions2))*np.sum((j[i,2] - j[0,2])**2 for j in WFsolutions2)

    displacement2 = np.zeros_like(t)
    for i in range(len(t)):
        displacement2[i] = np.sqrt(meanxseparation2[i] + meanyseparation2[i] + meanzseparation2[i]) #Note lack of 'squares' as separations are already distances squared

    fig3 = plt.figure()
    fig3.suptitle('Evolution of mean parcel displacement')
    ax9 = fig3.add_subplot(111)
    ax9.set_xlabel('time (T)')
    ax9.set_ylabel('displacement')
    ax9.plot(t[:-1], displacement[:-1], label='implicit midpoint')
    ax9.plot(t, displacement2, label='lsoda')
    ax9.legend(loc='upper left')

###############################################################################################

    #Plot of oscillation growth with time

    fig4 = plt.figure()
    fig4.suptitle('Growth of oscillations with time')

    ax10 = fig4.add_subplot(311)
    ax10.set_title('Zonal extent', fontsize=10)
    ax10.set_ylabel('(x-x0)^2')

    ax11 = fig4.add_subplot(312)
    ax11.set_title('Meridional extent', fontsize=10)
    ax11.set_ylabel('(y-y0)^2')

    ax12 = fig4.add_subplot(313)
    ax12.set_title('Vertical extent', fontsize=10)
    ax12.set_xlabel('time (T)')
    ax12.set_ylabel('(z-z0)^2')

    ax10.plot(t[:-1], meanxseparation[:-1], label='implicit midpoint')
    ax11.plot(t[:-1], meanyseparation[:-1], label='implicit midpoint')
    ax12.plot(t[:-1], meanzseparation[:-1], label='implicit midpoint')

    ax10.plot(t, meanxseparation2, label='lsoda')
    ax11.plot(t, meanyseparation2, label='lsoda')
    ax12.plot(t, meanzseparation2, label='lsoda')

    ax10.legend(loc='upper left')

###############################################################################################

    #WF isentropic surface projection figure

    fig5 = plt.figure()
    fig5.suptitle('Relative motion within surface at half-slope to isentropes')
    ax13 = fig5.add_subplot(111, projection = '3d')
    ax13.set_xlabel('x (L)')
    ax13.set_ylabel('y (L)')
    ax13.set_zlabel('z (H)')
    for j in projected_solutions_i:
        ax13.plot(j[:-1,0], j[:-1,1], j[:-1,2])

    times = np.arange(0,len(t),50)

    #Add a dot along trajectories every 50 timesteps to indicate time evolution
    for i in range(len(times)):
        for j in projected_solutions_i:
            ax13.scatter(j[times[i],0], j[times[i],1], j[times[i],2], marker='o', c='black', s=8)


###############################################################################################

    #Plotting total theta (i.e. background theta plus thetaprime)
    #For one parcel (a) with time

    absolutetheta_a = np.zeros_like(t)
    absolutetheta_a2 = np.zeros_like(t)

    for i in range(len(t)):
        absolutetheta_a[i] = theta(y=rel_sol_a[i,1], z=rel_sol_a[i,2]) + thetaprime(x=rel_sol_a[i,0], y=rel_sol_a[i,1], z=rel_sol_a[i,2], t=t[i])
        absolutetheta_a2[i] = theta(y=rel_sol_a2[i,1], z=rel_sol_a2[i,2]) + thetaprime(x=rel_sol_a2[i,0], y=rel_sol_a2[i,1], z=rel_sol_a2[i,2], t=t[i])

    fig6 = plt.figure()
    ax14 = fig6.add_subplot(111)
    ax14.set_xlabel('time (T)')
    ax14.set_ylabel('Absolute theta (K)')
    ax14.plot(t[:-1], absolutetheta_a[:-1], color='black', label='implicit midpoint')
    ax14.plot(t, absolutetheta_a2, color='blue', label='lsoda')
    ax14.axhline(y=absolutetheta_a[0], color='black',ls='dotted')
    ax14.legend(loc='lower left')

    plt.show()

###############################################################################################

#Run the programme

main(0,30,0.9)
