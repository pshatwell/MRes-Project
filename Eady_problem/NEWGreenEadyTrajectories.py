from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import numpy.linalg as la
import odespy

from Eadyinfo_oc import *

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

def main(start, stop, xpos=0, ypos=0.5):

    def velocity(s,t):
        x,y,z = s
        #dsdt = np.array((uprime(x=x,y=y,z=z,t=t), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)))
        dsdt = np.array((uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)))
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = start
    tmax = stop
    nt = 500

    time_points = np.linspace(tmin, tmax, nt)

    dt = (tmax - tmin)/nt

    print '\nnumber of time_points is', nt
    print 'timestep dt is', dt

###############################################################################################

    #Define numerical integrator and initial conditions

    def props(cls):
        return [i for i in cls.__dict__.keys() if i[:1] != '_']

    properties = props(odespy)
    print '\nodespy integrators are:', properties

    #Define initial positions of parcels
    #Note: zpos=0.5 determines steering level, i.e. where c.real matches mean flow speed

    #xshift = np.pi/(k*L) #nondimensional shift of half a wavelength
    xshift = 2*np.pi/(k*L) #nondimensional shift of a wavelength

    gamma = (dthetady*L)/(dthetadz*H)
    print '\ngamma is:', gamma

    def zpos(y): #enforce parcels to lie on line at half-slope to isentropes
        return -(gamma/2.)*y + 0.5

    def zpos2(y):
        return -gamma*y + 0.5

    print '(ypos, z(ypos)) is:', (ypos, zpos(ypos))
    print '(-ypos, z(-ypos)) is:', (-ypos, zpos(-ypos))

    #Set of 5 'warm' parcels (positive y, near top)
    s0_a = np.array((xpos, ypos, zpos(ypos)))
    s0_b = np.array((xpos+0.5, ypos, zpos(ypos)))
    s0_c = np.array((xpos+1, ypos, zpos(ypos)))
    s0_d = np.array((xpos+1.5, ypos, zpos(ypos)))
    s0_e = np.array((xpos+2, ypos, zpos(ypos)))

    #Set of 5 'cold' parcels (negative y, near bottom)
    s0_f = np.array((xpos+xshift, -ypos, zpos(-ypos)))
    s0_g = np.array((xpos+xshift+0.5, -ypos, zpos(-ypos)))
    s0_h = np.array((xpos+xshift+1, -ypos, zpos(-ypos)))
    s0_i = np.array((xpos+xshift+1.5, -ypos, zpos(-ypos)))
    s0_j = np.array((xpos+xshift+2, -ypos, zpos(-ypos)))


    solver_a = odespy.MidpointImplicit(velocity)
    solver_a.set_initial_condition(s0_a)

    solver_b = odespy.MidpointImplicit(velocity)
    solver_b.set_initial_condition(s0_b)
    solver_c = odespy.MidpointImplicit(velocity)
    solver_c.set_initial_condition(s0_c)
    solver_d = odespy.MidpointImplicit(velocity)
    solver_d.set_initial_condition(s0_d)
    solver_e = odespy.MidpointImplicit(velocity)
    solver_e.set_initial_condition(s0_e)


    solver_f = odespy.MidpointImplicit(velocity)
    solver_f.set_initial_condition(s0_f)

    solver_g = odespy.MidpointImplicit(velocity)
    solver_g.set_initial_condition(s0_g)
    solver_h = odespy.MidpointImplicit(velocity)
    solver_h.set_initial_condition(s0_h)
    solver_i = odespy.MidpointImplicit(velocity)
    solver_i.set_initial_condition(s0_i)
    solver_j = odespy.MidpointImplicit(velocity)
    solver_j.set_initial_condition(s0_j)


###############################################################################################

    #Solve for parcel trajectories

    sol_a, t = solver_a.solve(time_points)

    sol_b, t = solver_b.solve(time_points)
    sol_c, t = solver_c.solve(time_points)
    sol_d, t = solver_d.solve(time_points)
    sol_e, t = solver_e.solve(time_points)

    sol_f, t = solver_f.solve(time_points)

    sol_g, t = solver_g.solve(time_points)
    sol_h, t = solver_h.solve(time_points)
    sol_i, t = solver_i.solve(time_points)
    sol_j, t = solver_j.solve(time_points)

    #Create list of Earth frame trajectories
    EFsolutions = [sol_a, sol_b, sol_c, sol_d, sol_e, sol_f, sol_g, sol_h, sol_i, sol_j]
    #EFsolutions = [sol_a, sol_f]
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
    #WFsolutions = [rel_sol_a, rel_sol_f]
###############################################################################################


    #Solve for parcel trajectories AGAIN but using Scipy's odeint - 'lsoda' integrator

    sol_a2 = odeint(velocity, s0_a, t)
    sol_b2 = odeint(velocity, s0_b, t)
    sol_c2 = odeint(velocity, s0_c, t)
    sol_d2 = odeint(velocity, s0_d, t)
    sol_e2 = odeint(velocity, s0_e, t)

    sol_f2 = odeint(velocity, s0_f, t)
    sol_g2 = odeint(velocity, s0_g, t)
    sol_h2 = odeint(velocity, s0_h, t)
    sol_i2 = odeint(velocity, s0_i, t)
    sol_j2 = odeint(velocity, s0_j, t)

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

    #Create matrices for background potential temperature distribution
    #and streamfunction perturbation

    time = start #time at which to evaluate arrays for Eady solution

    #Find extents for arrays to display graphically behind trajectories
    xmins = []
    for trajectory in WFsolutions:
        xmins.append(min(trajectory[:,0]))

    xmaxes = []
    for trajectory in WFsolutions:
        xmaxes.append(max(trajectory[:,0]))

    ymins = []
    for trajectory in WFsolutions:
        ymins.append(min(trajectory[:,1]))

    ymaxes = []
    for trajectory in WFsolutions:
        ymaxes.append(max(trajectory[:,1]))

    zmins = []
    for trajectory in WFsolutions:
        zmins.append(min(trajectory[:,2]))

    zmaxes = []
    for trajectory in WFsolutions:
        zmaxes.append(max(trajectory[:,2]))

    xmin = min(xmins)
    xmax = max(xmaxes)

    ymin = min(ymins)
    ymax = max(ymaxes)

    zmin = min(zmins)
    zmax = max(zmaxes)

    xvalues = np.linspace(xmin, xmax, 20)
    yvalues = np.linspace(ymin, ymax, 20)
    zvalues = np.linspace(zmin, zmax, 20)

    xlength = len(xvalues)
    ylength = len(yvalues)
    zlength = len(zvalues)

    theta_matrix = np.zeros((zlength,ylength))
    psiprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate theta matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            theta_matrix[i,j] = theta(y=yvalues[j], z=zvalues[i])

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                psiprime_matrix[i,j,m] = psiprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    print 'dthetady is:', dthetady*L
    print 'dthetadz is:', dthetadz*H

    #Make array for full perturbation velocity field
    velocityfield = np.zeros((zlength, ylength, xlength, 3))

    print 'shape of velocity field is:', velocityfield.shape

    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                velocityfield[i,j,m,0] = uprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time) + umeanflow(z=zvalues[i]) - c.real*(T/L) #Subtracting phase speed, so we're in wave frame
                velocityfield[i,j,m,1] = vprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)
                velocityfield[i,j,m,2] = wprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

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
    #projected_solutions = [sol_a_p, sol_f_p]
    projected_solutions_i = [sol_a_pi, sol_b_pi, sol_c_pi, sol_d_pi, sol_e_pi, sol_f_pi, sol_g_pi, sol_h_pi, sol_i_pi, sol_j_pi]
    #projected_solutions_i = [sol_a_pi, sol_f_pi]
###############################################################################################

    #Plot the trajectories

    ypoints = np.linspace(-1,1,50)

    #Define slices for background distributions
    zslice = int(round(zlength/2.))
    yslice = int(round(ylength/2.))
    xslice = int(round(xlength/2.))

    #Plot in the Earth frame

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Earth frame', fontsize='18')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)', fontsize='16')
    ax1.set_ylabel('y (L)', fontsize='16')
    ax1.set_zlabel('z (H)', fontsize='16')
    for i in EFsolutions:
        ax1.plot(i[:,0], i[:,1], i[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    #ax2.set_title('x-z plane', fontsize='14')
    ax2.set_xlabel('x (L)', fontsize='16')
    ax2.set_ylabel('z (H)', fontsize='16')
    for i in EFsolutions:
        ax2.plot(i[:,0], i[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    #ax3.set_title('y-z plane', fontsize='14')
    ax3.set_xlabel('y (L)', fontsize='16')
    ax3.set_ylabel('z (H)', fontsize='16')
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    ax3.plot(ypoints, zpos(ypoints), color='white', ls='dotted')
    for i in EFsolutions:
        ax3.plot(i[:,1], i[:,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    #ax4.set_title('x-y plane', fontsize='14')
    ax4.set_xlabel('x (L)', fontsize='16')
    ax4.set_ylabel('y (L)', fontsize='16')
    for i in EFsolutions:
        ax4.plot(i[:,0], i[:,1])

###############################################################################################

    #Plot in the Wave frame

    #Define grids for quiver plots
    X, Y = np.meshgrid(xvalues, yvalues)
    X1, Z1 = np.meshgrid(xvalues, zvalues)
    Y2, Z2 = np.meshgrid(yvalues, zvalues)

    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('Wave frame', fontsize='18')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)', fontsize='16')
    ax5.set_ylabel('y (L)', fontsize='16')
    ax5.set_zlabel('z (H)', fontsize='16')
    for i in WFsolutions:
        ax5.plot(i[:,0], i[:,1], i[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    #ax6.set_title('x-z plane', fontsize='14')
    ax6.set_xlabel('x (L)', fontsize='16')
    ax6.set_ylabel('z (H)', fontsize='16')
    #psicontour1 = ax6.contourf(psiprime_matrix[:,yslice,:], origin='lower', extent=[xmin, xmax, zmin, zmax], aspect='auto', cmap='Greys')
    #plt.colorbar(psicontour1)
    Q1 = ax6.quiver(X1,Z1, velocityfield[:,yslice,:,0], 100*velocityfield[:,yslice,:,2], headlength=3, headaxislength=3)
    for i in WFsolutions:
        ax6.plot(i[:,0], i[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    #ax7.set_title('y-z plane', fontsize='14')
    ax7.set_xlabel('y (L)', fontsize='16')
    ax7.set_ylabel('z (H)', fontsize='16')
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    ax7.plot(ypoints, zpos(ypoints), color='white', ls='dotted')
    #ax7.plot(ypoints, zpos2(ypoints), color='white', ls='dotted')
    for i in WFsolutions:
        ax7.plot(i[:,1], i[:,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    #ax8.set_title('x-y plane', fontsize='14')
    ax8.set_xlabel('x (L)', fontsize='16')
    ax8.set_ylabel('y (L)', fontsize='16')
    Q2 = ax8.quiver(X,Y, velocityfield[zslice,:,:,0], 10*velocityfield[zslice,:,:,1], headlength=3, headaxislength=3)
    #psicontour2 = ax8.contourf(psiprime_matrix[zslice,:,:], origin='lower', extent=[xmin, xmax, ymin, ymax], aspect='auto', cmap='Greys')
    #plt.colorbar(psicontour2)
    for i in WFsolutions:
        ax8.plot(i[:,0], i[:,1])

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
    #fig3.suptitle('Evolution of mean parcel displacement', fontsize='16')
    ax9 = fig3.add_subplot(111)
    ax9.set_xlabel('time (T)', fontsize='16')
    ax9.set_ylabel('Mean parcel displacement', fontsize='16')
    ax9.plot(t[:], displacement[:], label='implicit midpoint')
    #ax9.plot(t, displacement2, label='lsoda')
    #ax9.legend(loc='upper left')

###############################################################################################

    #Plot of oscillation growth with time

    fig4 = plt.figure()
    #fig4.suptitle('Growth of oscillations with time', fontsize='16')

    ax10 = fig4.add_subplot(311)
    #ax10.set_title('Zonal extent', fontsize='14')
    ax10.set_ylabel('(x-x0)^2', fontsize='16')

    ax11 = fig4.add_subplot(312)
    #ax11.set_title('Meridional extent', fontsize='14')
    ax11.set_ylabel('(y-y0)^2', fontsize='16')

    ax12 = fig4.add_subplot(313)
    #ax12.set_title('Vertical extent', fontsize='14')
    ax12.set_xlabel('time (T)', fontsize='16')
    ax12.set_ylabel('(z-z0)^2', fontsize='16')

    ax10.plot(t[:], meanxseparation[:], label='implicit midpoint')
    ax11.plot(t[:], meanyseparation[:], label='implicit midpoint')
    ax12.plot(t[:], meanzseparation[:], label='implicit midpoint')

    #ax10.plot(t, meanxseparation2, label='lsoda')
    #ax11.plot(t, meanyseparation2, label='lsoda')
    #ax12.plot(t, meanzseparation2, label='lsoda')

    #ax10.legend(loc='upper left')

###############################################################################################

    #WF isentropic surface projection figure

    fig5 = plt.figure()
    #fig5.suptitle('Relative motion within surface at half-slope to isentropes', fontsize='16')
    ax13 = fig5.add_subplot(111, projection = '3d')
    ax13.set_xlabel('x (L)', fontsize='16')
    ax13.set_ylabel('y (L)', fontsize='16')
    ax13.set_zlabel('z (H)', fontsize='16')
    for j in projected_solutions_i:
        ax13.plot(j[:,0], j[:,1], j[:,2])

    times = np.arange(0,len(time_points),100)

    #Add a dot along trajectories every 100 timesteps to indicate time evolution
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
    ax14.set_xlabel('time (T)', fontsize='16')
    ax14.set_ylabel('Potential temperature (K)', fontsize='16')
    ax14.plot(t[:], absolutetheta_a[:], color='black', label='implicit midpoint')
    #ax14.plot(t, absolutetheta_a2, color='blue', label='lsoda')
    ax14.axhline(y=absolutetheta_a[0], color='black',ls='dotted')
    #ax14.legend(loc='lower left')

    plt.show()

###############################################################################################

#Run the programme

main(start=15, stop=60, xpos=0, ypos=0.2)
