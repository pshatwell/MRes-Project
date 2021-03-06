from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as la
import odespy

from Eadyinfo_oc import *

###############################################################################################

#Print physical parameters for simulation run

print 'f0 is (s^-1)', f0
print 'N is (s^-1)', N
print 'U_0 is (ms^-1)', U_0

print '\nH is (m)', H
print 'L is (m)', L
print 'T is (s)', T

print '\nk is (m^-1)', k
print 'l is (m^-1)', l
print 'mu is', mu

print '\nc is (ms^-1)', c
print 'sigma is (s^-1)', sigma

###############################################################################################

def main(start=0, stop=30, xpos=0, ypos=0.5, zpos=0.5, showvel=False, showgrowth=False, showgreen=False, showtheta=False):

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

    xshift = np.pi/(k*L) #nondimensional shift of half a wavelength

    s0_a = np.array((xpos, ypos, zpos))
    s0_b = np.array((xpos, 0, zpos))
    s0_c = np.array((xpos, -ypos, zpos))
    #s0_c = np.array((xpos + xshift, -ypos, zpos))

    #solver_a = odespy.ForwardEuler(velocity)
    solver_a = odespy.MidpointImplicit(velocity)
    #solver_a = odespy.lsoda_scipy(velocity)
    #solver_a = odespy.RK4(velocity)
    solver_a.set_initial_condition(s0_a)

    #solver_b = odespy.ForwardEuler(velocity)
    solver_b = odespy.MidpointImplicit(velocity)
    #solver_b = odespy.lsoda_scipy(velocity)
    #solver_b = odespy.RK4(velocity)
    solver_b.set_initial_condition(s0_b)

    #solver_c = odespy.ForwardEuler(velocity)
    solver_c = odespy.MidpointImplicit(velocity)
    #solver_c = odespy.lsoda_scipy(velocity)
    #solver_c = odespy.RK4(velocity)
    solver_c.set_initial_condition(s0_c)

    #Solve for parcel trajectories

    s_a, t = solver_a.solve(time_points)
    s_b, t = solver_b.solve(time_points)
    s_c, t = solver_c.solve(time_points)

    #Create list of Earth frame trajectories
    EFsolutions = np.array((s_a, s_b, s_c))

    print '\ntype of s_a is:', type(s_a)
    print 'shape of s_a is:', s_a.shape

###############################################################################################

    #Transform to wave frame by shifting x-coordinates

    shift = np.zeros_like(s_a)
    shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

    #Define new x-shifted trajectories
    rel_s_a = s_a + shift #rel for relative motion
    rel_s_b = s_b + shift
    rel_s_c = s_c + shift

    #Create list of Wave frame trajectories
    WFsolutions = np.array((rel_s_a, rel_s_b, rel_s_c))

###############################################################################################

    #Define normals to project trajectories onto surfaces

    #Normal of theta surface
    normal = np.array((0, dthetady*L, dthetadz*H))
    normalhat = normal/(la.norm(normal))

    #Normal to surface at half the slope of theta surface
    #Maximum energy release for baroclinic instability typically within this half-slope surface
    normal_i = np.array((0, dthetady*L, 2*dthetadz*H)) #i for instability
    normal_ihat = normal_i/(la.norm(normal_i))

    print '\nnormal to theta surfaces is:', normalhat
    print 'normal to half-slope surface is:', normal_ihat

    s_a_p = np.zeros_like(rel_s_a) #p for projected
    s_b_p = np.zeros_like(rel_s_a)
    s_c_p = np.zeros_like(rel_s_a)

    print 'shape of s_a_p is:', s_a_p.shape

    s_a_pi = np.zeros_like(rel_s_a) #i for instability
    s_b_pi = np.zeros_like(rel_s_a)
    s_c_pi = np.zeros_like(rel_s_a)

    #Project WF solutions onto theta surface
    for i in range(len(time_points)):
        s_a_p[i] = rel_s_a[i] - np.dot(rel_s_a[i], normalhat)*normalhat
        s_b_p[i] = rel_s_b[i] - np.dot(rel_s_b[i], normalhat)*normalhat
        s_c_p[i] = rel_s_c[i] - np.dot(rel_s_c[i], normalhat)*normalhat

    #Project WF solutions onto sheet at half-slope of theta surface
    for i in range(len(time_points)):
        s_a_pi[i] = rel_s_a[i] - np.dot(rel_s_a[i], normal_ihat)*normal_ihat
        s_b_pi[i] = rel_s_b[i] - np.dot(rel_s_b[i], normal_ihat)*normal_ihat
        s_c_pi[i] = rel_s_c[i] - np.dot(rel_s_c[i], normal_ihat)*normal_ihat

    projected_solutions = [s_a_p, s_b_p, s_c_p]
    projected_solutions_i = [s_a_pi, s_b_pi, s_c_pi]

###############################################################################################

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

    upmedian = np.median(rel_s_a[:,2]) #This should drift up
    downmedian = np.median(rel_s_c[:,2]) #This should drift down

    print 'upmedian is:', upmedian
    print 'downmedian is:', downmedian

    print 'upmedian - zpos is:', upmedian - zpos
    print 'downmedian - zpos is:', downmedian - zpos

    xvalues = np.linspace(xmin, xmax, 80)
    yvalues = np.linspace(ymin, ymax, 80)
    zvalues = np.linspace(zmin, zmax, 80)

    xlength = len(xvalues)
    ylength = len(yvalues)
    zlength = len(zvalues)

###############################################################################################

    #Define arrays for Eady solution
    #Solutions are functions in Eadyinfo script

    time = start #time at which to evaluate arrays for Eady solution

    theta_matrix = np.zeros((zlength, ylength))

    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            theta_matrix[i,j] = theta(y=yvalues[j], z=zvalues[i])

    print '\ndthetady is:', dthetady*L
    print 'dthetadz is:', dthetadz*H

    gradient = (dthetadz*H)/(dthetady*L)
    print 'gradient is:', gradient

    #Create empty matrix for thetaprime
    thetaprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate thetaprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                thetaprime_matrix[i,j,m] = thetaprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    #Create empty matrix for wprime
    wprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                wprime_matrix[i,j,m] = wprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    #Create empty matrix for streamfunction perturbation
    psiprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                psiprime_matrix[i,j,m] = psiprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)


###############################################################################################

    #Plot the trajectories

    #Plot in the Earth frame

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    fig1.suptitle('Earth frame', fontsize='16')
    ax1 = fig1.add_subplot(221, projection = '3d')
    ax1.set_xlabel('x (L)', fontsize='14')
    ax1.set_ylabel('y (L)', fontsize='14')
    ax1.set_zlabel('z (H)', fontsize='14')
    for i in EFsolutions:
        ax1.plot(i[:,0], i[:,1], i[:,2])

    #Projection in x-z plane
    ax2 = fig1.add_subplot(222)
    #ax2.set_title('x-z plane', fontsize='14')
    ax2.set_xlabel('x (L)', fontsize='14')
    ax2.set_ylabel('z (H)', fontsize='14')
    for i in EFsolutions:
        ax2.plot(i[:,0], i[:,2])

    #Projection in y-z plane
    ax3 = fig1.add_subplot(223)
    #ax3.set_title('y-z plane', fontsize='14')
    ax3.set_xlabel('y (L)', fontsize='14')
    ax3.set_ylabel('z (H)', fontsize='14')
    isentropes1 = ax3.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes1)
    for i in EFsolutions:
        ax3.plot(i[:,1], i[:,2])

    #Projection in x-y plane
    ax4 = fig1.add_subplot(224)
    #ax4.set_title('x-y plane', fontsize='14')
    ax4.set_xlabel('x (L)', fontsize='14')
    ax4.set_ylabel('y (L)', fontsize='14')
    for i in EFsolutions:
        ax4.plot(i[:,0], i[:,1])

###############################################################################################

    #Plot in the Wave frame

    fig2 = plt.figure()
    plt.set_cmap('inferno')
    fig2.suptitle('Wave frame', fontsize='16')
    ax5 = fig2.add_subplot(221, projection = '3d')
    ax5.set_xlabel('x (L)', fontsize='14')
    ax5.set_ylabel('y (L)', fontsize='14')
    ax5.set_zlabel('z (H)', fontsize='14')
    for i in WFsolutions:
        ax5.plot(i[:,0], i[:,1], i[:,2])

    #Projection in x-z plane
    ax6 = fig2.add_subplot(222)
    #ax6.set_title('x-z plane', fontsize='14')
    ax6.set_xlabel('x (L)', fontsize='14')
    ax6.set_ylabel('z (H)', fontsize='14')
    for i in WFsolutions:
        ax6.plot(i[:,0], i[:,2])

    #Projection in y-z plane
    ax7 = fig2.add_subplot(223)
    #ax7.set_title('y-z plane', fontsize='14')
    ax7.set_xlabel('y (L)', fontsize='14')
    ax7.set_ylabel('z (H)', fontsize='14')
    isentropes2 = ax7.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes2)
    for i in WFsolutions:
        ax7.plot(i[:,1], i[:,2])

    #Projection in x-y plane
    ax8 = fig2.add_subplot(224)
    #ax8.set_title('x-y plane', fontsize='14')
    ax8.set_xlabel('x (L)', fontsize='14')
    ax8.set_ylabel('y (L)', fontsize='14')
    for i in WFsolutions:
        ax8.plot(i[:,0], i[:,1])

###############################################################################################

    '''WALLACE 78 FIGURES'''

    if sigma != 0:
        efoldtime = 1./(sigma*T)
        print 'efoldtime is', efoldtime

    #Create times for plotting dots along trajectories every e-folding time
    times = np.arange(start, len(time_points), int((efoldtime*len(time_points)/stop)))

    #Define slices for background distributions
    zslice = int(round(zlength*zpos))
    yslice = int(round(ylength/2.))

    #Attempt at 'Figures 1 and 2' from Wallace 78 paper
    Wallacefig12 = plt.figure()

    W12ax = Wallacefig12.add_subplot(211)
    W12ax.set_ylabel('Latitude (L)', fontsize='15')
    #W12ax.set_xlabel('Longitude (L)', fontsize='15')
    plt.set_cmap('Reds')
    xy_contour_theta = W12ax.contourf(thetaprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_theta)
    #xy_contour_psi1 = W12ax.contourf(psiprime_matrix[zslice,:,:], cmap='viridis', origin='lower', aspect='auto', extent=[xmin,xmax,ymin,ymax])
    #plt.colorbar(xy_contour_psi1)
    #plt.clabel(xy_contour_psi1, inline=1, fontsize=10)
    for i in WFsolutions:
        W12ax.plot(i[:,0], i[:,1], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            W12ax.scatter(j[times[i],0], j[times[i],1], marker='o', c='black', s=9)


    W12ax2 = Wallacefig12.add_subplot(212)
    W12ax2.set_xlabel('Longitude (L)', fontsize='15')
    W12ax2.set_ylabel('Latitude (L)', fontsize='15')
    plt.set_cmap('Blues')
    xy_contour_w = W12ax2.contourf(wprime_matrix[zslice,:,:], aspect='auto', extent = [xmin,xmax,ymin,ymax])
    plt.colorbar(xy_contour_w)
    #xy_contour_psi2 = W12ax2.contour(psiprime_matrix[zslice,:,:], colors='white', origin='lower', aspect='auto', extent=[xmin,xmax,ymin,ymax])
    #plt.clabel(xy_contour_psi2, inline=1, fontsize=10)
    for i in WFsolutions:
        W12ax2.plot(i[:,0], i[:,1], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            W12ax2.scatter(j[times[i],0], j[times[i],1], marker='o', c='black', s=9)


    #'Figures 1 and 2' in the x-z plane

    Wallacefig12xz = plt.figure()

    W12ax3 = Wallacefig12xz.add_subplot(211)
    W12ax3.set_ylabel('Height (H)', fontsize='15')
    plt.set_cmap('Reds')
    xz_contour_theta = W12ax3.contourf(thetaprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    #xz_contour_theta = W12ax3.contourf(psiprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_theta)
    for i in WFsolutions:
        W12ax3.plot(i[:,0], i[:,2], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            W12ax3.scatter(j[times[i],0], j[times[i],2], marker='o', c='black', s=9)


    W12ax4 = Wallacefig12xz.add_subplot(212)
    W12ax4.set_xlabel('Longitude (L)', fontsize='15')
    W12ax4.set_ylabel('Height (H)', fontsize='15')
    plt.set_cmap('Blues')
    xz_contour_w = W12ax4.contourf(wprime_matrix[:,yslice,:], aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_w)
    for i in WFsolutions:
        W12ax4.plot(i[:,0], i[:,2], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            W12ax4.scatter(j[times[i],0], j[times[i],2], marker='o', c='black', s=9)


    #Attempt for 'Figure 3' from Wallace 78 paper

    gamma = (dthetady*L)/(dthetadz*H)
    print '\ngamma is:', gamma

    def zposition(y): #Line for half slope of isentropes
        return -(gamma/2.)*y + zpos

    ypoints = np.linspace(-1,1,50)

    Wallacefig3 = plt.figure()
    W3ax = Wallacefig3.add_subplot(111)
    W3ax.set_xlabel('Latitude (L)', fontsize='15')
    W3ax.set_ylabel('Height (H)', fontsize='15')
    isentropes = W3ax.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto', cmap='inferno')
    plt.colorbar(isentropes)
    W3ax.plot(ypoints, zposition(ypoints), color='white', ls='dotted')
    W3ax.axhline(y=zpos, color='white', ls='dashed')
    #W3ax.plot(rel_s_b[:,1], rel_s_b[:,2], lw='1.5', color='black')
    for i in WFsolutions:
        W3ax.plot(i[:,1], i[:,2], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        #W3ax.scatter(rel_s_b[times[i],1], rel_s_b[times[i],2], marker='o', c='black', s=9)
        for j in WFsolutions:
            W3ax.scatter(j[times[i],1], j[times[i],2], marker='o', c='black', s=9)

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

    disfig = plt.figure()
    #disfig.suptitle('Evolution of mean parcel displacement', fontsize='16')
    disax = disfig.add_subplot(111)
    disax.set_xlabel('time (T)', fontsize='15')
    disax.set_ylabel('Mean parcel displacement', fontsize='15')
    disax.plot(t[:], displacement[:], label='Implicit Midpoint')
    #disax.legend(loc='upper left')

###############################################################################################

    if showgrowth==True:
        #Plot of oscillation growth with time

        growthfig = plt.figure()
        #growthfig.suptitle('Growth of oscillations with time', fontsize='16')

        growthax1 = growthfig.add_subplot(311)
        #growthax1.set_title('Zonal extent', fontsize='14')
        growthax1.set_ylabel('(x-x0)^2', fontsize='15')

        growthax2 = growthfig.add_subplot(312)
        #growthax2.set_title('Meridional extent', fontsize='14')
        growthax2.set_ylabel('(y-y0)^2', fontsize='15')

        growthax3 = growthfig.add_subplot(313)
        #growthax3.set_title('Vertical extent', fontsize='14')
        growthax3.set_xlabel('time (T)', fontsize='15')
        growthax3.set_ylabel('(z-z0)^2', fontsize='15')

        growthax1.plot(t[:], meanxseparation[:], label='Implicit Midpoint')
        growthax2.plot(t[:], meanyseparation[:], label='Implicit Midpoint')
        growthax3.plot(t[:], meanzseparation[:], label='Implicit Midpoint')

        #growthax1.legend(loc='upper left')

###############################################################################################

    if showvel==True:

        #Plot of velocity perturbations with time

        uprime_a = np.zeros_like(time_points)
        uprime_b = np.zeros_like(time_points)
        uprime_c = np.zeros_like(time_points)

        vprime_a = np.zeros_like(time_points)
        vprime_b = np.zeros_like(time_points)
        vprime_c = np.zeros_like(time_points)

        wprime_a = np.zeros_like(time_points)
        wprime_b = np.zeros_like(time_points)
        wprime_c = np.zeros_like(time_points)

        for i in range(len(time_points)):
            uprime_a[i] = uprime(x=WFsolutions[0,i,0], y=WFsolutions[0,i,1], z=WFsolutions[0,i,2], t=time_points[i])
            uprime_b[i] = uprime(x=WFsolutions[1,i,0], y=WFsolutions[1,i,1], z=WFsolutions[1,i,2], t=time_points[i])
            uprime_c[i] = uprime(x=WFsolutions[2,i,0], y=WFsolutions[2,i,1], z=WFsolutions[2,i,2], t=time_points[i])

            vprime_a[i] = vprime(x=WFsolutions[0,i,0], y=WFsolutions[0,i,1], z=WFsolutions[0,i,2], t=time_points[i])
            vprime_b[i] = vprime(x=WFsolutions[1,i,0], y=WFsolutions[1,i,1], z=WFsolutions[1,i,2], t=time_points[i])
            vprime_c[i] = vprime(x=WFsolutions[2,i,0], y=WFsolutions[2,i,1], z=WFsolutions[2,i,2], t=time_points[i])

            wprime_a[i] = wprime(x=WFsolutions[0,i,0], y=WFsolutions[0,i,1], z=WFsolutions[0,i,2], t=time_points[i])
            wprime_b[i] = wprime(x=WFsolutions[1,i,0], y=WFsolutions[1,i,1], z=WFsolutions[1,i,2], t=time_points[i])
            wprime_c[i] = wprime(x=WFsolutions[2,i,0], y=WFsolutions[2,i,1], z=WFsolutions[2,i,2], t=time_points[i])

        velocityfig = plt.figure()
        uprimeax = velocityfig.add_subplot(311)
        uprimeax.set_ylabel('u prime', fontsize='15')
        uprimeax.plot(t[:], uprime_a[:], label='a')
        uprimeax.plot(t[:], uprime_b[:], label='b')
        uprimeax.plot(t[:], uprime_c[:], label='c')
        uprimeax.axhline(y = uprime_a[0], color='black', ls='dotted')
        uprimeax.legend(loc='upper left')

        vprimeax = velocityfig.add_subplot(312)
        vprimeax.set_ylabel('v prime', fontsize='15')
        vprimeax.plot(t[:], vprime_a[:], label='a')
        vprimeax.plot(t[:], vprime_b[:], label='b')
        vprimeax.plot(t[:], vprime_c[:], label='c')
        vprimeax.axhline(y = vprime_a[0], color='black', ls='dotted')

        wprimeax = velocityfig.add_subplot(313)
        wprimeax.set_xlabel('time (T)', fontsize='15')
        wprimeax.set_ylabel('w prime', fontsize='15')
        wprimeax.plot(t[:], wprime_a[:], label='a')
        wprimeax.plot(t[:], wprime_b[:], label='b')
        wprimeax.plot(t[:], wprime_c[:], label='c')
        wprimeax.axhline(y = wprime_a[0], color='black', ls='dotted')

###############################################################################################

    if showgreen==True:
        #WF isentropic surface projection figure

        greenfig = plt.figure()
        greenfig.suptitle('Relative motion within surface at half-slope to isentropes', fontsize='16')
        greenax = greenfig.add_subplot(111, projection = '3d')
        greenax.set_xlabel('x (L)', fontsize='14')
        greenax.set_ylabel('y (L)', fontsize='14')
        greenax.set_zlabel('z (H)', fontsize='14')
        for j in projected_solutions_i:
            greenax.plot(j[:,0], j[:,1], j[:,2])

        dot_times = np.arange(0,len(time_points),100)

        #Add a dot along trajectories every 100 timesteps to indicate time evolution
        for i in range(len(dot_times)):
            for j in projected_solutions_i:
                greenax.scatter(j[dot_times[i],0], j[dot_times[i],1], j[dot_times[i],2], marker='o', c='black', s=8)

###############################################################################################

    if showtheta==True:
        #Plotting total theta (i.e. background theta plus thetaprime)
        #For one parcel (a) with time

        absolutetheta_a = np.zeros_like(time_points)
        absolutetheta_b = np.zeros_like(time_points)
        absolutetheta_c = np.zeros_like(time_points)
        thetaprime_a = np.zeros_like(time_points)

        for i in range(len(time_points)):
            absolutetheta_a[i] = theta(y=rel_s_a[i,1], z=rel_s_a[i,2]) + thetaprime(x=rel_s_a[i,0], y=rel_s_a[i,1], z=rel_s_a[i,2], t=time_points[i])
            absolutetheta_b[i] = theta(y=rel_s_b[i,1], z=rel_s_b[i,2]) + thetaprime(x=rel_s_b[i,0], y=rel_s_b[i,1], z=rel_s_b[i,2], t=time_points[i])
            absolutetheta_c[i] = theta(y=rel_s_c[i,1], z=rel_s_c[i,2]) + thetaprime(x=rel_s_c[i,0], y=rel_s_c[i,1], z=rel_s_c[i,2], t=time_points[i])
            thetaprime_a[i] = thetaprime(x=rel_s_a[i,0], y=rel_s_a[i,1], z=rel_s_a[i,2], t=time_points[i])

        thetafig = plt.figure()
        thetaax1 = thetafig.add_subplot(111)
        thetaax1.set_xlabel('time (T)', fontsize='15')
        thetaax1.set_ylabel('Potential temperature (K)', fontsize='15')
        thetaax1.plot(t[:], absolutetheta_a[:], label='a')
        thetaax1.plot(t[:], absolutetheta_b[:], label='b')
        thetaax1.plot(t[:], absolutetheta_c[:], label='c')
        thetaax1.axhline(y=absolutetheta_a[0], color='black',ls='dotted')
        thetaax1.axhline(y=absolutetheta_b[0], color='black',ls='dotted')
        thetaax1.axhline(y=absolutetheta_c[0], color='black',ls='dotted')
        thetaax1.legend(loc='upper left')

        thetaprimefig = plt.figure()
        thetaax2 = thetaprimefig.add_subplot(111)
        thetaax2.set_xlabel('time (T)', fontsize='15')
        thetaax2.set_ylabel('Theta perturbation (K)', fontsize='15')
        thetaax2.plot(t[:], thetaprime_a[:], color='red')
        thetaax2.axhline(y=thetaprime_a[0], color='black',ls='dotted')

        if showvel==True:
            #Define vertical and meridional heat fluxes
            wheatflux = np.multiply(wprime_a, thetaprime_a)
            vheatflux = np.multiply(vprime_a, thetaprime_a)

            heatfluxfig = plt.figure()
            fluxax1 = heatfluxfig.add_subplot(211)
            fluxax1.set_ylabel('Vertical heat flux (K ms^-1)', fontsize='15')
            fluxax1.plot(t[:], wheatflux[:], color='blue')
            fluxax1.axhline(y=wheatflux[0], color='black', ls='dotted')

            fluxax2 = heatfluxfig.add_subplot(212)
            fluxax2.set_xlabel('time (T)', fontsize='14')
            fluxax2.set_ylabel('Meridional heat flux (K ms^-1)', fontsize='15')
            fluxax2.plot(t[:], vheatflux[:], color='green')
            fluxax2.axhline(y=vheatflux[0], color='black', ls='dotted')

###############################################################################################

    #Show all the Figures

    plt.show()

###############################################################################################

#Run the programme, choose stop to be close to a multiple of efoldtime

main(start=0, stop=31.6, xpos=0, ypos=0.4, zpos=0.95) #T is about a day
