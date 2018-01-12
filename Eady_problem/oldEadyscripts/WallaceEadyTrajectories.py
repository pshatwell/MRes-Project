from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import numpy.linalg as la

from Eadyinfo_oc import *
#from Eadyinfo_atm import *

'''WHY DO IMPLICIT MIDPOINT SOLUTION VALUES SEEM TO BE TWICE AS LARGE
AS THEY SHOULD BE??? E.G. MERIDIONAL CHANNEL EXTENT IS -2 TO + 2, VERTICAL HEIGHT IS +2'''

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

def main(start, stop, zpos=0.5, ypos=0.5):

    def velocity(t,s):
        x,y,z = s
        #dsdt = np.array((uprime(x=x,y=y,z=z,t=t), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)))
        dsdt = np.array((uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)))
        return dsdt

    def velocity2(s,t):
        x,y,z = s
        #dsdt = [uprime(x=x,y=y,z=z,t=t), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        dsdt = [uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
        return dsdt

###############################################################################################

    #Define timesteps for integration

    tmin = start
    tmax = stop
    nt = 500

    t = np.linspace(tmin, tmax, nt)

    dt = (tmax - tmin)/nt

    print 'length of t is', len(t)
    print 'dt is', dt

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

    '''CURRENTLY DIVIDING SOLUTIONS BY TWO TO CORRECT THEM. DON'T KNOW WHY.'''
    #Create list of Earth frame trajectories
    #EFsolutions = np.array((sol_a, sol_b, sol_c))/2.
    EFsolutions = np.array((sol_a, sol_b, sol_c))

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

    '''CURRENTLY DIVIDING SOLUTIONS BY TWO TO CORRECT THEM. DON'T KNOW WHY.'''
    #Create list of Wave frame trajectories
    #WFsolutions = np.array((rel_sol_a, rel_sol_b, rel_sol_c))/2.
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

    #Find extents for arrays to display graphically behind trajectories
    xmins = []
    for trajectory in WFsolutions:
        xmins.append(min(trajectory[:-1,0]))

    xmaxes = []
    for trajectory in WFsolutions:
        xmaxes.append(max(trajectory[:-1,0]))

    ymins = []
    for trajectory in WFsolutions:
        ymins.append(min(trajectory[:-1,1]))

    ymaxes = []
    for trajectory in WFsolutions:
        ymaxes.append(max(trajectory[:-1,1]))

    zmins = []
    for trajectory in WFsolutions:
        zmins.append(min(trajectory[:-1,2]))

    zmaxes = []
    for trajectory in WFsolutions:
        zmaxes.append(max(trajectory[:-1,2]))

    xmin = min(xmins)
    xmax = max(xmaxes)

    ymin = min(ymins)
    ymax = max(ymaxes)

    zmin = min(zmins)
    zmax = max(zmaxes)

    xvalues = np.linspace(xmin, xmax, 50)
    yvalues = np.linspace(ymin, ymax, 50)
    zvalues = np.linspace(zmin, zmax, 50)

    xlength = len(xvalues)
    ylength = len(yvalues)
    zlength = len(zvalues)

###############################################################################################

    #Define arrays for Eady solution
    #Solutions are functions in Eadyinfo script

    time = stop #time at which to evaluate arrays for Eady solution

    theta_matrix = np.zeros((zlength, ylength))

    for i in range(0, 50, 1):
        for j in range(0, 50, 1):
            theta_matrix[i,j] = theta(y=yvalues[j], z=zvalues[i])

    print 'dthetady is:', dthetady*L
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

    #Create empty matrix for uprime
    uprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate uprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                uprime_matrix[i,j,m] = uprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    #Create empty matrix for vprime
    vprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate vprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                vprime_matrix[i,j,m] = vprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

    #Create empty matrix for wprime
    wprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                wprime_matrix[i,j,m] = wprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)
                #wprime_matrix[i,j,m] = wprime(x=2*xvalues[m], y = 2*yvalues[j], z=2*zvalues[i], t=time)

    #Create empty matrix for streamfunction perturbation
    psiprime_matrix = np.zeros((zlength,ylength,xlength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, ylength, 1):
            for m in range(0, xlength, 1):
                psiprime_matrix[i,j,m] = psiprime(x=xvalues[m], y = yvalues[j], z=zvalues[i], t=time)

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

    sol_a_pi = np.zeros_like(rel_sol_a) #i for instability
    sol_b_pi = np.zeros_like(rel_sol_a)
    sol_c_pi = np.zeros_like(rel_sol_a)

    #Project WF solutions onto theta surface
    for i in range(len(t)):
        sol_a_p[i] = rel_sol_a[i] - np.dot(rel_sol_a[i], normalhat)*normalhat
        sol_b_p[i] = rel_sol_b[i] - np.dot(rel_sol_b[i], normalhat)*normalhat
        sol_c_p[i] = rel_sol_c[i] - np.dot(rel_sol_c[i], normalhat)*normalhat

    #Project WF solutions onto sheet at half-slope of theta surface
    for i in range(len(t)):
        sol_a_pi[i] = rel_sol_a[i] - np.dot(rel_sol_a[i], normal_ihat)*normal_ihat
        sol_b_pi[i] = rel_sol_b[i] - np.dot(rel_sol_b[i], normal_ihat)*normal_ihat
        sol_c_pi[i] = rel_sol_c[i] - np.dot(rel_sol_c[i], normal_ihat)*normal_ihat

    projected_solutions = [sol_a_p, sol_b_p, sol_c_p]
    projected_solutions_i = [sol_a_pi, sol_b_pi, sol_c_pi]

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


###############################################################################################

    #Plot of velocity perturbations with time

    uprime_a = np.zeros_like(t)
    uprime_b = np.zeros_like(t)
    uprime_c = np.zeros_like(t)

    vprime_a = np.zeros_like(t)
    vprime_b = np.zeros_like(t)
    vprime_c = np.zeros_like(t)

    wprime_a = np.zeros_like(t)
    wprime_b = np.zeros_like(t)
    wprime_c = np.zeros_like(t)

    for i in range(len(t)):
        uprime_a[i] = uprime(x=WFsolutions[0,i,0], y=WFsolutions[0,i,1], z=WFsolutions[0,i,2], t=t[i])
        uprime_b[i] = uprime(x=WFsolutions[1,i,0], y=WFsolutions[1,i,1], z=WFsolutions[1,i,2], t=t[i])
        uprime_c[i] = uprime(x=WFsolutions[2,i,0], y=WFsolutions[2,i,1], z=WFsolutions[2,i,2], t=t[i])

        vprime_a[i] = vprime(x=WFsolutions[0,i,0], y=WFsolutions[0,i,1], z=WFsolutions[0,i,2], t=t[i])
        vprime_b[i] = vprime(x=WFsolutions[1,i,0], y=WFsolutions[1,i,1], z=WFsolutions[1,i,2], t=t[i])
        vprime_c[i] = vprime(x=WFsolutions[2,i,0], y=WFsolutions[2,i,1], z=WFsolutions[2,i,2], t=t[i])

        wprime_a[i] = wprime(x=WFsolutions[0,i,0], y=WFsolutions[0,i,1], z=WFsolutions[0,i,2], t=t[i])
        wprime_b[i] = wprime(x=WFsolutions[1,i,0], y=WFsolutions[1,i,1], z=WFsolutions[1,i,2], t=t[i])
        wprime_c[i] = wprime(x=WFsolutions[2,i,0], y=WFsolutions[2,i,1], z=WFsolutions[2,i,2], t=t[i])

    velocityfig = plt.figure()
    uprimeax = velocityfig.add_subplot(311)
    uprimeax.set_xlabel('time (T)')
    uprimeax.set_ylabel('u prime')
    uprimeax.plot(t[:-1], uprime_a[:-1])
    uprimeax.plot(t[:-1], uprime_b[:-1])
    uprimeax.plot(t[:-1], uprime_c[:-1])
    uprimeax.axhline(y = uprime_a[0], color='black', ls='dotted')

    vprimeax = velocityfig.add_subplot(312)
    vprimeax.set_xlabel('time (T)')
    vprimeax.set_ylabel('v prime')
    vprimeax.plot(t[:-1], vprime_a[:-1])
    vprimeax.plot(t[:-1], vprime_b[:-1])
    vprimeax.plot(t[:-1], vprime_c[:-1])
    vprimeax.axhline(y = vprime_a[0], color='black', ls='dotted')

    wprimeax = velocityfig.add_subplot(313)
    wprimeax.set_xlabel('time (T)')
    wprimeax.set_ylabel('w prime')
    wprimeax.plot(t[:-1], wprime_a[:-1])
    wprimeax.plot(t[:-1], wprime_b[:-1])
    wprimeax.plot(t[:-1], wprime_c[:-1])
    wprimeax.axhline(y = wprime_a[0], color='black', ls='dotted')

###############################################################################################

    #Wallace 78 Figures

    '''These only really make sense ignoring the mean flow in the velocity function'''
    '''We only see the Wallace parcel motion below the steering level here, I think'''

    if sigma != 0:
        efoldtime = 1./(sigma*T)
        print 'efoldtime is', efoldtime

        times = np.arange(start, len(t), int((efoldtime*len(t)/stop)))


    zslice = int(round(zlength*zpos))
    yslice = int(round(ylength/2.))


    #Attempt for 'Figure 3' from Wallace 78 paper

    fig7 = plt.figure()
    ax15 = fig7.add_subplot(111)
    ax15.set_xlabel('Latitude (L)')
    ax15.set_ylabel('Height (H)')
    isentropes = ax15.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    plt.colorbar(isentropes)
    for i in WFsolutions:
        ax15.plot(i[:-1,1], i[:-1,2], lw='1.5', color='black')
    #for i in WFsolutions2:
    #    ax15.plot(i[:,1], i[:,2], lw='1.5', color='green')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            ax15.scatter(j[times[i],1], j[times[i],2], marker='o', c='black', s=9)


    #Attempt at 'Figures 1 and 2' from Wallace 78 paper

    fig8 = plt.figure()

    ax16 = fig8.add_subplot(211)
    #ax16.set_title('On top of theta perturbation')
    #ax16.set_xlabel('Longitude (L)')
    ax16.set_ylabel('Latitude (L)')
    plt.set_cmap('Reds')
    xy_contour_theta = ax16.contourf(thetaprime_matrix[zslice,:,:], origin='lower', aspect='auto', extent = [xmin,xmax,ymin,ymax])
    #plt.colorbar(xy_contour_theta)
    xy_contour_psi1 = ax16.contour(psiprime_matrix[zslice,:,:], colors='white', origin='lower', aspect='auto', extent=[xmin,xmax,ymin,ymax])
    #plt.clabel(xy_contour_psi1, inline=1, fontsize=10)
    for i in WFsolutions:
        ax16.plot(i[:-1,0], i[:-1,1], lw='1.5', color='black')
    #for i in WFsolutions2:
    #    ax16.plot(i[:,0], i[:,1], lw='1.5', color='green')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            ax16.scatter(j[times[i],0], j[times[i],1], marker='o', c='black', s=9)


    ax17 = fig8.add_subplot(212)
    #ax17.set_title('On top of w perturbation')
    ax17.set_xlabel('Longitude (L)')
    ax17.set_ylabel('Latitude (L)')
    plt.set_cmap('Blues')
    xy_contour_w = ax17.contourf(wprime_matrix[zslice,:,:], aspect='auto', extent = [xmin,xmax,ymin,ymax])
    #plt.colorbar(xy_contour_w)
    xy_contour_psi2 = ax17.contour(psiprime_matrix[zslice,:,:], colors='white', origin='lower', aspect='auto', extent=[xmin,xmax,ymin,ymax])
    #plt.clabel(xy_contour_psi2, inline=1, fontsize=10)
    for i in WFsolutions:
        ax17.plot(i[:-1,0], i[:-1,1], lw='1.5', color='black')
    #for i in WFsolutions2:
    #    ax17.plot(i[:,0], i[:,1], lw='1.5', color='green')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            ax17.scatter(j[times[i],0], j[times[i],1], marker='o', c='black', s=9)


    #'Figures 1 and 2' in the x-z plane

    fig9 = plt.figure()

    ax18 = fig9.add_subplot(211)
    ax18.set_ylabel('Height (H)')
    plt.set_cmap('Reds')
    xz_contour_theta = ax18.contourf(thetaprime_matrix[:,yslice,:], origin='lower', aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_theta)
    for i in WFsolutions:
        ax18.plot(i[:-1,0], i[:-1,2], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            ax18.scatter(j[times[i],0], j[times[i],2], marker='o', c='black', s=9)


    ax19 = fig9.add_subplot(212)
    ax19.set_xlabel('Longitude (L)')
    ax19.set_ylabel('Height (H)')
    plt.set_cmap('Blues')
    xz_contour_w = ax19.contourf(wprime_matrix[:,yslice,:], aspect='auto', extent = [xmin,xmax,zmin,zmax])
    plt.colorbar(xz_contour_w)
    for i in WFsolutions:
        ax19.plot(i[:-1,0], i[:-1,2], lw='1.5', color='black')

    #Add a dot along trajectories every e-folding time to indicate time evolution
    for i in range(len(times)):
        for j in WFsolutions:
            ax19.scatter(j[times[i],0], j[times[i],2], marker='o', c='black', s=9)


    plt.show()

###############################################################################################

#Run the programme

'''NOTE, THE STEERING LEVEL IS SOMEHOW Z=1 FOR THE IMPLICIT MIDPOINT SOLUTIONS?'''

main(start=0, stop=60, zpos=0.8, ypos=0.5)
