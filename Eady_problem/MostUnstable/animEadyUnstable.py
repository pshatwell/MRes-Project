import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint

import matplotlib.animation as anim
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

#Choose file to import to give trajectories in Earth frame or Wave frame
from EadyUnstableinfo_Earthframe import *
#from EadyUnstableinfo_Waveframe import *

###############################################################################################

#xslice defines meridional plane to view trajectories in
#parameter p defines position of initial cluster

def main(p, xslice):
    #Create coords for plotting

    xvalues = np.linspace(0,6,100)
    zvalues = np.linspace(-1,1,100)

    time = 0

    xlength = len(xvalues)
    zlength = len(zvalues)

###############################################################################################
    #Create arrays for quantities

    #Create empty matrix for streamfunction
    phiprime_matrix = np.zeros((zlength,xlength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            phiprime_matrix[i,j] = phiprime(x=xvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for meridional velocity perturbation
    vprime_matrix = np.zeros((zlength,xlength))

    #Evaluate vprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            vprime_matrix[i,j] = vprime(x=xvalues[j], z=zvalues[i], t=time)


    #Create empty matrix for vertical velocity perturbation
    wprime_matrix = np.zeros((zlength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            wprime_matrix[i,j] = wprime(x=xvalues[j], z=zvalues[i], t=time)

    print 'wprime at z=-H is:', wprime_matrix[0,10]
    print 'wprime at z=H is:', wprime_matrix[-1,10]


    #Create empty matrix for potential temperature perturbation
    thetaprime_matrix = np.zeros((zlength,xlength))

    #Evaluate thetaprime matrix values:
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            thetaprime_matrix[i,j] = thetaprime(x=xvalues[j], z=zvalues[i], t=time)

###############################################################################################

    #Define ODE for parcel trajectories
    def velocity(s,t):
        y,z = s
        dsdt = [vprime(x=xslice,z=z,t=t), wprime(x=xslice,z=z,t=t)]
        return dsdt

    t = np.linspace(0,15,1000) #Define timesteps to integrate over

    #Define line of intial positions for parcels
    s0_a = [p, p/5.]
    s0_b = [p, (p-0.2)/5.]
    s0_c = [p, (p+0.2)/5.]
    s0_d = [p, (p-0.4)/5.]
    s0_e = [p, (p+0.4)/5.]
    s0_f = [p, (p-0.6)/5.]
    s0_g = [p, (p+0.6)/5.]
    s0_h = [p, (p-0.8)/5.]
    s0_i = [p, (p+0.8)/5.]

    #Find solutions for trajectories
    sol_a = odeint(velocity, s0_a, t)
    sol_b = odeint(velocity, s0_b, t)
    sol_c = odeint(velocity, s0_c, t)
    sol_d = odeint(velocity, s0_d, t)
    sol_e = odeint(velocity, s0_e, t)
    sol_f = odeint(velocity, s0_f, t)
    sol_g = odeint(velocity, s0_g, t)
    sol_h = odeint(velocity, s0_h, t)
    sol_i = odeint(velocity, s0_i, t)

###############################################################################################

    #Plot the initial conditions

    fig1 = plt.figure()
    plt.set_cmap('inferno')
    ax1 = fig1.add_subplot(221)
    ax1.set_title('Streamfunction')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('z (H)')
    phicontour = ax1.contourf(phiprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(phicontour)

    ax2 = fig1.add_subplot(222)
    ax2.set_title('v')
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    vcontour = ax2.contourf(vprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(vcontour)

    ax3 = fig1.add_subplot(223)
    ax3.set_title('w')
    ax3.set_xlabel('x (L)')
    ax3.set_ylabel('z (H)')
    wcontour = ax3.contourf(wprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(wcontour)

    ax4 = fig1.add_subplot(224)
    ax4.set_title('theta')
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('z (H)')
    thetacontour = ax4.contourf(thetaprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(thetacontour)

    #plt.show()

###############################################################################################

    #Define figure to plot onto for animation
    fig2 = plt.figure()
    plt.set_cmap('inferno')
    ax5 = fig2.add_subplot(111)
    #ax5.contourf(thetaprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')

    ax5.set_xlabel('y (L)')
    ax5.set_ylabel('z (H)')
    #ax5.set_ylim([-1,1])
    #ax5.set_xlim([-50,50])
    #ax5.grid(axis='y')

    #Plot empty lines on axes to later add data to
    line_a, = ax5.plot([], [])
    line_b, = ax5.plot([], [])
    line_c, = ax5.plot([], [])
    line_d, = ax5.plot([], [])
    line_e, = ax5.plot([], [])
    line_f, = ax5.plot([], [])
    line_g, = ax5.plot([], [])
    line_h, = ax5.plot([], [])
    line_i, = ax5.plot([], [])

    #Define lists to store parcel positions
    xa=[]
    ya=[]

    xb=[]
    yb=[]

    xc=[]
    yc=[]

    xd=[]
    yd=[]

    xe=[]
    ye=[]

    xf=[]
    yf=[]

    xg=[]
    yg=[]

    xh=[]
    yh=[]

    xi=[]
    yi=[]

###############################################################################################

    #Define initialisation function for animation
    def init():
        line_a.set_data([],[])
        line_b.set_data([],[])
        line_c.set_data([],[])
        line_d.set_data([],[])
        line_e.set_data([],[])
        line_f.set_data([],[])
        line_g.set_data([],[])
        line_h.set_data([],[])
        line_i.set_data([],[])
        return line_a, line_b, line_c, line_d, line_e, line_f, line_g, line_h, line_i,

    #Define animation function
    def Eadyanimate(i):
        xa.append(sol_a[i,0])
        ya.append(sol_a[i,1])

        xb.append(sol_b[i,0])
        yb.append(sol_b[i,1])

        xc.append(sol_c[i,0])
        yc.append(sol_c[i,1])

        xd.append(sol_d[i,0])
        yd.append(sol_d[i,1])

        xe.append(sol_e[i,0])
        ye.append(sol_e[i,1])

        xf.append(sol_f[i,0])
        yf.append(sol_f[i,1])

        xg.append(sol_g[i,0])
        yg.append(sol_g[i,1])

        xh.append(sol_h[i,0])
        yh.append(sol_h[i,1])

        xi.append(sol_i[i,0])
        yi.append(sol_i[i,1])

        line_a.set_data(xa,ya)
        line_b.set_data(xb,yb)
        line_c.set_data(xc,yc)
        line_d.set_data(xd,yd)
        line_e.set_data(xe,ye)
        line_f.set_data(xf,yf)
        line_g.set_data(xg,yg)
        line_h.set_data(xh,yh)
        line_i.set_data(xi,yi)

        ax5.autoscale(tight=True)
        ax5.relim()

        return line_a, line_b, line_c, line_d, line_e, line_f, line_g, line_h, line_i,

###############################################################################################

    #Simpler animations, plotting only y or z with time
    #This is only for solution a from above

    fig3 = plt.figure()
    ax6 = fig3.add_subplot(111)
    ax6.set_xlabel('time (T)')
    ax6.set_ylabel('latitude (L)')
    ax6.set_xlim([0,15])
    #ax6.set_ylim([-100,100])
    #ax6.grid()

    fig4 = plt.figure()
    ax7 = fig4.add_subplot(111)
    ax7.set_xlabel('time (T)')
    ax7.set_ylabel('height (H)')
    ax7.set_xlim([0,15])
    #ax7.set_ylim([-1,1])
    #ax7.grid()

    yline, = ax6.plot([], [])
    zline, = ax7.plot([], [])

    ytimesteps = []
    ztimesteps = []

    latitude = []
    height = []

###############################################################################################

    #Define initialisation functions
    def init_y():
        yline.set_data([],[])
        return yline,

    def init_z():
        zline.set_data([],[])
        return zline,

    #Animate the latitude with time
    def y_Eadyanimate(i):
        ytimesteps.append(t[i])
        latitude.append(sol_a[i,0])

        yline.set_data(ytimesteps,latitude)

        ax6.autoscale(axis = 'y', tight=True)
        ax6.relim()

        return yline,

    #Animate the height with time
    def z_Eadyanimate(i):
        ztimesteps.append(t[i])
        height.append(sol_a[i,1])

        zline.set_data(ztimesteps,height)

        ax7.autoscale(axis = 'y', tight=True)
        ax7.relim()

        return zline,

###############################################################################################

    #Make and save animations as mp4 files

    Eadyanim = anim.FuncAnimation(fig2, Eadyanimate, init_func=init, frames=np.arange(0,len(t)-1), interval = 100, blit=True)

    Eadyanim.save('EarthFrameMovies/yztrajectory_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=40, bitrate=-1, codec='libx264', writer='ffmpeg')


    y_Eadyanim = anim.FuncAnimation(fig3, y_Eadyanimate, init_func=init_y, frames=np.arange(0,len(t)-1), interval=100, blit=True)

    y_Eadyanim.save('EarthFrameMovies/latitude_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=40, bitrate=-1, codec='libx264', writer='ffmpeg')


    z_Eadyanim = anim.FuncAnimation(fig4, z_Eadyanimate, init_func=init_z, frames=np.arange(0,len(t)-1), interval=100, blit=True)

    z_Eadyanim.save('EarthFrameMovies/height_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=40, bitrate=-1, codec='libx264', writer='ffmpeg')

###############################################################################################

#Run the program

#xslices = [i for i in np.arange(1,4,0.1)]

#for i in xslices:
#    main(p=0, xslice=i)

main(p=1, xslice=0)
