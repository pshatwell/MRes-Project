import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint

import matplotlib.animation as anim
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

from EadyUnstableinfo import *
from EadyUnstableSolution import *

###############################################################################################

#xslice defines meridional plane to view trajectories in
#parameter p defines position of initial cluster

def main(p, xslice):

###############################################################################################

    #Define ODE for parcel trajectories
    def velocity(s,t):
        y,z = s
        dsdt = [vprime(x=xslice,z=z,t=t), wprime(x=xslice,z=z,t=t)]
        return dsdt

    tmin = 0
    tmax = 5

    t = np.linspace(tmin, tmax, 500) #Define timesteps to integrate over

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

    #Define figure to plot onto for animation
    fig1 = plt.figure()
    plt.set_cmap('inferno')
    ax1 = fig1.add_subplot(111)

    ax1.set_xlabel('y (L)')
    ax1.set_ylabel('z (H)')
    #ax1.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')

    #Plot empty lines on axes to later add data to
    line_a, = ax1.plot([], [])
    line_b, = ax1.plot([], [])
    line_c, = ax1.plot([], [])
    line_d, = ax1.plot([], [])
    line_e, = ax1.plot([], [])
    line_f, = ax1.plot([], [])
    line_g, = ax1.plot([], [])
    line_h, = ax1.plot([], [])
    line_i, = ax1.plot([], [])

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

        ax1.autoscale(tight=True)
        ax1.relim()

        return line_a, line_b, line_c, line_d, line_e, line_f, line_g, line_h, line_i,

###############################################################################################

    #Simpler animations, plotting only y or z with time
    #This is only for solution a from above

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_xlabel('time (T)')
    ax2.set_ylabel('latitude (L)')
    ax2.set_xlim([tmin, tmax])
    ax2.grid()

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.set_xlabel('time (T)')
    ax3.set_ylabel('height (H)')
    ax3.set_xlim([tmin, tmax])
    ax3.set_ylim([zmin, zmax])
    ax3.grid()

    yline, = ax2.plot([], [])
    zline, = ax3.plot([], [])

    #Define lists to store positions
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

        ax2.autoscale(axis = 'y', tight=True)
        ax2.relim()

        return yline,

    #Animate the height with time
    def z_Eadyanimate(i):
        ztimesteps.append(t[i])
        height.append(sol_a[i,1])

        zline.set_data(ztimesteps,height)

        return zline,

###############################################################################################

    #Make and save animations as mp4 files
    #MAKE SURE TO CHOOSE CORRECT DIRECTORY DEPENDING ON REFERENCE FRAME

    Eadyanim = anim.FuncAnimation(fig1, Eadyanimate, init_func=init, frames=np.arange(0,len(t)-1), interval = 100, blit=True)
    #Eadyanim.save('movies/EarthFrameMovies/yztrajectory_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')
    Eadyanim.save('movies/WaveFrameMovies/yztrajectory_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')
    print 'Made trajectory animation'

    y_Eadyanim = anim.FuncAnimation(fig2, y_Eadyanimate, init_func=init_y, frames=np.arange(0,len(t)-1), interval=100, blit=True)
    #y_Eadyanim.save('movies/EarthFrameMovies/latitude_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')
    y_Eadyanim.save('movies/WaveFrameMovies/latitude_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')
    print 'Made latitude animation'

    z_Eadyanim = anim.FuncAnimation(fig3, z_Eadyanimate, init_func=init_z, frames=np.arange(0,len(t)-1), interval=100, blit=True)
    #z_Eadyanim.save('movies/EarthFrameMovies/height_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')
    z_Eadyanim.save('movies/WaveFrameMovies/height_c_%s_p_%s_xslice_%s.mp4' % (str(c), str(p), str(xslice)), fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')
    print 'Made height animation'

###############################################################################################

#Run the program

#xslices = [i for i in np.arange(1,4,0.1)]

#for i in xslices:
#    main(p=0, xslice=i)

main(p=4, xslice=0)
