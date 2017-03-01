import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from EadyUnstableinfo import *
import matplotlib.animation as anim
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

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

    t = np.linspace(0,15,100) #Define timesteps to integrate over

    #Define cluster of intial positions for parcels
    s0_a = [p, p/5.]
    s0_b = [p-0.1, (p-0.1)/5.]
    s0_c = [p+0.1, (p+0.1)/5.]
    s0_d = [p-0.2, (p-0.2)/5.]
    s0_e = [p+0.2, (p+0.2)/5.]
    s0_f = [p-0.3, (p-0.3)/5.]
    s0_g = [p+0.3, (p+0.3)/5.]
    s0_h = [p-0.4, (p-0.4)/5.]
    s0_i = [p+0.4, (p+0.4)/5.]

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

    ax5.set_xlabel('y (L)')
    ax5.set_ylabel('z (H)')

    ax5.autoscale_view(True,True,True)

    #Plot empty lines on axes to later add data to
    line_a, = ax5.plot([], [], lw=1.5)
    line_b, = ax5.plot([], [], lw=1.5)
    line_c, = ax5.plot([], [], lw=1.5)
    line_d, = ax5.plot([], [], lw=1.5)
    line_e, = ax5.plot([], [], lw=1.5)
    line_f, = ax5.plot([], [], lw=1.5)
    line_g, = ax5.plot([], [], lw=1.5)
    line_h, = ax5.plot([], [], lw=1.5)
    line_i, = ax5.plot([], [], lw=1.5)

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

        ax5.autoscale()
        ax5.relim()

        return line_a, line_b, line_c, line_d, line_e, line_f, line_g, line_h, line_i,

###############################################################################################

    #Make and save animation as mp4 file

    Eadyanim = anim.FuncAnimation(fig2, Eadyanimate, init_func=init, frames=np.arange(0,len(t)-1), interval = 100, blit=False)

    Eadyanim.save('movies/p=4/Eadybundle_p_%s_xslice_%s.mp4' % (str(p), str(xslice)), fps=20, bitrate=-1, codec='libx264', writer='ffmpeg')

xslices = [i for i in np.arange(0,6,0.1)]

for i in xslices:
    main(p=4, xslice=i)
