import numpy as np
from sympy import *
from RWinfo import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as anim
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

###############################################################################################

#CHOOSE THE POINT P TO DETERMINE CLUSTER OF INITIAL POINTS FOR TRAJECTORIES
p=0.9

#Define dimensions of arrays
xvalues = np.linspace(0,2,200)
yvalues = np.linspace(0,1,100) #Channel width in y direction is 1*L

time = 0 #Choose time step (in units of T)

xlength = len(xvalues)
ylength = len(yvalues)

###############################################################################################

#CREATE ARRAY FOR PSIPRIME VALUES
psiprime_matrix = np.zeros((ylength,xlength)) #Note x and y have switched so plotted output is oriented correctly

for i in range(0, ylength, 1):
    for j in range(0, xlength, 1):
        psiprime_matrix[i,j] = psiprime(x=xvalues[j], y=yvalues[i], t=time)

print 'Array for psiprime made.'

###############################################################################################

#DEFINE ODE TO INTEGRATE FOR PARCEL TRAJECTORIES
def velocity(s,t):
    x, y = s
    dsdt = [uprime(x=x,y=y,t=t), vprime(x=x,y=y,t=t)]
    return dsdt

t=np.linspace(0,2000,250)

#Create 3x3 grid of 9 points to evolve in time
s0_a = [p, p/2.]
s0_b = [p-0.1, p/2.]
s0_c = [p+0.1, p/2.]
s0_d = [p, (p/2.)+0.1]
s0_e = [p-0.1, (p/2.)+0.1]
s0_f = [p+0.1, (p/2.)+0.1]
s0_g = [p, (p/2.)-0.1]
s0_h = [p-0.1, (p/2.)-0.1]
s0_i = [p+0.1, (p/2.)-0.1]

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

#TRANSFORM TRAJECTORIES TO WAVE REFERENCE FRAME
shift = np.zeros_like(sol_a)
shift[:,0] = (c*t)/(L*f0) #to make non-dimensional

rel_sol_a = sol_a - shift #new parcel trajectory relative to Rossby wave
rel_sol_b = sol_b - shift
rel_sol_c = sol_c - shift
rel_sol_d = sol_d - shift
rel_sol_e = sol_e - shift
rel_sol_f = sol_f - shift
rel_sol_g = sol_g - shift
rel_sol_h = sol_h - shift
rel_sol_i = sol_i - shift

print 'Trajectory solutions found.'

###############################################################################################

#Animate the trajectories

fig = plt.figure()
plt.set_cmap('inferno')
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_xlim([0,2])
ax1.set_ylim([0,1])
ax1.set_xlabel('x (L)')
ax1.set_ylabel('y (L)')
#ax1.imshow(psiprime_matrix[:,:], origin='lower', extent=[0,2,0,1], aspect='auto')

ax2.set_xlim([0,2])
ax2.set_ylim([0,1])
ax2.set_xlabel('x (L)')
ax2.set_ylabel('y (L)')
ax2.imshow(psiprime_matrix[:,:], origin='lower', extent=[0,2,0,1], aspect='auto')

line_a1, = ax1.plot([], [], lw=1.5)
line_b1, = ax1.plot([], [], lw=1.5)
line_c1, = ax1.plot([], [], lw=1.5)
line_d1, = ax1.plot([], [], lw=1.5)
line_e1, = ax1.plot([], [], lw=1.5)
line_f1, = ax1.plot([], [], lw=1.5)
line_g1, = ax1.plot([], [], lw=1.5)
line_h1, = ax1.plot([], [], lw=1.5)
line_i1, = ax1.plot([], [], lw=1.5)

line_a2, = ax2.plot([], [], lw=1.5)
line_b2, = ax2.plot([], [], lw=1.5)
line_c2, = ax2.plot([], [], lw=1.5)
line_d2, = ax2.plot([], [], lw=1.5)
line_e2, = ax2.plot([], [], lw=1.5)
line_f2, = ax2.plot([], [], lw=1.5)
line_g2, = ax2.plot([], [], lw=1.5)
line_h2, = ax2.plot([], [], lw=1.5)
line_i2, = ax2.plot([], [], lw=1.5)

#LOOK INTO BETTER STRUCTURING ALL OF THESE LISTS
xa1=[]
ya1=[]
xa2=[]
ya2=[]

xb1=[]
yb1=[]
xb2=[]
yb2=[]

xc1=[]
yc1=[]
xc2=[]
yc2=[]

xd1=[]
yd1=[]
xd2=[]
yd2=[]

xe1=[]
ye1=[]
xe2=[]
ye2=[]

xf1=[]
yf1=[]
xf2=[]
yf2=[]

xg1=[]
yg1=[]
xg2=[]
yg2=[]

xh1=[]
yh1=[]
xh2=[]
yh2=[]

xi1=[]
yi1=[]
xi2=[]
yi2=[]

def init():
    line_a1.set_data([],[])
    line_a2.set_data([],[])

    line_b1.set_data([],[])
    line_b2.set_data([],[])

    line_c1.set_data([],[])
    line_c2.set_data([],[])

    line_d1.set_data([],[])
    line_d2.set_data([],[])

    line_e1.set_data([],[])
    line_e2.set_data([],[])

    line_f1.set_data([],[])
    line_f2.set_data([],[])

    line_g1.set_data([],[])
    line_g2.set_data([],[])

    line_h1.set_data([],[])
    line_h2.set_data([],[])

    line_i1.set_data([],[])
    line_i2.set_data([],[])
    return line_a1, line_a2, line_b1, line_b2, line_c1, line_c2, line_d1, line_d2, line_e1, line_e2, line_f1, line_f2, line_g1, line_g2, line_h1, line_h2, line_i1, line_i2,

def parcelanimate(i):
    xa1.append(sol_a[i,0])
    ya1.append(sol_a[i,1])
    xa2.append(rel_sol_a[i,0])
    ya2.append(rel_sol_a[i,1])

    xb1.append(sol_b[i,0])
    yb1.append(sol_b[i,1])
    xb2.append(rel_sol_b[i,0])
    yb2.append(rel_sol_b[i,1])

    xc1.append(sol_c[i,0])
    yc1.append(sol_c[i,1])
    xc2.append(rel_sol_c[i,0])
    yc2.append(rel_sol_c[i,1])

    xd1.append(sol_d[i,0])
    yd1.append(sol_d[i,1])
    xd2.append(rel_sol_d[i,0])
    yd2.append(rel_sol_d[i,1])

    xe1.append(sol_e[i,0])
    ye1.append(sol_e[i,1])
    xe2.append(rel_sol_e[i,0])
    ye2.append(rel_sol_e[i,1])

    xf1.append(sol_f[i,0])
    yf1.append(sol_f[i,1])
    xf2.append(rel_sol_f[i,0])
    yf2.append(rel_sol_f[i,1])

    xg1.append(sol_g[i,0])
    yg1.append(sol_g[i,1])
    xg2.append(rel_sol_g[i,0])
    yg2.append(rel_sol_g[i,1])

    xh1.append(sol_h[i,0])
    yh1.append(sol_h[i,1])
    xh2.append(rel_sol_h[i,0])
    yh2.append(rel_sol_h[i,1])

    xi1.append(sol_i[i,0])
    yi1.append(sol_i[i,1])
    xi2.append(rel_sol_i[i,0])
    yi2.append(rel_sol_i[i,1])

    line_a1.set_data(xa1,ya1)
    line_a2.set_data(xa2,ya2)
    line_b1.set_data(xb1,yb1)
    line_b2.set_data(xb2,yb2)
    line_c1.set_data(xc1,yc1)
    line_c2.set_data(xc2,yc2)
    line_d1.set_data(xd1,yd1)
    line_d2.set_data(xd2,yd2)
    line_e1.set_data(xe1,ye1)
    line_e2.set_data(xe2,ye2)
    line_f1.set_data(xf1,yf1)
    line_f2.set_data(xf2,yf2)
    line_g1.set_data(xg1,yg1)
    line_g2.set_data(xg2,yg2)
    line_h1.set_data(xh1,yh1)
    line_h2.set_data(xh2,yh2)
    line_i1.set_data(xi1,yi1)
    line_i2.set_data(xi2,yi2)

    #line1.axes.axis([0,2,0,1])
    #line2.axes.axis([0,2,0,1])

    return line_a1, line_a2, line_b1, line_b2, line_c1, line_c2, line_d1, line_d2, line_e1, line_e2, line_f1, line_f2, line_g1, line_g2, line_h1, line_h2, line_i1, line_i2,

print 'Start making animation:'

parcelanim = anim.FuncAnimation(fig, parcelanimate, init_func=init, frames=np.arange(0,len(t)-1), interval = 100, blit=True)

parcelanim.save('movies/parcelbundle_pvalue%s_cvalue%s.mp4' % (str(p), str(c)), fps=25, bitrate=-1, codec='libx264', writer='ffmpeg')

print 'Animation has been made.'
