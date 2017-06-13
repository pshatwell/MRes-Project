import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.integrate import odeint
import numpy.linalg as la

import matplotlib.animation as anim
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

from Eadyinfo import *

'''AUTOSCALING AXES CURRENTLY DOES NOT WORK. HAVE TO SET AXES LIMITS.'''

###############################################################################################

#Define velocity function for 3d parcel trajectories

def velocity3d(s,t):
    x,y,z = s
    dsdt = [uprime(x=x,y=y,z=z,t=t) + umeanflow(z=z), vprime(x=x,y=y,z=z,t=t), wprime(x=x,y=y,z=z,t=t)]
    return dsdt


###############################################################################################

#Define timesteps for integration

tmin = 0
tmax = 60

t = np.linspace(tmin, tmax, 500)

###############################################################################################

#Define initial positions of parcels

#setting both 'cold' and 'warm' parcels off at the same height but different y
#seems to give closest reproduction of Green picture
top = 0.5
bottom = 0.5
x1 = 0
xshift = np.pi/(k*L) #nondimensional shift of half a wavelength
y1 = 0.85

#Set of 5 'warm' parcels
s0_a = np.array((x1, y1, top))
s0_b = np.array((x1+0.5, y1, top))
s0_c = np.array((x1+1, y1, top))
s0_d = np.array((x1+1.5, y1, top))
s0_e = np.array((x1+2, y1, top))

#Set of 5 'cold' parcels
s0_f = np.array((x1+xshift, -y1, bottom))
s0_g = np.array((x1+xshift+0.5, -y1, bottom))
s0_h = np.array((x1+xshift+1, -y1, bottom))
s0_i = np.array((x1+xshift+1.5, -y1, bottom))
s0_j = np.array((x1+xshift+2, -y1, bottom))

###############################################################################################

#Solve for parcel trajectories

sol_a = odeint(velocity3d, s0_a, t)
sol_b = odeint(velocity3d, s0_b, t)
sol_c = odeint(velocity3d, s0_c, t)
sol_d = odeint(velocity3d, s0_d, t)
sol_e = odeint(velocity3d, s0_e, t)

sol_f = odeint(velocity3d, s0_f, t)
sol_g = odeint(velocity3d, s0_g, t)
sol_h = odeint(velocity3d, s0_h, t)
sol_i = odeint(velocity3d, s0_i, t)
sol_j = odeint(velocity3d, s0_j, t)

solutions = [sol_a, sol_b, sol_c, sol_d, sol_e, sol_f, sol_g, sol_h, sol_i, sol_j]

###############################################################################################

#Define figure to plot onto for animation

fig1 = plt.figure()
plt.set_cmap('inferno')
fig1.suptitle('Earth frame')
ax1 = p3.Axes3D(fig1)

ax1.set_xlim3d([0,30])
ax1.set_ylim3d([-1,1])
ax1.set_zlim3d([0,1])

ax1.set_xlabel('x (L)')
ax1.set_ylabel('y (L)')
ax1.set_zlabel('z (H)')

###############################################################################################

#Plot empty lines on axes to later add data to

line_a, = ax1.plot([], [], [])
line_b, = ax1.plot([], [], [])
line_c, = ax1.plot([], [], [])
line_d, = ax1.plot([], [], [])
line_e, = ax1.plot([], [], [])
line_f, = ax1.plot([], [], [])
line_g, = ax1.plot([], [], [])
line_h, = ax1.plot([], [], [])
line_i, = ax1.plot([], [], [])
line_j, = ax1.plot([], [], [])

lines = [line_a, line_b, line_c, line_d, line_e, line_f, line_g, line_h, line_i, line_j]

###############################################################################################

#Define lists to store parcel positions

a = [[],[],[]]
b = [[],[],[]]
c = [[],[],[]]
d = [[],[],[]]
e = [[],[],[]]
f = [[],[],[]]
g = [[],[],[]]
h = [[],[],[]]
i = [[],[],[]]
j = [[],[],[]]

positions = [a,b,c,d,e,f,g,h,i,j]

###############################################################################################

#Define initialisation function for animation
def init():
    for i in range(len(lines)):
        lines[i].set_data([],[])
        lines[i].set_3d_properties([])

    return lines

#Define animation function
def Eadyanimate(i):
    for m in range(len(positions)):
        positions[m][0].append(solutions[m][i,0])
        positions[m][1].append(solutions[m][i,1])
        positions[m][2].append(solutions[m][i,2])

        lines[m].set_data(positions[m][0], positions[m][1])
        lines[m].set_3d_properties(positions[m][2])

        #ax1.autoscale(tight=True, axis='x')
        #ax1.relim()

    return lines
###############################################################################################

#Make the animation
Eadyanimation = anim.FuncAnimation(fig1, Eadyanimate, init_func=init, frames=np.arange(0,len(t)-1), interval = 100, blit=True)
Eadyanimation.save('EFtest.mp4', fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')

print 'Made trajectory animation'
