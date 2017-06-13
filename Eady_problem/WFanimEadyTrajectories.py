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

#Transform to wave frame by shifting x-coordinates
shift = np.zeros_like(sol_a)
shift[:,0] = -c.real*t*(T/L) #Factor of (T/L) to make nondimensional

#Define new x-shifted trajectories
rel_sol_a = sol_a + shift
rel_sol_b = sol_b + shift
rel_sol_c = sol_c + shift
rel_sol_d = sol_d + shift
rel_sol_e = sol_e + shift

rel_sol_f = sol_f + shift
rel_sol_g = sol_g + shift
rel_sol_h = sol_h + shift
rel_sol_i = sol_i + shift
rel_sol_j = sol_j + shift

rel_solutions = [rel_sol_a, rel_sol_b, rel_sol_c, rel_sol_d, rel_sol_e, rel_sol_f, rel_sol_g, rel_sol_h, rel_sol_i, rel_sol_j]

###############################################################################################

#Define normal to theta surfaces

normal = np.array((0, dthetady*L, dthetadz*H))
normalhat = normal/(la.norm(normal))

print 'normal to theta surfaces is:', normalhat

sol_a_p = np.zeros_like(rel_sol_a)
sol_b_p = np.zeros_like(rel_sol_a)
sol_c_p = np.zeros_like(rel_sol_a)
sol_d_p = np.zeros_like(rel_sol_a)
sol_e_p = np.zeros_like(rel_sol_a)
sol_f_p = np.zeros_like(rel_sol_a)
sol_g_p = np.zeros_like(rel_sol_a)
sol_h_p = np.zeros_like(rel_sol_a)
sol_i_p = np.zeros_like(rel_sol_a)
sol_j_p = np.zeros_like(rel_sol_a)

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

projected_solutions = [sol_a_p, sol_b_p, sol_c_p, sol_d_p, sol_e_p, sol_f_p, sol_g_p, sol_h_p, sol_i_p, sol_j_p]

###############################################################################################

#Define figure to plot onto for animation

fig1 = plt.figure()
plt.set_cmap('inferno')
fig1.suptitle('Wave frame')
ax1 = p3.Axes3D(fig1)

ax1.set_xlim3d([0, 6])
ax1.set_ylim3d([-1, 1])
ax1.set_zlim3d([-0.25, 0.25])

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
        positions[m][0].append(projected_solutions[m][i,0])
        positions[m][1].append(projected_solutions[m][i,1])
        positions[m][2].append(projected_solutions[m][i,2])

        lines[m].set_data(positions[m][0], positions[m][1])
        lines[m].set_3d_properties(positions[m][2])

        #ax1.autoscale(tight=True, axis='x')
        #ax1.relim()

    return lines
###############################################################################################

#Make the animation
Eadyanimation = anim.FuncAnimation(fig1, Eadyanimate, init_func=init, frames=np.arange(0,len(t)-1), interval = 100, blit=True)
#Eadyanimation.save('WFtest.gif', fps=30, writer='imagemagick')
Eadyanimation.save('WFtest.mp4', fps=30, bitrate=-1, codec='libx264', writer='ffmpeg')

print 'Made trajectory animation'
