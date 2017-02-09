import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from Parameters import *
from scipy.integrate import odeint

xvalues = np.floor(np.arange(0,1e3*H_R, 5e5))
zvalues = np.floor(np.arange(-H, H, 5e2))

Xlength = len(xvalues)
height = len(zvalues)

x,z,t = symbols('x z t')

#Define streamfunction perturbation for most unstable solution, Gill 13.3.15
phiprime = (cos(k*x)*(sinh(z/H_R)/sinh(Hratio)) + sin(k*x)*(cosh(z/H_R)/cosh(Hratio)))*exp(sigma_max*t)

#use Sympy to find expressions for the partial derivatives
phiprime_x = diff(phiprime, x)
phiprime_z = diff(phiprime, z)
phiprime_t = diff(phiprime, t)

phiprime_zx = diff(phiprime_x, z)
phiprime_zt = diff(phiprime_t, z)

vprime = (1./f0)*phiprime_x #meridional perturbation velocity, Gill 12.9.3
wprime = (1./N**2)*(shear*phiprime_x - phiprime_zt + shear*z*phiprime_zx) #vertical perturbation velocity, Gill 12.9.6 BOUNDARY CONDITION? NEED TO SOLVE 12.9.5?

vprimef = vprime.subs(x, Xlength/2.).evalf()
wprimef = wprime.subs(x, Xlength/2.).evalf()

########################### FIND PARCEL TRAJECTORIES ###################################################################

#STILL UNSURE WHAT I'M ACTUALLY PLOTTING HERE.

def velocity(s,t): #vprime and wprime expressions found from Sympy, evaluated at Xlength/2
    y, z = s
    dsdt = [10000.0*(-8.50918128239322e-7*np.sin(1.0e-6*Xlength/2.)*np.sinh(0.0001*z) + 6.48054273663885e-7*np.cos(1.0e-6*Xlength/2.)*np.cosh(0.0001*z))*np.exp(1.549e-6*t), 5.0*z*(-8.50918128239322e-11*np.sin(1.0e-6*Xlength/2.)*np.cosh(0.0001*z) + 6.48054273663885e-11*np.cos(1.0e-6*Xlength/2.)*np.sinh(0.0001*z))*np.exp(1.549e-6*t) + 5.0*(-8.50918128239322e-7*np.sin(1.0e-6*Xlength/2.)*np.sinh(0.0001*z) + 6.48054273663885e-7*np.cos(1.0e-6*Xlength/2.)*np.cosh(0.0001*z))*np.exp(1.549e-6*t) - 0.01549*(6.48054273663885e-5*np.sin(1.0e-6*Xlength/2.)*np.sinh(0.0001*z) + 8.50918128239322e-5*np.cos(1.0e-6*Xlength/2.)*np.cosh(0.0001*z))*np.exp(1.549e-6*t)]
    #dsdt = [vprimef, wprimef]
    return dsdt

s0 = [0, 0] #initial position in y-z plane

t = np.linspace(0,1e5,101) #create array of timesteps

sol = odeint(velocity, s0, t) #solve for y and z positions

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(t, sol[:,0], label = 'y position')
#ax.plot(t,sol[:,1], label = 'z position')
ax.plot(sol[:,0],sol[:,1])
#ax.matshow(sol, origin='lower', aspect='auto')
plt.legend(loc='best')

plt.show()
