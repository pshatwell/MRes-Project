import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sympy import *
from scipy.integrate import odeint

from Parameters import *

'''
Program for plotting solution to Eady problem for most unstable wave - Peter Shatwell 2017
'''

############################## SET UP COORDS FOR PLOTTING ################################################################

#xvalues = np.floor(np.arange(0, 1e7, 5e5))
#zvalues = np.floor(np.arange(-5e3, 5e3, 5e2))

xvalues = np.floor(np.arange(0,1e3*H_R, 5e5))
zvalues = np.floor(np.arange(-H, H, 5e2))

Xlength = len(xvalues)
height = len(zvalues)

########################## DEFINE STREAMFUNCTION AND DERIVATIVES ####################################################################

#def phiprime(x,z,t): #Define function for streamfunction perturbation
#    return (np.cos(k*x)*np.divide(np.sinh(z/H_R), np.sinh(Hratio)) + np.sin(k*x)*np.divide(np.cosh(z/H_R), np.cosh(Hratio)))*np.exp(sigma_max*t)

x,z,t = symbols('x z t')

#Define streamfunction perturbation for most unstable solution, Gill 13.3.15
phiprime = (cos(k*x)*(sinh(z/H_R)/sinh(Hratio)) + sin(k*x)*(cosh(z/H_R)/cosh(Hratio)))*exp(sigma_max*t)

#use Sympy to find expressions for the partial derivatives
phiprime_x = diff(phiprime, x)
phiprime_z = diff(phiprime, z)
phiprime_t = diff(phiprime, t)

phiprime_zx = diff(phiprime_x, z)
phiprime_zt = diff(phiprime_t, z)
phiprime_xx = diff(phiprime_x, x)

phiprime_txx = diff(phiprime_xx, t)
phiprime_xxx = diff(phiprime_xx, x)

print 'phiprime_txx is:', phiprime_txx
print 'phiprime_xxx is:', phiprime_xxx

############################## CREATE MATRIX OF PHIPRIME VALUES ################################################################

#phiprime is the streamfunction perturbation

#Initialise matrix of zeros for x-z plane values for phiprime
phiprime_zxmatrix = np.zeros((height,Xlength))

#Evaluate phiprime matrix values
for i in range(0, height, 1):
    for j in range(0, Xlength, 1):
        phiprime_zxmatrix[i-1][j-1] = phiprime.subs(x,xvalues[j-1]).subs(z,zvalues[i-1]).subs(t,0).evalf()

########################### FIND VELOCITY PERTURBATIONS FROM STREAMFUNCTION ###################################################################

#THERE IS NO UPRIME. PHI SOLUTION IS INDEPENDENT OF Y

#CHECK EXPRESSION FOR WPRIME. SHOULD USE GILL 12.9.5? THAT'S A BIT HARD TO SOLVE.

vprime = (1./f0)*phiprime_x #meridional perturbation velocity, Gill 12.9.3
wprime = (1./N**2)*(shear*phiprime_x - phiprime_zt + shear*z*phiprime_zx) #vertical perturbation velocity, Gill 12.9.6 BOUNDARY CONDITION? NEED TO SOLVE 12.9.5?

print 'wprime at z=0 is:', wprime.subs(z,0).subs(t,1).subs(x, Xlength/2.).evalf()

#print 'vprime is:', vprime
#print 'wprime is:', wprime

#Initialise matrix of zeros for vprime values in x-z plane
vprime_matrix = np.zeros((height,Xlength))

#Evaluate values for vprime matrix
for i in range(0, height, 1):
    for j in range(0, Xlength, 1):
        vprime_matrix[i-1][j-1] = vprime.subs(x,xvalues[j-1]).subs(z,zvalues[i-1]).subs(t,0).evalf()

#Initialise matrix of zeros for wprime values in x-z plane
wprime_matrix = np.zeros((height,Xlength))

#Evaluate values for wprime matrix
for i in range(0, height, 1):
    for j in range(0, Xlength, 1):
        wprime_matrix[i-1][j-1] = wprime.subs(x,xvalues[j-1]).subs(z,zvalues[i-1]).subs(t,0).evalf()

woverv = np.divide(wprime_matrix, vprime_matrix)

########################### PLOTTING ###################################################################

########################### PLOT STREAMFUNCTION ###################################################################

#Plot using contourf
fig1 = plt.figure()

ax1 = fig1.add_subplot(221)
ax1.set_xlabel('x (zonal)')
ax1.set_ylabel('z (height)')
ax1.set_title('Streamfunction perturbation')
#ax1.imshow(phiprime_zxmatrix, origin = 'lower')
phicontour = plt.contourf(phiprime_zxmatrix, origin ='lower')
phicbar = plt.colorbar(phicontour)
#plt.clabel(phicontour, inline=1, fontsize=10)

########################### PLOT VELOCITY PERTURBATIONS ###################################################################

#Plot vprime using contourf
ax2 = fig1.add_subplot(222)
ax2.set_xlabel('x (zonal)')
ax2.set_ylabel('z (height)')
ax2.set_title('Meridional velocity perturbation')
#ax2.imshow(vprime_matrix, origin = 'lower')
vcontour = plt.contourf(vprime_matrix, origin = 'lower')
vcbar = plt.colorbar(vcontour)
#plt.clabel(vcontour, inline=1, fontsize=10)

#Plot wprime using contourf
ax3 = fig1.add_subplot(223)
ax3.set_xlabel('x (zonal)')
ax3.set_ylabel('z (height)')
ax3.set_title('Vertical velocity perturbation')
wcontour = plt.contourf(wprime_matrix, origin = 'lower')
wcbar = plt.colorbar(wcontour)
#plt.clabel(wcontour, inline=1, fontsize=10)

'''
#Plot wprime over vprime
ax4 = fig1.add_subplot(224)
ax4.set_xlabel('x (zonal)')
ax4.set_ylabel('z (height)')
ax4.set_title('w over v')
wovervcontour = plt.contour(woverv, origin = 'lower')
'''

########################### FIND PARCEL TRAJECTORIES ###################################################################

#STILL UNSURE WHAT I'M ACTUALLY PLOTTING HERE.

'''
def velocity(s,t): #vprime and wprime expressions found from Sympy, evaluated at Xlength/2
    y, z = s
    dsdt = [10000.0*(-8.50918128239322e-7*sin(1.0e-6*Xlength/2.)*sinh(0.0001*z) + 6.48054273663885e-7*cos(1.0e-6*Xlength/2.)*cosh(0.0001*z))*exp(1.549e-6*t), 5.0*z*(-8.50918128239322e-11*sin(1.0e-6*Xlength/2.)*cosh(0.0001*z) + 6.48054273663885e-11*cos(1.0e-6*Xlength/2.)*sinh(0.0001*z))*exp(1.549e-6*t) + 5.0*(-8.50918128239322e-7*sin(1.0e-6*Xlength/2.)*sinh(0.0001*z) + 6.48054273663885e-7*cos(1.0e-6*Xlength/2.)*cosh(0.0001*z))*exp(1.549e-6*t) - 0.01549*(6.48054273663885e-5*sin(1.0e-6*Xlength/2.)*sinh(0.0001*z) + 8.50918128239322e-5*cos(1.0e-6*Xlength/2.)*cosh(0.0001*z))*exp(1.549e-6*t)]
    return dsdt

s0 = [0, height/2.] #initial position in y-z plane

t = np.linspace(0,5e5,101) #create 101 timesteps up to ~6 days

sol = odeint(velocity, s0, t) #solve for y and z positions

ax4 = fig1.add_subplot(224)
ax4.plot(t, sol[:,0])
ax4.plot(t,sol[:,1])
ax4.contour(sol)
'''
plt.show()
