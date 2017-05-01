import numpy as np
import cmath

#I DON'T KNOW WHAT FRAME WE'RE IN. CURRENT FORMULATION MAKES c EITHER
#PURELY IMAGINARY OR PURELY REAL.

#RE-DO THIS WHOLE PROGRAM FOLLOWING VALLIS. GILL IS CONFUSING.

################################# LIST OF PHYSICAL PARAMETERS #############################################################

#Ocean parameters for Eady model

f0 = -1e-4 #(s^-1) f-plane approx. for Coriolis parameter (southern hemisphere)
N = 1e-2 #(s^-1) buoyancy frequency
H = 1e3 #(m) positions of boundaries at H and -H
U = 0.1 #(ms^(-1)) mean flow zonal velocity magnitude at boundaries H and -H
shear = U/H #(s^-1) velocity shear
L = np.abs((N*H)/f0) #(m) typical length scale given by deformation radius
T = L/U #Eady timescale
k = 1e-5 #(m^-1) zonal wavenumber
#k = 0.8031/L #for most unstable case, Gill 13.3.12
l = (np.pi)/(2*L) #(m^-1) meridional wavenumber of growing wave, defined for zonal channel
#l = 0 #for most unstable case, Gill 13.3.12
H_R = np.abs(f0/(N*k)) #(m) Rossby height
Hratio = np.abs(H/H_R)
alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration

#Phase speed of wave, Gill 13.3.5
csquared = (shear**2)*(H*np.tanh(Hratio) - H_R)*(H*(1./np.tanh(Hratio)) - H_R)
c = cmath.sqrt(csquared)

sigma = k*c.imag #growth rate

print 'c is:', c
print 'sigma is:', sigma

A_mag = 1.
B_mag = 1.

#Gill 13.3.4
def A(x,y,t):
    return A_mag*np.cos(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(sigma*T*t)

def B(x,y,t):
    return B_mag*np.cos(l*L*y)*np.sin(k*(L*x - c.real*T*t))*np.exp(sigma*T*t)

#Gill 13.3.2
def streamfunction(x,y,z,t): #Has dimensions m^(2)s^(-1) - HOW DOES THIS AFFECT THE NON-DIMENSIONALISATION??
    return A(x,y,t)*np.sinh(z*Hratio) + B(x,y,t)*np.cosh(z*Hratio)

def streamfunction_x(x,y,z,t):
    return (-k*A_mag*np.cos(l*L*y)*np.sin(k*(L*x - c.real*T*t))*np.exp(sigma*T*t))*np.sinh(z*Hratio) + (k*B_mag*np.cos(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(sigma*T*t))*np.cosh(z*Hratio)

def streamfunction_y(x,y,z,t):
    return (-l*A_mag*np.sin(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(sigma*T*t))*np.sinh(z*Hratio) + (-l*B_mag*np.sin(l*L*y)*np.sin(k*(L*x - c.real*T*t))*np.exp(sigma*T*t))*np.cosh(z*Hratio)

def streamfunction_z(x,y,z,t):
    return (1./H_R)*(A(x,y,t)*np.cosh(z*Hratio) + B(x,y,t)*np.sinh(z*Hratio))

def streamfunction_zt(x,y,z,t):
    return sigma*streamfunction_z(x,y,z,t)

def streamfunction_zx(x,y,z,t):
    return (1./H_R)*((-k*A_mag*np.cos(l*L*y)*np.sin(k*(L*x - c.real*T*t))*np.exp(sigma*T*t))*np.cosh(z*Hratio) + (k*B_mag*np.cos(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(sigma*T*t))*np.sinh(z*Hratio))

def umeanflow(z): #Non-dimensional zonal velocity mean shear flow
    return (1./(L*np.abs(f0)))*shear*H*z

def uprime(x,y,z,t): #Non-dimensional zonal velocity perturbation
    return np.abs(1./(L*f0))*(-1./f0)*streamfunction_y(x,y,z,t)

def vprime(x,y,z,t): #Non-dimensional meridional velocity perturbation
    return np.abs(1./(L*f0))*(1./f0)*streamfunction_x(x,y,z,t)

def wprime(x,y,z,t): #Non-dimensional vertical velocity perturbation, GILL 12.9.6
    return (1./(H*np.abs(f0)))*(1./(N*N))*(shear*streamfunction_x(x,y,z,t) - streamfunction_zt(x,y,z,t) - H*z*shear*streamfunction_zx(x,y,z,t))
