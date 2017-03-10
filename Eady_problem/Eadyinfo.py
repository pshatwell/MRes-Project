import numpy as np
import cmath

################################# LIST OF PHYSICAL PARAMETERS #############################################################

#Ocean parameters for Eady model

L = 1e4 #(m) typical length scale
f0 = -1e-4 #(s^-1) f-plane approx. for Coriolis parameter (southern hemisphere)
T = 10*np.abs(1./f0) #(s) typical time scale
N = 1e-3 #(s^-1) buoyancy frequency
k = 1e-4 #(m^-1) zonal wavenumber of growing wave
l = 1e-4 #(m^-1) meridional wavenumber of growing wave
H_R = np.abs(f0/(N*k)) #(m) Rossby height
H = 0.01*H_R #(m) positions of boundaries at H and -H for most unstable wave
Hratio = np.abs(H/H_R)
U = 0.1 #(ms^(-1)) mean flow zonal velocity magnitude at boundaries H and -H
shear = U/H #(s^-1) velocity shear
alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration

#phase speed of wave
csquared = (shear**2)*(H*np.tanh(Hratio) - H_R)*(H*(1./np.tanh(Hratio)) - H_R)
c = cmath.sqrt(csquared)

sigma = k*c.imag #growth rate

print 'c is:', c
print 'sigma is:', sigma

A_mag = 1.
B_mag = 1.

#change l wavenumber so it's zero at boundaries for channel piy/L etc.
def A(x,y,t):
    return A_mag*np.cos(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(k*c.imag*T*t)

def B(x,y,t):
    return B_mag*np.cos(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(k*c.imag*T*t)

def streamfunction(x,y,z,t):
    return A(x,y,t)*np.sinh(z*Hratio) + B(x,y,t)*np.cosh(z*Hratio)
