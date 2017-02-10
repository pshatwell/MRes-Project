import numpy as np
from sympy import *
from RossbyParams import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

print 'phase speed, c (ms^(-1)) is:', c
print 'amplitude of uprime (in ms^(-1)), psi0*pi/L, is:', (psi0*np.pi)/L #amplitude of uprime
print 'amplitude of vprime (in ms^(-1)), psi0*k, is:', psi0*k #amplitude of vprime

print 'nondimensional amplitude of uprime, (psi0*pi)/(L**2*f0), is:', (psi0*np.pi)/((L**2)*f0) #non-dimensional amplitude of uprime
print 'nondimensional amplitude of vprime, (psi0*k)/(L*f0), is:', (psi0*k)/(L*f0) #non-dimensional amplitude of vprime

print 'Scaling factor to dimensionalise equations (in ms^(-1)), L*f0 is:', L*f0

def main():

    #Define coords for plotting
    xvalues = np.floor(np.arange(0, L, 1e3))
    yvalues = np.floor(np.arange(0, L, 1e3))

    Xlength = len(xvalues)
    Ylength = len(yvalues)

    #Define streamfunction
    x, y, t = symbols('x y t') #all non-dimensional

    psiprime = psi0*cos(k*L*x - c*T*t)*sin(np.pi*y) #streamfunction perturbation (expressed with non-dimensional variables)

    psiprime_x = diff(psiprime, x) #partial derivatives found using Sympy
    psiprime_y = diff(psiprime, y)

    #velocity perturbations
    uprime_d = -(1./L)*psiprime_y
    vprime_d = (1./L)*psiprime_x

    #non-dimensional velocity perturbations
    uprime_nd = -(1./((L**2)*f0))*psiprime_y
    vprime_nd = (1./((L**2)*f0))*psiprime_x

    print 'type of psiprime:', type(psiprime)
    print 'psiprime is:', psiprime
    print 'psiprime_x is:', psiprime_x
    print 'psiprime_y is:', psiprime_y
    print 'uprime_d is:', uprime_d
    print 'vprime_d is:', vprime_d
    print 'uprime_nd is:', uprime_nd
    print 'vprime_nd is:', vprime_nd

    '''
    #Initialise matrix of zeros for uprime_nd values in x-y plane
    uprime_matrix = np.zeros((Ylength,Xlength))

    #Evaluate values for uprime matrix at time t=0
    for i in range(0, Ylength, 1):
        for j in range(0, Xlength, 1):
            uprime_matrix[i-1][j-1] = uprime_d.subs(x,xvalues[j-1]).subs(y,yvalues[i-1]).subs(t,0).evalf()

    #Initialise matrix of zeros for uprime_nd values in x-y plane
    vprime_matrix = np.zeros((Ylength,Xlength))

    #Evaluate values for vprime matrix at time t=0
    for i in range(0, Ylength, 1):
        for j in range(0, Xlength, 1):
            vprime_matrix[i-1][j-1] = vprime_d.subs(x,xvalues[j-1]).subs(y,yvalues[i-1]).subs(t,0).evalf()
    '''

    #Define ODE to plot parcel trajectories, expressions for uprime_nd and vprime_nd taken from printed output
    def nondimvelocity(s,t):
        x, y = s
        dsdt = [-0.00942477796076938*cos(3.14159265358979*y)*cos(270.0*t + 2.0*x), -0.006*sin(3.14159265358979*y)*sin(270.0*t + 2.0*x)]
        return dsdt

    t = np.linspace(0,0.2,1e3) #non-dimensional t, in units of f0

    s0_a = [0.2,0.5] #non-dimensional starting positions, in units of 1/L
    s0_b = [0.2,0.50001]
    s0_c = [0.2,0.50002]
    s0_d = [0.2,0.50003]
    #s0_e = [0.4005,0.5]

    sol_a = L*odeint(nondimvelocity, s0_a, t)
    sol_b = L*odeint(nondimvelocity, s0_b, t)
    sol_c = L*odeint(nondimvelocity, s0_c, t)
    sol_d = L*odeint(nondimvelocity, s0_d, t)
    #sol_e = odeint(nondimvelocity, s0_e, t)

    fig = plt.figure()
    #ax1 = fig.add_subplot(311)

    ax1 = fig.add_subplot(111)

    #ax1.plot(t, sol_a[:,0], label='trajectory a')
    ax1.plot(sol_a[:,0], sol_a[:,1], label='trajectory a')
    ax1.plot(sol_b[:,0], sol_b[:,1], label='trajectory b')
    ax1.plot(sol_c[:,0], sol_c[:,1], label='trajectory c')
    ax1.plot(sol_d[:,0], sol_d[:,1], label='trajectory d')
    #ax1.plot(sol_e[:,0],sol_e[:,1], label='trajectory e')
    plt.legend(loc='best')

    '''
    ax2 = fig.add_subplot(312)
    ax2.contour(uprime_matrix, origin='lower')

    ax3 = fig.add_subplot(313)
    ax3.contour(vprime_matrix, origin='lower')
    '''

    plt.show()

main()
