import numpy as np
from sympy import *
from RossbyInfo import *
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
    xvalues = np.floor(np.arange(0., 3*L, 5e3))
    yvalues = np.floor(np.arange(0., L, 5e3))

    Xlength = len(xvalues)
    Ylength = len(yvalues)

    print 'Xlength is:', Xlength
    print 'Ylength is:', Ylength

    #Define streamfunction
    x, y, t = symbols('x y t')

    #psiprime = psi0*cos(k*L*x - c*T*t)*sin(np.pi*y) #streamfunction perturbation (expressed with non-dimensional variables)
    psiprime = psi0*cos(k*x - c*t)*sin((np.pi*y)/L) #dimensional streamfunction perturbation

    psiprime_x = diff(psiprime, x) #partial derivatives found using Sympy
    psiprime_y = diff(psiprime, y)

    #non-dimensional velocity perturbations
    uprime_nd = -(1./(L*f0))*psiprime_y
    vprime_nd = (1./(L*f0))*psiprime_x

    print 'type of psiprime:', type(psiprime)
    print 'psiprime is:', psiprime
    print 'psiprime_x is:', psiprime_x
    print 'psiprime_y is:', psiprime_y

    print 'uprime_nd is:', uprime_nd
    print 'vprime_nd is:', vprime_nd
    print 'vprime_nd at y=0 is:', vprime_nd.subs(y,0).subs(x,1).subs(t,1).evalf()
    print 'vprime_nd at y=1 is:', vprime_nd.subs(y,1.).subs(x,1).subs(t,1).evalf()


    #Initialise matrix of zeros for uprime_nd values in x-y plane
    uprime_matrix = np.zeros((Ylength,Xlength))

    #Evaluate values for uprime matrix at time t=1
    for i in range(0, Ylength, 1):
        for j in range(0, Xlength, 1):
            uprime_matrix[i-1][j-1] = uprime_nd.subs(x,xvalues[j-1]).subs(y,yvalues[i-1]).subs(t,1.).evalf()

    print 'uprime matrix is:', uprime_matrix

    #Initialise matrix of zeros for vprime_nd values in x-y plane
    vprime_matrix = np.zeros((Ylength,Xlength))

    #Evaluate values for vprime_nd matrix at time t=1
    for i in range(0, Ylength, 1):
        for j in range(0, Xlength, 1):
            vprime_matrix[i-1][j-1] = vprime_nd.subs(x,xvalues[j-1]).subs(y,yvalues[i-1]).subs(t,1.).evalf()

    print 'vprime matrix is:', vprime_matrix

    #Define ODE to plot parcel trajectories, expressions for uprime_nd and vprime_nd taken from printed output
    def nondimvelocity(s,t):
        x, y = s
        #dsdt = [-0.009*cos(np.pi*y)*cos(270.0*t + 2.0*x), -0.006*sin(np.pi*y)*sin(270.0*t + 2.0*x)]
        dsdt = [-0.009*cos(3.14159265358979e-5*y)*cos(0.027*t + 2.0e-5*x), -0.006*sin(3.14159265358979e-5*y)*sin(0.027*t + 2.0e-5*x)]
        return dsdt

    t = np.linspace(0,1e3,1e2) #non-dimensional t, in units of f0

    s0_a = [0.1,0.5]
    s0_b = [0.1,0.50001]
    s0_c = [0.1,0.50002]
    s0_d = [0.1,0.50003]
    s0_e = [0.1,0.50004]

    #multiply solutions by L to get physical trajectories
    sol_a = L*odeint(nondimvelocity, s0_a, t)
    sol_b = L*odeint(nondimvelocity, s0_b, t)
    sol_c = L*odeint(nondimvelocity, s0_c, t)
    sol_d = L*odeint(nondimvelocity, s0_d, t)
    sol_e = L*odeint(nondimvelocity, s0_e, t)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    #ax1.plot(t, sol_a[:,0], label='trajectory a')
    ax1.plot(sol_a[:,0], sol_a[:,1], label='trajectory a')
    ax1.plot(sol_b[:,0], sol_b[:,1], label='trajectory b')
    ax1.plot(sol_c[:,0], sol_c[:,1], label='trajectory c')
    ax1.plot(sol_d[:,0], sol_d[:,1], label='trajectory d')
    ax1.plot(sol_e[:,0],sol_e[:,1], label='trajectory e')
    plt.legend(loc='best')


    fig2 = plt.figure()
    ax2 = fig2.add_subplot(211)
    #ax2.contour(uprime_matrix, origin='lower')
    ax2.matshow(uprime_matrix, origin='lower')

    ax3 = fig2.add_subplot(212)
    #ax3.contour(vprime_matrix, origin='lower')
    ax3.matshow(vprime_matrix, origin='lower')

    plt.show()

main()
