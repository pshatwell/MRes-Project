import numpy as np
from sympy import *
from RossbyInfo import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

###############################################################################################

#PRINT OUT PHYSICAL QUANTITIES FOR THE WAVE AND DIMENSIONAL SCALINGS

print 'phase speed, c (ms^(-1)) is:', c
print 'horizontal wavelength (m) is:', wavelength

print 'amplitude of uprime (in ms^(-1)), psi0*pi/L, is:', (psi0*np.pi)/L #amplitude of uprime
print 'amplitude of vprime (in ms^(-1)), psi0*k, is:', psi0*k #amplitude of vprime

print 'Scaling factor to dimensionalise velocity equations (in ms^(-1)), L*f0 is:', L*f0

print 'uprime_nd is: -((psi0*np.pi)/((L**2)*f0))*np.cos(k*L*x - c*T*t)*np.cos(np.pi*y)'
print 'vprime_nd is: -((psi0*k)/(L*f0))*np.sin(k*L*x - c*T*t)*np.sin(np.pi*y)'

print 'amplitude of uprime_nd, (psi0*pi)/((L**2)*f0), is:', A_und
print 'amplitude of vprime_nd, (psi0*k)/(L*f0), is:', A_vnd

print 'k is:', k
print 'L is:', L
print 'k*L is:', k*L

print 'c is:', c
print 'T is:', T
print 'k*c*T is:', k*c*T

###############################################################################################

def main():
    #Define dimensions of arrays
    xvalues = np.linspace(0,1,100)
    yvalues = np.linspace(0,1,100) #Channel width in y direction is 1*L

    time = 0 #Choose time step (in units of T)

    xlength = len(xvalues)
    ylength = len(yvalues)

###############################################################################################

    #CREATE ARRAY FOR PSIPRIME VALUES
    psiprime_nd = np.zeros((ylength,xlength)) #Note x and y have switched so plotted output is oriented correctly

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            psiprime_nd[i,j] = psiprime(x=xvalues[j], y=yvalues[i], t=time)

    print 'type of psiprime_nd is:', type(psiprime_nd)
    print 'shape of psiprime_nd is:', psiprime_nd.shape

###############################################################################################

    #CREATE ARRAY FOR UPRIME VALUES
    uprime_nd = np.zeros((ylength,xlength)) #Note x and y have switched so plotted output is oriented correctly

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            uprime_nd[i,j] = uprime(x=xvalues[j], y=yvalues[i], t=time)

    print 'type of uprime_nd is:', type(uprime_nd)
    print 'shape of uprime_nd is:', uprime_nd.shape

###############################################################################################

    #CREATE ARRAY FOR VPRIME VALUES
    vprime_nd = np.zeros((ylength,xlength))

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            vprime_nd[i,j] = vprime(x=xvalues[j], y=yvalues[i], t=time)


    print 'type of vprime_nd is:', type(vprime_nd)
    print 'shape of vprime_nd is:', vprime_nd.shape

    #Check boundary conditions for v
    print 'vprime_nd at y = 0 is:', vprime_nd[0,10]
    print 'vprime_nd at y = L is:', vprime_nd[-1,10]

###############################################################################################

    #DEFINE ODE TO INTEGRATE FOR PARCEL TRAJECTORIES

    #SOMETHING VERY WRONG WITH ORDERS OF MAGNITUDE SOMEWHERE? AMPLITUDES ARE ~1e-4, cT~1e2,
    #WHICH MEANS 1e-6 SCALING OF TRAJECTORIES. VERY SMALL LENGTHS, VERY FAST OSCILLATIONS.
    #WHAT PARAMETERS ARE INCORRECT HERE?

    #MAYBE ODEINT IS THE WRONG INTEGRATOR TO USE?? ARE U AND V ACTUALLY UPDATED WITH EACH TIMESTEP?

    def velocity(s,t):
        x, y = s
        dsdt = [uprime(x=x,y=y,t=t), vprime(x=x,y=y,t=t)]
        #dsdt = [testuprime(x=x,y=y,t=t), testvprime(x=x,y=y,t=t)]
        return dsdt

    t=np.linspace(0,100,500)

    #s0_a = [0.47, 0.1] #oscillatory trajectory
    #s0_e = [0.7857, 0.5] #oscillatory trajectory
    s0_a = [0.47, 0.1]
    s0_b = [0.5, 0.5]

    '''
    s0_c = [0.47, 0.12]
    s0_d = [0.47, 0.13]
    s0_e = [0.47, 0.14]
    s0_f = [0.47, 0.15]
    '''

    sol_a = odeint(velocity, s0_a, t, mxstep=1000000)
    sol_b = odeint(velocity, s0_b, t, mxstep=1000000)
    #sol_c = odeint(velocity, s0_c, t, mxstep=1000000)
    #sol_d = odeint(velocity, s0_d, t, mxstep=1000000)
    #sol_e = odeint(velocity, s0_e, t, mxstep=1000000)
    #sol_f = odeint(velocity, s0_f, t, mxstep=1000000)

    print 'type of trajectory solution is:', type(sol_a)
    print 'shape of trajectory solution is:', sol_a.shape


###############################################################################################

    #ATTEMPT TO TRANSFORM TRAJECTORIES TO WAVE REFERENCE FRAME
    shift = np.zeros_like(sol_a)
    shift[:,0] = c*t

    print 'shape of sol_a is:', sol_a.shape
    print 'shape of shift is:', shift.shape
    print 'shape of transposed shift is:', shift.transpose().shape

    rel_sol_a = sol_a - shift #new parcel trajectory relative to Rossby wave
    rel_sol_b = sol_b - shift

###############################################################################################

    #DO ALL THE PLOTTING
    fig1 = plt.figure()

    ax1=fig1.add_subplot(131)
    ax1.imshow(psiprime_nd[:,:], origin='lower')

    ax2 = fig1.add_subplot(132)
    ax2.imshow(uprime_nd[:,:], origin='lower')

    ax3 = fig1.add_subplot(133)
    ax3.imshow(vprime_nd[:,:], origin='lower')


    fig2 = plt.figure()
    ax4 = fig2.add_subplot(111)

    #Multiplying solutions by L to get physical lengths
    ax4.plot(L*sol_a[:,0], L*sol_a[:,1], label='a')
    ax4.plot(L*sol_b[:,0], L*sol_b[:,1], label='b')
    #ax4.plot(L*sol_c[:,0], L*sol_c[:,1], label='c')
    #ax4.plot(L*sol_d[:,0], L*sol_d[:,1], label='d')
    #ax4.plot(L*sol_e[:,0], L*sol_e[:,1], label='e')
    #ax4.plot(L*sol_f[:,0], L*sol_f[:,1], label='f')

    ax4.plot(L*rel_sol_a[:,0], L*rel_sol_a[:,1], color='black', label ='a in RW frame')
    ax4.plot(L*rel_sol_b[:,0], L*rel_sol_b[:,1], color='red', label ='b in RW frame')

    #Overlay initial streamfunction to see if trajectories make sense.
    #ax4.imshow(psiprime_nd[:,:], origin='lower', extent=[0,L,0,L])

    plt.legend(loc='best')

    plt.show()

main()
