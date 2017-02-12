import numpy as np
from sympy import *
from RossbyInfo import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

###############################################################################################

#PRINT OUT PHYSICAL QUANTITIES FOR THE WAVE AND DIMENSIONAL SCALINGS

print 'phase speed, c (ms^(-1)) is:', c
print 'amplitude of uprime (in ms^(-1)), psi0*pi/L, is:', (psi0*np.pi)/L #amplitude of uprime
print 'amplitude of vprime (in ms^(-1)), psi0*k, is:', psi0*k #amplitude of vprime

print 'nondimensional amplitude of uprime_nd, (psi0*pi)/(L**2*f0), is:', (psi0*np.pi)/((L**2)*f0) #non-dimensional amplitude of uprime
print 'nondimensional amplitude of vprime_nd, (psi0*k)/(L*f0), is:', (psi0*k)/(L*f0) #non-dimensional amplitude of vprime

print 'Scaling factor to dimensionalise velocity equations (in ms^(-1)), L*f0 is:', L*f0

###############################################################################################

def main():
    xvalues = np.linspace(0,3,100)
    yvalues = np.linspace(0,1,100) #Channel width in y direction is 1*L
    tvalues = np.linspace(0,2,100)

    #print 'xvalues are:', xvalues

    xlength = len(xvalues)
    ylength = len(yvalues)
    tlength = len(tvalues)

    #CREATE ARRAY FOR PSIPRIME VALUES
    psiprime_nd = np.zeros((ylength,xlength,tlength)) #Note x and y have switched so plotted output is oriented correctly

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            for k in range(0, tlength, 1):
                #print 'x is: ', xvalues[j]
                #print 'y is: ', yvalues[i]
                #print 't is: ', tvalues[k]

                psiprime_nd[i,j,k] = psiprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'uprime_nd is (from function):', uprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'uprime_nd is (from array):', uprime_nd[i,j,k]

    print 'type of psiprime_nd is:', type(psiprime_nd)
    print 'shape of psiprime_nd is:', psiprime_nd.shape
    #print 'uprime_nd x values are:', uprime_nd[:,0,0]

###############################################################################################

    #CREATE ARRAY FOR UPRIME VALUES
    uprime_nd = np.zeros((ylength,xlength,tlength)) #Note x and y have switched so plotted output is oriented correctly

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            for k in range(0, tlength, 1):
                #print 'x is: ', xvalues[j]
                #print 'y is: ', yvalues[i]
                #print 't is: ', tvalues[k]

                uprime_nd[i,j,k] = uprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'uprime_nd is (from function):', uprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'uprime_nd is (from array):', uprime_nd[i,j,k]

    print 'type of uprime_nd is:', type(uprime_nd)
    print 'shape of uprime_nd is:', uprime_nd.shape
    #print 'uprime_nd x values are:', uprime_nd[:,0,0]

###############################################################################################

    #CREATE ARRAY FOR VPRIME VALUES
    vprime_nd = np.zeros((ylength,xlength,tlength))

    for i in range(0, ylength, 1): #Note x and y have switched so plotted output is oriented correctly
        for j in range(0, xlength, 1):
            for k in range(0, tlength, 1):
                #print 'x is: ', xvalues[j]
                #print 'y is: ', yvalues[i]
                #print 't is: ', tvalues[k]

                vprime_nd[i,j,k] = vprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'vprime_nd is (from function):', vprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'vprime_nd is (from array):', vprime_nd[i,j,k]

    print 'type of vprime_nd is:', type(vprime_nd)
    print 'shape of vprime_nd is:', vprime_nd.shape
    #print 'vprime_nd x values are:', vprime_nd[:,0,0]

    print 'vprime_nd at y = 0 is:', vprime_nd[2,0,2]
    print 'vprime_nd at y = L is:', vprime_nd[2,-1,2]

###############################################################################################

    #DEFINE ODE TO INTEGRATE FOR PARCEL TRAJECTORIES
    #???WHY DON'T THE PARCELS MOVE VERY FAR AT ALL??

    def velocitynd(s,t):
        x, y = s
        dsdt = [uprime(x=x,y=y,t=t), vprime(x=x,y=y,t=t)]
        return dsdt

    t=np.linspace(0,2,200)
    s0_a = [0.3, 0.25]
    s0_b = [0.3, 0.25001]
    s0_c = [0.5,0.5]

    sol_a = odeint(velocitynd, s0_a, t)
    sol_b = odeint(velocitynd, s0_b, t)
    sol_c = odeint(velocitynd, s0_c, t)

###############################################################################################

    fig1 = plt.figure()

    ax1=fig1.add_subplot(311)
    ax1.matshow(psiprime_nd[:,:,0], origin='lower')

    ax2 = fig1.add_subplot(312)
    ax2.matshow(uprime_nd[:,:,0], origin='lower')

    #print 'uprime_nd[:,:,0] is:', uprime_nd[:, : ,0]
    #print 'uprime_nd[:,0,:] is:', uprime_nd[:, 0, :]
    #print 'uprime_nd[0,:,:] is:', uprime_nd[0, :, :]

    ax3 = fig1.add_subplot(313)
    ax3.matshow(vprime_nd[:,:,0], origin='lower')

    #print 'vprime_nd[:,:,0] is:', vprime_nd[:, : ,0]
    #print 'vprime_nd[:,0,:] is:', vprime_nd[:, 0, :]
    #print 'vprime_nd[0,:,:] is:', vprime_nd[0, :, :]

    fig2 = plt.figure()
    ax4 = fig2.add_subplot(111)
    ax4.plot(sol_a[:,0], sol_a[:,1], label='trajectory a')
    ax4.plot(sol_b[:,0], sol_b[:,1], label='trajectory b')
    ax4.plot(sol_c[:,0], sol_c[:,1], label='trajectory c')
    plt.legend(loc='best')

    plt.show()

main()
