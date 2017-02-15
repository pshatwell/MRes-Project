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

print 'Scaling factor to dimensionalise velocity equations (in ms^(-1)), L*f0 is:', L*f0

print 'uprime_nd is: -((psi0*np.pi)/((L**2)*f0))*np.cos(k*L*x - c*T*t)*np.cos(np.pi*y)'
print 'vprime_nd is: -((psi0*k)/(L*f0))*np.sin(k*L*x - c*T*t)*np.sin(np.pi*y)'

print 'amplitude of uprime_nd, (psi0*pi)/(L**2*f0), is:', A_und
print 'amplitude of vprime_nd, (psi0*k)/(L*f0), is:', A_vnd

print 'k is:', k
print 'L is:', L
print 'k*L is:', k*L

print 'c is:', c
print 'T is:', T
print 'c*T is:', c*T

###############################################################################################

def main():
    xvalues = np.linspace(0,3,100)
    yvalues = np.linspace(0,1,100) #Channel width in y direction is 1*L
    tvalues = np.linspace(0,1,100)

    #print 'xvalues are:', xvalues

    xlength = len(xvalues)
    ylength = len(yvalues)
    tlength = len(tvalues)

    #CREATE ARRAY FOR PSIPRIME VALUES
    psiprime_nd = np.zeros((ylength,xlength,tlength)) #Note x and y have switched so plotted output is oriented correctly

    #Indices in array go [y,x,t]
    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            for m in range(0, tlength, 1):
                #print 'x is: ', xvalues[j]
                #print 'y is: ', yvalues[i]
                #print 't is: ', tvalues[k]

                psiprime_nd[i,j,m] = psiprime(x=xvalues[j], y=yvalues[i], t=tvalues[m])
                #print 'uprime_nd is (from function):', uprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'uprime_nd is (from array):', uprime_nd[i,j,k]

    print 'type of psiprime_nd is:', type(psiprime_nd)
    print 'shape of psiprime_nd is:', psiprime_nd.shape
    #print 'uprime_nd x values are:', uprime_nd[:,0,0]

###############################################################################################

    #CREATE ARRAY FOR UPRIME VALUES
    uprime_nd = np.zeros((ylength,xlength,tlength)) #Note x and y have switched so plotted output is oriented correctly

    #Indices in array go [y,x,t]
    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            for m in range(0, tlength, 1):
                #print 'x is: ', xvalues[j]
                #print 'y is: ', yvalues[i]
                #print 't is: ', tvalues[k]

                uprime_nd[i,j,m] = uprime(x=xvalues[j], y=yvalues[i], t=tvalues[m])
                #print 'uprime_nd is (from function):', uprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'uprime_nd is (from array):', uprime_nd[i,j,k]

    print 'type of uprime_nd is:', type(uprime_nd)
    print 'shape of uprime_nd is:', uprime_nd.shape
    #print 'uprime_nd x values are:', uprime_nd[:,0,0]

###############################################################################################

    #CREATE ARRAY FOR VPRIME VALUES
    vprime_nd = np.zeros((ylength,xlength,tlength))

    #Indices in array go [y,x,t]
    for i in range(0, ylength, 1): #Note x and y have switched so plotted output is oriented correctly
        for j in range(0, xlength, 1):
            for m in range(0, tlength, 1):
                #print 'x is: ', xvalues[j]
                #print 'y is: ', yvalues[i]
                #print 't is: ', tvalues[k]

                vprime_nd[i,j,m] = vprime(x=xvalues[j], y=yvalues[i], t=tvalues[m])
                #print 'vprime_nd is (from function):', vprime(x=xvalues[j], y=yvalues[i], t=tvalues[k])
                #print 'vprime_nd is (from array):', vprime_nd[i,j,k]

    print 'type of vprime_nd is:', type(vprime_nd)
    print 'shape of vprime_nd is:', vprime_nd.shape
    #print 'vprime_nd x values are:', vprime_nd[:,0,0]

    print 'vprime_nd at y = 0 is:', vprime_nd[0,2,2]
    print 'vprime_nd at y = L is:', vprime_nd[-1,2,2]

    '''
    #Try to find streamlines by dividing v over u
    voveru = np.divide(vprime_nd, uprime_nd)
    print 'type of voveru is:', type(voveru)
    print 'shape of voveru is:', voveru.shape
    '''

###############################################################################################

    #DEFINE ODE TO INTEGRATE FOR PARCEL TRAJECTORIES

    #SOMETHING VERY WRONG WITH ORDERS OF MAGNITUDE SOMEWHERE. AMPLITUDES ARE ~1e-3, cT~1e2,
    #WHICH MEANS 1e-5 SCALING OF TRAJECTORIES. VERY SMALL LENGTHS, VERY FAST OSCILLATIONS.
    #WHAT PARAMETERS ARE INCORRECT HERE?

    #TRY SOLVING THIS ANALYTICALLY BY HAND AND PLOTTING THAT

    #CURRENTLY SOLVING UPRIME AND VPRIME WITH AMPLITUDES = 1.
    def velocitynd(s,t):
        x, y = s
        dsdt = [uprime(x=x,y=y,t=t), vprime(x=x,y=y,t=t)]
        #dsdt = [testuprime(x=x,y=y,t=t), testvprime(x=x,y=y,t=t)]
        return dsdt

    t=np.linspace(0,0.3,500)

    s0_a = [0, 0.5]
    s0_b = [0, 0.4999]

    print 'trajectory a initial x position:', A_und*L*s0_a[0]
    print 'trajectory a initial y position:', A_vnd*L*s0_a[1]

    print 'trajectory b initial x position:', A_und*L*s0_b[0]
    print 'trajectory b initial y position:', A_vnd*L*s0_b[1]

    print 'initial separation is:', np.sqrt((A_und*L*s0_a[0]-A_und*L*s0_b[0])**2 + (A_vnd*L*s0_a[1]-A_vnd*L*s0_b[1])**2)

    sol_a = odeint(velocitynd, s0_a, t, mxstep=1000000)
    sol_b = odeint(velocitynd, s0_b, t, mxstep=1000000)

###############################################################################################

    #ATTEMPT TO TRANSFORM TRAJECTORIES TO EDDY REFERENCE FRAME?

    shift_a = np.array(t*(uprime(x=sol_a[:,0], y=sol_a[:,1], t=t), vprime(x=sol_a[:,0], y=sol_a[:,1], t=t)))

    print 'shape of sol_a is:', sol_a.shape
    print 'shape of shift is:', shift_a.shape
    print 'shape of transposed shift is:', shift_a.transpose().shape

    #N.B. need to transpose to agree with shape of sol_a array
    rel_sol_a = sol_a - shift_a.transpose() #new parcel trajectory relative to WHAT???


    shift_b = np.array(t*(uprime(x=sol_b[:,0], y=sol_b[:,1], t=t), vprime(x=sol_b[:,0], y=sol_b[:,1], t=t)))
    rel_sol_b = sol_b - shift_b.transpose()

    print 'trajectory a final x position:', A_und*L*rel_sol_a[-1,0]
    print 'trajectory a final y position:', A_vnd*L*rel_sol_a[-1,1]

    print 'trajectory b final x position:', A_und*L*rel_sol_b[-1,0]
    print 'trajectory b final y position:', A_vnd*L*rel_sol_b[-1,1]

    print 'final separation is:', np.sqrt((A_und*L*rel_sol_a[-1,0]-A_und*L*rel_sol_b[-1,0])**2 + (A_vnd*L*rel_sol_a[-1,1]-A_vnd*L*rel_sol_b[-1,1])**2)

###############################################################################################

    #DO ALL THE PLOTTING

    fig1 = plt.figure()

    ax1=fig1.add_subplot(131)
    ax1.imshow(psiprime_nd[:,:,0], origin='lower')

    ax2 = fig1.add_subplot(132)
    ax2.imshow(uprime_nd[:,:,0], origin='lower')

    #print 'uprime_nd[:,:,0] is:', uprime_nd[:, : ,0]
    #print 'uprime_nd[:,0,:] is:', uprime_nd[:, 0, :]
    #print 'uprime_nd[0,:,:] is:', uprime_nd[0, :, :]

    ax3 = fig1.add_subplot(133)
    ax3.imshow(vprime_nd[:,:,0], origin='lower')

    #print 'vprime_nd[:,:,0] is:', vprime_nd[:, : ,0]
    #print 'vprime_nd[:,0,:] is:', vprime_nd[:, 0, :]
    #print 'vprime_nd[0,:,:] is:', vprime_nd[0, :, :]

    fig2 = plt.figure()
    ax4 = fig2.add_subplot(111)

    #MULTIPLYING SOLUTIONS BY THE VELOCITY AMPLITUDES AND BY L. MAKE SURE THIS IS CORRECT.
    ax4.plot(A_und*L*sol_a[:,0], A_vnd*L*sol_a[:,1], label='trajectory a')
    ax4.plot(A_und*L*sol_b[:,0], A_vnd*L*sol_b[:,1], label='trajectory b')

    ax4.plot(A_und*L*rel_sol_a[:,0], A_vnd*L*rel_sol_a[:,1], color='black', label ='trajectory a in different frame')
    ax4.plot(A_und*L*rel_sol_b[:,0], A_vnd*L*rel_sol_b[:,1], color='red', label ='trajectory b in different frame')

    #ax4.imshow(psiprime_nd[:,:,50], origin='lower', extent=[-3*L,0,0,L])

    plt.legend(loc='best')

    '''
    #Plot vprime over uprime
    fig3 = plt.figure()
    ax5 = fig3.add_subplot(111)
    ax5.set_xlabel('x')
    ax5.set_ylabel('y')
    ax5.set_title('v over u')
    voverucontour = plt.contour(voveru[:,:,0], origin = 'lower')
    '''

    plt.show()

main()
