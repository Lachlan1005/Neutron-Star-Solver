import numpy as np
import matplotlib.pyplot as plt
hbar = 3.33 * 10**(-67)   # solar mass * km^2 / s
mN   = 8.395 * 10**(-58)  # solar masses
c    = 3 * 10**5          # km/s
G    = 1.327 * 10**(11)   # km^3/s^2 solar masses
Rsun = 696340             # km
def WDrk4(f:callable, g:callable, x:float, y:float, dt:float, tmax:float, gamma, a, b, name:str="Test System", epsilon0:float=3*10**(-2),graph:bool=True,verbosity:bool=True)->list[list[float]]:
    """
    Implement rk4 to solve for the function x(t) in the following system of equations

    dx/dt = f(t,x,y)
    dy/dt = g(t,x,y)

    Store the output in a list sol=[ [tvalues] ,[xvalues] , [yvalues] ] 
    """

    t=0
    x0=1.5*x
    tvals=[t]
    xvals=[x]
    yvals=[y]
    sol=[tvals,xvals,yvals]
    print("Solving for: "+name+". See Below for Initial Conditions")
    print("r: ",tvals[-1], " || M:",yvals[-1], " || p:", xvals[-1],"\n")
    if not verbosity:
        print("Solving in progress... \nDo not close or interrupt the program. ")
    while isinstance(x,float) and 0<=x<=x0:
        if verbosity:
            print("Solving "+name+": Iteration "+str(int(t/dt+1)))
        t+=dt
        tvals.append(t)
        xvals.append(x)
        yvals.append(y)
        k1=f(t,x,y,gamma, a, b,epsilon0)
        l1=g(t,x,y,gamma, a, b,epsilon0)
        k2=f(t+dt/2 , x+dt*k1/2, y+dt*l1/2,gamma, a, b,epsilon0)
        l2=g(t+dt/2 , x+dt*k1/2, y+dt*l1/2,gamma, a, b,epsilon0)
        k3=f(t+dt/2 , x+dt*k2/2, y+dt*l2/2,gamma, a, b,epsilon0)
        l3=g(t+dt/2 , x+dt*k2/2, y+dt*l2/2,gamma, a, b,epsilon0)
        k4=f(t+dt, x+dt*k3, y+dt*l3,gamma, a, b,epsilon0)
        l4=g(t+dt, x+dt*k3, y+dt*l3,gamma, a, b,epsilon0)
        x+=dt/6 * (k1+2*k2+2*k3+k4)
        y+=dt/6 * (l1+2*l2+2*l3+l4)
        if verbosity:
            print("M: ", y, " || p: ", x, "\n")
    if graph:
        print("Solving complete. See output plot for result. See below for resultant parameters.")
        print("r: ",round(tvals[-1]*1000)/1000, " || M:",round(yvals[-1]*1000)/1000, " || p:", round(xvals[-1]*1000)/1000)
        sol=[tvals,xvals,yvals]
        plt.rcParams['text.usetex'] = True
        fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(15,8))
        fig.suptitle('Polytrope Equations of State, $R$='+str(round(tvals[-1]))+" km, $M/M_{\odot}=$"+str(round((yvals[-1]*1000))/1000), fontsize=18)
        fig.text(0.5, 0.92, r"$\epsilon_0=$"+str(epsilon0)+r" $M_{\odot}c^2/km^3,\ p_0/\epsilon_0=$"+str(xvals[0]),
        fontsize=14, ha='center')
        ax1.plot(tvals, xvals, color="black")
        ax1.set_title('Normalised Pressure vs Radial Distance', fontsize=16)
        ax1.set_xlabel('$r$ (km)', fontsize=16)
        ax1.set_ylabel(r'$p/\epsilon_0$ (1)', fontsize=16)

        ax2.plot(tvals, yvals, color="black")
        ax2.set_title("Mass vs Radial Distance", fontsize=16)
        ax2.set_xlabel('$r$ (km)', fontsize=16)
        ax2.set_ylabel(r'$M/M_{\odot}$ (1)', color='black', fontsize=16)
        plt.show()
    else:
        print("Solving complete. See below for resultant parameters.")
        print("r: ",round(tvals[-1]*1000)/1000, " || M:",round(yvals[-1]*1000)/1000, " || p:", round(xvals[-1]*1000)/1000)
    return sol 

#Reccomended Parameters: a=1.473, b=52.46, gamma=4/3, p0= 10**(-16), dr=1 

def pressureChange(t,p,M,gamma, a, b,epsilon0):
    r=t
    kappa=G/(a*c**2)
    sigma=4*np.pi*epsilon0/c**2
    R0=G/c**2
    if p==0 or M==0:
        return p
    return -(((a*M*p**(1/gamma))*(1+kappa*p**(1-1/gamma))*(1+sigma*(r**3)*(p/M))))/((1-2*R0*M/r)*r**2)

def massChange(t,p,M,gamma, a, b,epsilon0):
     r=t
     return b*(r**2)*(p**(1/gamma))

def runSimulation(a:float,b:float,gamma:float,epsilon0:float, p0, dr, customAB:bool=False, A:float=1, Z:float=2,graph:bool=True, verbosity:bool=True):
    if customAB: #This overrides any a and b selection
        K=(hbar*c)/(12*np.pi**2) * ( (Z * 3*np.pi**2) / (mN * A * c**2) )**(4/3)
        a=Rsun/( (K * epsilon0**(gamma -1))**(1/gamma) )
        b= (4 * np.pi * epsilon0) / (c **2 * (K * epsilon0**(gamma -1))**(1/gamma))
    return WDrk4(pressureChange,massChange,p0,0,dr,epsilon0, gamma, a, b, "White Dwarf",10**(-12),graph, verbosity)

def terminalSimulation():
    print("===================================================")
    preset=int(input("Welcome to TOV Polytrope Solver. Enter 0 to view preset solution or insert any other number to solve from custom parameters.\nA unit system in km and solar masses will be used.\nYour input: ")) 
    print("===================================================")
    if preset==0:
        graph=False
        graphNum=int(input("Enter 1 to plot result. Enter any number to omit the plot. "))
        if graphNum==1:
            graph=True
        verbosity=False 
        verbosityNum=int(input("Enter 0 to activate debug mode. Enter any number to continue. "))
        if verbosityNum==0:
            verbosity=True
        runSimulation(1.473, 52.46, 4/3, 10**(-12), 10**(-16), .003, False,1,1,graph,verbosity)
    else:
        print("Enter the parameters that characterise your star. See documentation for common parameters for red giants, neutron stars, and more.")
        p0=float(input("Enter central normalised pressure: "))
        a=float(input("Enter coefficient of pressure equation (parameter a): "))
        b=float(input("Enter coefficient of mass equation (parameter b): "))
        gamma=float(input("Enter polytropic exponent: "))
        epsilon0=10**(-12)
        dr=float(input("Enter step size: "))
        graph=False
        graphNum=int(input("Enter 1 to plot result. Enter any number to omit the plot. "))
        if graphNum==1:
            graph=True
        verbosity=False 
        verbosityNum=int(input("Enter 0 to activate debug mode. Enter any number to continue. "))
        if verbosityNum==0:
            verbosity=True
        WDrk4(pressureChange,massChange,p0,0,dr,epsilon0, gamma, a, b, "Polytrope",10**(-12),graph, verbosity)

#runSimulation(1.473, .0488, 2, .1, 0.01, .00001, False,1,1,True,True)

#terminalSimulation()
#dataset=runSimulation(1.473, 52.46, 4/3, 10**(-12), 10**(-16), .005, False,1,1,True,True)





