import numpy as np
import matplotlib.pyplot as plt
hbar = 3.33 * 10**(-71)   # solar mass * km^2 / s
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
    #print("Solving for: "+name+". See Below for Initial Conditions")
    #print("r: ",tvals[-1], " || M:",yvals[-1], " || p:", xvals[-1],"\n")
    if not verbosity:
        print("  ")
    while isinstance(x,float) and 10**(-23.69)<=x<=x0:
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

def RvsM(debug:int, a:float,b:float,gamma:float, epsilon0:float, p0min:float, p0max:float,dr:float):
    """
    debug=0 => debug 
    debug=1 => standard 
    debug=2 => fancy  
    """
    radii=[]
    masses=[]
    i=1
    maxmass=0
    maxmassRad=0
    totali=np.log(p0max/p0min)/np.log(1.1)
    while p0min<=p0max:
        if debug != 2:
            print("\n-------------------------------\nSimulating Polytrope No.",i, " || p0=", p0min)
            if i*100/totali>100:
                print("\nSimulation Progress: ",100,"%")
            elif i*100/totali<0.0001:
                print("\nSimulation Progress: ",0.001,"%")
            else:
                print("\nSimulation Progress: ",i*100/totali,"%")
        if debug ==2:
            print(30*"\n")
            print("Solver running...")
            if i*100/totali<0.0001:
                print("Simulation Progress: ",0.001,"%")
            if i*100/totali>100:
                print("\nSimulation Progress: ",100,"%")
            else:
                print("Simulation Progress: ",np.round(i*1000/totali)/10,"%")
            print(int(i*100/totali)*"|")
            print("0%"+95*" "+"100%")
        dataset=runSimulation(a,b,gamma,epsilon0,p0min, dr, False, 1, 2, False, False)
        mass=dataset[2][-1]
        radius=dataset[0][-1]
        radii.append(radius)
        masses.append(mass)
        if debug !=2:
            print("Local Results: R=",radius, "|| M=", mass)
        if debug==0:
            plt.plot(radius, mass,"o",color="black")
        if mass>maxmass:
            maxmass=mass 
            maxmassRad=radius
        p0min+=p0min*0.1
        i+=1
    print("Simulation Complete. See resultant plot for solution.",i,"stars simulated.")
    print("The specifications of the heaviest star is as follows: M=",maxmass, ", R=",maxmassRad)
    if debug != 0:
        plt.plot(radii, masses, color="black")
    plt.title("Mass vs Radius of Polytrope")
    plt.xlabel("Radii (km)")
    plt.ylabel("Mass in Solar Masses (1)")
    plt.show()

#White Dwarf Relativistic
#RvsM(0, 1.473, 52.46, 4/3, 10**(-8), 10**(-20), 10**(-12),.1)


#White Dwarf Nonrelativistic
#RvsM(2, 0.05, 0.005924, 5/3, 10**(-8), 10**(-5), 10**(3),.01)

#Neutron Star1


def terminalSimulation():
    print(100*"=")
    preset=int(input("Welcome to TOV Equation of State Solver. Enter 0 to view preset solution or insert any other number to solve from custom parameters.\nA unit system in km and solar masses will be used.\nYour input: ")) 
    print(100*"=")
    if preset==0:
        graph=False
        RvsM(2, 1.473, 0.0488, 2, 10**(-2),10**(-16.9), 10**(8),.001)
    else:
        print("Enter the parameters that characterise your sample of stars. See documentation for common parameters for red giants, neutron stars, and more.")
        p0min=float(input("Enter minimum central normalised pressure: "))
        p0max=float(input("Enter minimum central normalised pressure: "))
        a=float(input("Enter coefficient of pressure equation (parameter a): "))
        b=float(input("Enter coefficient of mass equation (parameter b): "))
        gamma=float(input("Enter polytropic exponent (parameter gamma): "))
        epsilon0=10**(-12)
        dr=float(input("Enter radial step size (parameter dr): "))
    return RvsM(2,a,b,gamma,epsilon0,p0min,p0max,dr)



