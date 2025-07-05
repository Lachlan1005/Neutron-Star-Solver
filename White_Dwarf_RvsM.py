import numpy as np
import matplotlib.pyplot as plt
hbar = 3.33 * 10**(-67)   # solar mass * km^2 / s
mN   = 8.395 * 10**(-58)  # solar masses
c    = 3 * 10**5          # km/s
Rsun = 696340             # km
def WDrk4(f:callable, g:callable, x:float, y:float, dt:float, tmax:float, gamma, a, b, name:str="Test System", epsilon0:float=10**(-12),graph:bool=True)->list[list[float]]:
    """
    Implement rk4 to solve for the function x(t) in the following system of equations

    dx/dt = f(t,x,y)
    dy/dt = g(t,x,y)

    Store the output in a list sol=[ [tvalues] ,[xvalues] , [yvalues] ] 
    """
    t=0
    tvals=[t]
    xvals=[x]
    yvals=[y]
   # print("M: ", y, " || p: ", x, "\n")
    while isinstance(x,float) and x>0:
       # print("Solving "+name+": Iteration "+str(int(t/dt+1)))
        t+=dt
        tvals.append(t)
        xvals.append(x)
        yvals.append(y)
        k1=f(t,x,y,gamma, a, b)
        l1=g(t,x,y,gamma, a, b)
        k2=f(t+dt/2 , x+dt*k1/2, y+dt*l1/2,gamma, a, b)
        l2=g(t+dt/2 , x+dt*k1/2, y+dt*l1/2,gamma, a, b)
        k3=f(t+dt/2 , x+dt*k2/2, y+dt*l2/2,gamma, a, b)
        l3=g(t+dt/2 , x+dt*k2/2, y+dt*l2/2,gamma, a, b)
        k4=f(t+dt, x+dt*k3, y+dt*l3,gamma, a, b)
        l4=g(t+dt, x+dt*k3, y+dt*l3,gamma, a, b)
        x+=dt/6 * (k1+2*k2+2*k3+k4)
        y+=dt/6 * (l1+2*l2+2*l3+l4)
     #   print("M: ", y, " || p: ", x, "\n")
    if graph:
        sol=[tvals,xvals,yvals]
        plt.rcParams['text.usetex'] = True
        fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(15,8))
        fig.suptitle('White Dwarf Equations of State, $R$='+str(round(tvals[-1]))+" km, $M/M_{\odot}=$"+str(round((yvals[-1]*1000))/1000), fontsize=18)
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
    sol=[tvals,xvals,yvals,t]
    return sol 

#Reccomended Parameters: a=1.473, b=52.46, gamma=4/3, p0= 10**(-16), dr=1 

def pressureChange(t,p,M,gamma, a, b):
    r=t
    return -(a*M*p**(1/gamma))/r**2

def massChange(t,p,M,gamma, a, b):
     r=t
     return b*(r**2)*(p**(1/gamma))

def runSimulation(a:float,b:float,gamma:float,epsilon0:float, p0, dr, customAB:bool=False, A:float=1, Z:float=2,graph:bool=True):
    if customAB: #This overrides any a and b selection
        K=(hbar*c)/(12*np.pi**2) * ( (Z * 3*np.pi**2) / (mN * A * c**2) )**(4/3)
        a=Rsun/( (K * epsilon0**(gamma -1))**(1/gamma) )
        b= (4 * np.pi * epsilon0) / (c **2 * (K * epsilon0**(gamma -1))**(1/gamma))
    return WDrk4(pressureChange,massChange,p0,0,dr,epsilon0, gamma, a, b, "White Dwarf",10**(-12),graph)

dataset=runSimulation(1.473, 52.46, 4/3, 10**(-12), 10**(-16), .1, False, 1,2,False)
def MassRadius(a,b,p0,p0max,dp=3*10**(-16)):
    masses=[]
    radii=[]
    i=1
    sol=[masses,radii]
    while p0<p0max:
        name="White Dwarf"+str(i) 
        print("Simulating White Dwarf No."+str(i), "|| p0="+str(p0))
        star=WDrk4(pressureChange,massChange,p0,0,10,10**(-12), 4/3, a, b, name,10**(-12),False)
        radius=star[-1]
        mass=star[-2][-1]
        print("r=", radius, " m=", mass)
        radii.append(radius)
        masses.append(mass)
        i+=1  
        p0+=dp        
    plt.plot(radii,masses,color="black")
    plt.title("Mass vs Radius of White Dwarves")
    plt.xlabel("Radii (km)")
    plt.ylabel("Mass in Solar Masses (1)")
    plt.show()
    return sol



def terminalSimulation():
    print(100*"=")
    preset=int(input("Welcome to Newtonian Equation of State Solver. Enter 0 to view preset solution or insert any other number to solve from custom parameters.\nA unit system in km and solar masses will be used.\nYour input: ")) 
    print(100*"=")
    if preset==0:
        graph=False
        MassRadius(1.473, 0.052, 10**(-5),100, 0.1)
    else:
        print("Enter the parameters that characterise your sample of stars. See documentation for common parameters for red giants, neutron stars, and more.")
        p0min=float(input("Enter minimum central normalised pressure: "))
        p0max=float(input("Enter minimum central normalised pressure: "))
        a=float(input("Enter coefficient of pressure equation (parameter a): "))
        b=float(input("Enter coefficient of mass equation (parameter b): "))
        dp=float(input("Enter pressure step size (parameter dp): "))
    return MassRadius(a,b,p0min,p0max,dp)


#print(MassRadius(1.473,52.46,+0.1*10**(-16)-0*10**(-12),10**(-16)+10**(-13.5), 10**(-18)))
#print(dataset[-2])