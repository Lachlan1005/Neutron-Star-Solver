import matplotlib.pyplot as plt
import numpy as np
hbar = 3.33 * 10**(-67)   # solar mass * km^2 / s
mN   = 8.395 * 10**(-58)  # solar masses
c    = 3 * 10**5          # km/s
Rsun = 696340             # km


def Equation_of_state(a:float, b:float, dr=10, customAB=False, p=10**(-16), M=0, r=0, gamma=4/3, epsilon0=10**(-12)):
    rList=[0]
    pList=[p]
    MList=[M]
    i=1
    plt.rcParams['text.usetex'] = True
    if customAB:
        Z=2
        A=1
        K=(hbar*c)/(12*np.pi**2) * ( (Z * 3*np.pi**2) / (mN * A * c**2) )**(4/3)
        a=Rsun/( (K * epsilon0**(gamma -1))**(1/gamma) )
        b= (4 * np.pi * epsilon0) / (c **2 * (K * epsilon0**(gamma -1))**(1/gamma))
    while p>0:
        print("iteration:", i, " r:", rList[-1], " M:", MList[-1], " p:", p)
        i+=1
        r+=dr
        M=((b* r**2 * p**(1/gamma))*dr+ M)
        p= (p- ( (a * M * p**(1/gamma))/( r**2) )*dr)
        rList.append(r)
        pList.append(p)
        MList.append(M)
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(15,8))
    fig.suptitle('White Dwarf Equations of State, $R$='+str(round(rList[-1]))+" km, $M/M_{\odot}=$"+str(round((MList[-1]*100))/100), fontsize=18)
    fig.text(0.5, 0.92, r"$\epsilon_0=$"+str(epsilon0)+r" $M_{\odot}c^2/km^3,\ p_0/\epsilon_0=$"+str(pList[0]),
    fontsize=14, ha='center')
    ax1.plot(rList, pList, color="black")
    ax1.set_title('Normalised Pressure vs Radial Distance', fontsize=16)
    ax1.set_xlabel('$r$ (km)', fontsize=16)
    ax1.set_ylabel(r'$p/\epsilon_0$ (1)', fontsize=16)

    ax2.plot(rList, MList, color="black")
    ax2.set_title("Mass vs Radial Distance", fontsize=16)
    ax2.set_xlabel('$r$ (km)', fontsize=16)
    ax2.set_ylabel(r'$M/M_{\odot}$ (1)', color='black', fontsize=16)

    plt.show()
    return r


#Equation_of_state(.013,.005924)

