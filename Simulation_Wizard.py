import Full_TOV 
import Neutron_Star_RvsM
import White_Dwarf_RK4
import White_Dwarf_RvsM

def terminal_Wizard():
    print(100*"=")
    print("Welcome to the Polytrope Modelling Program. \nA unit system in kilometers and solar masses will be used. Pressure and some other quantities will be dimensionless and scaled by the input variable epsilon0.\nSee README for more information.")
    print(100*"=")
    model=int(input("Enter 1 to use the Newtonian model or enter 2 to use the TOV model: "))
    if model==1:
        equation=int(input("Enter 1 to solve the structure equations once or enter 2 to solve fo the the full Equation of State: "))
        if equation==1:
            return White_Dwarf_RK4.terminalSimulation()
        elif equation==2:
            return White_Dwarf_RvsM.terminalSimulation()
        else:
            print("Error: Invalid input")
    if model==2:
        equation=int(input("Enter 1 to solve the structure equations once or enter 2 to solve fo the the full Equation of State: "))
        if equation==1:
            return Full_TOV.terminalSimulation()
        elif equation==2:
            return Neutron_Star_RvsM.terminalSimulation()
        else:
            print("Error: Invalid input")
    else:
        print("Error: Invalid input")

terminal_Wizard()