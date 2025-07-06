# Polytropic Neutron Star Solver
This project provides the tools to solve various polytropic models for neutron stars (main focus) and other celestial objects such as White Dwarves based on the Tolman-Oppenheimer-Volkov (TOV) Equations.
There is also tools for modelling using the classical Newtonian equations. Through this solver, you can obtain accurate solutions for 

- The structure equations 
- The equation of states (EoS)

For the object of interest, as long as you have the correct initial conditions for it. The units are in kilometers and solar masses. Pressure and several other quantities* are expressed in 
normalised units by the parameter $$\epsilon_0$$. More details about the choice of normalisation parameter and the choice of units can be found in "Neutron Stars for Undergraduates" by Reddy and Silbar - the paper which this 
project is based upon. To get started, run Simulation_Wizard.py and follow the prompts. For mac users, use 

```
python3 Simulation_Wizard.py
```

To run the script. The source code is written in several sections, explained below. 

*: The full list of normalised quantities used can be found in the paper

## Simulation_Wizard.py
Provides a terminal interface to run the code in this project. Run this script to get started. 

## Full_TOV.py
This is the Runge-Kutta 4 (RK4) solver for the structure equations for one single star, using the full TOV equations.

## White_Dwarf_RK4.py
This is the RK4 solver for the structure equations for one single star, using the Newtonian equations.

## White_Dwarf_RvsM.py
This is the EoS solver using the Newtonian model. Solves for the masses and radii of an array of stars using the Newtonian structure equations, varying central pressure.

## Neutron_Star_RvsM.py
This is the EoS solver using the TOV model. Solves for the masses and radii of an array of stars using the TOV structure equations, varying central pressure.

## Fermi_gas.py
This is an Euler solver for the Fermi gas model for a polytrope. Obselete in any precise modelling. Diverges quickly and often yields null results. Inaccessable by the simulation wizard. 

## Notes 
To view the fully rendered notes, open "derivations.pdf". These proofs are the derivations to the various equations presented in the "Neutron Stars for Undergraduates" paper mentioned above. 
