# Thermodynamic study of single-level Quantum Dots in out-or-equilibrium regimes
 This prioject is an assembly of some sections of my bachelor thesis on simulation of the thermodynamic evolution of quantum doots.

 Three cases are treated in two regimes
* Driven regime - single Quantum Dot in contact with one reservoir
* Stationary regime - sinlge Quantum Dot in contact with two reservoirs
* Sationary regime - two coupled Quantum Dots in contact with three reservoirs

The project is here for informative porposes, you can use the information you get from the code and from the explanations. However, keep in mind that here you see some of the sections of my actual thesis and because some of the mathematical formalism illustrated comes purely from me, I would apprechiate to be refenced by github link if you choose to use the information in your projects. Moreover, I do not allow a copy and pasted use of sections of the project without behing referenced.

# Theory
The current directory and the subdirectories are provided with explanations needed to understand the physics of the problems and to be able to reproduce the simulations.

In the current directory, you find a theoretical introduction to stochastic jump processes, stochastic themrodynamics and a explanation on how quantum dots are used in this project, with the corresponding assumptions.

In the subdirectories, for each simulation, you find an explainations on how the simulation is done in the particular case and what is the purpose of it. For each case you also find an analysis of results and the corresponding accordance with the stochastic-thermodynamic theory explained in the introduction_to_Stochastic-Thermodynamic_and_simulaions.pdf file.

# Possible Improvements
Of course this project is not perfect and has some issues. That I report here below.

## Descrete-time simulation
* Some issues occur in the C program, I get a memory error when I try to do a power calculation in the main function, but not if I do it inside an external function that is then called by the main. This is independent if I use the C pow() function or if I use a costume one. The problem is bypassed by calculating powers on external functions but this create a useless repetition of code.
* The C program is not able to recognise the type -t input by command line if compiled and executed from a Mac OS machine. Need to change the string characterising the type of energy protocol from the code. However, it works perfectly fine in Linux.

## Gillespie simulation
The program is written in Python and becomes exponentially slower when the increasing the number of iterations while it would be nice to see if the oscillation in the occupation probability decrease with more trajectories. Using a compiled language or maybe simply coding using more python libraries we would achieve bigger numbers of iterations.

I would apprichiate to have some feedback/comment and/or some advice on how to solve the issues or maybe on solutions for more general related problems.
