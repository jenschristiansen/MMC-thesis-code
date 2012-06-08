
#---------------------------------------------#

Code for Master's thesis: 

Topology optimization of flow domains 
using least-squares finite element methods

by Jens V. Christiansen

Technical University of Denmark

June, 2012

#---------------------------------------------#

The code is written in Python and uses FEniCS, see www.fenicsproject.org.

It is divided into two main parts:

1) Convergence analysis
2) Topology optimization

Part 2) makes use of the MMA algorithm. The algorithm is included in the repository as ksmma.so. The source code for the MMA algorithm was obtained on an academic license from Krister Svanberg, and is not included in this repository. ksmma.so was compiled in Ubuntu, and runs fine in version 11.10 and 12.04.


### Convergence analysis ###

The Python code setup for the convergence analysis in Chapter 3 is described here.

The 2D example is called CY1 and the 3D example is called CY2.

The convergence analysis was done using three scripts for each dimension. The setup for the 2D analysis is the following (the 3D case is analogous):

- CY1conv.py is the main script where the desired methods and grid dimensions are chosen. The results are saved to disk when done.
- CY1solver.py contains the actual FEniCS solver and is imported by CY1conv.py. This script imports the DOLFIN library.
- CY1plot.py loads the results saved by CY1conv.py and writes output to screen / plots the results.


### Topology optimization ###

The topology optimization problems in Chapter 4 were solved using two scripts for each problem: BP1,...,BP4. Using BP1 as an example:

BP1.py is the main script. It 
- contains the optimization loop;
- imports BP1import.py and calls functions from there;
- writes results to a log, saves them to disk, and plots them.

BP1import.py has multiple functions:
- it imports the DOLFIN} library;
- it imports the MMA algorithm from ksmma.so, which is compiled Fortran code;
- it contains the main 'solver' function, which is called by BP1.py in each iteration of the optimization process.

