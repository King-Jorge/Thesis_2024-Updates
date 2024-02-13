# Thesis_2024-Updates
This repository contains all project files and updates.

# File Update 4
In this file, Problem 2 (Travel restrictions only model) and Problem 3 (Mixed policy model) are presented. 
We applied the PMP to these problems and proved that the optimal solution to the models are bang-bang. One difficulty encountered here is determing the switching times of the controls.

We observe some nonlinearity in the system, hence analytic solutions becomes difficult. On the other hand, we have litte information on the adjoint varibles which form the equation of the switching function.
Some researchers tackles this kind of problem by suggesting some numerical schemes [Switching Time Compuation (STC), Switching-Time-Variation Method (SVTM) and Time Optimal Switching (TOS)] that can solve for the switching times. 
However, these methods are simply applied to linear systems, hence applying to nonlinear systems requires a lot of modications (linearization and parameterization) and assumptions, another work on it own.

The approach used by Hansen and Day takes the form of guessing some solutions and then testing if those solutions are optimal.
But there is little explanation to the whole approach and little is known about the adjoint variables in their model. Hence it is not easily understood.

# File Update 3:
This file contains the corrections from the previous files (files 1 and 2). Some sections in files 1 and 2 are still going through corrections and are not included in this update.

# File Update 2:
This document is a continuation of file 1. I have included the formulation of Problem 2 (The travel restrictions-only model) and Problem 3 (The mixed policy model).
The theorems for these problems follow that in Hansen and Day. The proofs are different due to the model modification and assumptions. 
However, I followed the procedure Hansen and Day used to achieve a conclusion on the proofs. 


# File Update 1:
This file contains the model formulations and PMP statements. I am not able to provide detailed proof of the PMP due to the complexity of proving it. 
Notwithstanding, I have provided an elementary theorem and proof of the PMP. Reference is given to the original work of Pontryagin for the statement and proof of the theorem.
The general problem of the control model is provided. Theorems and proofs of Problem 1 - (Isolation-only model) are included in this file.


