# Pseudo Particle Strong Coupling (ppsc)

**Copyright (C) 2017- by Hugo U.R. Strand**

This code is not public and not licensed under any open source license (yet).

Please do not distribute.


Author: Hugo U.R. Strand


(Possible incomplete) List of contributors:

Martin Eckstein, Philipp Werner, Denis Golez, Nikolai Bittner, ...


## Purpose

Aimed at solving real-time quantum-mechanical impurity problems using the strong coupling approach using the pseudo particle formalism.

OpenMP parallelized calculation of first and second order strong coupling diagrams, a.k.a. the Non-Crossing Approximation (NCA) and the One-Crossing Approximation (OCA).


## Theory

For an introduction please see: M. Eckstein, P. Werner, Phys. Rev. B 82, 115115 (2010) <https://doi.org/10.1103/PhysRevB.82.115115>


## Requirements

- C++, MPI, OpenMP
- hdf5, boost?
- `libcntr` from NESSi (<https://github.com/nessi-cntr/nessi>)


---------------------------------------------------------------------------------------------------
(Update by Eva)

## How to set up/compile

- create 'configure.sh' file and 'cbuild' directory in the same parent directory
- navigate to 'cbuild' directory and execute shell script using the shell command 'sh ../configure.sh'
- to compile all type 'make' and enter; compile only certain programs use the command 'make <program>' instead