# Latent3d
Julia module for latent symmetry in systems of high contrast scattering resonators

## Usage:
This module is intended to be downloaded and used with the Julia REPL.
To do this, ensure that the Julia packages `LinearAlgebra, Graphs, Polyhedra, GLPK, PlotlyJS` have been installed with `pkg`. Then, working in the directory containing `latent3d.jl` and `nautyinterface.so`, type `include("latent3d.jl"); import Latent3D` and you will have access to the data structures and functions that the module provides.

## Possible issues:
This module should work on any 64 bit machine. However on a 32 bit machine (or anything else, if you find yourself on something older or weirder) will require you to change the `WORDSIZE` constant and build `nautyinterface.so` and hence also build `nauty` with the `-fPIC` flag set in the makefile.

## Acknowledgements
This module was produced as a part of an undergraduate reaserach opportunity at Imperial College London and funded by an EPSRC vacation internship. 
I would like to thank Bryn Davies for supervising the project and Malte Röntgen for useful discussions.
