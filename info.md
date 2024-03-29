@mainpage Code-guide for phase-field solver

The solver has been successfully tested using OpenFOAM v8.

@section s1 Compiling the solver

* Following commands should create the executable of the solver
> cd ~/MicroSim/Grand_potential_OpenFOAM/phaseField

> wclean

> wmake

* The solver can be run by following the instructions in *userGuide*.

@section s2 Further details

The implementation, client and header files of the solver have been written following OpenFOAM conventions. These will be explained with flow charts generated from the source code using Doxygen. It must be noted that the solver is based on [laplacianFoam](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-basic-laplacianFoam.html "laplacianFoam") solver within OpenFOAM. Hence, it may be helpful for the user to become familiar with [OpenFOAM Programmer’s Guide](http://foam.sourceforge.net/docs/Guides-a4/ProgrammersGuide.pdf "OpenFOAM Programmer’s Guide") and [laplacianFoam](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-basic-laplacianFoam.html "laplacianFoam") beforehand.
