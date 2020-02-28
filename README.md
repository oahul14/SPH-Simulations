# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Compilation/Installation Guide

You need to install Cmake first:
    sudo apt-get -y install cmake

Clone this repositiory and create a folder called 'build'. Navigate to that folder 'cd build'. On the command line type 'cmake ..'. This will create c CMakeLists.txt file that contains the instructions for your compiler. 
Then type 'make' and the program will be compiled.


### User instructions

To execute the program the commands are : 
    ```
    SPH_main tmax interval-for-animation type-of-timestep
    ```
tmax: maximum time that the simulation will run 
interval-for-animation: every how many seconds the program will output files
type-of-timestep: Options are: forward_euler, improved_euler, AB2. 
AB2 is the most accurate and forward_euler is the fastest.

The files that we have output can be viewed using ParaView, a free program. Simply import them and press apply. A graph will be made.


### Testing


We have also added C++ unit tests that test the individual functions. To compile them navigate to the 'tests' subdirectory and type 'make'.
To execute them type './tests/test_SPH'. These tests are assert statements so if any of them fail, the exucutable will not be executed.
If all tests pass,  the following message will appear: 'All tests passed succesfully ! '
You can also type ctest and this will test the file writting and out timestepping.
Known issue: the python_tests will fail, but this is not a cause of concern because we are not using them.
